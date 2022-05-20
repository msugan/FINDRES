import logging
from itertools import combinations
from math import cos, hypot

import numpy as np
import obspy
import pandas as pd
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth

from obspy.signal.trigger import pk_baer

from . import custom_formats
from .utils import estimate_s_pick, p_picker


class InventoryLookupError(BaseException):
    pass


class MissingPhaseDataError(BaseException):
    pass


class SPickEstimationError(BaseException):
    pass


def zmap(catalog_path, extensions=None):
    names = ['longitude',
             'latitude',
             'year',
             'month',
             'day',
             'magnitude',
             'depth',
             'hour',
             'minute',
             'second']

    if extensions:
        names.extend(extensions)

    catalogue = pd.read_csv(catalog_path,
                            sep=r'\s+',
                            usecols=range(len(names)),
                            names=names,
                            parse_dates={'date': ['year', 'month', 'day', 'hour', 'minute', 'second']},
                            date_parser=lambda datestr: pd.to_datetime(datestr, format='%Y %m %d %H %M %S.%f',
                                                                       utc=True))
    return catalogue


# TODO: improve performance, reading of phase_file should be cached (also, only hypoinv has been thoroughly tested)
def errors(catalogue, phase_file, phase_type):
    if phase_type == 'hypoinv':
        error_function = _hypoinverse_errors
    elif phase_type == 'hypoel':
        error_function = _hypoel_errors
    elif phase_type in ('quakeml', 'nll'):
        error_function = _obspy_errors
    elif phase_type in custom_formats.registry:
        error_function = custom_formats.registry[phase_type]['errors']
    else:
        raise NotImplementedError()

    return {event.name: error_function(event.date, phase_file) for event in catalogue.itertuples()}


def _textfile_errors(origin_time, phase_file, columns_delimiters):
    error_values = {}
    with phase_file.open('r') as file:
        while line := file.readline():
            if line.startswith(origin_time.strftime('%Y%m%d%H%M%S')):
                for key, (a, b) in columns_delimiters.items():
                    try:
                        error_values[key] = _read_float(line, a, b)
                    except ValueError as err:
                        logging.warning("Found invalid string.", exc_info=err)
                        error_values[key] = None
    return error_values


# Column numbers taken from the documentation of the hypoinverse/hypoel formats
def _hypoinverse_errors(origin_time, phase_file):
    return _textfile_errors(origin_time, phase_file, {'time_uncertainty': (48, 52),
                                                      'horizontal_uncertainty': (85, 89),
                                                      'vertical_uncertainty': (89, 93)})


def _hypoel_errors(origin_time, phase_file):
    return _textfile_errors(origin_time, phase_file, {'time_uncertainty': (47, 50),
                                                      'horizontal_uncertainty': (56, 59),
                                                      'vertical_uncertainty': (74, 77)})


def _obspy_errors(origin_time, phase_file):
    catalogue = obspy.read_events(str(phase_file))
    error_values = {}
    for event in catalogue:
        if abs(event.time - obspy.UTCDateTime(origin_time)) < 1.0:
            origin, = event.origins
            error_values['time_uncertainty'] = origin.time_errors.uncertainty
            error_values['horizontal_uncertainty'] = _horizontal_error(origin)
            error_values['vertical_uncertainty'] = origin.depth_errors.uncertainty
    return error_values


def _horizontal_error(origin, earth_radius_km=6.371e3):
    lat = origin.latitude
    sigma_lat = origin.latitude_errors.uncertainty
    sigma_lon = origin.longitude_errors.uncertainty
    return (np.pi / 180.0) * earth_radius_km * hypot(sigma_lat, sigma_lon * cos(lat))


def _read_float(whole_str, start, end, decimals=2):
    if whole_str[start:end].isspace():
        return None
    float_str = whole_str[start:end - decimals] + '.' + whole_str[end - decimals:end]
    float_str = float_str.replace(' ', '0')
    return float(float_str)


def coordinates(inventory, network_code, station_code, time, tol=0.2):
    stations = []
    for network in inventory.select(network=network_code, station=station_code, time=time):
        for station in network:
            stations.append(station)
    if stations:
        if len(stations) > 1:
            maxdist = 1e-3 * np.max([gps2dist_azimuth(s.latitude, s.longitude, r.latitude, r.longitude)[0]
                                     for s, r in combinations(stations, 2)])
            if maxdist > tol:
                raise InventoryLookupError(f"Incompatible matches for {station_code} at {time} found in inventory.")
        return stations[0].latitude, stations[0].longitude
    else:
        raise InventoryLookupError(f"Station {station_code} at {time} cannot be found.")


def picks(event, event_coordinates, station_coordinates, trace: obspy.Trace, params, phase_file=None, phase_type=None, picker=False,
          picker_arguments=None, travel_times_function=None):
    station_latitude, station_longitude = station_coordinates

    if phase_file is None:
        if picker:
            if picker_arguments is None:
                raise MissingPhaseDataError("Picker arguments must be provided")
            p_pick = p_picker(trace, params)
            s_pick = None
        elif travel_times_function is not None:
            p_pick, s_pick = compute_picks(event, station_latitude, station_longitude, travel_times_function,
                                           p_model_corr=params['p_velocity_model_correction'],
                                           s_model_corr=params['s_velocity_model_correction'])
        else:
            raise MissingPhaseDataError("Either a phase file, picker or a model must be provided")
    else:
        p_pick, s_pick = _read_picks(event.date, trace.stats.station, phase_file, phase_type)
        if p_pick is None:
            if picker:
                if picker_arguments is None:
                    raise MissingPhaseDataError("Picker arguments must be provided")
                p_pick = p_picker(trace, params)
                s_pick = None
            elif travel_times_function is not None:
                p_pick, s_pick = compute_picks(event, station_latitude, station_longitude, travel_times_function,
                                               p_model_corr=params['p_velocity_model_correction'],
                                               s_model_corr=params['s_velocity_model_correction'])
            else:
                raise MissingPhaseDataError("Phase file with picking for P phase, picker or a model must be provided")
        if s_pick is None:
            event_latitude, event_longitude = event_coordinates
            epi_dist, _, _ = gps2dist_azimuth(event_latitude, event_longitude, station_latitude, station_longitude)
            try:
                s_pick = estimate_s_pick(trace, p_pick, epi_dist, params)
            except ValueError as err:
                raise SPickEstimationError(str(err))
    return p_pick, s_pick


# TODO: improve performance, phase_file should be cached (also, only hypoinv has been thoroughly tested)
def _read_picks(origin_time, station, phase_file, phase_type):
    if phase_type == 'hypoinv':
        phases_function = _hypoinverse_phases
    elif phase_type == 'hypoel':
        phases_function = _hypoel_phases
    elif phase_type in ('nll', 'quakeml'):
        phases_function = _ob_phases
    elif phase_type in custom_formats.registry:
        phases_function = custom_formats.registry[phase_type]['phases']
    else:
        raise NotImplementedError()

    return phases_function(origin_time, station, phase_file)


def _hypoinverse_phases(origin_time, station, phase_file):
    with phase_file.open('r') as file:
        lines = []
        start_section = False
        while line := file.readline():
            if line.startswith(origin_time.strftime('%Y%m%d%H%M%S')):
                start_section = True
                continue
            if start_section and line.startswith(30 * " "):
                break
            if start_section:
                lines.append(line)
    p_pick, s_pick = None, None
    for line in lines:
        if line.startswith(station) and int(line[16]) < 4:
            if not line[30:34].isspace():
                p_pick = UTCDateTime.strptime(line[17:34], '%Y%m%d%H%M %S%f')
            elif not line[42:46].isspace():
                s_pick = UTCDateTime.strptime(line[17:30] + line[42:46], '%Y%m%d%H%M %S%f')
    return p_pick, s_pick


def _hypoel_phases(origin_time, station, phase_file):
    with phase_file.open('r') as file:
        lines = []
        start_section = False
        while line := file.readline():
            if line.startswith(origin_time.strftime('%Y%m%d%H%M%S')):
                start_section = True
                continue
            if start_section and line.startswith(17 * " "):
                break
            if start_section:
                lines.append(line)
    p_pick, s_pick = None, None
    for line in lines:
        if line.startswith(station) and int(line[7]) < 4:
            if not line[20:24].isspace():
                p_pick = UTCDateTime.strptime(line[9:19] + line[20:24], '%y%m%d%H%M%S%f')
            if not line[32:36].isspace():
                s_pick = UTCDateTime.strptime(line[9:19] + line[32:36], '%y%m%d%H%M%S%f')
    return p_pick, s_pick


def _ob_phases(origin_time, station, phase_file):
    catalogue = obspy.read_events(str(phase_file))
    p_pick, s_pick = None, None
    for event in catalogue:
        if abs(event.time - obspy.UTCDateTime(origin_time)) < 1.0:
            for pick in event.picks:
                if pick.waveform_id.station_code == station and pick.phase_hint == 'P' \
                        and pick.time_errors.uncertainty < 1.0:
                    p_pick = pick.time
                if pick.waveform_id.station_code == station and pick.phase_hint == 'S' \
                        and pick.time_errors.uncertainty < 1.0:
                    s_pick = pick.time
    return p_pick, s_pick


def compute_picks(event, station_latitude, station_longitude, model_function, p_model_corr=None, s_model_corr=None,
                  earth_radius=6.371e6, min_event_depth=1.5e3):
    origin_time = UTCDateTime(event.date)
    event_latitude = event.latitude
    event_longitude = event.longitude
    event_depth = event.depth
    event_depth = max(event_depth, min_event_depth)

    epi_dist, _, _ = gps2dist_azimuth(event_latitude, event_longitude, station_latitude, station_longitude)
    deg = (180.0 / np.pi) * epi_dist / earth_radius

    arrivals_p = model_function(source_depth_in_km=event_depth * 1e-3, distance_in_degree=deg, phase_list=("p", "P"))
    arrivals_s = model_function(source_depth_in_km=event_depth * 1e-3, distance_in_degree=deg, phase_list=("s", "S"))
    pick_p = origin_time + arrivals_p[0].time
    pick_s = origin_time + arrivals_s[0].time
    if p_pick_correction := p_model_corr:
        pick_p -= p_pick_correction
    if s_pick_correction := s_model_corr:
        pick_s -= s_pick_correction
    return pick_p, pick_s
