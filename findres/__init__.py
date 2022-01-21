import logging
from collections import defaultdict
from itertools import combinations
from math import cos, hypot

import matplotlib.pyplot as plt
import numpy as np
import obspy
import pandas as pd
from mtspec import mt_coherence
from obspy import Stream
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.cross_correlation import correlate
from obspy.signal.cross_correlation import xcorr_max
from obspy.signal.filter import envelope
from sklearn.linear_model import LinearRegression


def read_zmap(catalog_path, extensions=None):
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

    zmap = pd.read_csv(catalog_path,
                       sep=r'\s+',
                       usecols=range(len(names)),
                       names=names,
                       parse_dates={'date': ['year', 'month', 'day', 'hour', 'minute', 'second']},
                       date_parser=lambda datestr: pd.to_datetime(datestr, format='%Y %m %d %H %M %S.%f', utc=True))

    return zmap


# TODO: improve performance, reading of phase_file should be cached (also, only hypoinv has been thoroughly tested)
def read_errors(catalogue, phase_file, phase_type):
    if phase_type == 'hypoinv':
        error_function = _hypoinverse_errors
    elif phase_type == 'hypoel':
        error_function = _hypoel_errors
    elif phase_type in ('quakeml', 'nll'):
        error_function = _obspy_errors
    else:
        raise NotImplementedError()

    return [error_function(event.date, phase_file) for event in catalogue.itertuples()]


def _textfile_errors(origin_time, phase_file, columns_delimiters):
    errors = {}
    with phase_file.open('r') as file:
        while line := file.readline():
            if line.startswith(origin_time.strftime('%Y%m%d%H%M%S')):
                for key, (a, b) in columns_delimiters.items():
                    errors[key] = _read_float(line, a, b)
    return errors


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
    errors = {}
    for event in catalogue:
        if abs(event.time - obspy.UTCDateTime(origin_time)) < 1.0:
            origin, = event.origins
            errors['time_uncertainty'] = origin.time_errors.uncertainty
            errors['horizontal_uncertainty'] = _horizontal_error(origin)
            errors['vertical_uncertainty'] = origin.depth_errors.uncertainty
    return errors


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


def piecewise_from_thresholds(x, thresholds):
    if thresholds:
        _, value = thresholds[0]
        for threshold, value in reversed(thresholds):
            if x > threshold:
                break
        return value
    else:
        raise RuntimeError()


def zip_streams(stream1, stream2):
    common_ids = set.intersection({trace.id for trace in stream1}, {trace.id for trace in stream2})
    return zip(Stream([trace for trace in stream1 if trace.id in common_ids]),
               Stream([trace for trace in stream2 if trace.id in common_ids]))


def get_coordinates(inventory, station_code):
    stations = []
    for network in inventory.select(station=station_code):
        for station in network:
            stations.append(station)
    if stations:
        logging.debug(f"More than one match for {station_code} found in inventory")
        return stations[0].latitude, stations[0].longitude
    else:
        raise RuntimeError(f"Station {station_code} cannot be found")


def get_picks(event, event_coordinates, station_coordinates, trace, params, phase_file=None, phase_type=None,
              travel_times_function=None):

    station_latitude, station_longitude = station_coordinates
    if phase_file is None:
        if travel_times_function is not None:
            p_pick, s_pick = _get_travel_times(event, station_latitude,
                                               station_longitude,
                                               travel_times_function)

            if p_pick_correction := params['p_velocity_correction']:
                p_pick -= p_pick_correction
            if s_pick_correction := params['s_velocity_correction']:
                s_pick -= s_pick_correction
        else:
            raise RuntimeError("Either a phase file or a model must be provided")
    else:
        p_pick, s_pick = _read_picks(event.date, trace.stats.station, phase_file, phase_type)
        if p_pick is None:
            raise RuntimeError("Phase file must provide the picking for P phase")
        if s_pick is None:
            event_latitude, event_longitude = event_coordinates
            epi_dist, _, _ = gps2dist_azimuth(event_latitude, event_longitude, station_latitude, station_longitude)
            s_pick = estimate_s_pick(trace, p_pick, epi_dist, params['Vp'], params['Vs'],
                                     params['s_pick_estimation_window'])
    return p_pick, s_pick


# TODO: improve performance, phase_file should be cached (also, only hypoinv has been thoroughly tested)
def _read_picks(origin_time, station, phase_file, phase_type):
    if phase_type == 'hypoinv':
        phases_function = _hypoinverse_phases
    elif phase_type == 'hypoel':
        phases_function = _hypoel_phases
    elif phase_type in ['nll', 'quakeml']:
        phases_function = _ob_phases
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


def _get_travel_times(event, station_latitude, station_longitude, model_function, earth_radius=6.371e6,
                      min_event_depth=1.5e3):
    origin_time = UTCDateTime(event.date)
    event_latitude = event.latitude
    event_longitude = event.longitude
    event_depth = event.depth
    event_depth = max(event_depth, min_event_depth)

    epi_dist, _, _ = gps2dist_azimuth(event_latitude, event_longitude, station_latitude, station_longitude)
    deg = (180.0 / np.pi) * epi_dist / earth_radius

    arrivals_p1 = model_function(source_depth_in_km=event_depth * 1e-3, distance_in_degree=deg, phase_list=("p", "P"))
    arrivals_s1 = model_function(source_depth_in_km=event_depth * 1e-3, distance_in_degree=deg, phase_list=("s", "S"))
    arr_p1 = origin_time + arrivals_p1[0].time
    arr_s1 = origin_time + arrivals_s1[0].time
    return arr_p1, arr_s1


def estimate_s_pick(trace, pick, distance, vp, vs, window, shift=1, freqmin=1, freqmax=10):
    trace_filt = trace.copy()
    trace_filt.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=2, zerophase=True)
    data_envelope = envelope(trace_filt.data)
    sample_shift = int(shift * trace.stats.sampling_rate)
    ind = np.argmax(data_envelope) - sample_shift

    pick_sample = int((pick - trace.stats.starttime) * trace.stats.sampling_rate)
    thumb_rule_delay_sample = int(distance * (vp - vs) / (vp * vs) * trace.stats.sampling_rate)
    reference_sample = pick_sample + thumb_rule_delay_sample

    delta1, delta2 = window
    start = reference_sample - delta1
    stop = reference_sample + delta2
    if not (pick_sample < ind and start <= ind <= stop):
        ind = start + np.argmax(data_envelope[start:stop]) - sample_shift

    return trace.stats.starttime + ind * trace.stats.delta


def correlate_waves(data1, data2, max_shift, **kwargs):
    cc = correlate(data1, data2, max_shift, method="direct", **kwargs)
    shift, value = xcorr_max(cc, abs_max=False)
    return shift, value


def relative_pick_time(stats1, stats2, pick1, pick2, shift):
    pick_time1 = pick1 - stats1.starttime
    pick_time2 = pick2 - stats2.starttime + shift * stats2.delta
    return (pick_time1 + pick_time2) / 2


def sync_traces(trace_1, trace_2, shift):
    trace_1.stats.update({'starttime': 0.0})
    trace_2.stats.update({'starttime': shift * trace_1.stats.delta})
    starttime = min(trace_1.stats.starttime, trace_2.stats.starttime)
    endtime = max(trace_1.stats.endtime, trace_2.stats.endtime)
    trace_1.trim(starttime=starttime, endtime=endtime, pad=True, fill_value=0.0)
    trace_2.trim(starttime=starttime, endtime=endtime, pad=True, fill_value=0.0)
    trace_1.stats.update({'starttime': 0.0})
    trace_2.stats.update({'starttime': 0.0})


def cc_preprocess(trace, pick_p_mean_delay, pick_s_mean_delay, freq_range, full_waveform_window, taper_length=0.15):
    new_trace = trace.copy()
    starttime = trace.stats.starttime + pick_p_mean_delay - full_waveform_window[0]
    endtime = trace.stats.starttime + pick_s_mean_delay + full_waveform_window[1]
    new_trace.trim(starttime=starttime, endtime=endtime)
    new_trace.detrend("constant")
    new_trace.filter("bandpass", freqmin=freq_range[0], freqmax=freq_range[1], corners=2, zerophase=True)
    new_trace.taper(taper_length)
    return new_trace


def normalize_cc(x, y, cc_value, shift):
    x_start = max(0, shift)
    x_end = min(len(x), len(y) + shift)
    n = x_end - x_start
    x = x[x_start:x_end]
    y_start = max(0, -shift)
    y_end = min(len(x), len(y) - shift)
    y = y[y_start:y_end]
    return (cc_value / n - np.mean(x) * np.mean(y)) / (np.std(x) * np.std(y))


def cross_spectrum_analysis(trace_1, trace_2, delay, time_window, frequency_range, max_shift, cs_params,
                            graphics_filename):
    trace_1_trimmed = trace_1.copy()
    trim_window(trace_1_trimmed, delay, time_window)
    trace_2_trimmed = trace_2.copy()
    trim_window(trace_2_trimmed, delay, time_window)

    cc_shift, _ = correlate_waves(trace_1_trimmed.data, trace_2_trimmed.data, max_shift)
    if cc_shift != 0:
        sync_traces(trace_1_trimmed, trace_2_trimmed, cc_shift)
    slope_p, at, bt, ct, out = get_slope(trace_1_trimmed, trace_2_trimmed, frequency_range, cs_params)
    if slope_p is None:
        raise RuntimeWarning(f"No coherence found")
    else:
        time_delay = cc_shift * trace_1_trimmed.stats.delta + slope_p / (2 * np.pi)
        if graphics_filename:
            plot_cross_spectrum(trace_1_trimmed, trace_2_trimmed, at, bt, ct, slope_p, time_delay, out,
                                graphics_filename)
        return time_delay


def trim_window(trace, shift, window):
    trace.trim(starttime=trace.stats.starttime + shift - window[0],
               endtime=trace.stats.starttime + shift + window[1],
               pad=True,
               fill_value=0.0)


def get_slope(trace1, trace2, frange, params):
    freq_range = (frange[0], frange[1] + params['coherence_frequency_upper_bound_correction'])
    out = mt_coherence(trace1.stats.delta, trace1.data, trace2.data,
                       params['mt_coherence_time_bandwidth_product'],
                       params['mt_coherence_num_taper'],
                       int((trace1.stats.npts // 2) + 1),
                       params['mt_coherence_confidence'],
                       freq=True, phase=True, cohe=True, iadapt=1)
    coherence = out['cohe']
    frequency = out['freq']
    phase = out['phase']
    index_a0, = np.nonzero(coherence > params['cs_coherence_threshold'])
    index_b, = np.nonzero((freq_range[0] < frequency) & (frequency < freq_range[1]))
    indices = np.intersect1d(index_a0, index_b)

    t = np.split(indices, np.where(np.diff(indices) != 1)[0] + 1)
    index_a_sel = max(t, key=lambda x: len(x))
    a_take = coherence[index_a_sel]
    b_take = frequency[index_a_sel]
    c_take = phase[index_a_sel] * np.pi / 180

    slope = None
    if len(index_a_sel) > params['coherence_minimum_num_points']:
        lr_fi_false = LinearRegression(fit_intercept=False)
        lr_fi_false.fit(b_take.reshape(-1, 1), c_take.reshape(-1, 1))
        slope = lr_fi_false.coef_[0, 0]
    return slope, a_take, b_take, c_take, out


def connected_components(edges):
    neighbors = defaultdict(set)
    for edge in edges:
        for node in edge:
            neighbors[node].update(edge)
    graph = set(neighbors.keys())
    visited = set()

    def component_from(starting_node):
        component = {starting_node}
        while component:
            current_node = component.pop()
            visited.add(current_node)
            component |= (neighbors[current_node] - visited)
            yield current_node

    components = []
    for node in graph:
        if node not in visited:
            components.append(sorted(component_from(node)))

    return components


def dump_hypodd(families, catalogue, errors, REs, parameters, output_dir):
    for i, family in enumerate(families):
        with (output_dir / f"{i}_RES_event.sel").open("w") as file:
            for n in family:
                e_t, e_h, e_v = map(lambda s: err if (err := errors[n].get(s)) else parameters['hypodd_default_' + s],
                                    ['time_uncertainty', 'horizontal_uncertainty', 'vertical_uncertainty'])
                file.write(f"{catalogue.loc[n, 'date'].strftime('%Y%m%d  %H%M%S%f')[:-4]}   "
                           f"{catalogue.loc[n, 'latitude']:.4f}     {catalogue.loc[n, 'longitude']:.4f}    "
                           f"{catalogue.loc[n, 'depth']:.3f}   {catalogue.loc[n, 'magnitude']:.1f}    "
                           f"{e_t:.2f}    {e_h:.2f}   {e_v:.2f}        "
                           f"{catalogue.loc[n, 'name']}\n")
        if len(family) > 2:
            with (output_dir / f"{i}_RES_dt.cc").open("w") as file:
                for (t1, t2) in combinations(family, 2):
                    if REs[(t1, t2)]:
                        file.write(f"#    {catalogue.loc[t1, 'name']}    {catalogue.loc[t2, 'name']}     0.0\n")
                        for station in sorted(REs[(t1, t2)]):
                            cc, delta_sp = REs[(t1, t2)][station]
                            delta_v = parameters['Vp'] - parameters['Vs']
                            ttp = parameters['Vs'] * delta_sp / delta_v
                            file.write(f"{station}     {ttp: 10.9f}    {cc:.2f}    P\n")
                            tts = -parameters['Vp'] * delta_sp / delta_v
                            file.write(f"{station}     {tts: 10.9f}    {cc:.2f}    S\n")


def plot_signals(trace1, trace2, mag_t1, mag_t2, p_wave_window_1, p_wave_window_2, s_wave_window_1, s_wave_window_2,
                 cc, frange,
                 filename):
    plt.figure(figsize=(6, 7))

    plt.subplot(211)
    plt.plot(np.arange(trace1.stats.npts) * trace1.stats.delta, trace1.data / abs(trace1.max()), linewidth=1,
             color="blue", alpha=0.5)
    plt.plot(np.arange(trace2.stats.npts) * trace2.stats.delta, trace2.data / abs(trace2.max()), linewidth=1,
             color="red", alpha=0.5)

    plt.xlabel("Times (sec)")
    plt.ylabel("Normalized Amplitude")

    plt.text(0.5, 1, f"M={mag_t1:.2f}-{mag_t2:.2f}     BP = {frange[0]}-{frange[1]}Hz     CC={cc:.2f}")
    plt.ylim(-1.3, 1.3)

    plt.hlines(-1, p_wave_window_1, p_wave_window_2, color="r", linewidth=1, linestyles="dashed")
    plt.hlines(-1, s_wave_window_1, s_wave_window_2, color="r", linewidth=1, linestyles="dashed")
    plt.savefig(filename)
    plt.close()


def plot_cross_spectrum(trace1, trace2, a_take, b_take, c_take, slope, time_delay, out, filename):
    plt.figure(figsize=(6, 7))

    plt.subplot(411)
    plt.plot(np.arange(trace1.stats.npts) * trace1.stats.delta, trace1.data / np.abs(trace1.max()), linewidth=1)
    plt.xlabel("Time [sec]")
    plt.ylabel("Amplitude")

    plt.subplot(412)
    plt.plot(np.arange(trace2.stats.npts) * trace2.stats.delta, trace2.data / np.abs(trace2.max()), linewidth=1)
    plt.xlabel("Time [sec]")
    plt.ylabel("Amplitude")

    plt.subplot(413)
    plt.plot(out["freq"], out["cohe"], linewidth=1)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Coherency")
    plt.ylim(0, 1.1)
    plt.plot(b_take, a_take, "rx", label="original data", markersize=4)

    plt.subplot(414)
    plt.plot(out["freq"], out['phase'] * np.pi / 180.0, linewidth=1)
    plt.plot(b_take, c_take, "rx", label="original data", markersize=4)
    plt.plot(b_take.reshape(-1, 1), slope * b_take.reshape(-1, 1), "g", label="fitted line")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Phase(rad)")
    plt.ylim(-3, 3)
    plt.tight_layout()
    plt.text(1, 1, f"time_delay = {time_delay:0.6f}")

    plt.savefig(filename)
    plt.close()


def plot_similarity(ev1, ev2, data, inventory, cc_threshold, delta_threshold, filename, **kwargs):
    event_latitude = (ev1.latitude + ev2.latitude) / 2
    event_longitude = (ev1.longitude + ev2.longitude) / 2
    plt.figure(figsize=(7, 5), **kwargs)
    xs = np.fromiter((x for x, _ in data[(ev1.Index, ev2.Index)].values()), dtype=np.float)
    ys = np.fromiter((y for _, y in data[(ev1.Index, ev2.Index)].values()), dtype=np.float)
    zs = np.fromiter(map(lambda s: gps2dist_azimuth(event_latitude, event_longitude, *get_coordinates(inventory, s))[0],
                         data[(ev1.Index, ev2.Index)]),
                     dtype=np.float)
    plt.scatter(xs, np.absolute(ys), c=1e-3 * zs, cmap="jet")
    plt.xlim(cc_threshold, 1.0)
    plt.ylim(0.0, delta_threshold)
    plt.xlabel("CC Z_component")
    plt.ylabel("|ΔS-P|")
    plt.title(f"{ev1.name}-{ev2.name}")
    cb = plt.colorbar()
    cb.set_label("epi station-event distance (km)")
    plt.text(cc_threshold * 1.01, delta_threshold * 0.95, f"M_{ev1.name}={ev1.magnitude}")
    plt.text(cc_threshold * 1.01, delta_threshold * 0.9, f"M_{ev2.name}={ev2.magnitude}")
    plt.savefig(filename)
    plt.close()
