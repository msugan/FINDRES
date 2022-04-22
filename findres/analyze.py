from collections import defaultdict

import numpy as np
from multitaper.mtcross import MTCross
from obspy.signal.cross_correlation import correlate_template
from sklearn.linear_model import LinearRegression

from .plot import cross_spectrum as plot_cross_spectrum
from .utils import trim_window, sync_traces


class EmptyOrConstantTraceError(BaseException):
    pass


def correlate_waves(trace1, trace2, max_shift, **kwargs):
    if trace1.stats.npts == 0 or trace1.stats.npts == 0:
        raise EmptyOrConstantTraceError("Cannot cross-correlate empty traces.")
    elif trace1.std() == 0 or trace2.std() == 0:
        raise EmptyOrConstantTraceError("Normalized cross-correlation of constant traces is undefined.")
    data1, data2 = (trace1.data, trace2.data) if trace1.stats.npts >= trace2.stats.npts else (trace2.data, trace1.data)
    data1 = np.concatenate([np.zeros(max_shift), data1, np.zeros(max_shift)])
    cc = correlate_template(data1, data2, **kwargs)
    argmax = np.argmax(cc)
    shift = argmax - max_shift
    return shift, cc[argmax]


def relative_pick_time(stats1, stats2, pick1, pick2, shift):
    pick_time1 = pick1 - stats1.starttime
    pick_time2 = pick2 - stats2.starttime + shift * stats2.delta
    return (pick_time1 + pick_time2) / 2


def cc_preprocess(trace, pick_p_mean_delay, pick_s_mean_delay, freq_range, full_waveform_window, taper_length=0.15):
    new_trace = trace.copy()
    starttime = trace.stats.starttime + pick_p_mean_delay - full_waveform_window[0]
    endtime = trace.stats.starttime + pick_s_mean_delay + full_waveform_window[1]
    new_trace.trim(starttime=starttime, endtime=endtime)
    new_trace.detrend("constant")
    new_trace.filter("bandpass", freqmin=freq_range[0], freqmax=freq_range[1], corners=2, zerophase=True)
    new_trace.taper(taper_length)
    return new_trace


def cross_spectrum(trace_1, trace_2, delay, time_window, frequency_range, max_shift, cs_params,
                   graphics_filename):
    trace_1_trimmed = trace_1.copy()
    trim_window(trace_1_trimmed, delay, time_window)
    trace_2_trimmed = trace_2.copy()
    trim_window(trace_2_trimmed, delay, time_window)

    cc_shift, _ = correlate_waves(trace_1_trimmed, trace_2_trimmed, max_shift)
    if cc_shift != 0:
        sync_traces(trace_1_trimmed, trace_2_trimmed, cc_shift)
    slope_p, at, bt, ct, out = get_slope(trace_1_trimmed, trace_2_trimmed, frequency_range, cs_params)
    if slope_p is None:
        return None
    else:
        time_delay = cc_shift * trace_1_trimmed.stats.delta + slope_p / (2 * np.pi)
        if graphics_filename:
            plot_cross_spectrum(trace_1_trimmed, trace_2_trimmed, at, bt, ct, slope_p, time_delay, out,
                                graphics_filename)
        return time_delay


def get_slope(trace1, trace2, frange, params):
    freq_range = (frange[0], frange[1] + params['coherence_frequency_upper_bound_correction'])
    out = MTCross(trace1.data, trace2.data,
                  dt=trace1.stats.delta,
                  nw=params['mt_coherence_time_bandwidth_product'],
                  kspec=params['mt_coherence_num_taper'],
                  nfft=int((trace1.stats.npts // 2) + 1), iadapt=1)
    coherence = out.cohe.flatten()
    frequency = out.freq.flatten()
    phase = out.phase.flatten()
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
