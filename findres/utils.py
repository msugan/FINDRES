from collections import defaultdict

from obspy import Stream
from obspy.signal.filter import envelope
import numpy as np


def piecewise_from_thresholds(x, thresholds):
    if thresholds:
        _, value = thresholds[0]
        for threshold, value in reversed(thresholds):
            if x >= threshold:
                break
        return value
    else:
        raise RuntimeError()


def zip_streams(stream1, stream2):
    common_ids = set.intersection({trace.id for trace in stream1}, {trace.id for trace in stream2})
    return list(zip(Stream([trace for trace in stream1 if trace.id in common_ids]),
                    Stream([trace for trace in stream2 if trace.id in common_ids])))


def relative_pick_time(stats1, stats2, pick1, pick2, shift):
    pick_time1 = pick1 - stats1.starttime
    pick_time2 = pick2 - stats2.starttime + shift * stats2.delta
    return (pick_time1 + pick_time2) / 2


def cc_preprocess(trace, pick_p_mean_delay, pick_s_mean_delay, freq_range, full_waveform_window, taper_length=0.15):
    new_trace = trace.copy()
    starttime = trace.stats.starttime + pick_p_mean_delay - full_waveform_window[0]
    endtime = trace.stats.starttime + pick_s_mean_delay + full_waveform_window[1]
    new_trace.detrend("constant")
    new_trace.filter("bandpass", freqmin=freq_range[0], freqmax=freq_range[1], corners=2, zerophase=True)
    new_trace.taper(taper_length)
    new_trace.trim(starttime=starttime, endtime=endtime)
    return new_trace


def sync_traces(trace_1, trace_2, shift):
    trace_1.stats.update({'starttime': 0.0})
    trace_2.stats.update({'starttime': shift * trace_1.stats.delta})
    starttime = min(trace_1.stats.starttime, trace_2.stats.starttime)
    endtime = max(trace_1.stats.endtime, trace_2.stats.endtime)
    trace_1.trim(starttime=starttime, endtime=endtime, pad=True, fill_value=0.0)
    trace_2.trim(starttime=starttime, endtime=endtime, pad=True, fill_value=0.0)
    trace_1.stats.update({'starttime': 0.0})
    trace_2.stats.update({'starttime': 0.0})


def trim_window(trace, center, window):
    trace.trim(starttime=trace.stats.starttime + center - window[0],
               endtime=trace.stats.starttime + center + window[1],
               pad=True,
               fill_value=0.0)


def estimate_s_pick(trace, pick, distance, vp, vs, window, shift=1, freqmin=1, freqmax=10):
    trace_filt = trace.copy()
    trace_filt.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=2, zerophase=True)
    data_envelope = envelope(trace_filt.data)
    sample_shift = int(shift * trace.stats.sampling_rate)

    pick_sample = int((pick - trace.stats.starttime) * trace.stats.sampling_rate)
    thumb_rule_delay_sample = int(distance * (vp - vs) / (vp * vs) * trace.stats.sampling_rate)
    reference_sample = pick_sample + thumb_rule_delay_sample

    ind = np.argmax(data_envelope) - sample_shift
    delta1, delta2 = window
    start = reference_sample - delta1
    stop = reference_sample + delta2
    if not (pick_sample < ind and start <= ind <= stop):
        ind = start + np.argmax(data_envelope[start:stop]) - sample_shift

    return trace.stats.starttime + ind * trace.stats.delta


def longest_subsequence_of_trues(bools):
    if len(bools) > 0:
        sequences = [None]
        for n, x in enumerate(bools):
            current_sequence = sequences[-1]
            if x:
                if current_sequence is None:
                    sequences[-1] = (n, n)
                else:
                    start, _ = current_sequence
                    sequences[-1] = (start, min(n + 1, len(bools)))
            else:
                if current_sequence is not None:
                    sequences.append(None)
        if sequences[-1] is None:
            sequences.pop()
        if sequences:
            return max(sequences, key=lambda pair: pair[1] - pair[0])
        else:
            return None
    else:
        return None


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
