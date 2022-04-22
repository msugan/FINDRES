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
