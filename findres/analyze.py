import numpy as np
from obspy.signal.cross_correlation import correlate_template

from findres.utils import longest_subsequence_of_trues


class EmptyOrConstantTraceError(BaseException):
    pass


class LowCoherenceError(BaseException):
    pass


def correlate_waves(trace1, trace2, max_shift, **kwargs):
    if trace1.stats.npts == 0 or trace1.stats.npts == 0:
        raise EmptyOrConstantTraceError("Cannot cross-correlate empty traces.")
    elif trace1.std() == 0 or trace2.std() == 0:
        raise EmptyOrConstantTraceError("Normalized cross-correlation of constant traces is undefined.")
    if trace1.stats.npts >= trace2.stats.npts:
        is_swapped = False
        data1, data2 = trace1.data, trace2.data
    else:
        is_swapped = True
        data1, data2 = trace2.data, trace1.data
    data1 = np.concatenate([np.zeros(max_shift), data1, np.zeros(max_shift)])
    cc = correlate_template(data1, data2, **kwargs)
    argmax = np.argmax(cc)
    shift = argmax - max_shift
    if is_swapped:
        shift = -shift
    return shift, cc[argmax]


def find_slope(cross_spectrum, frequency_range, params):
    coherence = cross_spectrum.cohe.flatten()
    frequency = cross_spectrum.freq.flatten()
    phase = cross_spectrum.phase.flatten()
    lower_frequency_bound = frequency_range[0],
    upper_frequency_bound = frequency_range[1] + params['coherence_frequency_upper_bound_correction']

    indices = longest_subsequence_of_trues((coherence > params['cs_coherence_threshold']) &
                                           (frequency > lower_frequency_bound) & (frequency < upper_frequency_bound))
    if indices:
        start, stop = indices
        if (stop - start) < params['coherence_minimum_num_points']:
            raise LowCoherenceError("Not enough points above coherence threshold.")
        else:
            x = frequency[start:stop]
            y = phase[start:stop] * np.pi / 180.0
            slope = np.mean(x * y) / np.mean(x * x)
            return slope, indices
    else:
        raise LowCoherenceError("No points above coherence threshold.")
