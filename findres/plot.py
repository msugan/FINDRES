from matplotlib import figure
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
from obspy.geodetics import gps2dist_azimuth

from .read import coordinates


def signals(trace1, trace2, p_wave_window_1, p_wave_window_2, s_wave_window_1, s_wave_window_2, text, filename):
    fig = figure.Figure()
    ax = fig.add_subplot()
    ax.plot(np.arange(trace1.stats.npts) * trace1.stats.delta, trace1.data / abs(trace1.max()), linewidth=1,
            color="blue", alpha=0.5)
    ax.plot(np.arange(trace2.stats.npts) * trace2.stats.delta, trace2.data / abs(trace2.max()), linewidth=1,
            color="red", alpha=0.5)
    ax.set_xlabel("Times (sec)")
    ax.set_ylabel("Normalized Amplitude")

    ax.text(0.5, 1, text)
    ax.set_ylim(-1.3, 1.3)
    ax.hlines(-1, p_wave_window_1, p_wave_window_2, color="r", linewidth=1, linestyles="dashed")
    ax.hlines(-1, s_wave_window_1, s_wave_window_2, color="r", linewidth=1, linestyles="dashed")
    fig.savefig(filename)


def cross_spectrum(trace1, trace2, cs, slope, indices, textline, filename):
    fig = figure.Figure(constrained_layout=True)
    waveforms, analysis = fig.subfigures(2, 1)

    waveform1, waveform2 = waveforms.subplots(2, 1, sharex=True)
    data1 = trace1.data
    data1 = 2 * (data1 - data1.min()) / (data1.max() - data1.min()) - 1
    waveform1.plot(np.arange(trace1.stats.npts) * trace1.stats.delta, data1, linewidth=1)
    waveform1.set_ylabel("Amplitude")
    data2 = trace2.data
    data2 = 2 * (data2 - data2.min()) / (data2.max() - data2.min()) - 1
    waveform2.plot(np.arange(trace2.stats.npts) * trace2.stats.delta, data2, linewidth=1)
    waveform2.set_ylabel("Amplitude")
    waveform2.set_xlabel("Time [sec]")

    coherence, phase = analysis.subplots(2, 1, sharex=True)
    positive_freq_indices = cs.freq.flatten() > 0.0
    x = cs.freq.flatten()
    y = cs.cohe.flatten()
    z = cs.phase.flatten() * np.pi / 180.0
    start, stop = indices
    coherence.plot(x[positive_freq_indices], y[positive_freq_indices], linewidth=1)
    coherence.plot(x[start:stop], y[start:stop], "rx", label="used data", markersize=4)
    coherence.set_ylabel("Coherency")
    coherence.set_ylim(0, 1.1)
    phase.plot(x[positive_freq_indices], z[positive_freq_indices], linewidth=1)
    phase.plot(x[start:stop], z[start:stop], "rx", label="used data", markersize=4)
    phase.plot(x[start:stop], slope * x[start:stop], "g", label="fitted line")
    phase.set_xlabel("Frequency [Hz]")
    phase.set_ylabel("Phase [rad]")
    phase.set_ylim(-3, 3)
    phase.text(1, 1, textline)
    fig.savefig(filename)


def similarity(ev1, ev2, data, inventory, cc_threshold, delta_threshold, filename):
    event_latitude = (ev1.latitude + ev2.latitude) / 2
    event_longitude = (ev1.longitude + ev2.longitude) / 2
    xs = np.fromiter((x for x, _ in data[(ev1.Index, ev2.Index)].values()), dtype=np.float)
    ys = np.fromiter((y for _, y in data[(ev1.Index, ev2.Index)].values()), dtype=np.float)
    zs = 1e-3 * np.fromiter(map(lambda s: gps2dist_azimuth(event_latitude, event_longitude,
                                                           *coordinates(inventory, s[0], s[1], ev1.date))[0],
                                data[(ev1.Index, ev2.Index)]),
                            dtype=np.float)
    fig = figure.Figure()
    ax = fig.add_subplot()
    ax.scatter(xs, np.absolute(ys), c=zs, cmap="jet")
    ax.set_xlim(cc_threshold, 1.0)
    ax.set_ylim(0.0, delta_threshold)
    ax.set_xlabel("CC Z_component")
    ax.set_ylabel("|ΔS-P|")
    ax.set_title(f"{ev1.name}-{ev2.name}")
    cb = fig.colorbar(ScalarMappable(norm=Normalize(vmin=zs.min(), vmax=zs.max()), cmap="jet"))
    cb.set_label("epi station-event distance (km)")
    ax.text(cc_threshold * 1.01, delta_threshold * 0.95, f"M_{ev1.name}={ev1.magnitude}")
    ax.text(cc_threshold * 1.01, delta_threshold * 0.9, f"M_{ev2.name}={ev2.magnitude}")
    fig.savefig(filename)
