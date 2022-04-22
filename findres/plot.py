import matplotlib.pyplot as plt
import numpy as np
from obspy.geodetics import gps2dist_azimuth

from .read import coordinates


def signals(trace1, trace2, mag_t1, mag_t2, p_wave_window_1, p_wave_window_2, s_wave_window_1, s_wave_window_2,
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


def cross_spectrum(trace1, trace2, a_take, b_take, c_take, slope, time_delay, out, filename):
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
    plt.plot(out.freq.flatten(), out.cohe.flatten(), linewidth=1)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Coherency")
    plt.ylim(0, 1.1)
    plt.plot(b_take, a_take, "rx", label="original data", markersize=4)

    plt.subplot(414)
    plt.plot(out.freq.flatten(), out.phase.flatten() * np.pi / 180.0, linewidth=1)
    plt.plot(b_take, c_take, "rx", label="original data", markersize=4)
    plt.plot(b_take.reshape(-1, 1), slope * b_take.reshape(-1, 1), "g", label="fitted line")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Phase(rad)")
    plt.ylim(-3, 3)
    plt.tight_layout()
    plt.text(1, 1, f"time_delay = {time_delay:0.6f}")

    plt.savefig(filename)
    plt.close()


def similarity(ev1, ev2, data, inventory, cc_threshold, delta_threshold, filename, **kwargs):
    event_latitude = (ev1.latitude + ev2.latitude) / 2
    event_longitude = (ev1.longitude + ev2.longitude) / 2
    plt.figure(figsize=(7, 5), **kwargs)
    xs = np.fromiter((x for x, _ in data[(ev1.Index, ev2.Index)].values()), dtype=np.float)
    ys = np.fromiter((y for _, y in data[(ev1.Index, ev2.Index)].values()), dtype=np.float)
    zs = np.fromiter(map(lambda s: gps2dist_azimuth(event_latitude, event_longitude, *coordinates(inventory, s, ev1.date))[0],
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
