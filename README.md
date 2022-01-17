# PyRES
Python code for repeating earthquakes

# Abstract
The software package PyRES is an open source seismological software. It consists of some utilities for data
preparation (utils directory) and the main program (xxxx.py). PyRES is designed to discriminate repeating earthquakes
starting from a family of candidate repeating earthquakes, based on the cross-correlation values and S-P time difference
between pairs of earthquakes (estimated using cross spectrum).

# Motivation and significance

The code PyRES is inspired to previously published methods that combine both seismic waveform similarity, using
cross-correlation function, and differential S-P travel time measured at each seismic station (Chen et al., 2008 and
Shakibay Senobari and Funning, 2019). The code is versatile and works with and without P and S-wave phase pickings. At
the moment the reading of phases in archive Hypoellipse and Hypoinverse formats are implemented. The code has been
tested using synthetic and real data, providing accurate results. It contributes to the implementations of open-source
Python packages in seismology aiming to support the activities of researchers and the reproducibility of scientific
results.

# Requirements

Python installed to run the program (version 3.0 or more) is required [http://python.org]. Dependancies include Obspy,
Numpy, Scipy, Scikit-learn, mtspec (©2009-2016, Lion Krischer, Moritz Beyreuther) and Matplotlib libraries. to do

# How to use it

The package has few dependencies; the recommended way of installing them is via the Conda package manager. You can
create a test environment using

```shell
conda create -n pyres-test python=3.8 numpy=1.21 scipy scikit-learn pandas tqdm obspy mtspec -c conda-forge
```

The package is registered on PyPi, you can install it using `pip`

```shell
pip install py-resw
```

After that, you'll have the `pyres` script in your path.

```shell
# $ pyres --help

usage: ./bin/pyres [-h] [--phase_file PHASE_FILE] [--phase_type {hypoinv,nll,quakeml,hypoel}] [--taup_model TAUP_MODEL] [--rebuild_model]
                   [--graphics_dir GRAPHICS_DIR] [--graphics_format GRAPHICS_FORMAT] [--stop] [--log LOG] [--progress]
                   catalogue inventory parameters output_dir

positional arguments:
  catalogue
  inventory
  parameters
  output_dir

optional arguments:
  -h, --help            show this help message and exit
  --phase_file PHASE_FILE
  --phase_type {hypoinv,nll,quakeml,hypoel}
  --taup_model TAUP_MODEL
  --rebuild_model
  --graphics_dir GRAPHICS_DIR
  --graphics_format GRAPHICS_FORMAT
  --stop
  --log LOG             Log level
  --progress            Show progress bar

```

You can run a test using the [data provided in this repository](data/california).

```shell
pyres cre.zmap inventory.xml parameters.json outputs --phase_file=phases_hypoinv.txt --phase_type=hypoinv --taup_model=ncmodel --progress
```

The numerical parameters are set using the `parameters.json` file. The name of the fields are self-explicative and more
extensive information can be found in the documentation.

# References

Sugan, M., Campanella, S., Vuan, A., Shakibay Senobari, N., (2022). A Python Code for Detecting Repeating Earthquakes
from Self-similar Waveforms (PyRES). Submitted to BSSA ?

Shakibay Senobari, N., Funning G. J., (2019), Widespread Fault Creep in the Northern San Francisco Bay Area Revealed by
Multistation Cluster Detection of Repeating Earthquakes, Geophysical Research Letters, 46(12),
6425-6434, https://doi.org/10.1029/2019GL082766.

Chen, K. H., Nadeau, R. M. and Rau, R.-J. (2008). Characteristic repeating microearthquakes on an arc-continent
collision boundary – the Chihshang fault of eastern Taiwan, Earth Planet. Sci. Lett., 276,
262–272, https://doi.org/10.1016/j.epsl.2008.09.021.

