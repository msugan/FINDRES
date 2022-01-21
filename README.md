# FINDRES
A Python code for detecting true Repeating Earthquakes from Self-similar Waveforms (FINDRES)

# Abstract
The software package FINDRES is an open-source seismological software. FINDRES is designed 
to discriminate repeating earthquakes starting from a family of candidate repeating earthquakes, 
based on the cross-correlation values and S-P time difference between pairs of earthquakes 
(estimated using cross-spectrum).

# Motivation and significance

The code FINDRES is inspired to previously published methods that combine both seismic waveform 
similarity, using cross-correlation function, and differential S-P travel time measured at each 
seismic station (Chen et al., 2008 and Shakibay Senobari and Funning, 2019). The code is versatile 
and works with and without P and S-wave phase pickings. At the moment the reading of phases in 
Hypoellipse (Lahr 1999), Hypoinverse (Klein, 2002), NonlinLoc (Lomax et al. 2000), 
and QuakeML (https://quake.ethz.ch/quakeml) format are implemented. The code has been tested using 
synthetic and real data, providing accurate results. It contributes to the implementations of open-source
Python packages in seismology aiming to support the activities of researchers and the reproducibility of scientific
results.

# Requirements

Python installed to run the program (version 3.8 or more) is required [http://python.org]. 

# How to use it

The package has few dependencies; the recommended way of installing them is via the Conda package manager. You can
create a test environment using

```shell
conda create -n findres-test python=3.8 numpy=1.21 scipy scikit-learn pandas tqdm pyyaml obspy mtspec -c conda-forge
```

Remember to activate the environment with

```shell
conda activate findres-test
```

The package is registered on PyPi, you can install it using `pip`

```shell
pip install findres
```

After that, you'll have the `findres` script in your path.

```shell
# findres --help

usage: ./bin/findres [-h] [--phase_file PHASE_FILE] [--phase_type {hypoinv,nll,quakeml,hypoel}] [--taup_model TAUP_MODEL] [--rebuild_model]
                     [--graphics_dir GRAPHICS_DIR] [--graphics_format GRAPHICS_FORMAT] [--hypodd] [--stop] [--log LOG] [--progress]
                     catalogue inventory parameters output

positional arguments:
  catalogue             Modified ZMAP catalogue containing repeater candidates, their MSEED location, and names
  inventory             Inventory of the stations data
  parameters            Numerical parameters file
  output                Output YAML file

optional arguments:
  -h, --help            show this help message and exit
  --phase_file PHASE_FILE
                        Catalogue containing the phase picking and other information if available (default: None)
  --phase_type {hypoinv,nll,quakeml,hypoel}
                        Type of PHASE_FILE (default: None)
  --taup_model TAUP_MODEL
                        Velocity model file without extension (assumed to be .tvel) if available (default: None)
  --rebuild_model       Force the rebuild of the velocity model (regenerate ObsPy local cache) (default: False)
  --graphics_dir GRAPHICS_DIR
                        Where to put the graphics relative to OUTPUT_DIR, if not specified the graphics is not generated (default: None)
  --graphics_format GRAPHICS_FORMAT
                        Graphics format, must be one of the extensions recognized by matplotlib (default: pdf)
  --hypodd              Whether to output hypodd input files (default: False)
  --stop                Stop if exceptions are raised during the analysis, otherwise skip to the next event station (default: False)
  --log LOG             Log level (default: info)
  --progress            Show progress bar (default: False)
```

You can run a test using the [data provided in this repository](data/california) and accessing the following folder 

```shell
cd data/california
```

For example the command to analyse the
data and produce graphics while showing a progress bar is

```shell
findres cre.zmap inventory.xml parameters.yaml results --phase_file=phases_hypoinv.txt --phase_type=hypoinv --taup_model=ncmodel --graphics_dir=figures --progress
```

The command to analyse the
data without graphics (speeding up the time computation) while showing a progress bar is

```shell
findres cre.zmap inventory.xml parameters.yaml results --phase_file=phases_hypoinv.txt --phase_type=hypoinv --taup_model=ncmodel --progress
```

The command to analyse the
data without graphics (speeding up the time computation) and producing the HypodDD output file while showing a progress bar is

```shell
findres cre.zmap inventory.xml parameters.yaml results --phase_file=phases_hypoinv.txt --phase_type=hypoinv --taup_model=ncmodel --progress --hypodd
```

The Hypodd output file can be used to locate the events. HypoDD software (Waldhauser, F., and W.L. Ellsworth 2000) must be installed.
To run an exmaple you can accesss the following directory:

```shell
cd /data/california/Hypodd_svd_1_RES
```

and type:

```shell
HypoDD hypoDD.inp
```

The numerical parameters are set using the `parameters.yaml` file. The name of the fields are self-explicative and more
extensive information can be found in the [provided example](data/california/parameters.yaml).

# References

Sugan, M., Campanella, S., Vuan, A., and Shakibay Senobari, N., (2022). A Python Code for detecting true Repeating Earthquakes
from Self-similar Waveforms (FINDRES). Submitted

Shakibay Senobari, N., and Funning, G. J., (2019). Widespread Fault Creep in the Northern San Francisco Bay Area Revealed by
Multistation Cluster Detection of Repeating Earthquakes, Geophysical Research Letters, 46(12),
6425-6434, https://doi.org/10.1029/2019GL082766.

Chen, K. H., Nadeau, R. M, and Rau, R.-J. (2008). Characteristic repeating microearthquakes on an arc-continent
collision boundary – the Chihshang fault of eastern Taiwan, Earth Planet. Sci. Lett., 276,
262–272, https://doi.org/10.1016/j.epsl.2008.09.021.

