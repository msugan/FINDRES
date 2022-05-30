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

```console
foo@bar ~ % conda create -n findres-test python=3.8 numpy=1.21 pandas tqdm pyyaml obspy -c conda-forge
```

Remember to activate the environment with

```console
foo@bar ~ % conda activate findres-test
```

The package multitaper is required, you can install it using `pip`

```console
foo@bar ~ % pip install multitaper
```

The package FINDRES is registered on PyPi, you can install it using `pip`

```console
(findres-test) foo@bar ~ % pip install findres
```

Otherwise, you can download or clone this repository in a specific path, then add `FINDRES/bin` to `PATH` and `FINDRES` to `PYTHONPATH`. This is useful if you want to modify the code, work with your own version of the code, get the latest one, or to extend [custom_formats.py](findres/custom_formats.py) to support your own formats. For example, if the repository is located at `/home/foo/dev/FINDRES` you can append the new paths to the relevant environment variables with

```console
(findres-test) foo@bar ~ % export PATH=$PATH:/home/foo/dev/FINDRES/bin
(findres-test) foo@bar ~ % export PYTHONPATH=$PYTHONPATH:/home/foo/dev/FINDRES
```

After that, you'll have the `findres` script in your path (you can check using `which` that you are calling the script located at the right path).

```console
(findres-test) foo@bar ~ % findres --help
usage: ./bin/findres [-h] [--phase_file PHASE_FILE] [--phase_type {hypoinv,nll,quakeml,hypoel,rise_custom,hyposynth}] [--taup_model TAUP_MODEL] [--rebuild_model] [--graphics_dir GRAPHICS_DIR]
                     [--graphics_format GRAPHICS_FORMAT] [--hypodd] [--include INCLUDE] [--log LOG] [--progress]
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
  --phase_type {hypoinv,nll,quakeml,hypoel,rise_custom,hyposynth}
                        Type of PHASE_FILE (default: None)
  --picker              Enable picker (default: False) 
  --taup_model TAUP_MODEL
                        Velocity model file without extension (assumed to be .tvel) if available (default: None)
  --rebuild_model       Force the rebuild of the velocity model (regenerate ObsPy local cache) (default: False)
  --graphics_dir GRAPHICS_DIR
                        Where to put the graphics relative to OUTPUT_DIR, if not specified the graphics is not generated (default: None)
  --graphics_format GRAPHICS_FORMAT
                        Graphics format, must be one of the extensions recognized by matplotlib (default: pdf)
  --hypodd              Whether to output hypodd input files (default: False)
  --include INCLUDE     Include only the listed events from catalogue (default: None)
  --log LOG             Log level (default: warning)
  --progress            Show progress bar (default: False)
```

You can run a test using the [data provided in this repository](data/california) and accessing the following folder 

```console
(findres-test) foo@bar ~ % cd /home/foo/dev/FINDRES/data/california
```

For example the command to analyse the
data and produce graphics while showing a progress bar is

```console
(findres-test) foo@bar california % findres cre.zmap inventory.xml parameters.yaml results_reference --phase_file=phases_hypoinv.txt --phase_type=hypoinv --taup_model=ncmodel --graphics_dir=figures --progress
```

The command to analyse the
data without graphics (speeding up the time computation) while showing a progress bar is

```console
(findres-test) foo@bar california % findres cre.zmap inventory.xml parameters.yaml results_reference --phase_file=phases_hypoinv.txt --phase_type=hypoinv --taup_model=ncmodel --progress
```

The command to analyse the
data without graphics (speeding up the time computation) and producing the HypodDD output file while showing a progress bar is

```console
(findres-test) foo@bar california % findres cre.zmap inventory.xml parameters.yaml results_reference --phase_file=phases_hypoinv.txt --phase_type=hypoinv --taup_model=ncmodel --progress --hypodd
```

The Hypodd output file can be used to locate the events. HypoDD software (Waldhauser, F., and W.L. Ellsworth 2000) must be installed.
To run an exmaple you can accesss the following directory:

```console
(findres-test) foo@bar california % cd Hypodd_svd_1_RES
```
and type:

```console
(findres-test) foo@bar Hypodd_svd_1_RES % HypoDD hypoDD.inp
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

