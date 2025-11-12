# UCL_PDR

UCL_PDR is a PDR code which solves the equilibrium chemistry and temperature of an arbitrary 1D cloud. It was written by [Bell et al. 2006](https://ui.adsabs.harvard.edu/abs/2006A%26A...459..805B/abstract) and modified by [Priestley et al. 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.4444P/abstract). The code is no longer maintained by Viti group, the code remains available online as a reference.

# Usage
First, get ahold of a copy of sundials version 2.5.0. A copy is provided within this reposistory. 
We thank the Sundials developers and acknowledge their license and terms. If you use this repository as a mirror,
please do so as well. You can install sundials 2.5.0 using a few easy steps:
#### Installing sundials 2.5.0
Go into the sundials directory and run the following commands with the path where you want to install
sundials
```
./configure --prefix $PWD/../sundials
make
make install
```
### Installing UCL_PDR

First go into the source directory, this allows you to install run the make command:
```
cd Source/
make
```
If you get an error with regards to sundials not being found, please edit the INCLUDES and
LIBRARIES to absolute paths on your system.
It can then be executed via (assure your path is no longer than 256 characters.)

```
./UCL-PDR INPUT_FILE_PATH
```
UCL_PDR will then use the contents of INPUT_FILE_PATH to determine the parameters of the model to be run. 

For example try:

```
./UCL-PDR Input/10_1e3.dat
```

## Inputs

A key input is the cloud file. The name of this file should be supplied in `Input/model-parameters.dat` and many cloud files exist in `Clouds/`. These are csv files containing the x co-ordinate and density of a series of points at which UCL_PDR will evaluate the equilibrium state. `Scripts/constant_density_cloud.py` and `Scripts/make_cloud.py` provide convenient ways to create new cloud files.

Chemical networks can be changed by copying the relevant files. All networks are stored in `Chemical-Networks/`, they contain .dat files which should be moved to `Datafiles/Chemical-Network/` and .c files which should be moved to `Source/` this mostly easily done by invoking

```
./Scripts/switch-network.sh full-network
```
where `full-network` is the name of the desired network in this example.


## Outputs
UCL_PDR produces a huge amount of output however, this can all be compressed into an easily managed .csv file using

```
python Scripts/make_results_df.py my_prefix
```
where `my_prefix` is the file prefix set in `Input/model-parameters.dat`. This csv file can also be plotted using `Scripts/plot_results.py`.
