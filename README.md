# UCL_PDR

UCL_PDR is a PDR code which solves the equilibrium chemistry and temperature of an arbitrary 1D cloud. It was written by [Bell et al. 2006](https://ui.adsabs.harvard.edu/abs/2006A%26A...459..805B/abstract) and modified by [Priestley et al. 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.4444P/abstract). It is currently maintained by Jon Holdship.

# Usage
UCL_PDR is currently command line only. After compiling the code with

```
cd Source/
make
```
it can be run via

```
./UCL_PDR
```

UCL_PDR will then use the contents of `Input/model-parameters.dat` to determine the parameters of the model to be run. 

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
