# Networks
We have several networks for UCL_PDR and (currently) no tool to produce Jacobian or ODE files from lists of species and reactions in UCL_PDR's required formats. 

It's therefore most likely that you will want to use one of the networks in this folder. We describe them below. To use a network in this folder, copy the csvs in the Datafiles/Chemical-Network and the code files into Source/

### benchmark

Benchmark is the network used in the Rollig et al. 2007 PDR model benchmarking paper.

### full-network

Full network is our most complete gas phase network including XDR reactions. It does a reasonable job of reproducing the Meijerink et al. 2005 XDR models.