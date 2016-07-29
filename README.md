# Hydroflow
Hydrological network model with linear reservoir elements, originally developed by Glen Liston (Liston and Mernild, 2012). This updated Fortran 2003 code by Arno C Hammann. Licensing information will be added. 

# Prerequisites
gfortran, netcdf fortran library

# Compilation
```
git clone https://github.com/betaplane/hydroflow.git
cd hydroflow
make
```
The executable file will be called ```driver```. A bash run script is provided, ```run.sh```, see comments therein.

# References
Liston, Glen E., and Sebastian H. Mernild. 2012. “Greenland Freshwater Runoff. Part I: A Runoff Routing Model for Glaciated and Nonglaciated Landscapes (HydroFlow).” Journal of Climate 25 (17): 5997–6014. doi:10.1175/JCLI-D-11-00591.1.


