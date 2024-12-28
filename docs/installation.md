---
title: Installation
---
# Installation

Clone the repository and add the path of the Cradar folder to your PYTHONPATH.

```{code} bash

>>> git clone https://github.com/stvnfrnk/Cradar.git
```

On Windows, go to *Environment Variables* and add the path as *Variable value* and PYTHONPATH as *Variable name*.

Mac or Linux add this to your ~/.bashrc

```{code} bash

>>> export PYTHONPATH="${PYTHONPATH}:/my/other/path"
```

:::{note}
It is probably best to create a new environment and install the following packages:
- numpy, pandas, matplotlib, 
- h5py, scipy, xarray,
- pyproj, geopy, geopandas, shapely, 
- obspy
:::

