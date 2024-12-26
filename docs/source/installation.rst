Installation
============


Required libraries
------------------

For radar data loading and processing:

- numpy, pandas, scipy, h5py, obspy, xarray, pyproj, geopy


For geo-file I/O:

- geopandas, fiona, shapely, pyproj, geopy, rioxarray


For plotting:

- matplotlib, geopandas, shapely, pyproj, rasterio,


Installation
------------

To use Cradar, clone the repository from GitHub:

.. code-block:: console

   $ git clone https://github.com/stvnfrnk/Cradar.git

Add the path to the Cradar folder to your $PYTHONPATH



Module import and Cradar class initialization
---------------------------------------------

>>> from Cradar import Cradar
>>> crd = Cradar()

imports the Cradar module and initiates an empty Cradar class in the "crd" object



