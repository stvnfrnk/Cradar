---
title: Usage
---
# Usage

## Getting started

Import the Cradar class from the Cradar folder and initialize the class.
```{code} python
from Cradar import Cradar

crd = Cradar()
```


## Loading data

Loading AWI nc radar data.
```{code} python
from Cradar import Cradar

awi_nc_file = "Data_20181013_01_016_standard.nc"
profile_id  = "20181013_01_016"

crd = Cradar().load_awi_nc(awi_nc_file, frame=profile_id, read_agc=False)
```
\
Loading CReSIS .mat radar data.
```{code} python
from Cradar import Cradar

cresis_mat_file = "Data_20181013_01_016.mat"

crd = Cradar().load_cresis_mat(cresis_mat_file, dB=False)
```

\
Loading SEGY radar data. Lontitude and Latitude must be provided externally, e.g., from a coordinate file.
```{code} python
from Cradar import Cradar
import pandas as pd

sgy_file   = "example.sgy"
coord_file = "coord_file.csv"

df       = pd.read_csv(coord_file, sep="\t")
Lon, Lat = df["Longitude"].values, df["Latitude"].values

crd = Cradar().load_segy(sgy_file, Longitude=Lon, Latitude=Lat, dB=False)
```


## Cradar attributes

**Key attributes**
| Attribute     | Type           | Definition                            |
|---------------|----------------|---------------------------------------|
| crd.Data      | 2d numpy array | Radar data matrix                     |
| crd.Time      | 1d numpy array | Vertical array of Two-way travel time |
| crd.Longitude | 1d numpy array | Logitude array of radar trace         |
| crd.Latitude  | 1d numpy array | Latitude array of radar trace         |
| crd.Frame     | string         | Profile ID                            |

\
**The Layer dictionary**
| Attribute                             | Type           | Definition                                             |
|---------------------------------------|----------------|--------------------------------------------------------|
| crd.Layer                             | Dictionary     | Dictionary with layer information                      |
| crd.Layer["Surface"]                  |                | Example of "Surface" Layer                             |
| crd.Longitude["Surface"]["trace"]     | 1d numpy array | Array containing trace numbers                         |
| crd.Longitude["Surface"]["value"]     | 1d numpy array | Array containing the twt of "Surface"                  |
| crd.Longitude["Surface"]["value_idx"] | 1d numpy array | Array containing the index (range sample) of "Surface" |

\
**Additional attributes**
| Attribute     | Type    | Definition                                     |
|---------------|---------|------------------------------------------------|
| crd.dB        | boolean | Is data in dB?                                 |
| crd.Reader    | string  | Reading library, e.g., scipy, h5py, obspy, ... |
| crd.Domain    | string  | TWT or elevation/depth                         |
| crd.Spacing   | float   | Trace spacing (add_distance() function)        |
| crd.Distance  | float   | Along-track distance (add_distance() function) |

\