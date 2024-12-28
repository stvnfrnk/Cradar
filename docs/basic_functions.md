---
title: Basic functions
---
# Basic functions

## ...






**Convert to and from dB**
```{code} python
# to dB (checks if crd.dB is False)
crd.to_dB()         # sets crd.dB=True

# inverse dB (checks if crd.dB is True)
crd.inverse_dB()    # sets crd.dB=False
```


**Get layer index**\
Adds the range bin index (*value_idx*) to a layer.
```{code} python
crd.get_layer_idx()    
```


**Add spacing and distance**
```{code} python
crd.add_distance()  
```



## Trimming and merging data

**Clip data in along-track (x-axis)**
```{code} python
crd.clip_along(start_val, end_val, mode='trace')       # via trace numbers
crd.clip_along(start_val, end_val, mode='distance')    # via distance (needs crd.add_distance()) 
```

**Clip data in range (y-axis)**
```{code} python
crd.clip_range(start_val, end_val, mode='range_bin')   # via range bin numbers
crd.clip_range(start_val, end_val, mode='twt')         # via twt limits
crd.clip_range(start_val, end_val, mode='elevation')   # via elevation limits
```

**Flip radargram on the y-axis**
```{code} python
crd.flip_lr()
```

**Concatenating frames**
Flip radargram on the y-axis
```{code} python
Cradar.concat_frames(added_objects)

## Example
crd1 = Cradar().load_awi_nc("Data_20181013_01_016_standard.nc", frame="20181013_01_016", read_agc=False)
crd2 = Cradar().load_awi_nc("Data_20181013_01_017_standard.nc", frame="20181013_01_017", read_agc=False)
crd3 = Cradar().load_awi_nc("Data_20181013_01_018_standard.nc", frame="20181013_01_018", read_agc=False)

crd_concat = Cradar.concat_frames([crd1, crd2, crd3])
```


## Loading layer data


**Track surface if *Surface* layer is missing**

track_surface(self, skip=100, llim='', gauss_factor=1, use_gradient=False, offset=2)

**Retrac surface**
retrack_surface(self, roll_factor=10, sigma=2, gate=[0, 0], differenciate=False, offset=1)


## Add data from Geo-Raster



## Transform between domains


## Gain

range_gain(self, gain_type='', b=2, n=2, f=2)


## Correcting for attenuation


## Filtering


## Test