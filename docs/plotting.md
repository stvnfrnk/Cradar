---
title: Plotting
---
# Plotting


## Plotting a radargram



| Argument       | Type          | Definition                                                         | Example               |
|----------------|---------------|--------------------------------------------------------------------|-----------------------|
| crd_object     | Cradar object |                                                                    | crd_object,           |
| ax             |               | axis label (None for single-panel plot)                            | ax=None,              |
| range_mode     | string        | "twt" or "elevation"                                               | range_mode='',        |
| every_km_dist  | int           | x-axis km along-track distance interval                            | every_km_dist=10,     |
| every_m_elev   | int           | y-axis m elevation intreval                                        | every_m_elev=1000,    |
| every_twt      | list          | y-axis tick interval ("ms"=µs, "ns"=nanoseconds)                   | every_twt=['ms', 10], |
| plot_layers    | bool          | show layers in plot                                                | plot_layers=False,    |
| markersize     | float         | markersize of layers                                               | markersize=0.2,       |
| show_legend    | bool          | show layer legend                                                  | show_legend=True,     |
| xlabels_as_int | bool          |                                                                    | xlabels_as_int=True,  |
| ylabels_as_int | bool          |                                                                    | ylabels_as_int=True,  |
| vline_list     | list          | plot additional vlines (takes list of integers with trace numbers) | vline_list=[],        |
| fontsize       | int           | fontsize of title                                                  | fontsize=12,          |
| show_cbar      | bool          | show colorbar                                                      | show_cbar=True,       |
| cmap           | string        | colormap (takes matplotlib names)                                  | cmap='binary',        |
| vmin           | string        | imshow vmin limit (is Data.min() by default)                       | vmin='',              |
| vmax           | string        | imshow vmax limit (is Data.max() by default)                       | vmax='',              |


## Plotting a map



## Multi-panel plots

Multi-panel plots can be embedded for example as plt.subplots where the axes need to be defined.

### Example with radargram and ice thickness below


### Example with profile location on a map and radargram

```{code} python

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from Cradar import Cradar
from lib.plot_functions import plot_radargram, plot_map

test_pol = "shapes/testpol.shp"
measures = "raster/measures_mag_ant_clip_yuting.tif"

crd1       = Cradar().load_awi_nc("Data_20181013_01_016_standard.nc", read_agc=False)
crd1.Frame = "20181013_01_016"
crd1.add_distance()
crd1.clip_range(0, 20000, mode="twt")
crd1.get_layer_idx()

crd2       = Cradar().load_awi_nc("Data_20181013_01_017_standard.nc", read_agc=False)
crd2.Frame = "20181013_01_017"
crd2.add_distance()
crd2.clip_range(0, 20000, mode="twt")
crd2.get_layer_idx()


# make figure with subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
                                    ncols=2,
                                    nrows=2,
                                    figsize=(15, 8))

plot_map( crd1, 
          ax=ax1,
          flight_lines="",
          geotif=measures,
          epsg=3031)

plot_radargram( crd1, 
                ax=ax2,
                range_mode='twt',
                every_km_dist=10,
                every_m_elev=1000,
                every_twt=['ms', 10],
                plot_layers=True,
                markersize=0.2,
                show_legend=True,
                xlabels_as_int=True,
                ylabels_as_int=True,
                vline_list=[200, 400],
                fontsize=12,
                show_figure=True, 
                show_cbar=False,
                cmap='binary',
                vmin='',
                vmax='')

plot_map( crd2, 
          ax=ax3,
          flight_lines="",
          geotif=measures,
          epsg=3031)

plot_radargram( crd2,
                ax=ax4,
                range_mode='twt',
                every_km_dist=10,
                every_m_elev=1000,
                every_twt=['ms', 10],
                plot_layers=True,
                markersize=0.2,
                show_legend=True,
                xlabels_as_int=True,
                ylabels_as_int=True,
                vline_list=[200, 400],
                fontsize=12,
                show_figure=True, 
                show_cbar=False,
                cmap='binary',
                vmin='',
                vmax='')

plt.tight_layout()
plt.subplots_adjust(left=0.2, bottom=None, right=None, top=None, wspace=None, hspace=None)

'''
left  = 0.125  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.2   # the amount of height reserved for white space between subplots
'''
```


```{figure} images/multi_panel_plot_with_map.png
:scale: 90 %
:alt: map to buried treasure

This is the caption of the figure (a simple paragraph).


## Saving plots