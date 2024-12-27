Plotting
=========



Plotting a radargram
--------------------

Function ``plot_radargram()``:

>>> from lib.plot_functions import plot_radargram

.. code-block:: python
    plot_radargram( crd_object, 
                    ax=None,
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
                    show_cbar=True,
                    cmap='binary',
                    vmin='',
                    vmax='')

``crd_object``: bla bla
``ax``        : next
``range_mode``: twt or elevation
   

Plotting a map with the radar profile
-------------------------------------

