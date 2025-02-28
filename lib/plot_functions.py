
###########################
# Method: plot_radargram()
###########################

def plot_radargram( crd_object, 
                    ax=None,
                    range_mode='',
                    every_dist=["km", 2],
                    every_m_elev=1000,
                    every_twt=['ms', 10],
                    plot_layers=False,
                    markersize=0.2,
                    show_legend=True,
                    xlabels_as_int=True,
                    ylabels_as_int=True,
                    vline_list=[],
                    title="",
                    fontsize=12,
                    show_cbar=True,
                    cmap='binary',
                    vmin='',
                    vmax='',
                    save_svg=False, 
                    save_raster=[False, "jpg"], 
                    suffix='',
                    out_folder='',
                    dpi=200):

    '''
    
        range_mode     : 'twt' or 'elevation'
        xlabels_as_int : 
        ylabels_as_int : 
        every_km_dist  :
        every_m_elev   : 
        every_twt_ms   : 


    '''

    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from pyproj import Transformer
    # import rasterio
    # import geopandas as gpd
    import pandas as pd
    import numpy as np 
    import os

    if ax is None:
        ax = plt.gca()

    if every_dist[0] == "m":
        # Build xticks every n meters
        distance   = crd_object.Distance
        num_xticks = int(distance[-1]/every_dist[1])
        tmp        = np.repeat(every_dist[1], num_xticks)
        tmp[0]     = 0
        xticks     = np.cumsum(tmp)
        xticks_idx = []

        for i in range(len(xticks)):
            idx = np.abs(distance - xticks[i]).argmin()
            xticks_idx.append(idx)
            
        xticks_idx    = np.array(xticks_idx)
        xtick_labels  = xticks.astype(int)
            
        xaxis_label   = 'Distance (m)'

    if every_dist[0] == "km":
        # Build xticks every n meters
        distance   = crd_object.Distance / 1000
        num_xticks = int(distance[-1]/every_dist[1])
        xticks     = np.linspace(0, num_xticks*every_dist[1], num_xticks + 1)
        xticks_idx = []

        for i in range(len(xticks)):
            idx = np.abs(distance - xticks[i]).argmin()
            xticks_idx.append(idx)
            
        xticks_idx    = np.array(xticks_idx)
        xtick_labels  = xticks.astype(int)
            
        xaxis_label   = 'Distance (km)'


    # Build yticks every n µs
    if range_mode == 'twt':

        if every_twt[0] == 'ns':
            twt        = crd_object.Time * 10e8
        if every_twt[0] == 'ms':
            twt        = crd_object.Time * 10e5


        num_yticks = int(twt[-1]/every_twt[1])

        yticks_ms     = np.linspace(0, num_yticks*every_twt[1], num_yticks + 1)
        yticks_ms_idx = []

        for i in range(len(yticks_ms)):
            idx = np.abs(twt - yticks_ms[i]).argmin()
            yticks_ms_idx.append(idx)

        yticks_ms_idx = np.array(yticks_ms_idx)

        yticks        = yticks_ms_idx
        ytick_labels  = yticks_ms
        
        if ylabels_as_int == True:
            ytick_labels  = yticks_ms.astype(int)
            
        if every_twt[0] == 'ns':
            yaxis_label   = 'TWT (ns)'
        if every_twt[0] == 'ms':
            yaxis_label   = 'TWT (µs)'
        

    # Build yticks every n meters
    if range_mode == 'elevation':
        
        elevation_m = crd_object.Z
        below_zero  = int(elevation_m.min() / every_m_elev)
        above_zero  = int(elevation_m.max() / every_m_elev)

        below_zero_steps = np.linspace(below_zero*every_m_elev, 0, np.abs(below_zero) + 1)
        above_zero_steps = np.linspace(every_m_elev, above_zero*every_m_elev, np.abs(above_zero))

        yticks_m     = np.concatenate((below_zero_steps, above_zero_steps))
        yticks_m_idx = []

        for i in range(len(yticks_m)):
            idx = np.abs(elevation_m - yticks_m[i]).argmin()
            yticks_m_idx.append(idx)

        yticks_m_idx = np.array(yticks_m_idx)

        yticks        = yticks_m_idx
        ytick_labels  = yticks_m
        
        if ylabels_as_int == True:
            ytick_labels  = yticks_m.astype(int)
            
        yaxis_label   = 'Elevation (m)'


    # Build yticks every n meters
    if range_mode == 'depth':
        
        depth_m     = crd_object.Depth
        num_yticks = int(depth_m[-1]/every_m_elev)

        yticks_m     = np.linspace(0, num_yticks*every_m_elev, num_yticks + 1)
        yticks_m_idx = []

        for i in range(len(yticks_m)):
            idx = np.abs(depth_m - yticks_m[i]).argmin()
            yticks_m_idx.append(idx)

        yticks_m_idx = np.array(yticks_m_idx)

        yticks        = yticks_m_idx
        ytick_labels  = yticks_m
        
        if ylabels_as_int == True:
            ytick_labels  = yticks_m.astype(int)
            
        yaxis_label   = 'Depth (m)'
        
    # set vmin/vmax depending if values are given
    if vmin == '':
        vmin = np.nanmin(crd_object.Data)
    else:
        vmin = vmin

    if vmax == '':
        vmax = np.nanmax(crd_object.Data)
    else:
        vmax = vmax    
        
    # imshow radargram    
    radargram = ax.imshow(crd_object.Data, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.8)

    # plot layers
    if plot_layers == True:
        layer_list = list(crd_object.Layer.keys())
        n          = len(layer_list) #- 2
        offset     = 0
        colors     = plt.cm.spring(np.linspace(0,1,n + offset))
        c = 0

        if range_mode == 'twt':
            for lr in layer_list:
                crd_object.reshape_layer_values(lr)

            for lr in layer_list:    
                if lr == 'Surface':
                    ax.scatter(x=crd_object.Layer[lr]['trace'], y=crd_object.Layer[lr]['value_idx'], 
                        s=markersize, label=lr)

                elif lr == 'Base':
                    ax.scatter(x=crd_object.Layer[lr]['trace'], y=crd_object.Layer[lr]['value_idx'], 
                        color=crd_object.Layer[lr]['color'], s=markersize, label=lr)

                elif lr != 'Surface' and lr != 'Base' and lr != 'Surface_m':
                    ax.scatter(x=crd_object.Layer[lr]['trace'], y=crd_object.Layer[lr]['value_idx'], 
                            color=crd_object.Layer[lr]['color'], s=markersize, label=lr)

        elif range_mode == 'elevation':
            for lr in layer_list:
                if '_m' not in lr:
                    pass
                else:
                    if lr == 'Surface_m':
                        ax.scatter(x=crd_object.Layer[lr]['trace'], y=crd_object.Layer[lr]['value_idx'], 
                            color=crd_object.Layer[lr]['color'], s=markersize, label=lr)

                    elif lr == 'Bed_m':
                        ax.scatter(x=crd_object.Layer[lr]['trace'], y=crd_object.Layer[lr]['value_idx'], 
                            color=crd_object.Layer[lr]['color'], s=markersize, linestyle='dashed', label=lr)

                    elif lr != 'Surface_m':
                        if lr != 'Bed_m':
                            if "_m" in lr:
                                ax.scatter(x=crd_object.Layer[lr]['trace'], y=crd_object.Layer[lr]['value_idx'], 
                                        color=crd_object.Layer[lr]['color'], s=markersize, label=lr)

    if vline_list != []:
        # print('vlines')
        for vline in vline_list:
            ax.vlines(x=vline, ymin=0, ymax=len(crd_object.Time), color='white')

    ax.set_xticks(xticks_idx, xtick_labels, fontsize=fontsize)
    ax.set_xlabel(xaxis_label, fontsize=fontsize)
    ax.set_yticks(yticks, ytick_labels, fontsize=fontsize)
    ax.set_ylabel(yaxis_label, fontsize=fontsize)
    ax.set_xlim(0, len(crd_object.Longitude))
    ax.set_ylim(crd_object.Data.shape[0], 0)
    ax.set_title(title, fontsize=fontsize)



    if show_cbar == True:
        cbr = plt.colorbar(radargram, ax=ax)
        cbr.set_label('dB')

    if show_legend == True:
        ax.legend()

    if save_raster[0] == True:
        if save_raster[1] == "jpg":
            raster_type = ".jpg"
        if save_raster[1] == "png":
            raster_type = ".png" 

        if out_folder == '':
            if not os.path.exists('figures'):
                os.makedirs('figures')

            figname = str(crd_object.Frame) + suffix + raster_type
            plt.savefig('figures/' + figname, dpi=dpi, bbox_inches='tight')
            print('==> Written: figures/{}'.format(figname))
        
        else:
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)

            figname = str(crd_object.Frame) + suffix + raster_type
            plt.savefig(out_folder + '/' + figname, dpi=dpi, bbox_inches='tight')
            print('==> Written: {}/{}'.format(out_folder, figname))


    if save_svg == True:
        if out_folder == '':
            if not os.path.exists('figures'):
                os.makedirs('figures')

            figname = str(crd_object.Frame) + suffix + '.svg'
            plt.savefig('figures/' + figname)
            print('==> Written: figures/{}'.format(figname))
        
        else:
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)

            figname = str(crd_object.Frame) + suffix + '.svg'
            plt.savefig(out_folder + '/' + figname)
            print('==> Written: {}/{}'.format(out_folder, figname))

    return radargram



#####################
# Method: plot_map()
#####################

def plot_map(crd_object, 
             ax=None,
             flight_lines="",
             geotif="",
             points="",
             lines="",
             polygons="",
             epsg=""):
                    
    '''


    '''

    import matplotlib.pyplot as plt
    from pyproj import Transformer
    import rasterio
    from rasterio.plot import show
    import geopandas as gpd
    import pandas as pd
    import numpy as np


    point_dict_list   = points
    line_dict_list    = lines
    polygon_dict_list = polygons
    

    if ax is None:
        ax = plt.gca()

    flight_lines = flight_lines
    Lon          = crd_object.Longitude
    Lat          = crd_object.Latitude

    if flight_lines != "":
        survey_lines = gpd.read_file(flight_lines)
        survey_lines = survey_lines.to_crs(epsg)

    # generate data frames for points
    df          = pd.DataFrame(Lon)
    df['Lat']   = pd.DataFrame(Lat)
    df.columns  = ['Lon', 'Lat']

    # reproject coordinates to target epsg
    transformer      = Transformer.from_crs(4326, epsg, always_xy=True)
    Lon, Lat         = crd_object.Longitude, crd_object.Latitude
    df["X"], df["Y"] = transformer.transform(Lon, Lat)

    df_first    = df[0:1]

    # create geopandas data frames
    frame  = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["X"], df["Y"]))
    first  = gpd.GeoDataFrame(df_first, geometry=gpd.points_from_xy(df_first["X"], df_first["Y"]))

    # set crs to EPSG:4326
    frame        = frame.set_crs(epsg=epsg)
    first        = first.set_crs(epsg=epsg)

    # plot geotif if available
    if geotif != "":
        src = rasterio.open(geotif)
        # extent = [src.bounds[0], src.bounds[2], src.bounds[1], src.bounds[3]]
        # rasterio.plot.show(src, extent=extent, ax=ax)
        show(src, ax=ax)

    # plot flight lines if available
    if flight_lines != "":
        survey_lines.plot(ax=ax, color='black', linewidth=0.5, zorder=1, alpha=1)

    # plot profile coordinates
    frame.plot(ax=ax, color='blue', markersize=15, zorder=2)
    first.plot(ax=ax, color='red', markersize=35, zorder=3)
    ax.set_xlabel("X (EPSG:{})".format(epsg))
    ax.set_ylabel("Y (EPSG:{})".format(epsg))

    # plot polygons if available
    if polygons == "":
        pass
    else:
        for polygon_dict in polygon_dict_list:
            pol_file  = polygon_dict["file"]
            edgecolor = polygon_dict["edgecolor"]
            linewidth = polygon_dict["linewidth"]
            alpha     = polygon_dict["alpha"]
            zorder    = polygon_dict["zorder"]
            try:
                facecolor = polygon_dict["facecolor"]
            except:
                pass

            gpd_polygon = gpd.read_file(pol_file)

            if polygon_dict["boundary"] == True:
                gpd_polygon.boundary.plot(ax=ax, edgecolor=edgecolor, linewidth=linewidth, alpha=alpha, zorder=zorder)
            else:
                gpd_polygon.plot(ax=ax, edgecolor=edgecolor, facecolor=facecolor, linewidth=linewidth, alpha=alpha, zorder=zorder)


    # plot lines if available
    if lines == "":
        pass
    else:
        for line_dict in line_dict_list:
            line_file  = line_dict["file"]
            color      = line_dict["color"]
            linewidth  = line_dict["linewidth"]
            alpha      = line_dict["alpha"]
            zorder     = line_dict["zorder"]

            gpd_line = gpd.read_file(line_file)
            gpd_line.plot(ax=ax, color=color, linewidth=linewidth, alpha=alpha, zorder=zorder)


    # plot points if available
    if points == "":
        pass
    else:
        for point_dict in point_dict_list:
            point_file = point_dict["file"]
            facecolor  = point_dict["facecolor"]
            edgecolor  = point_dict["edgecolor"]
            markersize = point_dict["markersize"]
            alpha      = point_dict["alpha"]
            zorder     = point_dict["zorder"]

            gpd_point = gpd.read_file(point_file)
            gpd_point.plot(ax=ax, facecolor=facecolor, edgecolor=edgecolor, markersize=markersize, alpha=alpha, zorder=zorder)




    return map

