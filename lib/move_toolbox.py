
######################################
# Write GOCAD (.pl) PolyLines

def layers2GoCAD(Line="", Layer_dict="", Lon="", Lat="", EPSG_out="", dir_out=""):

    """
    
    
    
    """

    import numpy as np
    import pandas as pd
    from pyproj import Transformer


    layer_list = list(Layer_dict.keys())

    for lr in layer_list:
        if "_m" not in lr or "Surface" in lr:
            pass
        else:
            print("... Layer2GoCAD PolyLine for: {}".format(lr))
            #try:
            trace = Layer_dict[lr]['trace']
            value = Layer_dict[lr]['value']
            color = Layer_dict[lr]['color']

            # Create sets where gaps are to avoid connecting lines
            nan_idxs = np.argwhere(np.isnan(value))
            trace    = np.delete(trace, nan_idxs, axis=0)
            value    = np.delete(value, nan_idxs, axis=0)

            mask       = np.diff(trace)
            split_vals = np.where(mask !=1)[0]

            trace_list = []
            value_list = []

            if len(split_vals) == 0:
                trace_list.append(trace)
                value_list.append(value)

            else:
                trace_list = []
                value_list = []

                for i in np.arange(len(split_vals) +1):
                    if i==0:
                        traces = trace[0:split_vals[i]+1]
                        values = value[0:split_vals[i]+1]
                    elif i>0 and i<len(split_vals):
                        traces = trace[split_vals[i-1]+1:split_vals[i]+1]
                        values = value[split_vals[i-1]+1:split_vals[i]+1]
                    elif i==len(split_vals):
                        traces = trace[split_vals[i-1]+1::]
                        values = value[split_vals[i-1]+1::]

                    trace_list.append(traces)
                    value_list.append(values)

            # Write Sets
            filename   = '{}\\{}_{}.pl'.format(dir_out, lr, Line)
            pline_name = '{}_{}'.format(lr, Line)
            f          = open(filename, 'w')

            for set in np.arange(len(split_vals)+1):

                idx_list  = trace_list[set] - 1
                #idx_list  = idx_list[5:-5]
                Longitude = np.array([Lon[i] for i in idx_list])
                Latitude  = np.array([Lat[i] for i in idx_list])
                value     = value_list[set]
                ids       = np.arange(len(idx_list)) + 1      

                transformer = Transformer.from_crs(4326, EPSG_out, always_xy=True)
                x, y        = transformer.transform(Longitude, Latitude)

                c1, c2, c3  = color[0], color[1], color[2]

                f.write('GOCAD PLine 1 \n')
                f.write('HEADER { \n')
                f.write('name:{}_{} \n'.format(pline_name, set))
                f.write('*visible:false \n')
                f.write('*line*color:{} {} {} 1 \n'.format(c1, c2, c3))
                f.write('} \n')
                f.write('GOCAD_ORIGINAL_COORDINATE_SYSTEM \n')
                f.write('NAME Default \n')
                f.write('AXIS_NAME "X" "Y" "Z" \n')
                f.write('AXIS_UNIT "m" "m" "m" \n')
                f.write('ZPOSITIVE Elevation \n')
                f.write('END_ORIGINAL_COORDINATE_SYSTEM \n')
                f.write('ILINE \n')

                val_mean    = np.array(pd.Series(value).rolling(window=10, center=True).mean())[5:-5]

                for i in range(int(len(val_mean)))[::5]:
                    try:
                        f.write('VRTX {}  {:.4f} {:.4f}   {:.4f} \n'.format(ids[i], x[i], y[i], value[i]))
                    except:
                        pass

                f.write('END \n')

            f.close()
            print('Written: {}'.format(filename))
            
            # except:
            #     print('Error with {}'.format(lr))