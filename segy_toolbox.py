# Toolbox to deal with SEGY Data


def plot_segy(segy_input, png=False):
    
    '''
    
    
    '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    from obspy.io.segy.segy import _read_segy
    #from obspy import Trace, Stream
    
    stream  = _read_segy(segy_input, headonly=True)
    data    = np.stack(t.data for t in list(stream.traces))
    
    fig = plt.subplots(figsize=(20,10))
    plt.imshow(data.T, cmap="bone_r", aspect='auto')
    plt.colorbar()
    plt.show()

    if png == True:
        plt.savefig('test_segy_plot.png', dpi=200, bbox_inches='tight')
        print('Saved Figure as ==> test_segy_plot.png')
        
    else:
        return fig




def mat2segy(matfile, elevation=True, region='', differenciate=False):

    '''  
    Usage:
            mat2segy(matfile, elevation=True, region='', differenciate=False)


    Creates a SEGY File from a CReSIS Radar .mat file
    
        - Parameters are:
            
            1) matfile 
            2) elevation=True/False
            3) region='Greenland' or 'Antarctica'
            4) differenciate=True/False (False is default)
    
                1) Provide a valid .mat file in TWT or 
                   true elevation (Elevation_WGS84)
            
                2) elevation is set True for Files that are already 
                   in true elevation and not as TWT
                
                3) if the parameter 'region' is set empty, you already
                   should have 'X' and 'Y' coordinates available for your region
                   
                   --> if you set it to 'Greenland'  they will be projected to EPSG:3413
                   --> if you set it to 'Antarctica' they will be projected to EPSG:3031
                   
                4) choose weather to take the data as it is or to differenciate
                   --> set differenciate=True (default is differenciate=False)
           
        
        - if elevation == True:    
                    
                mat             = scipy.io.loadmat(matfile)
                data            = mat['Data']
                X               = mat['X']
                Y               = mat['Y']
                Elevation_max   = mat['Elevation_WGS84'].max()
                num_of_samples  = mat['Elevation_WGS84'].shape[1]  
                gps_time        = mat['GPS_Time']
            
            
        - if elevation == False:
            
                mat             = scipy.io.loadmat(matfile)
            
                # Reproject
                EPSG_=pyproj.Proj(Projection)
                X, Y = EPSG_(mat['Longitude'], mat['Latitude'])
                
                data            = mat['Data']
                X               = X
                Y               = Y
                TWT             = mat['Time']
                num_of_samples  = mat['Time'].shape[0]  
                gps_time        = mat['GPS_Time']
                

    ''' 
    
    
    import numpy as np
    import scipy.io
    
    from obspy import Trace, Stream
    from obspy.core import AttribDict
    from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader
    
    import pyproj
    import time
    import sys
    
    
###############################################################################
######################        FOR Elevation DATA         ######################
###############################################################################
    
    
    print('')
    print('==> Processing {}'.format(matfile))
    
    # Check on TWT or Elevation Data
    if elevation == True:
        
        # Define variables
        mat             = scipy.io.loadmat(matfile)
        data            = mat['Data']
        X               = mat['X']
        Y               = mat['Y']
        Elevation       = mat['Elevation_WGS84'][0]
        Elevation_max   = mat['Elevation_WGS84'].max()
        num_of_samples  = mat['Elevation_WGS84'].shape[1]  
        gps_time        = mat['GPS_time']
    
        # Create empty stream
        stream = Stream()
        
        # Figure out TWT sampling interval
        diff_el = np.array([])
        
        for i in range(len(Elevation) - 1):
            diff_el = np.append(diff_el, Elevation[i] - Elevation[i+1])
    
        sample_interval_ms = diff_el.mean()
        
        
        # Create Traces and Trace Header
        for i in range(data.shape[1] - 1):
            
            year            = int(time.strftime("%Y", time.gmtime(gps_time[0][i])))
            day_of_year     = int(time.strftime("%j", time.gmtime(gps_time[0][i])))
            hour            = int(time.strftime("%H", time.gmtime(gps_time[0][i])))
            minute          = int(time.strftime("%M", time.gmtime(gps_time[0][i])))
            second          = int(time.strftime("%S", time.gmtime(gps_time[0][i])))
            
            # Create some random data.
            trace_  = data[:, i]
            trace_  = np.require(trace_, dtype=np.float32)
            trace_  = trace_.flatten()
            trace   = Trace(data=trace_)
        
            if differenciate == True:
                trace = trace.differentiate(method='gradient')
                
            else:
                pass
        
        
            trace.stats.delta = 0.01
            # SEGY does not support microsecond precision! Any microseconds will
            # be discarded.
            #trace.stats.starttime = UTCDateTime(2011,11,11,11,11,11)
        
            # If you want to set some additional attributes in the trace header,
            # add one and only set the attributes you want to be set. Otherwise the
            # header will be created for you with default values.
        
            if not hasattr(trace.stats, 'segy.trace_header'):
                trace.stats.segy = {}
            trace.stats.segy.trace_header = SEGYTraceHeader()
        
            trace.stats.segy.trace_header.trace_sequence_number_within_line = i + 1
            trace.stats.segy.trace_header.trace_sequence_number_within_segy_file = i + 1
            trace.stats.segy.trace_header.original_field_record_number = i + 1
            trace.stats.segy.trace_header.trace_number_within_the_original_field_record = i + 1
            trace.stats.segy.trace_header.energy_source_point_number = i + 1
        
            trace.stats.segy.trace_header.ensemble_number = 1
            trace.stats.segy.trace_header.trace_number_within_the_ensemble = 1
            trace.stats.segy.trace_header.trace_identification_code = 1
            trace.stats.segy.trace_header.number_of_vertically_summed_traces_yielding_this_trace = 1
            trace.stats.segy.trace_header.number_of_horizontally_stacked_traces_yielding_this_trace = 1
            trace.stats.segy.trace_header.data_use = 1
            trace.stats.segy.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group = 0
        
            trace.stats.segy.trace_header.receiver_group_elevation = Elevation_max
            trace.stats.segy.trace_header.surface_elevation_at_source = Elevation_max
        
            trace.stats.segy.trace_header.source_depth_below_surface = 0
            trace.stats.segy.trace_header.datum_elevation_at_receiver_group = 36167
            trace.stats.segy.trace_header.datum_elevation_at_source = 89275
            trace.stats.segy.trace_header.water_depth_at_source = 36292
            trace.stats.segy.trace_header.water_depth_at_group = 89581
            trace.stats.segy.trace_header.scalar_to_be_applied_to_all_elevations_and_depths = 1
            trace.stats.segy.trace_header.scalar_to_be_applied_to_all_coordinates = 1
        
            trace.stats.segy.trace_header.source_coordinate_x   = X[0][i].astype(int)
            trace.stats.segy.trace_header.source_coordinate_y   = Y[0][i].astype(int)
            trace.stats.segy.trace_header.group_coordinate_x    = X[0][i].astype(int)
            trace.stats.segy.trace_header.group_coordinate_y    = Y[0][i].astype(int)
        
            trace.stats.segy.trace_header.coordinate_units = 1
            trace.stats.segy.trace_header.weathering_velocity = 300
            trace.stats.segy.trace_header.subweathering_velocity = 168
            trace.stats.segy.trace_header.uphole_time_at_source_in_ms = 38
            trace.stats.segy.trace_header.uphole_time_at_group_in_ms = 0
            trace.stats.segy.trace_header.source_static_correction_in_ms = 0
            trace.stats.segy.trace_header.group_static_correction_in_ms = 0
            trace.stats.segy.trace_header.total_static_applied_in_ms = 0
            trace.stats.segy.trace_header.lag_time_A = 0
            trace.stats.segy.trace_header.lag_time_B = 0
            trace.stats.segy.trace_header.delay_recording_time = 0
            trace.stats.segy.trace_header.mute_time_start_time_in_ms = 0
            trace.stats.segy.trace_header.mute_time_end_time_in_ms = 0
            trace.stats.segy.trace_header.number_of_samples_in_this_trace = num_of_samples
            trace.stats.segy.trace_header.sample_interval_in_ms_for_this_trace = sample_interval_ms
            trace.stats.segy.trace_header.gain_type_of_field_instruments = 1
            trace.stats.segy.trace_header.instrument_gain_constant = 0
            trace.stats.segy.trace_header.instrument_early_or_initial_gain = 0
            trace.stats.segy.trace_header.correlated = 1
            trace.stats.segy.trace_header.sweep_frequency_at_start = 180
            trace.stats.segy.trace_header.sweep_frequency_at_end = 210
            trace.stats.segy.trace_header.sweep_length_in_ms = 10
            trace.stats.segy.trace_header.sweep_type = 4
            trace.stats.segy.trace_header.sweep_trace_taper_length_at_start_in_ms = 3
            trace.stats.segy.trace_header.sweep_trace_taper_length_at_end_in_ms = 0
            trace.stats.segy.trace_header.taper_type = 3
            trace.stats.segy.trace_header.alias_filter_frequency = 0
            trace.stats.segy.trace_header.alias_filter_slope = 0
            trace.stats.segy.trace_header.notch_filter_frequency = 0
            trace.stats.segy.trace_header.notch_filter_slope = 0
            trace.stats.segy.trace_header.low_cut_frequency = 0
            trace.stats.segy.trace_header.high_cut_frequency = 0
            trace.stats.segy.trace_header.low_cut_slope = 0
            trace.stats.segy.trace_header.high_cut_slope = 0
            trace.stats.segy.trace_header.year_data_recorded = year
            trace.stats.segy.trace_header.day_of_year = day_of_year
            trace.stats.segy.trace_header.hour_of_day = hour
            trace.stats.segy.trace_header.minute_of_hour = minute
            trace.stats.segy.trace_header.second_of_minute = second
            trace.stats.segy.trace_header.time_basis_code = 4
            trace.stats.segy.trace_header.trace_weighting_factor = 0
            trace.stats.segy.trace_header.geophone_group_number_of_roll_switch_position_one = 0
            trace.stats.segy.trace_header.geophone_group_number_of_trace_number_one = 1
            trace.stats.segy.trace_header.geophone_group_number_of_last_trace = 1
            trace.stats.segy.trace_header.gap_size = 0
            trace.stats.segy.trace_header.over_travel_associated_with_taper = 0
        
            trace.stats.segy.trace_header.x_coordinate_of_ensemble_position_of_this_trace = X[0][i].astype(int)
            trace.stats.segy.trace_header.y_coordinate_of_ensemble_position_of_this_trace = Y[0][i].astype(int)
        
            trace.stats.segy.trace_header.for_3d_poststack_data_this_field_is_for_in_line_number = 1
        
            trace.stats.segy.trace_header.for_3d_poststack_data_this_field_is_for_cross_line_number = i + 1
            trace.stats.segy.trace_header.shotpoint_number = i + 1
        
            trace.stats.segy.trace_header.scalar_to_be_applied_to_the_shotpoint_number = 1
            trace.stats.segy.trace_header.trace_value_measurement_unit = 1
            trace.stats.segy.trace_header.transduction_constant_mantissa = 1
            trace.stats.segy.trace_header.transduction_constant_exponent = 0
            trace.stats.segy.trace_header.transduction_units = 9
            trace.stats.segy.trace_header.device_trace_identifier = 0
            trace.stats.segy.trace_header.scalar_to_be_applied_to_times = 0
            trace.stats.segy.trace_header.source_type_orientation = 0
            trace.stats.segy.trace_header.source_energy_direction_mantissa = 0
            trace.stats.segy.trace_header.source_energy_direction_exponent = 0
            trace.stats.segy.trace_header.source_measurement_mantissa = 1
            trace.stats.segy.trace_header.source_measurement_exponent = 0
            trace.stats.segy.trace_header.source_measurement_unit = 2
        
            # Add trace to stream
            stream.append(trace)
        
        # A SEGY file has file wide headers. This can be attached to the stream
        # object.  If these are not set, they will be autocreated with default
        # values.
        stream.stats = AttribDict()
        stream.stats.textual_file_header = 'Textual Header!'
        stream.stats.binary_file_header = SEGYBinaryFileHeader()
        stream.stats.binary_file_header.trace_sorting_code = 5
        stream.stats.binary_file_header.number_of_samples_per_data_trace = num_of_samples
        stream.stats.binary_file_header.number_of_samples_per_data_trace_for_original_field_recording = num_of_samples
        stream.stats.binary_file_header.data_sample_format_code = 5
        
        # stream = stream.differentiate(method='gradient')
        
        
        if differenciate == True:
            stream.write(matfile.split('.')[0] + '_diff.sgy', format='SEGY', data_encoding=5, byteorder='>',textual_file_encoding='ASCII')
            print('==> Written: {}_diff.sgy'.format(matfile.split('.')[0]))
            
        else:    
            stream.write(matfile.split('.')[0] + '.sgy', format='SEGY', data_encoding=5, byteorder='>',textual_file_encoding='ASCII')
            print('==> Written: {}.sgy'.format(matfile.split('.')[0]))
        #return stream
        # Data Encoding = 4 byte IEEE floating points (float32)
        # Byte Order = big endian
        # print('Saved stream as : {}.sgy'.format(df))


###############################################################################
######################           FOR TWT DATA            ######################
###############################################################################
    
    elif elevation == False:

        mat = scipy.io.loadmat(matfile)
        
        # Reproject
        
        # Check on Region and if re-projection is required
        if region == '':
            print('')
            print('!!! Please prvide a region ==> Greenland or Antarctica')
            print('')
            sys.exit()
        
        elif region == 'Greenland':
            Projection = "+init=EPSG:3413"
            
        elif region == 'Antarctica':
            Projection = "+init=EPSG:3031"
        
        else:
            print('')
            print('==>> Something is Wrong with the Region, this might not end well...') 
            print("==>> Leave region='' if you have the correct X and Y coordinates")
            print("==>> or set region='Greenland' or 'Antarctica'")
            print('')
            sys.exit()
        
        EPSG_=pyproj.Proj(Projection)
        X, Y = EPSG_(mat['Longitude'], mat['Latitude'])
        
        data            = mat['Data']
        X               = X
        Y               = Y
        TWT             = mat['Time'].T[0]
        num_of_samples  = mat['Time'].shape[0]  
        gps_time        = mat['GPS_time']
     
        
        stream = Stream()
        
 
        # Figure out TWT sampling interval
        diff_TWT = np.array([])
        
        for i in range(len(TWT) - 1):
            diff_TWT = np.append(diff_TWT, TWT[i+1] - TWT[i])
    
        sample_interval_ms = diff_TWT.mean()
        
        for i in range(data.shape[1] - 1):
            
            year            = int(time.strftime("%Y", time.gmtime(gps_time[0][i])))
            day_of_year     = int(time.strftime("%j", time.gmtime(gps_time[0][i])))
            hour            = int(time.strftime("%H", time.gmtime(gps_time[0][i])))
            minute          = int(time.strftime("%M", time.gmtime(gps_time[0][i])))
            second          = int(time.strftime("%S", time.gmtime(gps_time[0][i])))
            
            # Create some random data.
            fish = data[:, i]
            fish = np.require(fish, dtype=np.float32)
            fish = fish.flatten()
            trace = Trace(data=fish)
        
            # trace = trace.differentiate(method='gradient')
        
        
            trace.stats.delta = 0.01
            # SEGY does not support microsecond precision! Any microseconds will
            # be discarded.
            #trace.stats.starttime = UTCDateTime(2011,11,11,11,11,11)
        
            # If you want to set some additional attributes in the trace header,
            # add one and only set the attributes you want to be set. Otherwise the
            # header will be created for you with default values.
        
            if not hasattr(trace.stats, 'segy.trace_header'):
                trace.stats.segy = {}
            trace.stats.segy.trace_header = SEGYTraceHeader()
        
            trace.stats.segy.trace_header.trace_sequence_number_within_line = i + 1
            trace.stats.segy.trace_header.trace_sequence_number_within_segy_file = i + 1
            trace.stats.segy.trace_header.original_field_record_number = i + 1
            trace.stats.segy.trace_header.trace_number_within_the_original_field_record = i + 1
            trace.stats.segy.trace_header.energy_source_point_number = i + 1
        
            trace.stats.segy.trace_header.ensemble_number = 1
            trace.stats.segy.trace_header.trace_number_within_the_ensemble = 1
            trace.stats.segy.trace_header.trace_identification_code = 1
            trace.stats.segy.trace_header.number_of_vertically_summed_traces_yielding_this_trace = 1
            trace.stats.segy.trace_header.number_of_horizontally_stacked_traces_yielding_this_trace = 1
            trace.stats.segy.trace_header.data_use = 1
            trace.stats.segy.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group = 0
        
            trace.stats.segy.trace_header.receiver_group_elevation = 0
            trace.stats.segy.trace_header.surface_elevation_at_source = 0
        
            trace.stats.segy.trace_header.source_depth_below_surface = 0
            trace.stats.segy.trace_header.datum_elevation_at_receiver_group = 36167
            trace.stats.segy.trace_header.datum_elevation_at_source = 89275
            trace.stats.segy.trace_header.water_depth_at_source = 36292
            trace.stats.segy.trace_header.water_depth_at_group = 89581
            trace.stats.segy.trace_header.scalar_to_be_applied_to_all_elevations_and_depths = 1
            trace.stats.segy.trace_header.scalar_to_be_applied_to_all_coordinates = 1
        
            trace.stats.segy.trace_header.source_coordinate_x   = X[0][i].astype(int)
            trace.stats.segy.trace_header.source_coordinate_y   = Y[0][i].astype(int)
            trace.stats.segy.trace_header.group_coordinate_x    = X[0][i].astype(int)
            trace.stats.segy.trace_header.group_coordinate_y    = Y[0][i].astype(int)
        
            trace.stats.segy.trace_header.coordinate_units = 1
            trace.stats.segy.trace_header.weathering_velocity = 300
            trace.stats.segy.trace_header.subweathering_velocity = 168
            trace.stats.segy.trace_header.uphole_time_at_source_in_ms = 38
            trace.stats.segy.trace_header.uphole_time_at_group_in_ms = 0
            trace.stats.segy.trace_header.source_static_correction_in_ms = 0
            trace.stats.segy.trace_header.group_static_correction_in_ms = 0
            trace.stats.segy.trace_header.total_static_applied_in_ms = 0
            trace.stats.segy.trace_header.lag_time_A = 0
            trace.stats.segy.trace_header.lag_time_B = 0
            trace.stats.segy.trace_header.delay_recording_time = 0
            trace.stats.segy.trace_header.mute_time_start_time_in_ms = 0
            trace.stats.segy.trace_header.mute_time_end_time_in_ms = 0
            trace.stats.segy.trace_header.number_of_samples_in_this_trace = num_of_samples
            trace.stats.segy.trace_header.sample_interval_in_ms_for_this_trace = sample_interval_ms
            trace.stats.segy.trace_header.gain_type_of_field_instruments = 1
            trace.stats.segy.trace_header.instrument_gain_constant = 0
            trace.stats.segy.trace_header.instrument_early_or_initial_gain = 0
            trace.stats.segy.trace_header.correlated = 1
            trace.stats.segy.trace_header.sweep_frequency_at_start = 180
            trace.stats.segy.trace_header.sweep_frequency_at_end = 210
            trace.stats.segy.trace_header.sweep_length_in_ms = 10
            trace.stats.segy.trace_header.sweep_type = 4
            trace.stats.segy.trace_header.sweep_trace_taper_length_at_start_in_ms = 3
            trace.stats.segy.trace_header.sweep_trace_taper_length_at_end_in_ms = 0
            trace.stats.segy.trace_header.taper_type = 3
            trace.stats.segy.trace_header.alias_filter_frequency = 0
            trace.stats.segy.trace_header.alias_filter_slope = 0
            trace.stats.segy.trace_header.notch_filter_frequency = 0
            trace.stats.segy.trace_header.notch_filter_slope = 0
            trace.stats.segy.trace_header.low_cut_frequency = 0
            trace.stats.segy.trace_header.high_cut_frequency = 0
            trace.stats.segy.trace_header.low_cut_slope = 0
            trace.stats.segy.trace_header.high_cut_slope = 0
            trace.stats.segy.trace_header.year_data_recorded = year
            trace.stats.segy.trace_header.day_of_year = day_of_year
            trace.stats.segy.trace_header.hour_of_day = hour
            trace.stats.segy.trace_header.minute_of_hour = minute
            trace.stats.segy.trace_header.second_of_minute = second
            trace.stats.segy.trace_header.time_basis_code = 4
            trace.stats.segy.trace_header.trace_weighting_factor = 0
            trace.stats.segy.trace_header.geophone_group_number_of_roll_switch_position_one = 0
            trace.stats.segy.trace_header.geophone_group_number_of_trace_number_one = 1
            trace.stats.segy.trace_header.geophone_group_number_of_last_trace = 1
            trace.stats.segy.trace_header.gap_size = 0
            trace.stats.segy.trace_header.over_travel_associated_with_taper = 0
        
            trace.stats.segy.trace_header.x_coordinate_of_ensemble_position_of_this_trace = X[0][i].astype(int)
            trace.stats.segy.trace_header.y_coordinate_of_ensemble_position_of_this_trace = Y[0][i].astype(int)
        
            trace.stats.segy.trace_header.for_3d_poststack_data_this_field_is_for_in_line_number = 1
        
            trace.stats.segy.trace_header.for_3d_poststack_data_this_field_is_for_cross_line_number = i + 1
            trace.stats.segy.trace_header.shotpoint_number = i + 1
        
            trace.stats.segy.trace_header.scalar_to_be_applied_to_the_shotpoint_number = 1
            trace.stats.segy.trace_header.trace_value_measurement_unit = 1
            trace.stats.segy.trace_header.transduction_constant_mantissa = 1
            trace.stats.segy.trace_header.transduction_constant_exponent = 0
            trace.stats.segy.trace_header.transduction_units = 9
            trace.stats.segy.trace_header.device_trace_identifier = 0
            trace.stats.segy.trace_header.scalar_to_be_applied_to_times = 0
            trace.stats.segy.trace_header.source_type_orientation = 0
            trace.stats.segy.trace_header.source_energy_direction_mantissa = 0
            trace.stats.segy.trace_header.source_energy_direction_exponent = 0
            trace.stats.segy.trace_header.source_measurement_mantissa = 1
            trace.stats.segy.trace_header.source_measurement_exponent = 0
            trace.stats.segy.trace_header.source_measurement_unit = 2
        
        
        
            # Add trace to stream
            stream.append(trace)
        
        # A SEGY file has file wide headers. This can be attached to the stream
        # object.  If these are not set, they will be autocreated with default
        # values.
        stream.stats = AttribDict()
        stream.stats.textual_file_header = 'Textual Header!'
        stream.stats.binary_file_header = SEGYBinaryFileHeader()
        stream.stats.binary_file_header.trace_sorting_code = 5
        stream.stats.binary_file_header.number_of_samples_per_data_trace = num_of_samples
        stream.stats.binary_file_header.number_of_samples_per_data_trace_for_original_field_recording = num_of_samples
        stream.stats.binary_file_header.data_sample_format_code = 5
        
        # stream = stream.differentiate(method='gradient')
        
        
        
        if differenciate == True:
            stream.write(matfile.split('.')[0] + '_diff.sgy', format='SEGY', data_encoding=5, byteorder='>',textual_file_encoding='ASCII')
            print('==> Written: {}_diff.sgy'.format(matfile.split('.')[0]))
            
        else:    
            stream.write(matfile.split('.')[0] + '.sgy', format='SEGY', data_encoding=5, byteorder='>',textual_file_encoding='ASCII')
            print('==> Written: {}.sgy'.format(matfile.split('.')[0]))
        #return stream
        # Data Encoding = 4 byte IEEE floating points (float32)
        # Byte Order = big endian
        # print('Saved stream as : {}.sgy'.format(df))


