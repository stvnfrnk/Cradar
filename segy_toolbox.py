# Toolbox to deal with SEGY Data

################################
# Quickplot a radar section
# from a SEGY file
################################

def plot_segy(segy_input, png=False, cmap='bone_r', log10=False):
    
    '''
    
    
    '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    from obspy.io.segy.segy import _read_segy
    
    stream  = _read_segy(segy_input, headonly=True)
    data    = np.stack(t.data for t in list(stream.traces))
    
    fig = plt.subplots(figsize=(20,10))
    if log10 == False:
        plt.imshow(data.T, cmap=cmap, aspect='auto')
    elif log10 == True:
        plt.imshow(np.log10(data).T, cmap=cmap, aspect='auto')
    plt.colorbar()
    plt.show()

    if png == True:
        plt.savefig('test_segy_plot.png', dpi=200, bbox_inches='tight')
        print('Saved Figure as ==> test_segy_plot.png')
        
    else:
        return fig






def radar2segy(data='', 
               receiver_elevation=0, 
               num_of_samples=0, 
               sample_interval=0,
               X='', 
               Y='',
               step=1,
               time_mode='gmtime',
               gps_time='',
               year='',
               day_of_year='',
               hour='',
               minute='',
               second='',
               differenciate=''
               ):

    '''  
    ==>    Writes Cradar Object as SEGY-Format File.
    Usage: write_segy(self, region='', differenciate=False, step=1)

    - Parameters are:
        1) Cradar Object
        2) region='Greenland' or 'Antarctica'
        3) differenciate=True/False (False is default)
        4) step


            1) Cradar Object, either in 'twt' or 'Z' Domain
            2) if the parameter 'region' is set empty, you already
               should have 'X' and 'Y' coordinates available for your region
               --> if you set it to 'Greenland'  they will be projected to EPSG:3413
               --> if you set it to 'Antarctica' they will be projected to EPSG:3031
            3) choose weather to take the data as it is or to differenciate
               --> set differenciate=True (default is differenciate=False)
            4) step=N, every N'th trace will be considered (to reduce data size if needed)
    ''' 

    import numpy as np

    from obspy import Trace, Stream
    from obspy.core import AttribDict
    from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader

    from pyproj import Transformer
    import time
    import sys


    data = data
    step = step

    gps_time = gps_time

    receiver_elevation = receiver_elevation
    num_of_samples     = num_of_samples
    sample_interval   = sample_interval

    X = X
    Y = Y

    time_array   = year
    year         = year
    day_of_year  = day_of_year
    hour         = hour
    minute       = minute
    second       = second

    differenciate = differenciate


    # Create empty stream
    stream = Stream()

    # Create Traces and Trace Header
    for i in range(0, len(X) - 1, step):

        if time_mode == 'gmtime':

            year            = int(time.strftime("%Y", time.gmtime(gps_time[i])))
            day_of_year     = int(time.strftime("%j", time.gmtime(gps_time[i])))
            hour            = int(time.strftime("%H", time.gmtime(gps_time[i])))
            minute          = int(time.strftime("%M", time.gmtime(gps_time[i])))
            second          = int(time.strftime("%S", time.gmtime(gps_time[i])))

        else:
            year         = year[i]
            day_of_year  = day_of_year[i]
            hour         = hour[i]
            minute       = minute[i]
            second       = second[i]

        # Create some random data.
        trace_  = np.array(data[:, i])
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

        trace.stats.segy.trace_header.receiver_group_elevation = int(receiver_elevation)
        trace.stats.segy.trace_header.surface_elevation_at_source = int(receiver_elevation)

        trace.stats.segy.trace_header.source_depth_below_surface = 0
        trace.stats.segy.trace_header.datum_elevation_at_receiver_group = 36167
        trace.stats.segy.trace_header.datum_elevation_at_source = 89275
        trace.stats.segy.trace_header.water_depth_at_source = 36292
        trace.stats.segy.trace_header.water_depth_at_group = 89581
        trace.stats.segy.trace_header.scalar_to_be_applied_to_all_elevations_and_depths = 1
        trace.stats.segy.trace_header.scalar_to_be_applied_to_all_coordinates = 1

        trace.stats.segy.trace_header.source_coordinate_x   = X[i].astype(int)
        trace.stats.segy.trace_header.source_coordinate_y   = Y[i].astype(int)
        trace.stats.segy.trace_header.group_coordinate_x    = X[i].astype(int)
        trace.stats.segy.trace_header.group_coordinate_y    = Y[i].astype(int)

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
        trace.stats.segy.trace_header.sample_interval_in_ms_for_this_trace = sample_interval
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

        trace.stats.segy.trace_header.x_coordinate_of_ensemble_position_of_this_trace = X[i].astype(int)
        trace.stats.segy.trace_header.y_coordinate_of_ensemble_position_of_this_trace = Y[i].astype(int)

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

    return stream







