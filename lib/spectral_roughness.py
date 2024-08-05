
def spectral_roughness( Line_ID='',
                        Longitude='',
                        Latitude='',
                        Bed_elevation='',
                        filename = '',
                        EPSG='',
                        make_Xi=True,
                        make_Eta=True,
                        make_RMS_h=True,
                        make_RMS_d=True,
                        N_list=[],
                        step=1,
                        tfl=2**13,
                        spacing_limit=50,
                        out_format=['csv']
                        ):

    ''' 

            filename                     : any filename (.csv, .shp, ... will be appended)
            region                       : 'Antarctica' (EPSG:3031) or 'Greenland' (EPSG:3413)
            make_[Xi, Eta, RMS_h, RMS_d] : can be set True or False
            N_list                       : any integer < 5 
            step                         : usually 1
            tfl                          : 2**Y (Y=13 usually)
            out_format                   : 'csv', 'shapefile', 'geojson'
    '''

    import pandas as pd 
    import numpy as np
    import geopandas
    from numpy import trapz
    from scipy.signal import detrend
    from scipy.fftpack import fft
    import geopy.distance
    from pyproj import Transformer
    import os
    
    print('Calculating the following parameters:')
    if make_Xi==True:
        print('==> Xi (spectral vertical roughness)')
    if make_Eta==True:
        print('==> Eta (spectral horizontal roughness)')
    if make_RMS_h==True:
        print('==> RMS height')
    if make_RMS_d==True:
        print('==> RMS deviation')
    print('')
    print("For the following N's: {}".format(N_list))
    print('')
    print('==> Reading Segment {}'.format(Line_ID))


    ####################################################################
    ##################
    #   Coordinates
    ##################

    EPSG = EPSG

    transformer = Transformer.from_crs(4326, EPSG)
    X, Y        = transformer.transform(Longitude, Latitude)

    #######################################################################
    #########################
    #  Distance and Spacing
    #########################

    spacing = np.array([])
    
    for i in range(1, len(Longitude)-1):
        coord_1    = (Latitude[i], Longitude[i])
        coord_2    = (Latitude[i + 1], Longitude[i + 1])
        spacing_i  = geopy.distance.geodesic(coord_1, coord_2).meters
        spacing    = np.append(spacing, spacing_i)

    distance       = np.cumsum(spacing)        # cumulative sum of spacing
    distance       = np.insert(distance, 0, 0) # insert zero at first location
    spacing        = np.insert(spacing, 0, spacing.mean()) # insert mean spacing at first position
    spacing        = np.insert(spacing, len(spacing), spacing.mean()) # insert mean spacing at last position

    # df['Spacing']  = pd.DataFrame(spacing)
    # df['Distance'] = pd.DataFrame(distance)


    if make_Xi == True:
        XI_list    = []
    if make_Eta == True:
        ETA_list   = []
    if make_RMS_h == True:
        RMS_h_list = []
    if make_RMS_d == True:
        RMS_d_list = []

    # The N Loop
    for N in N_list:

        print('')
        print('Loop N={}'.format(N))

        # moving window exponent
        N = N

        # steps (1 for full analysis)
        step = step

        # distance between data points
        #spacing = spacing.mean()

        # X coordinate
        X = X
        # Y coordinate
        Y = Y

        # Z coordinate (bed elevation)
        Z = Bed_elevation

        # Length of profile
        length_profile = len(Bed_elevation)

        # moving window (depends on exponent)
        MW = 2**N

        # half of moving window
        MW_half = MW/2

        # array for loop (depends on length of data frame and MW size)
        walk = np.arange(MW_half, length_profile - MW_half, 1).astype(int)


        if make_Xi == True:
            XI    = []

        if make_Eta == True:
            ETA   = []

        if make_RMS_h == True:
            RMS_h = []

        if make_RMS_d == True:
            RMS_d = []

        ###################################################################
        ###################
        # CALCULATE VALUES
        ###################

        for i in walk[0:len(walk)]:
            # define MW boundaries in data
            left, right = [int(i-MW_half), int(i+MW_half)]

            # slice accordingly
            Zi = Z[left:right]

            ##########
            ### XI ###
            ##########

            if make_Xi == True:
                try:
                    # define spacing between data points for the MW
                    ds = spacing[left:right].mean()

                    # detrend profile
                    Zdi = detrend(Zi)

                    # length of the MW in meters
                    Lmw = (MW - 1) * ds

                    # FFT
                    Zspec = abs(fft(Zdi, n=tfl)) / MW_half

                    # truncates array by half
                    Zspec = Zspec[0:int(tfl/2)]

                    # Power Zspec by ^2
                    PZspec = Zspec**2

                    # Define sample points for integration
                    SP = np.arange(0, tfl/2, 1) / Lmw / (tfl/MW)

                    # Integrate to calculate Xi (flip necessary to prevent neg. values)
                    xi = 2 * np.pi * trapz(np.flipud(SP), np.flipud(PZspec))

                except:
                    xi = np.nan

                XI.append(xi)



            ###########
            ### ETA ###
            ###########

            if make_Eta == True:
                try:
                    # gradient of the detrended Z profile
                    Zdi_grad = np.gradient(Zdi,ds);

                    # FFT
                    Zspec_grad = abs(fft(Zdi_grad, n=tfl)) / MW_half

                    # truncates array by half
                    Zspec_grad = Zspec_grad[0:int(tfl/2)]

                    # Power Zspec_grad by ^2
                    PZspec_grad = Zspec_grad**2 * Lmw / MW

                    # Define sample points for integration
                    SP_grad = np.arange(0, tfl/2, 1) / (tfl/MW)  # dfsi=(0:sN/2-1)/(sN/Ni);

                    # Integrate to calculate Xi (flip necessary to prevent neg. values)
                    xi_grad = 2 * np.pi * trapz(np.flipud(SP_grad), np.flipud(PZspec_grad))

                    # Calculate Eta & append
                    eta = xi / xi_grad
                except:
                    eta = np.nan

                ETA.append(eta)


            ##############
            # RMS Height
            ##############
            if make_RMS_h == True:
                try:
                    RMS_HI = []
                    for j in np.arange(0, MW):
                        rms_hi = ( Zdi[j] - Zdi.mean() )**2
                        RMS_HI.append(rms_hi)

                    RMS_HI_mean  = np.array(RMS_HI).mean()           
                    rms_h        = np.sqrt(RMS_HI_mean)
                except:
                    rms_h = np.nan

                RMS_h.append(rms_h)


            ################
            # RMS Deviation
            ################

            if make_RMS_d == True:
                try:
                    RMS_DI = []
                    RMS_DH = []

                    for j in np.arange(0, MW - 1):
                        rms_di = ( Zdi[j] - Zdi[j+1] )**2
                        rms_dh = Zdi[j] - Zdi[j+1]
                        RMS_DI.append(rms_di)
                        RMS_DH.append(rms_dh)

                    RMS_DI_mean  = np.array(RMS_DI).mean()           
                    rms_d        = np.sqrt(RMS_DI_mean)
                except:
                    rms_d = np.nan

                RMS_d.append(rms_d)


            '''
            ##################
            # Hurst Exponent
            ##################

            Hz      = np.cumsum(np.abs(RMS_DH))

            Hx      = np.array(np.log10(np.cumsum(df['Spacing'][left:right])))
            Hx      = Hx[2::]
            Hx      = np.insert(Hx, 0, 0)

            Hy      = np.log10(Hz)
            Hs, Hi  = np.polyfit(np.array(Hx), np.array(Hy), 1)
            try:
                Hy      = np.log10(Zdi)
                Hs, Hi  = np.polyfit(np.array(Hx), np.array(Hy), 1)
            except:    
                Hy      = np.log10(Zdi + np.abs(Zdi.min()) + 1)
                Hs, Hi  = np.polyfit(np.array(Hx), np.array(Hy), 1)


            HURST.append(Hs)
            '''
        ############################################################################
        ##################
        # GATHER DATA
        ##################

        if make_Xi == True:
            # fill with nan's where there is no data for xi and eta
            Xi      = np.array(XI)
            Xi      = np.append(Xi, np.repeat(np.nan, MW/2))
            Xi      = np.append(np.flipud(Xi), np.repeat(np.nan, MW/2))
            Xi      = np.flipud(Xi)
            xi_name = 'Xi_N_' + str(N)

            XI_list.append([Xi, xi_name])

        if make_Eta == True:
            Eta      = np.array(ETA)
            Eta      = np.append(Eta, np.repeat(np.nan, MW/2))
            Eta      = np.append(np.flipud(Eta), np.repeat(np.nan, MW/2))
            Eta      = np.flipud(Eta)
            eta_name = 'Eta_N_' + str(N)

            ETA_list.append([Eta, eta_name])

        if make_RMS_h == True:
            RMS_h      = np.array(RMS_h)
            RMS_h      = np.append(RMS_h, np.repeat(np.nan, MW/2))
            RMS_h      = np.append(np.flipud(RMS_h), np.repeat(np.nan, MW/2))
            RMS_h      = np.flipud(RMS_h)
            rms_h_name = 'RMS_h_N_' + str(N)

            RMS_h_list.append([RMS_h, rms_h_name])

        if make_RMS_d == True:
            RMS_d      = np.array(RMS_d)
            RMS_d      = np.append(RMS_d, np.repeat(np.nan, MW/2))
            RMS_d      = np.append(np.flipud(RMS_d), np.repeat(np.nan, MW/2))
            RMS_d      = np.flipud(RMS_d)
            rms_d_name = 'RMS_d_N_' + str(N)

            RMS_d_list.append([RMS_d, rms_d_name])


        print('... done.')
        
    #################################################################################
    ##############################
    # Build columns for Data Frame
    ##############################

    df             = pd.DataFrame(Longitude)
    df['Latitude'] = pd.DataFrame(Latitude)
    df.columns     = ['Longitude', 'Latitude']
    df['Line_ID']  = Line_ID
    df             = df[['Line_ID', 'Longitude', 'Latitude']]
    
    df['Spacing']      = pd.DataFrame(spacing)
    df['Distance']      = pd.DataFrame(distance)
    df['Bed_elevation'] = pd.DataFrame(Bed_elevation)
    
    if make_Xi == True:
        for N in np.arange(0, len(N_list)):
            data     = np.array(XI_list[N][0])
            name     = XI_list[N][1]
            df[name] = pd.DataFrame(data)

    if make_Eta == True:
        for N in np.arange(0, len(N_list)):
            data     = np.array(ETA_list[N][0])
            name     = ETA_list[N][1]
            df[name] = pd.DataFrame(data)

    if make_RMS_h == True:
        for N in np.arange(0, len(N_list)):
            data     = np.array(RMS_h_list[N][0])
            name     = RMS_h_list[N][1]
            df[name] = pd.DataFrame(data)

    if make_RMS_d == True:
        for N in np.arange(0, len(N_list)):
            data     = np.array(RMS_d_list[N][0])
            name     = RMS_d_list[N][1]
            df[name] = pd.DataFrame(data)


    ###############################################################################
    ###################
    # CLEAN AROUND GAPS
    ###################

    print("")
    print("==> Removing values around spacing gaps of more than {} m".format(spacing_limit))
    spacing_limit = 50
    gap_idxs      = df[df["Spacing"] > spacing_limit].index.values

    cleanup_list = []
    for N in N_list:
        cols = [col for col in df.columns if "N_{}".format(N) in col]
        gate = 2**N
        cleanup_list.append([N, gate, cols])

    for line in cleanup_list:
        win = int( (line[1] / 2) + 2 )
        for col_name in line[2]:
            for idx in gap_idxs:
                df[col_name][idx-win:idx+win] = np.nan

    print("... done.")




    ###############################################################################        
    ##################
    # WRITE DATA
    ##################
    print("")
    print("==> Writing {} data.".format(out_format))
    
    # Export to .csv
    if 'csv' in out_format:
        if not os.path.exists('csv'):
            os.makedirs('csv')

        # write output to csv file
        df.to_csv('csv/' + filename + '.csv', sep='\t', index=False)
        print('==> Saved: csv/{}.csv'.format(filename))


    #   Export as shape/geojson
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Longitude, df.Latitude))

    if 'shapefile' in out_format:
        if not os.path.exists('shapes'):
            os.makedirs('shapes')

        gdf.to_file('shapes/' + filename + '.shp')
        print('==> Saved: shapes/{}.shp'.format(filename))

    if 'geojson' in out_format: 
        if not os.path.exists('geojson'):
            os.makedirs('geojson')

        gdf.to_file('geojson/' + filename + '.geojson', driver='GeoJSON')
        print('==> Saved: geojson/{}.geojson'.format(filename))

    print("... done.")

    ######################################################################    
    
    ################
    # OVERVIEW PLOT
    ################

    '''if not os.path.exists('figures'):
        os.makedirs('figures')


    fig = plt.subplots(figsize=(12, 12))
    gs  = gridspec.GridSpec(6, 1,
                            width_ratios=[1],
                            height_ratios=[1, 1, 1, 1, 1, 1])

    gs.update(wspace=0.0, hspace=0.5)

    ax0 = plt.subplot(gs[0])
    plt.plot(df['Z'])
    plt.title('Elevation Profile')

    ax1 = plt.subplot(gs[1])
    plt.plot(df['Spacing'])
    plt.title('Spacing')

    ax2 = plt.subplot(gs[2])
    plt.plot(df['RMS_h'])
    plt.title('RMS Height')

    ax3 = plt.subplot(gs[3])
    plt.plot(df['RMS_d'])
    plt.title('RMS Deviation')

    #ax4 = plt.subplot(gs[4])
    #plt.plot(df['Hurst'])
    #plt.title('Hurst Exponent')

    ax4 = plt.subplot(gs[4])
    plt.plot(df['Xi'])
    plt.title('Vertical Roughness')

    ax5 = plt.subplot(gs[5])
    plt.plot(df['Eta'])
    plt.title('Horizontal Roughness')

    plt.savefig(filename + '_overview.png', dpi=100, bbox_inches='tight')'''

    return df
