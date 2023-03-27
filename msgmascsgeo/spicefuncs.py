"""Utilities function to calculate geoemtric variables via spiceypy."""
import spiceypy
import numpy as np
import pandas as pd


calculated_cols = ['trgepc','tarlon', 'tarlat','taralt', 'tardis', 'tarang ',
 'trgenpc','phase', 'incdnc', 'emissn','hr', 'mn', 'sc','sublon',
 'sublat', 'subalt', 'sunlon', 'sunlat', 'sunalt']

def geometry_intersect(method,
                       target,
                       et,  # calculated from MASCS data via spiceypy.utc2et 
                       pid, # index to aggretate multiple iFOV easily
                       target_frame, 
                       aberration_correction,
                       observer,
                       sensor_frame,
                       radii,
                       tarid,
                       vector):
    """Calculate instrument intersection defined in target_frame with target body at time et.
    """
    # try:

    re = radii[0]
    rp = radii[2]
    radii_relative_difference = (re - rp) / re # radii relative difference
    ######################################################################
    # intersection with the target, if the target is intersected we then 
    # compute the illumination angles using SPICE APIs: SINCPT and ILLUMF

    (spoint, trgepc, srfvec) = spiceypy.sincpt(
         method,target, et, target_frame, 
        'LT+S', observer, sensor_frame, vector)

    (tarlon, tarlat, taralt) = spiceypy.recgeo(spoint, re, radii_relative_difference)
    tardis = spiceypy.vnorm(srfvec)

    ######################################################################
    # Angular diameter

    tarang = np.degrees(2 * np.arctan( max(radii) / spiceypy.vnorm(spoint + srfvec) ))

    ######################################################################
    # Compute the illumination angles (phase, incidence, and emission) 
    # at a specified point on a target body

    (trgenpc, srfvec, phase, incdnc, emissn, visiblef, iluminatedf) = \
         spiceypy.illumf(method, target, 'SUN', et, target_frame,aberration_correction, observer, spoint)

    phase *= spiceypy.dpr()
    incdnc *= spiceypy.dpr()
    emissn *= spiceypy.dpr()

    ######################################################################
    # local solar time for an object on the surface of a body at 
    # a specified longitude
    (hr, mn, sc, ltime, ampm) = \
        spiceypy.et2lst(et, tarid, tarlon, 'PLANETOCENTRIC', 80, 80)

    ######################################################################
    # rectangular coordinates of the sub-observer point on a target body

    (subobserver_spoint, trgepc, srfev) = \
        spiceypy.subpnt(f'INTERCEPT/{method}', target, et, target_frame, aberration_correction, observer)

    (sublon, sublat, subalt) = spiceypy.recgeo(subobserver_spoint, re, radii_relative_difference)

    ######################################################################
    # Subsolar point
    (subsolar_spoint, trgepc, srfev) = \
        spiceypy.subslr(f'INTERCEPT/{method}', target, et, target_frame, aberration_correction, observer)

    (sunlon, sunlat, sunalt) = spiceypy.recgeo(subsolar_spoint, re, radii_relative_difference)

            # transform recgeo output from radians to degrees
    tarlon = np.rad2deg(tarlon)
    tarlat = np.rad2deg(tarlat)
    sublon = np.rad2deg(sublon)
    sublat = np.rad2deg(sublat)
    sunlon = np.rad2deg(sunlon)
    sunlat = np.rad2deg(sunlat)
    
    return {
    'method':method,
    'et':et,
    'pid':pid,
    'spoint':spoint,
    'trgepc':trgepc,
    'srfvec':srfvec,
    'tarlon':tarlon,
    'tarlat':tarlat,
    'taralt':taralt,
    'tardis':tardis,
    'tarang ':tarang ,
    'trgenpc':trgenpc,
    'phase':phase,
    'incdnc':incdnc,
    'emissn':emissn,
    'visiblef':visiblef,
    'iluminatedf':iluminatedf,
    'hr':hr,
    'mn':mn,
    'sc':sc,
    'ltime':ltime,
    'ampm':ampm,
    'sublon':sublon,
    'sublat':sublat,
    'subalt':subalt,
    'srfev':srfev,
    'sunlon':sunlon,
    'sunlat':sunlat,
    'sunalt':sunalt
    }
    # except Exception as e:
    #     return {'NotFound':1,
    #             'Exception':str(e)}



def row_to_shape(row, starting_fovs_n=3,radii=0.0,tarid='',fov_points=[],method='',target='',target_frame='',
                       aberration_correction='',
                       observer='',
                       sensor_frame='',
                 ):
    """Calculate the iFOVs for a given row in the MASCS meas list.
     Row must contain ET time, the pid from original data, used to aggreagate the iFOVs.
     Optional : int_time, if not given we assume 20, in 0.05s units, == 1 second integration.
    
    Input
    =====
    row: Pandas.Dataframe or dict type
         data for a single MASCS measurement.
    
    """
    import shapely
    ## initialize variables
    pid    = row['pid']
    shapes = {}
    ratio  = 0
    internal_points_ls = []
    # time variable used to calculate the FOVS
    # time_string = 'ET_to_UTC'
    time_string = 'ET'
    # this is used to store the calculated surface points and to calculate iFOV/unions(FOVs) 
    # inital time_steps
    if 'int_time' in row:
        int_time = row['int_time'] * 0.05 # convert to seconds
    else : 
        int_time = 20 * 0.05
    
    # time to calculate FOVs : shift based on MET
    time_steps = np.linspace(0,int_time,starting_fovs_n)
    for time_step_ind,time_step in enumerate(time_steps):
        # time = row[time_string]
        # et = spiceypy.utc2et(str(time))
        et = row[time_string]+time_step
        # print(f'{time_step_ind=} {time_step=:.3f} row[{time_string}]={row[time_string]:<21.20} {et=} row[MET]={row["MET(s.ms)"]}',)
        internal_points_df = pd.DataFrame.from_dict([geometry_intersect( method, target, et, pid, target_frame, aberration_correction, observer, sensor_frame, radii, tarid, v) for v in fov_points] )
        internal_points_df['time_step'] = time_step
        internal_points_ls.append(internal_points_df)
        try:
            shape = shapely.geometry.Polygon(internal_points_df.apply(lambda x : (x['tarlon'],x['tarlat']),axis=1).to_list() )    
            shapes[time_step] = [pid,shape]+internal_points_df[calculated_cols].std().to_list()+internal_points_df[calculated_cols].mean().to_list()
        except KeyError:
            rich.print(row)
            ipdb.set_trace()
    ##############################
    # calculate iFOV/unions(FOVs) 
    total_points = len(fov_points)
    # first FOV - 1/4 of a step (back point)
    first_FOV_back_point = internal_points_ls[0].loc[[3*total_points//4],'spoint'].values[0]
    # last FOV - 3/4 of a step (forwad point)
    last_FOV_forwad_point = internal_points_ls[-1].loc[[total_points//4],'spoint'].values[0]

    # calculate middle FOV angular size
    middle_FOV_points = internal_points_ls[1].loc[[internal_points_ls[1]['spoint'].shape[0]//4,3*internal_points_ls[1]['spoint'].shape[0]//4],'spoint'].values

    # calculate ratio (last - first) / middle
    ratio = 2*int(np.rad2deg(spiceypy.vsep(first_FOV_back_point,last_FOV_forwad_point))/np.rad2deg(spiceypy.vsep(*middle_FOV_points)))
    # round to closest bigger odd if ratio is even
    # print(f'{pid=},{len(internal_points_ls)=},{ratio=},D={np.rad2deg(spiceypy.vsep(first_FOV_back_point,last_FOV_forwad_point)):0.4f},d/2={np.rad2deg(spiceypy.vsep(*middle_FOV_points))/2:0.4f}',end=', ')
    ratio = ratio if ratio % 2 else ratio +1
    # print(f'{ratio=}')
    # calculate for other time steps if ratio > 3 => (last - first)/ (middle/2) > 3  => (last - first) / (middle/2)   
    if ratio > starting_fovs_n :
        time_steps = np.linspace(0,int_time,3)
        for time_step_ind,time_step in enumerate(time_steps):
            # chceck if we already calculated this step
            if not time_step in shapes :
                # time = row[time_string]
                # et = spiceypy.utc2et(str(time))
                et = row[time_string]+time_step
                # print(f'{time_step_ind=} {time_step=:5.0f} {time=:26s} {et=} ',)
                internal_points_df = pd.DataFrame.from_dict([geometry_intersect( method, target, et, pid, target_frame, aberration_correction, observer, sensor_frame, radii, tarid, v) for v in fov_points] )
                internal_points_df['time_step'] = time_step
                internal_points_ls.append(internal_points_df)
                shape = shapely.geometry.Polygon( internal_points_df.apply(lambda x : (x['tarlon'],x['tarlat']),axis=1).to_list() )    
                shapes[time_step] = [pid,shape]+internal_points_df[calculated_cols].std().to_list()+internal_points_df[calculated_cols].mean().to_list()
            # else:
                # print(f'{time_step_ind=} {time_step=:5.0f} ALREADY DONE - SKIPPING ',)

    tmp = pd.DataFrame(index=list(shapes.keys()),
                        data=list(shapes.values()),
                        columns=['pid','fov']+[f'std_{c}' for c in calculated_cols]+[f'mean_{c}' for c in calculated_cols]
                      )

    tmp.index.name = 'time_steps'

    return tmp.reset_index(drop=0).sort_values(['time_steps', 'pid']).set_index(['pid'])

