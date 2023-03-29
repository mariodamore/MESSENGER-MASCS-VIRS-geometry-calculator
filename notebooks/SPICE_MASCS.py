# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Import 

# %%
import matplotlib.pyplot as plt
import pathlib
import rich
import numpy as np
import pandas as pd

opts = {
    'display.max_columns': None,
    'display.max_colwidth': 50,
    'display.expand_frame_repr': False,  # Don't wrap to multiple pages
    'display.max_rows': 999, # maximum number of rows pandas should output 
    'display.max_seq_items': 50,         # Max length of printed sequence
    'display.precision': 4, # Floating point output precision
    'display.show_dimensions': True,
    'display.width':200, # Width of the display in characters. in a terminal this can be set to None, not in notebook or qtconsole.
    }

# set pandas options
[pd.set_option(n,o) for n,o in opts.items()]

# useful for displaying html, images and so on
from IPython.display import display

import geopandas as gpd

# %% [markdown]
# ## MASCS Data

# %% [markdown]
# ### direct osgeo/gdal read

# %%
from osgeo import ogr
import geopandas as gpd
import pandas as pd
import shapely
# useful for displaying html, images and so on
from IPython.display import display

# %%
file = "../data/raw/test/mess-e_v_h-mascs-3-virs-cdr-caldata-v1/messmas_2101/data/ddr/orb/virs/mascs20110812/vis/virsvd_orb_11224_214005.lbl"

# %%
df = gpd.read_file(file,engine='fiona').drop(columns=['SPARE_1', 'SPARE_2', 'SPARE_3','SPARE_4', 'SPARE_5',]) # default engine, no warnings

# %%
# convert SPECTRUM_UTC_TIME to python datetime object 
from datetime import datetime
df['SPECTRUM_UTC_TIME'] = pd.to_datetime(df['SPECTRUM_UTC_TIME'].map(lambda x : datetime.strptime(x.strip(), '%y%jT%H:%M:%S')))

# %%
df.iloc[[0,-1]].T

# %%
# drop the index in a column
df.index.name = 'pid' # primary ID
df = df.reset_index(drop=0)

# %%
inDataSource = ogr.Open(file)
table = inDataSource.GetLayer()

# %%
TARGET_LONGITUDE_SET = np.array([t.GetField('TARGET_LONGITUDE_SET') for t in table])
TARGET_LATITUDE_SET  = np.array([t.GetField('TARGET_LATITUDE_SET') for t in table])
stacked_latlon = np.stack((TARGET_LONGITUDE_SET,TARGET_LATITUDE_SET))
print(f'{TARGET_LATITUDE_SET.shape=}')
print(f'{TARGET_LONGITUDE_SET.shape=}')
print(f'{stacked_latlon.shape=}')

# %%
# this select only the polygons c1,c3,c2,c4 in this order, see Fig.4 in VIRS_CDRSIS.PDF
corners_index = [1,3,2,4]
coordinates = np.moveaxis(stacked_latlon[:,:,corners_index].T,1,0)
polygon_list = [shapely.geometry.Polygon(p)  for p in  coordinates]

# %%
corners_index = 0
points_list = [shapely.geometry.Point(p)  for p in  stacked_latlon[:,:,corners_index].T]

# %% [markdown]
# ## SPICE kernels use
#

# %%
import spiceypy

# %%
kernels_basepath = pathlib.Path().resolve() / '../data/raw/spice_kernels/NASA/messsp_1000/data'
# kernels_basepath = pathlib.Path('/home/kidpixo/Documents/work/MESSENGER/MASCS_mda_dlr_processed/spice_kernels_messenger/data')

print(f'{kernels_basepath=}')
print(f'{kernels_basepath.exists()=}')


# %%
# symlink SPICE path to /tmp/spice_kernels_messenger > SPICE Paths are limited to 255 chars!!!! 
original_spice_dir = pathlib.Path(kernels_basepath)

spice_dir = pathlib.Path('/tmp/spice_kernels_messenger')

if not spice_dir.exists():
    spice_dir.symlink_to(original_spice_dir)

# %%
# symlink SPICE path to /tmp/spice_kernels_messenger > SPICE Paths are limited to 255 chars!!!! 
esa_dsk_path = pathlib.Path('').resolve() / '../data/raw/spice_kernels/ESA/kernels/dsk'

new_dsk_path = original_spice_dir / 'dsk' 

if not new_dsk_path.exists():
    new_dsk_path.symlink_to(esa_dsk_path)

# %%
spiceypy.kclear()# Clear the Spice kernels pool

# %%
metakernels_paths = list(pathlib.Path(kernels_basepath / '../..').glob('*tm'))
# metakernels_paths = list(pathlib.Path(kernels_basepath / '../extras/mk_mod').glob('*tm'))
metakernels_paths_dict = {p.stem.split('_')[1]:p for p in metakernels_paths}
metakernels_paths_dict

# %%
print(f'Data year of aquisition : {df["SPECTRUM_UTC_TIME"].dt.year.unique()[0]}')

# %%
# # subset metakernel
# spiceypy.furnsh(str(kernels_basepath / 'msgr_2011_v10_110812_110812.tm'))  # load the kernel in the pool

# mission metakernel
spiceypy.furnsh(str(metakernels_paths_dict['2011']))  # load the kernels

# %%
sensor = 'MSGR_MASCS_VIRS'
sensor_id = spiceypy.bodn2c(sensor)
target = 'MERCURY'
tarid = spiceypy.bodn2c(target)

n, radii = spiceypy.bodvrd(target, 'RADII', 3)
re = radii[0]
rp = radii[2]
radii_relative_difference = (re - rp) / re # radii relative difference

rich.print({
    'sensor':sensor,
    'sensor_id':sensor_id,
    'target':target,
    'tarid':tarid,
    'n':n,
    'radii': radii,
    'radius equ.' : re,
    'radius pol.' : rp,
    'radii_relative_difference': radii_relative_difference
})

# %% [markdown]
# If the value of `shape` is "CIRCLE" the field of view of
# the instrument is a circular cone centered on the
# boresight vector. 
#
# The vertex of the cone is at the instrument focal point.
#
# A single vector will be returned in `bounds`.
#
# This vector will be parallel to a ray that lies in the cone that makes up the boundary
# of the field of view.
#

# %%
(shape, sensor_frame, bsight, vectors, bounds) = spiceypy.getfov(sensor_id,20)

# %%
rich.print({
    'shape':shape,
    'sensor_frame':sensor_frame,
    'bsight':bsight,
    'vectors':vectors,
    'bounds':bounds,
})

# %%
bound_norm = bounds[0]/spiceypy.spiceypy.unorm(bounds[0])[1]

# %%
# rotate bounds vector

# angular step in deg
steps = 12*10
angle_step = 360/steps
print(f'{steps=}')
print(f'{angle_step=}')
angle_step = 360//steps

# %%
# # print FOV steps
# for a in range(0,steps):
#     print(f'{a:2d} {angle_step*(a):3d}',spiceypy.vrotv([1,0,0],bsight, np.deg2rad(angle_step*(a))) )  

# %%
plt.figure(figsize=[8,8])

fov_points = np.asarray([spiceypy.vrotv(bound_norm, bsight, np.deg2rad(angle_step*(a))) for a in range(0,steps)])

plt.scatter(fov_points[:,0],fov_points[:,1],s=30)
plt.scatter(0,0,s=30)

plt.xlabel('X')
plt.ylabel('Y')
plt.title(f'MASCS VIS Istantaneous Field of View (iFOV) approximated as {steps} points')

# %%
observer = 'MESSENGER'
target_frame = 'IAU_MERCURY'

# %%
et = spiceypy.utc2et('2011-08-12 21:42:15')       # Convert the UTC datetime in Ephemris Time (ET)
rich.print(
{'computation method':'ELLIPSOID',
'target':target,
'et':et,
'target body frame':target_frame,
'aberration correction':'LT+S',
'observing body':observer,
'sensor_frame':sensor_frame,
})

# %%
# method = 'ELLIPSOID'
method = 'DSK/UNPRIORITIZED/SURFACES = \"199\"'
aberration_correction='LT+S'

# %%
# rich.print(geometry_intersect(method,target,et,target_frame, aberration_correction,observer,sensor_frame,bsight) )

# %%
from msgmascsgeo import spicefuncs

# %%
# calculate the Ephemeris Time (ET)
spacecraft_id = spiceypy.bodn2c('MESSENGER')
df['ET'] = df.apply(lambda x : spiceypy.scs2e(spacecraft_id,str(x.SPECTRUM_MET+x.SPECTRUM_SUBSECONDS/1000)) ,axis=1)

# %%
gdf_pds_poly  = gpd.GeoDataFrame(data=df.copy() ,geometry=polygon_list)
gdf_pds_point = gpd.GeoDataFrame(data=df.copy(),geometry= points_list)

# %%
gdf_pds_poly.plot()

# %%
ax = gdf_pds_poly.head(2).plot(figsize=[20,8],facecolor='none',column='SC_TIME',cmap=plt.cm.Spectral_r)
gdf_pds_point.head(2).plot(ax=ax)

# %%
(method,
target,
spiceypy.utc2et(str(df.loc[200,'SPECTRUM_UTC_TIME'])),
df.loc[200,'pid'],
target_frame,
aberration_correction,
observer,
sensor_frame,
radii,
tarid,
bsight)

# %%
# 1 fov point for fixed times 
spicefuncs.geometry_intersect(
                    method,
                    target,
                    spiceypy.utc2et(str(df.loc[200,'SPECTRUM_UTC_TIME'])),
                    df.loc[200,'pid'],
                    target_frame,
                    aberration_correction,
                    observer,
                    sensor_frame,
                    radii,
                    tarid,
                    bsight)

# %%
# fov points for fixed times
points_df = pd.DataFrame.from_dict([
    spicefuncs.geometry_intersect(
                    method,
                    target,
                    spiceypy.utc2et(str(df.loc[200,'SPECTRUM_UTC_TIME'])),
                    df.loc[200,'pid'],
                    target_frame,
                    aberration_correction,
                    observer,
                    sensor_frame,
                    radii,
                    tarid,
                    bsight) for v in fov_points] )
# points_df.T

# %%
from tqdm.notebook import trange, tqdm
from tqdm.gui import tqdm as tqdm_gui

tqdm.pandas()  # can use tqdm_gui, optional kwargs, etc

# %%
# calculate for all measurements in orbit at spectrum_utc_time (vector)
to_fovs_df = gdf_pds_poly.iloc[5:40]
print(to_fovs_df.shape)

# %%
to_fovs_df.plot(facecolor='none')

# %%
import functools

# %%
spicefuncs.row_to_shape(  to_fovs_df.iloc[0],
                          starting_fovs_n=5,
                          radii=radii,
                          tarid=tarid,
                          fov_points=fov_points, 
                          method=method,
                          target=target,
                          target_frame=target_frame,
                          aberration_correction=aberration_correction,
                          observer=observer,
                          sensor_frame=sensor_frame)

# %%
ifovs_df = pd.concat(to_fovs_df.reset_index().progress_apply(functools.partial(spicefuncs.row_to_shape,
                                                                              starting_fovs_n=5,
                                                                              radii=radii,
                                                                              tarid=tarid,
                                                                              fov_points=fov_points, 
                                                                              method=method,
                                                                              target=target,
                                                                              target_frame=target_frame,
                                                                              aberration_correction=aberration_correction,
                                                                              observer=observer,
                                                                              sensor_frame=sensor_frame,)
                                                            ,axis=1).values)

# %%
# convert iFOVS to Geopandas.GeoDataFrame
geofovs_df = gpd.GeoDataFrame(ifovs_df.reset_index(drop=0), geometry='fov')
# iFOVS unary_union.convex_hull : size == initial size, one geometry per measurement
geounion_gdf = gpd.GeoDataFrame(
    geofovs_df.drop(columns=['fov']).groupby('pid').agg(np.median),             # group per PID and average
    geometry=[v.unary_union.convex_hull for k,v in geofovs_df.groupby('pid')]# assemble union geometry
).reset_index(drop=0)


# %%
geofovs_df.columns

# %%
fig, axs = plt.subplots(ncols=2,nrows=2,figsize=[15,10])
axs = axs.flatten()

gdf_pds_poly.set_index('pid').loc[geounion_gdf.pid].plot(ax=axs[0],column='INCIDENCE_ANGLE',edgecolor='black',cmap=plt.cm.Spectral_r,alpha=0.75)
axs[0].set_title('PDS, color = incidence angle');

gdf_pds_poly.set_index('pid').loc[geounion_gdf.pid]['INCIDENCE_ANGLE'].plot(marker='.',ax=axs[2])
axs[2].set_title('PDS, incidence angle');


geofovs_df.plot(ax=axs[1],column='mean_incdnc',edgecolor='black',cmap=plt.cm.Spectral_r,alpha=0.75)
axs[1].set_title('NEW-Union, color = incidence angle');

geofovs_df['mean_incdnc'].plot(marker='.',ax=axs[3])
axs[3].set_title('NEW-Union, incidence angle');

_ = [a.set_xlabel('Sequential Meas. Index') for a in axs[-2:]]
_ = [(a.set_xlabel('Longitude'),a.set_ylabel('Latitude')) for a in axs[:2]]


# %% [markdown]
# ## Save output to csv , geojson, shapefiles etc
#
#
# Set `output` to a path and then run 

# %%
output = pathlib.Path('/tmp/')

# %%
gdf_pds_poly.to_csv( output / 'NEW-iFOV.csv.zip')

geounion_gdf.to_wkt().rename(columns={'geometry':'fov'}).to_csv( output / 'NEW-FOV-union-convex.csv.zip')


# %%
#  if yu wish, change the driver to shapefile 
gdf_pds_poly.to_file(output / 'NEW-iFOV.geojson', driver="GeoJSON")
geounion_gdf.to_file(output / 'NEW-FOV-union-convex.geojson', driver="GeoJSON")
