#Written by Arturo S. Leon, Last updated 04/03/2021

#Retrieve Bias-Corrected Global Forecast System (GFS) and convert to DSS for its use in HEC-HMS
import os, subprocess
from siphon.catalog import TDSCatalog
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pytz

from netCDF4 import num2date
#from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from pyproj import Proj
import pyproj
import rasterio
from rasterio.transform import from_origin
#from mpl_toolkits import mplot3d
import shutil, sys
import pandas as pd

# Define bounding box for Cypress Creek Watershed
North = 30.25
South = 29.75
West = -96
East = -95.25

#Bigger extent encompassing watershed and few states around it.
##North = 41
##South = 22
##West = -107
##East = -79

LeadTime = 5 # in days. This time will be added to the current time
res = 1000. #resampled resolution in meter UTM (1000). This needs to be of
            #the same resolution as your .MOD grid in HEC-HMS 

##################################################
#local = pytz.timezone ("America/Chicago") #Enter here the time zone to which the watershed belongs
#Uncomment one of the two lines below
#Initial_time_precipit = 'now' #use this for forecast data
#Initial_time_precipit = datetime(2018, 9, 15, 13, 39,0) #Enter initial local time of watershed in the following format
# Year, Month, Day, hour, minute, and second (use this for historical data)

##try: Initial_time_precipit
##except:
##    print("The variable Initial time was not defined.")
##    print("Make sure that one of the **Initial_time_precipit** lines is uncommented")
##    raise

##################################################
##if Initial_time_precipit == 'now':
##    initial_Date = datetime.utcnow()
##    print('Initial date for precipitation (UTM time) = ',initial_Date)
##    print('Forecast mode')
##else:
##    local_dt = local.localize(Initial_time_precipit, is_dst=None)
##    utc_dt = local_dt.astimezone(pytz.utc)
##    initial_Date = utc_dt
##    print('Initial date for precipitation (UTM time) = ',initial_Date)
##    print('historical data mode')


# Define directories
save_dir = 'Forecast_GFS/'
dir_openBATfile = 'util/ASCIIToDSS.bat'
dir_asc2dssGrid = 'util/asc2dssGrid.exe'
home_folder = os.getcwd()

if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

filesASC = os.listdir(save_dir)

for item in filesASC:
     if item.endswith(".asc"):
        os.remove(os.path.join(save_dir, item))
        if item.endswith(".prj"):
           os.remove(os.path.join(save_dir, item))


# use pyproj to convert lats lons to UTM
GRIB2_proj = Proj("+proj=longlat +ellps=WGS84 +pm=-360 +datum=WGS84 +no_defs") #put pm=-360 to indicate GRIB2 projection

#Projected = Proj(init="EPSG:26914") #26914 Cypress UTM 14N
Projected = Proj("EPSG:26914") #26914 Cypress UTM 14N
#Parameter = 'Total_precipitation_surface_Mixed_intervals_Accumulation'
Parameter = 'Precipitable_water_entire_atmosphere_single_layer'
   

# Get the dataset handle
##top_cat = TDSCatalog('http://thredds.ucar.edu/thredds/catalog.xml')
##models_cat = top_cat.catalog_refs[0].follow()
##gfs_cat = models_cat.catalog_refs['GFS Quarter Degree Forecast'].follow()
##ncss = gfs_cat.latest.subset()

top_cat = TDSCatalog('http://thredds.ucar.edu/thredds/catalog.xml')
ref = top_cat.catalog_refs['Forecast Model Data']
#ref = top_cat.catalog_refs['Radar Data']
print('ref',ref)
models_cat = ref.follow()
gfs_cat = models_cat.catalog_refs['GFS Quarter Degree Forecast'].follow()
#gfs_cat = models_cat.catalog_refs['NEXRAD Level III Radar'].follow()
ncss = gfs_cat.latest.subset()

#now = datetime.utcnow()
initial_Date = datetime.utcnow()

# Download a subset in box using NCSS
queryBox = ncss.query().lonlat_box(east=East, west=West, south=South, north=North)
queryBox.time_range(initial_Date, initial_Date + timedelta(days=LeadTime)).accept('netcdf4')
queryBox.variables(Parameter)
data = ncss.get_data(queryBox)

Precipitation_var = data.variables[Parameter] # pull time variable out of the coordinates attribute on Precipitation
time_name = Precipitation_var.coordinates.split()[0]
time_var = data.variables[time_name]

lat_var = data.variables['lat']
lon_var = data.variables['lon']
lat_vals = lat_var[:].squeeze()
lon_vals = lon_var[:].squeeze()

#convert coordinate
x,y = np.meshgrid(lon_vals, lat_vals)
UTMx, UTMy = pyproj.transform(GRIB2_proj, Projected, x.flatten(), y.flatten())
UTMx_grid = np.reshape(UTMx,(x.shape))
UTMy_grid = np.reshape(UTMy,(y.shape))
xmin, xmax, ymin, ymax = [UTMx.min(), UTMx.max(), UTMy.min(), UTMy.max()]
UTMx_res = np.arange(xmin, xmax, res)
UTMy_res = np.arange(ymin, ymax, res)
grid_x_res, grid_y_res = np.meshgrid(UTMx_res, UTMy_res)
nrows_res, ncols_res = np.shape(grid_x_res)

#Minimum value of precipitation matrix
Min_precip = np.array(Precipitation_var)
Min_single_precip = np.nanmin(Min_precip)

for timeId in range(0, np.shape(time_var)[0], 2):
    
    # Get the actual data values and remove any size 1 dimensions
    Precipitation_vals = Precipitation_var[timeId].squeeze()
    
    # Convert the number of hours since the reference time to an actual date, check the array index
    time_val = num2date(time_var[timeId].squeeze(), time_var.units)    
    
    #convert date to string
    filenameASC = save_dir + 'gfs.0p25.' + time_val.strftime("%Y%m%d%H") + '.asc'
       
    #resample raster
    Precipitation_res = interpolate.griddata((UTMx_grid.flatten(), UTMy_grid.flatten()), Precipitation_vals.flatten(), (grid_x_res, grid_y_res) , method='cubic')
    #Precipitation_res[np.isnan(Precipitation_res)]=-9999
    Precipitation_res[np.isnan(Precipitation_res)]=-9999
    Precip_plot = Precipitation_res
    #save resampled raster to ascii
    transform = from_origin(grid_x_res.min(), grid_y_res.max(), res, res)
    new_dataset = rasterio.open(filenameASC, 'w', nodata = -9999, driver='AAIGrid', decimal_precision=1,
                                height = nrows_res, width = ncols_res,
                                count=1, dtype=str(Precipitation_res.dtype),
                                transform=transform,
                                crs=Projected.srs)
    new_dataset.write(Precipitation_res, 1)
    new_dataset.close()

#plot
def Plot_Parameter(Precip_plot, grid_x_res, grid_y_res,Min_single_precip,time_val):       
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib import colors as mcolors
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import cartopy 
    import cartopy.crs as ccrs
    from cartopy.io import shapereader
    from cartopy.io.shapereader import Reader
    from cartopy.feature import ShapelyFeature
    import cartopy.io.shapereader as shpreader
    
    import numpy as np
    import pandas as pd
    import shapefile as shp    
    from cartopy.mpl.geoaxes import GeoAxes
    ZoneNo = "14" #Manually input or from other sources
    myProj = Proj("+proj=utm +zone="+\
    ZoneNo+", +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    Lon2, Lat2 = myProj(grid_x_res, grid_y_res,inverse=True)

# draw filled contours.
    clevs = [0, 5, 7.5, 10, 15, 20, 30, 40,
         50, 70, 100, 150, 200, 250, 300, 400, 500, 600, 750]
    cmap_data = [(1.0, 1.0, 1.0),
             (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
             (0.0, 1.0, 1.0),
             (0.0, 0.8784313797950745, 0.501960813999176),
             (0.0, 0.7529411911964417, 0.0),
             (0.501960813999176, 0.8784313797950745, 0.0),
             (1.0, 1.0, 0.0),
             (1.0, 0.6274510025978088, 0.0),
             (1.0, 0.0, 0.0),
             (1.0, 0.125490203499794, 0.501960813999176),
             (0.9411764740943909, 0.250980406999588, 1.0),
             (0.501960813999176, 0.125490203499794, 1.0),
             (0.250980406999588, 0.250980406999588, 1.0),
             (0.125490203499794, 0.125490203499794, 0.501960813999176),
             (0.125490203499794, 0.125490203499794, 0.125490203499794),
             (0.501960813999176, 0.501960813999176, 0.501960813999176),
             (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
             (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
             (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
             (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
             (0.4000000059604645, 0.20000000298023224, 0.0)]
    cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
    norm = mcolors.BoundaryNorm(clevs, cmap.N)

    useproj = ccrs.PlateCarree()
    f, ax = plt.subplots(1)
    ax = plt.axes(projection=useproj)    
    extent = [West, East, South, North]
    ax.set_extent(extent)
    
    Precip_plot[Precip_plot < 0]=Min_single_precip
##    plt.contourf(Lon2, Lat2, Precip_plot, 20, cmap='BuGn',
##             transform=ccrs.PlateCarree())
    plt.contourf(Lon2, Lat2, Precip_plot, 20, cmap=cmap, norm=norm,
             transform=ccrs.PlateCarree())   
    plt.colorbar()
    
    # Make a title with the time value
    ax.set_title(u'6-hour Cumul. Precip. Forecast (mm) at UTC time {0:%d %B %Y %H:%M}'.format(time_val), fontsize=20)
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    
    ax.add_feature(cfeature.STATES)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.legend()
    gl = ax.gridlines(crs=useproj, draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    #gl.xlines = False
    #gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 14, 'color': 'gray','weight': 'bold'}
    gl.ylabel_style = {'size': 14, 'color': 'gray','weight': 'bold'}
    #plt.show()
    fname = "Maps_shapefiles\Subbasin356.shp"
    assert os.path.exists(fname), "Input file does not exist."
    adm1_shapes = list(shpreader.Reader(fname).geometries())
    #Make sure shape file is in WGS 84 coordinates (latitude/Longitude)
    #Normally ArcMap creates shapefiles in UTM. You can convert the shapefile
    #to Lat/Long using ArcMAP or other convertors such as
    #https://mygeodata.cloud/converter/shp-to-latlong
    ax.add_geometries(adm1_shapes, useproj,
                  edgecolor='black', facecolor='red', alpha=0.90)
    #Show the plot during code execution
    #plt.show()

    
    #Save rather than show
    plt.savefig('Forecast_GFS\precip_plot.pdf')  

# Plot Results
Plot_Parameter(Precip_plot, grid_x_res, grid_y_res,Min_single_precip,time_val)

# convert to DSS
import re
import glob 
import platform
from datetime import datetime, timedelta, date

# Local variables:
# my_path = os.path.abspath(os.path.dirname(__file__))
# inDir = os.path.join(my_path, "..\Forecast_GFS")
# openBATfile = my_path  + "\ASCIIToDSS.bat"

# Local variables:
#inDir = 'Precipitation_Mongolia/ASC_2000_Dis/'
openBATfile = save_dir + 'ASCIIToDSS.bat'

# for generating asc2dss strings
#text_0 = set PATH = C:\temp\Gridded_Workshop_FINAL_MATLAB\utilitaries\hecexe\;%PATH%
text_0 = 'set PATH = ' + home_folder + '\\util\;%PATH%'
text_1 = "asc2dssGrid INPUT="
text_2 = " DSS=GFS.dss PATH=/UTM14/Cypress/Precip/"
text_3 = "/"
text_4 = "/PROJECTED/ GRIDTYPE=UTM ZONE=14N DUNITS=mm DTYPE=PER-CUM" 


#Copying BAT file and asc2dssGrid.dss
shutil.copy(dir_openBATfile, save_dir) # Copy BAT file
shutil.copy(dir_asc2dssGrid, save_dir) # Copy asc2dssGrid file

# Open and Write to file1  
file1 = open(openBATfile,"w")

# Write the first line of the openBATfile (set PATH = C:\...\util\;%PATH%)
file1.write(text_0 + "\n")

# date = [date for file in glob.glob(inDir + '*.asc') for date in re.findall("(\d{10})", file)]
date = [date for file in glob.glob(save_dir + '*.asc') for date in re.findall("(\d{10})", file)]
date.sort(key = lambda date: datetime.strptime(date, "%Y%m%d%H")) 
delt = datetime.strptime(date[1], "%Y%m%d%H") - datetime.strptime(date[0], "%Y%m%d%H")
# AscFiles = glob.glob(inDir + '*.asc')
AscFiles = glob.glob(save_dir + '*.asc')

for index in range(0, len(date)):

    namefile = 'gfs.0p25.' + date[index] + '.asc'
    objDate = datetime.strptime(date[index], '%Y%m%d%H')

    if (objDate.hour==0): 
        objDate0 = objDate - delt #timedelta(hours=1)
        objDate = objDate - delt #timedelta(days=1)
        converted_date = datetime.strftime(objDate,'%d%b%Y:2400') #HEC-DSS only accept 24:00 format!
    else:
        converted_date = datetime.strftime(objDate,'%d%b%Y:%H%M')
        objDate0 = objDate - delt#timedelta(hours=1)
    
    converted_date0 = datetime.strftime(objDate0,'%d%b%Y:%H%M')

    full_text = text_1 + namefile + text_2 + converted_date0 + text_3 + converted_date + text_4
    file1.write(full_text + "\n")   	 
file1.close() #to close file

