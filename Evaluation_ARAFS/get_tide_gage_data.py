import pandas as pd
import requests
import matplotlib.pyplot as plt

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Charleston, OR 
station_id = "9432780"

url = f"https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{station_id}.json"

response = requests.get(url)
data = response.json()["stations"][0]

lat_station = data["lat"]
lon_station = data["lng"]
name_station = data["name"]

# Setup NOAA CO-OPS API URL
url = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"

params = {
    "begin_date": "20230101",
    "end_date": "20230228",
    "station": station_id,
    "product": "hourly_height",  # Options: hourly_height, water_level (6-min), high_low
    "datum": "MSL",  # Mean Lower Low Water (Standard tidal datum)
    "units": "english",  # 'english' (feet) or 'metric' (meters)
    "time_zone": "gmt",  # 'gmt', 'lst' (local standard), or 'lst_ldt' (local day light)
    "format": "json",
    "application": "PythonScript",
}

response = requests.get(url, params=params)
data = response.json()

if "data" in data:
    df = pd.DataFrame(data["data"])
else:
    print("Error or no data returned:", data)

time_obs = pd.to_datetime(df['t'])
msl_obs = df['v'].astype(float)*0.3048

####################################################
plt.figure(figsize=(10,5))
plt.plot(time_obs,msl_obs,'.-')
plt.ylabel('Mean Sea level (m)')
plt.title(name_station,fontsize=18)

####################################################
plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.axis('scaled')
ax.plot(lon_station,lat_station,'*r',markersize=5)
# Add borders and coastlines
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
plt.xlim([-130,-110])
plt.ylim([10,60])
plt.title(name_station,fontsize=18)
