import numpy as np
import pyproj
import pandas as pd

# Generate transformer from latlon to polar stereographic
crs0 = pyproj.CRS('WGS84')
crs1 = pyproj.CRS('epsg:3413')
transformer_xy = pyproj.Transformer.from_crs(crs0, crs_to=crs1, always_xy=True)
transformer_ll = pyproj.Transformer.from_crs(crs1, crs_to=crs0, always_xy=True)


# Each of these scenes has distinct ice floes
cases = {
    'greenland_sea': {'center_lat': 77.1718,
                      'center_lon': -13.8062,
                      'date': '2018-07-12'},
    'barents_kara_seas': {'center_lat': 78.9808,
                       'center_lon': 48.8454,
                         'date': '2025-05-07'},
    'laptev_sea': {'center_lat': 76.5526,
                   'center_lon': 123.4811,
                   'date': '2023-06-13'},
    'sea_of_okhostk': {'center_lat': 57.0740,
                        'center_lon':  142.4997,
                        'date': '2018-04-12'},
    'east_siberian_sea': {'center_lat': 72.2032,
                          'center_lon': 172.7133,
                          'date': '2019-06-12'},
    'bering_chukchi_seas': {'center_lat': 71.2663,
                            'center_lon': -161.3376,
                            'date': '2021-03-14'},
    'beaufort_sea': {'center_lat': 72.4100,
                     'center_lon': -136.2380,
                     'date': '2007-06-14'},
    'hudson_bay': {'center_lat': 60.7905,
                   'center_lon': -84.6631,
                   'date': '2006-05-13'},
    'baffin_bay': {'center_lat': 68.8493,
                   'center_lon': -63.4067,
                   'date': '2006-07-04'}
}

columns = ['case_name', 'location', 'center_lat', 'center_lon', 'top_left_lat', 'top_left_lon',
'lower_right_lat', 'lower_right_lon', 'left_x', 'right_x', 'lower_y', 'top_y', 'startdate', 'enddate']

dx_list = [50, 100, 200, 50, 100, 200, 50, 100, 200]
case_list = []
for region, dx in zip(cases, dx_list):
    date = cases[region]['date']
    date2 = (pd.to_datetime(date) + pd.to_timedelta(1, 'day')).strftime('%Y-%m-%d')
    lat0, lon0 = cases[region]['center_lat'], cases[region]['center_lon']
    x0, y0 = transformer_xy.transform(lon0, lat0)
    left_x = x0 - dx/2*1e3
    right_x = x0 + dx/2*1e3
    lower_y = y0 - dx/2*1e3
    upper_y = y0 + dx/2*1e3
    top_left_lon, top_left_lat = transformer_ll.transform(left_x, upper_y)
    lower_right_lon, lower_right_lat = transformer_ll.transform(right_x, lower_y)
    
    case_list.append(
        [region + '-' + str(dx) + 'km' + '-' + date.replace('-', ''),
         region,
         lat0,
         lon0,
         np.round(top_left_lat, 5),
         np.round(top_left_lon, 5),
         np.round(lower_right_lat, 5),
         np.round(lower_right_lon, 5),
         np.round(left_x),
         np.round(right_x),
         np.round(lower_y),
         np.round(upper_y),
         date,
         date2]) # extra day for tracking
case_descriptions = pd.DataFrame(case_list, columns=columns)
case_descriptions.to_csv('test_case_descriptions.csv', index=False)