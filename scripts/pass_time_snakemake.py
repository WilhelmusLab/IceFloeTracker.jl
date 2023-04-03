# Code developed by Simon Hatcher (2022)
# Adapted for use in a snakemake pipeline for IceFloeTracker by Timothy Divoll (2023)

# ACTION REQUIRED FROM YOU:
# 1. Update the following parameters in the `snakemake-config.yaml`:

# startdate = YYYY-MM-DD
# enddate = YYYY-MM-DD
#
# centroid-x = DD.DDDD
# centoid-y = DD.DDDD

# Centroid is the approximate point in the middle of your bounding box area of interest
# Your www.space-track.org credentials (https://www.space-track.org/auth/createAccount for free account) need to be set as environment variables in .bash_profile or .zshrc or add to ENV VARS in Windows
# NOTE: PASSWORD FIELD IS NOT SECURE. DO NOT USE USER/PASSWORD DIRECTLY IN CONFIG FILES.

# 2. A stable internet connection is also required.

# Package imports.
import requests
import json
import datetime
from skyfield.api import wgs84, load, EarthSatellite
import numpy as np
import csv
import math
import yaml
import os


# URLs for space track login.
uriBase = "https://www.space-track.org"
requestLogin = "/ajaxauth/login"

# Define error.
class MyError(Exception):
    def __init___(self,args):
        Exception.__init__(self,"my exception was raised with arguments {0}".format(args))
        self.args = args

def main():

    # Read in username, password, start date, end dates from configuration file.
    
    def load_conf_file(config_file):
        with open(config_file, "r") as f:
            config = yaml.safe_load(f)
            startdate = config["startdate"]
            enddate = config["enddate"]
            csvoutpath = config["soit-outdir"]
            centroidx = config["centroid-x"]
            centroidy = config["centroid-y"]

        return startdate, enddate, csvoutpath, centroidx, centroidy
    
    startdate, enddate, csvoutpath, centroidx, centroidy = load_conf_file("hpc/snakemake-config.yaml")
    configUsr = os.environ["SPACEUSER"]
    configPwd = os.environ["SPACEPSWD"]
    siteCred = {'identity': configUsr, 'password': configPwd}
    end_date = datetime.datetime.strptime(enddate, "%Y-%m-%d").strftime("%m-%d-%Y").split('-')
    start_date = datetime.datetime.strptime(startdate, "%Y-%m-%d").strftime("%m-%d-%Y").split('-')
    lat = int(centroidx)
    long = int(centroidy)
    print(csvoutpath)
    print(f"Timeframe starts on {start_date}, and ends on {end_date}")
    print(centroidx, centroidy)

    end_date = getNextDay(end_date)

    # Get TLEs from space track.
    with requests.Session() as session:

        # Log in with username and password.
        resp = session.post(uriBase + requestLogin, data=siteCred)
        if resp.status_code != 200:
            raise MyError(resp, "POST fail on login. Your username/password may be incorrect.")

        # Retrieve Aqua TLEs from space track.
        resp = session.get(
            f'https://www.space-track.org/basicspacedata/query/class/gp_history/NORAD_CAT_ID/27424/orderby/TLE_LINE1%20ASC/EPOCH/{start_date[2]}-{start_date[0]}-{start_date[1]}--{end_date[2]}-{end_date[0]}-{end_date[1]}/format/json')
        if resp.status_code != 200:
            print(resp)
            raise MyError(resp, "GET fail on request")

        # Turn JSON into Python dict.
        aquaData = json.loads(resp.text)

        # Retrieve Terra TLEs from space track.
        resp = session.get(
            f'https://www.space-track.org/basicspacedata/query/class/gp_history/NORAD_CAT_ID/25994/orderby/TLE_LINE1%20ASC/EPOCH/{start_date[2]}-{start_date[0]}-{start_date[1]}--{end_date[2]}-{end_date[0]}-{end_date[1]}/format/json')
        if resp.status_code != 200:
            print(resp)
            raise MyError(resp, "GET fail on request")

        # Turn JSON into Python dict.
        terraData = json.loads(resp.text)

        # No more requests.
        session.close()

    # Define index variables for each satellite.
    aqua_i = 0
    terra_i = 0

    # Load in orbital mechanics tool timescale.
    ts = load.timescale()

    # Specify area of interest.
    aoi = wgs84.latlon(lat, long)

    # Define today and tomorrow.
    today = start_date
    tomorrow = getNextDay(start_date)

    # Define 2D array of values to be added to results CSV.
    rows = []

    # Loop through each day until the end date of interest is reached.
    while True:

        # Get UTC time values of the start of today and the start of tomorrow.
        # Passes between these times are considered.
        t0 = ts.utc(int(today[2]), int(today[0]), int(today[1]))
        t1 = ts.utc(int(tomorrow[2]), int(tomorrow[0]), int(tomorrow[1]))

        # Load in aqua TLE from space track data.
        aqua_tleline1 = aquaData[aqua_i]['TLE_LINE1']
        aqua_tleline2 = aquaData[aqua_i]['TLE_LINE2']

        # Create new aqua satellite with this TLE.
        aqua = EarthSatellite(aqua_tleline1, aqua_tleline2, 'AQUA', ts)

        # Load in terra TLE from space track data.
        terra_tleline1 = terraData[terra_i]['TLE_LINE1']
        terra_tleline2 = terraData[terra_i]['TLE_LINE2']

        # Create new terra satellite with this TLE.
        terra = EarthSatellite(terra_tleline1, terra_tleline2, 'TERRA', ts)

        # Update aqua index. If the epoch of the current TLE is more than a day away from the area of interest,
        # go to a newer TLE. This keeps the TLE epoch close to the date of the pass, decreasing error in the orbital
        # mechanics algorithm.
        while ((t0 - aqua.epoch) > 1):
            aqua_i+=1
            aqua_tleline1 = aquaData[aqua_i]['TLE_LINE1']
            aqua_tleline2 = aquaData[aqua_i]['TLE_LINE2']
            aqua = EarthSatellite(aqua_tleline1, aqua_tleline2, 'AQUA', ts)

        # Update terra index.
        while ((t0 - terra.epoch) > 1):
            terra_i+=1
            terra_tleline1 = terraData[terra_i]['TLE_LINE1']
            terra_tleline2 = terraData[terra_i]['TLE_LINE2']
            terra = EarthSatellite(terra_tleline1, terra_tleline2, 'TERRA', ts)

        # Find times in the time interval of interest where aqua's position in the sky passes over 30º for an observer
        # in the area of interest, i.e. when it's passing by relatively close.
        aqua_t, aqua_events = aqua.find_events(aoi, t0, t1, altitude_degrees=30.0)

        # Create empty arrays of today's passes and a dict for each pass.
        todays_aqua_passes = []
        aqua_passdict = {}

        # Loop through all Aqua events.
        for i in range(len(aqua_events)):

            # Get first event type and time.
            event = aqua_events[i]
            ti = aqua_t[i]

            # If this is the first event in the list and it's of type "overpass", this means that the rise of this
            # pass occurred on the previous day. This is an edge case when the area of interest is near the
            # international date line. The rise location values are set to nan, since these aren't necessary or easy to
            # retrieve in this case. The location of the satellite at the overpass point, the overpass time, and the
            # distance to the aoi at this point are recorded in the dictionary.
            if event == 1 and i == 0:
                aqua_passdict['rise_lat'] = float('nan')
                aqua_passdict['rise_lon'] = float('nan')
                difference = aqua - aoi
                topocentric = difference.at(ti)
                alt, az, distance = topocentric.altaz()
                aqua_passdict['distance'] = distance.km
                aqua_passdict['time'] = ti.utc_strftime('%Y %b %d %H:%M:%S')

                geocentric = aqua.at(ti)
                overlat, overlon = wgs84.latlon_of(geocentric)
                aqua_passdict['over_lat'] = overlat.degrees
                aqua_passdict['over_lon'] = overlon.degrees

            # This is another edge case when the first event in the list is of type "set". This occurs when the
            # satellite rises and passes over on the previous day, then sets today. If the overpass occurred yesterday,
            # it is already counted yesterday, so this event is ignored.
            elif event == 2 and i == 0:
                pass

            # This is an edge case where the final event of the day is an overpass. If so, this means that the satellite
            # sets on the next day, so the set location is set to nan.
            elif i == len(aqua_events) and event == 1:
                difference = aqua - aoi

                topocentric = difference.at(ti)

                alt, az, distance = topocentric.altaz()
                aqua_passdict['distance'] = distance.km
                aqua_passdict['time'] = ti.utc_strftime('%Y %b %d %H:%M:%S')

                geocentric = aqua.at(ti)
                overlat, overlon = wgs84.latlon_of(geocentric)
                aqua_passdict['over_lat'] = overlat.degrees
                aqua_passdict['over_lon'] = overlon.degrees

                aqua_passdict['set_lat'] = float('nan')
                aqua_passdict['set_lon'] = float('nan')
                todays_aqua_passes.append(aqua_passdict)

            # If the event is a 'rise', save the rise longitude and latitude to the dict representing the pass.
            elif event == 0:
                aqua_passdict = {}
                geocentric = aqua.at(ti)
                riselat, riselon = wgs84.latlon_of(geocentric)
                aqua_passdict['rise_lat'] = riselat.degrees
                aqua_passdict['rise_lon'] = riselon.degrees

            # If the event is an 'overpass', save the overpass longitude and latitude, the time, and the minimum
            # distance to the aoi to the dict representing the pass.
            elif event == 1:
                difference = aqua - aoi

                topocentric = difference.at(ti)

                alt, az, distance = topocentric.altaz()
                aqua_passdict['distance'] = distance.km
                aqua_passdict['time'] = ti.utc_strftime('%Y %b %d %H:%M:%S')

                geocentric = aqua.at(ti)
                overlat, overlon = wgs84.latlon_of(geocentric)
                aqua_passdict['over_lat'] = overlat.degrees
                aqua_passdict['over_lon'] = overlon.degrees

            # If the event is a 'set', save the set longitude and latitude to the dict representing the pass.
            else:
                geocentric = terra.at(ti)
                setlat, setlon = wgs84.latlon_of(geocentric)
                aqua_passdict['set_lat'] = setlat.degrees
                aqua_passdict['set_lon'] = setlon.degrees
                todays_aqua_passes.append(aqua_passdict)

        # Assume the minimum distance for a pass is very large until proven otherwise.
        least_distance = math.inf

        # Time of closest pass.
        closest_time = ""

        # Loop through each pass in the list of passes.
        for pass_dict in todays_aqua_passes:

            # If the rising latitude of the pass is less than the overpass latitude of the pass, this is a ascending
            # pass, which means a daytime pass for Aqua. These are the only passes of interest. If, furthermore, it's
            # less than the current minimum distance between the overpass and the aoi, this is the new closest pass.
            if not np.isnan(pass_dict['rise_lat']):
                if (pass_dict['rise_lat'] < pass_dict['over_lat']) and (pass_dict['distance'] < least_distance):
                    least_distance = pass_dict['distance']
                    closest_time = pass_dict['time']

            # For some edge case passes, rise latitude values are undefined. In this case, compare overpass latitude and
            # set latitude values. If the overpass latitude is less than the set latitude, it's an ascending pass.
            else:
                if (pass_dict['set_lat'] > pass_dict['over_lat']) and (pass_dict['distance'] < least_distance):
                    least_distance = pass_dict['distance']
                    closest_time = pass_dict['time']

        # Save this closest time as determined above.
        aqua_closest = closest_time.split(' ')[3]

        # Repeat for Terra.
        terra_t, terra_events = terra.find_events(aoi, t0, t1, altitude_degrees=30.0)

        todays_terra_passes = []
        terra_passdict = {}

        # Do the same as above for all terra events.
        for i in range(len(terra_events)):

            event = terra_events[i]
            ti = terra_t[i]

            if event == 1 and i == 0:
                terra_passdict['rise_lat'] = float('nan')
                terra_passdict['rise_lon'] = float('nan')
                difference = terra - aoi
                topocentric = difference.at(ti)
                alt, az, distance = topocentric.altaz()
                terra_passdict['distance'] = distance.km
                terra_passdict['time'] = ti.utc_strftime('%Y %b %d %H:%M:%S')

                geocentric = terra.at(ti)
                overlat, overlon = wgs84.latlon_of(geocentric)
                terra_passdict['over_lat'] = overlat.degrees
                terra_passdict['over_lon'] = overlon.degrees
            elif event == 2 and i == 0:
                pass
            elif i == len(terra_events) and event == 1:
                difference = terra - aoi

                topocentric = difference.at(ti)

                alt, az, distance = topocentric.altaz()
                terra_passdict['distance'] = distance.km
                terra_passdict['time'] = ti.utc_strftime('%Y %b %d %H:%M:%S')

                geocentric = terra.at(ti)
                overlat, overlon = wgs84.latlon_of(geocentric)
                terra_passdict['over_lat'] = overlat.degrees
                terra_passdict['over_lon'] = overlon.degrees

                terra_passdict['set_lat'] = float('nan')
                terra_passdict['set_lon'] = float('nan')
                todays_terra_passes.append(terra_passdict)

            elif event == 0:
                terra_passdict = {}
                geocentric = terra.at(ti)
                riselat, riselon = wgs84.latlon_of(geocentric)
                terra_passdict['rise_lat'] = riselat.degrees
                terra_passdict['rise_lon'] = riselon.degrees
            elif event == 1:
                difference = terra - aoi

                topocentric = difference.at(ti)

                alt, az, distance = topocentric.altaz()
                terra_passdict['distance'] = distance.km
                terra_passdict['time'] = ti.utc_strftime('%Y %b %d %H:%M:%S')

                geocentric = terra.at(ti)
                overlat, overlon = wgs84.latlon_of(geocentric)
                terra_passdict['over_lat'] = overlat.degrees
                terra_passdict['over_lon'] = overlon.degrees
            else:
                geocentric = terra.at(ti)
                setlat, setlon = wgs84.latlon_of(geocentric)
                terra_passdict['set_lat'] = setlat.degrees
                terra_passdict['set_lon'] = setlon.degrees
                todays_terra_passes.append(terra_passdict)

        least_distance = math.inf

        closest_time = ""

        for pass_dict in todays_terra_passes:

            # If the rising latitude of the pass is greater than the overpass latitude of the pass, this is a descending
            # pass, which means a daytime pass for Terra. These are the only passes of interest. If, furthermore, it's
            # less than the current minimum distance between the overpass and the aoi, this is the new closest pass.
            if not np.isnan(pass_dict['rise_lat']):
                if (pass_dict['rise_lat'] > pass_dict['over_lat']) and (pass_dict['distance'] < least_distance):
                    least_distance = pass_dict['distance']
                    closest_time = pass_dict['time']

            # For some edge case passes, rise latitude values are undefined. In this case, compare overpass latitude and
            # set latitude values. If the overpass latitude is greater than the set latitude, it's a descending pass.
            else:
                if (pass_dict['set_lat'] < pass_dict['over_lat']) and (pass_dict['distance'] < least_distance):
                    least_distance = pass_dict['distance']
                    closest_time = pass_dict['time']


        terra_closest = closest_time.split(' ')[3]

        # Add closest passes of the day to array of passes.
        rows.append(['-'.join(today), aqua_closest, terra_closest])

        # If today is the end date, break, otherwise set the current day to tomorrow.
        if np.array_equiv(tomorrow, end_date):
            break
        else:
            today = getNextDay(today)
            tomorrow = getNextDay(today)

    csvwrite(start_date, end_date, lat, long, rows, csvoutpath)

# Write CSV of all pass information.
def csvwrite(startdate, enddate, lat, lon, rows, outpath):
    fields = ['Date', 'Aqua pass time', 'Terra pass time']
    filename = f"{outpath}/passtimes_lat{lat}_lon{lon}_{''.join(startdate)}_{''.join(enddate)}.csv"

    with open(filename, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)

# Returns the date after a given date.
def getNextDay(date):
    month = int(date[0])
    year = int(date[2])
    day = int(date[1])
    monthDays = daysInMonth(month, year)

    nextmonth = month
    nextday = day
    nextyear = year

    if day == monthDays:
        nextday = 1
        nextmonth += 1
    else:
        nextday += 1
    if month == 12 and day == 31:
        nextyear += 1
        nextmonth = 1

    nextyearstr = str(nextyear)

    if nextday < 10:
        nextdaystr = '0' + str(nextday)
    else:
        nextdaystr = str(nextday)
    if nextmonth < 10:
        nextmonthstr = '0' + str(nextmonth)
    else:
        nextmonthstr = str(nextmonth)

    return [nextmonthstr, nextdaystr, nextyearstr]

# Returns the number of days in a certain month.
def daysInMonth(month, year):
    if month in {1, 3, 5, 7, 8, 10, 12}:
        return 31
    if month == 2:
        if year % 4 == 0:
            return 29
        return 28
    return 30

if __name__ == '__main__':
    main()