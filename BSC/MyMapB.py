################################################################################
################################## IMPORTS #####################################
################################################################################
import queue
import decimal
import scipy.optimize as optimize
from geopy.distance import great_circle
import math
import numpy as np
import time


################################################################################
################################# FUNCTIONS ####################################
################################################################################
# calculate distance of two points in kilometers
def distance(lat1,lon1,lat2,lon2):
    # first move the points to original space
    return great_circle((lat1, lon1), (lat2, lon2)).kilometers

def is_event_in_cell(event, cell):
    lat = event[LAT_IDX]
    lon = event[LON_IDX]
    # check if event falls inside the cell
    if lat >= cell[Y_STR_IDX] and lat <= cell[Y_END_IDX] and \
        lon >= cell[X_STR_IDX] and lon <= cell[X_END_IDX]:
        return True
    return False

# count number of events falls inside the input cell
def count_events(cell, related_events):
    cnt = ZERO
    for event in related_events: # loop over all the events
        if is_event_in_cell(event, cell):
            cnt += ONE
    return cnt

# read all the events from the input csv file and
# convert it to coordinate where (0,0) is the top
# left corner
def read_events():
    events = []# empty list of events
    with open(events_cvs_file) as f:# open the input file
        for line in f:# read file line by line
            record = line.replace('\n', '')# remove \n from end of the line
            # split the line by comma, separate lat and lon
            coordinate = record.split(',')
            events.append((float(coordinate[LAT_IDX]), float(coordinate[LON_IDX])))# add the (lat, lon) to the events list
    return events

# Backstrom's single-center algorithm
def find_center():
    global events
    print('Finding center based on Backstrom\'s single-center algorithm ...')
    cell_counts_map = {}
    cells = []
    x_start = x0
    while x_start < dx0:
        y_start = y0
        while y_start < dy0:
            key = str(x_start) + KEY_DELIMITER + str(y_start)
            cell_counts_map[key] = ZERO
            y_start += TWO
        x_start += TWO

    for event in events:
        y_floor = math.floor(event[LAT_IDX])
        if y_floor % TWO == ONE:
            y_floor -= ONE
        x_floor = math.floor(event[LON_IDX])
        if x_floor % TWO == ONE:
            x_floor -= ONE
        key = str(x_floor) + KEY_DELIMITER + str(y_floor)
        cell_counts_map[key] += ONE

    max_item = (None, -ONE)
    for ccm in cell_counts_map.items():
        if ccm[ONE] > max_item[ONE]:
            max_item = ccm

    max_item_x = int(max_item[ZERO].split(KEY_DELIMITER)[ZERO])
    max_item_y = int(max_item[ZERO].split(KEY_DELIMITER)[ONE])
    cells = []
    x_start = max_item_x
    while x_start < (max_item_x + TWO):
        y_start = max_item_y
        while y_start < (max_item_y + TWO):
            cells.append((x_start, round(x_start + ONETENTH, TWO), ONETENTH, y_start, round(y_start + ONETENTH, TWO), ONETENTH))
            y_start = round(y_start + ONETENTH, ROUND_DIGITS)
        x_start = round(x_start + ONETENTH, ROUND_DIGITS)

    final_cell = cells[ZERO]
    final_count = ZERO
    for i in range(len(cells)):
        count = count_events(cells[i], events)
        if count > final_count:
            final_count = count
            final_cell = cells[i]
    return (final_cell, final_count)

# find optimum values for alpha based on spatial variation paper
def find_optimum_params(center_cell_count):
    center_cell = center_cell_count[CELL_IDX]
    center_x = round((center_cell[X_STR_IDX] + center_cell[X_END_IDX]) / TWO, ROUND_DIGITS)
    center_y = round((center_cell[Y_STR_IDX] + center_cell[Y_END_IDX]) / TWO, ROUND_DIGITS)
    center_count = center_cell_count[COUNT_IDX]

    cell_counts_map = {}
    cell_dists_map = {}

    x_start = x0
    while x_start < dx0:
        y_start = y0
        while y_start < dy0:
            key = str(int(x_start * TEN)) + KEY_DELIMITER + str(int(y_start * TEN))
            cell_x = round(x_start + ONETWENTY, ROUND_DIGITS)
            cell_y = round(y_start + ONETWENTY, ROUND_DIGITS)
            cell_counts_map[key] = ZERO
            cell_dists_map[key] = 0#distance(center_y, center_x, cell_y, cell_x)
            y_start = round(y_start + ONETENTH, ROUND_DIGITS)
        x_start = round(x_start + ONETENTH, ROUND_DIGITS)

    for event in events:
        x_floor = math.floor(event[LON_IDX] * TEN)
        y_floor = math.floor(event[LAT_IDX] * TEN)
        key = str(x_floor) + KEY_DELIMITER + str(y_floor)
        cell_counts_map[key] += ONE

    # optimizing function which will be executed repeatedly to optimize x
    # x is alpha
    def f(x):
        result = ZERO
        logs = []
        one_minus_logs = []
        # summation over all the cells
        for cell_key in cell_counts_map.keys():
            if cell_dists_map[cell_key] == ZERO:
                continue
            cdx = center_count * (cell_dists_map[cell_key] ** (-x))
            # if there is at least one relevant event within the cell
            if cell_counts_map[cell_key] != ZERO:
                logs.append(cdx)
            else:
                one_minus_logs.append(cdx)

        # calculate sum for later normalization
        sum_all = np.sum(logs) + np.sum(one_minus_logs)

        # normalizing
        logs = np.divide(logs, sum_all)
        for ii in range(len(logs)):
            if logs[ii] == ZERO:
                logs[ii] = TEN**(-TEN)
        logs = np.log(logs)
        one_minus_logs = ONE - np.divide(one_minus_logs, sum_all)
        for ii in range(len(one_minus_logs)):
            if one_minus_logs[ii] == ZERO:
                one_minus_logs[ii] = TEN**(-TEN)
        one_minus_logs = np.log(one_minus_logs)

        # calculating the result
        result += np.sum(logs)
        result += np.sum(one_minus_logs)

        return (-ONE) * result

    # optimize alpha
    res = optimize.minimize_scalar(f, bounds=(0, 5), method="bounded")

    # res.x would be our alpha
    return (center_count, res.x)

################################################################################
################################## CONSTANTS ###################################
################################################################################
# x means lon
# y means lat
x0 = -180 # top left corner considered as (0,0) coordinate
y0 = -90 # top left corner considered as (0,0) coordinate
dx0 = 180 # from -180 to 180
dy0 = 90 # from -90 to 90 (for square cells, we suppose its 360 wide)

X_STR_IDX = 0
X_END_IDX = 1
X_LEN = 2
Y_STR_IDX = 3
Y_END_IDX = 4
Y_LEN = 5

CELL_IDX = 0
COUNT_IDX = 1
EVENTS_IDX = 2

LAT_IDX = 0
LON_IDX = 1

EIGHTEEN_HUNDRED = 180
NINETY = 90
TWELWE = 12
TEN = 10
FOUR = 4
TWO = 2
ONE = 1
ONETENTH = 0.1
ONETWENTY = 0.05
ZERO = 0
ROUND_DIGITS = 2

MAP_STR = 0
MAP_END = 1

KEY_DELIMITER = ','

################################################################################
################################# PARAMETERS ###################################
################################################################################
MAX_OPTIMIZATION_LOOP = 20
events_cvs_file = 'map.csv'
# events file should be a csv of lat,lon


################################################################################
################################### RUN ########################################
################################################################################
print('Maximum number of loops in optimization = ' + str(MAX_OPTIMIZATION_LOOP))
print('Input event file name = ' + str(events_cvs_file))
print('')

start_time = []
events = read_events()
results_celldists = {}
results_eventcells = {}

found_center = find_center()
# pritn some headlines
print('')
print('Center template is ((min x, max x, length, min y, max y, width, type), count)')
print('')
# print all centers and parameters
print('Center: ' + str(found_center))

params = find_optimum_params(found_center)
print('(C,alpha): ' + str(params))