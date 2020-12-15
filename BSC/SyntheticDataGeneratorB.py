import random
import math

def generate():
    # create a train of latitude
    lat = START_LAT + ONE
    lat_train = []
    while lat < (END_LAT - ONE):
        if (round(lat, TWO) + LARGE_CELL_SIZE) > END_LAT:
            lat_train.append((round(lat, TWO), END_LAT))
        else:
            lat_train.append((round(lat, TWO), round(lat + LARGE_CELL_SIZE, TWO)))
        lat += LARGE_CELL_SIZE

    # create a train of longitude
    lon = START_LON + ONE
    lon_train = []
    while lon < (END_LON - ONE):
        if (round(lon, TWO) + LARGE_CELL_SIZE) > END_LON:
            lon_train.append((round(lon, TWO), END_LON))
        else:
            lon_train.append((round(lon, TWO), round(lon + LARGE_CELL_SIZE, TWO)))
        lon += LARGE_CELL_SIZE

    # make NUMBER_OF_CENTERS distinct centers
    centers = []
    center_cells = {}
    count = ONE
    is_b = False
    while count <= NUMBER_OF_CENTERS + ONE:
        lat_index = round(random.random() * (len(lat_train) - ONE))
        lon_index = round(random.random() * (len(lon_train) - ONE))
        key = str(lat_index) + ',' + str(lon_index)
        if not center_cells.get(key):
            center_lat = round((lat_train[lat_index][STR_IDX] + lat_train[lat_index][END_IDX]) / TWO, TWO)
            center_lon = round((lon_train[lon_index][STR_IDX] + lon_train[lon_index][END_IDX]) / TWO, TWO)
            center_cells[key] = (lat_train[lat_index], lon_train[lon_index], is_b)
            centers.append((center_lat, center_lon, is_b))
            print('Center ' + str(count) + ': ' + str((center_lat, center_lon, is_b)))
            count += ONE
            if not is_b:
                is_b = True

    # create events for each center randomly inside it
    events = []
    for center in center_cells.items():
        lat_start = center[VAL_IDX][LAT_IDX][STR_IDX]
        lon_start = center[VAL_IDX][LON_IDX][STR_IDX]
        is_b = center[VAL_IDX][B_IDX]
        if not is_b:
            rand_lat = random.random() * int(LARGE_CELL_SIZE / TWO) + lat_start
            rand_lon = random.random() * int(LARGE_CELL_SIZE / TWO) + lon_start
            for i in range(math.floor(B_CENTER_FREQUENCY/TWO)):
                events.append((rand_lat, rand_lon))
            for i in range(math.floor(B_CENTER_FREQUENCY/TWO)):
                rand_lat = random.random() * int(LARGE_CELL_SIZE / TWO) + lat_start
                rand_lon = random.random() * int(LARGE_CELL_SIZE / TWO) + lon_start
                events.append((rand_lat, rand_lon))
        else:
            for i in range(math.floor(MIN_CENTER_FREQUENCY/FOUR)):
                rand_lat = random.random() * SMALL_CELL_SIZE + lat_start + \
                           (int(LARGE_CELL_SIZE/SMALL_CELL_SIZE) - ONE) * SMALL_CELL_SIZE
                rand_lon = random.random() * SMALL_CELL_SIZE + lon_start
                events.append((rand_lat, rand_lon))
            for i in range(math.floor(THREE*MIN_CENTER_FREQUENCY/FOUR)):
                rand_lat = random.random() * LARGE_CELL_SIZE + lat_start
                rand_lon = random.random() * LARGE_CELL_SIZE + lon_start
                events.append((rand_lat, rand_lon))

    # create remaining events randomly through all the map
    for i in range(NUMBER_OF_EVENTS - NUMBER_OF_CENTERS * MIN_CENTER_FREQUENCY - B_CENTER_FREQUENCY):
        rand_lat = random.random() * (END_LAT - START_LAT) + START_LAT
        rand_lon = random.random() * (END_LON - START_LON) + START_LON
        events.append((rand_lat, rand_lon))

    return centers, events

#############################################################
######################## CONSTANTS ##########################
#############################################################
START_LAT = -90
END_LAT = +90
START_LON = -180
END_LON = +180

ZERO = 0
ONE = 1
TWO = 2
THREE = 3
FOUR = 4

LAT_IDX = 0
LON_IDX = 1
B_IDX = 2
KEY_IDX = 0
VAL_IDX = 1
STR_IDX = 0
END_IDX = 1

LARGE_CELL_SIZE = 2
SMALL_CELL_SIZE = 0.1

#############################################################
####################### PARAMETERS ##########################
#############################################################
NUMBER_OF_EVENTS = 4000
NUMBER_OF_CENTERS = 4
MIN_CENTER_FREQUENCY = 800
B_CENTER_FREQUENCY = 600

#############################################################
########################### RUN #############################
#############################################################
if (NUMBER_OF_CENTERS * MIN_CENTER_FREQUENCY + B_CENTER_FREQUENCY) > NUMBER_OF_EVENTS:
    print('Constraints are not valid simultaneously. '
          'NUMBER_OF_CENTERS * MIN_CENTER_FREQUENCY + B_CENTER_FREQUENCY should be <= NUMBER_OF_EVENTS')
else:
    centers, events = generate()
    f = open('synthetic_events_b.txt', 'w')
    f.write('NUMBER_OF_EVENTS=' + str(NUMBER_OF_EVENTS) + '\n')
    f.write('NUMBER_OF_CENTERS=' + str(NUMBER_OF_CENTERS) + '\n')
    f.write('MIN_CENTER_FREQUENCY=' + str(MIN_CENTER_FREQUENCY) + '\n')
    f.write('LARGE_CELL_SIZE=' + str(LARGE_CELL_SIZE) + '\n')

    f.write('=================  CENTERS  ================\n')
    count = ONE
    for center in centers:
        f.write('Center ' + str(count) + ': ' + str(center) + '\n')
        count += ONE

    f.write('=================  EVENTS  ==================\n')
    for event in events:
        f.write(str(event[LAT_IDX]) + ',' + str(event[LON_IDX]) + '\n')
    f.close()