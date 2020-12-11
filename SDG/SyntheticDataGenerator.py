import random

def generate():
    # create a train of latitude
    lat = START_LAT
    lat_train = []
    while lat < END_LAT:
        if (round(lat, TWO) + CELL_SIZE) > END_LAT:
            lat_train.append((round(lat, TWO), END_LAT))
        else:
            lat_train.append((round(lat, TWO), round(lat + CELL_SIZE, TWO)))
        lat += CELL_SIZE

    # create a train of longitude
    lon = START_LON
    lon_train = []
    while lon < END_LON:
        if (round(lon, TWO) + CELL_SIZE) > END_LON:
            lon_train.append((round(lon, TWO), END_LON))
        else:
            lon_train.append((round(lon, TWO), round(lon + CELL_SIZE, TWO)))
        lon += CELL_SIZE

    # make NUMBER_OF_CENTERS distinct centers
    centers = []
    center_cells = {}
    count = ONE
    while count <= NUMBER_OF_CENTERS:
        lat_index = round(random.random() * (len(lat_train) - ONE))
        lon_index = round(random.random() * (len(lon_train) - ONE))
        key = str(lat_index) + ',' + str(lon_index)
        if not center_cells.get(key):
            center_cells[key] = (lat_train[lat_index], lon_train[lon_index])
            center_lat = round((lat_train[lat_index][STR_IDX] + lat_train[lat_index][END_IDX]) / TWO, TWO)
            center_lon = round((lon_train[lon_index][STR_IDX] + lon_train[lon_index][END_IDX]) / TWO, TWO)
            centers.append((center_lat, center_lon))
            print('Center ' + str(count) + ': ' + str((center_lat, center_lon)))
            count += ONE

    # create events for each center randomly inside it
    events = []
    for center in center_cells.items():
        lat_start = center[VAL_IDX][LAT_IDX][STR_IDX]
        lon_start = center[VAL_IDX][LON_IDX][STR_IDX]
        for i in range(MIN_CENTER_FREQUENCY):
            rand_lat = random.random() * CELL_SIZE + lat_start
            rand_lon = random.random() * CELL_SIZE + lon_start
            events.append((rand_lat, rand_lon))

    # create remaining events randomly through all the map
    for i in range(NUMBER_OF_EVENTS - NUMBER_OF_CENTERS*MIN_CENTER_FREQUENCY):
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

LAT_IDX = 0
LON_IDX = 1
KEY_IDX = 0
VAL_IDX = 1
STR_IDX = 0
END_IDX = 1

#############################################################
####################### PARAMETERS ##########################
#############################################################
NUMBER_OF_EVENTS = 10000
NUMBER_OF_CENTERS = 20
MIN_CENTER_FREQUENCY = 450
CELL_SIZE = 2 #precision for cell size is two floating point digits

#############################################################
########################### RUN #############################
#############################################################
if NUMBER_OF_CENTERS*MIN_CENTER_FREQUENCY > NUMBER_OF_EVENTS:
    print('Constraints are not valid simultaneously. '
          'NUMBER_OF_CENTERS*MIN_CENTER_FREQUENCY should be equal to NUMBER_OF_EVENTS')
else:
    centers, events = generate()
    f = open('synthetic_events.txt', 'w')
    f.write('NUMBER_OF_EVENTS=' + str(NUMBER_OF_EVENTS) + '\n')
    f.write('NUMBER_OF_CENTERS=' + str(NUMBER_OF_CENTERS) + '\n')
    f.write('MIN_CENTER_FREQUENCY=' + str(MIN_CENTER_FREQUENCY) + '\n')
    f.write('CELL_SIZE=' + str(CELL_SIZE) + '\n')

    f.write('=================  CENTERS  ================\n')
    count = ONE
    for center in centers:
        f.write('Center ' + str(count) + ': ' + str(center) + '\n')
        count += ONE

    f.write('=================  EVENTS  ==================\n')
    for event in events:
        f.write(str(event[LAT_IDX]) + ',' + str(event[LON_IDX]) + '\n')
    f.close()