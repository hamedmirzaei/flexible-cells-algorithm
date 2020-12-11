import random

def generate():
    # create a train of longitude
    lon = START_LON
    lon_train = []
    while lon <= END_LON:
        if round(lon, TWO) + CELL_SIZE > END_LON:
            lon_train.append((round(lon, TWO), END_LON))
        else:
            lon_train.append((round(lon, TWO), round(lon + CELL_SIZE, TWO)))
        lon += CELL_SIZE

    # create a train of latitude
    lat = START_LAT
    lat_train = []
    while lat <= END_LAT:
        if round(lat, TWO) + CELL_SIZE > END_LAT:
            lat_train.append((round(lat, TWO), END_LAT))
        else:
            lat_train.append((round(lat, TWO), round(lat + CELL_SIZE, TWO)))
        lat += CELL_SIZE

    # make NUMBER_OF_CENTERS distinct centers
    centers = []
    center_cells = {}
    count = ONE
    while count <= NUMBER_OF_CENTERS:
        rand_lon = round(random.random() * len(lon_train))
        rand_lat = round(random.random() * len(lat_train))
        key = str(rand_lon) + ',' + str(rand_lat)
        if not center_cells.get(key):
            center_cells[key] = (lon_train[rand_lon], lat_train[rand_lat])
            center_lon = round((lon_train[rand_lon][STR_IDX] + lon_train[rand_lon][END_IDX]) / 2, TWO)
            center_lat = round((lat_train[rand_lat][STR_IDX] + lat_train[rand_lat][END_IDX]) / 2, TWO)
            centers.append((center_lon, center_lat))
            print('Center ' + str(count) + ': ' + str((center_lon, center_lat)))
            count += ONE

    # create events for each center randomly inside it
    events = []
    for center in center_cells.items():
        lon_start = center[VAL_IDX][LON_IDX][STR_IDX]
        lat_start = center[VAL_IDX][LAT_IDX][STR_IDX]
        for i in range(MIN_CENTER_FREQUENCY):
            rand_lon = random.random() * CELL_SIZE + lon_start
            rand_lat = random.random() * CELL_SIZE + lat_start
            events.append((rand_lon, rand_lat))

    # create remaining events randomly through all the map
    for i in range(NUMBER_OF_EVENTS - NUMBER_OF_CENTERS*MIN_CENTER_FREQUENCY):
        rand_lon = random.random() * (END_LON - START_LON) + START_LON
        rand_lat = random.random() * (END_LAT - START_LAT) + START_LAT
        events.append((rand_lon, rand_lat))

    return centers, events

#############################################################
######################## CONSTANTS ##########################
#############################################################
START_LON = -180
END_LON = +180
START_LAT = -90
END_LAT = +90

ZERO = 0
ONE = 1
TWO = 2

LON_IDX = 0
LAT_IDX = 1
KEY_IDX = 0
VAL_IDX = 1
STR_IDX = 0
END_IDX = 1

#############################################################
####################### PARAMETERS ##########################
#############################################################
NUMBER_OF_EVENTS = 100000
NUMBER_OF_CENTERS = 20
MIN_CENTER_FREQUENCY = 1000
CELL_SIZE = 0.1 #precision for cell size is two floating point digits

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
        f.write(str(event[LON_IDX]) + ',' + str(event[LAT_IDX]) + '\n')
    f.close()