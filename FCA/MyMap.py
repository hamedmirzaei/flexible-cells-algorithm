################################################################################
################################## IMPORTS #####################################
################################################################################
import queue
import decimal
import scipy.optimize as optimize
from geopy.distance import great_circle
import math


################################################################################
################################# FUNCTIONS ####################################
################################################################################
# convert a point from real space to our space
def real_to_our_space_point(real_point):
    return (-real_point[LAT_IDX] + NINETY, real_point[LON_IDX] + EIGHTEEN_HUNDRED)

# convert a point from our space to real space
def our_to_real_space_point(our_point):
    return (-our_point[LAT_IDX] + NINETY, our_point[LON_IDX] - EIGHTEEN_HUNDRED)

# convert a result from our space to real space
def our_to_real_space_result(our_result):
    our_cell = our_result[CELL_IDX]
    top_left = our_to_real_space_point((our_cell[Y_STR_IDX], our_cell[X_STR_IDX]))
    bot_right = our_to_real_space_point((our_cell[Y_END_IDX], our_cell[X_END_IDX]))
    return ((top_left[LON_IDX], bot_right[LON_IDX],our_cell[X_LEN],
            top_left[LAT_IDX], bot_right[LAT_IDX], our_cell[Y_LEN],
            our_cell[TYPE_IDX]), our_result[COUNT_IDX])

# count number of events falls inside the input cell + min/max boundaries
def count_events_minmax(cell):
    #TODO: make it more efficient number of events
    # falling inside this specific cell
    min_max = [dx0,x0,ZERO,dy0,y0,ZERO]
    res = ZERO
    for event in events: # loop over all the events
        lat = event[LAT_IDX]
        lon = event[LON_IDX]
        # check if event falls inside the cell
        if lon >= cell[X_STR_IDX] and lon <= cell[X_END_IDX] and \
                lat >= cell[Y_STR_IDX] and lat <= cell[Y_END_IDX]:
            res += ONE

            #update min/max of lat/lon of events in this cell
            if lon < min_max[MIN_X_IDX]:
                min_max[MIN_X_IDX] = lon
            if lon > min_max[MAX_X_IDX]:
                min_max[MAX_X_IDX] = lon
            if lat < min_max[MIN_Y_IDX]:
                min_max[MIN_Y_IDX] = lat
            if lat > min_max[MAX_Y_IDX]:
                min_max[MAX_Y_IDX] = lat

    # make a small offset when all the events for this cell
    # are from one single coordinate
    if min_max[MIN_X_IDX] == min_max[MAX_X_IDX]:
        min_max[MIN_X_IDX] -= AREA_OFFSET
        min_max[MAX_X_IDX] += AREA_OFFSET
    if min_max[MIN_Y_IDX] == min_max[MAX_Y_IDX]:
        min_max[MIN_Y_IDX] -= AREA_OFFSET
        min_max[MAX_Y_IDX] += AREA_OFFSET

    # calculate diff of x and y
    min_max[DIFF_X_IDX] = round(min_max[MAX_X_IDX] - min_max[MIN_X_IDX], ROUND_DIGITS)
    min_max[DIFF_Y_IDX] = round(min_max[MAX_Y_IDX] - min_max[MIN_Y_IDX], ROUND_DIGITS)

    return res, min_max

# count number of events falls inside the input cell
def count_events(cell):
    #TODO: make it more efficient
    res = ZERO
    for event in events: # loop over all the events
        lat = event[LAT_IDX]
        lon = event[LON_IDX]
        # check if event falls inside the cell
        if lon >= cell[X_STR_IDX] and lon <= cell[X_END_IDX] and \
                lat >= cell[Y_STR_IDX] and lat <= cell[Y_END_IDX]:
            res += ONE
    return res

# read all the events from the input csv file and
# convert it to coordinate where (0,0) is the top
# left corner
def read_events():
    events = []# epty list of events
    with open(events_cvs_file) as f:# open the input file
        for line in f:# read file line by line
            record = line.replace('\n', '')# remove \n from end of the line
            # split the line by comma, separate lat and lon
            coordinate = record.split(',')
            # get lat and lon and convert it into our system
            our_point = real_to_our_space_point((float(coordinate[LAT_IDX]), float(coordinate[LON_IDX])))
            events.append((our_point[LAT_IDX], our_point[LON_IDX]))# add the (lat, lon) to the events list
    return events

# find all the cells met the provided thresholds and
# return them. It uses a BFS like structure to move
# deeper for cells
def find_centers():
    print('Finding centers based on flexible cells algorithm ...')
    bfsq = queue.Queue()# FIFO queue for processing cells
    # first cell which contains the full area
    cell = (x0, x0+dx0, dx0, y0, y0+dy0, dy0, CELL_TYPE_1)
    count = len(events)# all the events are inside the full area
    if count < TF_THR:# if there is not enough events
        return queue.Queue()# return empty queue
    bfsq.put((cell, count))# put the (full area cell, count) into the queue
    while bfsq:# do it til there is an item in the queue
        popItem = bfsq.get()# get an item from queue
        popCell = popItem[CELL_IDX]# get the cell info

        # if the cell met th size threshold
        if (popCell[X_LEN] <= SIZE_THR) and (popCell[Y_LEN] <= SIZE_THR):
            bfsq.put(popItem)# put the removed item back
            return bfsq# return final result

        if popCell[TYPE_IDX] == CELL_TYPE_1: #four corners, cells 1,3,7,9
            # move over nine different cells of size half of the original cell
            # move over horizontal cells
            for i in [ZERO, popCell[X_LEN] / QUARTER, popCell[X_LEN] / HALF]:
                # move over vertical cells
                for j in [ZERO, popCell[Y_LEN] / QUARTER, popCell[Y_LEN] / HALF]:
                    type = CELL_TYPE_1
                    if j == popCell[Y_LEN] / QUARTER:# middle vertical cells
                        type = CELL_TYPE_3
                    if i == popCell[X_LEN] / QUARTER:# middle horizontal cells
                        type = CELL_TYPE_2
                        # middle vertical/horizontal cells
                        if j == popCell[Y_LEN] / QUARTER:
                            type = CELL_TYPE_4
                    # make a cell of half size
                    cell = (popCell[X_STR_IDX] + i,
                            popCell[X_STR_IDX] + i + popCell[X_LEN] / HALF,
                            popCell[X_LEN] / HALF,
                            popCell[Y_STR_IDX] + j,
                            popCell[Y_STR_IDX] + j + popCell[Y_LEN] / HALF,
                            popCell[Y_LEN] / HALF,
                            type)
                    # count number of events inside this cell
                    count, min_max = count_events_minmax(cell)
                    # if count meets the threshold, add the cell to queue
                    if count >= TF_THR:
                        bfsq.put((cell, count, min_max))

        elif popCell[TYPE_IDX] == CELL_TYPE_2: #two cells 2,8
            # for type_2 cells, just need to move over middle horizontal
            # cells, others are covered already
            i = popCell[X_LEN] / QUARTER
            # move over vertical cells
            for j in [ZERO, popCell[Y_LEN] / QUARTER, popCell[Y_LEN] / HALF]:
                # make a cell of half size
                cell = (popCell[X_STR_IDX] + i,
                        popCell[X_STR_IDX] + i + popCell[X_LEN] / HALF,
                        popCell[X_LEN] / HALF,
                        popCell[Y_STR_IDX] + j,
                        popCell[Y_STR_IDX] + j + popCell[Y_LEN] / HALF,
                        popCell[Y_LEN] / HALF,
                        CELL_TYPE_2)
                # count number of events inside this cell
                count, min_max = count_events_minmax(cell)
                # if count meets the threshold, add the cell to queue
                if count >= TF_THR:
                    bfsq.put((cell, count, min_max))

        elif popCell[TYPE_IDX] == CELL_TYPE_3: #two cells 4,6
            # for type_3 cells, just need to move over middle vertical
            # cells, others are covered already
            j = popCell[Y_LEN] / QUARTER
            # move over horizontal cells
            for i in [ZERO, popCell[X_LEN] / QUARTER, popCell[X_LEN] / HALF]:
                # make a cell of half size
                cell = (popCell[X_STR_IDX] + i,
                        popCell[X_STR_IDX] + i + popCell[X_LEN] / HALF,
                        popCell[X_LEN] / HALF,
                        popCell[Y_STR_IDX] + j,
                        popCell[Y_STR_IDX] + j + popCell[Y_LEN] / HALF,
                        popCell[Y_LEN] / HALF,
                        CELL_TYPE_3)
                # count number of events inside this cell
                count, min_max = count_events_minmax(cell)
                # if count meets the threshold, add the cell to queue
                if count >= TF_THR:
                    bfsq.put((cell, count, min_max))

        elif popCell[TYPE_IDX] == CELL_TYPE_4: #one cell 5
            # for type_4 cells, just need to move over the middle
            # cell, others are covered already
            i = popCell[X_LEN] / QUARTER
            j = popCell[Y_LEN] / QUARTER
            # make a cell of half size
            cell = (popCell[X_STR_IDX] + i,
                    popCell[X_STR_IDX] + i + popCell[X_LEN] / HALF,
                    popCell[X_LEN] / HALF,
                    popCell[Y_STR_IDX] + j,
                    popCell[Y_STR_IDX] + j + popCell[Y_LEN] / HALF,
                    popCell[Y_LEN] / HALF,
                    CELL_TYPE_4)
            # count number of events inside this cell
            count, min_max = count_events_minmax(cell)
            # if count meets the threshold, add the cell to queue
            if count >= TF_THR:
                bfsq.put((cell, count, min_max))
    return bfsq

# checks confliction between two cells
def has_conflict(cell1, cell2):
    if cell1[X_STR_IDX] <= cell2[X_STR_IDX]:
        if cell1[X_END_IDX] <= cell2[X_STR_IDX]:
            return False
        if cell1[Y_STR_IDX] >= cell2[Y_STR_IDX] and \
                cell1[Y_STR_IDX] <= cell2[Y_END_IDX]:
            return True
        if cell2[Y_STR_IDX] >= cell1[Y_STR_IDX] and \
                cell2[Y_STR_IDX] <= cell1[Y_END_IDX]:
            return True
    else:
        if cell2[X_END_IDX] <= cell1[X_STR_IDX]:
            return False
        if cell1[Y_STR_IDX] >= cell2[Y_STR_IDX] and \
                cell1[Y_STR_IDX] <= cell2[Y_END_IDX]:
            return True
        if cell2[Y_STR_IDX] >= cell1[Y_STR_IDX] and \
                cell2[Y_STR_IDX] <= cell1[Y_END_IDX]:
            return True
    return False

# prune the results by removing conflict cells from find_centers output
def remove_conflicts(temp_result):
    print('Keep the higher_count/lower_x/lower_y cells and remove conflicting cells to them ...')
    # first try to sort current results by:
    # 1) count
    # 2) x of top-left corner
    # 3) y of top-left corner
    sort_res = []
    while not temp_result.empty():
        cell_count = temp_result.get()
        # making keys based on mentioned order above, later
        # sort will be done using this generated keys
        sort_res.append({KEY_NAME_1: cell_count[COUNT_IDX],
                         KEY_NAME_2: cell_count[CELL_IDX][X_STR_IDX],
                         KEY_NAME_3: cell_count[CELL_IDX][Y_STR_IDX],
                         VAL_NAME: cell_count})

    # sort based on keys:
    #  1) in reverse order of counts so that highest counts comes first
    #  2) in direct order of top-left x so that the lowest x's comes first
    #  3) in direct order of top-left y so that the lowest y's comes first
    sort_res = sorted(sorted(sorted(sort_res,
                                    key=lambda x: x[KEY_NAME_3]),
                             key=lambda x: x[KEY_NAME_2]),
                      key=lambda x: x[KEY_NAME_1], reverse=True)

    # remove conflicts
    pruned_res = []
    for item in sort_res:
        conflict = False
        for i in pruned_res:
            if has_conflict(item[VAL_NAME][CELL_IDX], i[CELL_IDX]):
                conflict = True
                break
        if not conflict:
            pruned_res.append(item[VAL_NAME])
    return pruned_res

# compress the boundaries for each cell to the min max of events in it
def compress_boundaries(temp_result):
    print('Compress boundaries to min/max of events within it ...')
    bounded = queue.Queue()
    if COMPRESS_BOUNDARIES_ACTIVE:
        while not temp_result.empty():
            # get the item from queue
            cell_count_minmax = temp_result.get()
            # put the cell and count in new queue
            bounded.put((cell_count_minmax[MINMAX_IDX], cell_count_minmax[COUNT_IDX]))
    else:
        while not temp_result.empty():
            # get the item from queue
            cell_count_minmax = temp_result.get()
            # put the cell and count in new queue
            bounded.put((cell_count_minmax[CELL_IDX], cell_count_minmax[COUNT_IDX]))
    return bounded

# calculate distance of two points in kilometers
def distance(lat1,lon1,lat2,lon2):
    # first move the points to original space
    return great_circle(our_to_real_space_point((lat1, lon1)), our_to_real_space_point((lat2, lon2))).kilometers

# find optimum values for alpha based on spatial variation paper
def find_optimum_params(cell_count):
    center_cell = cell_count[CELL_IDX]
    center_count = cell_count[COUNT_IDX]

    cell_center_x = center_cell[X_STR_IDX] + center_cell[X_LEN] / 2
    cell_center_y = center_cell[Y_STR_IDX] + center_cell[Y_LEN] / 2

    # first we should divide map into cells of size result cell
    # (which is the input of the method)
    i_indices = [] # to store all possible horizontal start/end of cells
    i = center_cell[X_STR_IDX]# start from cell start x backward
    while i > x0:
        if i - center_cell[X_LEN] < x0:
            i_indices.append((x0, i))
        else:
            i_indices.append((i - center_cell[X_LEN], i))
        i -= center_cell[X_LEN]
    i = center_cell[X_END_IDX]# start from cell end x forward
    while i < dx0:
        if i + center_cell[X_LEN] > dx0:
            i_indices.append((i, dx0))
        else:
            i_indices.append((i, i + center_cell[X_LEN]))
        i += center_cell[X_LEN]
    # sort i indices based on start x
    sorted(i_indices, key=lambda x: x[MAP_STR]);

    j_indices = [] # to store all possible vertical start/end of cells
    j = center_cell[Y_STR_IDX]# start from cell start y upward
    while j > y0:
        if j - center_cell[Y_LEN] < y0:
            j_indices.append((y0, j))
        else:
            j_indices.append((j - center_cell[Y_LEN], j))
        j -= center_cell[Y_LEN]
    j = center_cell[Y_END_IDX]# start from cell end y downward
    while j < dy0_exact:
        if j + center_cell[Y_LEN] > dy0_exact:
            j_indices.append((j, dy0_exact))
        else:
            j_indices.append((j, j + center_cell[Y_LEN]))
        j += center_cell[Y_LEN]
    # sort j indices based on start y
    sorted(j_indices, key=lambda y: y[MAP_STR]);

    # move over all the cells and count events for each one
    map = []
    for i in i_indices:
        for j in j_indices:
            new_cell = (i[MAP_STR], i[MAP_END], i[MAP_END] - i[MAP_STR],
                        j[MAP_STR], j[MAP_END], j[MAP_END] - j[MAP_STR])
            new_count = count_events(new_cell)
            new_center_x = (i[MAP_STR] + i[MAP_END]) / HALF
            new_center_y = (j[MAP_STR] + j[MAP_END]) / HALF
            new_dist = distance(new_center_y, new_center_x, cell_center_y, cell_center_x)
            map.append(('', new_count, new_dist))
    def f(x):
        result = 0
        cnt = 0
        for mapcell in map:
            if mapcell[MAPCELL_COUNT_IDX] != 0:
                cnt += 1
                result += math.log(center_count * (mapcell[MAPCELL_DIST_IDX] ** (-x)))
            else:
                result += math.log(1 - center_count * (mapcell[MAPCELL_DIST_IDX] ** (-x)))
        return (-1) * result

    res = optimize.minimize_scalar(f, bounds=(0, 5), method="bounded")
    # res.x would be our alpha
    del i_indices
    del j_indices
    del map
    return (center_count, res.x)

################################################################################
################################# PARAMETERS ###################################
################################################################################
TF_THR = 6 # threshold for tweet frequency
SIZE_THR = 0.8 # threshold for area size of interest
COMPRESS_BOUNDARIES_ACTIVE = False
events_cvs_file = 'CMPUT_692\map.csv'


################################################################################
################################## CONSTANTS ###################################
################################################################################
x0 = 0 # top left corner considered as (0,0) coordinate
y0 = 0 # top left corner considered as (0,0) coordinate
dx0 = 360 # from -180 to 180
dy0 = 360 # from -90 to 90 (for square cells, we suppose its 360 wide)
dy0_exact = 180 # from -90 to 90

X_STR_IDX = 0
X_END_IDX = 1
X_LEN = 2
Y_STR_IDX = 3
Y_END_IDX = 4
Y_LEN = 5
TYPE_IDX = 6

CELL_IDX = 0
COUNT_IDX = 1
MINMAX_IDX = 2

CELL_TYPE_1 = 1
CELL_TYPE_2 = 2
CELL_TYPE_3 = 3
CELL_TYPE_4 = 4

LAT_IDX = 0
LON_IDX = 1

MIN_X_IDX = 0
MAX_X_IDX = 1
DIFF_X_IDX = 2
MIN_Y_IDX = 3
MAX_Y_IDX = 4
DIFF_Y_IDX = 5

KEY_NAME_1 = 'key1'
KEY_NAME_2 = 'key2'
KEY_NAME_3 = 'key3'
VAL_NAME = 'val'

QUARTER = 4
HALF = 2
ZERO = 0
ONE = 1
NINETY = 90
EIGHTEEN_HUNDRED = 180
ROUND_DIGITS = 5
AREA_OFFSET = 0.001

MAP_STR = 0
MAP_END = 1

MAPCELL_CELL_IDX = 0
MAPCELL_COUNT_IDX = 1
MAPCELL_DIST_IDX = 2

################################################################################
################################### RUN ########################################
################################################################################
print('Tweet frequency threshold = ' + str(TF_THR))
print('Area size of interest threshold = ' + str(SIZE_THR))
print('Compress boundaries active status = ' + str(COMPRESS_BOUNDARIES_ACTIVE))
print('Input event file name = ' + str(events_cvs_file))
print('')

events = read_events()
results = remove_conflicts(compress_boundaries(find_centers()))
print('')
print('There are ' + str(len(results)) + ' centers - Calculating optimum parameters for each one:')
print('Center template is ((min x, max x, length, min y, max y, width, type), count)')
for result in results:
    print('Center: ' + str(our_to_real_space_result(result)))
    print('(C, Alpha) = ' + str(find_optimum_params(result)))