################################################################################
################################## IMPORTS #####################################
################################################################################
import queue
import decimal
import scipy.optimize as optimize
from geopy.distance import great_circle
import math
import numpy as np


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

# calculate distance of two points in kilometers
def distance(lat1,lon1,lat2,lon2):
    # first move the points to original space
    return great_circle(our_to_real_space_point((lat1, lon1)), our_to_real_space_point((lat2, lon2))).kilometers

# count number of events falls inside the input cell + min/max boundaries
def count_events_minmax(cell, related_events):
    #TODO: make it more efficient number of events
    # falling inside this specific cell
    min_max = [dx0,x0,ZERO,dy0,y0,ZERO]
    res = ZERO
    for event in related_events: # loop over all the events
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
def count_events(cell, related_events):
    #TODO: make it more efficient
    res = ZERO
    for event in related_events: # loop over all the events
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
# cell numbers for an square:
# [1  2  3
#  4  5  6
#  7  8  9]
def find_centers():
    print('Finding centers based on flexible cells algorithm ...')
    bfsq = queue.Queue()# FIFO queue for processing cells
    # first cell which contains the full area
    cell = (x0, x0+dx0, dx0, y0, y0+dy0, dy0, CELL_TYPE_1)
    count = len(events)# all the events are inside the full area
    if count < TF_THR:# if there is not enough events
        return queue.Queue()# return empty queue
    bfsq.put((cell, count))# put the (full area cell, count) into the queue
    while not bfsq.empty():# do it til there is an item in the queue
        popItem = bfsq.get()# get an item from queue
        popCell = popItem[CELL_IDX]# get the cell info

        # if the cell meets the size threshold
        if (popCell[X_LEN] <= SIZE_THR) and (popCell[Y_LEN] <= SIZE_THR):
            bfsq.put(popItem)# put the removed item back
            return bfsq# return final result

        if popCell[TYPE_IDX] == CELL_TYPE_1: #four corners, cells 1,3,7,9
            # move over nine different cells of size half of the original cell
            # move over horizontal cells
            for i in [ZERO, popCell[X_LEN] / FOUR, popCell[X_LEN] / TWO]:
                # move over vertical cells
                for j in [ZERO, popCell[Y_LEN] / FOUR, popCell[Y_LEN] / TWO]:
                    type = CELL_TYPE_1
                    if j == popCell[Y_LEN] / FOUR:# middle vertical cells
                        type = CELL_TYPE_3
                    if i == popCell[X_LEN] / FOUR:# middle horizontal cells
                        type = CELL_TYPE_2
                        # middle vertical/horizontal cells
                        if j == popCell[Y_LEN] / FOUR:
                            type = CELL_TYPE_4
                    # make a cell of half size
                    cell = (popCell[X_STR_IDX] + i,
                            popCell[X_STR_IDX] + i + popCell[X_LEN] / TWO,
                            popCell[X_LEN] / TWO,
                            popCell[Y_STR_IDX] + j,
                            popCell[Y_STR_IDX] + j + popCell[Y_LEN] / TWO,
                            popCell[Y_LEN] / TWO,
                            type)
                    # count number of events inside this cell
                    count, min_max = count_events_minmax(cell, events)
                    # if count meets the threshold, add the cell to queue
                    if count >= TF_THR:
                        bfsq.put((cell, count, min_max))

        elif popCell[TYPE_IDX] == CELL_TYPE_2: #two cells 2,8
            # for type_2 cells, just need to move over middle horizontal
            # cells, others are covered already
            i = popCell[X_LEN] / FOUR
            # move over vertical cells
            for j in [ZERO, popCell[Y_LEN] / FOUR, popCell[Y_LEN] / TWO]:
                # make a cell of half size
                cell = (popCell[X_STR_IDX] + i,
                        popCell[X_STR_IDX] + i + popCell[X_LEN] / TWO,
                        popCell[X_LEN] / TWO,
                        popCell[Y_STR_IDX] + j,
                        popCell[Y_STR_IDX] + j + popCell[Y_LEN] / TWO,
                        popCell[Y_LEN] / TWO,
                        CELL_TYPE_2)
                # count number of events inside this cell
                count, min_max = count_events_minmax(cell, events)
                # if count meets the threshold, add the cell to queue
                if count >= TF_THR:
                    bfsq.put((cell, count, min_max))

        elif popCell[TYPE_IDX] == CELL_TYPE_3: #two cells 4,6
            # for type_3 cells, just need to move over middle vertical
            # cells, others are covered already
            j = popCell[Y_LEN] / FOUR
            # move over horizontal cells
            for i in [ZERO, popCell[X_LEN] / FOUR, popCell[X_LEN] / TWO]:
                # make a cell of half size
                cell = (popCell[X_STR_IDX] + i,
                        popCell[X_STR_IDX] + i + popCell[X_LEN] / TWO,
                        popCell[X_LEN] / TWO,
                        popCell[Y_STR_IDX] + j,
                        popCell[Y_STR_IDX] + j + popCell[Y_LEN] / TWO,
                        popCell[Y_LEN] / TWO,
                        CELL_TYPE_3)
                # count number of events inside this cell
                count, min_max = count_events_minmax(cell, events)
                # if count meets the threshold, add the cell to queue
                if count >= TF_THR:
                    bfsq.put((cell, count, min_max))

        elif popCell[TYPE_IDX] == CELL_TYPE_4: #one cell 5
            # for type_4 cells, just need to move over the middle
            # cell, others are covered already
            i = popCell[X_LEN] / FOUR
            j = popCell[Y_LEN] / FOUR
            # make a cell of half size
            cell = (popCell[X_STR_IDX] + i,
                    popCell[X_STR_IDX] + i + popCell[X_LEN] / TWO,
                    popCell[X_LEN] / TWO,
                    popCell[Y_STR_IDX] + j,
                    popCell[Y_STR_IDX] + j + popCell[Y_LEN] / TWO,
                    popCell[Y_LEN] / TWO,
                    CELL_TYPE_4)
            # count number of events inside this cell
            count, min_max = count_events_minmax(cell, events)
            # if count meets the threshold, add the cell to queue
            if count >= TF_THR:
                bfsq.put((cell, count, min_max))

        # when there is no new cell, the only cell would be the last cell
        if bfsq.empty():
            bfsq.put(popItem)
            return bfsq
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

# find optimum values for alpha based on spatial variation paper
def find_optimum_params(cell_count, related_events):
    center_cell = cell_count[CELL_IDX]
    center_count = cell_count[COUNT_IDX]

    cell_center_x = center_cell[X_STR_IDX] + center_cell[X_LEN] / TWO
    cell_center_y = center_cell[Y_STR_IDX] + center_cell[Y_LEN] / TWO

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
    index = ZERO
    for i in i_indices:
        for j in j_indices:
            new_cell = (i[MAP_STR], i[MAP_END], i[MAP_END] - i[MAP_STR],
                        j[MAP_STR], j[MAP_END], j[MAP_END] - j[MAP_STR])
            new_count = count_events(new_cell, related_events)
            new_center_x = (i[MAP_STR] + i[MAP_END]) / TWO
            new_center_y = (j[MAP_STR] + j[MAP_END]) / TWO
            new_dist = distance(new_center_y, new_center_x, cell_center_y, cell_center_x)
            map.append((str(index), new_count, new_dist))
            index += ONE

    # optimizing function which will be executed repeatedly to optimize x
    # x is alpha
    def f(x):
        result = ZERO
        logs = []
        one_minus_logs = []
        # summation over all the cells
        for mapcell in map:
            cdx = center_count * (mapcell[MAPCELL_DIST_IDX] ** (-x))
            # if there is at least one relevant event within the cell
            if mapcell[MAPCELL_COUNT_IDX] != ZERO:
                logs.append(cdx)
                #result += math.log(cdx)
            else:
                one_minus_logs.append(cdx)
                #result += math.log(ONE - cdx)

        # calculate sum for later normalization
        sum_all = np.sum(logs) + np.sum(one_minus_logs)

        # normalizing
        logs = np.log(np.divide(logs, sum_all))
        one_minus_logs = np.log(ONE - np.divide(one_minus_logs, sum_all))

        # calculating the result
        result += np.sum(logs)
        result += np.sum(one_minus_logs)

        #del logs
        #del one_minus_logs

        return (-ONE) * result

    # optimize alpha
    res = optimize.minimize_scalar(f, bounds=(0, 5), method="bounded")

    #delete temp variables
    del i_indices
    del j_indices
    del map

    # res.x would be our alpha
    return (center_count, res.x)

# optimize parameters by first clustring events for each center
# and then optimizing the parameters based on those events
def optimize_params_by_clustring_events_repeatedly(result):
    print('Optimize parameters repeatedly by clustring events for ' + str(len(result)) + ' center(s)')
    # for storing which events are related to each result
    result_events = [[ONE] * len(events) for _ in range(len(results))]
    # for storing parameters C and alpha
    result_params = [(ZERO, ZERO)] * len(results)
    # for storing previous parameters and compare them with new parameters in each round
    previous_result_params = [(ZERO, ZERO)] * len(results)

    # loop maximum of MAX_OPTIMIZATION_LOOP rounds
    loop = ZERO
    while loop < MAX_OPTIMIZATION_LOOP:
        print("Loop " + str(loop+1) + ' of total ' + str(MAX_OPTIMIZATION_LOOP) + ' loop(s)')
        # do it for each single center
        result_index = ONE
        for ri in range(ZERO, len(results)):
            # get list of events that current result has the maximum probability for them
            print('Find related events for center ' + str(result_index))
            related_events = []
            for ei in range(ZERO, len(events)):
                if result_events[ri][ei] == ONE:
                    related_events.append(events[ei])

            print('Optimizing center ' + str(result_index))
            # find the optimum parameters for related events
            result_params[ri] = find_optimum_params(results[ri], related_events)
            result_index += ONE

        # check if we have any change in parameters of any center
        changed = False
        for ri in range(ZERO, len(results)):
            if previous_result_params[ri][PARAM_CENTER_IDX] != \
                    result_params[ri][PARAM_CENTER_IDX] or \
                    previous_result_params[ri][PARAM_ALPHA_IDX] != \
                    result_params[ri][PARAM_ALPHA_IDX]:
                changed = True

        if changed == False:
            break  # if there is no change, end of optimization
        else:
            # if there is some change, update value of previous parameteres for next round
            for ri in range(ZERO, len(results)):
                previous_result_params[ri] = result_params[ri]

        # update events probability for updated result parameters
        for ei in range(ZERO, len(events)):  # for each event
            max_prob = ZERO  # store maximum probability
            max_index = ZERO  # store index of result with maximum probability
            # calculate the probability for each result and find the max
            for ri in range(ZERO, len(results)):
                # prepare params involved in prob calculation
                C = result_params[ri][PARAM_CENTER_IDX]
                alpha = result_params[ri][PARAM_ALPHA_IDX]
                cell = results[ri][CELL_IDX]
                center_x = (cell[X_STR_IDX] + cell[X_END_IDX]) / TWO
                center_y = (cell[Y_STR_IDX] + cell[Y_END_IDX]) / TWO
                dist = distance(events[ei][LAT_IDX], events[ei][LON_IDX], center_y, center_x)

                # calculate the prob
                prob = C * (dist ** (-alpha))

                # check if we found a more suitable center for this event
                if prob > max_prob:
                    max_prob = prob
                    max_index = ri

                # set all to ZERO
                result_events[ri][ei] = ZERO
            # set the event for related center to one
            result_events[max_index][ei] = ONE

        # increment loops
        loop += ONE
    return result_params

################################################################################
################################# PARAMETERS ###################################
################################################################################
TF_THR = 30 # threshold for tweet frequency
SIZE_THR = 0.8 # threshold for area size of interest
COMPRESS_BOUNDARIES_ACTIVE = False
SINGLE_CENTERED = True
MAX_OPTIMIZATION_LOOP = 100
events_cvs_file = 'COVID_coords_1000.csv'
# events file should be a csv of lat,lon


################################################################################
################################## CONSTANTS ###################################
################################################################################
# x means lon
# y means lat
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

EIGHTEEN_HUNDRED = 180
NINETY = 90
FOUR = 4
TWO = 2
ONE = 1
ZERO = 0
ROUND_DIGITS = 5
AREA_OFFSET = 0.001

MAP_STR = 0
MAP_END = 1

MAPCELL_CELL_IDX = 0
MAPCELL_COUNT_IDX = 1
MAPCELL_DIST_IDX = 2

PARAM_CENTER_IDX = 0
PARAM_ALPHA_IDX = 1

################################################################################
################################### RUN ########################################
################################################################################
print('Tweet frequency threshold = ' + str(TF_THR))
print('Area size of interest threshold = ' + str(SIZE_THR))
print('Compress boundaries? = ' + str(COMPRESS_BOUNDARIES_ACTIVE))
print('Single center (True) or multiple center(False)? = ' + str(SINGLE_CENTERED))
print('Maximum number of loops in optimization = ' + str(MAX_OPTIMIZATION_LOOP))
print('Input event file name = ' + str(events_cvs_file))
print('')

events = read_events()
if not SINGLE_CENTERED:
    results = remove_conflicts(compress_boundaries(find_centers()))
else:
    first_center = queue.Queue()
    first_center.put(find_centers().get())
    results = remove_conflicts(compress_boundaries(first_center))
result_params = optimize_params_by_clustring_events_repeatedly(results)

# pritn some headlines
print('')
print('There are ' + str(len(results)) + ' centers - Calculating optimum parameters for each one:')
print('Center template is ((min x, max x, length, min y, max y, width, type), count)')
print('')
# print all centers and parameters
for ri in range(ZERO, len(results)):
    print('Center: ' + str(our_to_real_space_result(results[ri])))
    print('(C,alpha): ' + str(result_params[ri]))