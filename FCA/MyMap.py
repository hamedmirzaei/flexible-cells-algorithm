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
    return ((round(top_left[LON_IDX], ROUND_DIGITS),
             round(bot_right[LON_IDX], ROUND_DIGITS),
             round(our_cell[X_LEN], ROUND_DIGITS),
             round(top_left[LAT_IDX], ROUND_DIGITS),
             round(bot_right[LAT_IDX], ROUND_DIGITS),
             round(our_cell[Y_LEN], ROUND_DIGITS),
             our_cell[TYPE_IDX]),
            our_result[COUNT_IDX])

# calculate distance of two points in kilometers
def distance(lat1,lon1,lat2,lon2):
    # first move the points to original space
    return great_circle(our_to_real_space_point((lat1, lon1)), our_to_real_space_point((lat2, lon2))).kilometers

#storing the start time
def start():
    global start_time
    start_time.append(time.time())

#calculate time length and print it
def stop(location):
    global start_time
    if location == None:
        start_time.pop()
    else:
        print(location + ': done in ' + str(time.time() - start_time.pop()) + ' seconds')

#calculate time length and print it, then reset start_time
def stopAndStart(location):
    stop(location)
    start()

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
    events = []
    res = ZERO
    for event in related_events: # loop over all the events
        if is_event_in_cell(event, cell):
            res += ONE
            events.append(event)
    return res, events

# compress cell by looking over all possible squares with
# size equal to 1.2 * SIZE_THR (0.1 margin in each side)
def compressor(cell_count_events):
    cell = cell_count_events[CELL_IDX]
    events = cell_count_events[EVENTS_IDX]
    # to divide cell into smaller cells of 1/10 of cell size
    small_cell_size = round(SIZE_THR / TEN, ROUND_DIGITS)
    # number of events in each small cell => 2D array
    small_counts = []

    i = cell[X_STR_IDX]
    while i < cell[X_END_IDX]:
        i_end = i + small_cell_size
        if ((i + small_cell_size) > cell[X_END_IDX]):
            i_end = cell[X_END_IDX]

        # each array is a one row of the total small cells' event counts
        array = []
        j = cell[Y_STR_IDX]
        while j < cell[Y_END_IDX]:
            j_end = j + small_cell_size
            if ((j + small_cell_size) > cell[Y_END_IDX]):
                j_end = cell[Y_END_IDX]

            small_cell = (i, i_end, small_cell_size, j, j_end, small_cell_size)
            small_count,_ = count_events(small_cell, events)
            array.append(small_count)
            j += small_cell_size
        small_counts.append(array)
        i += small_cell_size

    # for storing maximum possible compressed cell of size 1.2 * SIZE_THR
    max_i = ZERO
    max_j = ZERO
    max_sum = ZERO
    # by i and j we find the top-left corner of cell,
    # by k and l we loop through small cells
    num_of_cells = int(SIZE_THR / small_cell_size + TWO)
    i = ZERO
    while i <= len(small_counts[ZERO]) - num_of_cells:
        j = ZERO
        while j <= len(small_counts[ZERO]) - num_of_cells:
            sum = ZERO
            for k in range(num_of_cells):
                for l in range(num_of_cells):
                    sum += small_counts[i + k][j + l]
            if sum > max_sum:
                max_sum = sum
                max_i = i
                max_j = j
            j += ONE
        i += ONE
    if max_sum >= TF_THR:
        max_cell_count = ((cell[X_STR_IDX] + max_i * small_cell_size,
                           cell[X_STR_IDX] + (max_i + num_of_cells) * small_cell_size,
                           TWELWE * small_cell_size,
                           cell[Y_STR_IDX] + max_j * small_cell_size,
                           cell[Y_STR_IDX] + (max_j + num_of_cells) * small_cell_size,
                           TWELWE * small_cell_size,
                           cell[TYPE_IDX]), max_sum)
        return max_cell_count
    return None

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
    global events
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

        # if the cell meets the size threshold,
        # we can conclude that all the others are ok
        if (popCell[X_LEN] >= TWO * SIZE_THR) and (popCell[Y_LEN] >= TWO * SIZE_THR) and \
           (popCell[X_LEN] <= FOUR * SIZE_THR) and (popCell[Y_LEN] <= FOUR * SIZE_THR):
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
                    count, cell_events = count_events(cell, events)
                    # if count meets the threshold, add the cell to queue
                    if count >= TF_THR:
                        bfsq.put((cell, count, cell_events))

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
                count, cell_events = count_events(cell, events)
                # if count meets the threshold, add the cell to queue
                if count >= TF_THR:
                    bfsq.put((cell, count, cell_events))

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
                count, cell_events = count_events(cell, events)
                # if count meets the threshold, add the cell to queue
                if count >= TF_THR:
                    bfsq.put((cell, count, cell_events))

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
            count, cell_events = count_events(cell, events)
            # if count meets the threshold, add the cell to queue
            if count >= TF_THR:
                bfsq.put((cell, count, cell_events))

        # when there is no new cell, the only cell would be the last cell
        #if bfsq.empty():
        #    bfsq.put(popItem)
        #    return bfsq
    return bfsq

# prune the results by removing conflict cells from find_centers output
def remove_duplicates(results):
    print('Remove duplicate centers ...')
    pruned_res = []
    for result in results:
        duplicate = False
        for res in pruned_res:
            if abs(res[CELL_IDX][X_STR_IDX] - result[CELL_IDX][X_STR_IDX]) < (SIZE_THR) \
                    and abs(res[CELL_IDX][X_END_IDX] - result[CELL_IDX][X_END_IDX]) < (SIZE_THR) \
                    and abs(res[CELL_IDX][Y_STR_IDX] - result[CELL_IDX][Y_STR_IDX]) < (SIZE_THR) \
                    and abs(res[CELL_IDX][Y_END_IDX] - result[CELL_IDX][Y_END_IDX]) < (SIZE_THR):
                duplicate = True
                break
        if not duplicate:
            pruned_res.append(result)
        else:
            if res[COUNT_IDX] < result[COUNT_IDX]:
                pruned_res.remove(res)
                pruned_res.append(result)
    return pruned_res

# find optimum values for alpha based on spatial variation paper
def find_optimum_params(cell_count, related_events):
    global results_celldists
    global results_eventcells
    start()

    center_cell = cell_count[CELL_IDX]
    center_count = cell_count[COUNT_IDX]

    cell_center_x = center_cell[X_STR_IDX] + center_cell[X_LEN] / TWO
    cell_center_y = center_cell[Y_STR_IDX] + center_cell[Y_LEN] / TWO

    # first we should divide map into cells of size result cell
    # (which is the input of the method)
    result_key = str(cell_center_x) + ',' + str(cell_center_y)
    stored_result_celldists = results_celldists.get(result_key)
    stored_result_eventcells = results_eventcells.get(result_key)
    result_celldists = {}
    result_eventcells = {}
    if not stored_result_celldists:
        i_indices = []  # to store all possible horizontal start/end of cells
        i_indices.append((center_cell[X_STR_IDX],center_cell[X_END_IDX]))
        j_indices = []  # to store all possible vertical start/end of cells
        j_indices.append((center_cell[Y_STR_IDX], center_cell[Y_END_IDX]))
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
        sorted(i_indices, key=lambda x: x[MAP_STR])

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
        sorted(j_indices, key=lambda y: y[MAP_STR])

        for i in i_indices:
            for j in j_indices:
                new_center_x = (i[MAP_STR] + i[MAP_END]) / TWO
                new_center_y = (j[MAP_STR] + j[MAP_END]) / TWO
                new_dist = distance(new_center_y, new_center_x, cell_center_y, cell_center_x)  # 7 s
                cell_key = str(i[MAP_STR]) + '_' + str(j[MAP_STR])
                result_celldists[cell_key] = new_dist
        results_celldists[result_key] = result_celldists

        for related_event in related_events:
            event_key = str(related_event[LAT_IDX]) + '_' + str(related_event[LON_IDX])
            if not result_eventcells.get(event_key):
                located = False
                for i in i_indices:
                    for j in j_indices:
                        new_cell = (i[MAP_STR], i[MAP_END], i[MAP_END] - i[MAP_STR],
                                    j[MAP_STR], j[MAP_END], j[MAP_END] - j[MAP_STR])
                        cell_key = str(i[MAP_STR]) + '_' + str(j[MAP_STR])
                        if is_event_in_cell(related_event, new_cell):
                            result_eventcells[event_key] = cell_key
                            located = True
                            break
                    if located:
                        break
        results_eventcells[result_key] = result_eventcells

    else:
        result_celldists = stored_result_celldists
        result_eventcells = stored_result_eventcells

    stopAndStart('In Optimization - Make, Sort and Dist')
    # move over all the cells and count events for each one
    result_cellcounts = {}
    for cell_key in result_celldists.keys():
        result_cellcounts[cell_key] = ZERO
    for related_event in related_events:
        event_key = str(related_event[LAT_IDX]) + '_' + str(related_event[LON_IDX])
        cell_key = result_eventcells[event_key]
        result_cellcounts[cell_key] = result_cellcounts[cell_key] + 1

    stopAndStart('In Optimization - Count')
    # optimizing function which will be executed repeatedly to optimize x
    # x is alpha
    def f(x):
        result = ZERO
        logs = []
        one_minus_logs = []
        # summation over all the cells
        for cell_key in result_celldists.keys():
            if result_celldists[cell_key] == ZERO:
                continue
            cdx = center_count * (result_celldists[cell_key] ** (-x))
            # if there is at least one relevant event within the cell
            if result_cellcounts[cell_key] != ZERO:
                logs.append(cdx)
                #result += math.log(cdx)
            else:
                one_minus_logs.append(cdx)
                #result += math.log(ONE - cdx)

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

        #del logs
        #del one_minus_logs

        return (-ONE) * result

    # optimize alpha
    res = optimize.minimize_scalar(f, bounds=(0, 5), method="bounded")
    stop('In Optimization - Min')

    # res.x would be our alpha
    return (center_count, res.x)

# optimize parameters by first clustring events for each center
# and then optimizing the parameters based on those events
def optimize_params_by_clustring_events_repeatedly():
    global results
    global events
    print('Optimize parameters repeatedly by clustring events for ' + str(len(results)) + ' center(s)')
    # for storing which events are related to each result
    result_events = [[ONE] * len(events) for _ in range(len(results))]
    # for storing parameters C and alpha
    result_params = [(ZERO, ZERO)] * len(results)
    # for storing previous parameters and compare them with new parameters in each round
    previous_result_params = [(ZERO, ZERO)] * len(results)

    # loop maximum of MAX_OPTIMIZATION_LOOP rounds
    loop = ZERO
    while loop < MAX_OPTIMIZATION_LOOP:
        start()
        print('========================================')
        print("Loop " + str(loop+1) + ' of total ' + str(MAX_OPTIMIZATION_LOOP) + ' loop(s)')
        # do it for each single center
        result_index = ONE
        for ri in range(ZERO, len(results)):
            start()
            # get list of events that current result has the maximum probability for them
            related_events = []
            for ei in range(ZERO, len(events)):
                if result_events[ri][ei] == ONE:
                    related_events.append(events[ei])
            stop('Group up related events for center ' + str(result_index))

            print('Optimizing center ' + str(result_index))
            # find the optimum parameters for related events
            result_params[ri] = find_optimum_params(results[ri], related_events)
            result_index += ONE

        stopAndStart('Overall optimization for loop ' + str(loop))
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

        stop('Updating event probabilities and clustering them for next round')
        # increment loops
        loop += ONE

    return result_params

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
EVENTS_IDX = 2

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
TWELWE = 12
TEN = 10
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
################################# PARAMETERS ###################################
################################################################################
TF_THR = 430 # threshold for tweet frequency
SIZE_THR = 3 # threshold for area size of interest
SINGLE_CENTERED = False
MAX_OPTIMIZATION_LOOP = 20
events_cvs_file = 'map.csv'
# events file should be a csv of lat,lon


################################################################################
################################### RUN ########################################
################################################################################
print('Tweet frequency threshold = ' + str(TF_THR))
print('Area size of interest threshold = ' + str(SIZE_THR))
print('Single center (True) or multiple center(False)? = ' + str(SINGLE_CENTERED))
print('Maximum number of loops in optimization = ' + str(MAX_OPTIMIZATION_LOOP))
print('Input event file name = ' + str(events_cvs_file))
print('')

start_time = []
events = read_events()
founded_centers = queue.Queue()
results = []
results_celldists = {}
results_eventcells = {}

if not SINGLE_CENTERED:
    founded_centers = find_centers()
else:
    temp_centers = find_centers()
    if not temp_centers.empty():
        founded_centers.put(temp_centers.get())

while not founded_centers.empty():
    res = founded_centers.get()
    result = compressor(res)
    if result != None:
        results.append(result)

results = remove_duplicates(results)

results_params = optimize_params_by_clustring_events_repeatedly()

# pritn some headlines
print('')
print('There are ' + str(len(results)) + ' centers - Calculating optimum parameters for each one:')
print('Center template is ((min x, max x, length, min y, max y, width, type), count)')
print('')
# print all centers and parameters
for ri in range(ZERO, len(results)):
    print('Center: ' + str(our_to_real_space_result(results[ri])))
    print('(C,alpha): ' + str(results_params[ri]))