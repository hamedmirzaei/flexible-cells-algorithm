################################################################################
################################## IMPORTS #####################################
################################################################################
import queue
import decimal


################################################################################
################################# FUNCTIONS ####################################
################################################################################
# count number of events falls inside the input cell
def count_events(cell):
    #TODO: make it more efficient number of events
    # falling inside this specific cell
    min_max = [dx0,x0,dy0,y0]
    res = 0
    for event in events: # loop over all the events
        lat = event[LAT_IDX]
        lon = event[LON_IDX]
        # check if event falls inside the cell
        if lon >= cell[X_STR_IDX] and lon <= cell[X_END_IDX] and \
                lat >= cell[Y_STR_IDX] and lat <= cell[Y_END_IDX]:
            res += 1

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
        min_max[MIN_X_IDX] = min_max[MIN_X_IDX] - 0.001
        min_max[MAX_X_IDX] = min_max[MAX_X_IDX] + 0.001
    if min_max[MIN_Y_IDX] == min_max[MAX_Y_IDX]:
        min_max[MIN_Y_IDX] = min_max[MIN_Y_IDX] - 0.001
        min_max[MAX_Y_IDX] = min_max[MAX_Y_IDX] + 0.001
    return res, min_max

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
            lat = -float(coordinate[LAT_IDX]) + 90
            lon = float(coordinate[LON_IDX]) + 180
            events.append((lat, lon))# add the (lat, lon) to the events list
    return events

# find all the cells met the provided thresholds and
# return them. It uses a BFS like structure to move
# deeper for cells
def find_centers():
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
                    count, min_max = count_events(cell)
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
                count, min_max = count_events(cell)
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
                count, min_max = count_events(cell)
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
            count, min_max = count_events(cell)
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
    bounded = queue.Queue()
    while not temp_result.empty():
        # get the item from queue
        cell_count_minmax = temp_result.get()

        # get parts of the item
        count = cell_count_minmax[COUNT_IDX]
        minmax = cell_count_minmax[MINMAX_IDX]

        # update cell boundaries
        cell = (minmax[MIN_X_IDX],
                minmax[MAX_X_IDX],
                round(minmax[MAX_X_IDX] - minmax[MIN_X_IDX], 5),
                minmax[MIN_Y_IDX],
                minmax[MAX_Y_IDX],
                round(minmax[MAX_Y_IDX] - minmax[MIN_Y_IDX], 5))

        #put the cell and count in new queue
        bounded.put((cell, count))
    return bounded


################################################################################
################################# PARAMETERS ###################################
################################################################################
TF_THR = 6 # threshold for tweet frequency
SIZE_THR = 0.8 # threshold for area size of interest
x0 = 0 # top left corner considered as (0,0) coordinate
y0 = 0 # top left corner considered as (0,0) coordinate
dx0 = 360 # from -180 to 180
dy0 = 360 # from -90 to 90 (for square cells, we suppose its 360 wide)
events_cvs_file = 'map.csv'


################################################################################
################################## CONSTANTS ###################################
################################################################################
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
MIN_Y_IDX = 2
MAX_Y_IDX = 3

KEY_NAME_1 = 'key1'
KEY_NAME_2 = 'key2'
KEY_NAME_3 = 'key3'
VAL_NAME = 'val'

QUARTER = 4
HALF = 2
ZERO = 0


################################################################################
################################### RUN ########################################
################################################################################
events = read_events()
print('events are: ')
print(events)
results = remove_conflicts(compress_boundaries(find_centers()))
print('results are: ')
for item in results:
    print(item)
