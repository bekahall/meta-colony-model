import matplotlib.pyplot as plt
import math
import random
from classes import *
from error_functions import *

def run_simulation(final_time, traits, size, cell_count, frequency, colony_no = 1, display = False):

    # check that traits argument is valid
    if not isinstance(traits, CellTraits):
        raise ValueError("traits argument must be a CellTraits object")
    
    # check that size argument is valid
    size = check_size(size)

    # check that display argument is valid
    if not isinstance(display, bool):
        raise ValueError("display argument must be a bool")

    [lattice, tracker] = initialise_subcolonies(frequency, cell_count, size, traits, colony_no=colony_no)

    t = 0

    if display:
        plt.ion()
        figure = lattice.return_colony_image()
        plt.show()
        plt.pause(1)
        colony_size = lattice.colony_size()
        prev_size = colony_size

    while (t < final_time):
        [lattice, deltat] = tracker.growth(lattice)

        t = t + deltat

        if display:
            colony_size = lattice.colony_size()
            if (colony_size - prev_size) > 1000:
                A = lattice.return_updated_masked_array()
                figure.set_array(A)
                plt.pause(0.005)
                print('t = {} - size = {}'.format(t, colony_size))
                prev_size = colony_size
    
    if display:
        plt.ioff()
        print('simulation complete...')
        figure = lattice.return_colony_image()
        plt.show()

    return lattice

def initialise_subcolonies(frequency, cell_count, size, traits, colony_no = 1):
    """
    Creates initial subcolonies to populate lattice.

    Parameters:
    frequency (int, float, list): Frequencies of the mutant allele in subcolonies. May be a list of different frequencies for each subcolony or a single value to be used for all subcolonies.
    cell_count (int, list): Number of cells in subcolonies. May be a list of different counts for each subcolony or a single value for all subcolonies.
    size (list): List of length two representing dimensions of the colony lattice.
    traits (CellTraits): Determines traits of the wild type and mutant cells.
    colony_no (int, default = 1): Determines the number of initial subcolonies. If both frequency and cell_count are lists, then this is not considered. If only one or the other is a list, colony_no must be the same length as that list.

    Returns:
    list: Contains a list of SubColonyCell objects which describe the initial subcolonies.
    """

    # get list of frequencies
    frequencies = return_frequency_list(frequency, colony_no)

    # get list of cell_counts for each subcolony
    cell_counts = return_cell_count_list(cell_count, colony_no)

    # check that frequencies and cell_counts are of the same length
    colony_no = get_colony_no(frequencies, cell_counts)

    # check that size argument is valid
    size = check_size(size)

    # check that traits argument is valid
    if not isinstance(traits, CellTraits):
        raise ValueError("traits argument must be a CellTraits object")

    # set coordinates for initial subcolonies
    x_mid = int(size[0]/2)
    y_mid = int(size[1]/2)

    radius = int(math.sqrt(colony_no/math.pi)) + 1
    i_list = list(range(-radius, radius + 1))
    j_list = list(range(-radius, radius + 1))
    pairs = []
    for i in i_list:
        for j in j_list:
            if math.sqrt(i**2 + j**2) <= radius:
                pairs.append([i,j])
    
    chosen_pairs = random.sample(pairs, colony_no)
    coords = []

    for pair in chosen_pairs:
        x = x_mid + pair[0]
        y = y_mid + pair[1]
        coords.append([x,y])

    # create initial subcolonies
    init_colonies = []
    for j in range(colony_no):
        init_colonies.append(SubColonyCell(coords[j], traits, cell_counts[j], frequencies[j]))

    # create lattice
    lattice = ColonyLattice(size, init_colonies)

    # define birth probabilities according to the neighbours
    updated_init_colonies = []
    for colony in init_colonies:
        colony.find_neighbours(lattice)
        updated_init_colonies.append(colony)

    # create tracker with the updated subcolonies
    tracker = ColonyTracker(updated_init_colonies)

    return lattice, tracker
