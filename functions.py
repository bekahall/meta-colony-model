import matplotlib.pyplot as plt
import math
import random
from classes import *
from error_functions import *

def run_simulation(final_time, traits, size, cell_count, frequency, deme_no = 1, display = False):
    """
    Runs a single simulation.

    Parameters:
    final_time (int, float): Maximum time (in hours) that the simulation should run for.
    traits (CellTraits): Contains the traits for wildtype and mutant cells.
    size (int): Size of the colony lattice. List of length two containing non-negative integers.
    cell_count (int, list): Number of cells in demes. May be a list of different counts for each deme or a single value for all demes.
    frequency (int, float, list): Frequencies of the mutant allele in demes. May be a list of different frequencies for each deme or a single value to be used for all demes.
    deme_no (int, default = 1): Determines the number of initial demes. If both frequency and cell_count are lists, then this is not considered. If only one or the other is a list, deme_no must be the same length as that list.
    display (bool): Display colony images or not. Default is False.

    Returns:
    states_list (list): List of tuples which hold the allele frequencies of each deme (NaN if deme is unoccupied) on the lattice at the end of each day of growth.
    counts_list (list): List of tuples which hold the cell counts of each deme on the lattice at the end of each day.
    cmaps_list (list): List of colour maps associated with the display of the colony image at the end of each day.
    norms_list (list): List of norms associated with the colour maps for the display of the colony image at the end of each day.
    time_list (list): List of floats giving the number of hours that have passed for each day.
    """

    # check that traits argument is valid
    if not isinstance(traits, CellTraits):
        raise ValueError("traits argument must be a CellTraits object")
    
    # check that size argument is valid
    size = check_size(size)

    # check that display argument is valid
    if not isinstance(display, bool):
        raise ValueError("display argument must be a bool")

    # initialise simulation
    [lattice, tracker] = initialise_demes(frequency, cell_count, size, traits, deme_no=deme_no)
    
    # save initial lattice conditions
    a = lattice.return_states_lattice()
    b = tuple(map(tuple, a))
    states_list = [b]
    a = lattice.return_counts_lattice()
    b = tuple(map(tuple, a))
    counts_list = [b]
    [a, b] = lattice.return_cmap()
    cmaps_list = [a]
    norms_list = [b]
    time_list = [0]

    t = 0
    prev_t = 0

    if display:
        plt.ion()
        figure = lattice.return_colony_image()
        plt.show()
        plt.pause(1)
        colony_size = lattice.colony_size()
        prev_size = colony_size

    # run simulation
    while (t < final_time):
        [lattice, deltat] = tracker.growth(lattice)

        t = t + deltat

        # save lattice conditions at the end of each day
        if t - prev_t > 24:
            a = lattice.return_states_lattice()
            b = tuple(map(tuple, a))
            states_list.append(b)
            a = lattice.return_counts_lattice()
            b = tuple(map(tuple, a))
            counts_list.append(b)
            [a, b] = lattice.return_cmap()
            cmaps_list.append(a)
            norms_list.append(b)
            time_list.append(t)
            prev_t = t

        if display:
            colony_size = lattice.colony_size()
            if (colony_size - prev_size) > 1000:
                A = lattice.return_updated_masked_array()
                figure.set_array(A)
                plt.pause(0.005)
                print('t = {} - size = {}'.format(t, colony_size))
                prev_size = colony_size
    
    # save final lattice conditions
    a = lattice.return_states_lattice()
    b = tuple(map(tuple, a))
    states_list.append(b)
    a = lattice.return_counts_lattice()
    b = tuple(map(tuple, a))
    counts_list.append(b)
    [a, b] = lattice.return_cmap()
    cmaps_list.append(a)
    norms_list.append(b)
    time_list.append((t+1))
    
    if display:
        plt.ioff()
        print('simulation complete...')
        figure = lattice.return_colony_image()
        plt.show()

    return states_list, counts_list, cmaps_list, norms_list, time_list

def initialise_demes(frequency, cell_count, size, traits, deme_no = 1):
    """
    Creates initial demes to populate lattice.

    Parameters:
    frequency (int, float, list): Frequencies of the mutant allele in demes. May be a list of different frequencies for each deme or a single value to be used for all demes.
    cell_count (int, list): Number of cells in demes. May be a list of different counts for each deme or a single value for all demes.
    size (list): List of length two representing dimensions of the colony lattice.
    traits (CellTraits): Determines traits of the wild type and mutant cells.
    deme_no (int, default = 1): Determines the number of initial demes. If both frequency and cell_count are lists, then this is not considered. If only one or the other is a list, deme_no must be the same length as that list.

    Returns:
    lattice (ColonyLattice): Holds information about the state of the colony lattice.
    tracker (ColonyTracker): Holds information about the state of the occupied demes and associated birth rates.
    """

    # get list of frequencies
    frequencies = return_frequency_list(frequency, deme_no)

    # get list of cell_counts for each deme
    cell_counts = return_cell_count_list(cell_count, deme_no)

    # check that frequencies and cell_counts are of the same length
    deme_no = get_deme_no(frequencies, cell_counts)

    # check that size argument is valid
    size = check_size(size)

    # check that traits argument is valid
    if not isinstance(traits, CellTraits):
        raise ValueError("traits argument must be a CellTraits object")

    # set coordinates for initial demes
    x_mid = int(size[0]/2)
    y_mid = int(size[1]/2)

    radius = int(math.sqrt(deme_no/math.pi)) + 1
    i_list = list(range(-radius, radius + 1))
    j_list = list(range(-radius, radius + 1))
    pairs = []
    for i in i_list:
        for j in j_list:
            if math.sqrt(i**2 + j**2) <= radius:
                pairs.append([i,j])
    
    chosen_pairs = random.sample(pairs, deme_no)
    coords = []

    for pair in chosen_pairs:
        x = x_mid + pair[0]
        y = y_mid + pair[1]
        coords.append([x,y])

    # create initial demes
    init_demes = []
    for j in range(deme_no):
        init_demes.append(Deme(coords[j], traits, cell_counts[j], frequencies[j]))

    # create lattice
    lattice = ColonyLattice(size, init_demes)

    # define birth probabilities according to the neighbours
    updated_init_demes = []
    for deme in init_demes:
        deme.find_neighbours(lattice)
        updated_init_demes.append(deme)

    # create tracker with the updated demes
    tracker = ColonyTracker(updated_init_demes)

    return lattice, tracker
