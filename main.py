import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from skimage import measure
from classes import *
from functions import *

def main():
    # Parameters
    run_no = int(sys.argv[1])
    
    KWT = int(sys.argv[2])
    KMT = int(sys.argv[3])
    mWT = float(sys.argv[4])
    mMT = float(sys.argv[5])
    p0 = float(sys.argv[6])

    r = 0.24
    mutation_rate = 4.104e-06

    total_days = int(sys.argv[7])
    size = [int(sys.argv[8]), int(sys.argv[8])]
    final_time = total_days * 24

    cell_count = 10
    colony_no = 1000

    # Save info
    image_folder = Path("data/KWT{}_KMT{}_mWT{}_mMT{}_p0{}/run{}/images".format(KWT,KMT,int(mWT*10),int(mMT*10),int(p0*10),run_no))
    data_folder = Path("data/KWT{}_KMT{}_mWT{}_mMT{}_p0{}/run{}/data/".format(KWT,KMT,int(mWT*10),int(mMT*10),int(p0*10),run_no))
    image_folder.mkdir(parents=True, exist_ok=True)
    data_folder.mkdir(parents=True, exist_ok=True)

    # Run simulation
    traits = CellTraits(r, KWT, mWT, mutation_rate, r, KMT, mMT, mutation_rate)
    [states_list, counts_list, cmaps_list, norms_list, time_list] = run_simulation(final_time, traits, size, cell_count, p0, colony_no = colony_no)

    # Analysis
    perimeter_list = []
    area_list = []
    p2a_list = []
    freq_list = []
    cell_count_list = []

    for i in range(len(time_list)):
        counts = np.array(counts_list[i])
        states = np.array(states_list[i])
        day = time_list[i]
        cmap = cmaps_list[i]
        norm = norms_list[i]
        masked_array = np.ma.array(states, mask = np.isnan(states))
        size = np.sum(counts)
        area = np.count_nonzero(counts)
        frequency = np.nanmean(states)

        counts[np.nonzero(counts)] = 1
        perimeter = measure.perimeter(counts, neighborhood=8)

        p2a = (perimeter**2)/(4 * math.pi * area)

        perimeter_list.append(perimeter)
        area_list.append(area)
        p2a_list.append(p2a)
        freq_list.append(frequency)
        cell_count_list.append(size)

        figure_name = 'hour{}.png'.format(str(int(day)))
        figure_location = image_folder / figure_name
        fig = plt.figure()
        plt.imshow(masked_array, cmap=cmap, norm=norm)
        plt.savefig(figure_location, dpi=600)
        plt.close(fig)

    data = {'area' : area_list, 'perimeter' : perimeter_list, 'P2A' : p2a_list, 'frequency' : freq_list, 'counts' : cell_count_list}
    df = pd.DataFrame.from_dict(data)
    file_name = 'results.csv'
    file_location = data_folder / file_name
    df.to_csv(file_location)

main()
