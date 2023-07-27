import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing
from joblib import Parallel, delayed
from pathlib import Path
from skimage import measure
from classes import *
from functions import *

def main():

    num_cores = 40

    r = 0.24
    mutation_rate = 4.104e-06
    KWT = 30
    mWT = 0.9

    run_list = [[1,90,0.5,0],[2,90,0.5,0],[3,90,0.5,0],[4,90,0.5,0],[5,90,0.5,0],[1,90,0.5,1],[2,90,0.5,1],[3,90,0.5,1],[4,90,0.5,1],[5,90,0.5,1],
            [1,90,0.7,0],[2,90,0.7,0],[3,90,0.7,0],[4,90,0.7,0],[5,90,0.7,0],[1,90,0.7,1],[2,90,0.7,1],[3,90,0.7,1],[4,90,0.7,1],[5,90,0.7,1],
            [1,120,0.5,0],[2,120,0.5,0],[3,120,0.5,0],[4,120,0.5,0],[5,120,0.5,0],[1,120,0.5,1],[2,120,0.5,1],[3,120,0.5,1],[4,120,0.5,1],[5,120,0.5,1],
            [1,120,0.7,0],[2,120,0.7,0],[3,120,0.7,0],[4,120,0.7,0],[5,120,0.7,0],[1,120,0.7,1],[2,120,0.7,1],[3,120,0.7,1],[4,120,0.7,1],[5,120,0.7,1]]

    total_days = 30
    final_time = total_days * 24
    size = [1000, 1000]
    cell_count = 10
    colony_no = 1000

    image_folder = Path("competition_data/images/")
    data_folder = Path("competition_data/data/")
    image_folder.mkdir(parents=True, exist_ok=True)
    data_folder.mkdir(parents=True, exist_ok=True)

    if __name__ == "__main__":
        Parallel(n_jobs=num_cores)(delayed(run_function)(pars, r, KWT, mWT, mutation_rate, final_time, size, cell_count, colony_no, image_folder, data_folder) for pars in run_list)

def run_function(pars, r, KWT, mWT, mutation_rate, final_time, size, cell_count, colony_no, image_folder, data_folder):
    # runs to be run in parallel
    run = pars[0]
    KMT = pars[1]
    mMT = pars[2]
    p0 = pars[3]
    traits = CellTraits(r, KWT, mWT, mutation_rate, r, KMT, mMT, mutation_rate)
    [states_list, counts_list, cmaps_list, norms_list, time_list] = run_simulation(final_time, traits, size, cell_count, p0, colony_no = colony_no)

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

        figure_name = 'KMT' + str(KMT) + 'mMT' + str(int(mMT*10)) + 'p0' + str(int(p0*10)) + 'run' + str(run) + 'day' + str(int(day)) + '.png'
        figure_location = image_folder / figure_name
        fig = plt.figure()
        plt.imshow(masked_array, cmap=cmap, norm=norm)
        plt.savefig(figure_location, dpi=600)
        plt.close(fig)

    data = {'area' : area_list, 'perimeter' : perimeter_list, 'P2A' : p2a_list, 'frequency' : freq_list, 'counts' : cell_count_list}
    df = pd.DataFrame.from_dict(data)
    file_name = 'KMT' + str(KMT) + 'mMT' + str(int(mMT*10)) + 'p0' + str(int(p0*10)) + 'run' + str(run) + '.csv'
    file_location = data_folder / file_name
    df.to_csv(file_location)

main()