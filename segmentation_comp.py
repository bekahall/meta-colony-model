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

    num_cores = multiprocessing.cpu_count()

    runs = 10

    r = 0.24
    mutation_rate = 0

    Kvals = [30, 60, 90, 120]
    mvals = [0.3, 0.5, 0.7, 0.9]

    total_days = 5
    final_time = total_days * 24
    size = [1000, 1000]
    cell_count = 10
    frequency = 0.5
    colony_no = 1000

    image_folder = Path("segmentation_data/images/")
    data_folder = Path("segmentation_data/data/")
    image_folder.mkdir(parents=True, exist_ok=True)
    data_folder.mkdir(parents=True, exist_ok=True)

    for K in Kvals:
        for m in mvals:
            traits = CellTraits(r, K, m, mutation_rate, r, K, m, mutation_rate)
            run_list = list(range(1, runs + 1))

            if __name__ == "__main__":
                data = Parallel(n_jobs=num_cores)(delayed(run_function)(run, final_time, traits, size, cell_count, frequency, colony_no, image_folder, K, m) for run in run_list)

            df = pd.DataFrame.from_dict(data)
            file_name = 'K' + str(K) + 'm' + str(int(m*10)) + '.csv'
            file_location = data_folder / file_name
            df.to_csv(file_location)

def run_function(run, final_time, traits, size, cell_count, frequency, colony_no, image_folder, K, m):
    # runs to be run in parallel
    lattice = run_simulation(final_time, traits, size, cell_count, frequency, colony_no = colony_no)

    image = lattice.counts
    image[np.nonzero(image)] = 1
    perimeter = measure.perimeter(image, neighborhood=8)
    area = np.count_nonzero(lattice.counts)
    p2a = (perimeter**2)/(4 * math.pi * area)

    figure_name = 'K' + str(K) + 'm' + str(int(m*10)) + 'run' + str(run) + '.png'
    figure_location = image_folder / figure_name
    fig = plt.figure()
    lattice.return_colony_image()
    plt.savefig(figure_location)
    plt.close(fig)

    return {'area' : area, 'perimeter' : perimeter, 'P2A' : p2a}

main()