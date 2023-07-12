import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from skimage import measure
from classes import *
from functions import *

def main():

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
            data = []
            for run in range(1, runs + 1):
                lattice = run_simulation(final_time, traits, size, cell_count, frequency, colony_no = colony_no)

                image = lattice.counts
                image[np.nonzero(image)] = 1
                perimeter = measure.perimeter(image, neighborhood=8)
                area = np.count_nonzero(lattice.counts)
                p2a = (perimeter**2)/(4 * math.pi * area)

                data.append({'area' : area, 'perimeter' : perimeter, 'P2A' : p2a})

                figure_name = 'K' + str(K) + 'm' + str(int(m*10)) + 'run' + str(run) + '.png'
                figure_location = image_folder / figure_name
                fig = plt.figure()
                lattice.return_colony_image()
                plt.savefig(figure_location)
                plt.close(fig)

            df = pd.DataFrame.from_dict(data)
            file_name = 'K' + str(K) + 'm' + str(int(m*10)) + '.csv'
            file_location = data_folder / file_name
            df.to_csv(file_location)

main()