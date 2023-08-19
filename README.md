# Metacolony Model
This repo contains the code used to perform the simulations presented in Bandara et al. (202_) (insert DOI). 

The metacolony model implements a modified stochastic stepping-stone model to simulate growth of *Saccharomyces cerevisiae* colonies, using different intrinsic growth rates, carrying capacities, and expansion rates to mimic different phenotypic growth behaviours. A complete description of the model can be found in Bandara et al. (202_). 

## Requirements

This code requires the following packages:

    matplotlib
    pandas
    numpy
    scikit-image

Version details are found in the `requirements.txt` file.

## Usage

A single simulation can be perfomed by running the file `main.py` from the terminal:

    main.py run_no KWT KMT mWT mMT p0 total_days size

`K_` and `m_` assign the carrying capacities and expansion rates of the wildtype (WT) and mutant (MT). `p0` assigns the initial frequency of the mutant, `total_days` is the amount of time the simulation with run for, and `size` is a single integer assigning the size of the square lattice. 

`main.py` creates a directory named KWT`KWT`_KMT`KMT`_mWT`mWT`_mMT`mMT`_p0`p0`, which contains a directory with the run number (set by `run_no`). This contains two subdirectories: data, which contains a CSV file `results.csv` and images, which contains images of the colony every twenty-four hours, named by the number of hours which have passed. 

The mutation rate and intrinsic growth rate (shared by both the wildtype and mutant cells) are 0 and 0.24 by default and can be changed in `main.py`. By default, the colony lattice is seed with 1000 demes each containing 10 cells. This may also be edited within `main.py`.

The role of `main.py` is to run the function `run_simulation` and process and output the resulting data. To run simulations with more customisation, the user may use the function `run_simulation` directly. More detailed documentation for all functions can be found below.

## Files

### `classes.py`

`classes.py` contains the following user-defined classes:

- `ColonyLattice` manages two attributes: an array containing the allele frequencies of each deme in the colony (where unoccupied demes are indicated by NaN) and an array containing the cell counts in each deme of the colony. This class contains methods to update a single deme and methods to return information about the state of the colony or individual demes.

- `Deme` is an object describing the state and actions of a single deme. It holds the coordinates of the given deme, the traits of the wildtype and mutant cells in the deme, the cell counts of each cell type found in the deme, information about its neighbouring demes, and the birth rates all reproductive events which may occur in the deme, namely the rates at which wildtype and mutant cells will bud new cells either within their own deme or into one of their neighbours. This class contains methods to update the state of the deme, calculate birth rates and choose random events to occur within the deme, and return information about the deme.

- `ColonyTracker` manages two attributes: a dictionary of `Deme` objects representing all currently occupied demes, keyed by their coordinates, and a dictionary of the non-zero birth rates associated with these demes, also keyed by their coordinates. This class has a single method: to randomly select a colony birth event and update the status of the demes, birth rates, and colony lattice accordingly.

- `CellTraits` is built on the dictionary class and serves to store the life-history and dispersal traits as well as mutation rates for wildtype and mutant cells.

- `DemeDict` is built on the dictionary class and serves to store the `Deme` objects for `ColonyTracker`.

- `DemeBirthsDict` is built on the dictionary class and serves the non-zero birth rates throughout the colony for `ColonyTracker`.

### `functions.py`

`functions.py` contains the following functions to run the simulation:

- `run_simulation` runs a single simulation. It takes parameters defining the amount of time to run the simulation for in simulated hours (`final_time`), a `CellTraits` object assigning traits and mutation rates for the wildtype and mutant cells (`traits`), and a list of length two assigning the size of the lattice (`size`). It also takes parameters in several forms to seed the lattice:
    - To set the number of cells in each deme, `cell_count` can take the value of a single integer assigning the same number of cells to each deme, or a list of integers where each integer is number of cells in a deme.
    - To set the initial frequency in each deme, `frequency` can take the value of a single float assignming the same frequency to each deme, or a list of floats where each float is the frequency in a deme. If both `cell_count` and `frequency` are lists, they must be the same length and the order of both lists corresponds to the same order of the demes, i.e. if `cell_count = [10, 20]` and `frequency = [0.5, 0.75]`, then there will be two seed demes. The first will have 10 cells with 5 mutant cells and 5 wildtype cells, and the second will have 20 cells with 15 mutant cells and 5 wildtype cells.
    - To set the number of demes, `deme_no` assigns the number of demes. `deme_no` has a default value of 1 and will be ignored if either `cell_count` or `frequency` are lists.

- `initialise_demes` is called within `run_simulation` and creates the `ColonyTracker` and `ColonyLattice` objects for the simulation, seeded with a circle of demes. 

### `helper_ftns.py`

`helper_ftns.py` contains the functions for calculating birth rates.

### `error_functions.py`

`error_functions.py` contains functions which will call error messages in certain parameters are not in the appropriate form.
