import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import random
import ast
from helper_ftns import *

class ColonyLattice:
    # Attributes
    
    def __init__(self, size, init_colonies):
        self.size = size
        self.states = np.empty(size)   # create lattice of states
        self.states[:] = np.NaN
        self.counts = np.zeros(size)   # create lattice of total cell counts
        for colony in init_colonies:
            x = colony.coords[0]
            y = colony.coords[1]
            self.states[x,y] = colony.MT_frequency()
            self.counts[x,y] = colony.cell_total()

        self.cmap = plt.get_cmap('Reds', 256)

    # Methods

    def update_state(self, coords, new_state):
        x = coords[0]
        y = coords[1]
        self.states[x,y] = new_state

    def update_count(self, coords, new_count):
        x = coords[0]
        y = coords[1]
        self.counts[x,y] = new_count

    def return_colony_image(self):
        masked_array = np.ma.array(self.states, mask = np.isnan(self.states))
        cmap = self.cmap
        cmap.set_bad('black', 1.)
        norm = mpl.colors.Normalize(vmin=0, vmax=1)
        return plt.imshow(masked_array, cmap=cmap, norm=norm)
    
    def return_updated_masked_array(self):
        return np.ma.array(self.states, mask = np.isnan(self.states))
    
    def return_state(self, coords):
        x = coords[0]
        y = coords[1]

        return self.states[x,y]
    
    def return_count(self, coords):
        x = coords[0]
        y = coords[1]

        return self.counts[x,y]
    
    def colony_size(self):
        return np.sum(self.counts)
    
    def colony_frequency(self):
        return np.nanmean(self.states)

    def change_cmap(self, newcmap):
        self.cmap = newcmap

    def __str__(self):
        return str(self.states) + '\n' + str(self.counts)
    
class SubColonyCell:
    def __init__(self, coords, cell_traits, cell_count, frequency):
        self.coords = coords   # set lattice coordinates 

        self.cell_traits = cell_traits   # set parameters for WT and MT cells
        self.mutationWT = cell_traits['mutationWT']
        self.rWT = cell_traits['rWT']
        self.KWT = cell_traits['KWT']
        self.mWT = cell_traits['mWT']
        self.mutationMT = cell_traits['mutationMT']
        self.rMT = cell_traits['rMT']
        self.KMT = cell_traits['KMT']
        self.mMT = cell_traits['mMT']
        self.max = max(self.KWT, self.KMT)

        self.MT_count = int(frequency * cell_count)   # set cell counts
        self.WT_count = cell_count - self.MT_count

        Nt = self.WT_count + self.MT_count
        self.birthWT = in_birth_rate(self.rWT, self.KWT, Nt, self.WT_count)
        self.birthMT = in_birth_rate(self.rMT, self.KMT, Nt, self.MT_count)

    def find_neighbours(self, states):
        # identifies coordinates of neighbouring sites and associated birth rates

        size = states.size
        neighboursWT = {}
        neighboursMT = {}
        x = self.coords[0]
        y = self.coords[1]
        for i in [0,1,-1]:
            for j in [-1,0,1]:
                if (i != 0) | (j != 0):
                    ncoords = [x + i, y + j]
                    if (x + i >= 0) & (x + i < size[0]) & (y + j >= 0) & (y + j < size[1]):   # check if coords are within lattice
                        Nt = states.return_count(ncoords)
                        WT_birth_rate = out_birth_rate(self.rWT, self.KWT, self.mWT, Nt, self.WT_count)
                        MT_birth_rate = out_birth_rate(self.rMT, self.KMT, self.mMT, Nt, self.MT_count)
                    else:
                        WT_birth_rate = 0
                        MT_birth_rate = 0
                    neighboursWT[str(ncoords)] = WT_birth_rate
                    neighboursMT[str(ncoords)] = MT_birth_rate
        self.neighboursWT = neighboursWT
        self.neighboursMT = neighboursMT

    def return_neighbour_coords(self):
        return list(self.neighboursWT)

    def update_neighbour_birth_rate(self, coords, states):
        # updates birth rates associated with the neighbour at the given coordinates

        Nt = states.return_count(coords)
        self.neighboursWT[str(coords)] = out_birth_rate(self.rWT, self.KWT, self.mWT, Nt, self.WT_count)
        self.neighboursMT[str(coords)] = out_birth_rate(self.rMT, self.KMT, self.mMT, Nt, self.MT_count)

    def update_birth_rate(self):
        # updates rate for birth within the current lattice site

        Nt = self.WT_count + self.MT_count
        self.birthWT = in_birth_rate(self.rWT, self.KWT, Nt, self.WT_count)
        self.birthMT = in_birth_rate(self.rMT, self.KMT, Nt, self.MT_count)

    def return_total_birth(self):
        # returns the sum of all the birth rates

        self.total_birth = self.birthWT + self.birthMT + sum(list(self.neighboursWT.values())) + sum(list(self.neighboursMT.values()))

        return self.total_birth
    
    def update_counts(self, cell_type):
        if cell_type == 'WT':
            self.WT_count = self.WT_count + 1
        elif cell_type == 'MT':
            self.MT_count = self.MT_count + 1

    def cell_total(self):
        return self.WT_count + self.MT_count
    
    def MT_frequency(self):
        return self.MT_count/(self.WT_count + self.MT_count)
    
    def birth_event(self):
        WT_birth_total = (self.birthWT + sum(self.neighboursWT.values()))
        MT_birth_total = (self.birthMT + sum(self.neighboursMT.values()))
        WT_birth_prob = WT_birth_total/self.total_birth

        rand1 = np.random.uniform()
        rand2 = np.random.uniform()
        if rand1 < WT_birth_prob:
            coords = [str(self.coords)] + list(self.neighboursWT.keys())
            p = np.array(([self.birthWT] + list(self.neighboursWT.values())))/WT_birth_total
            destination = np.random.choice(coords, p = p)
            if rand2 < self.mutationWT:
                print('a mutation!')
                cell_type = 'MT'
            else:
                cell_type = 'WT'
        else:
            coords = [str(self.coords)] + list(self.neighboursMT.keys())
            p = np.array(([self.birthMT] + list(self.neighboursMT.values())))/MT_birth_total
            destination = np.random.choice(coords, p = p)
            if rand2 < self.mutationMT:
                print('a back mutation')
                cell_type = 'WT'
            else:
                cell_type = 'MT'

        return destination, cell_type
    
    def __str__(self):
        return 'this is a subcolony at ' + str(self.coords) + ' with allele frequency ' + str(self.MT_frequency()) + ' and cell count ' + str(self.cell_total())
    
class ColonyTracker:
    
    def __init__(self, init_colonies):
        self.epsilon = 16
        self.subcolonies = SubColonyDict(init_colonies)   # dictionary of all existing subcolonies keyed by coordinates
        self.birth_rates = SubColonyBirthsDict(init_colonies, 0)
        # self.high_birth_rates = SubColonyBirthsDict(init_colonies, self.epsilon)   # dictionary of all high subcolony birth rates keyed by coordinates
        # self.low_birth_rates = SubColonyBirthsDict(init_colonies, self.epsilon, above=False)   # dictionary of all low subcolony birth rates keyed by coordinates

    def growth(self, states):
        # high_rates = list(self.high_birth_rates.values())
        # low_rates = list(self.low_birth_rates.values())
        # total_high_rate = sum(high_rates)
        # total_low_rate = sum(low_rates)
        # total_rate = total_high_rate + total_low_rate
        # rand1 = np.random.uniform()
        # if rand1 < (total_high_rate/total_rate):
        #     coords = list(self.high_birth_rates.keys())
        #     p = np.array(high_rates)/total_high_rate
        # else:
        #     coords = list(self.low_birth_rates.keys())
        #     p = np.array(low_rates)/total_low_rate
        # chosen_coords = np.random.choice(coords, p = p)
        # chosen_colony = self.subcolonies[chosen_coords]

        coords = list(self.birth_rates.keys())
        temp = list(self.birth_rates.values())
        total_rate = sum(temp)
        p = np.array(temp)/total_rate
        chosen_coords = np.random.choice(coords, p = p)
        chosen_colony = self.subcolonies[chosen_coords]

        [destination, cell_type] = chosen_colony.birth_event()
        coord_list = ast.literal_eval(destination)   # get destination as a list

        # update subcolony experiencing growth
        if destination in self.subcolonies:
            subcolony = self.subcolonies[destination]
            subcolony.update_counts(cell_type)
            subcolony.update_birth_rate()
            total_birth_rate = subcolony.return_total_birth()
            if total_birth_rate > 0:
                self.birth_rates[destination] = total_birth_rate
            elif destination in self.birth_rates:
                del self.birth_rates[destination]

            # if total_birth_rate > self.epsilon:
            #     self.high_birth_rates[destination] = total_birth_rate
            #     if destination in self.low_birth_rates:
            #         del self.low_birth_rates[destination]
            # elif total_birth_rate > 0:
            #     self.low_birth_rates[destination] = total_birth_rate
            #     if destination in self.high_birth_rates:
            #         del self.high_birth_rates[destination]
            # elif destination in self.low_birth_rates:
            #     del self.low_birth_rates[destination]
            # elif destination in self.high_birth_rates:
            #     del self.high_birth_rates[destination]

            self.subcolonies[destination] = subcolony
        else:
            if cell_type == 'WT':
                frequency = 0
            elif cell_type == 'MT':
                frequency = 1
            subcolony = SubColonyCell(coord_list, chosen_colony.cell_traits, 1, frequency)
            subcolony.find_neighbours(states)
            total_birth_rate = subcolony.return_total_birth()
            self.birth_rates[destination] = total_birth_rate
            self.subcolonies[destination] = subcolony

            # if total_birth_rate > self.epsilon:
            #     self.high_birth_rates[destination] = total_birth_rate
            # elif total_birth_rate > 0:
            #     self.low_birth_rates[destination] = total_birth_rate

        # update states lattice
        states.update_state(subcolony.coords, subcolony.MT_frequency())
        states.update_count(subcolony.coords, subcolony.cell_total())
        
        # update neighbouring subcolonies
        neighbour_list = self.subcolonies[destination].return_neighbour_coords()
        for coord in neighbour_list:
            if coord in self.subcolonies:
                colony = self.subcolonies[coord]
                colony.update_neighbour_birth_rate(coord_list, states)
                total_birth_rate = colony.return_total_birth()

                if total_birth_rate > 0:
                    self.birth_rates[coord] = total_birth_rate
                elif coord in self.birth_rates:
                    del self.birth_rates[coord]

                # if total_birth_rate > self.epsilon:
                #     self.high_birth_rates[coord] = total_birth_rate
                #     if coord in self.low_birth_rates:
                #         del self.low_birth_rates[coord]
                # elif total_birth_rate > 0:
                #     self.low_birth_rates[coord] = total_birth_rate
                #     if coord in self.high_birth_rates:
                #         del self.high_birth_rates[coord]
                # elif coord in self.low_birth_rates:
                #     del self.low_birth_rates[coord]
                # elif coord in self.high_birth_rates:
                #     del self.high_birth_rates[coord]

                self.subcolonies[coord] = colony

        deltat = 1/total_rate

        return states, deltat

class CellTraits(dict):

    def __init__(self, rWT, KWT, mWT, mutationWT, rMT, KMT, mMT, mutationMT):
        super().__init__([
            ('rWT', rWT),
            ('KWT', KWT),
            ('mWT', mWT),
            ('mutationWT', mutationWT),
            ('rMT', rMT),
            ('KMT', KMT),
            ('mMT', mMT),
            ('mutationMT', mutationMT)
        ])

class SubColonyDict(dict):
    def __init__(self, init_colonies):
        mapping = {}
        for i in range(len(init_colonies)):
            mapping[str(init_colonies[i].coords)] = init_colonies[i]
        super().__init__(mapping)

class SubColonyBirthsDict(dict):
    def __init__(self, init_colonies, epsilon, above = True):
        mapping = {}
        for i in range(len(init_colonies)):
            rate = init_colonies[i].return_total_birth()
            if above:
                if rate > epsilon:
                    mapping[str(init_colonies[i].coords)] = rate
            else:
                if rate <= epsilon:
                    mapping[str(init_colonies[i].coords)] = rate
        super().__init__(mapping)

#####################


# class SubColonyCell:
#     # Attributes

#     def __init__(self, coords, cell_traits, cell_count, frequency):
#         self.coords = coords

#         self.cell_traits = cell_traits
#         self.mutationWT = cell_traits['mutationWT']
#         self.rWT = cell_traits['rWT']
#         self.KWT = cell_traits['KWT']
#         self.mWT = cell_traits['mWT']
#         self.mutationMT = cell_traits['mutationMT']
#         self.rMT = cell_traits['rMT']
#         self.KMT = cell_traits['KMT']
#         self.mMT = cell_traits['mMT']
#         self.max = max(self.KWT, self.KMT)

#         self.MT_count = int(frequency * cell_count)
#         self.WT_count = cell_count - self.MT_count

#         self.growing = True
#         self.migrating = True
#         self.neighbours = []

#     def update_neighbours(self, states):
#         # identifies coordinates of neighbouring sites which do not contain a full subcolony
#         size = states.size
#         neighbours = []
#         x = self.coords[0]
#         y = self.coords[1]
#         for i in [0,1,-1]:
#             for j in [-1,0,1]:
#                 if (i != 0) | (j != 0):
#                     ncoords = [x + i, y + j]
#                     if (x + i >= 0) & (x + i < size[0]) & (y + j >= 0) & (y + j < size[1]):    # check if coords are within the lattice
#                         if states.return_count(ncoords) < self.max:
#                             neighbours.append(ncoords)
#         self.neighbours = neighbours
    
#     def update_growing(self):
#         if self.wtbirth() + self.mtbirth() > 0:
#             self.growing = True
#         else:
#             self.growing = False

#     def update_migrating(self):
#         if len(self.neighbours) == 0:
#             self.migrating = False

#     def cell_total(self):
#         return self.WT_count + self.MT_count
    
#     def MT_frequency(self):
#         return self.MT_count/(self.WT_count + self.MT_count)
    
#     def update_counts(self, cell_type):
#         if cell_type == 'WT':
#             self.WT_count = self.WT_count + 1
#         elif cell_type == 'MT':
#             self.MT_count = self.MT_count + 1

#     def wtbirth(self):
#         Nt = self.cell_total()
#         if Nt > self.KWT:
#             return 0
#         else:
#             return (self.rWT * self.WT_count) / (self.max*2)
        
#         # return max(0, self.rWT * (1 - (Nt/self.KWT)) * self.WT_count)/Nt
    
#     def mtbirth(self):
#         Nt = self.cell_total()
#         if Nt > self.KMT:
#             return 0
#         else:
#             return (self.rMT * self.MT_count) / (self.max*2)
#         # return max(0, self.rMT * (1 - (Nt/self.KMT)) * self.MT_count)/Nt

#     def grow(self):
#         Nt = self.cell_total()
#         wtbirth = self.wtbirth()
#         mtbirth = self.mtbirth()
        
#         rand1 = np.random.uniform()
#         if rand1 < wtbirth:
#             rand2 = np.random.uniform()
#             if rand2 < self.mutationWT:
#                 print('a mutation!')
#                 self.MT_count = self.MT_count + 1
#             else:
#                 self.WT_count = self.WT_count + 1
#         elif rand1 < wtbirth + mtbirth:
#             rand2 = np.random.uniform()
#             if rand2 < self.mutationMT:
#                 print(' a back mutation!')
#                 self.WT_count = self.WT_count + 1
#             else:
#                 self.MT_count = self.MT_count + 1
        
#         return self.max

#     def migrate(self):
#         Nt = self.cell_total()
#         neighbour_count = len(self.neighbours)
#         random.shuffle(self.neighbours)
#         wtmigration = (self.rWT * self.WT_count * self.mWT * (neighbour_count/8))/(self.max*2)
#         mtmigration = (self.rMT * self.MT_count * self.mMT * (neighbour_count/8))/(self.max*2)

#         rand1 = np.random.uniform()
#         if rand1 < wtmigration:
#             rand2 = np.random.uniform()
#             if rand2 < self.mutationWT:
#                 cell_type = 'MT'
#                 print('a migration mutation!')
#             else:
#                 cell_type = 'WT'
#             destination = random.choice(self.neighbours)
#             migration = True
#         elif rand1 < wtmigration + mtmigration:
#             rand2 = np.random.uniform()
#             if rand2 < self.mutationMT:
#                 cell_type = 'WT'
#                 print('a migration back mutation!')
#             else:
#                 cell_type = 'MT'
#             destination = random.choice(self.neighbours)
#             migration = True
#         else:
#             cell_type = None
#             destination = None
#             migration = False

#         return migration, destination, cell_type, self.max

#     def __str__(self):
#         return 'this is a subcolony at ' + str(self.coords) + ' with allele frequency ' + str(self.MT_frequency()) + ' and cell count ' + str(self.cell_total())
        

# class ColonyTracker:

#     def __init__(self, init_colonies):
#         self.subcolonies = SubColonyDict(init_colonies)   # list of coordinates for all existing subcolonies
#         self.growing = list(self.subcolonies.keys())   # list of coordinates for growing subcolonies
#         self.migrating = list(self.subcolonies.keys())   # list of coordinates for migrating subcolonies

#     def growth(self, states):
                                    
#         to_delete = []   # will collect keys of the subcolonies no longer growing
#         delta_t = 0   # total time spent in growth stage
#         update_no = 0

#         for k in self.growing:
#             # fetch subcolony
#             subcolony = self.subcolonies[k]
#             # update subcolony
#             Nt = subcolony.grow()
#             subcolony.update_growing()
#             if not subcolony.growing:   # subcolony is no longer growing and can be removed from list of growing subcolonies
#                 to_delete.append(k)
#             # update state lattice
#             states.update_state(subcolony.coords, subcolony.MT_frequency())
#             states.update_count(subcolony.coords, subcolony.cell_total())
#             # update subcolony in list of subcolonies
#             self.subcolonies[k] = subcolony
#             # increase delta_t and update count
#             delta_t = delta_t + (1/Nt)
#             update_no = update_no + 1
    

#         # remove coordinates of subcolonies which are no longer growing
#         if len(to_delete) > 0:
#             for k in to_delete:
#                 self.growing.remove(k)

#         # average time step
#         if update_no > 0:
#             delta_t = delta_t/update_no
#         else:
#             print('dt: ' + str(delta_t))

#         # return updated states lattice
#         return states, delta_t
    
#     def migration(self, states):
#         to_delete = []   # will collect keys of the subcolonies no longer migrating
#         delta_t = 0   # total time spend in migration stage
#         update_no = 0

#         random.shuffle(self.migrating)
#         for k in self.migrating:
#             # fetch subcolony
#             subcolony = self.subcolonies[k]
#             # update subcolony migration status
#             subcolony.update_neighbours(states)
#             subcolony.update_migrating()
#             if subcolony.migrating:
#                 # determine if and where migration occurs
#                 [migration, destination, cell_type, Nt] = subcolony.migrate()
#                 # if migration occurs, update subcolonies
#                 if migration:
#                     dest_key = str(destination)
#                     # check if destination subcolony already exists, if not, create new subcolony
#                     if dest_key in self.subcolonies:
#                         colony = self.subcolonies[dest_key]
#                         colony.update_counts(cell_type)
#                         self.subcolonies[dest_key] = colony
#                         self.growing.append(dest_key)
#                     else:
#                         if cell_type == 'WT':
#                             frequency = 0
#                         elif cell_type == 'MT':
#                             frequency = 1
#                         new_colony = SubColonyCell(destination, subcolony.cell_traits, 1, frequency)
#                         self.subcolonies.update({dest_key : new_colony})
#                         self.growing.append(dest_key)
#                         self.migrating.append(dest_key)
#                     # update state lattice
#                     states.update_state(self.subcolonies[dest_key].coords, self.subcolonies[dest_key].MT_frequency())
#                     states.update_count(self.subcolonies[dest_key].coords, self.subcolonies[dest_key].cell_total())
#                 # update delta_t and migration count
#                 delta_t = delta_t + (1/Nt)
#                 update_no = update_no + 1
#             else:
#                 to_delete.append(k)
            
#         # remove coordinates of subcolonies which are no longer migrating
#         if len(to_delete) > 0:
#             for k in to_delete:
#                 self.migrating.remove(k)

#         # average time step
#         delta_t = delta_t/update_no

#         # return updated states lattice
#         return states, delta_t