#!python3

from itertools import chain
from collections import defaultdict

def load_targets(filename, levels=None, limit=None):
    """Loads the target coordinates from a CSV file.  This function will now infer
       the number of levels represented in the file, but you can opt to load just a
       subset.
        filename: File to open
        levels: Number of levels to load inclusive of the centre,
                else all levels in the file will be loaded.
    """

    all_targets = AllTargets()

    with open(filename, 'r') as coord_fh:

        #This used to read the file N lines at a time but now it just keeps reading until
        #it finds a single number which it assumes must be the centre of a new target.
        targ_lines = None

        #Silly way to get all the lines in a file followed by a blank,
        #without reading the whole file into a list.
        for aline in chain( (x.rstrip() for x in coord_fh), ['']):

            if ',' not in aline:
                #We've either hit EOF or the start of the next record.
                if targ_lines:
                    all_targets.add_target([
                            [int(x) for x in l.split(',')] for l in targ_lines[:levels]
                        ])
                    if limit and len(all_targets) == limit:
                        break

                targ_lines = []

            targ_lines.append(aline)

    return all_targets

class AllTargets:
    """A holder for a bunch of Target objects.  Normally produced by a
       call to load_targets.
    """

    def __init__(self):
        self._target_dict = dict()

        # This will contain [ index : [ target, target, ... ] ]
        # So a natural use for a defaultdict of lists
        self._reverse_lookup = defaultdict(list)

        self.levels = None

    def __len__(self):
        return len(self._target_dict)

    def __iter__(self):
        """Iteration yields a list of target objects"""
        return self._target_dict.values().__iter__()

    def get_target_by_centre(self, centre):

        return self._target_dict[centre]

    def add_target(self, coords):

        new_target = Target(coords)

        #You shouldn't add the same target twice
        assert new_target.get_centre() not in self._target_dict

        #All targets should be the same size
        if self.levels is None:
            self.levels = new_target.get_levels()
        else:
            assert self.levels == new_target.get_levels()

        self._target_dict[new_target.get_centre()] = new_target

        #As well as indexing by centre, hash all indices
        for idx in new_target.get_indices():
                self._reverse_lookup[idx].append(new_target)

    def get_all_indices(self, level=None):
        """Returns all the indices held in all targets.
           Optionally limit to level.
        """
        if level == 0:
            #Do it the quick way
            return list(self._target_dict.keys())
        elif level is None:
            #Use the reverse lookup dict
            return list(self._reverse_lookup.keys())
        else:
            #Scan all targets and flatten the list (standard Python-ism)
            return [ y for x in self for y in x.get_indices(level) ]

    def get_from_index(self, index):

        res = []
        for target in self._reverse_lookup[index]:
            res.append( (target, target.get_level_from_index(index)) )

        return res

class Target:
    def __init__(self, coords):
        """Parse the data which is in the format of
           [[n],[n,n,...],[n,n,n,n,n,...]]
        """
        assert len(coords[0]) == 1, "Centre of target must be a single int, not " + str(coords)
        #centre = coords[0][0]

        self.coords = coords

    def get_indices(self, level = None):

        if level is None:
            #Flatten the list (standard Python-ism)
            return [ y for x in self.coords for y in x ]
        else:
            return self.coords[level]

    def get_centre(self):

        return self.get_indices(0)[0]

    def get_levels(self):
        """Get the target size
        """
        return len(self.coords)

    def get_level_from_index(self, index):

        for lev, arr in enumerate(self.coords):
            if index in arr:
                return lev

        #Do we want this??
        #raise Exception("No such index")
        # or this?
        return None
