#!python
from __future__ import print_function, division, absolute_import

def load_targets(filename, levels):

    all_targets = AllTargets()

    with open(filename, 'r') as coord_fh:

        #There are various smart-arse ways to read a file N lines
        #at a time, but I'll just do it the obvious dumb way.
        while True:

            targ_lines = []
            for l in range(levels):
                aline = coord_fh.readline().rstrip()

                if aline == '' and len(targ_lines) == 0:
                    #EOF - we're done
                    return all_targets
                elif aline == '':
                    raise Exception("Found blank line or EOF in the middle of a record")

                targ_lines.append(aline)

            all_targets.add_target([ map(int,x.split(',')) for x in targ_lines])

    #This function returns an all_targets object, but it only exits in the middle!


class AllTargets:

    def __init__(self):
        self._target_dict = dict()

        # This will contain [ index : [ target, target, ... ] ]
        self._reverse_lookup = dict()

    def get_all_targets(self):

        return self._target_dict.values()

    def get_target_by_centre(self, centre):

        return self._target_dict[centre]

    def add_target(self, coords):

        new_target = Target(coords)

        #You shouldn't add the same target twice
        assert new_target.get_centre() not in self._target_dict

        self._target_dict[new_target.get_centre()] = new_target

        #As well as indexing by centre, hash all indices
        for idx in new_target.get_indices():
            if idx in self._reverse_lookup:
                self._reverse_lookup[idx].append(new_target)
            else:
                self._reverse_lookup[idx] = [new_target]


    def get_all_indices(self, level = None):

        if level == 0:
            #Do it the quick way
            return list(self._target_dict.keys())
        if level is None:
            #Flatten the list (standard Python-ism)
            return [ y for x in self.get_all_targets() for y in x.get_indices(level) ]

    def get_from_index(self, index):

        res = []
        for target in self._reverse_lookup[index]:
            res.append( (target, target.get_level_from_index(index)) )

        return res

class Target:
    def __init__(self, coords):
        """Parse the data which is in the format of ...
        """
        assert len(coords[0]) == 1
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

    def get_level_from_index(self, index):

        for lev, arr in enumerate(self.coords):
            if index in arr:
                return lev

        #Do we want this??
        #raise Exception("No such index")
        # or this?
        return None
