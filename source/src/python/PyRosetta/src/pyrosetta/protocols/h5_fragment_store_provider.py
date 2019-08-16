from __future__ import print_function

import logging

try:
    import h5py

    h5py_is_present = True
except ImportError:
    h5py_is_present = False

from pyrosetta.rosetta.protocols.indexed_structure_store import (
    FragmentSpecification,
    FragmentStoreProvider,
    FragmentStoreManager,
    FragmentStore,
)
from pyrosetta.rosetta.std import (
    vector_unsigned_long,
    vector_double,
    vector_std_vector_double_std_allocator_double_t,
    vector_std_string,
    vector_numeric_xyzVector_double_t,
)
from pyrosetta.rosetta.numeric import xyzVector_double_t


def requires_h5py(func):
    def hard_stop(*args, **kwargs):
        import sys

        logging.error("The `h5py` module is missing. Please install it and try again.")
        sys.exit(1)

    if not h5py_is_present:
        return hard_stop
    return func


class H5PyFragmentStoreProvider(FragmentStoreProvider):

    """
    str_group = {"aa" : [str,str,str,str, ... ]}
    int64_group = {"ss_bin" : [1,1,1,1,4,7,34,76,34, ... ]}
    real_group = {}
    realVector_group = {
        "phi" : [[9],[9],[9],[9],[9], ... ],
        "psi" : [[9],[9],[9],[9],[9], ... ],
        "omega" : [[9],[9],[9],[9],[9], ... ],
        "cen" : [[9],[9],[9],[9],[9], ... ]
    }
    fragment_coordinates = [xyzVector, xyzVector, xyzVector, xyzVector, ... ]
    """

    instance = None

    def __init__(self, *args, **kwargs):
        FragmentStoreProvider.__init__(self, *args, **kwargs)
        self._target_filename = ""

    @requires_h5py
    def append_to_fragment_store(
        self, fragment_store, store_name, group_field, group_type
    ):
        with h5py.File(self._target_filename, "r") as dfile:
            print(
                "Reopening: ",
                self._target_filename,
                "/fragments/",
                store_name,
                " for group ",
                group_field,
                sep="",
            )
            indb = dfile["fragments"][store_name]
            if group_type == "int64":
                # Make vector of type pyrosetta.rosetta.std.vector_unsigned_long() with values extracted from H5file
                # Put it into the int64_groups map
                numbers = vector_unsigned_long()
                for number in indb[:][group_field]:
                    numbers.append(number)
                fragment_store.int64_groups[group_field] = numbers

            elif group_type == "real":
                numbers = vector_double()
                for number in indb[:][group_field]:
                    numbers.append(number)
                fragment_store.real_groups[group_field] = numbers

            elif group_type == "real_per_residue":
                numbers_vector = vector_std_vector_double_std_allocator_double_t()
                numbers = vector_double()
                for set_of_numbers in indb[:][group_field]:
                    numbers.clear()
                    for number in set_of_numbers:
                        numbers.append(number)
                    numbers_vector.append(numbers)
                fragment_store.realVector_groups[group_field] = numbers_vector

            elif group_type == "char_per_residue":
                strings = vector_std_string()
                for list_of_char in indb[:][group_field]:
                    string = ""
                    for char in list_of_char:
                        string = string + str(char.decode("UTF-8"))
                    strings.append(string)
                fragment_store.string_groups[group_field] = strings

            elif group_type == "five_char_per_residue":
                strings = vector_std_string()
                for list_of_char in indb[:][group_field]:
                    string = ""
                    for char in list_of_char:
                        string = string + str(char.decode("UTF-8"))
                    strings.append(string)
                fragment_store.string_groups[group_field] = strings

            else:
                quit(
                    group_type
                    + " is not a valid entry in the fragment store. Currently only int64,real,char and 5char are implemented"
                )

    @requires_h5py
    def get_fragment_store(self, store_name):
        fragment_atoms = vector_std_string()
        fragment_length = 0
        fragment_coords = vector_numeric_xyzVector_double_t()

        # Open the HDF5 File
        with h5py.File(self._target_filename, "r") as dfile:
            print("Opened file: ", self._target_filename)
            # extract the relevent dataset and attributes
            indb = dfile["fragments"][store_name]
            for atom in indb.attrs["fragment_atoms"].decode("UTF-8").split(","):
                fragment_atoms.append(str(atom))
            fragment_length = int(indb.attrs["fragment_length"])

            # Iterate over all fragments and append all coordinates
            # to fragment_coords to initialize the FragmentStore
            for frag in indb[:]["coordinates"]:
                for coord in frag:
                    xyzVector = xyzVector_double_t(coord[0], coord[1], coord[2])
                    fragment_coords.append(xyzVector)

        print("read coordinates")
        fragment_specification = FragmentSpecification(fragment_length, fragment_atoms)
        result = FragmentStore(
            fragment_specification, int(len(fragment_coords) / fragment_length)
        )
        result.fragment_coordinates = fragment_coords
        print("put cooordinates into FragmentStore")

        return result

    def set_target_filename(self, filename):
        self._target_filename = filename


def init_H5FragmentStoreProvider():
    if H5PyFragmentStoreProvider.instance is None:
        H5PyFragmentStoreProvider.instance = H5PyFragmentStoreProvider()
        FragmentStoreManager.get_instance().register_store_provider(
            -1, "h5py", H5PyFragmentStoreProvider.instance
        )
