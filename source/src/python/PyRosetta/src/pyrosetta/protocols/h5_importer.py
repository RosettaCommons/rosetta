import logging

try:
    import h5py
    from pyrosetta.utility.array import assign_by_field_names

    h5py_is_present = True
except ImportError:
    h5py = None
    assign_by_field_names = None
    h5py_is_present = False


def requires_h5py(func):
    def hard_stop(*args, **kwargs):
        import sys

        logging.error("The `h5py` module is missing. Please install it and try again.")
        sys.exit(1)

    if not h5py_is_present:
        return hard_stop
    return func

