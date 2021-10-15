from .h5_importer import h5py, requires_h5py, assign_by_field_names

from pyrosetta.rosetta.protocols.indexed_structure_store import (
    StructureStoreProvider,
    StructureStoreManager,
)


class H5StructureStoreProvider(StructureStoreProvider):

    instance = None

    @requires_h5py
    def can_load(self, store_path):
        return h5py.is_hdf5(store_path)

    def can_write(self, store_path):
        return False

    @requires_h5py
    def load_store(self, store_path):
        with h5py.File(store_path, "r") as indb:
            structures = indb["structures"][:]
            residues = indb["residues"][:]
            self.last_residues = residues.copy()

        from pyrosetta.rosetta.protocols.indexed_structure_store import StructureStore
        result = StructureStore()
        result.resize(len(structures), len(residues))
        result.structure_entries[:] = structures
        assign_by_field_names(result.residue_entries, residues)
        result.residue_entries["structure_id"][:] = residues["id"]
        result.residue_entries["residue_id"][:] = residues["resn"]
        result.residue_entries["ss"][:] = residues["ss"]
        result.residue_entries["cen"][:] = residues["cen"]

        return result


#  called in pyrosetta/__init__.py
def init_H5StructureStoreProvider():
    if H5StructureStoreProvider.instance is None:
        H5StructureStoreProvider.instance = H5StructureStoreProvider()
        StructureStoreManager.get_instance().register_store_provider(
            -1, "h5py", H5StructureStoreProvider.instance
        )
