import h5py

from pyrosetta.rosetta.protocols.indexed_structure_store import (
    StructureStoreProvider,
    StructureStoreManager,
    StructureStore,
)

from pyrosetta.utility.array import assign_by_field_names


class H5StructureStoreProvider(StructureStoreProvider):
    def can_load(self, store_path):
        return h5py.is_hdf5(store_path)

    def can_write(self, store_path):
        return False

    def load_store(self, store_path):
        with h5py.File(store_path, "r") as indb:
            structures = indb["structures"][:]
            residues = indb["residues"][:]
            self.last_residues = residues.copy()

        result = StructureStore()
        result.resize(len(structures), len(residues))
        result.structure_entries[:] = structures
        assign_by_field_names(result.residue_entries, residues)
        result.residue_entries["structure_id"][:] = residues["id"]
        result.residue_entries["residue_id"][:] = residues["resn"]

        return result


_h5py_provider = H5StructureStoreProvider()
StructureStoreManager.get_instance().register_store_provider(-1, "h5py", _h5py_provider)
