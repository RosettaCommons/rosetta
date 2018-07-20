import pyrosetta.distributed

__all__ = ["movers", "filters"]

class ComponentDoc(object):
    def __init__(self, name, doc):
        self.name = name
        self.doc = doc
        self.__doc__ = doc

    def __str__(self):
        return self.doc

    def __repr__(self):
        return "%r(name=%s)" % (self.__class__, self.name)

    def _repr_pretty_(self, p, cycle):
        p.text(self.doc)


class InlineDocs(object):

    @property
    def _component_names(self):
        raise NotImplementedError

    @pyrosetta.distributed.requires_init
    def __dir__(self):
        return self._component_names

    def __getattr__(self, name):
        if name in self._component_names:
            return self.get_component_doc(name)
        else:
            raise AttributeError("No defined component: %r" % name)

    @pyrosetta.distributed.requires_init
    def get_component_doc(self, name):
        os = pyrosetta.rosetta.std.stringstream()
        pyrosetta.rosetta.protocols.rosetta_scripts.print_information(name, os)
        return ComponentDoc(name, os.bytes().decode())


class MoverDocs(InlineDocs):
    @property
    def _component_names(self):
        import pyrosetta.rosetta.protocols.moves as moves
        return list(moves.MoverFactory.get_instance().mover_creator_map())


class FilterDocs(InlineDocs):
    @property
    def _component_names(self):
        import pyrosetta.rosetta.protocols.filters as filters
        return list(filters.FilterFactory.get_instance().filter_creator_map())

movers = MoverDocs()
filters = FilterDocs()
