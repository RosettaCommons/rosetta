import traitlets

import pyrosetta.distributed.utility.log


class TaskBase(pyrosetta.distributed.utility.log.LoggerMixin, traitlets.HasTraits):
    """Base class for traitlets-based callable, pickle-able task objects.

    Callable-base parallelization frameworks rely on two core properties to define a task:
        - The task must be callable, defining a sensible __call__ attribute.
        - The task must be pickle-able.

    Function objects fulfill these requirements, however parallelization frameworks differ in their
    support of variable closures. Additionally, tasks may require the execution of a relatively costly
    initialization procedure to support rapid execution against a large number of input arguments.

    The TaskBase object is designed to capture the behavior of a functional
    closure, support single-run setup, and defining two types of member attributes:
        state    - pickle-able configuration attributes defining the behavior of
                   the task object.
        run data - attributes required to execute the task that can be
                   generated from configuration state

    Tasks partition between these classes by annotating state as traitlets, and
    run data as arbitrary attributes. Tasks are then pickled by capturing
    all named traites.

    """

    def __repr__(self):
        return "%s(%s)" % (
            self.__class__.__name__.split(".")[-1],
            ", ".join(
                "%s = %r" % (k, getattr(self, k))
                for k in self.trait_names()
            ))

    @property
    def __name__(self):
        return self.__repr__()

    def __getstate__(self):
        return {k: getattr(self, k) for k in self.trait_names()}

    def __setstate__(self, state):
        self.__init__(**state)

    def maybe_setup(self):
        if not getattr(self, "is_setup", False):
            self.logger.debug("Performing setup.")
            self.setup()
            self.is_setup = True

        return self

    def setup(self):
        pass

    def execute(self):
        raise NotImplemented("TaskBase.execute")

    def __call__(self, *args, **kwargs):
        return self.maybe_setup().execute(*args, **kwargs)
