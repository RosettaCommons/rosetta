#! /usr/bin/python
# List of commands used in PyRosetts Workshop #9
import rosetta
import math
#init()

class MyShape:
    """A base class for a generic shape."""
    def __init__(self):
        self.color = "black" # default value for color

    def __str__(self):
        return self.__doc__

    def area(self):
        """Return the area of the shape."""
        return # Code to calculate the area goes here.

    def draw(self, line_width = 1):
        """Draw the shape on the screen."""
        pass # Code to draw the shape goes here.

class MyCircle(MyShape):
    """A subclass of MyShape for a circle."""
    def __init__(self):
        # This overrides the __init__() method inherited
        # from MyShape.
        MyShape.__init__(self)
        self.radius = 1.0 # default value

    def area(self):
        """Return the area of the circle."""
        # This overrides the area() method inherited from
        # MyShape.
        return math.pi * self.radius**2

class MySquare(MyShape):
    """A subclass of MyShape for a square."""
    def __init__(self):
        # This overrides the __init__() method inherited
        # from MyShape.
        MyShape.__init__(self)
        self.side_length = 1

    def area(self):
        """Return the area of the square."""
        # This overrides the area() method inherited from
        # MyShape.
        return self.side_length**2

class PhiNByXDegreesMover(rosetta.protocols.moves.Mover):
    """A mover that increments the phi angle of residue N
    By X degrees.

    Default values are residue 1 and 15 degrees.

    """
    def __init__(self, N_in=1, X_in=15):
        """Construct PhiNByXDegreesMover."""
        rosetta.protocols.moves.Mover.__init__(self)

        self.N = N_in
        self.X = X_in

    def __str__(self):
        return "residue: " + str(self.N) + \
               " phi increment: " + str(self.X) + \
               " degrees"

    def get_name(self):
        """Return name of class."""
        return self.__class__.__name__

    def apply(self, pose):
        """Applies move to pose."""
        print "Incrementing phi of residue", self.N, "by",
        print self.X, "degrees...."
        pose.set_phi(self.N, pose.phi(self.N) + self.X)

from rosetta.core.scoring.methods import ContextIndependentOneBodyEnergy

@rosetta.EnergyMethod()
class LengthScoreMethod(ContextIndependentOneBodyEnergy):
    """A scoring method that favors longer peptides by
    assigning negative one Rosetta energy unit per
    residue.
    """
    def __init__(self):
        """Construct LengthScoreMethod."""
        ContextIndependentOneBodyEnergy.__init__(self, self.creator())

    def residue_energy(self, res, pose, emap):
        """Calculate energy of res of pose and set emap"""
        # 1 energy unit per residue
        emap.set(self.scoreType, -1.0)

from rosetta.core.scoring.methods import ContextIndependentTwoBodyEnergy

@rosetta.EnergyMethod()
class CI2BScoreMethod(ContextIndependentTwoBodyEnergy):
    """A scoring method that depends on pairs of residues.
    """
    def __init__(self):
        """Construct CI2BScoreMethod."""
        ContextIndependentTwoBodyEnergy.__init__(self,
        self.creator())

    def residue_pair_energy(self, res1, res2, pose, sf, emap):
        """Calculate energy of each pair of res1 and res2
        of pose and set emap."""
        # A real method would calculate a value for score
        # from res1 and res2.
        score = 1.0
        emap.set(self.scoreType, score)

    def atomic_interaction_cutoff(self):
        """Get the cutoff."""
        return 0.0 # Change this value to set the cutoff.

    def defines_intrares_energy(self, weights):
        """Return True if intra-residue energy is
        Defined."""
        return True

    def eval_intrares_energy(self, res, pose, sf, emap):
        """Calculate intra-residue energy if defined."""
        pass
