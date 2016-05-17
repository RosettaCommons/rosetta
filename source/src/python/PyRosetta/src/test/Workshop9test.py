#! /usr/bin/python
# List of commands used in PyRosetts Workshop #9

from rosetta import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

from Workshop9_my_shapes import MyCircle, MySquare, PhiNByXDegreesMover, LengthScoreMethod

def hollow(shape_in):
    """Modify the draw() method of an input shape class to
    output a hollow shape."""
    pass # Code to modify the draw() method goes here.

circle = MyCircle()
square = MySquare()
print circle
print square
print circle.color
print square.color
square.side_length = 2
print square.area()
circle.radius, circle.color = 1.5, "pink"
print circle.area()
circle.draw(2)

N_in_2 = 10
X_in_2 = 45
test_pose = pose_from_file("../test/data/workshops/1YY8.clean.pdb")
PhiNByX_1 = PhiNByXDegreesMover()
PhiNByX_2 = PhiNByXDegreesMover()
PhiNByX_2.N = N_in_2
PhiNByX_2.X = X_in_2
PhiNByX_1.apply(test_pose)
PhiNByX_2.apply(test_pose)

circle2 = MyCircle()
circle2.draw() # Draws a filled circle.
#hollow_circle = hollow(circle2)
#hollow_circle.draw() # Draws a hollow circle.

#MyCircle = hollow(MyCircle)
#hollow_circle = MyCircle()


pose = pose_from_sequence("ACDEFGHIKLMNPQRSTVWY")
sf = ScoreFunction()
print "Score of the pose:", sf(pose)
len_score = LengthScoreMethod.scoreType
sf.set_weight(len_score, 1.0)
print "New score of the pose:", sf(pose)
