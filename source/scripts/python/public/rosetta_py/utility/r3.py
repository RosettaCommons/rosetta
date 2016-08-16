#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

'''Functions for dealing with vectors and points in 3-space (R3).

Adapted from IWD's driftwood.r3.Triple class in Java.
The object is to be correct, not necessarily fast.

Functions work on anything with .x, .y, and .z properties defined.
Later I may include support for numeric indexing instead [0], [1], [2]
and/or support for NumPy array objects.

Author: Ian W. Davis
'''
import math

'''A simple class for XYZ triples as the result of computations.'''
class Triple:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "(%8.3f, %8.3f, %8.3f)" % (self.x, self.y, self.z)

def is_nan(f):
    '''Stupid implementation so we're not dependent on outside libs.
    (NumPy has one that's probably better, isnan().)
    Python 2.3 has a bug where NaN1 == NaN2 is true, so we do it this way.'''
    return str(f) == 'nan'

def add(v1, v2, vOut=None):
    if vOut is None: vOut = Triple()
    # Do this in two steps in case v1 == vOut or v2 == vOut
    x = v1.x + v2.x;
    y = v1.y + v2.y;
    z = v1.z + v2.z;
    vOut.x = x
    vOut.y = y
    vOut.z = z
    return vOut

def sub(v1, v2, vOut=None):
    if vOut is None: vOut = Triple()
    # Do this in two steps in case v1 == vOut or v2 == vOut
    x = v1.x - v2.x;
    y = v1.y - v2.y;
    z = v1.z - v2.z;
    vOut.x = x
    vOut.y = y
    vOut.z = z
    return vOut

def mult(v, k, vOut=None):
    '''Multiplies this vector by k.'''
    if vOut is None: vOut = Triple()
    # Do this in two steps in case v == vOut
    x = v.x * k
    y = v.y * k
    z = v.z * k
    vOut.x = x
    vOut.y = y
    vOut.z = z
    return vOut

def div(v, k, vOut=None):
    '''Divides this vector by k.'''
    if vOut is None: vOut = Triple()
    # Do this in two steps in case v == vOut
    x = v.x / k
    y = v.y / k
    z = v.z / k
    vOut.x = x
    vOut.y = y
    vOut.z = z
    return vOut

def midpoint(v1, v2, vOut=None):
    if vOut is None: vOut = Triple()
    # Do this in two steps in case v1 == vOut or v2 == vOut
    x = v1.x + v2.x;
    y = v1.y + v2.y;
    z = v1.z + v2.z;
    vOut.x = x / 2.0
    vOut.y = y / 2.0
    vOut.z = z / 2.0
    return vOut

def from_to(v1, v2, vOut=None):
    '''The vector originating at v1 and pointing to v2 (v2 - v1)'''
    return sub(v2, v1, vOut)

def mag2(v):
    '''Returns the squared maginitude of a vector from the origin to this point.
    This is equivalent to the dot product of the vector with itself.'''
    return (v.x**2) + (v.y**2) + (v.z**2)

def mag(v):
    '''Returns the maginitude of a vector from the origin to this point.'''
    return math.sqrt(mag2(v))

def unit(v, vOut=None):
    '''Makes this vector one unit in length (magnitude) with the same directionality.
    If this vector is (0,0,0), no change is made.'''
    if vOut is None: vOut = Triple()
    m = mag(v)
    if m != 0: div(v, m, vOut)
    return vOut

def dot(v1, v2):
    '''Returns the vector dot product of these two vectors.
    The dot product of A and B, A.B, is equal to |A||B|cos(theta),
    where theta is the angle between vectors from the origin to A and B.'''
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z

def cross(v1, v2, vOut=None):
    '''Returns a new vector equal to the cross product of v1 and v2.
    The cross product of A and B, AxB, is orthogonal to the plane defined by vectors
    from the origin to A and B. Its direction (sign) is given by the right-hand rule.'''
    if vOut is None: vOut = Triple()
    # Do this in two steps in case v1 == vOut or v2 == vOut
    x = v1.y*v2.z - v1.z*v2.y
    y = v1.z*v2.x - v1.x*v2.z
    z = v1.x*v2.y - v1.y*v2.x
    vOut.x = x
    vOut.y = y
    vOut.z = z
    return vOut

def distance(a, b):
    '''Distance between two points.'''
    x = a.x - b.x
    y = a.y - b.y
    z = a.z - b.z
    return math.sqrt( x**2 + y**2 + z**2 )

def angle(a, b, c=None):
    '''Returns an angle between two vectors or three points, in degrees from 0 to 180.'''
    if c is not None:
        u = Triple(); v = Triple()
        sub(a, b, u)
        sub(c, b, v)
        a = u
        b = v
    dt = dot(a, b)
    # for Java at least, acos returns NaN sometimes when we're
    # too close to an angle of 0 or 180 (if |dot/mag| > 1)
    # Python behavior varies by platform:  may be NaN or throw an exception
    try:
        ang = math.acos( dt / (mag(a)*mag(b)) )
    except:
        if dt >= 0: return 0.0
        else: return 180.0
    if is_nan(ang):
        if dot >= 0: return 0.0
        else: return 180.0
    else: return math.degrees(ang)


def dihedral(a, b, c, d):
    '''Returns the dihedral ABCD, in degrees, from -180 to 180.'''
    # Schematic for naming points (uppercase) and vectors (lowercase).
    # Shown is a dihedral of 0 degress, with imaginary vectors u and v.
    # ABCD are all in a plane; u, v, and f are in a plane orthogonal to it.
    #
    #     ^        ^
    #   v :        : u
    #     :        :
    #     C<-------B
    #  g /    f    /\
    #   /            \ e
    # \/              \
    # D                A
    e = Triple(); f = Triple(); g = Triple()
    sub(b, a, e)
    sub(c, b, f)
    sub(d, c, g)
    u = cross(e, f, e) # overwrites e
    v = cross(f, g, f) # overwrites f
    dihe = angle(u, v)
    # Solve handedness problem:
    if angle(u, g) > 90.0: dihe = -dihe
    return dihe

def centroid(vs, vOut=None):
    '''Returns an unweighted center-of-mass for the list of points.'''
    if vOut is None: vOut = Triple()
    # Do this in two steps in case vOut is part of vs
    x = 0.0; y = 0.0; z = 0.0
    for v in vs:
        x += v.x
        y += v.y
        z += v.z
    l = len(vs)
    vOut.x = x/l
    vOut.y = y/l
    vOut.z = z/l
    return vOut

