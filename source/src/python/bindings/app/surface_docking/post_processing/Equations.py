#!/usr/bin/env python

from math import sqrt, pi, acos, cos, sin

def vector(pointA, pointB):
    vector = [0,0,0]
    vector[0] = pointB[0] - pointA[0]
    vector[1] = pointB[1] - pointA[1]
    vector[2] = pointB[2] - pointA[2]

    return vector

def length(vector):
    length = sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

    return length

def unit_vector(vector, length):
    unit_vector = [vector[0]/length, vector[1]/length, vector[2]/length]

    return unit_vector

def dot_product(vector1, vector2):
    dot_product = vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2]

    return dot_product

def cross_product(vector1, vector2):
    cross_product = [0,0,0]
    cross_product[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1]
    cross_product[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2]
    cross_product[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0]

    return cross_product

def angle(vector1, vector2):
    unit_vector1 = unit_vector(vector1, length(vector1))
    unit_vector2 = unit_vector(vector2, length(vector2))
    angle = acos(dot_product(unit_vector1, unit_vector2))*180/pi

    return angle

def scalar_triple_product(vector1, vector2, vector3):
    scalar_triple_product = dot_product( vector1, cross_product(vector2, vector3) )

    return scalar_triple_product
