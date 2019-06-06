#This file creates a very small network to be used for the integration test
import os

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras import metrics

from tensorflow.keras.models import load_model
import tensorflow.keras.backend as K
import tensorflow.keras.callbacks
import tensorflow.keras
import numpy

import sys
import h5py

import argparse
import random
import time
import subprocess

import tensorflow as tf

########
# INIT #
########

parser = argparse.ArgumentParser()
parser.add_argument( "--model", help="filename for output file", default="very_small_blank_model", required=False )
args = parser.parse_args()

# Create Layers

num_input_dimensions = 10
num_neurons_in_middle_layer = 5
num_output_dimensions = 2
model = Sequential()

model.add( Dense( num_neurons_in_middle_layer, input_dim=num_input_dimensions, activation='relu') )
model.add( Dense( num_output_dimensions, activation='linear') )

# Compile Model

metrics_to_output=[ 'accuracy' ]
model.compile( loss='mean_squared_error', optimizer='adam', metrics=metrics_to_output )

# Save Model

model.save( args.model + ".h5" )

# Run model so that we know what output to expect

dummy_test_input = numpy.zeros( shape=( 1, 10 ) )
for x in range( 0, 9 ):
    dummy_test_input[ 0 ][ x ] = 3
print( "Network Input:" )
print( dummy_test_input )
    
predictions = model.predict( x=dummy_test_input[:] )
print( "Network Output:" )
print( predictions )
