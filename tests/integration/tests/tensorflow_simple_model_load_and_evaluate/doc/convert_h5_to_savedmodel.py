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

from tensorflow.keras import backend as K
import tensorflow as tf

from tensorflow.python.util.tf_export import tf_export

# Tag for the `serving` graph.
SERVING = "serve"
tf_export(
    "saved_model.SERVING",
    v1=["saved_model.SERVING",
        "saved_model.tag_constants.SERVING"]).export_constant(
            __name__, "SERVING")


########
# INIT #
########

parser = argparse.ArgumentParser()
parser.add_argument( "--model", help="Most recent model file", required=True )
args = parser.parse_args()

if os.path.isfile( args.model ):
    model = load_model( args.model )
else:
    print( "Model " + args.model + " is not a file" )
    exit( 1 )


import time
saved_model_path = "./saved_models/{}".format(int(time.time()))

tf.keras.experimental.export_saved_model(model, saved_model_path)
saved_model_path

