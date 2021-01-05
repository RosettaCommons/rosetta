import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = ""


import tensorflow

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Reshape
from tensorflow.keras.layers import Input
from tensorflow.keras.layers import Flatten
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import MaxPooling2D
from tensorflow.keras.layers import LocallyConnected2D
from tensorflow.keras.layers import LeakyReLU
from tensorflow.keras import metrics
from tensorflow.keras import optimizers
from tensorflow.keras.layers import UpSampling2D
from tensorflow.keras.layers import UpSampling3D

from tensorflow.keras.models import load_model
from tensorflow.keras.models import Model
import tensorflow.keras.backend as K

import numpy as np

model = load_model( "test_model.h5" )

test_in1 = [[ -1, 0, 1 ]]
test_in2 = [[ -2, 2 ]]
print( test_in1, test_in2, model.predict( [test_in1,test_in2] ) )


test_in1 = [[ 1, 2, 3 ]]
test_in2 = [[ -4, 2 ]]
print( test_in1, test_in2, model.predict( [test_in1,test_in2] ) )


test_in1 = [[ -1, 5, 1 ]]
test_in2 = [[ -5, 2 ]]
print( test_in1, test_in2, model.predict( [test_in1,test_in2] ) )


test_in1 = [[ -1, 0, -1 ]]
test_in2 = [[ -2, -2 ]]
print( test_in1, test_in2, model.predict( [test_in1,test_in2] ) )

test_in1 = [[ -11, 0, 11 ]]
test_in2 = [[ -22, 22 ]]
print( test_in1, test_in2, model.predict( [test_in1,test_in2] ) )

''' Output:
[[-1, 0, 1]] [[-2, 2]] [[0.500397]]
[[1, 2, 3]] [[-4, 2]] [[0.40970826]]
[[-1, 5, 1]] [[-5, 2]] [[0.40430167]]
[[-1, 0, -1]] [[-2, -2]] [[0.3356494]]
[[-11, 0, 11]] [[-22, 22]] [[0.71121025]]
'''
