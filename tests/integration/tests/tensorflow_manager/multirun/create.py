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

input = Input(shape=(3,), name="in1", dtype="float32" )
dense = Dense( name="dense", units=5, activation="relu" )( input )
output = Dense( name="out", units=1, activation="sigmoid" )( dense )

model = Model(inputs=input, outputs=output )


metrics_to_output=[ 'accuracy' ]
model.compile( loss='mean_squared_error', optimizer='adam', metrics=metrics_to_output )


for _ in range( 0, 1000 ):
    in1 = [np.random.random_sample((3,))]
    out = [np.random.random_sample((1,))]
    model.fit( x=[in1], y=out )

print( model.predict( [[[ -1, 0, 1 ]]] ) )
print( model.predict( [[[  1, 2, 3 ]]] ) )
print( model.predict( [[[  0, 0, 0 ]]] ) )

model.save( "test_model.h5" )
