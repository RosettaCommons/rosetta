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

input1 = Input(shape=(3,), name="in1", dtype="float32" )
input2 = Input(shape=(2,), name="in2", dtype="float32" )

in1_dense = Dense( name="in1_dense", units=5, activation="relu" )( input1 )
in2_dense = Dense( name="in2_dense", units=10, activation="relu" )( input2 )

merge = tensorflow.keras.layers.concatenate( [in1_dense,in2_dense], name="merge", axis=-1 )

output = Dense( name="dense3", units=1, activation="sigmoid" )( merge )

model = Model(inputs=[input1, input2], outputs=output )


metrics_to_output=[ 'accuracy' ]
model.compile( loss='mean_squared_error', optimizer='adam', metrics=metrics_to_output )


for _ in range( 0, 1000 ):
    in1 = [np.random.random_sample((3,))]
    in2 = [np.random.random_sample((2,))]
    out = [np.random.random_sample((1,))]
    model.fit( x=[in1,in2], y=out )

test_in1 = [[ -1, 0, 1 ]]
test_in2 = [[ -2, 2 ]]
print( model.predict( [test_in1,test_in2] ) )
model.save( "test_model.h5" )
