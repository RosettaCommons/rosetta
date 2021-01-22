################################################################################
##
## @file integration/tests/basic_gcn_tensorflow_test/setup_code/create_gcn.py
##
## @brief Python3 code to create a simple graph convolutional neural network.
##
## @details This requires the Spektral package ("pip install spektral") and
## Tensorflow 2.  Note that this GCN doesn't _do_ anything interesting.
##
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
##
## @note This script also requires pydot ("pip install pydot"), pydotplus
## ("pip install pydotplus"), and graphviz ("brew install graphviz" on MacOS).
## For some reason, pip installation of graphviz doesn't work on MacOS.
##
## Also note that the model produced should be copied to the database (to
## database/protocol_data/tensorflow_graphs.  The model's input and output
## names can be interrogated with the saved_model_cli program (part of the
## Tensorflow 2 installation) as follows:
##
## saved_model_cli show --dir <path>/gcn_test_model/ --all
##
## In the above, <path> should be replaced with the path to the model.  The
## start of the output looks like this:
##
## MetaGraphDef with tag-set: 'serve' contains the following SignatureDefs:
## 
## signature_def['__saved_model_init_op']:
##   The given SavedModel SignatureDef contains the following input(s):
##   The given SavedModel SignatureDef contains the following output(s):
##     outputs['__saved_model_init_op'] tensor_info:
##         dtype: DT_INVALID
##         shape: unknown_rank
##         name: NoOp
##   Method name is: 
## 
## signature_def['serving_default']:
##   The given SavedModel SignatureDef contains the following input(s):
##     inputs['A_in'] tensor_info:
##         dtype: DT_FLOAT
##         shape: (-1, 7, 7)
##         name: serving_default_A_in:0
##     inputs['E_in'] tensor_info:
##         dtype: DT_FLOAT
##         shape: (-1, 7, 7, 2)
##         name: serving_default_E_in:0
##     inputs['X_in'] tensor_info:
##         dtype: DT_FLOAT
##         shape: (-1, 7, 2)
##         name: serving_default_X_in:0
##   The given SavedModel SignatureDef contains the following output(s):
##     outputs['OutputLayer'] tensor_info:
##         dtype: DT_FLOAT
##         shape: (-1, 7, 1)
##         name: StatefulPartitionedCall:0
##   Method name is: tensorflow/serving/predict
## 
## Given the above output, the input names are "serving_default_A_in",
## "serving_default_E_in", and "serving_default_X_in", and the output name is
## "StatefulPartitionedCall".  These are needed when writing the C++ code to
## run the model.
##
################################################################################

# Include Tensorflow and Spektral:
from tensorflow import bool, float32, uint8, constant
from tensorflow.keras import Input, Model
from tensorflow.keras.layers import Dense
from tensorflow.keras.utils import plot_model
from spektral.layers import EdgeConditionedConv
print( "Initialization complete." )

# Build a simiple graph with seven nodes and eight edges.
# Connectivity should be 1-2, 1-3, 2-4, 3-4, 2-5, 4-5, 4-6, 6-7:
A = constant(
    [[
        [0, 1, 1, 0, 0, 0, 0],
        [1, 0, 0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0, 0, 0],
        [0, 1, 1, 0, 1, 1, 0],
        [0, 1, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 1],
        [0, 0, 0, 0, 0, 1, 0]
    ]],
    dtype=float32
)
print(A)

# Add some values to the nodes.
# Nodes will have 2-vectors assigned to each.
X = constant(
    [[
        [0.2,  0.5],
        [0.8,  -0.5],
        [0.4,  0.2],
        [0.7,  -0.8],
        [-0.1,  0.5],
        [0.0,  0.7],
        [0.1,  0.8]
    ]],
    dtype=float32
)
print(X)

# Add some values to the edges.
# Edges will have a 2-vector assigned to each.
E = constant(
    [[
        [ [0,0], [7,3], [6,1], [0,0], [0,0], [0,0], [0,0] ],
        [ [7,3], [0,0], [0,0], [8,9], [5,5], [0,0], [0,0] ],
        [ [6,1], [0,0], [0,0], [6,6], [0,0], [0,0], [0,0] ],
        [ [0,0], [8,9], [6,6], [0,0], [8,3], [7,2], [0,0] ],
        [ [0,0], [5,5], [0,0], [8,3], [0,0], [0,0], [0,0] ],
        [ [0,0], [0,0], [0,0], [7,2], [0,0], [0,0], [1,9] ],
        [ [0,0], [0,0], [0,0], [0,0], [0,0], [1,9], [0,0] ]
    ]],
    dtype=float32
)
print(E)

# Build an edge-conditioned graph convolutional neural network.  We'll randomly initialize the weights.
A_in = Input( name="A_in", dtype=float32, shape=(7,7) )
X_in = Input( name="X_in", dtype=float32, shape=(7,2))
E_in = Input( name="E_in", dtype=float32, shape=(7,7,2))
Conv1 = EdgeConditionedConv( 3, name="Conv1Layer", kernel_network=(2,2), activation='elu', kernel_initializer='glorot_uniform', bias_initializer='zeros', dtype=float32 )([X_in,A_in,E_in])
output = Dense( 1, name="OutputLayer", activation='elu', kernel_initializer='glorot_uniform', bias_initializer='zeros', dtype=float32 )(Conv1)

# Package these layers into a model:
model = Model( name="MyModel", inputs=[X_in,A_in,E_in], outputs=output )
model.compile()
model._layers = [layer for layer in model._layers if not isinstance(layer, dict)] #For some reason, an empty dict somewhere impairs visualization.
model.summary()
print( model.get_layer(name="Conv1Layer").variables )

# Visualize the model:
plot_model(model, to_file="gcn_test_model_plot.png", show_shapes=True, show_layer_names=True, expand_nested=True)

# Run the model twice (to confirm stable output):
result1 = model( [X,A,E] )
result2 = model( [X,A,E] )

# Print the output:
print( "***RESULT 1***" )
print(result1)
print( "***RESULT 2***" )
print(result2)

# Save the model:
model.save( "gcn_test_model", save_format='tf', overwrite=True )
print( "Wrote model to gcn_test_model." )
