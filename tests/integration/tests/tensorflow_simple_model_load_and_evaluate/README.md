# Integration test "tensorflow\_simple\_model\_load\_and\_evaluate"

##Author

Jack Maguire

## Description

This integration test confirms that Rosetta can load in a simple TensorFlow model and evaluate it correctly.

The model was created by using tensorflow 2.0 alpha and keras, using the scripts in the docs folder:
```
cd docs
python3 create_very_small_model.py
mkdir saved_models
python3 ../utilities/convert_h5_to_pb.py --model very_small_blank_model.h5
```

## Expected output

The model that we are using should output
`-0.40661177 -0.5512837`
when we input
`3. 3. 3. 3. 3. 3. 3. 3. 3. 0.`

Note that the weights are randomly assigned so you will want to update the expected values in the c++ app with the output of the new run.
