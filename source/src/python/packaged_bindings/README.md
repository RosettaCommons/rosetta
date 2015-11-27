#Preparation

1) Install `boost`, `gcc-xml`, `python`, `pygccxml`, and `numpy`. For example, on osx via homebrew:

````
brew install python
brew install boost-python

# Legacy gcc needed to support gccxml
brew install apple-gcc42
brew install --HEAD homebrew/head-only/gccxml

#Ensure targeting homebrew python
which pip

pip install numpy 
pip install pygccxml
````


#Build Instructions


1) Execute `source/src/python/packaged_bindings/BuildPackagedBindings.py` optionally specifying your source `python` and `boost_python` libraries as well as desired result path.

````
source/src/python/packaged_bindings/BuildPackagedBindings.py 
````

2) Execute binding tests.

```
pushd pyrosetta
python setup.py nosetests
popd
```

3) Build installation wheel.

````
pushd pyrosetta
python setup.py wheel
popd
````

4) Install wheel

````
pip install -f pyrosetta/dist --pre pyrosetta
````

To install into a "standalone" directory, suitable for inclusion via PYTHONPATH:

````
pip install -f pyrosetta/dist --target=<target_dir> --pre pyrosetta
````

Wheel install requires unpacking the target wheel into a temporary "build" directory. To customize the build dir:

````
pip install -f pyrosetta/dist --build=<target_build_dir> --pre pyrosetta
````
