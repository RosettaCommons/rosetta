#Preparation

1) Install `boost`, `gcc-xml`, `python`, `pygccxml`, and `numpy`.

#Build Instructions

_Example commands assume digs environment._

1) Execute `source/src/python/packaged_bindings/BuildPackagedBindings.py` specifying your source `python` and `boost_python` libraries as your desired result path.

````
source/src/python/packaged_bindings/BuildPackagedBindings.py --python_lib=python2.7 --boost_lib=boost_python --boost_path=/work/buildbot/opt --python_path=/usr/local --compiler=gcc --jobs=20 --bindings_path=pyrosetta/rosetta
````

2) Symlink database into pyrosetta destination folder.

````
pushd pyrosetta
ln -s ../database
popd
````

3) Build installation egg.

````
cp source/src/python/packaged_bindings/packaging/* pyrosetta

pushd pyrosetta
python setup.py bdist_egg
popd
````

4) Install egg.

````
pushd pyrosetta
easy_install -U -H None -f dist/*.egg pyrosetta
popd
````
