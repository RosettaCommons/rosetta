P=`pwd`
#export APBS_SRC=${P}/apbs-1.4-dev
#export APBS_PREFIX=${P}/apbs-1.4-dev
export APBS_SRC=${P}
export APBS_PREFIX=${P}

export MCSH_HOME=/dev/null
	
# get the source
#git clone git://git.code.sf.net/p/apbs/code apbs-1.4-dev

# create a build directory and cd to it, then:
cd ${APBS_SRC}/build
sh cleanup.sh

cmake -DCMAKE_INSTALL_PREFIX:PATH=${APBS_PREFIX} -DENABLE_QUIET=ON \
-DBUILD_TOOLS=OFF -DENABLE_PYTHON=OFF -DENABLE_OPENMP=OFF \
-DENABLE_iAPBS=ON -DBUILD_WRAPPER=ON -DENABLE_MPI=OFF \
-DENABLE_DEBUG=OFF \
-DBUILD_SHARED_LIBS=ON ..

make

cp ${APBS_SRC}/src/routines.h ${APBS_SRC}/include/apbs
cp ${APBS_SRC}/contrib/iapbs/src/apbs_driver.h ${APBS_SRC}/include/iapbs
