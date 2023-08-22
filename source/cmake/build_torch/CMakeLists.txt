CMAKE_MINIMUM_REQUIRED(VERSION 3.18)
PROJECT(RosettaTorch)
SET(CMAKE_CXX_STANDARD 14)

##############################################################################
# This file only needs to specify which mode (i.e., release, debug),
# extras (i.e., graphics), and compiler (i.e., gcc, clang) should be used. The
# appropriate settings for each mode, extras, and compiler should be
# implemented in ../build/build.settings.cmake.
##############################################################################

############ Build Type, Mode, Link Rules, Compiler, Extras
SET(MODE release)
SET(LINK_RULES shared)
SET(CMAKE_EXPORT_COMPILE_COMMANDS 1)

############ Build Rules
INCLUDE(../build/compiler.cmake)
INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)
INCLUDE(../build/build.settings.cmake)
INCLUDE(../build/build.post.cmake)

############ PyTorch Rules
# The path is hard coded to follow the SCons compile scheme
# https://github.com/RosettaCommons/main/pull/6281#issuecomment-1302466301

# Download PyTorch source files:
# cd source/external/
# curl https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.13.0%2Bcpu.zip -o libtorch.zip
# unzip libtorch.zip && rm libtorch.zip
LIST(APPEND CMAKE_PREFIX_PATH "../../external/libtorch/share/cmake/")

# The cuDNN library is not available with standard libtorch
#SET(CAFFE2_USE_CUDNN 1)
FIND_PACKAGE(Torch REQUIRED)
ADD_DEFINITIONS(-DUSE_PYTORCH)
SET(COMPILE_FLAGS "${COMPILE_FLAGS} ${TORCH_CXX_FLAGS}")
SET(LINK_ALL_LIBS ${LINK_ALL_LIBS} ${TORCH_LIBRARIES})
INCLUDE_DIRECTORIES("../../external/libtorch/include")

########### Uncomment the following if you only want to build RosettaScripts:
#INCLUDE( ../build/apps/rosetta_scripts.cmake )

############ Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )

############ Unit Tests
#INCLUDE(../build/build.test.cmake)

############ Installation rules
#INCLUDE( ../build/install.cmake )