#!/bin/bash

rm -rf ../bin/* ../tools/bin/* CMakeFiles CMakeCache.txt cmake_install.cmake Makefile src install_manifest.txt apbs.cbp maloc-prefix tools fontconfig ../tests/*.pyc ../tests/io.mc ../tests/*~ ../tests/test.log doc ../doc/programmer/html ../doc/programmer/programmer.html ../doc/programmer/mainpage.h ../doc/programmer/Doxyfile ../doc/programmer/latex ../doc/programmer/programmer.pdf

find .. -name *~ | xargs rm

