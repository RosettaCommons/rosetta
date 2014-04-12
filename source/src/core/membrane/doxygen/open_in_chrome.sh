#!/usr/bin/env sh

absolute_path=$(pwd)/$(dirname $0)/html/index.html
clean_path=$(realpath -s $absolute_path)
chrome file://$clean_path
