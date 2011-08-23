#!/bin/sh

#THIS IS NOT GUARUNTEED TO NOT DESTROY YOUR CHANGES - it just deletes everything svn marks with ?

svn st | grep "?" | awk '{print $2 }' | xargs rm -r

exit 0
