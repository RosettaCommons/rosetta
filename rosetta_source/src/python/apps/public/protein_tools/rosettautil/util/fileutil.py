import gzip
import sys

def universal_open(path,mode):
    """given a path and a mode, return a filehandle.
    allows you to open files regardless of whether they are gzipped"""
    extension = path.split(".")[-1]
    if(extension == "gz"):
        try:
            return gzip.open(path,mode)
        except IOError:
            sys.exit("could not open "+path+" in mode "+mode)
    else:
        try:
            return open(path,mode)
        except IOError:
            sys.exit("could not open "+path+" in mode "+mode)
