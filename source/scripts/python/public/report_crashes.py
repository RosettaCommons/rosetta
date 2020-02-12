#!/usr/bin/env python
'''A script to submit Rosetta crash report files to the RosettaCommons.

This script will run through the current working directory and any subdirectories,
looking for Rosetta crash report files. If any are found, it will then submit
the files to the Rosetta crash report server at https://crash.rosettacommons.org/

Note that this script is just a convenience.
You can manually submit the crash report files at https://crash.rosettacommons.org/

Also, feel free to manually redact any sensitive or company proprietary information
that might be inadvertantly included in the crash report log before sending.
(For example, by using a text editor to replace it with "[REDACTED]".)
While the submitted contents of the crash report are not publically viewable by default,
RosettaCommons can make no gurantees about the long-term security of submitted information.
'''

from __future__ import print_function

import os, sys
import argparse
import socket

#Python 2/3 compatibility. Prefer raw_input if we've got it, but fall back to Python3's input if we don't.
try:
    get_input = raw_input
except NameError:
    get_input = input

if sys.version_info > (3, 0):
    def sending_encode(val):
        return val.encode(errors="replace")
    def recieve_decode(val):
        return val.decode(errors="replace")
else:
    def sending_encode(val):
        return val # Should already be a sendable string.
    def recieve_decode(val):
        return val # Should already be a sendable string.


##########################################################################

def collect_files(dir):
    crash_files = []
    for (root, dirs, files) in os.walk(dir):
        for filename in files:
            if filename.startswith("ROSETTA_CRASH"):
                crash_files.append( os.path.join( root, filename ) )

    return crash_files

HEADER = '''POST /report HTTP/1.0
Host: crash.rosettacommons.org
Content-length: {}

'''

def submit(filename, keep=False, verbose=False):
    with open(filename) as f:
        contents = f.read()

    body = HEADER.format(len(contents)) + contents

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        s.connect(("crash.rosettacommons.org", 80))
        s.sendall( sending_encode(body) )
        response = recieve_decode( s.recv(2048) ) # Should be more than enough for the type of responses we need.
    finally:
        s.close()

    split_response = response.split()

    if split_response[1] != "200": # This should be the server response code.
        print("Problem contacting crash report server. Please try again later.")
        if verbose:
            print("SERVER RESPONSE:")
            print(response)
        sys.exit() # Don't bother to continue processing if there's a server issue.

    if "FAIL" in split_response or "OK:" not in split_response: # Specifically with colon - we get an OK from the HTTP response
        print("Trouble sending file", filename, "to crash report server. Is it empty or malformed?")
        if verbose:
            print("SERVER RESPONSE:")
            print(response)
        return 0 # No reports this file, but we can continue to try other files.

    try:
        nreports = int( split_response[ split_response.index("OK:") + 1 ] )
    except IndexError:
        if verbose:
            print("Server response malformed - can't interpret the number of reports:")
            print(response)
        nreports = 0
    except ValueError:
        if verbose:
            print("Server response malformed - can't interpret the number of reports:")
            print(response)
        nreports = 0

    if keep: # Skip deletion
        return nreports

    if verbose:
        print("Deleting file", filename, "after a successful transmission.")

    os.remove(filename)

    return nreports

def main(args):
    crash_files = collect_files(args.dir)
    if args.verbose or not args.y:
        print( "Found {} crash file(s):".format(len(crash_files)) )
        for f in crash_files:
            print( "    " + f )
        print()
    if not args.y:
        while True:
            print("\nDo you want to submit these files? (y/n) ",end='')
            response = get_input().strip()
            if response[:1].lower() == "n":
                print("Cancelling crash report submission.")
                sys.exit()
            elif response[:1].lower() == "y":
                break
            else:
                print("Could not understand '{}'. Please enter 'y' or 'n'.".format(response) )

    print( "Attempting to submit {} crash file(s) to https://crash.rosettacommons.org/".format(len(crash_files)) )

    total_files, total_reports = 0, 0
    for filename in crash_files:
        nreports = submit(filename, keep=args.keep, verbose=args.verbose) # Will also delete if keep is false
        total_reports += nreports
        if nreports > 0:
            total_files += 1

    print("Submitted {} crash reports from {} files.".format(total_reports,total_files))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-y', action="store_true", help="Skip the confirmation dialog, and just send the reports.")
    parser.add_argument('-k','--keep',action="store_true", help="Keep the crash reports locally. By default, the script will delete the crash reports after sending to reduce double submission.")
    parser.add_argument('-d','--dir', default=".", help="The directory to work on. Defaults to the current working directory.")
    parser.add_argument('-v','--verbose',action="store_true",help="Print extra diagnostic information.")

    args = parser.parse_args()

    main(args)
