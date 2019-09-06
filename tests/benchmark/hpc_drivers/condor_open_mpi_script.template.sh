#!/bin/sh
##**************************************************************
##
## Copyright (C) 1990-2010, Condor Team, Computer Sciences Department,
## University of Wisconsin-Madison, WI.
##
## Licensed under the Apache License, Version 2.0 (the "License"); you
## may not use this file except in compliance with the License.  You may
## obtain a copy of the License at
##
##    http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##
##**************************************************************

# YOU MUST CHANGE THIS TO THE PREFIX DIR OF OPENMPI
MPDIR=/usr/lib/openmpi

if `uname -m | grep "64" 1>/dev/null 2>&1`
then
    MPDIR=/usr/lib64/openmpi
fi

PATH=$MPDIR/bin:$MPDIR/1.4-gcc/bin:.:$PATH
export PATH

# This is a script to run openmpi jobs under the Condor parallel universe
# openmpi assumes that a full openmpi is pre-installed on all execute
# machines
# A sample submit file might look like...
#

#universe = parallel
#executable = openmpiscript
#getenv=true
#arguments = actual_mpi_job arg1 arg2 arg3
#
#should_transfer_files = yes
#when_to_transfer_output = on_exit_or_evict
#
#output = o.$(NODE)
#error  = e.$(NODE)
#log    = l
#
#notification = never
#machine_count = 8
#queue

_CONDOR_PROCNO=$_CONDOR_PROCNO
_CONDOR_NPROCS=$_CONDOR_NPROCS

CONDOR_SSH=`condor_config_val libexec`
CONDOR_SSH=$CONDOR_SSH/condor_ssh

SSHD_SH=`condor_config_val libexec`
SSHD_SH=$SSHD_SH/sshd.sh

. $SSHD_SH $_CONDOR_PROCNO $_CONDOR_NPROCS

# If not the head node, just sleep forever, to let the
# sshds run
if [ $_CONDOR_PROCNO -ne 0 ]
then
		wait
		sshd_cleanup
		exit 0
fi

EXECUTABLE=$1
shift

# the binary is copied but the executable flag is cleared.
# so the script have to take care of this
chmod +x $EXECUTABLE

CONDOR_CONTACT_FILE=$_CONDOR_SCRATCH_DIR/contact
export CONDOR_CONTACT_FILE

# The second field in the contact file is the machine name
# that condor_ssh knows how to use
sort -n -k 1 < $CONDOR_CONTACT_FILE | awk '{{print $2}}' > machines


let effective_number_of_process=$[_CONDOR_NPROCS * {process_coefficient}]

## run the actual mpijob
if `ompi_info --param all all | grep orte_rsh_agent 1>/dev/null 2>&1`
then
    mpirun -v --prefix $MPDIR --mca orte_rsh_agent $CONDOR_SSH -n $effective_number_of_process -hostfile machines $EXECUTABLE $@
else
    ########## For mpi versions 1.1 & 1.2 use the line below
    mpirun -v --mca plm_rsh_agent $CONDOR_SSH -n $effective_number_of_process -hostfile machines $EXECUTABLE $@
fi

sshd_cleanup
rm -f machines

exit $?
