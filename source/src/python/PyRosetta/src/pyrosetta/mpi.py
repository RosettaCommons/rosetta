# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# MPI support for PyRosetta
# @author Sergey Lyskov

from __future__ import print_function

from mpi4py import MPI

import pyrosetta


# MPI version of init function, use it instead of init(...)
def mpi_init(*args, **kargs):

    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()

    kargs['extra_options'] = kargs.get('extra_options', '') + ' -seed_offset {}'.format(rank*10000)

    pyrosetta.init(*args, **kargs)


def MPIJobDistributor(njobs, fun):
    ''' Execute fun across MPI nodes to complete njobs.
        arguments:
          - int: njobs
          - fun(int job) function to execute N'th job

        Example: suppose we have 64 decoys to generate using 4 MPI nodes.
        In this case we will need to call MPIJobDistributor(64, foo) and the following foo(...) calls will be made on each node:

        Node 0: for i in [0, 4,  8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60]: foo(i)
        Node 1: for i in [1, 5,  9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61]: foo(i)
        Node 2: for i in [2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62]: foo(i)
        Node 3: for i in [3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63]: foo(i)
    '''

    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()

    myjobs = []

    if rank == 0:
        jobs = list( range(njobs) )
        jobs.extend( [None]*(size - njobs % size) )
        n = len(jobs) // size
        for i in range(size):
            queue = []  # list of jobs for individual cpu
            for j in range(n):
                queue.append(jobs[j*size+i])

            if( i == 0 ):
                myjobs = queue
            else:
                # now sending the queue to the process
                logger.info('Sending {} to node {}'.format(queue, i) )
                comm.send(queue, dest=i)
    else:
        # getting decoy lists
        myjobs = comm.recv(source=0)

    logger.info('Node {}, got queue:{}'.format(rank, myjobs) )

    for j in myjobs:
        if j is not None: fun(j)
