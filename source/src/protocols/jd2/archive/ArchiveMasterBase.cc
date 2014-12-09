// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MpiFileBuffer.hh
/// @brief  header file for MPISilentFileJobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @detail this outputter will send silentstructs via MPI to dedicated node that will collect all structures
/// @author Oliver Lange olange@u.washington.edu


#ifdef USEMPI
#include <mpi.h>
#endif


//unit headers
#include <protocols/mpi/ArchiveMasterBase.hh>

//project headers
#include <core/types.hh>

//utility headers
#include <utility/vector1.hh>
#include <utility/io/mpistream.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>

//C++ headers

#include <string>
#include <map>

namespace protocols {
namespace mpi {

ArchiveMasterBase::get_new_decoys( SilentStructOPs& new_structures ) {
  std::string filename = protocols::jd2::JobDistributor::get_instance()->job_outputter()->filename();
  MpiFileBuffer mpi_file_buffer( file_buffer_rank_ );
  mpi_file_buffer.block_file( filename );

  //.... do stuff

  mpi_file_buffer.release_file( filename );
}


}
}

#endif
