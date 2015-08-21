// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/wum2/WUQueueBuffer.cc
/// @brief  Foward decls for WUQueueBuffer, memory aware structure that links irecv and isend
///     with their mpi::request status, so buffers are maintained until messages are sent/recvd
/// @author Ken Jung

#ifdef USEBOOSTMPI

#include <protocols/wum2/WUQueueBuffer.hh>

namespace protocols {
namespace wum2 {

// mpi-aware, links buffer for an irecv/isend with the mpi::request status


WUQueueBuffer::riterator 
WUQueueBuffer::allocate_buffer(boost::uint64_t size) {
  boost::mpi::request req;
	boost::shared_ptr< std::vector < WorkUnitSP > > vec( new std::vector< WorkUnitSP > ) ;
  buffer_.push_back( boost::make_tuple(size,req, vec) );  
  current_mem_ += size;
  return buffer_.rbegin();
}

std::vector< WorkUnitSP > 
WUQueueBuffer::cleanup_reqs() {
  std::vector< WorkUnitSP > vec;
	iterator itr = buffer_.begin();
	while( itr != buffer_.end() ) {
     if( itr->get<1>().test().is_initialized() ){
       current_mem_ -= itr->get<0>();
       vec.insert( vec.end(), itr->get<2>()->begin(), itr->get<2>()->end() );
       itr = buffer_.erase(itr);
     } else {
			 ++itr;
		 }
  }
  return vec;
}

} // wum2
} // protocols
#endif
