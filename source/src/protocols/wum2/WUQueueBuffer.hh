// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/wum2/WUQueueBuffer.hh
/// @brief  memory aware structure that links irecv and isend
/// 				with their mpi::request status, so buffers are maintained until messages are sent/recvd
/// @author Ken Jung

#ifndef INCLUDED_protocols_wum2_WUQueueBuffer_hh
#define INCLUDED_protocols_wum2_WUQueueBuffer_hh


#ifdef USEBOOSTMPI
// this is useless without mpi
#include <boost/mpi.hpp>

#include <protocols/wum2/WUQueueBuffer.fwd.hh>
#include <protocols/wum2/WorkUnit.fwd.hh>
#include <boost/tuple/tuple.hpp>
#include <boost/cstdint.hpp>
#include <list>
#include <vector>

namespace protocols {
namespace wum2 {

class WUQueueBuffer {


  public:
		typedef boost::tuple< boost::uint64_t, boost::mpi::request, boost::shared_ptr< std::vector< WorkUnitSP > >  >  mem_req_data_tuple;
		typedef std::list< mem_req_data_tuple>::reverse_iterator riterator;
		typedef std::list< mem_req_data_tuple>::iterator iterator;

    WUQueueBuffer():current_mem_(0){}
    ~WUQueueBuffer(){}

    boost::uint64_t current_mem() { return current_mem_; }

    // allocates an entry in buffer for delayed linking to a mpi::request
    // size isn't used in actually allocating the buffer, but used to track memory usage
    // if buffer isn't linked to real mpi request, can have psuedo memory leak
    riterator allocate_buffer( boost::uint64_t size );

    // pops off all requests that are completed
    std::vector< WorkUnitSP > cleanup_reqs();

  private:
    boost::uint64_t current_mem_;
    std::list< mem_req_data_tuple> buffer_;
};

}
}

#endif
#endif

