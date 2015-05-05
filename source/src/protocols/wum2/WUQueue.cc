// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/wum2/WUQueue.hh
/// @brief  deque of WU with memory tracking
/// @author Ken Jung

#include <protocols/wum2/WUQueue.hh>
#include <protocols/wum2/WorkUnit.hh>
#include <sstream>
#include <iostream>

namespace protocols{
namespace wum2{

boost::uint64_t WUQueue::serialized_size( WorkUnitSP /*wu*/ ) {
#ifdef SERIALIZE
/*
  std::stringstream s;
  core::io::serialization::toBinary(s, wu);
  return s.str().length();
*/
  return 0;
#else
	std::cerr << "Memory usage tracked only if compiled against boost::serialize" << std::endl;
  return 0;
#endif
}

WorkUnitSP WUQueue::pop_front() {
  WorkUnitSP tmp;
  if( empty() )
    return tmp;
  current_mem_ -= deque_.front().first;
  tmp = deque_.front().second;
  deque_.pop_front();
  return tmp;
}

// so much unnecessary copying
std::vector<WorkUnitSP> WUQueue::pop_all() {
	std::vector<WorkUnitSP> tmp;
	for( std::deque<wu_mem_pair>::iterator itr=deque_.begin(); itr != deque_.end(); itr++) {
		tmp.push_back( itr->second );
	}
	deque_.clear();
	return tmp;
}

void WUQueue::push_back( std::vector<WorkUnitSP> wulist ) {
	for( std::vector<WorkUnitSP>::iterator itr=wulist.begin(); itr != wulist.end(); itr++) {
		push_back( *itr );
	}
}

void WUQueue::push_front( WorkUnitSP wu ) {
  boost::uint64_t mem_size = serialized_size( wu );
  deque_.push_front( wu_mem_pair( mem_size, wu ) );
  current_mem_ += mem_size;
}

void WUQueue::push_back( WorkUnitSP wu ) {
  boost::uint64_t mem_size = serialized_size( wu );
  deque_.push_back( wu_mem_pair( mem_size, wu ) );
  current_mem_ += mem_size;
}

}// wum2
}// protocols
