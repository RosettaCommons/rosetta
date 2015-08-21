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

#ifndef INCLUDED_protocols_wum2_WUQueue_hh
#define INCLUDED_protocols_wum2_WUQueue_hh

#include <protocols/wum2/WUQueue.fwd.hh>
#include <protocols/wum2/WorkUnit.fwd.hh>
#include <boost/cstdint.hpp>
#include <deque>
#include <vector>

namespace protocols {
namespace wum2 {

class WUQueue {
	typedef std::pair< boost::uint64_t, WorkUnitSP > wu_mem_pair;
	typedef std::deque< wu_mem_pair >::iterator iterator;

public:
	WUQueue(): current_mem_(0){}
	~WUQueue(){}

	boost::uint64_t current_mem() { return current_mem_; }

	boost::uint64_t size_front() { return empty() ? 0 : deque_.front().first; }

	void push_front( WorkUnitSP wu );
	void push_back( WorkUnitSP wu );
	void push_back( std::vector<WorkUnitSP> wulist );
	WorkUnitSP pop_front();
	std::vector<WorkUnitSP> pop_all();
	bool empty() { return deque_.empty(); }
	void clear() { deque_.clear(); }
	int size() { return deque_.size(); }

private:
	boost::uint64_t serialized_size( WorkUnitSP wu );

	boost::uint64_t current_mem_;
	std::deque< wu_mem_pair > deque_;
};

}
}

#endif

