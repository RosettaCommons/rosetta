// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/loops/filters/LoopAnalyzerFilterCreator.hh
/// @brief intense analysis of loop quality (creator)
/// @author Steven Lewis (smlewi@gmail.com)

#ifndef INCLUDED_protocols_loops_filters_LoopAnalyzerFilterCreator_hh
#define INCLUDED_protocols_loops_filters_LoopAnalyzerFilterCreator_hh

#include <protocols/filters/FilterCreator.hh>

namespace protocols {
namespace loops {
namespace filters {

class LoopAnalyzerFilterCreator : public protocols::filters::FilterCreator {
public:
	virtual protocols::filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
	static std::string filter_name();
};

} //protocols
} //loops
} //filters

#endif //INCLUDED_protocols_loops_filters_LoopAnalyzerFilter_fwd_hh



