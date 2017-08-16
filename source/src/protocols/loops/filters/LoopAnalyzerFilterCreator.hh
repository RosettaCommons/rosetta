// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

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
	protocols::filters::FilterOP create_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //protocols
} //loops
} //filters

#endif //INCLUDED_protocols_loops_filters_LoopAnalyzerFilter_fwd_hh



