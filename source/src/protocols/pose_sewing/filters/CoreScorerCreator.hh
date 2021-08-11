// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/filters/CoreScorerCreator.hh
/// @brief a filter that evaluates pairwise MotifScores
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_filters_CoreScorerCreator_hh
#define INCLUDED_protocols_pose_sewing_filters_CoreScorerCreator_hh

#include <protocols/filters/FilterCreator.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from FilterCreator.hh

namespace protocols {
namespace pose_sewing {
namespace filters {

class CoreScorerCreator : public protocols::filters::FilterCreator {
public:
	protocols::filters::FilterOP create_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} //protocols
} //pose_sewing
} //filters

#endif //INCLUDED_protocols_pose_sewing_filters_CoreScorer_fwd_hh
