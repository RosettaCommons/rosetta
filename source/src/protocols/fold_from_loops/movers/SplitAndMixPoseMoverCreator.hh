// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   MoveMapFromSelectorMover.hh
/// @brief  Generates a movemap from a ResidueSelector
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

#ifndef INCLUDED_protocols_fold_from_loops_movers_SplitAndMixPoseMoverCreator_hh
#define INCLUDED_protocols_fold_from_loops_movers_SplitAndMixPoseMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace fold_from_loops {
namespace movers {

/// @brief RosettaScripts factory for NubInitio.
class SplitAndMixPoseMoverCreator : public protocols::moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}
}
}

#endif
