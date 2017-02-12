// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/LocalBackboneMoverCreator.hh
/// @brief LocalBackboneMover moves a stretch of backbone locally.
/// @author xingjiepan (xingjiepan@gmail.com)

#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_LocalBackboneMoverCreator_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_LocalBackboneMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from MoverCreator

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {

class LocalBackboneMoverCreator : public protocols::moves::MoverCreator {

public:

	protocols::moves::MoverOP
	create_mover() const override;

	std::string
	keyname() const override;

	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;


};

} //protocols
} //backbone_moves
} //local_backbone_mover

#endif //INCLUDED_protocols/backbone_moves/local_backbone_mover_LocalBackboneMover_fwd_hh
