// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loop_mover/refine/LoopMover_BackrubCreator.hh
/// @brief  Header for LoopMover_BackrubCreator
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_loops_loop_mover_refine_LoopMover_BackrubCreator_hh
#define INCLUDED_protocols_loops_loop_mover_refine_LoopMover_BackrubCreator_hh

// Unit Headers
#include <protocols/moves/MoverCreator.hh>
namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

/// @brief creator for the LoopMover_Refine_BackrubCreator class
class LoopMover_Refine_BackrubCreator : public moves::MoverCreator
{
public:

	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} //namespace refine
} //namepace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_refine_LoopMover_BackrubCreator_hh
