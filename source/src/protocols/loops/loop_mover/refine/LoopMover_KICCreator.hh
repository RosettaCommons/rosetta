// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/LoopMover_KICCreator.hh
/// @brief  Header for LoopMover_KICCreator
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_loops_loop_mover_refine_LoopMover_KICCreator_hh
#define INCLUDED_protocols_loops_loop_mover_refine_LoopMover_KICCreator_hh

// Unit Headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

/// @brief creator for the LoopMover_Perturb_KICCreator class
class LoopMover_Refine_KICCreator : public moves::MoverCreator {
public:
	// XRW TEMP  LoopMover_Refine_KICCreator() {}
	// XRW TEMP  virtual ~LoopMover_Refine_KICCreator();

	// XRW TEMP  virtual moves::MoverOP create_mover() const;

	// XRW TEMP  virtual std::string keyname() const;
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} //namespace refine
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_refine_LoopMover_KICCreator_hh
