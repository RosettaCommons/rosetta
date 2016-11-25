// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/symmetric_docking/membrane/MPSymDockMover.fwd.hh
///
/// @brief      Membrane symmetric docking protocol
/// @details    Given an asymmetric starting pose, setup the pose for symmetry,
///             add a membrane representation for the full symmetric complex,
///             position the pose at the center of mass of all transmembrane spans,
///             and proceed with the standard symmetric docking protocol.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 2/9/15


#ifndef INCLUDED_protocols_symmetric_docking_MPSymDockMoverCreator_hh
#define INCLUDED_protocols_symmetric_docking_MPSymDockMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace symmetric_docking {
namespace membrane {

/// @brief Mover Creator
class MPSymDockMoverCreator : public protocols::moves::MoverCreator {

public:

	// XRW TEMP  virtual protocols::moves::MoverOP create_mover() const;
	// XRW TEMP  virtual std::string keyname() const;
	// XRW TEMP  static std::string mover_name();
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} // membrane
} // symmetric_docking
} // protocols

#endif // INCLUDED_protocols_symmetric_docking_MPSymDockMoverCreator_hh
