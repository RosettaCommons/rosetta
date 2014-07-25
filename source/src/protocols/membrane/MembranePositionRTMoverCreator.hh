// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		protocols/membrane/MembranePositionRTMoverCreator.hh
///
/// @brief		Membrane Position Rotation/Translation Mover Creator
///	@details	Apply a uniform rigid translation & rotation of the
///				membrane defined by the center/normal coordinates
///				stored in the membrane.
///				Last Modified: 6/28/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_MembranePositionRTMoverCreator_hh
#define INCLUDED_protocols_membrane_MembranePositionRTMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {

/// @brief Mover Creator for composed RT mover
class MembranePositionRTMoverCreator : public protocols::moves::MoverCreator {
	
public:
	
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
	
};

/// @brief Mover creator for membrane rotation mover
class MembranePositionRotationMoverCreator : public protocols::moves::MoverCreator {
	
public:
	
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
	
};

/// @brief Mover Creator for membrane translation mover
class MembranePositionTranslationMoverCreator : public protocols::moves::MoverCreator {
	
public:
	
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
	
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MembranePositionRTMoverCreator_hh
