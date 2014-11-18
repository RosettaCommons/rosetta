// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief      Does translation and rotation of a pose or the membrane
/// @details	Translates along a vector and rotates a pose or membrane around
///				the axis perpendicular to both vectors. Works for both
///				fixed protein/flexible membrane or fixed membrane/flexible protein
///				Takes the jump number as input, default is the membrane jump
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_TranslationRotationMoverCreator_hh
#define INCLUDED_protocols_membrane_TranslationRotationMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {

/// @brief Mover Creator
class TranslationMoverCreator : public protocols::moves::MoverCreator {
public:
	
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
	
};

/// @brief Mover Creator
class RotationMoverCreator : public protocols::moves::MoverCreator {
public:
	
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
	
};

/// @brief Mover Creator
class TranslationRotationMoverCreator : public protocols::moves::MoverCreator {	
public:
	
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
	
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TranslationRotationMoverCreator_hh
