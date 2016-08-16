// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/dunbrack/DunbrackConstraintCreator.hh
/// @brief  Base class for ConstraintCreators for the Constraint load-time factory registration scheme
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_dunbrack_DunbrackConstraintCreator_hh
#define INCLUDED_core_pack_dunbrack_DunbrackConstraintCreator_hh

// Unit Headers
#include <core/scoring/constraints/ConstraintCreator.hh>

#include <core/scoring/constraints/Constraint.fwd.hh>

// c++ headers

namespace core {
namespace pack {
namespace dunbrack {

/// @brief Mover creator for the DunbrackConstraint constraint
class DunbrackConstraintCreator : public scoring::constraints::ConstraintCreator
{
public:
	DunbrackConstraintCreator();
	virtual ~DunbrackConstraintCreator();

	virtual scoring::constraints::ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};


} //namespace dunbrack
} //namespace pack
} //namespace core

#endif
