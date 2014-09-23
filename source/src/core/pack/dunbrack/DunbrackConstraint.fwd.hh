// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/dunbrack/DunbrackConstraint.fwd.hh
/// @brief
/// @author James Thompson


#ifndef INCLUDED_core_pack_dunbrack_DunbrackConstraint_fwd_hh
#define INCLUDED_core_pack_dunbrack_DunbrackConstraint_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace dunbrack {


class DunbrackConstraint; // fwd declaration
typedef utility::pointer::shared_ptr< DunbrackConstraint > DunbrackConstraintOP;
typedef utility::pointer::shared_ptr< DunbrackConstraint const > DunbrackConstraintCOP;


} // namespace constraints
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_DunbrackConstraint_FWD_HH
