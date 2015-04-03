// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ResidualDipolarCoupling_Rohl.fwd.hh
/// @brief
/// @author Srivatsan Raman


#ifndef INCLUDED_core_scoring_ResidualDipolarCoupling_Rohl_fwd_hh
#define INCLUDED_core_scoring_ResidualDipolarCoupling_Rohl_fwd_hh

#include <utility/pointer/owning_ptr.hh>

//#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace scoring {

class ResidualDipolarCoupling_Rohl;
typedef utility::pointer::shared_ptr< ResidualDipolarCoupling_Rohl > ResidualDipolarCoupling_RohlOP;
typedef utility::pointer::shared_ptr< ResidualDipolarCoupling_Rohl const > ResidualDipolarCoupling_RohlCOP;


class RDC_Rohl;


} // scoring
} // core

#endif
