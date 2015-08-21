// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/CentroidRelax.fwd.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_relax_CENTROIDRELAX_FWD_HH
#define INCLUDED_protocols_relax_CENTROIDRELAX_FWD_HH

// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace relax {

// Forward
class CentroidRelax;

// Types
typedef  utility::pointer::shared_ptr< CentroidRelax >  CentroidRelaxOP;
typedef  utility::pointer::shared_ptr< CentroidRelax const >  CentroidRelaxCOP;

typedef  utility::pointer::weak_ptr< CentroidRelax >  CentroidRelaxAP;
typedef  utility::pointer::weak_ptr< CentroidRelax const >  CentroidRelaxCAP;


} // namespace kinematics
} // namespace core

#endif //#ifndef INCLUDED_protocols/relax_CENTROIDRELAX_FWD_HH
