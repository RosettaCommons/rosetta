// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: LoopModeler
///
///
/// @author Vatsan Raman


#ifndef INCLUDED_devel_residual_dipolar_coupling_RDC_energies_hh
#define INCLUDED_devel_residual_dipolar_coupling_RDC_energies_hh

// Package headers
#include <devel/residual_dipolar_coupling/RDC_main.hh>

//Core
#include <core/types.hh>

//Objexx headers
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

//// C++ headers
#include <string>
#include <map>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace devel {
namespace residual_dipolar_coupling {

void eval_dipolar(
	core::pose::Pose const & pose,
	utility::vector1<devel::residual_dipolar_coupling::RDC> const & All_RDC_lines
);

void assemble_datamatrix(
	core::pose::Pose const & pose,
	utility::vector1<devel::residual_dipolar_coupling::RDC> const & All_RDC_lines,
	ObjexxFCL::FArray2D< core::Real > & A,
	ObjexxFCL::FArray1D< core::Real > & b
);

void calc_ordermatrix(
	core::Size const & nrow,
	core::Size const & ORDERSIZE,
	ObjexxFCL::FArray2D< core::Real > & A,
	ObjexxFCL::FArray1D< core::Real > & b,
	ObjexxFCL::FArray1D< core::Real > & x,
	bool & reject
);

void svdcmp(
	ObjexxFCL::FArray2D< core::Real > & a,
	core::Size const & m,
	core::Size const & n,
	ObjexxFCL::FArray1D< core::Real > & w,
	ObjexxFCL::FArray2D< core::Real > & v
);

core::Real pythag(
	core::Real const & a,
	core::Real const & b
);

void svbksb(
	ObjexxFCL::FArray2D< core::Real > const & u,
	ObjexxFCL::FArray1D< core::Real > const & w,
	ObjexxFCL::FArray2D< core::Real > const & v,
	core::Size const & m,
	core::Size const & n,
	ObjexxFCL::FArray1D< core::Real > const & b,
	ObjexxFCL::FArray1D< core::Real > & x
);

void calc_orderparam(
	ObjexxFCL::FArray1D< core::Real > x,
	ObjexxFCL::FArray2D< core::Real > vec,
	core::Real & Azz,
	core::Real & eta
);

void calc_dipscore(
	ObjexxFCL::FArray2D< core::Real > const & A,
	ObjexxFCL::FArray1D< core::Real > const & x,
	ObjexxFCL::FArray1D< core::Real > const & b,
	utility::vector1<devel::residual_dipolar_coupling::RDC> const & All_RDC_lines,
	core::Size const & ORDERSIZE,
	core::Real const & Azz
);

} //ResidualDipolarCoupling
} //devel

#endif
