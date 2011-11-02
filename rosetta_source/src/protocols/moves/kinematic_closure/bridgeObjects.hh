// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   bridgeObjects.hh
/// @brief  Header file for bridgeObjects code for ALC
/// @author Evangelos A. Coutsias
/// @author Daniel J. Mandell

#ifndef INCLUDE_src_protocols_moves_kinematic_closure_bridgeObjects_hh
#define INCLUDE_src_protocols_moves_kinematic_closure_bridgeObjects_hh

// Rosetta Headers
#include <core/types.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>

#include <utility/vector1.fwd.hh>


namespace protocols {
namespace moves {
namespace kinematic_closure {

void bridgeObjects (const utility::vector1<utility::vector1<core::Real> >& atoms,
					const utility::vector1<core::Real> & dt, const utility::vector1<core::Real> & da,
					const utility::vector1<core::Real> & db, const utility::vector1<int>& pivots,
					const utility::vector1<int>& order, utility::vector1<utility::vector1<core::Real> >& t_ang,
					utility::vector1<utility::vector1<core::Real> >& b_ang,
					utility::vector1<utility::vector1<core::Real> >& b_len, int& nsol);
void chainTORS (const int& n, const utility::vector1<utility::vector1<core::Real> >& atoms,
				utility::vector1<core::Real>& t_ang, utility::vector1<core::Real>& b_ang,
				utility::vector1<core::Real>& b_len, utility::vector1<core::Real>& R0,
				utility::vector1<utility::vector1<core::Real> >& Q);
core::Real torsion(const utility::vector1<core::Real>& a, const utility::vector1<core::Real>& b,
				   const utility::vector1<core::Real>& c, const utility::vector1<core::Real>& d);
void chainXYZ  (const int& n,
				const utility::vector1<core::Real>& b_len,
				const utility::vector1<core::Real>& b_ang,
				const utility::vector1<core::Real>& t_ang,
				const bool space,
				const utility::vector1<core::Real>& R0,
				const utility::vector1<utility::vector1<core::Real> >& Q,
				utility::vector1<utility::vector1<core::Real> >& atoms
				);
core::Real bondangle(const utility::vector1<core::Real>& a, const utility::vector1<core::Real>& b,
					 const utility::vector1<core::Real>& c);

} // end namespace kinematic_closure
} // end namespace moves
} // end namespace protocols

#endif
