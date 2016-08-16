// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_dna_base_geometry_hh
#define INCLUDED_core_scoring_dna_base_geometry_hh

#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <iosfwd>

#include <utility/vector1.hh>


#ifdef WIN32
#include <string>
#endif

namespace core {
namespace scoring {
namespace dna {

void
get_y_axis_atoms(
	chemical::ResidueType const & rsd_type,
	int const strand, // 1 or 2
	std::string & a1,
	std::string & a2
);


bool
is_orthonormal(
	numeric::xyzMatrix< Real > const & M,
	Real const tol
);


void
get_sugar_pucker(
	conformation::Residue const & rsd,
	std::pair< std::string, int > & pucker
);

void
get_sugar_torsions(
	conformation::Residue const & rsd,
	utility::vector1< Real > & torsions
);

Vector
get_y_axis(
	conformation::Residue const & rsd,
	int const strand
);


Vector
get_z_axis(
	conformation::Residue const & rsd,
	Vector const & y_axis
);


kinematics::Stub
get_base_stub(
	conformation::Residue const & rsd,
	int const strand
);


void
get_base_pair_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	utility::vector1< Real > & params // output
);


void
show_base_pair_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
);


void
get_base_step_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	utility::vector1< Real > & params // output
);


void
show_new_base_step_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
);


void
show_dna_geometry( pose::Pose const & pose, std::ostream & out );

void
get_base_step_params(
	conformation::Residue const & rsd11, // pair1 strand I
	conformation::Residue const & rsd12, // pair1 strand II
	conformation::Residue const & rsd21, // pair2 strand I
	conformation::Residue const & rsd22, // pair2 strand II
	utility::vector1< Real > & params // output
);

void
get_base_pucker(
	conformation::Residue const & rsd,
	std::pair< std::string, int > & pucker,
	Real & pseudorotation,
	Real & amplitude
);

kinematics::Stub
get_base_pair_stub_slow(
	conformation::Residue const & rsd1, // on strand I
	conformation::Residue const & rsd2  // on strand II
);


kinematics::Stub
get_midstep_stub(
	kinematics::Stub const & in_stub1,
	kinematics::Stub const & in_stub2
);


} // namespace dna
}} // scoring core

#endif
