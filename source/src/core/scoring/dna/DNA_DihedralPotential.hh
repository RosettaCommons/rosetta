// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/dna/DNA_DihedralPotential.hh
/// @brief  dna scoring
/// @author Phil Bradley

#ifndef INCLUDED_core_scoring_dna_DNA_DihedralPotential_HH
#define INCLUDED_core_scoring_dna_DNA_DihedralPotential_HH

//#include <core/scoring/dna/DNABFormPotential.hh>

#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.fwd.hh>
#include <utility/vector0.fwd.hh>
#include <utility/exit.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <iosfwd>
#include <map>

namespace core {
namespace scoring {
namespace dna {

class DNA_DihedralPotential {
public:
	DNA_DihedralPotential( std::string const & filename );

	DNA_DihedralPotential();

	void
	eval_harmonic_backbone_torsion_score_and_deriv(
		Size const tor,
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		Real & score,
		Real & dscore_dtor
	) const;

	void
	eval_harmonic_sugar_pucker_dependent_chi_torsion_score_and_deriv(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		Size const pucker,
		Real & score,
		Real & dscore_dchi
	) const;

	Real
	get_mean_sugar_pucker_dependent_chi( conformation::Residue const & rsd ) const;

	void
	eval_sugar_torsion_score_and_deriv(
		Real const torsion,
		Size const tor,
		conformation::Residue const & rsd,
		Size const pucker,
		Real & score,
		Real & dscore_dtor
	) const;

	void
	get_sugar_torsion_mean_and_sdev(
		Size const tor,
		conformation::Residue const & rsd,
		Size const pucker,
		Real & mean,
		Real & sdev
	) const;

private:

	bool
	skip_torsion( conformation::Residue const & rsd, Size const tor ) const;

	void
	read_dna_geometry_log_file(
		std::string const & filename
	);

	void
	read_dna_geometry_log_file_from_database(
		std::string const & database_file
	);

	void
	parse_dna_geometry_log( std::istream & data );

	// this is not still around; to see the code (including smoothed torsion probs) use git, pre-2016: phbradley/mergetest
	// void
	// compute_torsion_probabilities(); // after reading the database_file

private:

	// utility::vector1< utility::vector1< utility::vector1< Real > > > torsion_probability_;
	// utility::vector1< utility::vector1< utility::vector1< Real > > > dtorsion_probability_;

	utility::vector1< utility::vector1< Real > > mean_backbone_torsion_;

	//utility::vector1< utility::vector1< utility::vector1< Real > > > mean_sugar_torsion_;
	utility::vector1< utility::vector0< utility::vector1< Real > > > mean_sugar_torsion_;

	//DNABFormPotential jjh_potential_;
};



} // ns dna
} // ns scoring
} // ns core

#endif
