// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rna/RNA_SuiteName.hh
/// @brief RNA suite assignment ported from suitename program
/// @author Fang-Chieh Chou

// NOTE: The definition of suite numbering here follows the original
// definition of the Richardson's. Suite i consists of the delta, epsilon,
// zeta of residue (i-1) and alpha, beta, gamma, delta of residue i. This is
// different to the one being used in current SWA codes.

#ifndef INCLUDED_core_pose_rna_RNA_SuiteName_HH
#define INCLUDED_core_pose_rna_RNA_SuiteName_HH

// Unit headers
#include <core/pose/rna/RNA_SuiteName.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace core {
namespace pose {
namespace rna {
//////////////////////////////////////////////////
class RNA_SuiteAssignment {
public:
	std::string name;
	Real suiteness;
	Real dist;

	RNA_SuiteAssignment(
		std::string const & name_in,
		Real const & suiteness_in,
		Real const & dist_in
	):
		name( name_in ),
		suiteness( suiteness_in ),
		dist( dist_in )
	{}

	RNA_SuiteAssignment():
		name( "" ),
		suiteness( -1 ),
		dist( -1 )
	{}
};
/////////////////////////////////////////
class RNA_SuiteInfo {
public:
	std::string name;
	Size classifier;
	utility::vector1<Real> torsion;

	RNA_SuiteInfo(
		std::string const & name_in,
		Size const & classifier_in,
		utility::vector1<Real> const & torsion_in
	):
		name( name_in ),
		classifier( classifier_in ),
		torsion( torsion_in )
	{}

	RNA_SuiteInfo():
		name( "" ),
		classifier( 0 )
	{}
};

////////////////////////////////////////
class RNA_SuiteName : public utility::pointer::ReferenceCount {
public:
	RNA_SuiteName();

	~RNA_SuiteName();

	RNA_SuiteInfo name2suite( std::string const & name ) const;

	RNA_SuiteAssignment assign(utility::vector1<Real> const & torsions_in) const;

	RNA_SuiteAssignment assign(utility::vector1<Real> const & torsions_in,
		utility::vector1<Real> & deriv ) const;

	RNA_SuiteAssignment assign(Pose const & pose, Size const res) const;

	void
	update_centers( utility::vector1< utility::vector1< Real > > const & centers,
		utility::vector1< std::string > const & tags );

	Real const
		epsilonmin, epsilonmax,
		delta3min , delta3max,
		delta2min , delta2max,
		gammapmin , gammapmax,
		gammatmin , gammatmax,
		gammammin , gammammax,
		alphamin  , alphamax,
		betamin   , betamax,
		zetamin   , zetamax;

	Real const delta_cutoff, gamma_cutoff1, gamma_cutoff2;

private:
	//Disable copy constructor and assignment
	RNA_SuiteName( const RNA_SuiteName & );
	void operator=( const RNA_SuiteName & );

	void init();

	Size
	get_classifier( utility::vector1< Real > const & torsions, bool & is_outlier ) const;

	Size
	get_classifier( utility::vector1< Real > const & torsions ) const;

	Real distance_4d(utility::vector1<Real> const &torsion1, utility::vector1<Real> const &torsion2,
		utility::vector1<Size> const & half_width) const;

	Real distance_7d(utility::vector1<Real> const &torsion1, utility::vector1 <Real> const &torsion2,
		utility::vector1<Size> const & half_width) const;

	Real distance_7d(
		utility::vector1<Real> const & torsion1,
		utility::vector1<Real> const & torsion2,
		utility::vector1<Size> const & half_width,
		utility::vector1<Real> & deriv // fill if non-empty
	) const;

	Real
	get_suiteness( Real const & dist_7d ) const;

	Real
	get_suiteness_derivative( Real const & dist_7d ) const;

	void
	fill_suiteness_derivative_7d (
		utility::vector1< Real  > & deriv,
		utility::vector1< Real > const & torsions,
		utility::vector1< Real > const & torsions_center,
		utility::vector1< Real > const & half_width ) const;

	bool is_in_between( utility::vector1<Real> const & target,
		utility::vector1<Real> const & dominant,
		utility::vector1<Real> const & satellite ) const;

	//////////////////////////////////
	RNA_SuiteAssignment const suite_undefined;
	Real const dist_pow_;
	utility::vector1<RNA_SuiteInfo> all_suites_;
	utility::vector1<Size> regular_half_width;
	utility::vector1<std::string> dominant_suites;
	utility::vector1<std::string> satellite_suites;
	utility::vector1< utility::vector1<Size> > half_width_sat;
	utility::vector1< utility::vector1<Size> > half_width_dom;
};

}
}
}

#endif
