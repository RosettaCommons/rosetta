// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/RNA_SuitePotential.cc
/// @brief  RNA_SuitePotential potential class implementation
/// @author Fang-Chieh Chou

// Unit Headers
#include <core/scoring/rna/RNA_SuitePotential.hh>

// Package Headers

// Project Headers
#include <core/chemical/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/conformation/Residue.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>

// Utility Headers
#include <utility/io/izstream.hh>

using namespace boost::numeric;

namespace core {
namespace scoring {
namespace rna {

RNA_SuitePotential::~RNA_SuitePotential() {}

RNA_SuitePotential::RNA_SuitePotential():
	n_torsions_( 7 ),
	inv_cov_( n_torsions_, n_torsions_ ),
	score_( 0 ),
	deriv_( n_torsions_, 0 )
{
	using namespace basic::options;

	std::string path;
	if ( option[ OptionKeys::score::rna_suite_potential ].user() ){
		path = "scoring/rna/suite_potentials/" +
			option[ OptionKeys::score::rna_suite_potential ]();
	} else {
		path = "scoring/rna/suite_potentials/Richardson";
	}

	utility::io::izstream stream;
	std::string line;
	Real weight, angle;

	// Load the suite centers and the corresponding weights
	std::string const filename_centers =
			basic::database::full_name( path + "/centers.txt" );
	stream.open( filename_centers );

	while ( getline( stream, line ) ) {
		std::istringstream line_stream( line );
		line_stream >> weight;
		weights_.push_back( weight );
		ublas::vector<Real> center( n_torsions_ );
		for ( Size i = 0; i != n_torsions_; ++i ) {
			line_stream >> angle;
			center( i ) = angle;
		}
		centers_.push_back( center );
	}

	// Set the offset to make the minimal score to be roughly 0.
	Real max( 0 );
	for ( Size i = 1; i <= weights_.size(); ++i ) {
		if ( weights_[i] > max ) max = weights_[i];
	}
	offset_ = log( max );

	// Load the inverse covariance
	std::string const filename_invcov =
			basic::database::full_name( path + "/inv_cov.txt" );
	stream.open( filename_invcov );

	for ( Size i = 0; i != n_torsions_; ++i ) {
		getline( stream, line );
		std::istringstream line_stream( line );
		for ( Size j = 0; j != n_torsions_; ++j ) {
			line_stream >> angle;
			inv_cov_( i, j ) = angle;
		}
	}
}

// Compute the score and derivatives
// Score = -log( w1 * exp(-sq_d1) + w2 * exp(-sq_d2) + ... )
// Where w1, w2... are weights_[i]
// sq_d1, sq_d2... are
// (torsion - centers_[i]).T.dot(inv_cov_).dot(torsion - centers_[i])
// which is the mahalanobis distance between the torsions and the centers
// Return false if there is no valid suite.
bool RNA_SuitePotential::eval_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose
) const {
	using namespace core::chemical::rna;
	using namespace core::pose::rna;
	using namespace core::id;
	// Default derivatives to 0
	score_ = 0;
	for ( Size i = 1; i <= n_torsions_; ++i ) deriv_[i] = 0;

	// Only counts consecutive RNA residues.
	if ( !rsd1.is_RNA() ) return false;
	if ( !rsd2.is_RNA() ) return false;
	Size rsdnum1( rsd1.seqpos() ), rsdnum2( rsd2.seqpos() );
	if ( rsdnum1 + 1 == rsdnum2 ) {
		;  // Good: this is what's assumed in the following calculation
	} else if ( rsdnum1 == rsdnum2 + 1 ) {
		std::swap( rsdnum1, rsdnum2 );  // Swap rsdnum1 and rsdnum2
	} else {
		return false;  // No valid suite
	}

	// Find the suite torsions
	torsion_ids_.clear();
	torsion_ids_.push_back( TorsionID( rsdnum1, BB, DELTA ) );
	torsion_ids_.push_back( TorsionID( rsdnum1, BB, EPSILON ) );
	torsion_ids_.push_back( TorsionID( rsdnum1, BB, ZETA ) );
	torsion_ids_.push_back( TorsionID( rsdnum2, BB, ALPHA ) );
	torsion_ids_.push_back( TorsionID( rsdnum2, BB, BETA ) );
	torsion_ids_.push_back( TorsionID( rsdnum2, BB, GAMMA ) );
	torsion_ids_.push_back( TorsionID( rsdnum2, BB, DELTA ) );

	ublas::vector<Real> torsions( n_torsions_ );
	for ( Size i = 1; i <= torsion_ids_.size(); ++i ) {
		// A suite must have all its torsions being valid
		if ( !is_torsion_valid( pose, torsion_ids_[i] ) ) return false;
		Real const tor( pose.torsion( torsion_ids_[i] ) );
		torsions( i-1 ) = tor; // ublas vector starts with 0
	}

	eval_score( torsions );
	return true;
}

void RNA_SuitePotential::eval_score(
	ublas::vector<Real> const & torsions
) const {
	utility::vector1<Real> unweight_likelyhood;
	utility::vector1< ublas::vector<Real> > neg_likelyhood_deriv;

	// Compute the likelyhood and derivative for each center
	for ( Size i = 1; i <= centers_.size(); ++i ) {
		ublas::vector<Real> diff( torsions - centers_[i] );
		regularize_torsions( diff );
		ublas::vector<Real> const temp( ublas::prod( inv_cov_, diff ) );
		Real const likelyhood0( exp( -ublas::inner_prod( diff, temp ) ) );
		ublas::vector<Real> const deriv0( 2 * likelyhood0 * temp );
		unweight_likelyhood.push_back( likelyhood0 );
		neg_likelyhood_deriv.push_back( deriv0 );
	}

	Real weighted_sum_likelyhood( 0 );
	for ( Size i = 1; i <= unweight_likelyhood.size(); ++i ) {
		weighted_sum_likelyhood += weights_[i] * unweight_likelyhood[i];
	}
	score_ = -log( weighted_sum_likelyhood ) + offset_;

	ublas::vector<Real> deriv_vec( n_torsions_ );
	deriv_vec *= 0;  // Initialize to 0

	for ( Size i = 1; i <= neg_likelyhood_deriv.size(); ++i ) {
		deriv_vec += weights_[i] * neg_likelyhood_deriv[i];
	}
	deriv_vec /= weighted_sum_likelyhood;
	for ( Size i = 1; i <= deriv_vec.size(); ++i ) {
		deriv_[i] = deriv_vec( i-1 ); // Convert to vector1
	}
}

// Folds the torsions into (-180, 180]
void RNA_SuitePotential::regularize_torsions(
	ublas::vector<Real> & torsions
) const {
	for ( Size i = 0; i != torsions.size(); ++i ) {
		torsions(i) = numeric::principal_angle_degrees( torsions(i) );
	}
}

}
}
}
