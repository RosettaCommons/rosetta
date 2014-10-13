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
/// @author Rhiju Das

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
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>

// Utility Headers
#include <utility/io/izstream.hh>

static basic::Tracer TR( "core.scoring.rna.RNA_SuitePotential", basic::t_info );

using namespace boost::numeric;

///////////////////////////////////////////////////////////////////////////////////////////
//
// This now has two 'modes' for computing energies based on seven
//  torsions t = (delta, epsilon, zeta, alpha, beta, gamma, delta):
//
// (1) Suite_potential ('likelihood based') -- computes:
//
//    E( t ) =  -log[  w_1a * G_1a( t ) + w_1c * G_1c( t ) + ... ]
//
//    See below ('Mahalanobis distance') for functional form of each G( t ).
//    By default E( t ) is zero when torsions are right 'on rotamer' and
//    can jump up to a very high value off-rotamer. Note that the 'weights',
//    which can be specified in files through -rna_suite_potential, enter by
//    their logarithms, as they are assumed to reflect statistical populations.
//
// (2) Suiteness bonus  -- computes:
//
//    E( t ) = w_1a * suiteness_1a( t ) + w_1c * suiteness_1c( t ) + ...
//
//   Here w_1a, etc. are simple linear coefficients on suiteness, which is calculated with
//    the cosine of 3-euclidean distance form due to the Richardsons (which is nice and flat
//    when on rotamer). To help reduce code copying, this makes use of the
//    RNA_SuiteName framework.
//
// The choice for putting both modes into here was to avoid copying a bunch of code for
//  packing torsions, deciding which torsions to score, and computing derivatives.
//
///////////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {

RNA_SuitePotential::~RNA_SuitePotential() {}

RNA_SuitePotential::RNA_SuitePotential( bool const calculate_suiteness_bonus /* = false */ ):
	n_torsions_( 7 ),
	inv_cov_( n_torsions_, n_torsions_ ),
	offset_( 0.0 ),
	score_( 0 ),
	deriv_( n_torsions_, 0 ),
	calculate_suiteness_bonus_( calculate_suiteness_bonus )
{
	using namespace basic::options;

	std::string path;
	if ( calculate_suiteness_bonus_ ) {
		path = "scoring/rna/suiteness_bonus/" +	option[ OptionKeys::score::suiteness_bonus ]() /* default is Richardson*/;
	} else {
		path = "scoring/rna/suite_potentials/" +	option[ OptionKeys::score::rna_suite_potential ]() /* default is Richardson*/;
	}
	utility::io::izstream stream;
	std::string line;
	Real weight, angle;

	// Load the suite centers and the corresponding weights
	std::string const filename_centers =
			basic::database::full_name( path + "/centers.txt" );
	stream.open( filename_centers );
	if ( !stream.good() ) utility_exit_with_message( "Trouble with path: " + path + "/centers.txt" );

	Size new_rotamer( 0 ); // just for anything without a tag.
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

		std::string tag( "" );
		if ( !line_stream.fail() ){
			line_stream >> tag;
		}
		if ( tag.size() == 0 ) tag = "new"+ObjexxFCL::string_of( ++new_rotamer);
		runtime_assert( !tags_.has_value( tag ) );
		tags_.push_back( tag );
	}

	// offset_ adjustment was too confusing, as we optimize potentials. -- rhiju, 2014
	//figure_out_offset();

	// Load information on how 'wide' the basins are.
	if (calculate_suiteness_bonus_ ) {
		// RNA_SuiteName contains info on half_widths for regular, dominant, and satellite half widths.
		// It will be used to actually calculate suiteness.
		rna_suite_name_ = pose::rna::RNA_SuiteNameOP( new pose::rna::RNA_SuiteName );

		// But need to update Suite centers if user has specified differently in suite_bonus file. Do this based on *tags*.
		// have to convert from ublas::vector to utility::vector1
		utility::vector1< utility::vector1< Real > >  centers_vector1;
		for ( Size n = 1; n <= centers_.size(); n++ ){
			utility::vector1< Real > center_vector1;
			for ( Size k = 1; k <= n_torsions_; k++ )	center_vector1.push_back( centers_[n](k-1) );
			centers_vector1.push_back( center_vector1 );
		}
		rna_suite_name_->update_centers( centers_vector1, tags_); // could relax this. just need some way to input information.

	} else {
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

	utility::vector1<Real> torsions;
	for ( Size i = 1; i <= torsion_ids_.size(); ++i ) {
		// A suite must have all its torsions being valid
		if ( !is_torsion_valid( pose, torsion_ids_[i] ) ) return false;
		Real const tor( pose.torsion( torsion_ids_[i] ) );
		torsions.push_back( tor );
	}


	eval_score( torsions );
	//	TR << rsdnum1 << "--" << rsdnum2 << ": " << score_ << std::endl;
	return true;
}

////////////////////////////////////////////////////////
void RNA_SuitePotential::eval_score(
	utility::vector1<Real> const & torsions
) const {

	if ( calculate_suiteness_bonus_ ) {
		eval_suiteness_bonus( torsions );
	} else {
		eval_likelihood_potential( torsions );
	}

}

////////////////////////////////////////////////////////
void RNA_SuitePotential::eval_suiteness_bonus(
	utility::vector1<Real> const & torsions
) const {

	pose::rna::RNA_SuiteAssignment assignment( rna_suite_name_->assign( torsions, deriv_ ) );

	if ( !tags_.has_value( assignment.name ) ) return; // outlier
	Real const bonus_weight = weights_[ tags_.index( assignment.name ) ];

	score_ = bonus_weight * assignment.suiteness;
	for ( Size n = 1; n <= deriv_.size(); n++ ) deriv_[ n ] *= bonus_weight;

}

////////////////////////////////////////////////////////
void RNA_SuitePotential::eval_likelihood_potential(
	utility::vector1<Real> const & torsions_in
) const {

	ublas::vector<Real> torsions( n_torsions_);
	for ( Size n = 1; n <= n_torsions_; n++ ) {
		torsions( n-1 ) = torsions_in[ n ]; // ublas vector starts with 0
	}
	utility::vector1<Real> unweight_likelihood;
	utility::vector1< ublas::vector<Real> > neg_likelihood_deriv;

	// Compute the likelihood and derivative for each center
	for ( Size i = 1; i <= centers_.size(); ++i ) {
		ublas::vector<Real> diff( torsions - centers_[i] );
		regularize_torsions( diff );
		ublas::vector<Real> const temp( ublas::prod( inv_cov_, diff ) );
		Real const likelihood0( exp( -ublas::inner_prod( diff, temp ) ) );
		ublas::vector<Real> const deriv0( 2 * likelihood0 * temp );
		unweight_likelihood.push_back( likelihood0 );
		neg_likelihood_deriv.push_back( deriv0 );
	}

	Real weighted_sum_likelihood( 0 );
	for ( Size i = 1; i <= unweight_likelihood.size(); ++i ) {
		//		TR << "comp: " << i << " " << weights_[i] << " " << unweight_likelihood[i] << std::endl;
		weighted_sum_likelihood += weights_[i] * unweight_likelihood[i];
	}
	score_ = -log( weighted_sum_likelihood ) + offset_;

	ublas::vector<Real> deriv_vec( n_torsions_ );
	deriv_vec *= 0;  // Initialize to 0

	for ( Size i = 1; i <= neg_likelihood_deriv.size(); ++i ) {
		deriv_vec += weights_[i] * neg_likelihood_deriv[i];
	}
	deriv_vec /= weighted_sum_likelihood;
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

////////////////////////////////////////////////
// Set the offset to make the minimal score to be roughly 0.
// Note: this is no longer activated. OK to delete this function
// and offset_ variable if still not in use in 2015.
void
RNA_SuitePotential::figure_out_offset() {
	Real max( 0 );
	for ( Size i = 1; i <= weights_.size(); ++i ) {
		if ( weights_[i] > max ) max = weights_[i];
	}
	offset_ = log( max );
}

} //rna
} //scoring
} //core
