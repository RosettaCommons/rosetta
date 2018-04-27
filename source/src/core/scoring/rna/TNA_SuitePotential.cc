// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/TNA_SuitePotential.cc
/// @brief  TNA_SuitePotential potential class implementation
/// @author Andy Watkins

// Unit Headers
#include <core/scoring/rna/TNA_SuitePotential.hh>

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

static basic::Tracer TR( "core.scoring.rna.TNA_SuitePotential", basic::t_info );

using namespace boost::numeric;
using namespace core::chemical;

///////////////////////////////////////////////////////////////////////////////////////////
//
// Unlike the RNA suite energy, this only has one 'mode' for computing
// energies and RIGHT NOW only supports four torsions. While an RNA suite is
//     (delta, epsilon, zeta, alpha, beta, gamma, delta)
// we have little understanding of TNA and will momentarily assume that ring
// flips -- which would be ludicrous polymerically because of the way the
// polymeric connection is 2,3-anti -- don't happen. So delta/'ring torsion'
// doesn't matter. Furthermore, there is one fewer dihedral anyway. So we
// just have:
//        t = (epsilon, zeta, alpha, beta)
//
// suite_potential ('likelihood based') -- computes:
//
//    E( t ) =  -log[  w_1a * G_1a( t ) + w_1c * G_1c( t ) + ... ]
//
//    See below ('Mahalanobis distance') for functional form of each G( t ).
//    By default E( t ) is zero when torsions are right 'on rotamer' and
//    can jump up to a very high value off-rotamer. Note that the 'weights',
//    which can be specified in files through -rna_suite_potential, enter by
//    their logarithms, as they are assumed to reflect statistical populations.
//
///////////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {

TNA_SuitePotential::~TNA_SuitePotential() = default;

/// @details This constructor reads in data from disk and should only 
/// be called from the ScoringManager
TNA_SuitePotential::TNA_SuitePotential() // TODO move to in class init
{
	using namespace basic::options;

	std::string path;
	path = "scoring/rna/suite_potentials/tna_suite_potential";
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
		if ( !line_stream.fail() ) {
			line_stream >> tag;
		}
		if ( tag.size() == 0 ) tag = "new"+ObjexxFCL::string_of( ++new_rotamer);
		runtime_assert( !tags_.has_value( tag ) );
		tags_.push_back( tag );
	}

	// offset_ adjustment was too confusing, as we optimize potentials. -- rhiju, 2014

	// Load information on how 'wide' the basins are.

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
bool
TNA_SuitePotential::eval_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose
) const {
#ifdef MULTI_THREADED
	utility_exit_with_message( "The TNA_SuitePotential caches pose-specific scoring data in the global instance of the object.  As such, it is fundamentally non-threadsafe, and cannot be used in a multi-threaded context.  Please contact Rhiju Das if you need to use this scoring term in a multi-threaded context." );
#endif

	using namespace core::chemical::rna;
	using namespace core::pose::rna;
	using namespace core::id;

	if ( rsd1.has_variant_type( REPLONLY ) ) return false;
	if ( rsd2.has_variant_type( REPLONLY ) ) return false;

	// Default derivatives to 0
	score_ = 0;
	for ( Size i = 1; i <= n_torsions_; ++i ) deriv_[i] = 0;

	// Only counts consecutive TNA residues.
	if ( !rsd1.is_TNA() ) return false;
	if ( !rsd2.is_TNA() ) return false;

	TR.Debug << "Computing suite potential for " << rsd1.seqpos() << " " << rsd2.seqpos() << std::endl;
	// AMW TODO:
	// 1. Compute based on bonding.
	// 2. Don't cache anything.
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
	//torsion_ids_.emplace_back( rsdnum1, BB, DELTA );
	torsion_ids_.emplace_back( rsdnum1, BB, DELTA/*EPSILON*/ );
	torsion_ids_.emplace_back( rsdnum1, BB, EPSILON/*ZETA*/ );
	torsion_ids_.emplace_back( rsdnum2, BB, ALPHA );
	torsion_ids_.emplace_back( rsdnum2, BB, BETA );
	//torsion_ids_.emplace_back( rsdnum2, BB, GAMMA );
	//torsion_ids_.emplace_back( rsdnum2, BB, DELTA );

	utility::fixedsizearray1<Real, 4> torsions;
	for ( Size i = 1; i <= torsion_ids_.size(); ++i ) {
		// A suite must have all its torsions being valid
		// Getting NaNs if I don't skip chainbreak torsions; this is kinda
		// fine with me really.
		if ( !is_torsion_valid( pose, torsion_ids_[i], false, true ) ) return false;
		torsions[ i ] = pose.torsion( torsion_ids_[i] );
	}

	eval_likelihood_potential( torsions );
	// TR << rsdnum1 << "--" << rsdnum2 << ": " << score_ << std::endl;
	return true;
}

////////////////////////////////////////////////////////
void TNA_SuitePotential::eval_likelihood_potential(
	utility::fixedsizearray1<Real, 4> const & torsions_in
) const {
#ifdef MULTI_THREADED
	utility_exit_with_message( "The TNA_SuitePotential caches pose-specific scoring data in the global instance of the object.  As such, it is fundamentally non-threadsafe, and cannot be used in a multi-threaded context.  Please contact Rhiju Das if you need to use this scoring term in a multi-threaded context." );
#endif
	TR.Debug << "Looking at torsions: " << torsions_in << std::endl;
	ublas::vector<Real> torsions( n_torsions_);
	for ( Size n = 1; n <= n_torsions_; n++ ) {
		torsions( n-1 ) = torsions_in[ n ]; // ublas vector starts with 0
	}
	utility::vector1<Real> unweight_likelihood;
	utility::vector1< ublas::vector<Real> > neg_likelihood_deriv;

	// Compute the likelihood and derivative for each center
	for ( auto const & center : centers_ ) {
		ublas::vector<Real> diff( torsions - center );
		regularize_torsions( diff );
		TR.Trace << "Torsions center diff: " << diff(0) << " " << diff(1) << " " << diff(2) << " " << diff(3) << std::endl;
		ublas::vector<Real> const temp( ublas::prod( inv_cov_, diff ) );
		Real const likelihood0( exp( -ublas::inner_prod( diff, temp ) ) );
		ublas::vector<Real> const deriv0( 2 * likelihood0 * temp );
		unweight_likelihood.push_back( likelihood0 );
		neg_likelihood_deriv.push_back( deriv0 );
	}

	Real weighted_sum_likelihood( 0 );
	for ( Size i = 1; i <= unweight_likelihood.size(); ++i ) {
		//  TR << "comp: " << i << " " << weights_[i] << " " << unweight_likelihood[i] << std::endl;
		weighted_sum_likelihood += weights_[i] * unweight_likelihood[i];
	}
	score_ = -log( weighted_sum_likelihood ) + offset_;
	TR.Debug << "This residue pair score: " << score_ << std::endl;
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
void TNA_SuitePotential::regularize_torsions(
	ublas::vector<Real> & torsions
) const {
#ifdef MULTI_THREADED
	utility_exit_with_message( "The TNA_SuitePotential caches pose-specific scoring data in the global instance of the object.  As such, it is fundamentally non-threadsafe, and cannot be used in a multi-threaded context.  Please contact Rhiju Das if you need to use this scoring term in a multi-threaded context." );
#endif

	for ( Size i = 0; i != torsions.size(); ++i ) {
		torsions(i) = numeric::principal_angle_degrees( torsions(i) );
	}
}

////////////////////////////////////////////////
// Set the offset to make the minimal score to be roughly 0.
// Note: this is no longer activated. OK to delete this function
// and offset_ variable if still not in use in 2015.
void
TNA_SuitePotential::figure_out_offset() {
#ifdef MULTI_THREADED
	utility_exit_with_message( "The TNA_SuitePotential caches pose-specific scoring data in the global instance of the object.  As such, it is fundamentally non-threadsafe, and cannot be used in a multi-threaded context.  Please contact Rhiju Das if you need to use this scoring term in a multi-threaded context." );
#endif

	Real max( 0 );
	for ( Size i = 1; i <= weights_.size(); ++i ) {
		if ( weights_[i] > max ) max = weights_[i];
	}
	offset_ = log( max );
}

} //rna
} //scoring
} //core
