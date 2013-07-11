// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_LoopEnergy.cc
/// @brief  Cost of bringing two chains together.
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_LoopEnergy.hh>
#include <core/scoring/rna/RNA_LoopEnergyCreator.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++

static basic::Tracer tr("core.scoring.rna.RNA_LoopEnergy");

namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the RNA_LoopEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_LoopEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new RNA_LoopEnergy;
}

ScoreTypes
RNA_LoopEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_loop_fixed );
	sts.push_back( rna_loop_logN );
	sts.push_back( rna_loop_harmonic );
	return sts;
}


/// c-tor
RNA_LoopEnergy::RNA_LoopEnergy() :
	parent( new RNA_LoopEnergyCreator ),
	persistence_length2_( 8.0 * 8.0 )
{}

/// clone
methods::EnergyMethodOP
RNA_LoopEnergy::clone() const
{
	return new RNA_LoopEnergy;
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void
RNA_LoopEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	update_rna_loop_atoms_and_lengths( pose );
}

/////////////////////////////////////////////////////////////////////////////
void
RNA_LoopEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	update_rna_loop_atoms_and_lengths( pose );
}

/////////////////////////////////////////////////////////////////////////////
void
RNA_LoopEnergy::update_rna_loop_atoms_and_lengths( pose::Pose & pose ) const {

	using namespace core::pose::full_model_info;
	using namespace core::id;

	FullModelInfo const & full_model_info = const_full_model_info_from_pose( pose );
	utility::vector1< Size > const & sub_to_full = full_model_info.sub_to_full();
	utility::vector1< Size > const chains = figure_out_chains_from_full_model_info( pose );

	loop_takeoff_atoms_.clear();
	loop_landing_atoms_.clear();
	rna_loop_lengths_.clear();

	for ( Size n = 1; n < pose.total_residue(); n++ ){

		tr.Debug << "chain at " << n << " is " <<  chains[ n ] << std::endl;

		if ( pose.residue( n ).is_RNA() &&
				 pose.residue( n+1 ).is_RNA() &&
				 sub_to_full[ n+1 ] > sub_to_full[ n ] + 1 &&
				 chains[ n ] == chains[ n+1 ] ){

			// better not be a closed chainbreak!
			runtime_assert( ! pose.residue( n   ).has_variant_type( chemical::CUTPOINT_LOWER ) );
			runtime_assert( ! pose.residue( n+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) );

			/////////////////////////////////////////////////////////////////////
			// a single residue bulge has two suites
			Size num_loop_suites = sub_to_full[ n+1 ] - sub_to_full[ n ];
			runtime_assert( num_loop_suites > 1 );

			rna_loop_lengths_.push_back( num_loop_suites );
			loop_takeoff_atoms_.push_back( AtomID( pose.residue( n ).atom_index( " O3*" )  , n ) );
			loop_landing_atoms_.push_back( AtomID( pose.residue( n+1 ).atom_index( " C5*" ), n+1 ) );

			tr.Debug << "Found RNA loop: " << AtomID( pose.residue( n ).atom_index( " O3*" )  , n ) << " -- " <<  AtomID( pose.residue( n+1 ).atom_index( " C5*" ), n+1 ) << std::endl << "  loop length: " << num_loop_suites;

		}

	}

	runtime_assert( loop_takeoff_atoms_.size() == loop_landing_atoms_.size() );
	runtime_assert( loop_takeoff_atoms_.size() == rna_loop_lengths_.size() );

}

///////////////////////////////////////////////////////////////////////////////
void
RNA_LoopEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {


	Real intermol_( 0.0 );

	totals[ rna_loop_fixed    ] = 0.0;
	totals[ rna_loop_logN     ] = 0.0;
	totals[ rna_loop_harmonic ] = 0.0;

	for ( Size n = 1; n <= loop_takeoff_atoms_.size(); n++ ){

		/////////////////////////////////////////////////////////////////////
		// score function weighting will determine actual penalty.
		totals[ rna_loop_fixed ] += 1.0;

		// prefactor should be something like (3/2) * k_B T.
		// score function weighting will determine actual penalty.
		totals[ rna_loop_logN  ] += log( static_cast<Real>( rna_loop_lengths_[ n ] ) );

		/////////////////////////////////////////////////////////////////////
		// harmonic term -- assuming a random Gaussian chain.
		Real const gaussian_variance = rna_loop_lengths_[ n ] * persistence_length2_;
		Real const loop_distance2 = get_loop_distance2( pose.xyz( loop_takeoff_atoms_[ n ] ), pose.xyz( loop_landing_atoms_[ n ] ) );
		totals[ rna_loop_harmonic ] += loop_distance2 / ( 2 * gaussian_variance );

	}

} // finalize_total_energy

///////////////////////////////////////////////////////////////////////////////
Real
RNA_LoopEnergy::get_loop_distance2( Vector const & v_takeoff, Vector const & v_landing ) const {

	Vector dist_vec = v_takeoff - v_landing;

	return dist_vec.length_squared();

}


///////////////////////////////////////////////////////////////////////////////
Vector
RNA_LoopEnergy::get_loop_distance2_deriv( Vector const & v_takeoff, Vector const & v_landing ) const {

	Vector dist_vec = 2 * (v_takeoff - v_landing);

	return dist_vec;

}

////////////////////////////////////////////////////////////////////////////////////
// silly helper function. Perhaps this should be made available as part of vector1.
Size
find_index( utility::vector1< id::AtomID > const & vec, id::AtomID const & value ){

	for ( Size n = 1; n <= vec.size(); n++ ) if ( vec[ n ] == value ) return n;

	utility_exit_with_message( "find_index failed" );

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
void
RNA_LoopEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
 	) const
{

	using namespace conformation;

	runtime_assert( loop_takeoff_atoms_.size() == loop_landing_atoms_.size() );
	runtime_assert( loop_takeoff_atoms_.size() == rna_loop_lengths_.size() );

	// following assumes that setup_for_scoring or setup_for_derivatives has been calculated, which
	// fills loop_takeoff_atoms_, loop_landing_atoms_, and gaussian_variances_

	if ( loop_takeoff_atoms_.has_value( atom_id ) ){

		//		tr.Debug << "Found loop takeoff atom! " << atom_id << std::endl;

		runtime_assert( ! loop_landing_atoms_.has_value( atom_id ) );

		Size const n = find_index( loop_takeoff_atoms_, atom_id );
		id::AtomID const & loop_landing_atom = loop_landing_atoms_[ n ];
		Real const gaussian_variance = rna_loop_lengths_[ n ] * persistence_length2_;

		Vector f2 = get_loop_distance2_deriv( pose.xyz( atom_id ), pose.xyz( loop_landing_atom ) ) / ( 2 * gaussian_variance );
		Vector f1 = cross( f2, pose.xyz( atom_id ) );

		F1 += weights[ rna_loop_harmonic ] * f1;
		F2 += weights[ rna_loop_harmonic ] * f2;

	}

	if ( loop_landing_atoms_.has_value( atom_id ) ){

		//		tr.Debug << "Found loop landing atom! " << atom_id << std::endl;

		runtime_assert( ! loop_takeoff_atoms_.has_value( atom_id ) );

		Size const n = find_index( loop_landing_atoms_, atom_id );
		id::AtomID const & loop_takeoff_atom = loop_takeoff_atoms_[ n ];
		Real const gaussian_variance = rna_loop_lengths_[ n ] * persistence_length2_;

		//		tr.Debug << "Found corresponding loop takeoff atom! " << loop_takeoff_atom << std::endl;

		Vector f2 = get_loop_distance2_deriv( pose.xyz( atom_id ), pose.xyz( loop_takeoff_atom ) ) / ( 2 * gaussian_variance );
		Vector f1 = cross( f2, pose.xyz( atom_id ) );

		F1 += weights[ rna_loop_harmonic ] * f1;
		F2 += weights[ rna_loop_harmonic ] * f2;

	}


} // eval atom derivative

core::Size
RNA_LoopEnergy::version() const
{
	return 1; // Initial versioning
}



} // rna
} // scoring
} // core
