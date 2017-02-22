// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNA_FullAtomStacking.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_StubCoordinateEnergy.hh>
#include <core/scoring/rna/RNA_StubCoordinateEnergyCreator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/DenseEnergyContainer.hh>
#include <basic/Tracer.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/kinematics/Stub.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

// Utility headers
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

//Auto Headers
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Computes penalty for having an O5' atom far from a target location, as defined by the
///                  O3' stub of a previous residue.
///
///  Using this for umbrella sampling of end-to-end geometries of RNA loops with rna_denovo--
///   the resulting six-dimensional potential will go into an improved loop_close term for stepwise
///   and other RNA modeling applications
///
///  TODO: Ideally, would make a clean version that lives in StubCoordinateConstraint and
///           handles derivatives. Probably could sub-class from CoordinateConstraint.
///
///    -- rhiju, 2017
///////////////////////////////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "core.scoring.rna.RNA_StubCoordinateEnergy" );

using namespace basic::options;
using namespace basic::options::OptionKeys::rna;

namespace core {
namespace scoring {
namespace rna {

/// @details This must return a fresh instance of the RNA_StubCoordinateEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_StubCoordinateEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RNA_StubCoordinateEnergy );
}

ScoreTypes
RNA_StubCoordinateEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_stub_coord_hack );
	return sts;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// c-tor
RNA_StubCoordinateEnergy::RNA_StubCoordinateEnergy() :
	parent( methods::EnergyMethodCreatorOP( new RNA_StubCoordinateEnergyCreator ) )
{
	// meant for use with RNA_FragmentMonteCarlo in mode where it simuilates an RNA loop and
	//  outputs translation & rotation information from a takeoff to a landing residue.
	runtime_assert( option[ farna::out::output_jump_o3p_to_o5p ]() );
	runtime_assert( option[ farna::out::output_jump_res ]().size() == 2 );
	takeoff_res_ = option[ farna::out::output_jump_res ]()[ 1 ];
	landing_res_ = option[ farna::out::output_jump_res ]()[ 2 ];

	runtime_assert( option[ farna::out::target_xyz ]().size() == 3 );
	utility::vector1< Real > const & xyz = option[ farna::out::target_xyz ]();
	target_xyz_in_takeoff_frame_ = Vector( xyz[1], xyz[2], xyz[3] );

	// Following defines: ( target_xyz - xyz )**2   [note: no 1/2]
	using namespace core::scoring::func;
	func_ = FuncOP( new HarmonicFunc( 0.0 /*target distance*/, 1.0 /*standard deviation*/ ) );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
methods::EnergyMethodOP
RNA_StubCoordinateEnergy::clone() const
{
	return methods::EnergyMethodOP( new RNA_StubCoordinateEnergy );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// copied from SASAEnergy.cc
void
RNA_StubCoordinateEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	using namespace core::scoring::methods;
	LongRangeEnergyType const & lr_type( long_range_type() );

	// create a container
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast< DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.size() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		LREnergyContainerOP new_dec = LREnergyContainerOP( new DenseEnergyContainer( pose.size(), rna_stub_coord_hack ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_StubCoordinateEnergy::defines_residue_pair_energy(
														pose::Pose const &,
														Size res1,
														Size res2
														) const
{
	if ( res1 == takeoff_res_ && res2 == landing_res_ ) return true;
	if ( res1 == landing_res_ && res2 == takeoff_res_ ) return true;
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StubCoordinateEnergy::residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & scorefxn,
		EnergyMap & emap
) const
{
	using namespace core::kinematics;
	using namespace core::chemical::rna;
	using namespace core::conformation;

	if ( rsd1.seqpos() == landing_res_ &&
			 rsd2.seqpos() == takeoff_res_ ) {
		residue_pair_energy( rsd2, rsd1, pose, scorefxn, emap );
		return;
	}

	if ( rsd1.seqpos() != takeoff_res_ ) return;
	if ( rsd2.seqpos() != landing_res_ ) return;

	// Code copied from protocols::farna::RNA_FragmentMonteCarlo.cc:
	// takeoff
	Stub stub1( rsd1.xyz( " O3'") /* center */,
							rsd1.xyz( " O3'") /* a */,
							rsd1.xyz( " C3'") /* b  [b->a defines x] */,
							rsd1.xyz( " C4'") /* c  [c->b defines y] */ );
	stub1.M = Matrix::cols( stub1.M.col_y(), stub1.M.col_z(), stub1.M.col_x() ); // Prefer to have C3'->O3' (takeoff vector) along z

	// could easily generalize this function to also compute penalties for having rotation different from target rotation
	//   using, e.g., a von-Mises function.
	// landing
	// stub2 = Stub( rsd2.xyz( " O5'") /* center */,
	// 							rsd2.xyz( " C5'") /* a */,
	// 							rsd2.xyz( " O5'") /* b  [b->a defines x] */,
	// 							rsd2.xyz( " C4'") /* c  [c->b defines y] */ );
	// stub2.M = Matrix::cols( stub2.M.col_y(), stub2.M.col_z(), stub2.M.col_x() ); // Prefer to have O5'->C5' (landing vector) along z

	Vector const & landing_xyz = rsd2.xyz( " O5'" );
	Vector const landing_xyz_in_takeoff_frame = stub1.global2local( landing_xyz );
	emap[ rna_stub_coord_hack ] += func_->func( landing_xyz_in_takeoff_frame.distance( target_xyz_in_takeoff_frame_ ) );

}


} //rna
} //scoring
} //core
