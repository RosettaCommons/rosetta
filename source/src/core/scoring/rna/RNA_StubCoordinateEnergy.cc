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
#include <core/pose/rna/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>
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
///  TODO: Could also use this framework to align rotations of coordinate frames. See notes below.
///
///  TODO: Setup still copies a ton of code from RNA_FragmentMonteCarlo. Could be
///           unified if we use Constraint above.
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
	using namespace core::pose::rna;
	// meant for use with RNA_FragmentMonteCarlo in mode where it simuilates an RNA loop and
	//  outputs translation & rotation information from a takeoff to a landing residue.
	stub_stub_type_ = BASE_CENTROID;
	if ( option[ farna::out::output_jump_o3p_to_o5p ]() ) stub_stub_type_ = O3P_TO_O5P;
	if ( option[ farna::out::output_jump_chainbreak ]() ) stub_stub_type_ = CHAINBREAK;
	if ( option[ farna::out::output_jump_reference_RT ].user() ){
		reference_RT_ = kinematics::RTOP( new kinematics::RT );
		std::stringstream rt_stream( option[ farna::out::output_jump_reference_RT ]() );
		rt_stream >> *reference_RT_;
	}

	jump_resnum_and_chain_ = option[ farna::out::output_jump_res ].resnum_and_chain();
	runtime_assert( jump_resnum_and_chain_.first.size() == 2 );
	runtime_assert( jump_resnum_and_chain_.second.size() == 2 );

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

	using namespace core::pose::full_model_info;
	if ( full_model_info_defined( pose ) ) {
		res1_ = const_full_model_info( pose ).full_model_parameters()->conventional_to_full( jump_resnum_and_chain_.first[1],
																																												 jump_resnum_and_chain_.second[1] );
		res2_ = const_full_model_info( pose ).full_model_parameters()->conventional_to_full( jump_resnum_and_chain_.first[2],
																																												 jump_resnum_and_chain_.second[2] );
	} else {
		res1_ = jump_resnum_and_chain_.first[ 1 ];
		res2_ = jump_resnum_and_chain_.first[ 2 ];
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
	if ( res1 == res1_ && res2 == res2_ ) return true;
	if ( res1 == res2_ && res2 == res1_ ) return true;
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

	if ( rsd1.seqpos() == res2_ &&
			 rsd2.seqpos() == res1_ ) {
		residue_pair_energy( rsd2, rsd1, pose, scorefxn, emap );
		return;
	}

	if ( rsd1.seqpos() != res1_ ) return;
	if ( rsd2.seqpos() != res2_ ) return;

	Stub stub1, stub2;
	core::pose::rna::get_stub_stub( rsd1, rsd2, stub1, stub2, stub_stub_type_ );

	if ( reference_RT_ != 0 ) reference_RT_->make_jump( stub1 /*start*/, stub1 /*end*/ );

	// could easily generalize this function to also compute penalties for having rotation different from target rotation
	//   using, e.g., a von-Mises function -- compute Jump( stub1, stub2 ) and then check out, e.g. axis-angle of rotation.
	Vector const & landing_xyz = stub2.center();
	Vector const landing_xyz_in_takeoff_frame = stub1.global2local( landing_xyz );
	emap[ rna_stub_coord_hack ] += func_->func( landing_xyz_in_takeoff_frame.distance( target_xyz_in_takeoff_frame_ ) );

}


} //rna
} //scoring
} //core
