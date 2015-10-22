// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/MgWaterHydrogenPacker.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/magnesium/MgWaterHydrogenPacker.hh>
#include <protocols/magnesium/util.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/UniformRotationSampler.hh>
#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.magnesium.MgWaterHydrogenPacker" );

using namespace core;
using utility::vector1;
typedef  numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace magnesium {


//Constructor
MgWaterHydrogenPacker::MgWaterHydrogenPacker():
	excise_mini_pose_( true ),
	use_fast_heuristic_( true )
{}

//Constructor
MgWaterHydrogenPacker::MgWaterHydrogenPacker( utility::vector1< Size > const & mg_res_list ):
	mg_res_list_( mg_res_list ),
	excise_mini_pose_( true ),
	use_fast_heuristic_( true )
{}

//Destructor
MgWaterHydrogenPacker::~MgWaterHydrogenPacker()
{}

///////////////////////////////////////////
// go over each Mg(2+)
///////////////////////////////////////////
void
MgWaterHydrogenPacker::apply( pose::Pose & pose )
{
	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::chemical::rna;

	// go over all Mg(2+)-bound waters:
	vector1< Size > mg_res_list( mg_res_list_ );
	if ( mg_res_list.size() == 0 ) mg_res_list = get_mg_res( pose );
	mg_water_pairs_ = get_mg_water_pairs( pose, mg_res_list );
	for ( Size i = 1; i <= mg_water_pairs_.size(); i++ ) {

		clock_t start_time = clock();
		apply( pose, mg_water_pairs_[ i ] );
		clock_t end_time = clock();

		std::cout << "pack_mg_water finished in " << double(end_time - start_time) / CLOCKS_PER_SEC << " seconds." << std::endl;
	}

	std::cout << "Processed " << mg_water_pairs_.size() <<  std::endl;
}

///////////////////////////////////////////
void
MgWaterHydrogenPacker::apply( pose::Pose & pose,
	std::pair< Size, Size > const & mg_water ) {
	using namespace core::id;
	using namespace core::pose;

	///////////////////////////////////////
	///////////////////////////////////////
	///////////////////////////////////////
	///////////////////////////////////////
	// UNIFY THIS WITH THE OTHER EXCISE_MINI_POSE
	//  IN HYDRATER.
	///////////////////////////////////////
	///////////////////////////////////////
	///////////////////////////////////////
	///////////////////////////////////////
	if ( excise_mini_pose_ ) {

		Pose pose_full = pose;

		Size const mg_res_in_full    = mg_water.first;
		Size const water_res_in_full = mg_water.second;
		vector1< Size > slice_res = pdbslice( pose, water_res_in_full );
		runtime_assert( slice_res.has_value( mg_res_in_full ) );

		Size const mg_res_in_mini    = slice_res.index( mg_res_in_full );
		Size const water_res_in_mini = slice_res.index( water_res_in_full );
		pack_mg_water_hydrogens_in_pose( pose, std::make_pair( mg_res_in_mini, water_res_in_mini ) );

		pose_full.set_xyz( NamedAtomID( " H1 ", water_res_in_full ),
			pose.xyz( NamedAtomID( " H1 ", water_res_in_mini ) ) );
		pose_full.set_xyz( NamedAtomID( " H2 ", water_res_in_full ),
			pose.xyz( NamedAtomID( " H2 ", water_res_in_mini ) ) );

		pose = pose_full;
	} else {
		pack_mg_water_hydrogens_in_pose( pose, mg_water );
	}
}

///////////////////////////////////////////////////////////
void
MgWaterHydrogenPacker::pack_mg_water_hydrogens_in_pose( pose::Pose & pose,
	std::pair< Size, Size > const & mg_water_pair ) {
	///////////////////////////////////////////
	// Some of the following was copied from dimaio/waterstuff branch:
	//
	//  src/core/pack/rotamer_set/rotamer_building_functions.cc
	//
	//  Function build_rotated_water_rotamers.
	//
	///////////////////////////////////////////
	using namespace core::conformation;
	using namespace core::scoring;


	Size const mg_res = mg_water_pair.first;
	Size const water_res = mg_water_pair.second;
	TR.Debug << "Doing water ligand with PDB numbering " << pose.pdb_info()->number( water_res ) << " near Mg " << pose.pdb_info()->number( mg_res ) << std::endl;

	core::Vector const MG( pose.residue( mg_res ).xyz("MG") );
	core::Vector const O ( pose.residue(water_res).xyz("O") );

	// uniformly sample rotations
	if ( urs_ == 0 ) urs_ = get_water_uniform_rotation_sampler();

	chemical::ResidueTypeSet const & residue_set( pose.residue( water_res ).residue_type_set() );
	Vector Oc, OH1c, OH2c;
	ResidueOP rsd_canonic = get_useful_HOH_coords( Oc, OH1c, OH2c, residue_set );

	// full enumeration!
	numeric::xyzMatrix< core::Real > best_R( numeric::xyzMatrix<Real>::identity() );
	numeric::xyzMatrix< core::Real > R( best_R );

	vector1< Vector > acc_vecs, rep_vecs;
	Vector mg_vec( 0.0 );
	ScoreFunctionOP scorefxn = get_mg_scorefxn();
	Real best_score( 0.0 );
	if ( use_fast_heuristic_ ) {
		find_water_neighbor_vecs( pose, water_res, mg_res, acc_vecs, rep_vecs, mg_vec );
		best_score = get_heuristic_water_hydrogen_score( acc_vecs, rep_vecs, mg_vec, OH1c, OH2c, R );
	} else {
		rotate_water_away_from_magnesium( pose, water_res, O, OH1c, OH2c, MG, R );
		best_score = ( *scorefxn )( pose );
	}

	for ( core::Size ctr=1; ctr<=urs_->nrots() ; ++ctr ) {
		urs_->get(ctr, R);
		Real score( 0.0 );
		if ( use_fast_heuristic_ ) {
			score = get_heuristic_water_hydrogen_score( acc_vecs, rep_vecs, mg_vec, OH1c, OH2c, R );
		} else {
			bool pointed_away_from_mg = rotate_water_away_from_magnesium( pose, water_res, O, OH1c, OH2c, MG, R );
			if ( !pointed_away_from_mg ) continue;
			score = ( *scorefxn )( pose );
		}
		//  std::cout << "score: " << F(8,3,score) << " for rotation " << I(4,ctr) << " out of " << I(6,urs_->nrots()) << std::endl;
		if ( score < best_score ) {
			best_R = R;
			best_score = score;
		}
	}
	TR.Debug << "Number of acceptors " << acc_vecs.size() << "   number of rep " << rep_vecs.size() << " -- best_score: " << best_score << std::endl;
	rotate_water_away_from_magnesium( pose, water_res, O, OH1c, OH2c, MG, best_R );

	return;
}

////////////////////////////////////////////////////////////////
// point towards acceptors, away from 'rep' atoms (Mg2+, etc.)
////////////////////////////////////////////////////////////////
Real
MgWaterHydrogenPacker::get_heuristic_water_hydrogen_score( utility::vector1< Vector > const & acc_vecs,
	utility::vector1< Vector > const & rep_vecs,
	Vector const & mg_vec,
	Vector const & OH1c,
	Vector const & OH2c,
	Matrix const & R ) const {
	vector1< Vector > OH_vecs = utility::tools::make_vector1( ( R * OH1c ).normalized(),
		( R * OH2c ).normalized() );
	Real score( 0.0 );

	// try to point towards acceptors -- but allow acceptor to "H-bond" to only one OH
	for ( Size j = 1; j <= acc_vecs.size(); j++ ) {
		Real acc_score( 0.0 );
		for ( Size i = 1; i <= OH_vecs.size(); i++ ) {
			Real possible_acc_score = -1.0 * dot( OH_vecs[ i ], acc_vecs[ j ] );
			if ( possible_acc_score < acc_score ) acc_score = possible_acc_score;
		}
		score += acc_score;
	}

	// weak preference to point away from non-acceptors
	Real const rep_weight( 0.25 );
	for ( Size j = 1; j <= rep_vecs.size(); j++ ) {
		for ( Size i = 1; i <= OH_vecs.size(); i++ ) {
			score += (  rep_weight ) * dot( OH_vecs[ i ], rep_vecs[ j ] );
		}
	}

	// strong preference to point away from Mg2+;  any angle better than 120.0 is fine.
	if ( mg_vec.length() > 0.0 ) {
		for ( Size i = 1; i <= OH_vecs.size(); i++ ) {
			score += (  1.0 ) * std::max( dot( OH_vecs[ i ], mg_vec ), -0.5 );
		}
	}

	return score;
}

///////////////////////////////////////////////////////////
// used for fast (heuristic) hydrogen packing.
void
MgWaterHydrogenPacker::find_water_neighbor_vecs( pose::Pose const & pose,
	Size const water_res,
	Size const mg_res,
	utility::vector1< Vector > & acc_vecs,
	utility::vector1< Vector > & rep_vecs,
	Vector & mg_vec ) const {
	using namespace core::conformation;
	using namespace core::chemical;
	Distance NBR_CUTOFF( 3.2 );
	Vector xyz_water = pose.residue( water_res ).xyz( "O" );
	Vector xyz_mg = pose.residue( mg_res ).xyz( "MG" );
	acc_vecs.clear();
	rep_vecs.clear();
	for ( Size j = 1; j <= pose.total_residue(); j++ ) {
		Residue const & rsd_j = pose.residue( j );
		if ( j == water_res ) continue;
		for ( Size jj = 1; jj <= rsd_j.natoms(); jj++ ) {
			Vector const & xyz_nbr = rsd_j.xyz( jj );
			if ( rsd_j.is_virtual( jj ) ) continue;
			// don't include other waters around mg2+
			if ( j != mg_res && ( xyz_nbr - xyz_mg ).length() < NBR_CUTOFF ) continue;

			AtomType atom_type = pose.residue( j ).atom_type( jj );
			if ( ( xyz_nbr - xyz_water ).length() < NBR_CUTOFF ) {
				Vector const vec = ( xyz_nbr - xyz_water ).normalized();
				if ( j == mg_res && atom_type.name() == "Mg2p" ) {
					mg_vec = vec;
					continue;
				}
				if ( atom_type.is_acceptor() ) {
					TR.Debug << "Found acceptor  neighbor for " << pose.pdb_info()->number( water_res ) << ": " << pose.pdb_info()->number( j ) << " " << rsd_j.atom_name( jj ) << std::endl;
					acc_vecs.push_back( vec );
				} else {
					TR.Debug << "Found repulsive neighbor for " << pose.pdb_info()->number( water_res ) << ": " << pose.pdb_info()->number( j ) << " " << rsd_j.atom_name( jj ) << std::endl;
					rep_vecs.push_back( vec );
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// checks if the two OH's are pointed away from Mg(2+).
// If so, makes the rotation. Otherwise, no change.
///////////////////////////////////////////////////////////////////////////////
bool
MgWaterHydrogenPacker::rotate_water_away_from_magnesium( pose::Pose & pose,
	Size const seqpos,
	Vector const & O,
	Vector const & OH1c,
	Vector const & OH2c,
	Vector const & MG,
	numeric::xyzMatrix< core::Real > const & R ) const {

	using namespace core::id;

	// rotated OH vectors, relative to O
	core::Vector newOH1 = R * OH1c;
	core::Vector newOH2 = R * OH2c;

	// Mg(2+) relative to O
	core::Vector mg_o = ( MG - O );

	Real align_cutoff( 0.0 ); // could make -0.3 for more stringent directionality...
	if ( dot( mg_o, newOH1 ) / ( mg_o.length() * newOH1.length() ) > align_cutoff ) return false;
	if ( dot( mg_o, newOH2 ) / ( mg_o.length() * newOH2.length() ) > align_cutoff ) return false;

	// add <x,y,z> to 'O'; compute H's from this
	pose.set_xyz( NamedAtomID( " H1 ", seqpos), O+newOH1);
	pose.set_xyz( NamedAtomID( " H2 ", seqpos), O+newOH2);
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
MgWaterHydrogenPacker::remove_waters_except_mg_bound( pose::Pose & pose ) const {
	vector1< std::pair< Size, Size > > mg_water_pairs = get_mg_water_pairs( pose );
	protocols::magnesium::remove_waters_except_mg_bound( pose, mg_water_pairs /*except these*/ );
}


} //magnesium
} //protocols
