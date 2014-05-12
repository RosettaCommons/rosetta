// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/phosphate/PhosphateMover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/rna/phosphate/PhosphateMover.hh>
#include <protocols/stepwise/sampling/rna/phosphate/util.hh>
#include <core/chemical/rna/RNA_SamplerUtil.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.rna.phosphate.PhosphateMover" );

using namespace core;
typedef  utility::vector1< Real >  TorsionList;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace phosphate {

	//Constructor
	PhosphateMover::PhosphateMover( Size const sample_res,
																	PhosphateTerminus const which_terminus,
																	scoring::ScoreFunctionCOP scorefxn ):
		phosphate_move_( PhosphateMove( sample_res, which_terminus ) ),
		scorefxn_( scorefxn )
	{
		initialize_variables();
	}

	//Constructor
	PhosphateMover::PhosphateMover( PhosphateMove const phosphate_move,
																	scoring::ScoreFunctionCOP scorefxn ):
		phosphate_move_( phosphate_move ),
		scorefxn_( scorefxn )
	{
		initialize_variables();
	}

	/////////////////////
	void
	PhosphateMover::initialize_variables(){
		do_screening_ = true;
		screen_for_donor_contact_ = true;
		instantiated_phosphate_ = false;
		force_phosphate_instantiation_ = false; // do not check in as true, just for testing
		number_score_calls_ = 0;
	}

	//Destructor
	PhosphateMover::~PhosphateMover()
	{}

	/////////////////////
	std::string
	PhosphateMover::get_name() const {
		return "PhosphateMover";
	}


	/////////////////////
	void
	PhosphateMover::apply( core::pose::Pose & pose ) {
		setup_variants_and_free_pose_for_terminal_phosphate( pose );
		if ( do_screening_ )	screen_phosphate( pose );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	PhosphateMover::screen_phosphate( pose::Pose & pose ) {
		instantiated_phosphate_ = false;
		number_score_calls_ = 0;
		setup_atom_and_neighbor_list( pose );
		if ( phosphate_move_.terminus() == FIVE_PRIME_PHOSPHATE ){
			screen_five_prime_phosphate(  pose );
		} else {
			runtime_assert( phosphate_move_.terminus() == THREE_PRIME_PHOSPHATE );
			screen_three_prime_phosphate( pose );
		}
		//TR << "After " << number_score_calls_ << " did we instantiate phosphate?: " << instantiated_phosphate_ << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	PhosphateMover::setup_variants_and_free_pose_for_terminal_phosphate( pose::Pose & pose  ){

		if ( phosphate_move_.terminus() == FIVE_PRIME_PHOSPHATE ){
			setup_variants_and_free_pose_for_five_prime_phosphate( pose );
		} else {
			runtime_assert( phosphate_move_.terminus() == THREE_PRIME_PHOSPHATE );
			setup_variants_and_free_pose_for_three_prime_phosphate( pose );
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	PhosphateMover::setup_variants_and_free_pose_for_five_prime_phosphate( pose::Pose & pose ){

		Size const sample_res = phosphate_move_.rsd();

		if ( pose.residue( sample_res ).has_variant_type( "FIVE_PRIME_PHOSPHATE" ) ){

			runtime_assert( !pose.residue( sample_res ).has_variant_type( "VIRTUAL_PHOSPHATE" ) );
			if ( pose.residue( sample_res ).has_variant_type( "VIRTUAL_RIBOSE" ) ){
				std::cerr << pose.fold_tree() << std::endl;
				std::cerr << pose.annotated_sequence() << std::endl;
				std::cerr << phosphate_move_ << std::endl;
			}

			runtime_assert( !pose.residue( sample_res ).has_variant_type( "VIRTUAL_RIBOSE" ) );

			pose_free_ = pose.clone();
			remove_variant_type_from_pose_residue( *pose_free_, "FIVE_PRIME_PHOSPHATE", sample_res );
			add_variant_type_to_pose_residue( *pose_free_, "VIRTUAL_PHOSPHATE", sample_res );

		} else {
			runtime_assert( !pose.residue( sample_res ).has_variant_type( "FIVE_PRIME_PHOSPHATE" ) );
			runtime_assert( pose.residue( sample_res ).has_variant_type( "VIRTUAL_PHOSPHATE" ) );

			pose_free_ = pose.clone();

			remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", sample_res );
			remove_variant_type_from_pose_residue( pose, "LOWER_TERMINUS",    sample_res );
			remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE",    sample_res );
			add_variant_type_to_pose_residue( pose, "FIVE_PRIME_PHOSPHATE", sample_res );
			correctly_position_five_prime_phosphate( pose, sample_res );

		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	PhosphateMover::setup_variants_and_free_pose_for_three_prime_phosphate( pose::Pose & pose ){

		Size const sample_res = phosphate_move_.rsd();

		if ( pose.residue( phosphate_move_.rsd() ).has_variant_type( "THREE_PRIME_PHOSPHATE" ) ){
			pose_free_ = pose.clone();
			remove_variant_type_from_pose_residue( *pose_free_, "THREE_PRIME_PHOSPHATE", sample_res );
		} else {
			pose_free_ = pose.clone();
			remove_variant_type_from_pose_residue( pose, "UPPER_TERMINUS",    sample_res );
			remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", sample_res );
			add_variant_type_to_pose_residue( pose, "THREE_PRIME_PHOSPHATE", sample_res );
		}

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	PhosphateMover::setup_atom_and_neighbor_list( pose::Pose & pose ) {

		static Real const  phosphate_takeoff_donor_distance_cutoff2_( 8.0 * 8.0 );
		static Real const  phosphate_nbr_distance_cutoff2_( 12.0 * 12.0 );
		Size const & n = phosphate_move_.rsd();
		Vector const phosphate_takeoff_xyz = (phosphate_move_.terminus() == FIVE_PRIME_PHOSPHATE) ?
			pose.residue( n ).xyz( " C5'" ) :
			pose.residue( n ).xyz( " O3'" );

		donor_atom_xyz_list_.clear();
		donor_base_atom_xyz_list_.clear();
		neighbor_copy_dofs_.clear();

		for ( Size m = 1; m <= pose.total_residue(); m++ ){

			if ( ( pose.residue( m ).nbr_atom_xyz()  - phosphate_takeoff_xyz ).length_squared() > phosphate_nbr_distance_cutoff2_ ) continue;
			neighbor_copy_dofs_.push_back( m );

			core::chemical::AtomIndices Hpos_polar = pose.residue_type( m ).Hpos_polar();
			for ( Size ii = 1; ii <= Hpos_polar.size(); ii++ ){
				Size const & q = Hpos_polar[ ii ];
				Vector const & donor_atom_xyz = pose.residue( m ).xyz( q );
				if ( ( donor_atom_xyz - phosphate_takeoff_xyz ).length_squared() < phosphate_takeoff_donor_distance_cutoff2_ ){
					donor_atom_xyz_list_.push_back( donor_atom_xyz );
					Size  q_base = pose.residue_type(m).atom_base( q );
					Vector const & donor_base_atom_xyz = pose.residue( m ).xyz( q_base );
					donor_base_atom_xyz_list_.push_back( donor_base_atom_xyz );
				}
			}

		}

		Size const & terminus = phosphate_move_.terminus();
		if ( terminus == FIVE_PRIME_PHOSPHATE ){
			op1_atom_idx_ = pose.residue( n ).atom_index( " OP1" );
			op2_atom_idx_ = pose.residue( n ).atom_index( " OP2" );
		} else {
			op1_atom_idx_ = pose.residue( n ).atom_index( "YOP1" );
			op2_atom_idx_ = pose.residue( n ).atom_index( "YOP2" );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	PhosphateMover::screen_five_prime_phosphate( pose::Pose & pose ) {

		using namespace core::chemical::rna;
		using namespace core::id;
		using namespace core::scoring;

		Real const bin_size_( 20.0 );
		TorsionList full_torsions = get_full_torsions( bin_size_ );

		Real const score_free = ( *scorefxn_ )( *pose_free_ );
		Size const & n = phosphate_move_.rsd();

		Pose pose_best = pose;
		Real score_best = ( *scorefxn_ )( pose );
		for ( Size i = 1; i <= full_torsions.size(); i++ ){
			pose.set_torsion( TorsionID( n, BB, GAMMA ), full_torsions[i] );
			if ( !pass_clash_check( " O5'", n, pose ) ) continue;

			for ( Size j = 1; j <= full_torsions.size(); j++ ){
				pose.set_torsion( TorsionID( n, BB, BETA ), full_torsions[j] );
				if ( !pass_clash_check( " P  ", n, pose ) ) continue;

				for ( Size k = 1; k <= full_torsions.size(); k++ ){
					pose.set_torsion( TorsionID( n, BB, ALPHA ), full_torsions[k] );
					//					if ( !pass_clash_check( " OP1", n, pose ) ) continue;
					//					if ( !pass_clash_check( " OP2", n, pose ) ) continue;
					//					if ( !pass_clash_check( "XO3'", n, pose ) ) continue;
					if ( screen_for_donor_contact_ && !check_phosphate_contacts_donor( pose ) ) continue;
					Real const score = ( *scorefxn_ )( pose );
					//					TR << "Score: " << score << " compared to best score " << score_best << std::endl;
					//					scorefxn_->show( pose );
					number_score_calls_++;
					if ( score <= score_best ) {
						score_best = score;
						pose_best = pose;
					}
				}
			}
		}
		if ( score_best < score_free || force_phosphate_instantiation_ ){
			pose = pose_best;
			instantiated_phosphate_ = true;
		} else {
			pose = *pose_free_;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	PhosphateMover::screen_three_prime_phosphate( pose::Pose & pose ){

		using namespace core::chemical::rna;
		using namespace core::id;
		using namespace core::scoring;

		Size const & n = phosphate_move_.rsd();

		Real const bin_size_( 20.0 );
		TorsionList full_torsions = get_full_torsions( bin_size_ );
		bool const extra_epsilon_( true );
		TorsionList epsilon_torsions = get_epsilon_torsions( pose.delta( n ), extra_epsilon_, bin_size_ );

		Real const score_free = ( *scorefxn_ )( *pose_free_ );
		//		TR << "score for start without 3' phosphate" << score_free << std::endl;

		Pose pose_best = pose;
		Real score_best = ( *scorefxn_ )( pose );
		for ( Size i = 1; i <= epsilon_torsions.size(); i++ ){
			pose.set_torsion( TorsionID( n, BB, EPSILON ), epsilon_torsions[i] );
			if ( !pass_clash_check( "YP  ", n, pose ) ) continue;

			for ( Size j = 1; j <= full_torsions.size(); j++ ){
				pose.set_torsion( TorsionID( n, BB, ZETA ), full_torsions[j] );
				if ( !pass_clash_check( "YOP1", n, pose ) ) continue;
				if ( !pass_clash_check( "YOP2", n, pose ) ) continue;
				if ( !pass_clash_check( "YO5'", n, pose ) ) continue;

				if ( screen_for_donor_contact_ && !check_phosphate_contacts_donor( pose ) ) continue;
				Real const score = ( *scorefxn_ )( pose );
				number_score_calls_++;
				if ( score <= score_best ) {
					score_best = score;
					pose_best = pose;
				}
			}
		}
		//TR << "score for sample with 3' phosphate " << score_best << " vs free " << score_free;

		if ( score_best < score_free || force_phosphate_instantiation_ ){
			pose = pose_best;
			instantiated_phosphate_ = true;
		} else {
			pose = *pose_free_;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	PhosphateMover::pass_clash_check( std::string const atom_name,
																		Size const n,
																		pose::Pose & pose ){

		Size const atom_idx = pose.residue_type( n ).atom_index( atom_name );
		Vector const & atom_xyz = pose.residue( n ).xyz( atom_idx );
		Real const VDW_radius_1 = pose.residue( n ).atom_type( atom_idx  ).lj_radius();
		static Real const clash_dist_cutoff = 0.8; //Fail van der Waals repulsion screen if two atoms radius within 0.5 Angstrom of each other

		for ( Size k = 1; k <= neighbor_copy_dofs_.size(); k++ ){
			Size const & m = neighbor_copy_dofs_[ k ];
			if ( m == n ) continue;
			for ( Size j = 1; j <= pose.residue( m ).natoms(); j++ ){
				if ( pose.residue( m ).is_virtual( j )  ) continue;
				Real const VDW_radius_2 = pose.residue( m ).atom_type( j ).lj_radius();
				Real const clash_radius = VDW_radius_1 + VDW_radius_2 - clash_dist_cutoff;
				if ( ( atom_xyz  - pose.residue( m ).xyz( j ) ).length_squared() < clash_radius*clash_radius ){
					return false;
				}
			}
		}
		return true;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	PhosphateMover::check_phosphate_contacts_donor( pose::Pose & pose ) const {
		static Real const max_hbond_dist_cutoff2_( 2.2 * 2.2 );
		static Real const max_hbond_base_dist_cutoff2_( 3.2 * 3.2 );
		static Real const cos_angle_cutoff = -0.5;

		Size const & n = phosphate_move_.rsd();
		utility::vector1< Vector > op_xyz_list;
		op_xyz_list.push_back( pose.residue( n ).xyz( op1_atom_idx_ ) );
		op_xyz_list.push_back( pose.residue( n ).xyz( op2_atom_idx_ ) );

		for ( Size i = 1; i <= donor_atom_xyz_list_.size(); i++ ){
			for ( Size j = 1; j <= op_xyz_list.size(); j++ ){
				Real const dist2 =  ( donor_atom_xyz_list_[ i ] - op_xyz_list[ j ] ).length_squared();
				if ( dist2 > max_hbond_dist_cutoff2_ ) continue;

				Real const dist_to_base2 =  ( donor_base_atom_xyz_list_[ i ] - op_xyz_list[ j ] ).length_squared();
				if ( dist_to_base2 > max_hbond_base_dist_cutoff2_  ) continue;

				Real const costheta = 	cos_of( donor_base_atom_xyz_list_[i], donor_atom_xyz_list_[i], op_xyz_list[j]  );
				if ( costheta > cos_angle_cutoff ) continue;
				//				TR << "dist! " << std::sqrt( dist2 ) << " to " << i << " out of " << donor_atom_xyz_list_.size() << " and cos angle is " <<	costheta <<	std::endl;
				return true;
			}
		}

		return false;
	}

} //phosphate
} //rna
} //sampling
} //stepwise
} //protocols
