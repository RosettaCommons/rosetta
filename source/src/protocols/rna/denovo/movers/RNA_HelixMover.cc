// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RNA_HelixMover.cc
/// @brief protocols that are specific to RNA_HelixMover
/// @details RNA helix rotations and translations along helical axis
/// @author Kalli Kappel


#include <protocols/rna/denovo/movers/RNA_HelixMover.hh>
#include <core/import_pose/RNA_BasePairHandler.hh>

#include <protocols/rna/denovo/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/conformation/Residue.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <basic/options/option.hh>


//Relaxer stuff

// ObjexxFCL Headers

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>

// External library headers


//C++ headers
#include <utility>
#include <vector>
#include <string>
#include <sstream>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


using namespace core;


static basic::Tracer TR( "protocols.rna.denovo.movers.RNA_HelixMover" );

namespace protocols {
namespace rna {
namespace denovo {
namespace movers {

//////////////////////////////////////////////////////////////////////////////////////////
RNA_HelixMover::RNA_HelixMover(  utility::vector1< utility::vector1< Size > > const & helix_regions,
	core::import_pose::RNA_BasePairHandlerCOP rna_base_pair_handler,
	bool const & move_first_rigid_body ):
	Mover(),
	helix_regions_( helix_regions ),
	rna_base_pair_handler_(std::move( rna_base_pair_handler )),
	pose_is_set_( false ),
	move_first_rigid_body_( move_first_rigid_body ),
	rot_mag_( 10.0 ),
	trans_mag_( 2.0 )
{
	Mover::type("RNA_HelixMover");
	get_helix_ends();
}

RNA_HelixMover::~RNA_HelixMover() = default;

/// @details  Apply the RNA helix mover
///
///////////////////////////////////////////////////////////////////////////////////////////
void RNA_HelixMover::apply( core::pose::Pose & pose )
{
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !pose_is_set_ ) {
		set_pose( pose );
	}

	if ( pose.fold_tree() != pose_fold_tree_ ) {
		set_pose( pose );
	}

	if ( helix_regions_with_jumps_and_ends_.size() == 0 ) {
		TR.Warning << "No helix regions that can be moved! Doing nothing." << std::endl;
		return;
	}

	// randomly choose a helical region
	core::Size region = static_cast<Size>( numeric::random::rg().uniform() * helix_regions_with_jumps_and_ends_.size() ) + 1;

	// get the helical axis
	std::pair< core::Vector, core::Vector > helix_axis_and_rot_center = get_helical_axis_and_center( pose, region );

	//choose rotate or translate
	bool rotate = false;
	if ( numeric::random::rg().uniform() < 0.5 ) {
		rotate = true;
	}

	if ( rotate ) {
		// set the helical axis and center of rotation
		spin_movers_[region]->spin_axis( helix_axis_and_rot_center.first );
		spin_movers_[region]->rot_center( helix_axis_and_rot_center.second );
		spin_movers_[region]->spin_mag( rot_mag_ );
		// apply the mover
		spin_movers_[ region ]->apply( pose );

	} else { //translate
		// set the helical axis
		trans_movers_[ region ]->trans_axis( helix_axis_and_rot_center.first );
		// pick a random step size
		//  core::Real step_size = numeric::random::rg().uniform()*trans_mag_;
		core::Real step_size = trans_mag_ * numeric::random::rg().gaussian();
		trans_movers_[ region ]->step_size( step_size );
		trans_movers_[ region ]->apply( pose );
	}

}

///////////////////////////////////////////////////////////////////////////////////////////
std::pair< core::Vector, core::Vector >
RNA_HelixMover::get_helical_axis_and_center( core::pose::Pose const & pose,
	Size const & region ) const
{

	core::Vector helix_axis( 0.0, 0.0, 0.0 );
	core::Vector helix_start, helix_end, rot_center;

	// first base pair
	helix_start = get_bp_center( pose, helix_ends_final_[ region ].first );

	// last base pair
	helix_end = get_bp_center( pose, helix_ends_final_[ region ].second );

	helix_axis = helix_start - helix_end;
	rot_center = ( helix_start + helix_end ) / 2.0;

	return std::make_pair( helix_axis, rot_center );

}
///////////////////////////////////////////////////////////////////////////////////////////
core::Vector
RNA_HelixMover::get_bp_center( core::pose::Pose const & pose,
	std::pair< core::Size, core::Size > const & bp_res ) const
{
	core::Vector res1 = get_bb_pos( pose, bp_res.first );
	core::Vector res2 = get_bb_pos( pose, bp_res.second );

	return ( res1 + res2 ) / 2.0;
}
///////////////////////////////////////////////////////////////////////////////////////////
core::Vector
RNA_HelixMover::get_bb_pos( core::pose::Pose const & pose,
	core::Size const & res ) const
{

	core::Vector bp_center( 0.0, 0.0, 0.0 );
	core::Vector backbone_centroid( 0.0, 0.0, 0.0 );
	bool backbone_centroid_calculated = false;

	// check that the residue has the atoms, if not average the backbone pos?
	if ( pose.residue( res ).has( " C1'" ) ) {
		bp_center += pose.residue( res ).xyz( " C1'" );
	} else {
		backbone_centroid = get_backbone_centroid( pose, res );
		bp_center += backbone_centroid;
		backbone_centroid_calculated = true;
	}

	if ( pose.residue( res ).is_purine() && pose.residue( res ).has( " N9 " ) ) {
		bp_center += pose.residue( res ).xyz( " N9 " );
	} else if ( pose.residue( res ).is_pyrimidine() && pose.residue( res ).has( " N1 " ) ) {
		bp_center += pose.residue( res ).xyz( " N1 " );
	} else {
		if ( !backbone_centroid_calculated ) {
			backbone_centroid = get_backbone_centroid( pose, res );
		}
		bp_center += backbone_centroid;
	}

	bp_center /= 2.0;

	return bp_center;

}
///////////////////////////////////////////////////////////////////////////////////////////
core::Vector
RNA_HelixMover::get_backbone_centroid( core::pose::Pose const & pose,
	core::Size const & res ) const
{
	core::Vector backbone_centroid( 0.0, 0.0, 0.0 );

	for ( Size j=1; j<=pose.residue( res ).last_backbone_atom(); ++j ) {
		backbone_centroid += pose.residue( res ).atom( j ).xyz();
	}
	backbone_centroid /= pose.residue( res ).last_backbone_atom();

	return backbone_centroid;
}
///////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixMover::set_pose( core::pose::Pose const & pose ) {

	pose_is_set_ = true;

	// save the fold tree - so we can check in apply that it hasn't changed
	pose_fold_tree_ = pose.fold_tree();

	// figure out the jump for each helix region
	utility::vector1< core::Size > rigid_body_jumps = core::pose::rna::get_rigid_body_jumps( pose );

	if ( !move_first_rigid_body_ ) {
		if ( rigid_body_jumps.size() >= 1 ) {
			rigid_body_jumps.erase( rigid_body_jumps.begin() );
		}
	}

	//helix_regions_jumps_ vector of the jumps
	std::map< Size, Size > jump_res_to_jump_num;

	// get the jump residues for each rigid body jump
	for ( Size i =1; i<=rigid_body_jumps.size(); ++i ) {
		// exactly one of these is guaranteed to be true
		// (b/c rigid body jumps need to be to the vrt res which is at the end of the pose
		// - get_rigid_body_jumps checks for this
		core::Size upstream_res = pose.fold_tree().upstream_jump_residue( rigid_body_jumps[i] );
		core::Size downstream_res = pose.fold_tree().downstream_jump_residue( rigid_body_jumps[i] );
		if ( upstream_res != pose.size() ) {
			jump_res_to_jump_num[ upstream_res ] = rigid_body_jumps[i];
		} else if ( downstream_res != pose.size() ) {
			jump_res_to_jump_num[ downstream_res ] = rigid_body_jumps[i];
		}
	}

	for ( Size i=1; i<=helix_regions_.size(); ++i ) {
		bool found_jump = false;
		for ( Size j=1; j<=helix_regions_[i].size(); ++j ) {
			// check whether it's a jump residue
			if ( jump_res_to_jump_num.count( helix_regions_[i][j] ) ) {
				helix_regions_jumps_.push_back( jump_res_to_jump_num[ helix_regions_[i][j] ] );
				found_jump = true;
				break; // only one jump per region
			}

		}
		if ( !found_jump ) {
			helix_regions_jumps_.push_back( 0 );
		}

	}

	for ( Size i=1; i<= helix_regions_.size(); ++i ) {

		if ( helix_regions_jumps_[i] == 0 ) continue;
		if ( helix_ends_[i].first.first == 0 ) continue;
		if ( helix_ends_[i].first.second == 0 ) continue;
		if ( helix_ends_[i].second.first == 0 ) continue;
		if ( helix_ends_[i].second.second == 0 ) continue;

		helix_regions_with_jumps_and_ends_.push_back( helix_regions_[i] );
		helix_ends_final_.push_back( helix_ends_[i] );
		helix_regions_jumps_final_.push_back( helix_regions_jumps_[i] );
	}

	// set up the movers
	for ( Size i=1; i<= helix_regions_with_jumps_and_ends_.size(); ++i ) {

		protocols::rigid::RigidBodySpinMoverOP spin_mover = protocols::rigid::RigidBodySpinMoverOP(
			new protocols::rigid::RigidBodySpinMover(  helix_regions_jumps_final_[i] /*jump*/) );

		spin_movers_.push_back( spin_mover );

		core::Vector zero_spin_axis(0.0, 0.0, 0.0);
		protocols::rigid::RigidBodyTransMoverOP trans_mover = protocols::rigid::RigidBodyTransMoverOP(
			new protocols::rigid::RigidBodyTransMover( zero_spin_axis, helix_regions_jumps_final_[i], false ));

		trans_movers_.push_back( trans_mover );

	}

	//std::cout << "HELIX REGIONS FINAL: " << helix_regions_with_jumps_and_ends_ << std::endl;
	//std::cout << "HELIX JUMPS FINAL: " << helix_regions_jumps_final_ << std::endl;
	//std::cout << "HELIX ENDS FINAL: " << helix_ends_final_ << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixMover::get_helix_ends() {

	std::map< core::Size, core::Size > base_pairs_first_second;
	std::map< core::Size, core::Size > base_pairs_second_first;
	for ( core::Size i=1; i<=rna_base_pair_handler_->rna_pairing_list().size(); ++i ) {
		//std::cout << "Res1: " << rna_base_pair_handler_->rna_pairing_list()[i].res1() <<
		// "Res2: " << rna_base_pair_handler_->rna_pairing_list()[i].res2() << std::endl;
		base_pairs_first_second[ rna_base_pair_handler_->rna_pairing_list()[i].res1() ] = rna_base_pair_handler_->rna_pairing_list()[i].res2();
		base_pairs_second_first[ rna_base_pair_handler_->rna_pairing_list()[i].res2() ] = rna_base_pair_handler_->rna_pairing_list()[i].res1();
	}

	for ( core::Size i = 1; i<= helix_regions_.size(); ++i ) {
		std::pair< std::pair< Size, Size>, std::pair< Size, Size > > helix_end;
		std::pair< Size, Size> end1 = std::make_pair( 0, 0 );
		std::pair< Size, Size> end2 = std::make_pair( 0, 0 );
		// don't bother if there's only a single base pair
		if ( helix_regions_[i].size() <= 2 ) {
			helix_end = std::make_pair( end1, end2 );
			helix_ends_.push_back( helix_end );
			continue;
		}
		std::sort( helix_regions_[i].begin(), helix_regions_[i].end() );
		//base_pair_ends_
		for ( core::Size j = 1; j<=helix_regions_[i].size(); ++j ) {
			// find the first residue with a partner
			if ( base_pairs_first_second.count( helix_regions_[i][j] ) && end1.first == 0 ) {
				// this is the first base pair
				end1 = std::make_pair( helix_regions_[i][j], base_pairs_first_second[ helix_regions_[i][j] ] );
			}
			if ( base_pairs_second_first.count( helix_regions_[i][j] ) && end2.first == 0 ) {
				// this is the last base pair
				end2 = std::make_pair( base_pairs_second_first[ helix_regions_[i][j] ], helix_regions_[i][j] );

			}
		}
		helix_end = std::make_pair( end1, end2 );
		helix_ends_.push_back( helix_end );

	}

	//std::cout << "HELIX ENDS " << helix_ends_ << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////
std::string
RNA_HelixMover::get_name() const {
	return "RNA_HelixMover";
}
///////////////////////////////////////////////////////////////////////////////////////////

} //movers
} //denovo
} //rna
} //protocols
