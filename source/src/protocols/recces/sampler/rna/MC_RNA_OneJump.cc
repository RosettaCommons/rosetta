// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/rna/MC_RNA_OneJump.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/recces/sampler/rna/MC_RNA_OneJump.hh>
#include <core/chemical/rna/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

#include <numeric/random/random_xyz.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/constants.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.recces.sampler.rna.MC_RNA_OneJump" );

using namespace core;

namespace protocols {
namespace recces {
namespace sampler {
namespace rna {

	//Constructor
	MC_RNA_OneJump::MC_RNA_OneJump( core::pose::Pose const & pose,
																	core::Size const & jump_num ):
		MC_Sampler(),
		jump_num_( jump_num ),
		rmsd_cutoff_( 2.0 )
	{
		using namespace core::pose;
		using namespace core::kinematics;

		// just get out the two residues connected by the jump from the pose...
		scratch_pose_ = PoseOP( new Pose );
		FoldTree const & f( pose.fold_tree() );
		scratch_pose_->append_residue_by_bond( *( pose.residue( f.upstream_jump_residue( jump_num ) ).clone() ) );
		scratch_pose_->append_residue_by_jump( *( pose.residue( f.downstream_jump_residue( jump_num ) ).clone() ), 1,
																					 f.upstream_atom( jump_num ),
																					 f.downstream_atom( jump_num ),
																					 true /* new chain */ );
		ref_scratch_pose_ = scratch_pose_->clone();

		stored_upstream_stub_ = scratch_pose_->conformation().upstream_jump_stub( 1 ); // only 1 jump in 2-residue scratch pose
		active_jump_ = scratch_pose_->jump( 1 );
		update();

		set_name( "MC_RNA_OneJump" );
	}

	//Destructor
	MC_RNA_OneJump::~MC_RNA_OneJump()
	{}

	///////////////////////////////////////////////////////////////////////////
	void
	MC_RNA_OneJump::operator++()
	{
		using namespace numeric;
		using namespace numeric::random;

		// need to separately define uniform_modeler() vs. markov-chain -- currently doing markov chain:
		active_jump_ = stored_jump_;

		// following few lines are from protocols/rigid/UniformRigidBodyMover.cc -- I'm not sure if
		//  its "kosher", however. Would be better to write myself in terms of, say, Euler angles or quaternions.
		core::Real theta = random_rotation_angle<core::Real>( rotation_mag_, numeric::random::rg() );
		xyzVector<core::Real> axis = random_point_on_unit_sphere<core::Real>( numeric::random::rg() );
		xyzMatrix<core::Real> delta_rot = rotation_matrix_radians( axis, theta );
		active_jump_.rotation_by_matrix( stored_upstream_stub_, stored_base_centroid_, delta_rot );

		// translation
		xyzVector<core::Real> delta_trans = random_translation( translation_mag_, numeric::random::rg() );
		active_jump_.set_translation( active_jump_.get_translation() +  delta_trans );

		found_move_ = true;
		if ( !check_jump_in_range() ) {
			active_jump_ = stored_jump_;
			found_move_ = false;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	void
	MC_RNA_OneJump::apply( core::pose::Pose & pose )
	{
		pose.set_jump( jump_num_, active_jump_ );
	}
	///////////////////////////////////////////////////////////////////////////
	void
	MC_RNA_OneJump::update()
	{
		if ( update_pose_ != 0 ) {
			active_jump_ = update_pose_->jump( jump_num_ );
			scratch_pose_->set_jump( 1, active_jump_ );
		}
		stored_jump_ = active_jump_;
		stored_base_centroid_ = core::chemical::rna::get_rna_base_centroid( scratch_pose_->residue( 2 /*moving_rsd_*/ ) );
	}


	// copy/pasted from rb_entropy -- unify into util
	core::Real
	calc_base_centroid_rmsd( core::conformation::Residue const & rsd1, core::conformation::Residue const & rsd2 )
	{
		Real rmsd( 0.0 ); Size numatoms( 0 );
		for ( Size i = rsd1.first_sidechain_atom() + 1; i <= rsd1.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2prime oxygen.
			if ( rsd1.is_virtual( i ) ) continue;
			if ( rsd1.is_repulsive( i ) ) continue;
			Vector dist = ( rsd1.xyz(i) - rsd2.xyz(i) );
			rmsd += dist.length_squared();
			numatoms++;
		}
		rmsd = sqrt( rmsd/numatoms );
		return rmsd;
	}


	///////////////////////////////////////////////////////////////////////////
	// Note -- this could be made *EVEN FASTER* if we use moments-of-inertia
	//   trick, and also, we don't need to catch ref_scratch_pose -- just its atoms.
	bool
	MC_RNA_OneJump::check_jump_in_range()
	{
		Size const moving_rsd( 2 ); // always second residue in scratch poses
		scratch_pose_->set_jump( 1, active_jump_ );
		Real rmsd( calc_base_centroid_rmsd( scratch_pose_->residue(moving_rsd),
																				ref_scratch_pose_->residue(moving_rsd) ) );
		return ( rmsd <= rmsd_cutoff_ );
	}

	///////////////////////////////////////////////////////////////////////////
	void
	MC_RNA_OneJump::show( std::ostream & out, Size const indent ) const {
		for ( Size n = 1; n <= indent; n++ ) out << ' ';
		out << get_name() << " jump_number:" << jump_num_ << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////
	MC_SamplerOP
	MC_RNA_OneJump::find( core::id::TorsionID const & torsion_id ) {
		if ( torsion_id.rsd() == jump_num_  &&
				 torsion_id.type() == core::id::TorsionType::JUMP ) {
			return std::dynamic_pointer_cast< MC_RNA_OneJump >( shared_from_this() );
		}
		return 0;
	}


} //rna
} //sampler
} //recces
} //protocols
