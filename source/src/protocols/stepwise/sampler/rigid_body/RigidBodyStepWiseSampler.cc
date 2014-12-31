// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh>
#include <protocols/toolbox/rigid_body/util.hh>
#include <protocols/stepwise/sampler/rigid_body/EulerAngles.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerOneValue.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <basic/Tracer.hh>
#include <numeric/NumericTraits.hh>

////////////////////////////////////////////////////////////////////////////////////////////
//
// Iterates through x, y, z, euler_alpha, euler_z [ = cos(beta)], and euler_gamma.
//
// Note that actual order has gamma first -- making that the inner loop allows
//  for convenient 'fast forward' past gamma or to next translation.
//
// This StepWiseSampler also gives out a Stub of the moving residues, which can be useful for screening.
//
// Originally, the class requires a 'template' residue that is centered at the origin,
//  and a reference stub, at which rotations & translations are centered.
// More recently, can just supply a pose and the moving_residue.
//
// Developed for protocols/swa/rna/StepWiseRNA_FloatingBaseSampler.
//
////////////////////////////////////////////////////////////////////////////////////////////


static thread_local basic::Tracer TR( "protocols.sampler.rigid_body.RigidBodyStepWiseSampler" );

typedef  numeric::xyzMatrix< core::Real > Matrix;
static Real const RADS_PER_DEG = numeric::NumericTraits < Real > ::pi() / 180.;

using namespace protocols::toolbox::rigid_body;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rigid_body {


	//Constructor
	// example pose must have fold tree defined, with a jump connecting a reference res to the moving res.
	RigidBodyStepWiseSampler::RigidBodyStepWiseSampler( pose::Pose const & pose,
																			Size const moving_res ):
		moving_res_( moving_res ),
	 	reference_res_( figure_out_reference_res_for_jump( pose, moving_res ) )
	{
		moving_partition_res_ = figure_out_moving_partition_res( pose, moving_res, reference_res_ );
	 	pose_at_origin_ = transform_moving_partition_to_origin( pose, moving_res_, moving_partition_res_ );
		moving_residue_at_origin_ = pose_at_origin_->residue( moving_res_ ).clone();
		reference_stub_ = initialize_stub( pose_at_origin_->residue( reference_res_ ) );
	}

	//Constructor -- old-style.
	RigidBodyStepWiseSampler::RigidBodyStepWiseSampler( Size const moving_res,
																			core::conformation::Residue const & template_moving_residue,
																			core::kinematics::Stub const & reference_stub ):
		moving_res_( moving_res ),
		moving_residue_at_origin_( template_moving_residue.clone() ),
		reference_stub_( reference_stub ),
		reference_res_( 0 ) // unknown with this old style of initialization
	{}

	//Destructor
	RigidBodyStepWiseSampler::~RigidBodyStepWiseSampler()
	{}

	////////////////////////////////////////////////////////////
	void RigidBodyStepWiseSampler::init() {

		value_range_.init();
		clear_rotamer();

		StepWiseSamplerOneValueOP x_rotamer( new StepWiseSamplerOneValue( value_range_.x_values(), "x"   ) );
		StepWiseSamplerOneValueOP y_rotamer( new StepWiseSamplerOneValue( value_range_.y_values(), "y"   ) );
		StepWiseSamplerOneValueOP z_rotamer( new StepWiseSamplerOneValue( value_range_.z_values(), "z"   ) );
		StepWiseSamplerOneValueOP euler_alpha_rotamer( new StepWiseSamplerOneValue( value_range_.euler_alpha_values(), "euler_alpha"   ) );
		StepWiseSamplerOneValueOP euler_z_rotamer( new StepWiseSamplerOneValue( value_range_.euler_z_values(), "euler_beta" ) );
		StepWiseSamplerOneValueOP euler_gamma_rotamer( new StepWiseSamplerOneValue( value_range_.euler_gamma_values(), "euler_gamma"  ) );

		// first rotamer is the inner-most loop.
		add_external_loop_rotamer( euler_gamma_rotamer );
		add_external_loop_rotamer( euler_z_rotamer );
		add_external_loop_rotamer( euler_alpha_rotamer );
		add_external_loop_rotamer( z_rotamer );
		add_external_loop_rotamer( y_rotamer );
		add_external_loop_rotamer( x_rotamer );

		StepWiseSamplerOneValueComb::init();

	}

	///////////////////////////////////////////////////////////////////////////
	void RigidBodyStepWiseSampler::apply( core::pose::Pose & pose ){
		apply( pose, id_ );
	}

	///////////////////////////////////////////////////////////////////////////
	void RigidBodyStepWiseSampler::apply( core::pose::Pose & pose,
																Size const id ) {
		runtime_assert( is_init() );
		apply_by_jump( pose, moving_res_, get_stub( id ) );
	}

	///////////////////////////////////////////////////////////////////////////
	void RigidBodyStepWiseSampler::apply( core::pose::Pose & pose,
																core::conformation::Residue const & template_moving_residue ){
		// note that this may be deprecated soon.
		apply( pose, template_moving_residue, id_ );
	}

	///////////////////////////////////////////////////////////////////////////
	void RigidBodyStepWiseSampler::apply( core::pose::Pose & pose,
																core::conformation::Residue const & template_moving_residue,
																Size const id ) {
		runtime_assert( is_init() );
		// note that this may be deprecated soon.
		transform_single_residue( pose, moving_res_,
															template_moving_residue, get_stub( id ) );
	}


	///////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::apply( core::conformation::Residue & residue_initially_at_origin ){
		Size const seqpos = residue_initially_at_origin.seqpos();

		// this function only works when RigidBodyStepWiseSampler is initialized with a pose. Don't call this otherwise.
		runtime_assert( moving_partition_res_.size() > 0 );

		if ( moving_partition_res_.has_value( seqpos  ) ){
			utility::vector1< std::pair < id::AtomID, numeric::xyzVector< core::Real > > > xyz_list;
			get_atom_coordinates( xyz_list, residue_initially_at_origin.seqpos(),
														residue_initially_at_origin, get_stub() );
			for ( Size n = 1; n <= xyz_list.size(); n++ )	residue_initially_at_origin.set_xyz( xyz_list[n].first.atomno(), xyz_list[n].second );
		}
	}

	///////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::apply( Vector & xyz_initially_at_origin, Size const seqpos ){

		// this function only works when RigidBodyStepWiseSampler is initialized with a pose. Don't call this otherwise.
		runtime_assert( moving_partition_res_.size() > 0 );
		if ( moving_partition_res_.has_value( seqpos  ) ){
			get_specific_atom_coordinate( xyz_initially_at_origin, get_stub() );
		}
	}

	///////////////////////////////////////////////////////////////////////////
	core::kinematics::Stub const &
	RigidBodyStepWiseSampler::get_stub() {
		return get_stub( id_list_ );
	}

	///////////////////////////////////////////////////////////////////////////
	core::kinematics::Stub const &
	RigidBodyStepWiseSampler::get_stub( utility::vector1< Size > const & id_list ) {

		ValueList const & rotamer_values = get_value_list( id_list );
		runtime_assert( rotamer_values.size() == 6 ); // rotations & translations

		Vector O_frame_centroid;
		O_frame_centroid.x() = ( rotamer_values[6] );
		O_frame_centroid.y() = ( rotamer_values[5] );
		O_frame_centroid.z() = ( rotamer_values[4] );
		moving_res_stub_.v = ( reference_stub_.M * O_frame_centroid ) + reference_stub_.v;

		Matrix O_frame_rotation;
		EulerAngles euler_angles;
		euler_angles.set_alpha( rotamer_values[3] );
		euler_angles.set_z(     rotamer_values[2] );
		euler_angles.set_gamma( rotamer_values[1] );
		euler_angles.convert_to_rotation_matrix( O_frame_rotation );

		moving_res_stub_.M = reference_stub_.M * O_frame_rotation;

		return moving_res_stub_;
	}


	///////////////////////////////////////////////////////////////////////////
	core::kinematics::Stub const &
	RigidBodyStepWiseSampler::get_stub( Size const id ) {
		utility::vector1< Size > id_list = id2list( id );
		return get_stub( id_list );
	}

	///////////////////////////////////////////////////////////////////////////
	utility::vector1< Real > const &
	RigidBodyStepWiseSampler::get_values() {
		return get_value_list( id_list_ );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::transform_single_residue( pose::Pose & pose, Size const & seq_num,
																							core::conformation::Residue const & rsd_at_origin,
																							core::kinematics::Stub const & moving_res_stub ){
		// [OLD] option 1 -- old school, just works for single residue -- note that it replaces
		//  that residue with the input residue.
		utility::vector1< std::pair < id::AtomID, numeric::xyzVector< core::Real > > > xyz_list;
		get_atom_coordinates( xyz_list, seq_num, rsd_at_origin, moving_res_stub );
		for ( Size n = 1; n <= xyz_list.size(); n++ )	pose.set_xyz( xyz_list[n].first, xyz_list[n].second );
	}


	///////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::apply_by_jump( pose::Pose & pose, Size const & seq_num,
																	 core::kinematics::Stub const & ){
		// [NEW] option 2 -- new, generalizable to full moving partition, but using Jump setter.
		//  requires that any changes to conformation of Residue occur /*outside*/.
		calculate_jump( pose, seq_num, moving_res_stub_ );
		pose.set_jump( jump_atom_id_, jump_ );
	}

	///////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::calculate_jump( pose::Pose & pose, Size const seq_num, kinematics::Stub const & moving_res_stub ){

		using namespace core::kinematics;
		using namespace core::id;
		using namespace core::pose;

		Size const jump_no = stepwise::modeler::look_for_unique_jump_to_moving_res( pose.fold_tree(), seq_num );
		std::string downstream_atom_name = pose.fold_tree().downstream_atom( jump_no );
		Size const i = seq_num;
		Size const j = pose.residue_type( i ).atom_index( downstream_atom_name );
		kinematics::tree::AtomCOP current_atom ( pose.atom_tree().atom_dont_do_update( AtomID(j,i) ).get_self_ptr() );
		runtime_assert( current_atom->is_jump() );

		core::kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );
		core::kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
		core::kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );

		Stub input_stub( pose.xyz( input_stub_atom1->id() ),
										 pose.xyz( input_stub_atom2->id() ),
										 pose.xyz( input_stub_atom3->id() ) );

		core::kinematics::tree::AtomCOP stub_atom1( current_atom->stub_atom1() );
		core::kinematics::tree::AtomCOP stub_atom2( current_atom->stub_atom2() );
		core::kinematics::tree::AtomCOP stub_atom3( current_atom->stub_atom3() );

		Vector stub_atom1_xyz, stub_atom2_xyz, stub_atom3_xyz;
		core::conformation::Residue const & template_rsd = moving_residue_at_origin();
		get_specific_atom_coordinate( pose.residue_type(i).atom_name( stub_atom1->id().atomno() ), stub_atom1_xyz,
																	template_rsd, moving_res_stub );
		get_specific_atom_coordinate( pose.residue_type(i).atom_name( stub_atom2->id().atomno() ), stub_atom2_xyz,
																	template_rsd, moving_res_stub );
		get_specific_atom_coordinate( pose.residue_type(i).atom_name( stub_atom3->id().atomno() ), stub_atom3_xyz,
																	template_rsd, moving_res_stub );

		Stub const stub( stub_atom1_xyz, stub_atom2_xyz, stub_atom3_xyz );

		jump_ = Jump( input_stub, stub );
		jump_atom_id_ = AtomID( j, i );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::fast_forward_to_end(){
		fast_forward( 6 );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::fast_forward_to_next_translation(){
		fast_forward( 3 ); // go to end of euler angles, so that next increment will update translation.
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::fast_forward_to_next_euler_gamma(){
		fast_forward( 1 ); // go to end of gamma
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::set_x_values( Real const centroid_x_min, Real const centroid_x_max, Real const centroid_x_bin ){
		value_range_.set_x_values( centroid_x_min, centroid_x_max, centroid_x_bin );
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::set_y_values( Real const centroid_y_min, Real const centroid_y_max, Real const centroid_y_bin ){
		value_range_.set_y_values( centroid_y_min, centroid_y_max, centroid_y_bin );
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::set_z_values( Real const centroid_z_min, Real const centroid_z_max, Real const centroid_z_bin ){
		value_range_.set_z_values( centroid_z_min, centroid_z_max, centroid_z_bin );
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::set_euler_alpha_values( Real const centroid_euler_alpha_min, Real const centroid_euler_alpha_max, Real const centroid_euler_alpha_bin ){
		value_range_.set_euler_alpha_values( centroid_euler_alpha_min, centroid_euler_alpha_max, centroid_euler_alpha_bin );
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::set_euler_z_values( Real const centroid_euler_z_min, Real const centroid_euler_z_max, Real const centroid_euler_z_bin ){
		value_range_.set_euler_z_values( centroid_euler_z_min, centroid_euler_z_max, centroid_euler_z_bin );
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyStepWiseSampler::set_euler_gamma_values( Real const centroid_euler_gamma_min, Real const centroid_euler_gamma_max, Real const centroid_euler_gamma_bin ){
		value_range_.set_euler_gamma_values( centroid_euler_gamma_min, centroid_euler_gamma_max, centroid_euler_gamma_bin );
	}

	core::pose::PoseCOP
	RigidBodyStepWiseSampler::pose_at_origin(){ return pose_at_origin_; }

	core::conformation::Residue const &
	RigidBodyStepWiseSampler::get_residue_at_origin( Size const seqpos ){ return pose_at_origin_->residue( seqpos ) ; }

	//////////////////////////////////////////////////////////////////////////
	/// @brief Name of the class
	std::string
	RigidBodyStepWiseSampler::get_name() const {
		return "RigidBodyStepWiseSampler residue:" + utility::to_string(moving_res_);
	}


} //rigid_body
} //sampler
} //stepwise
} //protocols
