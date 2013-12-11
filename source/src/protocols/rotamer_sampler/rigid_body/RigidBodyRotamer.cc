// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyUtil.hh>
#include <protocols/rotamer_sampler/rigid_body/EulerAngles.hh>
#include <protocols/rotamer_sampler/RotamerOneValue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <numeric/NumericTraits.hh>
#include <basic/Tracer.hh>

////////////////////////////////////////////////////////////////////////////////////////////
//
// Iterates through x, y, z, euler_alpha, euler_z [ = cos(beta)], and euler_gamma.
//
// Note that actual order has gamma first -- making that the inner loop allows
//  for convenient 'fast forward' past gamma or to next translation.
//
// This Rotamer also gives out a Stub of the moving residues, which can be useful for screening.
// The class requires a 'template' residue that is centered at the origin,
//  and a reference stub, at which rotations & translations are centered.
// Developed for protocols/swa/rna/StepWiseRNA_FloatingBaseSampler.
//
////////////////////////////////////////////////////////////////////////////////////////////


static basic::Tracer TR( "protocols.rotamer_sampler.rigid_body.RigidBodyRotamer" );

typedef  numeric::xyzMatrix< core::Real > Matrix;
static Real const RADS_PER_DEG = numeric::NumericTraits < Real > ::pi() / 180.;

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {

	//Constructor
	RigidBodyRotamer::RigidBodyRotamer( Size const moving_res,
																			core::conformation::Residue const & template_moving_residue,
																			core::kinematics::Stub const & reference_stub ):
		moving_res_( moving_res ),
		template_moving_residue_( template_moving_residue ),
		reference_stub_( reference_stub )
	{}

	//Destructor
	RigidBodyRotamer::~RigidBodyRotamer()
	{}

	////////////////////////////////////////////////////////////
	void RigidBodyRotamer::init() {

		clear_rotamer();

		RotamerOneValueOP x_rotamer = new RotamerOneValue( x_values_   );
		RotamerOneValueOP y_rotamer = new RotamerOneValue( y_values_   );
		RotamerOneValueOP z_rotamer = new RotamerOneValue( z_values_   );
		RotamerOneValueOP euler_alpha_rotamer = new RotamerOneValue( euler_alpha_values_   );
		RotamerOneValueOP euler_z_rotamer     = new RotamerOneValue( euler_z_values_ );
		RotamerOneValueOP euler_gamma_rotamer = new RotamerOneValue( euler_gamma_values_   );

		// first rotamer is the inner-most loop.
		add_rotamer( euler_gamma_rotamer );
		add_rotamer( euler_z_rotamer );
		add_rotamer( euler_alpha_rotamer );
		add_rotamer( z_rotamer );
		add_rotamer( y_rotamer );
		add_rotamer( x_rotamer );

		RotamerOneValueComb::init();

	}

	///////////////////////////////////////////////////////////////////////////
	void RigidBodyRotamer::apply( core::pose::Pose & pose ){
		apply( pose, template_moving_residue_, id_ );
	}

	///////////////////////////////////////////////////////////////////////////
	void RigidBodyRotamer::apply( core::pose::Pose & pose,
																core::conformation::Residue const & template_moving_residue ){
		apply( pose, template_moving_residue, id_ );
	}

	///////////////////////////////////////////////////////////////////////////
	void RigidBodyRotamer::apply( core::pose::Pose & pose,
																Size const id ) {
		apply( pose, template_moving_residue_, id );
	}

	///////////////////////////////////////////////////////////////////////////
	void RigidBodyRotamer::apply( core::pose::Pose & pose,
																core::conformation::Residue const & template_moving_residue,
																Size const id ) {
		runtime_assert( is_init() );
		set_coordinate_frame( pose, moving_res_, template_moving_residue, get_stub( id ) );
	}

	///////////////////////////////////////////////////////////////////////////
	core::kinematics::Stub const &
	RigidBodyRotamer::get_stub() {
		return get_stub( id_list_ );
	}

	///////////////////////////////////////////////////////////////////////////
	core::kinematics::Stub const &
	RigidBodyRotamer::get_stub( utility::vector1< Size > const & id_list ) {

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
	RigidBodyRotamer::get_stub( Size const id ) {
		utility::vector1< Size > id_list = id2list( id );
		return get_stub( id_list );
	}

	///////////////////////////////////////////////////////////////////////////
	utility::vector1< Real > const &
	RigidBodyRotamer::get_values() {
		return get_value_list( id_list_ );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyRotamer::set_coordinate_frame( pose::Pose & pose, Size const & seq_num, core::conformation::Residue const & rsd_at_origin, core::kinematics::Stub const & moving_res_stub ){

		utility::vector1< std::pair < id::AtomID, numeric::xyzVector< core::Real > > > xyz_list;
		get_atom_coordinates( xyz_list, seq_num, rsd_at_origin, moving_res_stub );

		for ( Size n = 1; n <= xyz_list.size(); n++ ){
			pose.set_xyz( xyz_list[n].first, xyz_list[n].second );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyRotamer::fast_forward_to_next_translation(){
		fast_forward( 3 ); // go to end of translations.
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyRotamer::fast_forward_to_next_euler_gamma(){
		fast_forward( 1 ); // go to end of gamma
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyRotamer::set_sampler_values( Real const & val_min, Real const & val_max, Real const & val_bin, utility::vector1< Real > & values ){
		values.clear();
		runtime_assert( val_min < val_max );
		for ( Real val = val_min; val <= val_max + 1.0e-6 /*floating point error*/; val += val_bin ) values.push_back( val );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyRotamer::set_x_values( Real const centroid_x_min, Real const centroid_x_max, Real const centroid_x_bin ){
		set_sampler_values( centroid_x_min, centroid_x_max, centroid_x_bin, x_values_ );
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyRotamer::set_y_values( Real const centroid_y_min, Real const centroid_y_max, Real const centroid_y_bin ){
		set_sampler_values( centroid_y_min, centroid_y_max, centroid_y_bin, y_values_ );
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyRotamer::set_z_values( Real const centroid_z_min, Real const centroid_z_max, Real const centroid_z_bin ){
		set_sampler_values( centroid_z_min, centroid_z_max, centroid_z_bin, z_values_ );
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyRotamer::set_euler_alpha_values( Real const centroid_euler_alpha_min, Real const centroid_euler_alpha_max, Real const centroid_euler_alpha_bin ){
		set_sampler_values( RADS_PER_DEG * centroid_euler_alpha_min, RADS_PER_DEG * centroid_euler_alpha_max, RADS_PER_DEG * centroid_euler_alpha_bin, euler_alpha_values_ );
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyRotamer::set_euler_z_values( Real const centroid_euler_z_min, Real const centroid_euler_z_max, Real const centroid_euler_z_bin ){
		set_sampler_values( centroid_euler_z_min, centroid_euler_z_max, centroid_euler_z_bin, euler_z_values_ );
		// need to really get these in bounds...
		for ( Size n = 1; n <= euler_z_values_.size(); n++ ){
			euler_z_values_[n] = std::max( euler_z_values_[n], -1.0 );
			euler_z_values_[n] = std::min( euler_z_values_[n], +1.0 );
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodyRotamer::set_euler_gamma_values( Real const centroid_euler_gamma_min, Real const centroid_euler_gamma_max, Real const centroid_euler_gamma_bin ){
		set_sampler_values( RADS_PER_DEG * centroid_euler_gamma_min, RADS_PER_DEG * centroid_euler_gamma_max, RADS_PER_DEG * centroid_euler_gamma_bin, euler_gamma_values_ );
	}


} //rigid_body
} //rotamer_sampler
} //protocols
