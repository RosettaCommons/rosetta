// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_rigid_body_RigidBodyRotamer_HH
#define INCLUDED_protocols_rotamer_sampler_rigid_body_RigidBodyRotamer_HH

#include <protocols/rotamer_sampler/RotamerOneValueComb.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {

	class RigidBodyRotamer: public RotamerOneValueComb {

	public:

		//constructor
		RigidBodyRotamer( Size const moving_res,
																				core::conformation::Residue const & template_moving_residue,
																				core::kinematics::Stub const & reference_stub );

		//destructor
		~RigidBodyRotamer();

	public:

		void init();

		virtual void apply( core::pose::Pose & pose );

		virtual void apply( core::pose::Pose & pose,
												Size const id );

		void apply( core::pose::Pose & pose,
								core::conformation::Residue const & template_moving_residue );

		void apply( core::pose::Pose & pose,
								core::conformation::Residue const & template_moving_residue,
								Size const id );

		core::kinematics::Stub const &
		get_stub();

		core::kinematics::Stub const &
		get_stub( utility::vector1< Size > const & id_list);

		core::kinematics::Stub const &
		get_stub( Size const id );

		void
		fast_forward_to_next_translation();

		void
		fast_forward_to_next_euler_gamma();

		// @brief set & get functions
		void set_x_values( Real const centroid_x_min, Real const centroid_x_max, Real const centroid_x_bin );
		void set_y_values( Real const centroid_y_min, Real const centroid_y_max, Real const centroid_y_bin );
		void set_z_values( Real const centroid_z_min, Real const centroid_z_max, Real const centroid_z_bin );
		void set_euler_alpha_values( Real const centroid_euler_alpha_min, Real const centroid_euler_alpha_max, Real const centroid_euler_alpha_bin );
		void set_euler_z_values( Real const centroid_euler_z_min, Real const centroid_euler_z_max, Real const centroid_euler_z_bin );
		void set_euler_gamma_values( Real const centroid_euler_gamma_min, Real const centroid_euler_gamma_max, Real const centroid_euler_gamma_bin );

		ValueList const & get_values();

	private:

		void
		set_sampler_values( Real const & val_min, Real const & val_max, Real const & val_bin, ValueList & values );

		void
		set_coordinate_frame( pose::Pose & pose, Size const & seq_num, core::conformation::Residue const & rsd_at_origin, core::kinematics::Stub const & moving_res_stub );

	private:

		Size const moving_res_;
		conformation::Residue const & template_moving_residue_;
		kinematics::Stub const & reference_stub_;
		kinematics::Stub moving_res_stub_;
		int centroid_bin_min_;
		int centroid_bin_max_;
		Real centroid_bin_size_;
		int euler_angle_bin_min_;
		int euler_angle_bin_max_;
		Real euler_angle_bin_size_;
		int euler_z_bin_min_;
		int euler_z_bin_max_;
		Real euler_z_bin_size_;

		ValueList x_values_, y_values_, z_values_;
		ValueList euler_alpha_values_, euler_z_values_, euler_gamma_values_;
	};

} //rigid_body
} //rotamer_sampler
} //protocols

#endif
