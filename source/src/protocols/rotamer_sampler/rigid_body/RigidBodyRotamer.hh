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
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerValueRange.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {

	class RigidBodyRotamer: public RotamerOneValueComb {

	public:

		//constructor
		RigidBodyRotamer( pose::Pose const & pose,
											Size const moving_res );

		//constructor
		RigidBodyRotamer( Size const moving_res,
											core::conformation::Residue const & moving_residue_at_origin,
											core::kinematics::Stub const & reference_stub );

		//destructor
		~RigidBodyRotamer();

	public:

		void init();

		virtual void apply( core::pose::Pose & pose );

		virtual void apply( core::pose::Pose & pose,
												Size const id );

		void apply( core::pose::Pose & pose,
								core::conformation::Residue const & moving_residue_at_origin );

		void apply( core::pose::Pose & pose,
								core::conformation::Residue const & moving_residue_at_origin,
								Size const id );


		void apply( core::conformation::Residue & residue_initially_at_origin );

		void
		apply( Vector & xyz_initially_at_origin, Size const seqpos );

		core::kinematics::Stub const &
		get_stub();

		core::kinematics::Stub const &
		get_stub( utility::vector1< Size > const & id_list);

		core::kinematics::Stub const &
		get_stub( Size const id );

		core::conformation::Residue const &
		get_residue_at_origin( Size const seqpos );

		core::pose::PoseCOP pose_at_origin();
		void
		fast_forward_to_end();

		void
		fast_forward_to_next_translation();

		void
		fast_forward_to_next_euler_gamma();

		ValueList const & get_values();

		Size const & moving_res() const { return moving_res_; }

		core::conformation::Residue const & moving_residue_at_origin() const { return *moving_residue_at_origin_; }

		/// @brief Type of class (see enum in RotamerTypes.hh)
		virtual RotamerType type() const { return RIGID_BODY; }

		Size const & reference_res() const { return reference_res_; }
		utility::vector1< Size > const & moving_partition_res() const { return moving_partition_res_;}

		RigidBodyRotamerValueRange & value_range(){ return value_range_; }

		void set_x_values( Real const centroid_x_min, Real const centroid_x_max, Real const centroid_x_bin );
		void set_y_values( Real const centroid_y_min, Real const centroid_y_max, Real const centroid_y_bin );
		void set_z_values( Real const centroid_z_min, Real const centroid_z_max, Real const centroid_z_bin );
		void set_euler_alpha_values( Real const centroid_euler_alpha_min, Real const centroid_euler_alpha_max, Real const centroid_euler_alpha_bin );
		void set_euler_z_values( Real const centroid_euler_z_min, Real const centroid_euler_z_max, Real const centroid_euler_z_bin );
		void set_euler_gamma_values( Real const centroid_euler_gamma_min, Real const centroid_euler_gamma_max, Real const centroid_euler_gamma_bin );

	private:

		void
		calculate_jump( pose::Pose & pose, Size const seq_num, kinematics::Stub const & moving_res_stub );

		void
		transform_single_residue( pose::Pose & pose, Size const & seq_num,
															core::conformation::Residue const & rsd_at_origin, core::kinematics::Stub const & moving_res_stub );

		void
		apply_by_jump( pose::Pose & pose, Size const & seq_num,
									 core::kinematics::Stub const & moving_res_stub );

	private:

		Size const moving_res_;
		conformation::ResidueCOP moving_residue_at_origin_;
		kinematics::Stub reference_stub_;
		kinematics::Stub moving_res_stub_;

		RigidBodyRotamerValueRange value_range_;

		kinematics::Jump jump_;
		id::AtomID jump_atom_id_;

		// following is only set when RigidBodyRotamer is initialized with a pose.
		Size const reference_res_;
		utility::vector1< Size > moving_partition_res_;
		pose::PoseCOP pose_at_origin_;

	};

} //rigid_body
} //rotamer_sampler
} //protocols

#endif
