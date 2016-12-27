// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_rigid_body_RigidBodyStepWiseSampler_HH
#define INCLUDED_protocols_sampler_rigid_body_RigidBodyStepWiseSampler_HH

#include <protocols/stepwise/sampler/StepWiseSamplerOneValueComb.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.fwd.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerValueRange.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rigid_body {

class RigidBodyStepWiseSampler: public StepWiseSamplerOneValueComb {

public:

	//constructor
	RigidBodyStepWiseSampler( core::pose::Pose const & pose,
		core::Size const moving_res );

	//constructor
	RigidBodyStepWiseSampler( core::Size const moving_res,
		core::conformation::Residue const & moving_residue_at_origin,
		core::kinematics::Stub const & reference_stub );

	//destructor
	~RigidBodyStepWiseSampler();

public:

	void init();

	virtual void apply( core::pose::Pose & pose );

	virtual void apply( core::pose::Pose & pose,
		core::Size const id );

	void apply( core::pose::Pose & pose,
		core::conformation::Residue const & moving_residue_at_origin );

	void apply( core::pose::Pose & pose,
		core::conformation::Residue const & moving_residue_at_origin,
		core::Size const id );


	void apply( core::conformation::Residue & residue_initially_at_origin );

	void
	apply( core::Vector & xyz_initially_at_origin, core::Size const seqpos );

	core::kinematics::Stub const &
	get_stub();

	core::kinematics::Stub const &
	get_stub( utility::vector1< core::Size > const & id_list);


	core::kinematics::Stub const &
	get_stub( core::Size const id );

	core::kinematics::Stub const &
	reference_stub(){ return reference_stub_; }

	core::conformation::Residue const &
	get_residue_at_origin( core::Size const seqpos );

	core::pose::PoseCOP pose_at_origin();
	void
	fast_forward_to_end();

	void
	fast_forward_to_next_translation();

	void
	fast_forward_to_next_euler_gamma();

	ValueList const & get_values();

	core::Size const & moving_res() const { return moving_res_; }

	core::conformation::Residue const & moving_residue_at_origin() const { return *moving_residue_at_origin_; }

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::RIGID_BODY; }

	core::Size const & reference_res() const { return reference_res_; }
	utility::vector1< core::Size > const & moving_partition_res() const { return moving_partition_res_;}

	RigidBodyStepWiseSamplerValueRange & value_range(){ return value_range_; }

	void set_x_values( core::Real const centroid_x_min, core::Real const centroid_x_max, core::Real const centroid_x_bin );
	void set_y_values( core::Real const centroid_y_min, core::Real const centroid_y_max, core::Real const centroid_y_bin );
	void set_z_values( core::Real const centroid_z_min, core::Real const centroid_z_max, core::Real const centroid_z_bin );
	void set_euler_alpha_values( core::Real const centroid_euler_alpha_min, core::Real const centroid_euler_alpha_max, core::Real const centroid_euler_alpha_bin );
	void set_euler_z_values( core::Real const centroid_euler_z_min, core::Real const centroid_euler_z_max, core::Real const centroid_euler_z_bin );
	void set_euler_gamma_values( core::Real const centroid_euler_gamma_min, core::Real const centroid_euler_gamma_max, core::Real const centroid_euler_gamma_bin );

	/// @brief Name of the class
	virtual std::string get_name() const;

private:

	void
	calculate_jump( core::pose::Pose & pose, core::Size const seq_num, core::kinematics::Stub const & moving_res_stub );

	void
	transform_single_residue( core::pose::Pose & pose, core::Size const & seq_num,
		core::conformation::Residue const & rsd_at_origin, core::kinematics::Stub const & moving_res_stub );

	void
	apply_by_jump( core::pose::Pose & pose, core::Size const & seq_num,
		core::kinematics::Stub const & moving_res_stub );

private:

	core::Size const moving_res_;
	core::conformation::ResidueCOP moving_residue_at_origin_;
	core::kinematics::Stub reference_stub_;
	core::kinematics::Stub moving_res_stub_;

	RigidBodyStepWiseSamplerValueRange value_range_;

	core::kinematics::Jump jump_;
	core::id::AtomID jump_atom_id_;

	// following is only set when RigidBodyStepWiseSampler is initialized with a pose.
	core::Size const reference_res_;
	utility::vector1< core::Size > moving_partition_res_;
	core::pose::PoseCOP pose_at_origin_;

};

} //rigid_body
} //sampler
} //stepwise
} //protocols

#endif
