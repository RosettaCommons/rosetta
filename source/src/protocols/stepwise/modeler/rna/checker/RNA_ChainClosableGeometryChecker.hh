// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_checker_RNA_ChainClosableGeometryChecker_HH
#define INCLUDED_protocols_stepwise_rna_checker_RNA_ChainClosableGeometryChecker_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {

class RNA_ChainClosableGeometryChecker: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_ChainClosableGeometryChecker( core::Size const five_prime_chain_break_res, core::Size const gap_size );

	RNA_ChainClosableGeometryChecker( core::Size const five_prime_chain_break_res, core::Size const three_prime_chain_break_res, core::Size const gap_size );

	//destructor
	~RNA_ChainClosableGeometryChecker();

public:

	bool
	check_screen( core::pose::Pose const & pose,
		bool const strict = false ) const;

	bool
	check_screen( core::pose::Pose const & moving_pose,
		core::pose::Pose const & reference_pose,
		bool const is_prepend,
		bool const strict = false ) const;

	bool
	check_screen(  core::pose::Pose const & pose,
		utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
		core::kinematics::Stub const & moving_res_base_stub,
		core::Size const & reference_res ) const;

	bool
	check_screen(  utility::vector1< core::pose::PoseOP > const & pose_data_list,
		utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
		core::kinematics::Stub const & moving_res_base_stub,
		core::Size const & reference_res ) const;

	bool
	check_chain_closable_geometry( numeric::xyzVector< core::Real > const & xyz_1, numeric::xyzVector< core::Real > const & xyz_2 ) const;

	bool
	check_chain_closable_geometry( core::conformation::Residue const & five_prime_residue,
		core::conformation::Residue const & three_prime_residue ) const;

	core::Size const & five_prime_chain_break_res() const { return five_prime_chain_break_res_; }
	core::Size const & three_prime_chain_break_res() const { return three_prime_chain_break_res_; }

	core::Real const & dist_squared() const { return dist_squared_; }
	core::Real const & max_dist_squared() const { return max_dist_squared_; }

private:

	bool
	check_chain_closable_geometry( core::pose::Pose const & five_prime_pose, core::pose::Pose const & three_prime_pose ) const;

	bool
	check_chain_closable_geometry( core::pose::Pose const & five_prime_pose, core::pose::Pose const & three_prime_pose, bool const strict ) const;

	bool
	check_chain_closable_geometry( core::Size const & reference_res,
		utility::vector1< core::pose::PoseOP > const & pose_data_list,
		utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
		core::kinematics::Stub const & moving_res_base_stub,
		bool const is_prepend ) const;
	bool
	check_chain_closable_geometry( core::Size const & reference_res,
		core::pose::Pose const & pose,
		utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list, //this one correspond to the moving_base
		core::kinematics::Stub const & moving_res_base_stub,
		bool const is_prepend ) const;

	bool
	check_chain_closable_geometry_strict( core::pose::Pose const & five_prime_pose, core::pose::Pose const & three_prime_pose ) const;

	void
	get_C4_C3_distance_range( core::conformation::Residue const & five_prime_rsd,
		core::conformation::Residue const & three_prime_rsd,
		core::Distance & C4_C3_dist_min,
		core::Distance & C4_C3_dist_max ) const;

	void
	initialize_distance_range();

private:

	core::Size const five_prime_chain_break_res_;
	core::Size const three_prime_chain_break_res_;
	core::Size const gap_size_;
	core::Real min_dist_squared_;
	core::Real max_dist_squared_;
	mutable core::Real dist_squared_;

};

} //checker
} //rna
} //modeler
} //stepwise
} //protocols

#endif
