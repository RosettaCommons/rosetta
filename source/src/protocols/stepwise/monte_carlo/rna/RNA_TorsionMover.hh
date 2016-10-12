// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RNA_TorsionMover.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_monte_carlo_RNA_TorsionMover_hh
#define INCLUDED_protocols_stepwise_monte_carlo_RNA_TorsionMover_hh

#include <core/pose/Pose.fwd.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_TorsionMover.fwd.hh>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class RNA_TorsionMover: public protocols::moves::Mover {
public:


	//destructor
	~RNA_TorsionMover();

	// Undefinded, commenting out to fix PyRosetta build  void apply( core::pose::Pose & pose, core::Size const res_to_torsion, protocols::stepwise::monte_carlo::MovingResidueCase const moving_residue_case  );

	/// @brief Apply the minimizer to one pose
	using protocols::moves::Mover::apply;
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	void
	apply( core::pose::Pose & pose, std::string & move_type, core::Real const & sample_range );

	void
	random_torsion_move( core::pose::Pose & pose,
		utility::vector1< core::Size > const & moving_res_list,
		std::string & move_type,
		core::Real const & sample_range );

	void
	sample_near_suite_torsion(utility::vector1< core::Real > & torsion_list, core::Real const stddev);

	void
	sample_near_nucleoside_torsion(utility::vector1< core::Real > & torsion_list, core::Real const stddev);

	void
	apply_random_nucleoside_torsion( core::pose::Pose & pose,
		core::Size const moving_res );

	void
	apply_random_suite_torsion( core::pose::Pose & pose,
		core::Size const moving_suite );

	void
	apply_nucleoside_torsion_Aform(
		core::pose::Pose & pose,
		core::Size const moving_res );

	void
	apply_suite_torsion_Aform(
		core::pose::Pose & pose,
		core::Size const moving_suite );

	void
	sample_near_suite_torsion( core::pose::Pose & pose, core::Size const moving_suite, core::Real const sample_range);

	void
	sample_near_nucleoside_torsion( core::pose::Pose & pose, core::Size const moving_res, core::Real const sample_range);

	void
	crankshaft_alpha_gamma( core::pose::Pose & pose, core::Size const moving_suite, core::Real const sample_range);


private:

	void
	apply_nucleoside_torsion( utility::vector1< core::Real > const & torsion_set,
		core::pose::Pose & pose,
		core::Size const moving_res);


	void
	apply_suite_torsion( utility::vector1< core::Real > const & torsion_set,
		core::pose::Pose & pose,
		core::Size const moving_suite );

	utility::vector1< core::Real>
	get_suite_torsion( core::pose::Pose const & pose, core::Size const moving_suite );

	utility::vector1< core::Real>
	get_nucleoside_torsion( core::pose::Pose const & pose, core::Size const moving_nucleoside );

private:

	core::Size default_sample_range_;
	core::chemical::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info_;

};

} //rna
} //monte_carlo
} //stepwise
} //protocols

#endif
