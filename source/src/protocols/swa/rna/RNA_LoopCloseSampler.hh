// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_LoopCloseSampler
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_rna_RNA_LoopCloseSampler_HH
#define INCLUDED_protocols_swa_rna_RNA_LoopCloseSampler_HH

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>


namespace protocols {
namespace swa {
namespace rna {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class RNA_LoopCloseSampler: public protocols::moves::Mover {
public:

	//constructor!
	RNA_LoopCloseSampler ( core::Size const moving_suite, core::Size const chainbreak_suite );

	//destructor -- necessary?
	~RNA_LoopCloseSampler();

	/// @brief Apply the minimizer to one pose
	virtual void apply ( core::pose::Pose & pose_to_visualize );

	/// @brief Each derived class must specify its name.  The class name.
	virtual std::string get_name() const {
		return "RNA_LoopCloseSampler";
	}
	utility::vector1< utility::vector1< core::Real > >  all_torsion_info() {
		return all_torsion_info_;
	}
	Size n_construct() {
		return n_construct_;
	}
	void fill_pose ( core::pose::Pose & pose, core::Size construct_number );
	void clear_all();

	void set_bin_size ( core::Size const setting ) {
		bin_size_ = setting;
	}

	void set_epsilon_range ( core::Real const setting ) {
		epsilon_range_ = setting;
	}

	void set_rep_cutoff ( core::Real const rep_cutoff ) {
		rep_cutoff_ = rep_cutoff;
	}
	void set_torsion_range ( core::Real const torsion_range ) {
		torsion_range_ = torsion_range;
	}
	void set_torsion_increment ( core::Real const torsion_increment ) {
		torsion_increment_ = torsion_increment;
	}

	void set_just_output_score ( bool const setting ) {
		just_output_score_ = setting;
	}
	void set_sample_only ( bool const setting ) {
		sample_only_ = setting;
	}

	void set_sample_native_torsion ( bool const setting ) {
		sample_native_torsion_ = setting;
	}

	void set_include_current ( bool const setting ) {
		include_current_ = setting;
	}

	void set_choose_random ( bool const setting ) {
		choose_random_ = setting;
	}

	void set_rbs_new_pair ( utility::vector1< core::Real > const & rbs_new_pair ) {
		rbs_new_pair_ = rbs_new_pair;
	}

	void set_scorefxn ( core::scoring::ScoreFunctionOP const & scorefxn );
	// Undefined, commenting out to fix PyRosetta build  Size get_n_construct();


private:

	bool
	torsion_angles_within_cutoffs ( core::pose::Pose const & pose,
	                                core::Size const moving_suite,
	                                core::Size const chainbreak_suite );

	void initialize_rep_scorefxn();

	core::Real
	initialize_fa_rep ( core::pose::Pose const & pose,
	                    utility::vector1< Size > const & moving_suites,
	                    core::scoring::ScoreFunctionOP rep_scorefxn );

	bool
	check_clash ( core::pose::Pose & pose,
	              core::Real const & fa_rep_score_baseline,
	              core::Real const & rep_cutoff_,
	              core::scoring::ScoreFunctionOP rep_scorefxn );

	void
	add_driver_torsion_set(
												 utility::vector1< utility::vector1< core::Real > > & driver_torsion_sets,
												 utility::vector1< core::id::TorsionID > const & driver_torsion_IDs,
												 core::pose::Pose const & pose );

private:

	core::Size const moving_suite_, chainbreak_suite_;

	core::scoring::ScoreFunctionOP scorefxn_, rep_scorefxn_;
	core::Size bin_size_, n_construct_;
	core::Real rep_cutoff_, torsion_range_, torsion_increment_, epsilon_range_;
	bool just_output_score_, sample_only_, include_current_, sample_native_torsion_;
	bool choose_random_; // new 2013
	utility::vector1< core::Real > rbs_new_pair_;

	utility::vector1< core::Real > torsion_info_;
	utility::vector1< utility::vector1< core::Real > > all_torsion_info_;
};

}
} //swa
} // protocols

#endif
