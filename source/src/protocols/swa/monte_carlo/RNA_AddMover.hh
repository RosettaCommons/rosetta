// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_AddMover.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_monte_carlo_RNA_AddMover_hh
#define INCLUDED_protocols_swa_monte_carlo_RNA_AddMover_hh

#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/swa/monte_carlo/SWA_Move.hh>
#include <protocols/swa/monte_carlo/RNA_TorsionMover.fwd.hh>
#include <protocols/swa/monte_carlo/RNA_AddMover.fwd.hh>


namespace protocols {
namespace swa {
namespace monte_carlo {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class RNA_AddMover: public protocols::moves::Mover {
public:


	RNA_AddMover( core::scoring::ScoreFunctionOP scorefxn );

	//destructor
	~RNA_AddMover();

	using protocols::moves::Mover::apply;
  void
	apply( core::pose::Pose & pose, Size const res_to_add, protocols::swa::monte_carlo::MovingResidueCase const moving_residue_case  );

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	void set_kT( core::Real const & setting ){ kT_ = setting; }

	void set_sample_range_small( core::Real const setting ){ sample_range_small_ = setting; }

	void set_sample_range_large( core::Real const setting ){ sample_range_large_ = setting; }

	void set_internal_cycles( Size const setting ){ internal_cycles_ = setting; }

	void set_presample_added_residue( Size const setting ){ presample_added_residue_ = setting; }

	void set_presample_by_swa( Size const setting ){ presample_by_swa_ = setting; }

	void set_start_added_residue_in_aform( Size const setting ){ start_added_residue_in_aform_ = setting; }

	void set_minimize_single_res( Size const setting ){ minimize_single_res_ = setting; }

	void set_num_random_samples( Size const & setting ){ num_random_samples_ = setting; }

	void set_use_phenix_geo( Size const & setting ){ use_phenix_geo_ = setting; };

	void set_erraser( Size const & setting ){ erraser_ = setting; };

private:

	void sample_by_swa( core::pose::Pose & pose, Size const res_to_add  ) const;

	void sample_by_monte_carlo_internal( core::pose::Pose &  pose, Size const nucleoside_num, Size const suite_num ) const;

private:

	core::chemical::ResidueTypeSetCAP rsd_set_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool presample_added_residue_;
	bool presample_by_swa_;
	bool minimize_single_res_;
	bool start_added_residue_in_aform_;
	bool use_phenix_geo_;
	bool erraser_;
	Size internal_cycles_;
	RNA_TorsionMoverOP rna_torsion_mover_;
	core::Real sample_range_small_;
	core::Real sample_range_large_;
	core::Real kT_;
	Size num_random_samples_;
};

} // monte_carlo
} // swa
} // protocols

#endif
