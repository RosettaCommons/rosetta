// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/mpi_refinement/StructAvrgMover.hh
/// @brief
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_mpi_refinement_StructAvrgMover_hh
#define INCLUDED_protocols_mpi_refinement_StructAvrgMover_hh

// Package headers
#include <protocols/wum/SilentStructStore.hh>

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace mpi_refinement {

using namespace core;

class StructAvrgMover : public protocols::moves::Mover
{
public:
  // Constructor
  StructAvrgMover( pose::Pose const &pose,
									 protocols::wum::SilentStructStore const &decoys,
									 bool const minimize );

  ~StructAvrgMover();

	virtual std::string get_name() const{ return "StructAvrgMover"; }
	virtual void apply( pose::Pose &ref_pose );

  void
  set_default();

	void set_mulfactor( core::Real const mulfactor ){ mulfactor_ = mulfactor; }

	void report_dev( pose::Pose const&ref_pose ) const;

private:

  utility::vector1< Real >
  calculate_variations( utility::vector1< utility::vector1< Real > > const deviation );

  void
  add_deviations( pose::Pose ref_pose, 
		  pose::Pose pose,
		  utility::vector1< utility::vector1< Real > > &deviation
		  );

	utility::vector1< std::pair< Size, Size > >
	predict_region( utility::vector1< core::Real > const CAvar, 
									utility::vector1< bool > &is_region,
									core::Real const min_fluc = 0.2,
									core::Real const frac_base = 0.4
									) const;

  pose::Pose 
  weighted_average( utility::vector1< pose::Pose > &poses,
		    scoring::ScoreFunctionCOP sfxn,
		    pose::Pose const &pose_ref,
				utility::vector1< Real > const, //CAvar,
		    bool const weighted = true
		    );

  void
  shave_poses( utility::vector1< pose::Pose > &poses,
	       pose::Pose const &avrg_pose,
	       Real const frac );

private:
	utility::vector1< pose::Pose > poses_;
	utility::vector1< core::Real > CAvar_;

	Real mulfactor_;
	Real kT_; // denominator for score-weighted averaging
	Real shave_frac_;
	bool minimize_;
	core::Size niter_; // Num iter for shaving outlier poses
	scoring::ScoreFunctionOP sfxn_;

};

} //namespace mpi_refinement
} //namespace protocols

#endif
