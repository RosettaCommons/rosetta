// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_DeleteMover.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_monte_carlo_RNA_DeleteMover_hh
#define INCLUDED_protocols_swa_monte_carlo_RNA_DeleteMover_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/swa/monte_carlo/types.hh>
#include <protocols/swa/monte_carlo/RNA_DeleteMover.fwd.hh>


namespace protocols {
namespace swa {
namespace monte_carlo {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class RNA_DeleteMover: public protocols::moves::Mover {
public:


	//destructor -- necessary? -- YES destructors are necessary.
	~RNA_DeleteMover();

	using protocols::moves::Mover::apply;
  void
	apply( core::pose::Pose & pose, Size const res_to_delete, protocols::swa::monte_carlo::MovingResidueCase const moving_residue_case  );

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	void
	wipe_out_moving_residues( core::pose::Pose & pose );

private:

  void
	remove_cutpoint_variants_at_res_to_delete( core::pose::Pose & pose, Size const & res_to_delete );

};

} // monte_carlo
} // swa
} // protocols

#endif
