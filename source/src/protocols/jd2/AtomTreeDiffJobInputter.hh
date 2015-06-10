// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/AtomTreeDiffJobInputter.hh
/// @brief  header file for AtomTreeDiffJobInputter class
/// @author Gordon Lemmon (gordon.h.lemmon@vanderbilt.edu); Rocco Moretti (rmoretti@u.washington.edu)


#ifndef INCLUDED_protocols_jd2_AtomTreeDiffJobInputter_hh
#define INCLUDED_protocols_jd2_AtomTreeDiffJobInputter_hh

//unit headers
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/AtomTreeDiffJobInputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobsContainer.hh>

//project headers
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

#include <utility/vector1.hh>


//utility headers

namespace protocols {
namespace jd2 {

/// @details This is the simplest implementation of JobInputter, which reads from -s/-l and AtomTree files.
class AtomTreeDiffJobInputter : public protocols::jd2::JobInputter
{
public:

	AtomTreeDiffJobInputter();

	virtual ~AtomTreeDiffJobInputter();

	/// @brief this function is responsible for filling the pose reference with
	/// the pose indicated by the job.  The Job object (within its InnerJob)
	/// contains a PoseCOP.  This function needs to either fill the pose
	/// reference from the InnerJob or, on first demand of a pose from that
	/// InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill
	/// the reference.  This implementation uses pose_from_pdb
	virtual void pose_from_job( core::pose::Pose & pose, JobOP job );

	/// @brief this function determines what jobs exist from -in::file::silent and
	/// -in::file::tags.
	virtual void fill_jobs( JobsContainer & jobs );

	/// @brief Return the type of input source that the AtomTreeDiffJobInputter is currently
	///  using.
	/// @return Always <em>ATOM_TREE_FILE</em>.
	virtual JobInputterInputSource::Enum input_source() const;

	/// @brief Utility function to allow mutatable (!) access to the reference
	/// poses of the undelying AtomTreeDiff object. Use with caution.
	///
	/// @details This exists to allow setup on stored reference poses for properties that
	/// don't get saved/restored in PDB format, like covalent constraints for enzyme design.
	core::pose::PoseOPs const & all_ref_poses() const;

private:
	core::import_pose::atom_tree_diffs::AtomTreeDiff atom_tree_diff_;

}; // AtomTreeDiffJobInputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_AtomTreeDiffJobInputter_HH
