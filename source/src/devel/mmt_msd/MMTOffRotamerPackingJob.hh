// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/mmt_msd/MMTOffRotamerPackingJob.hh
/// @brief  declaration for class MMTOffRotamerPackingJob
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_devel_mmt_msd_MMTOffRotamerPackingJob_HH
#define INCLUDED_devel_mmt_msd_MMTOffRotamerPackingJob_HH

// Unit headers
#include <devel/mmt_msd/MMTOffRotamerPackingJob.fwd.hh>
#include <devel/mmt_msd/MMTPackingJob.hh>

// Core headers
#include <core/pack/interaction_graph/SimpleInteractionGraph.fwd.hh>
#include <core/pack/rotamer_set/ContinuousRotamerSet.fwd.hh>
#include <core/pack/scmin/AtomTreeCollection.fwd.hh>
#include <core/pack/scmin/SidechainStateAssignment.fwd.hh>


namespace devel {
namespace mmt_msd {

class MMTOffRotamerPackingJob : public MMTPackingJob
{
public:
	typedef core::pack::scmin::SidechainStateAssignment SidechainStateAssignment;
	typedef core::pack::scmin::SidechainStateAssignmentOP SidechainStateAssignmentOP;

public:
	MMTOffRotamerPackingJob();
	virtual ~MMTOffRotamerPackingJob();

	virtual void setup();
	virtual void optimize();
	virtual void update_pose( core::pose::Pose & pose );

	bool best_assignment_exists() const;

	SidechainStateAssignment const &
	get_best_assignment() const;

	virtual core::Real final_energy() const;

protected:
	virtual void clean_up();

private:

	core::pack::rotamer_set::ContinuousRotamerSetsOP        rotsets_;
	core::pack::scmin::AtomTreeCollectionOP                 atc_;
	core::pack::interaction_graph::SimpleInteractionGraphOP ig_;
	core::pack::scmin::SidechainStateAssignmentOP           best_assignment_;

};

}
}

#endif
