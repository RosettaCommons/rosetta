// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file devel/constrianed_sequence_design/ConstaintManager.cc
/// @brief
/// @author Javier Castellanos	(javiercv@uw.edu)

// unit header
#include <devel/constrained_sequence_design/ConstraintManager.hh>

// package headers
#include <devel/constrained_sequence_design/SequenceConstraint.hh>
#include <devel/constrained_sequence_design/SequenceConstraintSet.hh>

// types header
#include <core/types.hh>

// project headers
#include <core/pose/Pose.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>

// STL headers
#include <algorithm>

static basic::Tracer TR("devel.constrained_sequence_design.ConstraintManager");

namespace devel {
namespace constrained_sequence_design {

using namespace core;
using namespace utility;

/// @brief default constructor
ConstraintManager::ConstraintManager( Pose& p, SequenceConstraintSetOP constraints ):
    pose_(p)
{
	pick_positions();
	for(SequenceConstraintSet::iterator it = constraints->begin(); it != constraints->end(); ++it) {
		add_constraint(*it);
	}
}

void
ConstraintManager::add_constraint( SequenceConstraintOP c ) {
  seq_constraints_[ c->name() ] =  c;
}

Real
ConstraintManager::score( const std::string& constraint_name ) {
    return seq_constraints_[constraint_name]->score();
}

Real
ConstraintManager::raw_score(const std::string& constraint_name) {
    return seq_constraints_[constraint_name]->raw_score();
}

/// @brief returns the sum of the weighted scores for all the constraints
Real
ConstraintManager::score() const 
{
    Real sum = 0.0;
    for(SequenceConstraintMap::const_iterator it = seq_constraints_.begin(); it != seq_constraints_.end(); ++it) {
        sum += (it->second)->score();
    }
    return sum;
}


/// @ brief return the sum of the weighed scores for one type of AA in a particular position.
Real
ConstraintManager::score(Size position, core::chemical::AA aa) const 
{
    Real sum = 0.0;
    for(SequenceConstraintMap::const_iterator it = seq_constraints_.begin(); it != seq_constraints_.end(); ++it) {
        sum += (it->second)->apply(aa, position);
    }
		return sum;
}

// @brief updates all the positions in the packertask that have not been specified in
// the initial resfile trying to satisfy as many constraints as possible.
void 
ConstraintManager::apply(const Pose& pose, PackerTask& task)  const
{
	using core::pack::task::ResidueLevelTask;
	typedef core::pack::task::ResidueLevelTask::ResidueTypeCAPListConstIter ResidueTypeCAPListConstIter;

	// Update the constraints to reflect the state of the pose.
  for(SequenceConstraintMap::const_iterator cntr = seq_constraints_.begin(); cntr != seq_constraints_.end(); ++cntr) 
		( cntr->second )->update( pose );
	// pick the positions of the packertask that need to be updated
  VecSize positions = shuffled_positions();

	// iterate over the positions of the packer task
	for(VecSize::const_iterator p = positions.begin(); p != positions.end(); ++p) {
		// and for each position iterate over the residue types 
		ResidueTypeCAPListConstIter restype_begin  = task.nonconst_residue_task(*p).allowed_residue_types_begin();
		ResidueTypeCAPListConstIter restype_end    = task.nonconst_residue_task(*p).allowed_residue_types_end();
		for(ResidueTypeCAPListConstIter it = restype_begin; it != restype_end; ++it) {
			core::chemical::AA  aa = (*it)->aa();
			Real score_p = score(*p, aa);
			if( score_p >= 0.0)
				TR << "ALLOW " << (*it)->name() << " at position " << *p << std::endl; // Add residue to the task
			else
				TR << "DISALLOW " << (*it)->name() << " at position " << *p << std::endl; // remove residue from the task	
		}
	}
} // apply

utility::vector1<Size>
ConstraintManager::shuffled_positions() const {
				utility::vector1<Size>  pos(included_positions_);
		std::random_shuffle(pos.begin(), pos.end());
    return pos; 
} // shuffled_positions

/// @brief updates the list of included positions to reflect the size of the pose
/// and the excluded positions. 
void 
ConstraintManager::pick_positions()
{
	included_positions_.clear();
	for(Size i = 1; i <= pose_->total_residue(); ++i)
		if( ! exclude_positions_.count(i) )
			included_positions_.push_back(i);
} // pick_positions


} // constrained_sequence_design
} // devel

