// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file devel/constrianed_sequence_design/ConstaintManager.hh
/// @brief
/// @author Javier Castellanos	( javiercv@uw.edu )
///

#ifndef INCLUDED_devel_constrained_sequence_ConstraintManager_HH
#define INCLUDED_devel_constrained_sequence_ConstraintManager_HH

// Unit header
#include <devel/constrained_sequence_design/ConstraintManager.fwd.hh>
// Package headers
#include <devel/constrained_sequence_design/SequenceConstraint.hh>
#include <devel/constrained_sequence_design/SequenceConstraintSet.hh>

// types header
#include <core/types.hh>
// Project headers
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// STL headers
#include <set>
#include <map>


namespace devel {
namespace constrained_sequence_design {

class ConstraintManager : public core::pack::task::operation::TaskOperation {
public:
  typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::PoseOP PoseOP;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef std::map< std::string, SequenceConstraintOP > SequenceConstraintMap;
	typedef utility::vector1<Size> VecSize;
public:
		//	---	constructor / destructor ---
			
		// Default constructor (only declared, always should be called with a working pose)
		ConstraintManager();
		
		/// @brief constructor	
		ConstraintManager(Pose& p, SequenceConstraintSetOP constraints);

		//	---	accessors ---
		
		///@brief returns the weighted score of a specific constraint
		Real score(const std::string& constraint_name) ;
		/// @brief returns the raw score of the specific constraint
		Real raw_score(const std::string& constraint_name) ;
		/// @ brief return the sum of the weighed scores for one type of AA in a particular position.
		Real score(Size position, core::chemical::AA aa) const ;
		/// @brief returns the sum of the weighted scores for all the constraints
		Real score() const;

		// --- setters --- 
		/// @brief add a sequence constraint
		void add_constraint(SequenceConstraintOP c);

		// --- inherit methods from Task Operations ---
		/// @brief updates all the positions in the packertask accoring to the constraints
		virtual void apply(const Pose& pose, PackerTask& task) const;
		
		/// @brief clone method
		virtual TaskOperationOP clone() const { return  new ConstraintManager( *this ); }



private:
		// --- private methods ---

		/// @brief this method is called by update() to get the positions in the packer task
		/// than are modifiable. It excludes any position specified in the initial blueprint.
		/// the positions are going to be modified in the same order as returned by pick_positions
		/// so they have to be shuffled before returning them.
		VecSize shuffled_positions() const;

		/// @brief updates the list of included positions to reflect the size of the pose
		/// and the excluded positions. 
		void pick_positions();

		/// @brief intializes the packertask to fit to the pose, loads the information
		/// from a resfile if given and set the included positions.
		/// Sets pose_ as the working for all the constraints in seq_constriants.
		void initialize();

private:
		// --- private data ---
		SequenceConstraintMap seq_constraints_;
		PoseOP pose_;
		std::set<Size> exclude_positions_;	
		VecSize				included_positions_;				
		std::string resfile_;
		bool initialized_;
};

} // constrained_sequence_design
} // devel

#endif
