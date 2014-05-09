// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  Restrict design to residues matching user-specified SASA criteria in the monomeric, bound, or unbound state.
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_SelectByDeltaScoreOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_SelectByDeltaScoreOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/SelectByDeltaScoreOperation.fwd.hh>
#include <protocols/toolbox/task_operations/SelectByDeltaScoreOperationCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility Headers
#include <string>

// C++ Headers

namespace protocols { 
namespace toolbox {
namespace task_operations {

class SelectByDeltaScoreOperation : public core::pack::task::operation::TaskOperation {
public:

	/// @brief default constructor
	SelectByDeltaScoreOperation();
	/// @brief destructor
	virtual ~SelectByDeltaScoreOperation();
	/// @brief make clone
	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual void apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;
	void parse_tag( TagCOP tag , DataMap & );
	
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	core::scoring::ScoreType score_type() const;
	void score_type( core::scoring::ScoreType const st );
	std::string score_type_name() const;
	void score_type_name( std::string const name );
	core::Real threshold() const;
	void threshold( core::Real const threshold );
	bool lower() const;
	void lower( bool const lower );
	core::pose::PoseOP reference_pose() const;
	void reference_pose( core::pose::PoseOP const reference_pose );
	bool individual_hbonds() const;
	void individual_hbonds( bool const individual_hbonds );

private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreType score_type_; // scoretype of interest, defaults to total_energy 
	std::string score_type_name_;
 	core::Real threshold_; 
 	bool lower_, individual_hbonds_; 
	core::pose::PoseOP reference_pose_;
};

} //namespace task_operations
} //namespace toolbox
} //namespace protocols

#endif 
