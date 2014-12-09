// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/task_operations/RestrictToMoveMapChiOperation.hh
/// @brief Makes it easier to use TF and MM in the same class
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictToMoveMapChiOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictToMoveMapChiOperation_hh

#include <protocols/toolbox/task_operations/RestrictToMoveMapChiOperationCreator.hh>
#include <protocols/toolbox/task_operations/RestrictToMoveMapChiOperation.fwd.hh>
#include <basic/datacache/DataMap.hh>

#include <basic/datacache/DataMap.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/tag/Tag.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {
	
	using namespace core::pack::task::operation;
	using core::kinematics::MoveMapCOP;
	using core::pose::Pose;
	using core::pack::task::PackerTask;
	
	
///@brief A TaskOperation that accepts a movemap and restricts chi that are false to either packing or design.
/// Does not turn anything on, just like the rest of the RestrictTo operations.
class RestrictToMoveMapChiOperation : public core::pack::task::operation::TaskOperation {

public:
	
	typedef core::pack::task::operation::TaskOperation parent;
	
public:
	
	RestrictToMoveMapChiOperation();
	
	RestrictToMoveMapChiOperation( MoveMapCOP movemap );

	virtual ~RestrictToMoveMapChiOperation();
	

public:
	
	void
	set_movemap( MoveMapCOP movemap );
	 
	///@brief Set residues from movemap to designable.  Default false.
	void
	set_design( bool setting );

	///@brief Set to use neighbor residues in vacinity of movemap chi residues for packing.  Default False.
	void
	set_include_neighbors( bool setting );
	
	///@brief Cutoff distance for neighbor detection. Default is 10 A.
	void
	set_cutoff_distance( core::Real cutoff );
	
	
public:
	
	//@brief.  Yes, I don't know wtf I'm doing with RosettaScripts yet.
	//virtual void parse_tag( TagCOP, basic::datacache::DataMap & data, core);
	
	virtual TaskOperationOP clone() const;
	
	virtual
	void
	apply( Pose const & pose, PackerTask & task ) const;
	
	
private:
	
	void
	init();
	
	MoveMapCOP movemap_;
	bool design_;
	bool include_neighbors_;
	bool movemap_set_; //Since I can't assign a COP to null it seems.
	core::Real cutoff_;
	
};
	
	
	
	
	
}
}
}
#endif	//INCLUDED_protocols_toolbox_task_operations_RestrictToMoveMapChiOperation_hh


