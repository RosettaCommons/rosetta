// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/ImportUnboundRotamersOperation.hh
/// @brief  eliminate aromatic rotamers, of which chi2 are around 0, 180 degree.
/// @author Dave La ( davela@uw.edu )


#ifndef INCLUDED_protocols_toolbox_task_operations_ImportUnboundRotamersOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_ImportUnboundRotamersOperation_hh


// Unit Headers
#include <protocols/toolbox/task_operations/ImportUnboundRotamersOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

#include <core/pack/rotamer_set/UnboundRotamersOperation.fwd.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {


class ImportUnboundRotamersOperation : public core::pack::task::operation::TaskOperation {
public:


	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
	typedef utility::tag::TagCOP TagCOP;


public:


	/// @brief default constructor
	ImportUnboundRotamersOperation();

	/// @brief destructor
	virtual ~ImportUnboundRotamersOperation();

	/// @brief make clone
	virtual TaskOperationOP clone() const;


public:


	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task ) const;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "ImportUnboundRotamers"; }

public:


	void parse_tag( TagCOP tag , DataMap & );

private: // data

	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot_;

};


} // TaskOperations
} // toolbox
} // protocols


#endif
