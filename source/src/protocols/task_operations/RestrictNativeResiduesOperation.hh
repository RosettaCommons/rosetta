// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/RestrictNativeResiduesOperation.hh
/// @brief Restrict every residue in the current pose that is native to repacking. ie, only allow mutated positions to be designed.
/// @author Jacob Bale, balej@u.washington.edu

#ifndef INCLUDED_protocols_task_operations_RestrictNativeResiduesOperation_hh
#define INCLUDED_protocols_task_operations_RestrictNativeResiduesOperation_hh

// unit headers
#include <protocols/task_operations/RestrictNativeResiduesOperation.fwd.hh>

//package headers
#include <core/pack/task/operation/TaskOperation.hh>

//project headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace protocols {
namespace task_operations {

class RestrictNativeResiduesOperation : public core::pack::task::operation::TaskOperation {
public:

	typedef std::string String;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
	typedef utility::tag::TagCOP TagCOP;

public:

	/// @brief default constructor
	RestrictNativeResiduesOperation();

	/// @brief destructor
	~RestrictNativeResiduesOperation();

	/// @brief make clone
	virtual TaskOperationOP clone() const;

public:

	void parse_tag( TagCOP tag , DataMap & );

	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task ) const;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "RestrictNativeResidues"; }

	core::pose::PoseCOP reference_pose() const;
	void reference_pose( core::pose::PoseCOP reference_pose );
	void reference_pose( core::pose::Pose const & pose );
	bool verbose() const;
	void verbose( bool const verb );
	bool prevent_repacking() const;
	void prevent_repacking( bool const prev );

private:
	core::pose::PoseCOP reference_pose_;
	bool verbose_,prevent_repacking_;
	bool invert_;
};


} // TaskOperations
} // protocols
#endif
