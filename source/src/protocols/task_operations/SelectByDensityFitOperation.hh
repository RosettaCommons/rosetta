// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Restrict design to residues with electron density correlation above threshold value
/// @author Patrick Conway

#ifndef INCLUDED_protocols_task_operations_SelectByDensityFitOperation_hh
#define INCLUDED_protocols_task_operations_SelectByDensityFitOperation_hh

// Unit Headers
#include <protocols/task_operations/SelectByDensityFitOperation.fwd.hh>
#include <protocols/task_operations/SelectByDensityFitOperationCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Utility Headers
#include <string>

// C++ Headers

namespace protocols {
namespace task_operations {

class SelectByDensityFitOperation : public core::pack::task::operation::TaskOperation {
public:
	SelectByDensityFitOperation( core::Real threshold = 0.72, bool invert = 0 );

	virtual ~SelectByDensityFitOperation();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	void parse_tag( TagCOP tag, DataMap & );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "SelectByDensityFit"; }

private:
	core::Real threshold_;
	bool invert_;
};

} //namespace task_operations
} //namespace protocols

#endif
