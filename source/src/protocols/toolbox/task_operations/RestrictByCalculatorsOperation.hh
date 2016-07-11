// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh
/// @brief  A class that applies arbitrary calculators (whose calculations return std::set< core::Size >) to restrict a PackerTask
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictByCalculatorsOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictByCalculatorsOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Utility Headers

// C++ Headers

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

/// @details this class is a TaskOperation to prevent repacking of residues according to std::set< core::Size >-returning PoseMetricCalculators
class RestrictByCalculatorsOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;
	typedef std::pair< std::string, std::string> calc_calcn; //calculator and calculation

	RestrictByCalculatorsOperation();
	RestrictByCalculatorsOperation( utility::vector1< calc_calcn > const & calcs_and_calcns );

	virtual ~RestrictByCalculatorsOperation();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "RestrictByCalculators"; }

private:

	utility::vector1< calc_calcn > const calcs_and_calcns_;

};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictByCalculatorsOperation_HH
