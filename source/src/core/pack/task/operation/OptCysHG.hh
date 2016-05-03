// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/OptCysHG.hh
/// @brief  run optH on non-disulfided bonded CYS only; meant to relieve
///         any clashes caused by swapping of CYD->CYS after calling
///         Conformation::detect_disulfides()
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pack_task_operation_OptCysHG_hh
#define INCLUDED_core_pack_task_operation_OptCysHG_hh

// unit headers
#include <core/pack/task/operation/OptCysHG.fwd.hh>

// project headers
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.fwd.hh>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {


/// @brief run optH on non-disulfided bonded CYS only; meant to relieve any
///  clashes caused by swapping of CYD->CYS after calling Conformation::detect_disulfides()
class OptCysHG : public core::pack::task::operation::TaskOperation {

private: // typedefs

	typedef core::pack::task::operation::TaskOperation Super;

public: // typedefs

	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pose::Pose Pose;

public: // construct/destruct

	/// @brief default constructor
	OptCysHG();

	/// @brief copy constructor
	OptCysHG( OptCysHG const & rval );

	/// @brief default destructor
	virtual ~OptCysHG();

public: // virtual constructors

	/// @brief clone this object
	virtual TaskOperationOP clone() const;

public: // methods

	/// @brief apply operations to PackerTask
	virtual void apply( Pose const & pose, PackerTask & task ) const;

	static std::string keyname();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core


#endif /* INCLUDED_core_pack_task_operation_OptCysHG_HH */
