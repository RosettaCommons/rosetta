// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/RestrictInteractionGraphThreadsOperation
/// @brief A task operation that restricts the number of threads allowed for interaction graph computation.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_core_pack_task_operation_RestrictInteractionGraphThreadsOperation_hh
#define INCLUDED_core_pack_task_operation_RestrictInteractionGraphThreadsOperation_hh

#include <core/pack/task/operation/RestrictInteractionGraphThreadsOperation.fwd.hh>

// Core headers
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/types.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

///@brief A task operation that restricts the number of threads allowed for interaction graph computation.
class RestrictInteractionGraphThreadsOperation: public core::pack::task::operation::TaskOperation {

public:

	/// @brief Default constructor.
	RestrictInteractionGraphThreadsOperation();

	/// @brief Initialization constructor.
	/// @param[in] thread_limit_in The number of threads to limit the packer to using. If this is set to zero,
	/// this operation does nothing.  Otherwise, it *reduces* the allowed number of threads to the specified value.
	/// If the allowed number of threads is less than the specified value, it does nothing.
	RestrictInteractionGraphThreadsOperation( core::Size const thread_limit_in );

	/// @brief Default copy constructor.
	RestrictInteractionGraphThreadsOperation(RestrictInteractionGraphThreadsOperation const & src);

	/// @brief Destructor.
	~RestrictInteractionGraphThreadsOperation() override;

	/// @brief Clone operation: make a copy of this object, and return a smart pointer to the copy.
	core::pack::task::operation::TaskOperationOP clone() const override;

public:

	/// @brief Configure from a RosettaScripts XML tag.
	void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & ) override;

	/// @brief Alter a PackerTask by reducing the number of allowed threads for packer setup to the number
	/// specified by this TaskOperation.
	/// @details Does nothing if this TaskOperation's allowed threads are set to zero or if the number already
	/// allowed in the PackerTask is less than the number allowed by the TaskOperation.
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const override;

	/// @brief Return the name used to construct this TaskOperation from an XML file.
	static std::string keyname();

	/// @brief Describe the format of XML file used to initialize this TaskOperation.
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Set the maximum number of threads that this TaskOperation will allow for packer setup.
	/// @details Set this to zero to indicate no limit.
	void set_thread_limit( core::Size const setting );

	/// @brief Get the maximum number of threads that this TaskOperation will allow for packer setup.
	/// @details Zero indicates no limit.
	inline core::Size thread_limit() const { return num_threads_; }

private:

	/// @brief The number of threads to limit the packer to using.  If this is set to zero,
	/// this operation does nothing.  Otherwise, it *reduces* the allowed number of threads
	/// to the specified value.  If the allowed number of threads is less than the specified
	/// value, it does nothing.
	core::Size num_threads_ = 0;


};

} //operation
} //task
} //pack
} //core

#endif //INCLUDED_core/pack/task/operation_RestrictInteractionGraphThreadsOperation_fwd_hh
