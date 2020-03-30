// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowSessionContainer.cc
/// @brief A container for Rosetta Tensorflow sessions, allowing sessions to be loaded once and stored in the global Tensorflow session manager.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.hh>

// Basic headers:
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// Utility headers:
#include <utility/pointer/memory.hh>
#include <utility/file/file_sys_util.hh>

static basic::Tracer TR( "basic.tensorflow_manager.RosettaTensorflowSessionContainer" );


namespace basic {
namespace tensorflow_manager {

#ifdef USE_TENSORFLOW

/// @brief Initialization constructor.
/// @details Takes name of stored session file.  Triggers read from disk.  Only accessible to TensorflowManager.
/// @note If session options are passed in, then the RosettaTensorflowSessionContainer takes control of them and manages their
/// lifetime and destruction.  They must be allocated on the heap with `new`!
RosettaTensorflowSessionContainer::RosettaTensorflowSessionContainer(
	std::string const & filename,
	std::string const & tag,
	TF_SessionOptions* sess_options /*=nullptr*/
) :
	utility::VirtualBase(),
	status_(TF_NewStatus()),
	graph_(TF_NewGraph()),
	sess_options_( sess_options == nullptr ? TF_NewSessionOptions() : sess_options ),
	has_custom_options_( sess_options != nullptr ),
	session_(nullptr) //Initialized below
{
	static const std::string errmsg( "Error in constructor for RosettaTensorflowSessionContainer: " );
	runtime_assert_string_msg( !filename.empty(), errmsg + "The filename for the session cannot be empty." );
	std::string full_filename( filename );
	if ( !utility::file::file_exists(full_filename) ) {
		full_filename = basic::database::full_name( full_filename );
	}
	runtime_assert_string_msg( utility::file::file_exists( full_filename ), errmsg + "Could not find file \"" + filename + "\"." );
	const char *tagchars[] = { &tag.c_str()[0] };
	session_ = TF_LoadSessionFromSavedModel( sess_options_, nullptr, full_filename.c_str(), tagchars, 1, graph_, nullptr, status_ );
	runtime_assert_string_msg( session_ != nullptr, errmsg + "Failed to initialize Tensorflow session from file \"" + full_filename + "\" with tag \"" + tag + "\"." );
}

/// @brief Destructor.
/// @details It is necessary to clean up Tensorflow objects here.
RosettaTensorflowSessionContainer::~RosettaTensorflowSessionContainer(){
	TF_DeleteSession(session_, status_);
	TF_DeleteSessionOptions( sess_options_ );
	TF_DeleteGraph( graph_ );
	TF_DeleteStatus( status_ );
}

/// @brief List all of the operations available in a loaded model.
void
RosettaTensorflowSessionContainer::list_all_operations(
	std::ostream & outstream
) const {
	debug_assert( graph_ != nullptr );
	size_t pos = 0;
	TF_Operation* oper;
	outstream << "All Operations: " << std::endl;
	while ((oper = TF_GraphNextOperation(graph_, &pos)) != nullptr) {
		outstream << TF_OperationName( oper ) << std::endl;
	}
}

#endif //USE_TENSORFLOW

} //tensorflow_manager
} //basic
