// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaProtocolBase.cc
/// @brief A pure virtual class, derived from RosettaTensorflowProtocolBase, which will serve as a base for protocols
/// for predicting peptide fold propensity given peptides of different lengths and different Tensorflow models.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <protocols/trRosetta/trRosettaProtocolBase.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:

#ifdef USE_TENSORFLOW
#include <protocols/trRosetta/trRosettaMultipleSequenceAlignment.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#endif

static basic::Tracer TR( "protocols.trRosetta.trRosettaProtocolBase" );


namespace protocols {
namespace trRosetta {

/// @brief Destructor.
trRosettaProtocolBase::~trRosettaProtocolBase() = default;

////////////////////////////////
// PUBLIC MEMBER FUNCTIONS
////////////////////////////////

#ifdef USE_TENSORFLOW
/// @brief Set the input multiple sequence alignment file.
/// @details Triggers read from disk!  At the end of this operation, input_msa_file_contents_
/// is populated with the file's contents.  Throws if read has already occurred, if file can't
/// be found, or if file is empty.
/// @param[in] input_msa_file The file to read.
void
trRosettaProtocolBase::set_input_msa_file(
	std::string const & input_msa_file
) {
	std::string const errmsg( "Error in trRosettaProtocolBase::set_input_msa_file(): " );
	runtime_assert_string_msg( !input_msa_file.empty(), errmsg + "An empty input filename was provided to this function." );
	runtime_assert_string_msg( input_msa_ == nullptr, errmsg + "A multiple sequence alignment file has already been loaded!" );
	runtime_assert_string_msg( utility::file::file_exists( input_msa_file ), errmsg + "The file \"" + input_msa_file + "\" does not exist!" );
	std::string const filecontents( utility::file_contents( input_msa_file ) );
	runtime_assert_string_msg( !filecontents.empty(), errmsg + "The file \"" + input_msa_file + "\" was empty!" );
	input_msa_ = utility::pointer::make_shared< trRosettaMultipleSequenceAlignment >(filecontents);
	debug_assert( input_msa_ != nullptr );
}
#endif

#ifdef USE_TENSORFLOW
/// @brief Allow derived classes to set the Tensorflow session.
/// @details Should only be called from derived class constructors!
void
trRosettaProtocolBase::set_tensorflow_session(
	basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP session_in
) {
	debug_assert( session_in != nullptr );
	tensorflow_session_ = session_in;
}

/// @brief Construct the inputs into the trRosetta model, as Tensorflow tensors.
/// @details The set_input_msa_file() function must be called first!
utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< int32_t > >
trRosettaProtocolBase::construct_input_tensors() const {
	runtime_assert_string_msg( input_msa_ != nullptr, "Error in trRosettaProtocolBase::construct_input_tensors(): The set_input_msa_file() function must be called before calling this function!" );
	return input_msa_->construct_input_tensors();
}
#endif //USE_TENSORFLOW

} //trRosetta
} //protocols
