// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaProtocolBase.hh
/// @brief A pure virtual base class for trRosetta protocols, derived from RosettaTensorflowProtocolBase.
/// @details Subclasses will be for particular trRosetta versions.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_trRosetta_trRosettaProtocolBase_hh
#define INCLUDED_protocols_trRosetta_trRosettaProtocolBase_hh

#include <protocols/trRosetta/trRosettaProtocolBase.fwd.hh>
#include <basic/tensorflow_manager/RosettaTensorflowProtocolBase.hh>

#ifdef USE_TENSORFLOW
#include <protocols/trRosetta/trRosettaMultipleSequenceAlignment.fwd.hh>
#include <protocols/trRosetta/trRosettaOutputsBase.fwd.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.fwd.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.fwd.hh>
#include <utility/vector1.hh>
#endif //USE_TENSORFLOW


namespace protocols {
namespace trRosetta {

/// @brief A pure virtual base class for trRosetta protocols, derived from RosettaTensorflowProtocolBase.
/// @details Subclasses will be for particular trRosetta versions.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class trRosettaProtocolBase : public basic::tensorflow_manager::RosettaTensorflowProtocolBase {

public:

	/// @brief Default constructor.
	trRosettaProtocolBase() = default;

	/// @brief Copy constructor.
	trRosettaProtocolBase( trRosettaProtocolBase const & ) = default;

	/// @brief Destructor.
	~trRosettaProtocolBase() override;

public:

#ifdef USE_TENSORFLOW
	/// @brief Set the input multiple sequence alignment file.
	/// @details Triggers read from disk!  At the end of this operation, input_msa_file_contents_
	/// is populated with the file's contents.  Throws if read has already occurred, if file can't
	/// be found, or if file is empty.
	/// @param[in] input_msa_file The file to read.
	void set_input_msa_file( std::string const & input_msa_file );

	/// @brief Run the protocol and generate outputs.
	/// @returns Whatever the outputs are for the current model version.
	/// @details Must be implemented by derived classes.
	virtual trRosettaOutputsBaseCOP run() const = 0;

#endif


protected:

#ifdef USE_TENSORFLOW
	/// @brief Initialize the Tensorflow session used by this class.
	/// @details Implemented by derived classes.  Triggers read from disk!
	virtual void init_tensorflow_session() = 0;

	/// @brief Allow derived classes to access the Tensorflow session.
	inline basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP tensorflow_session() const { return tensorflow_session_; }

	/// @brief Allow derived classes to set the Tensorflow session.
	/// @details Should only be called from derived class constructors!
	void set_tensorflow_session( basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP session_in );

	/// @brief Construct the inputs into the trRosetta model, as Tensorflow tensors.
	/// @details The set_input_msa_file() function must be called first!
	utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< int32_t > > construct_input_tensors() const;
#endif //USE_TENSORFLOW

private:

#ifdef USE_TENSORFLOW
	/// @brief The tensorflow session.
	basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP tensorflow_session_;

	/// @brief The contents of the input multiple sequence alignment file.
	trRosettaMultipleSequenceAlignmentOP input_msa_;
#endif //USE_TENSORFLOW

};

} //trRosetta
} //protocols

#endif //INCLUDED_protocols_trRosetta_trRosettaProtocolBase_hh
