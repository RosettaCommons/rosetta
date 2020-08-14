// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/tensorflow_manager/RosettaTensorflowProtocolBase.hh
/// @brief A pure virtual base class for storing Tensorflow sessions and the code for executing them and returning a result.
/// Derived classes will be protocol-specific, and will accept a RosettaTensorflowInput and produce a RosettaTensorflowOutput.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
/// @note   This is interface: it has no fields, and only pure virtual methods.  No further constructors should be defined.


#ifndef INCLUDED_basic_tensorflow_manager_RosettaTensorflowProtocolBase_hh
#define INCLUDED_basic_tensorflow_manager_RosettaTensorflowProtocolBase_hh


// Project forward headers
#include <basic/tensorflow_manager/RosettaTensorflowProtocolBase.fwd.hh>

#include <utility/VirtualBase.hh>

#if WIN32
#include <string>
#endif


namespace basic {
namespace tensorflow_manager {


/// @brief A pure virtual base class for storing Tensorflow sessions and the code for executing them and returning a result.
/// Derived classes will be protocol-specific, and will accept a RosettaTensorflowInput and produce a RosettaTensorflowOutput.
class RosettaTensorflowProtocolBase : public utility::VirtualBase {

public: // Creation

	/// @brief Destructor
	~RosettaTensorflowProtocolBase() override;

	/// @brief Clone operation.
	virtual RosettaTensorflowProtocolBaseOP clone() const = 0;

protected: // Creation

	/// @brief Prevent direct instantiation.
	RosettaTensorflowProtocolBase() = default;

	/// @brief Copy constructor.
	RosettaTensorflowProtocolBase( RosettaTensorflowProtocolBase const & ) = default;

public: // Methods
	// Further subsections of methods allowed

	/// @brief Get the name of this RosettaTensorflowProtocol.
	/// @details Must be implemented by derived class.
	virtual std::string name() const = 0;

}; // RosettaTensorflowProtocolBase


} //tensorflow_manager
} //basic


#endif //INCLUDED_basic_tensorflow_manager_RosettaTensorflowProtocolBase_hh



