// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaProtocol_v1.hh
/// @brief Version 1 model for trRosetta.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_trRosetta_trRosettaProtocol_v1_hh
#define INCLUDED_protocols_trRosetta_trRosettaProtocol_v1_hh

#include <protocols/trRosetta/trRosettaProtocol_v1.fwd.hh>
#include <protocols/trRosetta/trRosettaProtocolBase.hh>
#include <basic/tensorflow_manager/RosettaTensorflowProtocolBase.hh>
#include <basic/citation_manager/CitationCollectionBase.fwd.hh>

// Basic headers
#ifdef USE_TENSORFLOW
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.fwd.hh>
#endif //USE_TENSORFLOW

// Core headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

// C++ headers
#include <utility>

namespace protocols {
namespace trRosetta {

/// @brief Version 1 model for trRosetta.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class trRosettaProtocol_v1 : public protocols::trRosetta::trRosettaProtocolBase {

public:

	/// @brief Constructor.
	trRosettaProtocol_v1();

	/// @brief Copy constructor.
	trRosettaProtocol_v1( trRosettaProtocol_v1 const & ) = default;

	/// @brief Destructor.
	~trRosettaProtocol_v1() override;

	/// @brief Copy this ojbect and return an owning pointer to the copy.
	basic::tensorflow_manager::RosettaTensorflowProtocolBaseOP clone() const override;

public:

#ifdef USE_TENSORFLOW

	/// @brief Run the protocol and generate outputs.
	/// @returns Outputs are for the model version 1: distance, phi, theta,
	/// and omega probability distributions for every pair of amino acids.
	/// @details The set_input_msa_file() method must be called first.
	trRosettaOutputsBaseCOP run() const override;

#endif

	/// @brief Get the name of this RosettaTensorflowProtocol.
	/// @details Must be implemented by derived class.
	std::string name() const override;

	/// @brief Get the citation for trRosetta.
	/// @details ONLY provides the citation for the Yang et al. paper.  Includes
	/// no unpublished module info.
	static
	basic::citation_manager::CitationCollectionBaseCOP
	get_trRosetta_neural_net_citation();

protected:

#ifdef USE_TENSORFLOW
	/// @brief Initialize the Tensorflow model used by this class.
	/// @details Called ONCE by class constructor.  Triggers read from disk!
	void init_tensorflow_session() override;
#endif //USE_TENSORFLOW

private:



};

} //trRosetta
} //protocols

#endif //INCLUDED_protocols_trRosetta_trRosettaProtocol_v1_hh
