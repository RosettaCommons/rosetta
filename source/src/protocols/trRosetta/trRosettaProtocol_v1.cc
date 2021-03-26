// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaProtocol_v1.cc
/// @brief A pure virtual class, derived from RosettaTensorflowProtocolBase, which will serve as a base for protocols
/// for predicting peptide fold propensity given peptides of different lengths and different Tensorflow models.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <protocols/trRosetta/trRosettaProtocol_v1.hh>
#include <protocols/trRosetta/trRosettaOutputs_v1.hh>

// Basic headers:
#include <basic/Tracer.hh>

#ifdef USE_TENSORFLOW
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.tmpl.hh>
#include <basic/tensorflow_manager/RosettaTensorflowManager.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.tmpl.hh>
#else //Not USE_TENSORFLOW
#include <basic/tensorflow_manager/util.hh>
#endif //USE_TENSORFLOW

#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/CitationManager.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.trRosetta.trRosettaProtocol_v1" );

namespace protocols {
namespace trRosetta {

/// @brief Constructor.
trRosettaProtocol_v1::trRosettaProtocol_v1() :
	trRosettaProtocolBase()
{
#ifdef USE_TENSORFLOW
	init_tensorflow_session();
#else //Not USE_TENSORFLOW
	utility_exit_with_message( basic::tensorflow_manager::get_tensorflow_compilation_instructions( "trRosettaProtocol_v1 protocol" ) );
#endif //USE_TENSORFLOW
}

/// @brief Destructor.
trRosettaProtocol_v1::~trRosettaProtocol_v1() = default;

/// @brief Copy this ojbect and return an owning pointer to the copy.
basic::tensorflow_manager::RosettaTensorflowProtocolBaseOP
trRosettaProtocol_v1::clone() const {
	return utility::pointer::make_shared< trRosettaProtocol_v1 >( *this );
}

////////////////////////////////
// PUBLIC MEMBER FUNCTIONS:
////////////////////////////////


#ifdef USE_TENSORFLOW
/// @brief Run the protocol and generate outputs.
/// @returns Outputs are for the model version 1: distance, phi, theta,
/// and omega probability distributions for every pair of amino acids.
/// @details The set_input_msa_file() method must be called first.
trRosettaOutputsBaseCOP
trRosettaProtocol_v1::run() const {
	using namespace basic::tensorflow_manager;

	// Deliberately using double and float and int32_t below here.  Do not change to core::Real!
	std::chrono::duration< double, std::micro > runtime(0.0);
	utility::vector1< RosettaTensorflowTensorContainer< int32_t > > const input_tensors( construct_input_tensors() );
	utility::vector1< RosettaTensorflowTensorContainer< float > > output_tensors{
		RosettaTensorflowTensorContainer<float>( utility::vector1< int64_t >{ static_cast<int64_t>(input_tensors[2](1)), static_cast<int64_t>(input_tensors[2](1)), static_cast<int64_t>(trRosetta_v1_expected_dist_bins) }, 40.0 ),
		RosettaTensorflowTensorContainer<float>( utility::vector1< int64_t >{ static_cast<int64_t>(input_tensors[2](1)), static_cast<int64_t>(input_tensors[2](1)), static_cast<int64_t>(trRosetta_v1_expected_omega_bins) }, 40.0 ),
		RosettaTensorflowTensorContainer<float>( utility::vector1< int64_t >{ static_cast<int64_t>(input_tensors[2](1)), static_cast<int64_t>(input_tensors[2](1)), static_cast<int64_t>(trRosetta_v1_expected_phi_bins) }, 40.0 ),
		RosettaTensorflowTensorContainer<float>( utility::vector1< int64_t >{ static_cast<int64_t>(input_tensors[2](1)), static_cast<int64_t>(input_tensors[2](1)), static_cast<int64_t>(trRosetta_v1_expected_theta_bins) }, 40.0 )
	};

	RosettaTensorflowSessionContainerCOP sess( tensorflow_session() );

	TR << "Running version 1 trRosetta neural network..." << std::endl;

	sess->run_multiinput_multioutput_session(
		utility::vector1< std::string >{ "input/XXX_nrow", "input/XXX_ncol", "input/XXX_msa" },
		utility::vector1< std::string >{ "YYY_prob_dist", "YYY_prob_omega", "YYY_prob_phi", "YYY_prob_theta" },
		input_tensors,
		output_tensors,
		runtime
	);

	TR << "Ran version 1 trRosetta protocol in " << static_cast<core::Real>(runtime.count())/1.0e3 << " milliseconds." << std::endl;

	return utility::pointer::make_shared< trRosettaOutputs_v1 >( input_tensors[2](1) /*The sequence length*/, output_tensors /*distance, phi, theta, omega*/ );
}
#endif

/// @brief Get the name of this RosettaTensorflowProtocol.
/// @details Must be implemented by derived class.
std::string
trRosettaProtocol_v1::name() const {
	return "trRosettaProtocol_v1";
}

/// @brief Get the citation for trRosetta.
/// @details ONLY provides the citation for the Yang et al. paper.  Includes
/// no unpublished module info.
/*static*/
basic::citation_manager::CitationCollectionBaseCOP
trRosettaProtocol_v1::get_trRosetta_neural_net_citation() {
	using namespace basic::citation_manager;
	CitationCollectionOP citation(
		utility::pointer::make_shared< CitationCollection >(
		"trRosetta", CitedModuleType::NeuralNetwork
		)
	);
	citation->add_citation( CitationManager::get_instance()->get_citation_by_doi( "10.1073/pnas.1914677117" ) );
	return citation;
}


////////////////////////////////
// PROTECTED MEMBER FUNCTIONS:
////////////////////////////////

#ifdef USE_TENSORFLOW
/// @brief Initialize the Tensorflow model used by this class.
/// @details Called ONCE by class constructor.  Triggers read from disk!
void
trRosettaProtocol_v1::init_tensorflow_session() {
	using namespace basic::tensorflow_manager;
	debug_assert( tensorflow_session() == nullptr );
	set_tensorflow_session( RosettaTensorflowManager::get_instance()->get_session( "protocol_data/tensorflow_graphs/tensorflow_graph_repo_submodule/trRosetta/model_v1/", "serve" ) );
}
#endif //USE_TENSORFLOW

} //trRosetta
} //protocols
