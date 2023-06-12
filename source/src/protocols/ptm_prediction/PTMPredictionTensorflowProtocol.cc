// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ptm_prediction/PTMPredictionTensorflowProtocol.cc
/// @brief A class for predicting post-translational modifications using neural networks.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adapted from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)

// project headers
#include <protocols/ptm_prediction/PTMPredictionTensorflowProtocol.hh>


// // utility headers
#include <utility/pointer/owning_ptr.hh>
#include <basic/database/open.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowManager.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.tmpl.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.tmpl.hh>

// core headers
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/conformation/Residue.hh>

// protocol headers
#include <protocols/moves/DsspMover.hh>
#include <protocols/ptm_prediction/PTMPredictionMetric.hh>

// STL headers:
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>

static basic::Tracer TR( "protocols.ptm_prediction.PTMPredictionTensorflowProtocol" );

namespace protocols {
namespace ptm_prediction {

/// @brief Default constructor.
PTMPredictionTensorflowProtocol::PTMPredictionTensorflowProtocol() :
	PTMPredictionTensorflowProtocolBase()

{

}

/// @brief constructor with file path
PTMPredictionTensorflowProtocol::PTMPredictionTensorflowProtocol( std::string path ) : path_to_model( path ) {

#ifdef USE_TENSORFLOW
		set_tensorflow_session( get_tensorflow_session() );
#endif // USE_TENSORFLOW
}

/// @brief Copy constructor.  Keep default unless deep copying is needed (and in that case,
/// consider using DeepCopyOPs.)
PTMPredictionTensorflowProtocol::PTMPredictionTensorflowProtocol( PTMPredictionTensorflowProtocol const & )=default;

/// @brief Destructor.
PTMPredictionTensorflowProtocol::~PTMPredictionTensorflowProtocol(){}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
basic::tensorflow_manager::RosettaTensorflowProtocolBaseOP
PTMPredictionTensorflowProtocol::clone() const {
	return utility::pointer::make_shared< PTMPredictionTensorflowProtocol >( *this );
}

/// @brief Get the name of this protocol.
std::string
PTMPredictionTensorflowProtocol::name() const {
	return "PTMPredictionTensorflowProtocol";
}

#ifdef USE_TENSORFLOW
/// @brief Analyze a given Asparagine in a pose, returning an deamidation probability (ranging from 0 to 1).
/// @details Deamidation is the spontaneous reaction of Asn to (Iso-)Asp.
std::map< core::Size, core::Real >
PTMPredictionTensorflowProtocol::compute_deamidation_probability(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const &selection
) const {

	std::chrono::time_point<std::chrono::system_clock> const overall_starttime( basic::tensorflow_manager::ROSETTA_TENSORFLOW_CLOCK::now() );  

	//Selected pose residues:
	utility::vector1< core::Size > const selected_residues( get_selected_residues(selection) );

 	// Defining the name of the first and last layer of our NN 
 	std::string const input_name = "input";
 	std::string const output_name = "output/Sigmoid";
	// Loading in the deamidation rate table
	std::vector< std::vector < float > > deamidation_rates = parseCSV();
	// setting up return_map
	std::map< core::Size, core::Real > return_map; 
	for ( auto const position : selected_residues ) {
      // Creating an tensorflow container for the input & output tensor
    	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > input_tensor;
    	input_tensor.initialize( TF_FLOAT, { 1, 26 } );

    	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > output_tensor;
    	output_tensor.initialize( TF_FLOAT, {1, 2 } );
     	// Calculating the features
      copy_deamidation_features_from_pose_to_tensor( pose, position, input_tensor, deamidation_rates );

    	TR  << "================== Starting prediction ===============" << std::endl;
    	std::chrono::duration< double, std::micro > runtime;
    	// predicting based on the input tensor 
			basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP session = PTMPredictionTensorflowProtocolBase::tensorflow_session();
    	session->run_session( input_name, output_name, input_tensor, output_tensor, runtime );		
    	TR << "Predicted probability at position " << std::fixed << std::setprecision( 1 ) << position << ": " << output_tensor( 1 ) << std::endl;
    	//Copy outputs out of tensors, and delete tensors:
    	core::Real const returnval( static_cast< float >( output_tensor( 1 ) ) );
    	return_map[position] = returnval;

    	std::chrono::time_point<std::chrono::system_clock> const overall_endtime( basic::tensorflow_manager::ROSETTA_TENSORFLOW_CLOCK::now() );
    	std::chrono::duration< double, std::micro > const overall_diff( overall_endtime - overall_starttime );

      TR << "Overall time for preparing Tensorflow inputs, running Tensorflow session, and retrieving outputs: " << overall_diff.count() << " microseconds." << std::endl;
	} 

	return return_map;
}
#endif //USE_TENSORFLOW

#ifdef USE_TENSORFLOW
/// @brief Get the Tensorflow session used by this class.
/// @details Implemented by derived classes.
basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP
PTMPredictionTensorflowProtocol::get_tensorflow_session(
				) const {
	return basic::tensorflow_manager::RosettaTensorflowManager::get_instance()->get_session( path_to_model, "serve" );
}
/// @brief Given a pose and an already-allocated (but empty) input tensor, store the relevant pose data for the
/// deamidation site prediction in the tensor.
/// @details The input tensor is a 1D matrix. Following features are stored:
/// 1: SASA of the selected Asn.
/// 2: Phi angle of the selected Asn.
/// 3: Psi angle of the selected Asn.
/// 4: Chi_1 angle of the selected Asn.
/// 5: Chi_2 angle of the selected Asn.
/// 6: Attack distance between C-Beta of Asn and next residue N atom
/// 7: Deamidation rates from Robinson, N. E. et al. Journal of Peptide Research (2001)
/// 8-26: one-hot encoded amino acid
void
PTMPredictionTensorflowProtocol::copy_deamidation_features_from_pose_to_tensor(
	core::pose::Pose const & pose,
  core::Size const position,
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & input_tensor,
  std::vector< std::vector< float > > & deamidation_rates
) const {

	// get pose sequence
	std::string pose_sequence = pose.sequence();

	// make sure that an Asparagine is selected
	std::string aa_select = pose.residue_type( position ).base_name();
	runtime_assert_string_msg( aa_select == "ASN", "Error in PTMPredictionTensorflowProtocol::copy_deamidation_features_from_pose_to_tensor(): The " + name() + " module requires a single Asparagine, but an " + aa_select + " was selected instead!" );

	using namespace core::scoring::sasa; 
	SasaCalc sasa_calc = SasaCalc();
	sasa_calc.calculate(pose);
	utility::vector1< core::Real > residue_sasa_vec = sasa_calc.get_residue_sasa();
	core::Real residue_sasa = residue_sasa_vec[ position ];

	core::Real phi = pose.phi( position ); 
	core::Real psi = pose.psi( position );
	core::Real chi_1 = pose.chi( 1, position );
	core::Real chi_2 = pose.chi( 2, position );

	// Alphabet for one-hot encoding 
	std::string aaalphabet = "ADEFGHIKLMNPQRSTVWY"; // does not include Cysteine
	// make sure its all canonical amin acids
	runtime_assert_string_msg( pose.residue( position + 1 ).has_property( "CANONICAL_AA" ) && pose.residue( position - 1 ).has_property( "CANONICAL_AA" ), "Residues around modification must be a canonical amino acid!" );	

  // residue following Asparagine (N+1)
	char aa_plus1 = pose_sequence[ position ];
	// residue before Asparagine (N-1)
	char aa_min1 = pose_sequence[ position-2 ];

	core::Real attack_distance = pose.residue( position ).xyz( "CG" ).distance( pose.residue( position+1 ).xyz( "N" ) );
	// Order of the AAs in the rate table:
	std::string rate_alphabet =  "ARDCEGHILKMFPSTWYVNQ";

	// Loading in the deamidation rate table
	// std::vector< std::vector < float > > deamidation_rates = parseCSV();
  // Assigning the different features to the tensor
	input_tensor( 1 ) = residue_sasa;
	input_tensor( 2 ) = phi;
	input_tensor( 3 ) = psi;
	input_tensor( 4 ) = chi_1;
	input_tensor( 5 ) = chi_2;
	input_tensor( 6 ) = attack_distance;
	// half lifes for Q/N weren't reported so we use the average half life time like wid did in training
	if ( aa_min1 == 'N' || aa_min1 == 'Q' || aa_plus1 == 'N' || aa_plus1 == 'Q' ) {
			input_tensor( 7 ) = 5.7;
	} else {
			input_tensor( 7 ) = deamidation_rates[ rate_alphabet.find( aa_min1 ) ][ rate_alphabet.find( aa_plus1 ) ];
	}
	for ( int i = 8; i < 27; i++ ) {
			input_tensor( i ) = 0;
	}
  input_tensor( aaalphabet.find( aa_plus1 ) + 8 ) = 1.0;
}	
#endif //USE_TENSORFLOW
#ifdef USE_TENSORFLOW
/// @brief Analyze a given position in a pose, returning PTM probabilities (ranging from 0 to 1).
/// @details Uses a neural network to predict the probability of a particular modification.
std::map< core::Size, core::Real >
PTMPredictionTensorflowProtocol::compute_ptm_probability(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	PTMPredictionMetric::PTMPredictionMetricModes const ptm_mode
) const {
#ifdef USE_TENSORFLOW_GPU
 	TR << "Starting PTM prediction using Tensorflow-GPU" << std::endl;
#else //USE_TENSORFLOW_GPU
#ifdef USE_TENSORFLOW_CPU
 	TR << "Starting PTM prediction using Tensorflow-CPU" << std::endl;
#else
	utility_exit_with_message( "Either USE_TENSORFLOW_GPU or USE_TENSORFLOW_CPU should be defined! If not you are probably in deep trouble." );
#endif //USE_TENSORFLOW_CPU
#endif //USE_TENSORFLOW_GPU


	std::chrono::time_point<std::chrono::system_clock> const overall_starttime( basic::tensorflow_manager::ROSETTA_TENSORFLOW_CLOCK::now() );  
// Calculating the features
	//Selected pose residues:
	utility::vector1< core::Size > const selected_residues( get_selected_residues(selection) );
  std::string const output_name = "StatefulPartitionedCall";
  // get pose sequence
	std::map< core::Size, core::Real > return_map;
  std::string pose_sequence = pose.sequence();
	for ( auto const position : selected_residues ) {
      // Creating an tensorflow container for the input & output tensor
     	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > input_tensor_seq;
     	input_tensor_seq.initialize( TF_FLOAT, { 1, 9 } );

     	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > input_tensor_struc;
     	input_tensor_struc.initialize( TF_FLOAT, { 1, 10 } );    	

     	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > output_tensor;

		  size_t index;
		  utility::vector1< std::string > input_names;
			std::string aa_select = pose.residue_type( position ).base_name();

		  // setup the shape of the output tensor depending on the model used
			setup_and_check( aa_select, ptm_mode, output_tensor, index, input_names );

			// std::ostream& os = std::cout;
			// session->list_all_operations( os ); // useful for debugging graph signatures

    	copy_feat_to_tensor( pose, position, input_tensor_seq, input_tensor_struc );
			utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > > input_vec {
					input_tensor_seq,
				  input_tensor_struc
			};
    	TR  << "================== Running tensorflow graph  ===============" << std::endl;

      // Chrono timing required by run_session
    	std::chrono::duration< double, std::micro > runtime;
    	// predicting based on the input tensor 
			basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP session = PTMPredictionTensorflowProtocolBase::tensorflow_session();
    	session->run_multiinput_session(
							input_names,
							output_name, 
							input_vec, 
							output_tensor, 
							runtime 
			);
    	TR << "Computed probability in " << runtime.count() << " microseconds." << std::endl;
		  core::Real const returnval( static_cast< float >( output_tensor( index ) ) );
		   // NlinkedGlycosylation was just trained on sequons that can be glycosylated, easy cases were not part of training
			 // thats why we set those to 0 by hand, as they won't get glycosylated anyway
			core::Size const next_pose_index = pose.residue( position ).connected_residue_at_upper();
			core::Size const next_next_pose_index = pose.residue( next_pose_index ).connected_residue_at_upper();
       if ( ptm_mode == PTMPredictionMetric::PTMPredictionMetricModes::N_LINKED_GLYCOSYLATION && next_pose_index != 0 && pose.residue_type( next_pose_index ).aa() == core::chemical::AA::aa_pro ) {
					return_map[ position ] = 0.0;
					TR << "WARNING: Proline after the potential site, extremely unlikely to be glycosylated." << std::endl;
			 }
			 else if ( ptm_mode == PTMPredictionMetric::PTMPredictionMetricModes::N_LINKED_GLYCOSYLATION && next_pose_index != 0 && next_next_pose_index != 0 && pose.residue_type( next_next_pose_index ).aa() != core::chemical::AA::aa_ser  && pose.residue_type( next_next_pose_index ).aa() != core::chemical::AA::aa_thr ) {
					 return_map[ position ] = 0.0;
					 TR << "WARNING: Potential glycosylation site has no sequon, extremely unlikely to be glycosylated." << std::endl;
			 }
			 else {
				  return_map[ position ] = returnval;
			 }
    	std::chrono::time_point<std::chrono::system_clock> const overall_endtime( basic::tensorflow_manager::ROSETTA_TENSORFLOW_CLOCK::now() );
    	std::chrono::duration< double, std::micro > const overall_diff( overall_endtime - overall_starttime );

			TR << "Predicted Probabiltiy at position " << position << ": " << std::fixed << std::setprecision( 1 ) << return_map[ position ] << std::endl;

    	TR << "Overall time for preparing Tensorflow inputs, running Tensorflow session, and retrieving outputs: " << overall_diff.count() << " microseconds." << std::endl;
    	TR << "===========================================================" << std::endl;
	} 
	return return_map;

}
#endif //USE_TENSORFLOW

#ifdef USE_TENSORFLOW
/// @brief Given a pose and an already-allocated (but empty) input tensors, store the relevant pose data for the
/// PTM prediction in the tensors.
/// @details The input tensor is a 1D matrix. Following features are stored:
// input_tensor_seq: Sequence Window of 4 residues before/after potential modification site
// input_tensor_struc: Phi Angles (-1/+1), Psi Angles (-1/+1), sasa, secondary structure E/H/L
void
PTMPredictionTensorflowProtocol::copy_feat_to_tensor(
	core::pose::Pose const & pose,
	core::Size const & position,
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & input_tensor_seq,
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & input_tensor_struc
) const {
	// get pose sequence
	std::string const pose_sequence = pose.sequence();
	// define an AA alphabet, include X placeholder for N/C-termini
	std::string alphabet =  "ACDEFGHIKLMNPQRSTVWY";
	// loop through the 21mer sequence window of the potential modification site
	// and input the amino acid identity in the tensor
	// check that we are at least 4 residues inside N/C-Termini
	runtime_assert_string_msg( 
					position > 4 and position + 4 < pose_sequence.length(),
					"Position must be 4 residues away from N/C Terminus!"
	);
	bool residues_correct_connect;
	for ( int i = -5; i < 4; i++ ) {
			int window_position = position + i + 1; // 1-indexed position
			size_t window_position_usigned = window_position; // to make comparisons to functions returning usigned integers
			// check that the residues are correctly connected to residue before/after
		  if ( window_position_usigned < pose.total_residue() && pose.residue( window_position_usigned ).connected_residue_at_upper() == window_position_usigned+1 && pose.residue( window_position_usigned+1 ).connected_residue_at_lower() == window_position_usigned ) {
					residues_correct_connect = true;
      // check if first residue is N-terminus 
			} else if ( i == -5 && pose.residue( window_position_usigned+1 ).connected_residue_at_lower() == window_position_usigned && pose.residue( window_position_usigned ).is_lower_terminus() == true ) {
					residues_correct_connect = true;
			// check if last residue is C-terminus		 
			} else if ( i == 3 && pose.residue( window_position_usigned ).connected_residue_at_upper() == window_position_usigned+1 && pose.residue( window_position_usigned ).is_upper_terminus() == true ) {
					residues_correct_connect = true;
      // if none of this applies we don't have a correctly connected residue window
			} else {
					residues_correct_connect = false;
			}
			runtime_assert_string_msg( residues_correct_connect == true, "Residues around modification site are either too close to a termini (<4 residues, or are not correctly connected!" );
	    // check that the current residue is canonical amino acid
			runtime_assert_string_msg( pose.residue( window_position ).has_property( "CANONICAL_AA" ), "Residues around modification must be a canonical amino acid!" );	
			// get amino acid identity and input into the tensor
			char aa_name = pose_sequence[ window_position - 1 ]; // -1 since sequence isn't 1 indexed
			input_tensor_seq( i + 6 ) = alphabet.find( aa_name ); 
			if ( i >= -1 and i <= 1 ) {
					// add phi angle
					core::Real phi_angle = pose.phi( window_position );
					input_tensor_struc( i + 2 ) = phi_angle;
					// add psi angle
					core::Real psi_angle = pose.psi( window_position );
					input_tensor_struc( i + 5 ) = psi_angle;
			}
	}

  // Caculating SASA for potential modified site
	using namespace core::scoring::sasa; 
	SasaCalc sasa_calc = SasaCalc();
	sasa_calc.calculate( pose );
	utility::vector1< core::Real > residue_sasa_vec = sasa_calc.get_residue_sasa();
	core::Real residue_sasa = residue_sasa_vec[ position ];
  // Adding SASA to input tensor
	input_tensor_struc( 7 ) = residue_sasa;
  // Checking secondary structure, have to clone pose so its not const
	moves::DsspMover dssp;
	core::pose::Pose dssp_pose = *pose.clone();
	dssp.apply( dssp_pose );
	runtime_assert_string_msg( dssp.get_last_move_status() == moves::MoverStatus::MS_SUCCESS, "DSSP mover failed, check your input pose" );
	char secstruct = dssp_pose.secstruct( position );
  // Adding coded sec. structure to tensor
	if ( secstruct == 'E' ) {
			input_tensor_struc( 8 ) = 1;
			input_tensor_struc( 9 ) = 0;
			input_tensor_struc( 10 ) = 0;
	}
	else if ( secstruct == 'H' ) {
			input_tensor_struc( 8 ) = 0;
			input_tensor_struc( 9 ) = 1;
			input_tensor_struc( 10 ) = 0;
	}
	else if ( secstruct == 'L' ) {
			input_tensor_struc( 8 ) = 0;
			input_tensor_struc( 9 ) = 0;
			input_tensor_struc ( 10 ) = 1;
	}
	else {
			utility_exit_with_message("Secstruct doesn't match E/H/L, sth with the DSSP mover went horribly wrong!");
	}
}

/// @brief Helper function to load the deamidation rates from a csv file
std::vector< std::vector< float > > 
PTMPredictionTensorflowProtocol::parseCSV() const {
	std::string filename = basic::database::full_name( "protocol_data/tensorflow_graphs/tensorflow_graph_repo_submodule/predict_PTMs/deamidation_rates.csv" );
	std::ifstream data( filename );
	std::string line;
	std::vector< std::vector< float > > parsedCsv;
	while( std::getline( data, line ) ) 
	{
			std::stringstream lineStream( line );
			std::string cell;
			std::vector< float > parsedRow;
			while( std::getline( lineStream, cell, ',' ) )
			{
					parsedRow.push_back( std::stof ( cell ) );
			}
			parsedCsv.push_back(parsedRow);
	}
	return parsedCsv;
}

/// @brief Helper function to setup in-/output names/shapes for different PTMs
void 
PTMPredictionTensorflowProtocol::setup_and_check( 
  std::string const & aa_select,
  PTMPredictionMetric::PTMPredictionMetricModes const & ptm_mode,
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & output_tensor,
	size_t & index,
	utility::vector1< std::string > &input_names
) const {
		bool is_multi = true; // default set to true to keep compiler happy in switch statement
		std::string expected_aa;

		switch( ptm_mode ) {

				case PTMPredictionMetric::PTMPredictionMetricModes::N_LINKED_GLYCOSYLATION:
						index = 5;
						input_names.push_back( "serving_default_input_5" );
						input_names.push_back( "serving_default_input_6" );
						is_multi = true;
						expected_aa = "ASN";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::HYDROXYLATION:
						index = 7;
						input_names.push_back( "serving_default_input_3" );
						input_names.push_back( "serving_default_input_4" );
						is_multi = true;
						expected_aa = "PRO";
						break;

			case PTMPredictionMetric::PTMPredictionMetricModes::GAMMA_CARBOXY_GLUTAMIC_ACID:
						index = 3;
						input_names.push_back( "serving_default_input_17" );
						input_names.push_back( "serving_default_input_18" );
						is_multi = true;
						expected_aa = "GLU";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::LYS_METHYLATION:
						index = 4;
						input_names.push_back( "serving_default_input_19" );
						input_names.push_back( "serving_default_input_20" );
						is_multi = true;
						expected_aa = "LYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::MALONYLATION:
						index = 1;
						input_names.push_back( "serving_default_input_17" );
						input_names.push_back( "serving_default_input_18" );
						is_multi = false;
						expected_aa = "LYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::ARG_METHYLATION:
						index = 1;
						input_names.push_back( "serving_default_input_3" );
						input_names.push_back( "serving_default_input_4" );
						is_multi = true;
						expected_aa = "LYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::CROTONYLATION:
						index = 2;
						input_names.push_back( "serving_default_input_11" );
						input_names.push_back( "serving_default_input_12" );
						is_multi = true;
						expected_aa = "LYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::UBIQUITINATION:
						index = 7;
						input_names.push_back( "serving_default_input_19" );
						input_names.push_back( "serving_default_input_20" );
						is_multi = true;
						expected_aa = "LYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::SUCCINYLATION:
						index = 1;
						input_names.push_back( "serving_default_input_19" );
						input_names.push_back( "serving_default_input_20" );
						is_multi = false;
						expected_aa = "LYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::GLUTATHIONYLATION:
						index = 1;
						input_names.push_back( "serving_default_input_17" );
						input_names.push_back( "serving_default_input_18" );
						is_multi = false;
						expected_aa = "CYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::SUMOYLATION:
						index = 7;
						input_names.push_back( "serving_default_input_17" );
						input_names.push_back( "serving_default_input_18" );
						is_multi = true;
						expected_aa = "LYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::S_NITROSYLATION:
						index = 1;
						input_names.push_back( "serving_default_input_11" );
						input_names.push_back( "serving_default_input_12" );
						is_multi = false;
						expected_aa = "CYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::ACETYLATION:
						index = 1;
						input_names.push_back( "serving_default_input_1" );
						input_names.push_back( "serving_default_input_2" );
						is_multi = true;
						expected_aa = "LYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::O_LINKED_GLYCOSYLATION:
						index = 6;
						input_names.push_back( "serving_default_input_13" );
						input_names.push_back( "serving_default_input_14" );
						is_multi = true;
						if ( aa_select == "SER" || aa_select == "THR" ) {
								expected_aa = aa_select;
						} else {
								expected_aa = "SER/THR";
						}
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::PHOSPHORYLATION:
						index = 7;
						input_names.push_back( "serving_default_input_15" );
						input_names.push_back( "serving_default_input_16" );
						is_multi = true;
						if ( aa_select == "SER" || aa_select == "THR" ) {
								expected_aa = aa_select;
						} else {
								expected_aa = "SER/THR";
						}
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::GLUTARYLATION:
						index = 3;
						input_names.push_back( "serving_default_input_5" );
						input_names.push_back( "serving_default_input_6" );
						is_multi = true;
						expected_aa = "LYS";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::CITRULLINATION:
						index = 1;
						input_names.push_back( "serving_default_input_17" );
						input_names.push_back( "serving_default_input_18" );
						is_multi = true;
						expected_aa = "ARG";
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::DEAMIDATION:
						utility_exit_with_message( "Deamidation should automatically lead to another function before this, if you see this message sth is very wrong, its just here to make the switch statement complete and the compiler happy" );
						break;

				case PTMPredictionMetric::PTMPredictionMetricModes::INVALID_MODE:
						utility_exit_with_message( "Chosen modification does not match any of the available options!" );
						break;

		}

		if ( is_multi ) {
				output_tensor.initialize( TF_FLOAT, {1, 8 } );
		}
		else {
				output_tensor.initialize( TF_FLOAT, {1, 2 } );
		}
		runtime_assert_string_msg( aa_select == expected_aa, "Error in PTMPredictionTensorflowProtocol::setup_and_check(): The " + name() + " module expects an " + expected_aa + ", but an " + aa_select + " was selected instead!" );
}
#endif // USE_TENSORFLOW
} // ptm_prediction
} // protocols
