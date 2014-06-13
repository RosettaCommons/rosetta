// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/UnfoldedStateEnergyCalculator/UnfoldedStateEnergyCalculatorJobDistributor.cc
/// @brief  Job distributor for UnfoldedStateEnergyCalculator
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit headers
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor.hh>
#ifdef USEMPI
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorUtil.hh>
#endif
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>

// Package headers
// AUTO-REMOVED #include <protocols/moves/Mover.hh>

// Project headers
#include <core/scoring/ScoreType.hh>

// Utility headers
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
#ifdef USEMPI
#include <utility/exit.hh>
#endif
#include <utility/assert.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>


static basic::Tracer TR("protocols.UnfoldedStateEnergyCalculator.UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor");

namespace protocols {
namespace unfolded_state_energy_calculator {

using namespace basic::options;
using namespace basic::options::OptionKeys;

///@brief ctor
UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor::UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor() :
  MPIWorkPoolJobDistributor()
{}

///@brief dtor (don't put anything in here)
UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor::~UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor()
{}

///@brief unforntunatly this is pretty much copied from the MPIWorkPoolJobDistributor, I should make that more compartmentalized
void
UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor::master_go( protocols::moves::MoverOP /*mover*/ )
{
#ifdef USEMPI

	using namespace core;
	using namespace core::scoring;
	using namespace protocols::jd2;

	runtime_assert( rank_ == 0 );

	// local arrays to store incomming messages
	int slave_data( 0 );
	double slave_data_vector[ n_score_types ];
	EMapVector slave_data_emap;
	char slave_tlc_vector[LENGTH_TLC];
	std::string slave_tlc_string;


	MPI_Status status;

	// set first job to assign
	master_get_new_job_id();

	// Job Distribution Loop
	while ( next_job_to_assign_ != 0 ) {
		TR << "Master Node: Waiting for job requests..." << std::endl;
		MPI_Recv( &slave_data, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		TR << "Master Node: Received message from " << status.MPI_SOURCE << " with tag " << status.MPI_TAG << std::endl;

		// decide what to do based on message tag
		switch ( status.MPI_TAG ) {
		case NEW_JOB_ID_TAG:
			TR << "Master Node: Sending new job id " << next_job_to_assign_ << " to node " << status.MPI_SOURCE << " with tag " << NEW_JOB_ID_TAG << std::endl;
			MPI_Send( &next_job_to_assign_, 1, MPI_INT, status.MPI_SOURCE, NEW_JOB_ID_TAG, MPI_COMM_WORLD );
			master_get_new_job_id();
			break;
		case BAD_INPUT_TAG:
			TR << "Master Node: Received job failure message for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
			bad_job_id_ = slave_data;
			master_remove_bad_inputs_from_job_list();
			break;
		case JOB_SUCCESS_TAG:
			TR << "Master Node: Received job success message for job id " << slave_data << " from node " << status.MPI_SOURCE << " blocking till output is done " << std::endl;
			MPI_Send( &next_job_to_assign_, 1, MPI_INT, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD );
			MPI_Recv( &slave_data, 1, MPI_INT, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD, &status);
			TR << "Master Node: Received job output finish messege for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
			break;
		case UNFOLDED_ENERGY_DATA_TAG:
			TR << "Master Node: Received unfolded energy data message for job id " << slave_data << " from node " << status.MPI_SOURCE << " waiting for data " << std::endl;
			// clear local data
			slave_data_emap.clear(); slave_tlc_string.clear();

			// send responce to send data
			MPI_Send( &slave_data, 1, MPI_INT, status.MPI_SOURCE, UNFOLDED_ENERGY_DATA_TAG, MPI_COMM_WORLD );

			// get vector of tlc data
			MPI_Recv( &slave_data_vector, LENGTH_TLC, MPI_CHAR, status.MPI_SOURCE, UNFOLDED_ENERGY_DATA_TAG, MPI_COMM_WORLD, &status);

			// convert to string
			for( int i(0); i < LENGTH_TLC; ++i ) {
				slave_tlc_string += slave_tlc_vector[i];
			}

			// get vector of energy data
			MPI_Recv( &slave_data_vector, n_score_types, MPI_DOUBLE, status.MPI_SOURCE, UNFOLDED_ENERGY_DATA_TAG, MPI_COMM_WORLD, &status);

			// convert to emap
			for ( Size i( 1 ); i <= n_score_types; ++i ) {
				slave_data_emap[ ScoreType( i ) ] = slave_data_vector[ i - 1 ];
			}

			// push into vector
			master_add_unfolded_energy_data( slave_tlc_string, slave_data_emap );
			TR << "Master Node: Finished receiving energy data message for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
			break;
		case UNFOLDED_ENERGY_TERMS_TAG:
			TR << "Master Node: Received unfolded energy terms message for job id " << slave_data << " from node " << status.MPI_SOURCE << " waiting for data " << std::endl;
			// clear emap
			slave_data_emap.clear();

			// send responce to send data
			MPI_Send( &slave_data, 1, MPI_INT, status.MPI_SOURCE, UNFOLDED_ENERGY_TERMS_TAG, MPI_COMM_WORLD );

			// get vector of data
			MPI_Recv( &slave_data_vector, n_score_types, MPI_DOUBLE, status.MPI_SOURCE, UNFOLDED_ENERGY_TERMS_TAG, MPI_COMM_WORLD, &status);

			// convert to emap
			for ( Size i( 1 ); i <= n_score_types; ++i ) {
				slave_data_emap[ ScoreType( i ) ] = slave_data_vector[ i - 1 ];
			}

			// push into vector
			master_set_energy_terms( slave_data_emap );
			TR << "Master Node: Finished receiving energy terms message for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
			break;
		}
	}
	TR << "Master Node: Finished handing out jobs" << std::endl;

	core::Size n_nodes_left_to_spin_down( npes_ - 1 ); // don't have to spin down self

	// Node Spin Down loop
	while ( n_nodes_left_to_spin_down > 0 ) {
		TR << "Master Node: Waiting for " << n_nodes_left_to_spin_down << " slaves to finish jobs" << std::endl;
		MPI_Recv( &slave_data, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		TR << "Master Node: Recieved message from  " << status.MPI_SOURCE << " with tag " << status.MPI_TAG << std::endl;

		// decide what to do based on message tag
		switch ( status.MPI_TAG ) {
		case NEW_JOB_ID_TAG:
			TR << "Master Node: Sending spin down signal to node " << status.MPI_SOURCE << std::endl;
			MPI_Send( &next_job_to_assign_, 1, MPI_INT, status.MPI_SOURCE, NEW_JOB_ID_TAG, MPI_COMM_WORLD );
			n_nodes_left_to_spin_down--;
			break;
		case BAD_INPUT_TAG:
			break;
		case JOB_SUCCESS_TAG:
			TR << "Master Node: Received job success message for job id " << slave_data << " from node " << status.MPI_SOURCE << " blocking till output is done " << std::endl;
			MPI_Send( &next_job_to_assign_, 1, MPI_INT, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD );
			MPI_Recv( &slave_data, 1, MPI_INT, status.MPI_SOURCE, JOB_SUCCESS_TAG, MPI_COMM_WORLD, &status);
			TR << "Master Node: Received job output finish messege for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
			break;
		case UNFOLDED_ENERGY_DATA_TAG:
			TR << "Master Node: Received unfolded energy data message for job id " << slave_data << " from node " << status.MPI_SOURCE << " waiting for data " << std::endl;
			// clear local data
			slave_data_emap.clear(); slave_tlc_string.clear();

			// send responce to send data
			MPI_Send( &slave_data, 1, MPI_INT, status.MPI_SOURCE, UNFOLDED_ENERGY_DATA_TAG, MPI_COMM_WORLD );

			// get vector of tlc data
			MPI_Recv( &slave_data_vector, LENGTH_TLC, MPI_CHAR, status.MPI_SOURCE, UNFOLDED_ENERGY_DATA_TAG, MPI_COMM_WORLD, &status);

			// convert to string
			for( int i(0); i < LENGTH_TLC; ++i ) {
				slave_tlc_string += slave_tlc_vector[i];
			}

			// get vector of data
			MPI_Recv( &slave_data_vector, n_score_types, MPI_DOUBLE, status.MPI_SOURCE, UNFOLDED_ENERGY_DATA_TAG, MPI_COMM_WORLD, &status);

			// convert to emap
			for ( Size i( 1 ); i <= n_score_types; ++i ) {
				slave_data_emap[ ScoreType( i ) ] = slave_data_vector[ i - 1 ];
			}

			// push into vector
			master_add_unfolded_energy_data( slave_tlc_string, slave_data_emap );
			TR << "Master Node: Finished receiving energy data message for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
			break;
		case UNFOLDED_ENERGY_TERMS_TAG:
			TR << "Master Node: Received unfolded energy terms message for job id " << slave_data << " from node " << status.MPI_SOURCE << " waiting for data " << std::endl;
			// clear emap
			slave_data_emap.clear();

			// send responce to send data
			MPI_Send( &slave_data, 1, MPI_INT, status.MPI_SOURCE, UNFOLDED_ENERGY_TERMS_TAG, MPI_COMM_WORLD );

			// get vector of data
			MPI_Recv( &slave_data_vector, n_score_types, MPI_DOUBLE, status.MPI_SOURCE, UNFOLDED_ENERGY_TERMS_TAG, MPI_COMM_WORLD, &status);

			// convert to emap
			for ( Size i( 1 ); i <= n_score_types; ++i ) {
				slave_data_emap[ ScoreType( i ) ] = slave_data_vector[ i - 1 ];
			}

			// push into vector
			master_set_energy_terms( slave_data_emap );
			TR << "Master Node: Finished receiving energy terms message for job id " << slave_data << " from node " << status.MPI_SOURCE << std::endl;
			break;
		}
	}
	TR << "Master Node: Finished sending spin down signals to slaves" << std::endl;

	// calc average unweighted energies for all amino acids in the map
	for ( uem_iter i( unweighted_energies_map_.begin() ), e( unweighted_energies_map_.end() ); i != e; ++i ) {
		TR << "Calculating averages for " << i->first << std::endl;
		calc_all_averages( i->second, energy_terms_ );
	}

#endif
}

void
UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor::add_unfolded_energy_data( std::string tlc, core::scoring::EMapVector const & scores )
{
	if ( rank_ == 0 ) {
    master_add_unfolded_energy_data( tlc, scores );
  } else {
    slave_add_unfolded_energy_data( tlc, scores );
  }
}

void
UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor::master_add_unfolded_energy_data( std::string tlc, core::scoring::EMapVector const & scores )
{
	unweighted_energies_map_[tlc].push_back( scores );
}

void
UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor::slave_add_unfolded_energy_data( std::string MPI_ONLY( tlc ), core::scoring::EMapVector const & MPI_ONLY( scores ) )
{
#ifdef USEMPI
	using namespace core;
	using namespace core::scoring;

	TR << "Slave Node " << rank_ << ": Requesting to send unfolded energy data to master" <<std::endl;

	int empty_data( 0 );
	MPI_Status status;
	double data_vector[ n_score_types ];
	char tlc_vector[LENGTH_TLC];

	// send message that we want to send a vector of energies
	MPI_Send( &current_job_id_, 1, MPI_INT, 0, UNFOLDED_ENERGY_DATA_TAG, MPI_COMM_WORLD );

	// get responce from master
	MPI_Recv( &empty_data, 1, MPI_INT, 0, UNFOLDED_ENERGY_DATA_TAG, MPI_COMM_WORLD, &status );
	TR << "Slave Node " << rank_ << ": Sending unfolded energy data to master" <<std::endl;

	// convert string to vector of chars
	for ( int i(0); i < LENGTH_TLC; ++i ) {
		tlc_vector[i] = tlc[i];
	}

	// send vector to master
	MPI_Send( &tlc_vector, 3, MPI_CHAR, 0, UNFOLDED_ENERGY_DATA_TAG, MPI_COMM_WORLD );

	// convert EMapVector to vector of doubles
	for ( Size i( 1 ); i <= n_score_types; ++i ) {
		data_vector[ i - 1 ] = scores[ ScoreType( i ) ];
	}

	// send vector to master
	MPI_Send( &data_vector, n_score_types, MPI_DOUBLE, 0, UNFOLDED_ENERGY_DATA_TAG, MPI_COMM_WORLD );

	TR << "Slave Node " << rank_ << ": Sent unfolded energy data to master" <<std::endl;
#endif
}

void
UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor::set_energy_terms(core::scoring::EMapVector const & weights )
{
	using namespace protocols::jd2;

	if ( rank_ == 0 ) {
    master_set_energy_terms( weights );
  } else {
		slave_set_energy_terms( weights );
  }
}

void
UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor::master_set_energy_terms(core::scoring::EMapVector const & weights )
{
	using namespace core;
	using namespace core::scoring;

	// get weights
	energy_terms_ = weights;

	// for each energy term in the EMapVector
	for ( Size i( 1 ); i <= n_score_types; ++i ) {

		// if the energy term has a non-zero weight set it to one
		if ( energy_terms_[ ScoreType( i ) ] > 0 ) {
			energy_terms_[ ScoreType( i ) ] = 1;
		} else if ( energy_terms_[ ScoreType( i ) ] < 0 ) {
			energy_terms_[ ScoreType( i ) ] = -1;
		}
	}
}

void
UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor::slave_set_energy_terms(core::scoring::EMapVector const & MPI_ONLY( weights ) )
{
#ifdef USEMPI
	using namespace core;
	using namespace core::scoring;

	TR << "Slave Node " << rank_ << ": Requesting to send unfolded energy terms to master" <<std::endl;

	int empty_data( 0 );
	MPI_Status status;
	double data_vector[ n_score_types ];

	// send message that we want to send a vector of energies
	MPI_Send( &current_job_id_, 1, MPI_INT, 0, UNFOLDED_ENERGY_TERMS_TAG, MPI_COMM_WORLD );

	// get responce from master
	MPI_Recv( &empty_data, 1, MPI_INT, 0, UNFOLDED_ENERGY_TERMS_TAG, MPI_COMM_WORLD, &status );
	TR << "Slave Node " << rank_ << ": Sending unfolded energy terms to master" <<std::endl;

	// convert EMapVector to vector of doubles
	for ( Size i( 1 ); i <= n_score_types; ++i ) {
		data_vector[ i - 1 ] = weights[ ScoreType( i ) ];
	}

	// send vector to master
	MPI_Send( &data_vector, n_score_types, MPI_DOUBLE, 0, UNFOLDED_ENERGY_TERMS_TAG, MPI_COMM_WORLD );

	TR << "Slave Node " << rank_ << ": Sent unfolded energy terms to master" <<std::endl;
#endif
}


} // UnfoldedStateEnergyCalculator
} // protocols
