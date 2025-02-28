// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/WorkManager.hh
/// @brief  Class declaration for %WorkManager
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_WorkManager_HH
#define INCLUDED_protocols_ligand_evolution_WorkManager_HH

// unit headers
#include <protocols/ligand_evolution/WorkManager.fwd.hh>

#ifdef USEMPI

// package headers
#include <protocols/ligand_evolution/Scorer.hh>

// C/C++ headers
#include <mpi.h>
#include <stack>

namespace protocols {
namespace ligand_evolution {

	/// @brief Manages parallel processing of ligand scoring and all associated steps
	class WorkManager {
	public:

		explicit WorkManager( ScorerOP& scorer, core::Size id_length, FragmentLibrary const& library );
		~WorkManager();

		/// @brief Scans all individuals for unscored ligands and starts working on them
		void score( Population& pop );

		/// @brief If this process is not rank_ == 0, it waits for tasks or a termination signal with non-blocking receives
		void work_loop();

		/// @brief Returns rank of the current process
		core::Size get_rank() const;

		/// @brief Deletes all workers to terminate their processes
		void clean_up();

	private:

		/// @brief Distributes the work on all unscored ligands to the workers
		void distribute_work();

		/// @brief Checks in on all workers if they have scores ready and retrieves them if possible
		void retrieve_scores();

		/// @brief Returns the index of the next idle process or 0 if non are idle
		core::Size next_idle_process() const;

		/// @brief Returns true if all processes are idle
		bool all_processes_idle() const;

	private:

		/// @brief Stores the total number of processes used
		int n_processes_ = 1;

		/// @brief Stores the rank of this process
		int rank_ = 0;

		/// @brief The element at top should be the next worked on ligand. Expects all unscored and unique ligands
		std::stack< LigandIdentifier > unscored_ligands_;

		ScorerOP scorer_;

		utility::vector1< WorkerOP > worker_;

		FragmentLibrary const& library_;

	};

	/// @brief Keeps track of a single task send out to another process and retrieves result
	class Worker {
	public:

		/// @brief Sets rank, number of scores and length of ids
		Worker( int processor_rank, core::Size n_scores, core::Size id_length );

		/// @brief Terminates the associated process when deleting object
		~Worker();

		// deleted these because copying is 100% undesired and will crash the system
		Worker( Worker const& ) = delete;
		Worker& operator=( Worker const& ) = delete;

		/// @brief Blocking sends to the associated process and starts a non-blocking receive for scores
		void send_ligand();

		/// @brief Retrieves scores if arrived
		double* retrieve_scores();

		/// @brief Checks if scores are received
		bool has_scores() const;

		/// @brief Sets the ligand that should be worked on
		void set_ligand( LigandIdentifier const& ligand );

		/// @brief Sets the ligand that should be worked on
		void set_smiles( std::string const& ligand );

		/// @brief checks if a termination message has arrived
		bool check_termination() const;

		/// @brief checks if a task has arrived
		bool check_task() const;

		/// @brief checks if a smile has arrived
		bool check_smiles() const;

		/// @brief Sets the scores that should be send
		void set_scores( utility::vector1< core::Real > raw_scores );

		/// @brief Sends the calculated scores to
		void send_scores() const;

		/// @brief If a task is available, this sets the current ligand and starts a new non-blocking receive
		void receive_ligand();

		/// @brief Returns true if the associated process is currently idle
		bool is_idle() const;

		/// @brief Returns the ligand set for this process
		LigandIdentifier ligand() const;

		/// @brief Returns the smiles set for this process
		std::string smiles() const;

	private:

		/// @brief Sends a termination signal to the associated process
		void terminate() const;

	private:

		/// @brief Tells the Manager if this worker accepts tasks
		bool idle_ = true;

		/// @brief The rank of the processor with which this Worker communicates
		int processor_rank_;

		/// @brief The ligand this worker is currently working on
		LigandIdentifier ligand_;

		/// @brief The smiles for the currently worked ligand
		std::string ligand_smiles_;

		/// @brief Stores how many score terms will be calculated
		core::Size n_score_terms_;

		/// @brief The scores calculated by the other process
		std::unique_ptr< double[] > raw_scores_;

		/// @brief Request handle to check on non-blocking communication for termination
		std::unique_ptr< MPI_Request > termination_handle_ = nullptr;

		/// @brief Request handle to check on non-blocking communication for tasks
		std::unique_ptr< MPI_Request > task_handle_ = nullptr;

		/// @brief Request handle to check on non-blocking communication for scores
		std::unique_ptr< MPI_Request > scores_handle_ = nullptr;

		/// @brief Request handle to check on non-blocking communication for smiles
		std::unique_ptr< MPI_Request> smiles_handle_ = nullptr;

		/// @brief Stores the result of termination message. Unused.
		int termination_in_buffer_ = 0;

		/// @brief A reference how long an expected ligand identifier will be
		core::Size id_length_;

		/// @brief Stores the result of a task message
		std::unique_ptr< int[] > raw_ligand_;

		/// @brief Stores the result of a task message sending a ligand smiles
		std::unique_ptr< char[] > raw_ligand_smiles_;

		static int const max_smiles_length_ = 2500;

	};

	enum Tag {
		TERMINATION,
		LIGAND,
		SCORES,
		SMILES
	};

}
}

#endif  // include guards
#endif  // USEMPI
