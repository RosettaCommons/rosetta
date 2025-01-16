// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/WorkManager.cc
/// @brief  The definition for %WorkManager class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifdef USEMPI

// unit headers
#include <protocols/ligand_evolution/WorkManager.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/stream_util.hh>

static basic::Tracer TR( "protocols.ligand_evolution.WorkManager" ); // NOLINT(cert-err58-cpp)

namespace protocols {
namespace ligand_evolution {

    // #########################################################################
    // WorkManager
    // #########################################################################

    WorkManager::WorkManager( ScorerOP& scorer, core::Size id_length, FragmentLibrary const& library )
    :
    scorer_( scorer ),
    library_( library )
    {
        MPI_Comm_size (MPI_COMM_WORLD, &n_processes_);
        MPI_Comm_rank (MPI_COMM_WORLD, &rank_);
        if( rank_ != 0 ) {
            worker_.emplace_back( new Worker( 0, scorer_->n_score_terms(), id_length ) );
        } else {
            for( int ii = 1; ii < n_processes_; ++ii ) {
                worker_.emplace_back( new Worker( ii, scorer_->n_score_terms(), id_length ) );
            }
        }
    }

    WorkManager::~WorkManager() {
        clean_up();
    }

    void WorkManager::score( Population& pop ) {

        std::set< LigandIdentifier > added_ligands;
        std::set< LigandIdentifier > encountered_ligands;

        core::Size ligand_duplicates = 0;
        core::Size ligands_in_memory = 0;
        core::Size ligands_scored = 0;

        for( Individual const& individual : pop.individuals() ) {
            LigandIdentifier const& identifier = individual.identifier();

            bool ligand_is_duplicate = !encountered_ligands.count( identifier ) == 0;
            bool ligand_is_scored = scorer_->is_scored( identifier );
            bool score_in_memory = scorer_->check_memory( identifier );

            encountered_ligands.insert( identifier );

            if ( score_in_memory && !ligand_is_duplicate && !ligand_is_scored ) ligands_in_memory++;
            if ( ligand_is_scored && !ligand_is_duplicate ) ligands_scored++;
            if ( ligand_is_duplicate ) ligand_duplicates++;

            if( !score_in_memory && !ligand_is_scored && !ligand_is_duplicate ) {
                added_ligands.insert( identifier );
                unscored_ligands_.push( identifier );
            }
        }

        TR << "Received " << pop.size() << " ligands. " << ligands_in_memory << " scores are taken from file, " << ligands_scored << " were encountered before and " << ligand_duplicates << " are duplicates in this population. Start working on " << added_ligands.size() << " remaining ligands." << std::endl;

        distribute_work();

        TR.Debug << "Finished all work." << std::endl;

        scorer_->score_population( pop );

    }

    void WorkManager::distribute_work() {

        // work until all ligands are scored
        while( !unscored_ligands_.empty() ) {

            LigandIdentifier ligand = unscored_ligands_.top();

            // TODO this can be optimized when I give the main process only short ligands so that it does not need to calculate rotamers for 30 seconds

            // 1 check if a process is done calculating and delivered scores
            retrieve_scores();
            // 2 check for idle processes and give them some work
            // 3 start working yourself on a score if there is nothing to do
            core::Size idle_process = next_idle_process();
            if( idle_process != 0 ) {
                // there is an idle process so give him something to work on
                worker_[ idle_process ]->set_ligand( ligand );
                worker_[ idle_process ]->set_smiles( library_.run_reaction( ligand ) );
                unscored_ligands_.pop();
                worker_[ idle_process ]->send_ligand();
            } else {
                // check if the scorer has a ligand set
                if( !scorer_->has_ligand() ) {
                    scorer_->set_ligand( ligand );
                    unscored_ligands_.pop();
                }
                scorer_->next_step();
            }

            TR.Debug << unscored_ligands_.size() << " unscored ligands not distributed." << std::endl;

        }

        TR.Debug << "Distributed all tasks." << std::endl;

        // finish your own work if something is left undone
        if( scorer_->has_ligand() ) {
            while( !scorer_->next_step ());
        }

        TR.Debug << "Finished own work." << std::endl;

        // check that all processes are idle again
        while( !all_processes_idle() ) {
            retrieve_scores();
        }

        TR.Debug << "All processes answered." << std::endl;

    }

    void WorkManager::work_loop() {

        if( rank_ == 0 ) {
            TR.Error << "'work_loop()' should not be called by process with rank 0." << std::endl;
            utility_exit_with_message( "Mixed up processor tasks" );
        }

        TR.Debug << "Select worker" << std::endl;
        // This WorkManager should only have one worker to handle communication with the Manager process
        Worker& worker = *worker_.front();
        TR.Debug << "Start working" << std::endl;

        while( !worker.check_termination() ) {
            if( worker.check_task() && worker.check_smiles() ) {
                // the worker has received a ligand to work on
                worker.receive_ligand();
                LigandIdentifier ligand( worker.ligand() );
                std::string smiles( worker.smiles() );

                // score the ligand
                scorer_->score_ligand( ligand, smiles );
                worker.set_scores( scorer_->get_raw_scores( ligand ) );

                TR.Debug << "Send scores..." << std::endl;
                // send the scores back
                worker.send_scores();
                TR.Debug << "Wait for new task." << std::endl;
            }
        }


    }

    core::Size WorkManager::next_idle_process() const {
        for( core::Size ii = 1; ii <= worker_.size(); ++ii ) {
            if( worker_[ ii ]->is_idle() ) {
                return ii;
            }
        }
        return 0;
    }

    void WorkManager::retrieve_scores() {
        for( core::Size ii = 1; ii <= worker_.size(); ++ii ) {
            if( !worker_[ ii ]->is_idle() ) {
                if( worker_[ ii ]->has_scores() ) {
                    double* raw_scores = worker_[ ii ]->retrieve_scores();
                    scorer_->set_raw_scores( worker_[ ii ]->ligand(), raw_scores );
                }
            }
        }
    }

    core::Size WorkManager::get_rank() const {
        return rank_;
    }

    bool WorkManager::all_processes_idle() const {
        for( core::Size ii = 1; ii <= worker_.size(); ++ii ) {
            if( !worker_[ ii ]->is_idle() ) {
                return false;
            }
        }
        return true;
    }

    void WorkManager::clean_up() {
        worker_.clear();
    }

    // #########################################################################
    // Worker
    // #########################################################################

    Worker::Worker( int processor_rank, core::Size n_scores, core::Size id_length )
    :
            processor_rank_( processor_rank ),
            n_score_terms_( n_scores ),
            raw_scores_( new double[ n_score_terms_ ] ),
            termination_handle_( new MPI_Request ),
            task_handle_( new MPI_Request ),
            scores_handle_( new MPI_Request ),
            smiles_handle_( new MPI_Request ),
            id_length_( id_length ),
            raw_ligand_( new int[ id_length_ ] ),
            raw_ligand_smiles_( new char[ max_smiles_length_ ] )
    {
        if( processor_rank_ == 0 ) {
            // starts a non blocking receiving for the termination signal
            MPI_Irecv( &termination_in_buffer_, 1, MPI_INT, processor_rank_, TERMINATION, MPI_COMM_WORLD, termination_handle_.get() );
            // starts a non blocking receiving for a ligand
            MPI_Irecv( raw_ligand_.get(), int( id_length_ ), MPI_INT, processor_rank_, LIGAND, MPI_COMM_WORLD, task_handle_.get() );
            // starts a non blocking receiving for a ligand smiles
            MPI_Irecv( raw_ligand_smiles_.get(), max_smiles_length_, MPI_CHAR, processor_rank_, SMILES, MPI_COMM_WORLD, smiles_handle_.get() );
        }
    }

    Worker::~Worker() {
        if( processor_rank_ == 0 ) {
            // maybe I need to call mpi cancel here to stop receiving, but so far I have to yet encounter a problem
        } else {
            terminate();
        }
    }

    // only called by rank != 0
    bool Worker::check_termination() const {
        int status = 0;
        MPI_Test( termination_handle_.get(), &status, MPI_STATUS_IGNORE );
        // 0 values of basic types map to false, all others to true
        return status;
    }

    // only called by rank != 0
    bool Worker::check_task() const {
        int status = 0;
        MPI_Test( task_handle_.get(), &status, MPI_STATUS_IGNORE );
        // 0 values of basic types map to false, all others to true
        return status;
    }

    bool Worker::check_smiles() const {
        int status = 0;
        MPI_Test( smiles_handle_.get(), &status, MPI_STATUS_IGNORE );
        // 0 values of basic types map to false, all others to true
        return status;
    }

    // only called by rank != 0
    void Worker::receive_ligand() {
        if( check_task() && check_smiles() ) {
            ligand_.clear();
            ligand_.resize( id_length_ );
            for( core::Size ii = 1; ii <= id_length_; ++ii ) {
                ligand_[ ii ]  = core::Size( raw_ligand_[ ii - 1 ] );
            }
            ligand_smiles_ = utility::strip( std::string( raw_ligand_smiles_.get() ).substr( 0, 990 ), ' ' );
            TR.Debug << "Received ligand " << ligand_ << " and its smiles " << ligand_smiles_ << std::endl;
            // starts a new non-blocking receiving for a ligand
            MPI_Irecv( raw_ligand_.get(), int( id_length_ ), MPI_INT, processor_rank_, LIGAND, MPI_COMM_WORLD, task_handle_.get() );
            // starts a new non-blocking receiving for a ligand smiles
            MPI_Irecv( raw_ligand_smiles_.get(), max_smiles_length_, MPI_CHAR, processor_rank_, SMILES, MPI_COMM_WORLD, smiles_handle_.get() );
        }
    }

    bool Worker::is_idle() const {
        return idle_;
    }

    LigandIdentifier Worker::ligand() const {
        return ligand_;
    }

    // only called by rank != 0
    void Worker::set_scores( utility::vector1< core::Real > raw_scores ) {
        for( core::Size ii = 1; ii <= n_score_terms_; ++ii ) {
            raw_scores_[ ii - 1 ] = raw_scores[ ii ];
        }
    }

    // only called by rank != 0
    void Worker::send_scores() const {
        MPI_Send( raw_scores_.get(), int( n_score_terms_ ), MPI_DOUBLE, processor_rank_, SCORES, MPI_COMM_WORLD );
    }

    // only called by rank == 0
    void Worker::terminate() const {
        int terminate = 1;
        MPI_Send( &terminate, 1, MPI_INT, processor_rank_, TERMINATION, MPI_COMM_WORLD );
    }

    // only called by rank == 0
    void Worker::set_ligand( LigandIdentifier const& ligand ) {
        TR.Debug << "Set ligand " << ligand << " for process " << processor_rank_ << std::endl;
        ligand_ = ligand;
    }

    // only called by rank == 0
    void Worker::set_smiles(std::string const& smiles) {
        TR.Debug << "Set smiles " << smiles << " for process " << processor_rank_ << std::endl;
        ligand_smiles_ = smiles;
    }

    // only called by rank == 0
    void Worker::send_ligand() {

        if( ligand_.empty() ) {
            TR.Error << "Cannot send empty ligand." << std::endl;
            utility_exit_with_message( "No ligand was set" );
        }

        if( ligand_smiles_.empty() ) {
            TR.Error << "Cannot send empty smiles." << std::endl;
            utility_exit_with_message( "No smiles was set" );
        }

        for( core::Size ii = 1; ii <= id_length_; ++ii ) {
            raw_ligand_[ ii - 1 ] = int( ligand_[ ii ] );
        }

        // add padding to always send fixed size smiles
        if ( ligand_smiles_.size() > max_smiles_length_ ) {
            TR.Error << "Ligand " << ligand_ << " with smiles " << ligand_smiles_ << " exceeds maximum length of " << max_smiles_length_ << " for MPI_Send." << std::endl;
            utility_exit_with_message( "Smiles to long to process with MPI send protocol." );
        }
        ligand_smiles_ = utility::pad_right( ligand_smiles_, max_smiles_length_, ' ' );

        MPI_Send( raw_ligand_.get(), int( id_length_ ), MPI_INT, processor_rank_, LIGAND, MPI_COMM_WORLD );
        MPI_Send( ligand_smiles_.c_str(), int( ligand_smiles_.size() ), MPI_CHAR, processor_rank_, SMILES, MPI_COMM_WORLD );
        MPI_Irecv( raw_scores_.get(), int( n_score_terms_ ), MPI_DOUBLE, processor_rank_, SCORES, MPI_COMM_WORLD, scores_handle_.get() );
        idle_ = false;
    }

    // only called by rank == 0
    bool Worker::has_scores() const {
        int status = 0;
        MPI_Test( scores_handle_.get(), &status, MPI_STATUS_IGNORE );
        // 0 values of basic types map to false, all others to true
        return status;
    }

    // only called by rank == 0
    double* Worker::retrieve_scores() {
        if( has_scores() ) {
            idle_ = true;
            return raw_scores_.get();
        } else {
            TR.Error << "No scores are ready." << std::endl;
            utility_exit_with_message( "Tried to access scores before receiving" );
        }
    }

    std::string Worker::smiles() const {
        return ligand_smiles_;
    }
}
}

#endif
