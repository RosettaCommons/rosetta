// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MPIWorkPoolJobDistributor.hh
/// @brief  header for MPIWorkPoolJobDistributor - intended for continuous resamplig jobs  that spawn new jobs based on a pool/archive of
///         structures
/// @author Oliver Lange olange@u.washington.edu

#ifndef INCLUDED_protocols_jd2_archive_ArchiveManager_hh
#define INCLUDED_protocols_jd2_archive_ArchiveManager_hh

// Unit headers
#include <protocols/jd2/archive/ArchiveManager.fwd.hh>
#include <protocols/jd2/archive/ArchiveBase.hh>
#include <protocols/jd2/archive/MPIArchiveJobDistributor.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <core/io/silent/silent.fwd.hh>


#include <protocols/moves/Mover.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/excn/Exceptions.hh>

// Utility headers
#include <core/types.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {
namespace archive {

class EXCN_Archive : public utility::excn::EXCN_Msg_Exception {
public:
	EXCN_Archive( std::string const& msg ) : EXCN_Msg_Exception( msg ) {};
};

extern void report_batch_inconsistency( Batch& new_batch, std::string const &tag );

/// @brief a Batch represents a directory "batch_000xxx" that contains flags, broker-setup input-files and output-files
/// @detail the Batch-class helps to get the correct file- and directory names,
///          and has some knowledge about its status: finished, unfinished ... decoys already processed by Archive
class Batch {
private:
	//do not use the default constructor
	Batch() {};
public:

	//c'stor
	Batch( utility::options::OptionCollection const& options, bool intermediate_structs, bool has_silent_in, core::Size nstruct )
	: nstruct_( nstruct ),
		intermediate_structs_( intermediate_structs ),
		has_silent_in_( has_silent_in ),
		has_finished_( false ),
		is_cancelled_( false ),
		allow_reading_cancelled_decoys_( true ),
		invalid_( false ),
		options_( options ),
		decoys_returned_to_archive_( 0 )
	{};

	///c'stor
	Batch( core::Size id ) : batch_id_( id ),
		nstruct_( 0 ),
		intermediate_structs_( false ),
		has_silent_in_( false ),
		has_finished_( false ),
		is_cancelled_( false ),
		allow_reading_cancelled_decoys_( true ),
		invalid_( false ),
		decoys_returned_to_archive_( 0 )

	{};

	///some useful file- and directory names
	std::string batch() const;
	std::string dir() const;

	std::string silent_in() const; // Input
	std::string flag_file() const;
	std::string broker_file() const; // the broker file written by the Archive

	// extra broker files that might have been supplied by cmd-line
	std::string extra_broker_files() const; // borker files deliminiated by " "
	// extra broker files that might have been supplied by cmd-line

	std::string all_broker_files() const {
		return broker_file() + " " + extra_broker_files(); // borker files deliminiated by " "
	}

	std::string silent_out() const; // Output
	std::string alternative_decoys_out() const;

	std::string score_file() const;

	// the options defined in this batch, that are not controlled directly by Batch, as, e.g., silent-in/out , or nstruct
	utility::options::OptionCollection const& user_options() const { return options_; };
	utility::options::OptionCollection& user_options() { return options_; };

	///Getters

	/// has input decoys
	bool has_silent_in() const { return has_silent_in_; }

	/// writes out intermediate decoys
	bool intermediate_structs() const { return intermediate_structs_; };

	/// batch has finished
	bool has_finished() const { return has_finished_; }

	/// batch has finished
	bool is_cancelled() const { return is_cancelled_; }

	/// batch has finished
	bool allow_reading_cancelled_decoys() const { return allow_reading_cancelled_decoys_; }

	/// nstruct ...
	core::Size& nstruct() { return nstruct_; }

	core::Size nstruct() const { return nstruct_; }

	/// batch id
	core::Size id() const { return batch_id_; }

	/// how many structures have been processed by archive already
	core::Size decoys_returned() const { return decoys_returned_to_archive_; }

	///Setters
	void set_has_silent_in( bool setting = true ) {
		has_silent_in_ = setting;
	}

	void set_intermediate_structs( bool setting = true ) {
		intermediate_structs_ = setting;
	}

	void mark_as_finished() {
		has_finished_ = true;
	}

	void mark_as_cancelled( bool allow_reading_of_decoys = true ) {
		is_cancelled_ = true;
		allow_reading_cancelled_decoys_ = allow_reading_cancelled_decoys_ && allow_reading_of_decoys;
	}

	void mark_as_invalid() {
		invalid_ = true;
	}

	bool valid() const {
		return !invalid_;
	}

	void set_id( core::Size id ) {
		batch_id_ = id;
	}

	void set_decoys_returned( core::Size setting ) {
		decoys_returned_to_archive_ = setting;
	}

	//input-output: batch_id, nstruct, has_finished, decoys_returned...
	void show( std::ostream&, bool single_line = false ) const;
	friend std::ostream& operator<< (std::ostream&, Batch const& );
	friend std::istream& operator>> (std::istream&, Batch & );

	/// read and write BATCH_INFO ( decoys_returned/ finished etc..  )
	void write_info_file() const;
	void read_info_file();
private:
	core::Size batch_id_;
	core::Size nstruct_;
	bool intermediate_structs_;
	bool has_silent_in_;
	bool has_finished_;
	bool is_cancelled_;
	bool allow_reading_cancelled_decoys_;
	bool invalid_;
	utility::options::OptionCollection options_;
	core::Size decoys_returned_to_archive_;
	std::string rosetta_script_file_;
};


/// @brief ArchiveManager is responsible for communication with JobDistributor and organization of Batches and returning decoys
/// @detail he owns an Archive (AbstractArchiveBase) that will be handed the decoys and is asked to generate_batch() if the QUEUE_EMPTY .
class BaseArchiveManager {
public:
	/// @brief ctor is protected; singleton pattern
	BaseArchiveManager() : theArchive_( /* NULL */ ) {};
	virtual ~BaseArchiveManager() {};  //virtual destructor because we have virtual functions

public:
	typedef utility::vector1< Batch > BatchList;
	//static void register_options();

	virtual Batch& start_new_batch();
	virtual void finalize_batch( Batch&, bool reread = false );
	virtual void cancel_batches_previous_to( core::Size batch_id, bool allow_reading_of_decoys = true );

	virtual void save_archive() = 0;
	BatchList const& batches() const {
		return batches_;
	}
	virtual bool restore_archive() = 0;

	core::Size last_batch_id() const {
		return batches_.size();
	}

protected:
	virtual void queue_batch( Batch const& batch ) = 0;
	virtual void cancel_batch( Batch& batch, bool allow_reading_of_decoys = true );

	void set_archive( AbstractArchiveBaseOP );
	AbstractArchiveBase& the_archive() { return *theArchive_; }
	virtual void unlock_file( Batch const&, bool ) {};
	void read_returning_decoys( Batch& batch, bool final );

	static bool options_registered_;

	utility::vector1< Batch > batches_;

private:
	AbstractArchiveBaseOP theArchive_;
};

/// @brief ArchiveManager is responsible for communication with JobDistributor and organization of Batches and returning decoys
/// @detail he owns an Archive (AbstractArchiveBase) that will be handed the decoys and is asked to generate_batch() if the QUEUE_EMPTY .
class ArchiveManager : public BaseArchiveManager {
	typedef BaseArchiveManager Parent;
public:
	/// @brief ctor is protected; singleton pattern
	ArchiveManager( core::Size archive_rank, core::Size jd_master_rank, core::Size file_buf_rank );
	virtual ~ArchiveManager() {}; //virtual destructor because we have virtual functions


public:
	static void register_options();

	void go( ArchiveBaseOP );

	//this will read options to check
	//it will read the broker_file and check
	//register and queue the option...


	core::Size unfinished_batches() const;

	void save_archive();
	virtual bool restore_archive();

protected:
	/// @brief triggered in slave if new batch_ID comes in.
	virtual void batch_underflow() {};

	//  /// @brief return true if message was understood
	//  virtual bool process_message( int msg_tag, int slave_rank, int slave_job_id );

	void idle();

	void jobs_completed();// core::Size batch_id, bool final, core::Size bad );
	void queue_batch( Batch const& batch );
	void cancel_batch( Batch& batch, bool allow_reading_of_decoys = true );
	void read_existing_batches();
	void register_batch( Batch new_batch );
	void send_stop_to_jobdistributor();


	friend class JobDistributorFactory; //ctor access

	virtual void unlock_file( Batch const& batch, bool final );
private:


	core::Size const archive_rank_;
	core::Size const jd_master_rank_;
	core::Size const file_buf_rank_;

	typedef MPIArchiveJobDistributor::CompletionMessage CompletionMessage;
	typedef std::map< core::Size, CompletionMessage > CompletionMessages;
	CompletionMessages jobs_completed_;
	/// @brief specify seconds between automatic saves to filesystem
	core::Size save_archive_time_interval_;

	static bool options_registered_;
};

}//archive
}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_ArchiveManager_HH
