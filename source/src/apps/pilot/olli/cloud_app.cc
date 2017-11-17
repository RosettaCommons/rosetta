// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange


#include <devel/init.hh>

#include <core/types.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>

#include <protocols/jd2/archive/ArchiveManager.hh>
#include <protocols/abinitio/IterativeAbrelax.hh>
#include <protocols/abinitio/BrokerMain.hh>

// Utility headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

static basic::Tracer tr( "main" );

using namespace core;
using namespace protocols;

//using namespace jumping;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( Boolean, new_batch )
OPT_KEY( Integer, read_batch )

class CloudArchiveManagerWrapper : public protocols::jd2::archive::BaseArchiveManager {
public:
	CloudArchiveManagerWrapper();
	virtual ~CloudArchiveManagerWrapper() {}; //virtual destructor because we have virtual functions
	typedef utility::vector1< protocols::jd2::archive::Batch > BatchList;
	// static void register_options();
	virtual void save_archive();
	virtual bool restore_archive();

	void read_batch( core::Size batch_id );
	void generate_batch();

protected:
	void queue_batch( protocols::jd2::archive::Batch const& batch ) {
		std::cout << "START BATCH" << batch << std::endl;
	};


	// static bool options_registered_;

	void save_batchlist();
	void restore_batchlist();

private:
	// protocols::jd2::archive::AbstractArchiveBaseOP theArchive_;
	// core::Size cancel_batches_previous_to_;
	// bool allow_reading_of_decoys_;
};

CloudArchiveManagerWrapper::CloudArchiveManagerWrapper() {
	set_archive( protocols::jd2::archive::AbstractArchiveBaseOP( new protocols::abinitio::IterativeAbrelax ) );
}

void CloudArchiveManagerWrapper::save_archive() {
	//write batches
	save_batchlist();
	the_archive().save_to_file();
}

bool CloudArchiveManagerWrapper::restore_archive() {
	restore_batchlist();
	return the_archive().restore_from_file();
}

void CloudArchiveManagerWrapper::save_batchlist() {
	std::string const BATCH_DB( "batch_db.txt" );
	utility::io::ozstream batch_info( BATCH_DB );
	for ( BatchList::const_iterator it = batches_.begin(); it != batches_.end(); ++it ) {
		batch_info << "BATCH " << *it << std::endl;
	}
}

void CloudArchiveManagerWrapper::restore_batchlist() {
	std::string const BATCH_DB( "batch_db.txt" );
	utility::io::izstream batch_info( BATCH_DB );
	std::string tag( "NONE" );
	while ( batch_info.good() ) {
		while ( batch_info.good() && tag != "BATCH" ) batch_info >> tag;
		if ( tag == "BATCH" ) {
			protocols::jd2::archive::Batch new_batch(0);
			batch_info >> new_batch;
			std::cerr << "just read batch " << new_batch << std::endl;
			batches_.push_back( new_batch );
			tag = "NONE";
		}
	}
}

void CloudArchiveManagerWrapper::read_batch( core::Size batch_id ) {
	if ( batch_id > batches_.size() ) {
		tr.Error << "invalid batch_id " << batch_id << " ! only " << batches_.size() << " have been registered " << std::endl;
		return ;
	}
	read_returning_decoys( batches_[ batch_id ], true /*final*/ );
}

void CloudArchiveManagerWrapper::generate_batch() {
	the_archive().generate_batch();
}

int main( int argc, char** argv ) {
	try{
		protocols::abinitio::IterativeAbrelax::register_options();
		protocols::jd2::archive::ArchiveManager::register_options();
		protocols::abinitio::register_options_broker();
		NEW_OPT( new_batch, "generate new batch", false );
		NEW_OPT( read_batch, "read decoys from batch", 0 );

		devel::init( argc, argv );

		using namespace basic::options;

		CloudArchiveManagerWrapper theArchiveManager;

		theArchiveManager.restore_archive();

		if ( option[ OptionKeys::read_batch ].user() ) {
			core::Size const batch_id = option[ OptionKeys::read_batch ]();
			theArchiveManager.read_batch( batch_id );
		}

		if ( option[ OptionKeys::new_batch ]() ) {
			theArchiveManager.generate_batch();
		}

		theArchiveManager.save_archive();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		std::exit( 1 );
	}
// finish with 0 when there's no error
	return 0;
}
