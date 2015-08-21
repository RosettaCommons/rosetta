// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/jd2/ShuffleJobDistributor.hh>
#include <protocols/jd2/FileSystemJobDistributor.hh>
// Package headers
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>


// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>
#include <numeric/random/random.hh>

// Option headers
#include <basic/options/keys/out.OptionKeys.gen.hh>


// C++ headers
#include <string>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.jd2.ShuffleFileSystemJobDistributor" );


namespace protocols {
namespace jd2 {

/// @details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land.
ShuffleFileSystemJobDistributor::ShuffleFileSystemJobDistributor() :
	FileSystemJobDistributor()
{
	next_random_job_ = 0; // indicate that list needs generating

}

/// @brief dtor
///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
ShuffleFileSystemJobDistributor::~ShuffleFileSystemJobDistributor()
{ }


/// @brief dummy for master/slave version
core::Size
ShuffleFileSystemJobDistributor::get_new_job_id()
{
	JobsContainer const & jobs( get_jobs() );
	if ( jobs.size() == 0 ) return 0; // no jobs;

	JobOutputterOP outputter = job_outputter();

	if ( next_random_job_ == 0 ) {
		scrambled_job_order_.clear();
		for ( int i=0; i < basic::options::option[ basic::options::OptionKeys::out::nstruct ].value() ; i++ ) {
			scrambled_job_order_.push_back( numeric::random::rg().random_range(1,jobs.size() ) );
		}


	}
	next_random_job_++;

	// for( core::Size ijob = 1; ijob <= jobs.size(); ijob ++ ){
	//  std::cerr  << "JL: " << jobs[ ijob ]->input_tag() << "  " << jobs[ ijob ]->nstruct_index()  << " Done? " << outputter->job_has_completed( jobs[ ijob ] ) << std::endl;
	// }
	//
	core::Size choice = scrambled_job_order_[next_random_job_];
	if ( (int)next_random_job_ > (int)basic::options::option[ basic::options::OptionKeys::out::nstruct ].value() ) return 0; // we're done (nstruct decoys done)
	if ( next_random_job_ > scrambled_job_order_.size() ) next_random_job_ = 0; // re-scramble the list

	// Now go backwards to find the last nstruct that hasn't been done yet.
	// If we're not done
	if ( choice <= 0 ) choice = 1;
	if ( choice > jobs.size() ) choice = jobs.size();
	if ( !outputter->job_has_completed( jobs[ choice ] ) ) {
		while ( ( choice > 1 )  &&
				( !outputter->job_has_completed( jobs[ choice-1 ]) ) &&
				( jobs[ choice - 1 ]->input_tag() == jobs[ choice ]->input_tag() )
				) {
			choice -- ;
		}

	}

	if ( choice <= 0 ) choice = 1;
	if ( choice > jobs.size() ) choice = jobs.size();
	if ( outputter->job_has_completed( jobs[ choice ] ) ) {
		while ( ( choice <= (jobs.size()-1) )  &&
				( outputter->job_has_completed( jobs[ choice ]) )
				) {
			choice ++ ;
		}
	}

	// still no luck ? then go back wards again looking for anything non-completed.
	if ( choice <= 0 ) choice = 1;
	if ( choice > jobs.size() ) choice = jobs.size();
	if ( outputter->job_has_completed( jobs[ choice ] ) ) {
		while ( ( choice > 1 )  &&
				( outputter->job_has_completed( jobs[ choice ]) )
				) {
			choice -- ;
		}
	}

	// if at 0 and still completed status is true, we must have done all the jobs - return 0;
	if ( choice <= 0 ) choice = 1;
	if ( choice > jobs.size() ) choice = jobs.size();
	if ( outputter->job_has_completed( jobs[ choice ] ) ) {
		TR << "No more jobs found." << std::endl;
		return 0;
	}

	// Now go backwards to find the last nstruct that hasn't been done yet.
	// If we're not done
	if ( choice <= 0 ) choice = 1;
	if ( choice > jobs.size() ) choice = jobs.size();
	if ( !outputter->job_has_completed( jobs[ choice ] ) ) {
		while ( ( choice > 1 )  &&
				( !outputter->job_has_completed( jobs[ choice-1 ]) ) &&
				( jobs[ choice - 1 ]->input_tag() == jobs[ choice ]->input_tag() )
				) {
			choice -- ;
		}

	}


	return choice;
}

void
ShuffleFileSystemJobDistributor::mark_current_job_id_for_repetition()
{
	// do nothing - no repetitions allowed. Behave as if job had succeeded (handled by job_failed)
}

}//jd2
}//protocols
