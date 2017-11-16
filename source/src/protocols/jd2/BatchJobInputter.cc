// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/BatchJobInputter.cc
/// @brief
/// @author Oliver Lange

///Unit headers
#include <protocols/jd2/BatchJobInputter.hh>

///Project headers
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>


///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

///C++ headers
#include <string>

// option key includes

static basic::Tracer tr( "protocols.jd2.BatchJobInputter" );

namespace protocols {
namespace jd2 {

std::string const BatchJobInputter::BOGUS_BATCH_ID( "NO_BATCH" );

BatchJobInputter::BatchJobInputter( std::string batch1 ) :
	current_batch_( batch1 ),
	vanilla_options_( basic::options::option )
{
	if ( batch1 != BOGUS_BATCH_ID ) {
		tr.Debug << "Instantiate BatchJobInputter with batch" << current_batch_ << std::endl;
		read_batch();
	} else {
		this_batch_job_inputter_ = JobDistributorFactory::create_job_inputter();
	}
}

BatchJobInputter::~BatchJobInputter() {
	basic::options::option=vanilla_options_;
}

/// @brief Return the type of input source that the BatchJobInputter is currently
///  using.
/// @return The input source for the current batch.
JobInputterInputSource::Enum BatchJobInputter::input_source() const {
	return this_batch_job_inputter_->input_source();
}

void BatchJobInputter::check_batch() {
	JobDistributor* jd( protocols::jd2::JobDistributor::get_instance() );
	if ( jd != nullptr && current_batch_ != jd->get_current_batch() ) {
		current_batch_ = jd->get_current_batch();
		read_batch();
	}
}

void BatchJobInputter::read_batch() {
	using namespace basic::options;
	option=vanilla_options_;
	option.load_options_from_file(current_batch_);
	//unfortunately I have to copy code from the JobDistributor factory unless I can remove "in:file:batch" from options...
	this_batch_job_inputter_ = JobDistributorFactory::create_job_inputter();
}

} // jd2
} // protocols
