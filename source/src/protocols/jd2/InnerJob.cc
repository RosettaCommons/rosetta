// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/Job.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Job classes
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/InnerJob.hh>

///Project headers


#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
///Utility headers
#include <basic/Tracer.hh>

///C++ headers
#include <string>

#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.InnerJob" );

namespace protocols {
namespace jd2 {

/////////////////////////////InnerJob/////////////////////////////
InnerJob::InnerJob( std::string  input_tag, core::Size nstruct_max ) :
	input_tag_(std::move(input_tag)),
	nstruct_max_(nstruct_max),
	pose_(/* NULL */),
	bad_( false )
{
	//commented these out... somehow they don't respond to -mute
	//TR.Trace << "Using InnerJob (base class) for JobDistributor" << std::endl;
}

InnerJob::InnerJob( core::pose::PoseCOP pose, std::string  input_tag, core::Size nstruct_max ) :
	input_tag_(std::move(input_tag)),
	nstruct_max_(nstruct_max),
	pose_(std::move( pose )),
	bad_( false )
{
	//TR.Trace << "Using InnerJob (base class) for JobDistributor" << std::endl;
}

/// @brief Copy constructor.
///
InnerJob::InnerJob( InnerJob const &src ) :
	utility::pointer::ReferenceCount(),
	input_tag_(src.input_tag_),
	nstruct_max_(src.nstruct_max_),
	pose_(src.pose_->clone()),
	bad_( src.bad_ )
{
}

InnerJob::~InnerJob()= default;

/// @brief Return an owning pointer to a copy of this object.
///
InnerJobOP InnerJob::clone() const {
	return InnerJobOP( new InnerJob(*this) );
}

std::string const & InnerJob::input_tag() const { return input_tag_; }

core::Size InnerJob::nstruct_max() const { return nstruct_max_; }

core::pose::PoseCOP InnerJob::get_pose() const { return pose_; }

void InnerJob::set_pose( core::pose::PoseCOP pose ) { pose_ = pose; }

/// @details Only compare the pointer value of the pose the inner
///job is referencing.
bool
InnerJob::operator == (
	InnerJob const & other
) const {
	return same( other ) && other.same( *this );
}


bool
InnerJob::operator!=(
	InnerJob const & other
) const {
	return !(*this == other);
}

bool
InnerJob::same( InnerJob const & other ) const
{
	return !( input_tag_.compare( other.input_tag_) ) &&
		nstruct_max_ == other.nstruct_max_ &&
		pose_.get() == other.pose_.get() &&
		bad_ == other.bad_;
}


void
InnerJob::show(
	std::ostream & out
) const {
	out
		<< "input_tag:" << input_tag_ << std::endl
		<< "nstruct max: " << nstruct_max_ << std::endl
		<< "bad: " << (bad_ ? "true" : "false" ) << std::endl
		<< "Pose:";
	if ( pose_ ) {
		out << std::endl << *pose_ << std::endl;
	} else {
		out << " NULL" << std::endl;
	}
}

std::ostream &
operator<< (
	std::ostream & out,
	const InnerJob & inner_job
) {
	inner_job.show(out);
	return out;
}


} // jd2
} // protocols
