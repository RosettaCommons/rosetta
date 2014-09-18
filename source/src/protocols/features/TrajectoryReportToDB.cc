// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/features/TrajectoryReportToDB.cc
///
/// @brief  Report features data to database multiple times per structure, creating a trajectory
/// @author Kyle Barlow (kb@kylebarlow.com)

// #ifdef USEMPI
// #include <mpi.h>
// #endif

#include <protocols/features/TrajectoryReportToDB.hh>
#include <protocols/features/TrajectoryMapFeatures.hh>

// Setup Mover
#include <protocols/features/TrajectoryReportToDBCreator.hh>
#include <basic/database/sql_utils.hh>

// Platform Headers
#include <basic/Tracer.hh>

namespace protocols{
namespace features{

// Constructors

TrajectoryReportToDB::TrajectoryReportToDB() :
	ReportToDB("TrajectoryReportToDB"),
	stride_(1),
	trajectory_map_features_reporter_()
{
	initialize_trajectory_reporter();
	ReportToDB::set_batch_name("trajectory_features");
}

TrajectoryReportToDB::TrajectoryReportToDB(
	core::Size stride
) :
	ReportToDB("TrajectoryReportToDB"),
	stride_(stride),
	trajectory_map_features_reporter_()
{
	initialize_trajectory_reporter();
	ReportToDB::set_batch_name("trajectory_features");
}

TrajectoryReportToDB::TrajectoryReportToDB(
	std::string const & name
) :
	ReportToDB(name),
	stride_(1),
	trajectory_map_features_reporter_()
{
	initialize_trajectory_reporter();
	ReportToDB::set_batch_name("trajectory_features");
}

TrajectoryReportToDB::TrajectoryReportToDB(
	utility::sql_database::sessionOP db_session,
	std::string const & batch_name,
	std::string const & batch_description,
	bool use_transactions,
	core::Size cache_size
) :
	ReportToDB(
		"TrajectoryReportToDB",
		db_session, batch_name, batch_description,
		use_transactions, cache_size ),
	stride_(1),
	trajectory_map_features_reporter_()
{
	initialize_trajectory_reporter();
}

TrajectoryReportToDB::TrajectoryReportToDB(
	TrajectoryReportToDB const & src ) :
	ReportToDB(src),
	stride_(src.stride_),
	trajectory_map_features_reporter_(src.trajectory_map_features_reporter_)
{}

TrajectoryReportToDB::~TrajectoryReportToDB(){}

// Functions

void
TrajectoryReportToDB::set_stride(
	core::Size stride
) {
	stride_ = stride;
}

core::Size
TrajectoryReportToDB::get_stride () const {
	return stride_;
}

std::map<std::string, core::Size>
TrajectoryReportToDB::get_cycle_counts() const {
	return cycle_counts_;
}

std::string
TrajectoryReportToDBCreator::keyname() const
{
	return TrajectoryReportToDBCreator::mover_name();
}

moves::MoverOP
TrajectoryReportToDBCreator::create_mover() const {
	return new TrajectoryReportToDB;
}

std::string
TrajectoryReportToDBCreator::mover_name()
{
	return "TrajectoryReportToDB";
}

static thread_local basic::Tracer TR( "protocols.features.TrajectoryReportToDB" );

moves::MoverOP
TrajectoryReportToDB::fresh_instance() const { return new TrajectoryReportToDB; }

moves::MoverOP
TrajectoryReportToDB::clone() const
{
	return new TrajectoryReportToDB( *this );
}

void
TrajectoryReportToDB::initialize_trajectory_reporter()
{
	trajectory_map_features_reporter_ = new TrajectoryMapFeatures();
	ReportToDB::add_features_reporter( trajectory_map_features_reporter_ );
}

void
TrajectoryReportToDB::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & pose )
{
	ReportToDB::parse_my_tag(
		tag, data, filters, movers, pose
	);

	parse_stride_tag_item(tag);
}

void
TrajectoryReportToDB::parse_stride_tag_item(
	TagCOP const tag) {
	if(tag->hasOption("stride")){
	  set_stride( tag->getOption<core::Size>("stride") );
	}
}

void
TrajectoryReportToDB::apply( Pose& pose )
{
	std::map<std::string, core::Size>::iterator tag_iter;
	std::string structure_tag;
	core::Size cycle_count;

	ReportToDB::ensure_structure_tags_are_ready();
	structure_tag = ReportToDB::get_structure_tag();

	tag_iter = cycle_counts_.find( structure_tag );
	if ( tag_iter == cycle_counts_.end() ) {
		// New job output tag - initialize cycle count at 0
		cycle_counts_[structure_tag] = 0;
		cycle_count = 0;
	}
	else {
		cycle_count = tag_iter->second;
	}

	trajectory_map_features_reporter_->set_current_cycle( cycle_count );

	if ( cycle_count % stride_ == 0 ) {
		ReportToDB::apply( pose );
	}

	// Increase cycle count
	cycle_counts_[structure_tag] += 1;
}

} // namespace
} // namespace
