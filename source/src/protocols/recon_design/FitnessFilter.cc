// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/recon_design/FitnessFilter.cc
/// @brief Returns the sum of energy of input poses. Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#include <protocols/recon_design/FitnessFilter.hh>
#include <protocols/recon_design/FitnessFilterCreator.hh>

#include <protocols/filters/VectorPoseFilter.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/Energies.hh>
#include <utility/tag/Tag.hh>
#include <utility/mpi_util.hh>
#include <basic/datacache/DataMap.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace recon_design {

static basic::Tracer TR( "protocols.recon_design.FitnessFilter" );

filters::FilterOP
FitnessFilterCreator::create_filter() const {
	return FitnessFilterOP( new FitnessFilter );
}

std::string
FitnessFilterCreator::keyname() const {
	return "FitnessFilter";
}


FitnessFilter::FitnessFilter() :
	filters::VectorPoseFilter( "FitnessFilter" ),
	sfxn_( core::scoring::get_score_function() ),
	output_to_scorefile_( false ),
	threshold_( std::numeric_limits<float>::max() )
{}

FitnessFilter::~FitnessFilter() {}

filters::FilterOP FitnessFilter::clone() const {
	return filters::FilterOP( new FitnessFilter( *this ) );
}

filters::FilterOP FitnessFilter::fresh_instance() const {
	return filters::FilterOP( new FitnessFilter() );
}


void
FitnessFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const & ) {

	if ( tag->hasOption("scorefxn") ) {
		std::string const scorefxn_key( tag->getOption<std::string>("scorefxn") );
		if ( datamap.has( "scorefxns", scorefxn_key ) ) {
			sfxn_ = datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_key );
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,"ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
		}
	}

	output_to_scorefile_ = tag->getOption<bool>( "output_to_scorefile", false );
	threshold_= tag->getOption<core::Real>( "threshold", std::numeric_limits<float>::max() );
}

void
FitnessFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "scorefxn", xs_string, "Score function to use when evaluating fitness", "talaris2014" )
		+ XMLSchemaAttribute::attribute_w_default( "output_to_scorefile", xsct_rosetta_bool, "Output the fitness to score file?", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Threshold to pass or fail structure based on fitness", std::to_string(std::numeric_limits<float>::max() ));

	protocols::filters::xsd_type_definition_w_attributes( xsd, "FitnessFilter", "Calculates the fitness over all input poses. To be used with the VectorPoseMover job distributor. Fitness is defined as the sum of energy of all input states", attlist );
}

void FitnessFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FitnessFilter::provide_xml_schema( xsd );
}

bool
FitnessFilter::apply( core::pose::Pose const & pose ) const {

	// Add scores to fitness
	core::Real fitness = calculate_fitness( pose );

	if ( output_to_scorefile_ ) {
		protocols::jd2::JobOP current_job = protocols::jd2::JobDistributor::get_instance()->current_job();
		current_job->add_string_real_pair( "fitness", fitness );
	}

	return fitness < threshold_;
}

bool
FitnessFilter::apply_mpi( core::pose::Pose const & pose ) const {

	core::Real overall_fitness = calculate_fitness_mpi( pose );

	if ( output_to_scorefile_ ) {
		protocols::jd2::JobOP current_job = protocols::jd2::JobDistributor::get_instance()->current_job();
		current_job->add_string_real_pair( "fitness", overall_fitness );
	}

	return overall_fitness < threshold_;
}

/// NOTE: report and report_sm are not implemented for use with the VectorPoseJobDistributor and the recon application
void
FitnessFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {

	core::Real fitness = calculate_fitness( pose );

	out << "fitness: " << fitness << std::endl;
}

/// NOTE: report and report_sm are not implemented for use with the VectorPoseJobDistributor and the recon application
core::Real
FitnessFilter::report_sm( core::pose::Pose const & pose ) const {

	return calculate_fitness( pose );
}

core::Real FitnessFilter::calculate_fitness( core::pose::Pose const & ) const {
	core::Real fitness = 0;
	for ( core::Size i = 1; i <= poses_.size(); ++i ) {
		core::pose::PoseOP pose_copy = poses_[ i ]->clone();
		core::Real pose_score = sfxn_->score( *pose_copy );

		fitness += pose_score;
	}
	return fitness;
}

core::Real FitnessFilter::calculate_fitness_mpi( core::pose::Pose const & pose ) const {

	core::pose::PoseOP pose_copy = pose.clone();
	core::Size n_procs = utility::mpi_nprocs();
	core::Size rank = utility::mpi_rank();

	// Step 1: calculate my pose's fitness
	core::Real pose_score = sfxn_->score( *pose_copy );


	// Step 2: have the master compile all scores into one fitness
	core::Real overall_fitness;
	if ( rank == 0 ) {
		overall_fitness = pose_score;
		for ( core::Size rank_to_receive = 1; rank_to_receive < n_procs; ++rank_to_receive ) {
			core::Real other_pose_fitness = utility::receive_double_from_node( rank_to_receive );
			overall_fitness += other_pose_fitness;
		}

	} else {
		utility::send_double_to_node( 0, pose_score );
	}


	// Step 3: now that I know my overall fitness send to all the other nodes
	if ( rank == 0 ) {
		for ( core::Size rank_to_send = 1; rank_to_send < n_procs; ++rank_to_send ) {
			utility::send_double_to_node( rank_to_send, overall_fitness );
		}
	} else {
		overall_fitness = utility::receive_double_from_node( 0 );
	}

	return overall_fitness;
}

core::Real FitnessFilter::threshold() const { return threshold_; }
void FitnessFilter::threshold( core::Real thresh ) { threshold_ = thresh; }

core::scoring::ScoreFunctionOP FitnessFilter::sfxn() const { return sfxn_; }
void FitnessFilter::sfxn( core::scoring::ScoreFunctionOP sfxn) { sfxn_ = sfxn; }

} //recon_design
} //protocols
