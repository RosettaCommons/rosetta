// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SigmoidFilter.cc
/// @brief
/// @author Gabi Pszolla & Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/SigmoidFilter.hh>
#include <protocols/simple_filters/SigmoidFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <fstream>
#include <utility/io/izstream.hh>
#include <sstream>

namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.Sigmoid" );

protocols::filters::FilterOP
SigmoidFilterCreator::create_filter() const { return protocols::filters::FilterOP( new Sigmoid ); }

std::string
SigmoidFilterCreator::keyname() const { return "Sigmoid"; }

//default ctor
Sigmoid::Sigmoid() :
	protocols::filters::Filter( "Sigmoid" ),
	filter_( /* NULL */ ),
	steepness_( 1.0 ),
	offset_( 0.0 ),
	baseline_( 0.0 ),
	negate_( false ),
	threshold_( 0 ),
	baseline_checkpointing_filename_( "" )
{
}

Sigmoid::~Sigmoid() {}

/// @brief The first MC trajectory should not read the baseline from the checkpoint file, instead, it should set the baseline
/// attempt_read_from_checkpoint determines whether a read attempt from checkpoint should be attempted
void
Sigmoid::reset_baseline( core::pose::Pose const & pose, bool const attempt_read_from_checkpoint ){
	using namespace std;

	//bool compute_new_baseline( false );
	if ( attempt_read_from_checkpoint ) {
		TR<<"Reading baseline from checkpoint file, if one exists"<<std::endl;
	} else {
		TR<<"Not reading from checkpoint file"<<std::endl;
	}
	if ( attempt_read_from_checkpoint && baseline_checkpointing_filename_ != "" ) {
		ifstream f( baseline_checkpointing_filename_.c_str(), ios::in );
		bool compute_new_baseline( false );
		if ( !f.good() ) {
			compute_new_baseline = true;
		} else {
			core::Size const begin = f.tellg();
			f.seekg( 0, ios::end );
			core::Size const end = f.tellg();
			f.seekg( 0, ios::beg );
			if ( end - begin == 0 ) { //file size == 0
				compute_new_baseline = true;
			}
		}
		if ( !compute_new_baseline ) {
			std::string line;
			getline( f, line );
			std::istringstream line_stream( line );
			line_stream >> baseline_;
			TR<<"Loading Sigmoid baseline from checkpoint. Loaded baseline: "<<baseline_<<std::endl;
			f.close();
			return;
		}
		f.close();
	}

	baseline_ = filter()->report_sm( pose );
	TR<<"Computed new baseline and set to: "<<baseline_<<std::endl;
	if ( baseline_checkpointing_filename_ != "" ) {
		ofstream f;
		f.open( baseline_checkpointing_filename_.c_str(), ios::out );
		if ( !f.good() ) {
			utility_exit_with_message( "Unable to open Sigmoid checkpointing file: " + baseline_checkpointing_filename_ );
		}
		f << baseline_;
		f.close();
		TR<<"Wrote baseline "<<baseline_<<" to checkpointing file "<<baseline_checkpointing_filename_<<std::endl;
	}//fi baseline_checkpointing_filename_
}

void
Sigmoid::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & )
{
	steepness( tag->getOption< core::Real >( "steepness", 1.0 ) );
	offset( tag->getOption< core::Real >( "offset", 0 ));
	negate( tag->getOption< bool >( "negate", false ) );
	threshold( tag->getOption< core::Real >( "threshold", 0 ) );
	if ( tag->hasOption( "filter" ) ) {
		filter( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "filter" ),filters  ) );
		TR<<" filter: "<<tag->getOption< std::string >( "filter" )<<std::endl;
	} else {
		TR<<"Filter not defined. I'm expecting another filter/mover to set my filter, else I will crash!!!"<<std::endl;
	}
	baseline_checkpointing_filename_ = tag->getOption< std::string >( "baseline_checkpoint", "" );
	TR<<"Sigmoid with options: steepness "<<steepness()<<" offset "<<offset()<<" negate "<<negate()<<" threshold "<<threshold()<<" baseline checkpointing file: "<<baseline_checkpointing_filename_<<std::endl;
}

bool
Sigmoid::apply( core::pose::Pose const & pose ) const {
	core::Real const val ( compute( pose ) );
	return( val >= threshold() );
}

void
Sigmoid::report( std::ostream &o, core::pose::Pose const & pose ) const {
	core::Real const val = compute( pose );
	o << "Sigmoid returns "<<val<<std::endl;
}

core::Real
Sigmoid::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
Sigmoid::compute(
	core::pose::Pose const & pose
) const {
	runtime_assert( filter().get() ); /// if I'm null the filter has not been set
	core::Real const val( filter()->report_sm( pose ) - baseline_ );
	core::Real const transform( 1.0 / ( ( 1.0 + std::exp( ( val - offset_ ) * steepness_ ) ) ) );
	core::Real const complement( negate() ? 1.0 - transform : transform ); // negate means to take the complement of the transform
	TR<<"filter val/transform: "<<val<<" "<<complement<<std::endl;
	TR<<"returning: "<<complement<<std::endl;
	return( complement );
}

protocols::filters::FilterOP
Sigmoid::filter() const{ return filter_; }

void
Sigmoid::filter( protocols::filters::FilterOP f ){ filter_ = f; }
}
}
