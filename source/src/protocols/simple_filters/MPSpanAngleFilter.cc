// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/MPSpanAngleFilter.cc
/// @brief calculates the angle between the TM span and the membrane normal
/// @details naturla TM span angles are generally distributed between 0-60 degrees. this filter claculates the angle
/// of the speicifed TM span and either returns the angle, or a score corresponding to how different that angle is
/// from the natural distribution
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)


// Project Headers

#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <map>
#include <numeric/random/random.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_filters/MPSpanAngleFilter.hh>
#include <protocols/simple_filters/MPSpanAngleFilterCreator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <iostream>
#include <fstream>
#include <core/scoring/membrane/MPSpanAngleEnergy.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.MPSpanAngleFilter" );

protocols::filters::FilterOP
MPSpanAngleFilterCreator::create_filter() const { return protocols::filters::FilterOP( new MPSpanAngleFilter ); }

std::string
MPSpanAngleFilterCreator::keyname() const { return "MPSpanAngle"; }

MPSpanAngleFilter::~MPSpanAngleFilter(){}

void
MPSpanAngleFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	threshold_ = tag->getOption< core::Real >( "threshold", -10.0 );
	output_ = tag->getOption< std::string >( "output_file", "TR" );
	tm_num_ = tag->getOption< core::Size >( "tm", 0 );
	angle_or_score_ = tag->getOption< std::string >( "angle_or_score", "angle" );
	ang_max_ = tag->getOption< core::Real >( "ang_max", 90 );
	ang_min_ = tag->getOption< core::Real >( "ang_min", 0 );

	TR << "testing MPSpanAngle with threshold " << threshold_ << std::endl;
}

bool
MPSpanAngleFilter::apply( core::pose::Pose const & pose ) const {

	utility::vector1< core::Real > total_mp_span_angle = compute( pose );
	if ( total_mp_span_angle.size() < tm_num_ ) {
		TR.Warning << "not enough angles for TM num " << tm_num_ << " returning 0" << std::endl;
		return( 0 );
	}

	bool status = 1;
	bool temp_status = 0;
	if ( angle_or_score_ == "score" ) {
		for ( core::Size i=1; i<=total_mp_span_angle.size(); ++i ) {
			if ( tm_num_ == 0 || tm_num_ == i ) {
				temp_status = total_mp_span_angle[ i ] <= threshold_;
				status = status && temp_status ? true : false;
				TR << "MPSpanAngle threshold is at " << threshold_ << " angle for span " << i << " is " << total_mp_span_angle[ i ] << " scoring: " << total_mp_span_angle[ i ] << std::endl;
			}
		}
	} else {
		for ( core::Size i=1; i<=total_mp_span_angle.size(); ++i ) {
			if ( tm_num_ == 0 || tm_num_ == i ) {
				temp_status = total_mp_span_angle[ i ] <= ang_max_ && total_mp_span_angle[ i ] >= ang_min_;
				status = status && temp_status ? true : false;
				TR << "the angle for TM " << i << " is " << total_mp_span_angle[ i ] << " and its status is " << temp_status << " the overall status is " << status << std::endl;

			}
		}
	}
	return status;
}

void
MPSpanAngleFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {

	utility::vector1< core::Real > total_mp_span_angle = compute( pose );
	if ( total_mp_span_angle.size() < tm_num_ ) {
		TR.Warning << "not enough angles for TM num " << tm_num_ << " returning 0" << std::endl;
		return;
	}

	bool status = 1;
	bool temp_status = 0;
	if ( angle_or_score_ == "score" ) {
		for ( core::Size i=1; i<=total_mp_span_angle.size(); ++i ) {
			if ( tm_num_ == 0 || tm_num_ == i ) {
				temp_status = total_mp_span_angle[ i ] <= threshold_;
				status = status && temp_status ? true : false;
				out << "MPSpanAngle threshold is at " << threshold_ << " angle for span " << i << " is " << total_mp_span_angle[ i ] << " scoring: " << total_mp_span_angle[ i ] << std::endl;
			}
		}
	} else {
		for ( core::Size i=1; i<=total_mp_span_angle.size(); ++i ) {
			if ( tm_num_ == 0 || tm_num_ == i ) {
				temp_status = total_mp_span_angle[ i ] <= ang_max_ && total_mp_span_angle[ i ] >= ang_min_;
				status = status && temp_status ? true : false;
				out << "the angle for TM " << i << " is " << total_mp_span_angle[ i ] << " and its status is " << temp_status << " the overall status is " << status << std::endl;

			}
		}
	}
}

core::Real
MPSpanAngleFilter::report_sm( core::pose::Pose const & pose ) const {
	utility::vector1< core::Real > total_mp_span_angle = compute( pose );
	if ( total_mp_span_angle.size() < tm_num_ ) {
		TR.Warning << "not enough angles for TM num " << tm_num_ << " returning 0" << std::endl;
		return( 0 );
	}
	core::Real result = 0;
	bool status;
	utility::vector1< core::Real > results;
	if ( angle_or_score_ == "score" ) {
		for ( core::Size i=1; i<=total_mp_span_angle.size(); ++i ) {
			if ( tm_num_ == 0 || tm_num_ == i ) {
				status = total_mp_span_angle[ i ] <= threshold_;
				TR << "MPSpanAngle threshold is at " << threshold_ << " angle for span " << i << " is " << total_mp_span_angle[ i ] << " scoring: " << total_mp_span_angle[ i ] << " status: " << status << std::endl;
				results.push_back( total_mp_span_angle[ i ] );
			}
		}
	} else {
		for ( core::Size i=1; i<=total_mp_span_angle.size(); ++i ) {
			if ( tm_num_ == 0 || tm_num_ == i ) {
				status = total_mp_span_angle[ i ] <= ang_max_ && total_mp_span_angle[ i ] >= ang_min_;
				TR << "the angle for TM " << i << " is " << total_mp_span_angle[ i ] << " and the status is " << status << std::endl;
				results.push_back( total_mp_span_angle[ i ] );
			}
		}
	}
	for ( auto const & i : results ) result += i;
	return( result / results.size() );
}

utility::vector1 < core::Real >
MPSpanAngleFilter::compute( core::pose::Pose const & pose ) const {
	core::scoring::membrane::MPSpanAngleEnergy mp_span_angle = core::scoring::membrane::MPSpanAngleEnergy();
	utility::vector1< core::Real > total_mp_span_angle;

	if ( output_ == "TR" ) {
		total_mp_span_angle = mp_span_angle.compute( pose, TR, true );
	} else if ( output_ == "auto" ) {
		std::ofstream outfile;
		std::string file_name = protocols::jd2::JobDistributor::get_instance()->current_output_name() + "_MPSpanAngle.txt";
		TR << "printing report to " << file_name << std::endl;
		outfile.open( file_name, std::ios::out  );
		total_mp_span_angle = mp_span_angle.compute( pose, outfile, true );
	} else {
		std::ofstream outfile;
		outfile.open( output_.c_str(), std::ios::out );
		total_mp_span_angle = mp_span_angle.compute( pose, outfile, true );
	}
	return( total_mp_span_angle );
}

void MPSpanAngleFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "if MPSpanAngle score is lower than threshold, than pass." , "5" )
		+ XMLSchemaAttribute::attribute_w_default( "tm" , xsct_non_negative_integer, "which TM to calculate angle for. if zero will test all tms" , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "ang_max" , xsct_real , "maxiaml angle to pass" , "90" )
		+ XMLSchemaAttribute::attribute_w_default( "ang_min" , xsct_real , "minimal angle to pass" , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "output_file" , xs_string , "where to print the result to. default to TRACER. use auto for job_name_MPSpanAngle.txt" , "TR" )
		+ XMLSchemaAttribute::attribute_w_default( "angle_or_score" , xs_string , "whether to return the angle or the score. either angle or score" , "angle" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the angle (in degrees) with the membrane normal for each TM span, and its relevant score.", attlist );
}

void MPSpanAngleFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MPSpanAngleFilter::provide_xml_schema( xsd );
}

std::string MPSpanAngleFilter::name() const {
	return class_name();
}

std::string MPSpanAngleFilter::class_name() {
	return "MPSpanAngle";
}

}
}
