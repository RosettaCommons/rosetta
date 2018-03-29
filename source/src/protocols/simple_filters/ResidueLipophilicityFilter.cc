// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueLipophilicityFilter.cc
/// @brief returns the value calculated by the energy term mp_res_lipo.
/// useful for understanding and debugging this energy term. use the output_file tag
/// to print a full table of parameters from within the energy term
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
#include <protocols/simple_filters/ResidueLipophilicityFilter.hh>
#include <protocols/simple_filters/ResidueLipophilicityFilterCreator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <iostream>
#include <fstream>
#include <core/scoring/membrane/MPResidueLipophilicityEnergy.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.ResidueLipophilicityFilter" );

protocols::filters::FilterOP
ResidueLipophilicityFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ResidueLipophilicityFilter ); }

std::string
ResidueLipophilicityFilterCreator::keyname() const { return "ResidueLipophilicity"; }

ResidueLipophilicityFilter::~ResidueLipophilicityFilter(){}

void
ResidueLipophilicityFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	threshold_ = tag->getOption< core::Real >( "threshold", -10.0 );
	output_ = tag->getOption< std::string >( "output_file", "TR" );
	print_splines_ = tag->getOption< bool >( "print_splines", false );

	TR << "testing ResidueLipophilicity with threshold " << threshold_ << std::endl;
}

bool
ResidueLipophilicityFilter::apply( core::pose::Pose const & pose ) const {

	core::Real total_mp_res_lipo = compute( pose );

	bool const status = total_mp_res_lipo <= threshold_;

	TR << "ResidueLipophilicity threshold is at " << threshold_ << " total ResidueLipophilicity is " << total_mp_res_lipo << std::endl;

	return status;
}

void
ResidueLipophilicityFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real total_mp_res_lipo = compute( pose );
	out << "total res solv: " << total_mp_res_lipo << std::endl;
}

core::Real
ResidueLipophilicityFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real total_mp_res_lipo = compute( pose );
	return( total_mp_res_lipo );
}

core::Real
ResidueLipophilicityFilter::compute( core::pose::Pose const & pose ) const {
	core::scoring::membrane::MPResidueLipophilicityEnergy mp_mp_res_lipo = core::scoring::membrane::MPResidueLipophilicityEnergy();
	core::Real total_mp_res_lipo = 0.0;

	if ( output_ == "TR" ) {
		total_mp_res_lipo = mp_mp_res_lipo.report_ressolv( TR, pose, print_splines_ );
	} else if ( output_ == "auto" ) {
		std::ofstream outfile;
		std::string file_name = protocols::jd2::JobDistributor::get_instance()->current_output_name() + "_RS.txt";
		TR << "printing report to " << file_name << std::endl;
		outfile.open( file_name, std::ios::out  );
		total_mp_res_lipo = mp_mp_res_lipo.report_ressolv( outfile, pose, print_splines_ );
	} else {
		std::ofstream outfile;
		outfile.open( output_.c_str(), std::ios::out );
		total_mp_res_lipo = mp_mp_res_lipo.report_ressolv( outfile, pose, print_splines_ );
	}
	return( total_mp_res_lipo );
}

void ResidueLipophilicityFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "if ResidueLipophilicity score is lower than threshold, than pass." , "-10" )
		+ XMLSchemaAttribute::attribute_w_default( "output_file" , xs_string , "where to print the result to. default to TRACER. use auto for job_name_RS.txt" , "TR" )
		+ XMLSchemaAttribute::attribute_w_default( "print_splines" , xsct_rosetta_bool , "whehter to print the splines, default=false" , "false" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the overall ResidueLipophilicity of the pose.", attlist );
}

void ResidueLipophilicityFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueLipophilicityFilter::provide_xml_schema( xsd );
}

std::string ResidueLipophilicityFilter::name() const {
	return class_name();
}

std::string ResidueLipophilicityFilter::class_name() {
	return "ResidueLipophilicity";
}

}
}
