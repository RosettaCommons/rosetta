// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SpanTopologyMatchPoseFilter.cc
/// @brief return the matching (%) between span file determined topology and of the actual pose
/// @details this filter checks whether the topology descriobed by the use (i.e in the span file on tags)
/// is actually kept in the pose. it dose this by creating a topology string, made of i (in), o (out) and M
/// (membrnae) corresponding to each amino acid in the pose. such a topology string is made both for the user
/// specified topology, and the actual pose in the current conformation. the percentage of identity between
/// the two strings is reported.
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
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/SpanTopologyMatchPoseFilter.hh>
#include <protocols/simple_filters/SpanTopologyMatchPoseFilterCreator.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <fstream>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <protocols/membrane/util.hh>
#include <core/conformation/membrane/Exceptions.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.SpanTopologyMatchPoseFilter" );

protocols::filters::FilterOP
SpanTopologyMatchPoseFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SpanTopologyMatchPoseFilter ); }

std::string
SpanTopologyMatchPoseFilterCreator::keyname() const { return "SpanTopologyMatchPose"; }

SpanTopologyMatchPoseFilter::~SpanTopologyMatchPoseFilter(){}

void
SpanTopologyMatchPoseFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	threshold_ = tag->getOption< core::Real >( "threshold", 0.3 );
	output_ = tag->getOption< std::string >( "output_file", "" );

	TR << "testing SpanTopologyMatchPose with threshold " << threshold_ << std::endl;
}

bool
SpanTopologyMatchPoseFilter::apply( core::pose::Pose const & pose ) const {
	core::Real ts_match = compute( pose );
	bool const status = ts_match <= threshold_;
	return status;
}

void
SpanTopologyMatchPoseFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real ts_match = compute( pose );
	out << "calculated a topology match at " << ts_match << std::endl;
}

core::Real
SpanTopologyMatchPoseFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real ts_match = compute( pose );
	return( ts_match );
}

core::Real
SpanTopologyMatchPoseFilter::compute( core::pose::Pose const & pose ) const {
	std::string span_file_topo = span_file_topology( pose );
	std::string actual_topo = actual_topology( pose );

	TR << "topology string in span file " << span_file_topo << std::endl;
	TR << "topology string in real pose " << actual_topo << std::endl;

	core::Real score = compare_topo_strings( span_file_topo, actual_topo );
	TR << "the matching score is " << score << std::endl;
	return( score );
}

core::Real
SpanTopologyMatchPoseFilter::compare_topo_strings( std::string const ts_span, std::string const ts_actual ) const {
	core::Real total = 0.0;
	core::Real match = 0.0;

	for ( core::Size i=0; i < ts_span.size(); ++i ) {
		total += 1.0;
		// this is what it should be:
		if ( ts_span[ i ] == ts_actual[ i ] ) {
			match += 1.0;
		}
	}

	return( match / total );
}

std::string
SpanTopologyMatchPoseFilter::span_file_topology( core::pose::Pose const & pose ) const {
	std::string result = "";
	utility::vector1< utility::pointer::shared_ptr< core::conformation::membrane::Span > > const spans = pose.conformation().membrane_info()->spanning_topology()->get_spans();

	char current_in_out;
	char char_i = 'i';
	char char_o = 'o';
	if ( spans[ 1 ]->orientation() == 1 ) {
		current_in_out = char_i;
	} else {
		current_in_out = char_o;
	}

	int end_of_last = 1;
	for ( core::Size span_i=1; span_i <= spans.size(); span_i++ ) {
		if ( end_of_last > (int) spans[ span_i   ]->start() ) continue;
		result += std::string( (int) spans[ span_i ]->start() - (int) end_of_last, current_in_out);
		result += std::string( (int) spans[ span_i ]->end() - (int) spans[ span_i ]->start() + 1, 'M');
		end_of_last = spans[ span_i ]->end() + 1;
		if ( current_in_out == char_i ) {
			current_in_out = char_o;
		} else {
			current_in_out = char_i;
		}
	}

	// find total number of protein residues
	core::Size total_residues = 0;
	for ( core::Size i=1; i <= pose.size(); i++ ) {
		if ( ! pose.residue(i).is_protein() ) continue;
		total_residues += 1;
	}

	if ( !core::pose::symmetry::is_symmetric( pose ) ) {
		// if pose is not symmetric, just add the end..
		result += std::string( (int) total_residues - (int) spans[ spans.size() ]->end(), current_in_out);
	} else {
		// if pose symmettric create string for single component,
		// and muliply by the number of compoenents in the total pose
		core::conformation::symmetry::SymmetryInfoOP symm_info;
		core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
		core::Size const nres_component = sym_info->get_nres_subunit();

		result += std::string( (int) nres_component - (int) spans[ spans.size() ]->end(), current_in_out);

		core::Size const num_components = sym_info->get_num_components();
		for ( core::Size i=1; i <= num_components; i++ ) {
			result += result;
		}
	}
	return( result );
}

std::string
SpanTopologyMatchPoseFilter::actual_topology( core::pose::Pose const & pose ) const {
	std::string result = "";

	for ( core::Size res_i=1; res_i <= pose.size(); res_i++ ) {
		core::conformation::Residue const & residue_i = pose.residue(res_i);
		if ( ! residue_i.is_protein() ) continue;

		core::Real z = pose.conformation().membrane_info()->residue_z_position( pose.conformation(), residue_i.seqpos() );

		if ( z < -15 ) {
			result += 'i';
		} else if ( -15 < z && z < 15 ) {
			result += 'M';
		} else if ( 15 < z ) {
			result += 'o';
		} else {
			TR << "residue " << residue_i << " has a wierd z " << z << std::endl;
		}
	}

	return( result );
}

std::string SpanTopologyMatchPoseFilter::name() const {
	return class_name();
}

std::string SpanTopologyMatchPoseFilter::class_name() {
	return "SpanTopologyMatchPose";
}

void SpanTopologyMatchPoseFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "threhsold underneath to fail pose. in 0-1" , "0.3" )
		+ XMLSchemaAttribute::attribute_w_default( "outputi_file" , xs_string , "file name to print output to" , "" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the percentage of positions that are topologically where they should be by the provided span.", attlist );
}

void SpanTopologyMatchPoseFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SpanTopologyMatchPoseFilter::provide_xml_schema( xsd );
}


}
}
