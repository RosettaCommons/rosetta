// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/TMsAACompFilter.cc
/// @brief return the matching (%) between span file determined topology and of the actual pose
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
#include <protocols/simple_filters/TMsAACompFilter.hh>
#include <protocols/simple_filters/TMsAACompFilterCreator.hh>
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

#include <core/scoring/membrane/MPSpanInsertionEnergy.hh>
#include <map>
#include <string>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.TMsAACompFilter" );

protocols::filters::FilterOP
TMsAACompFilterCreator::create_filter() const { return protocols::filters::FilterOP( new TMsAACompFilter ); }

std::string
TMsAACompFilterCreator::keyname() const { return "TMsAAComp"; }

TMsAACompFilter::~TMsAACompFilter(){}

void
TMsAACompFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	threshold_ = tag->getOption< core::Real >( "threshold", 0.5 );
}

bool
TMsAACompFilter::apply( core::pose::Pose const & pose ) const {
	core::Real aa_comp_agree = compute( pose );
	bool const status = aa_comp_agree <= threshold_;
	return status;
}

void
TMsAACompFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real aa_comp_agree = compute( pose );
	out << "TMs amino acid composition RMSD to natural TMs is " << aa_comp_agree << std::endl;
}

core::Real
TMsAACompFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real aa_comp_agree = compute( pose );
	return( aa_comp_agree );
}

core::Real
TMsAACompFilter::compute( core::pose::Pose const & pose ) const {
	std::map < const char, core::Real > wt_aa_comp = {{'L', 0.16}, {'I', 0.1}, {'A', 0.09},
		{'F', 0.08}, {'G', 0.08}, {'V', 0.08}, {'S', 0.07}, {'T', 0.07}, {'M', 0.04}, {'P', 0.04},
		{'Y', 0.03}, {'W', 0.03}, {'N', 0.03}, {'H', 0.02}, {'C', 0.02}, {'Q', 0.01}, {'E', 0.01},
		{'R', 0.01}, {'K', 0.01}, {'D', 0.01}};

	core::scoring::membrane::MPSpanInsertionEnergy mp_span_ins = core::scoring::membrane::MPSpanInsertionEnergy();
	utility::vector1< core::conformation::membrane::Span > spans = mp_span_ins.create_updated_span( pose );
	core::Real dist = 0;
	core::Real num = 0;

	if ( spans.size() == 0 ) {
		TR << "the pose has NO spans, returning -1" << std::endl;
		core::Real minus_one = -1;
		return minus_one;
	}

	for ( core::Size span_i = 1; span_i <= spans.size(); ++span_i   ) {
		//TR << "for span " << spans[span_i].start() << "-" << spans[span_i].end() << std::endl;
		std::map< const char, core::Real > aa_comp = find_aa_composition( pose, spans[span_i].start(),
			spans[span_i].end() );
		std::string all_aas = "ACDEFGHIKLMNPQRSTVWY";
		for ( char & c : all_aas ) {
			TR << "span " << span_i <<  " AA " << c << " pose composition " << aa_comp.at( c ) << " WT composition " << wt_aa_comp.at( c ) << std::endl;
			dist +=  std::pow( aa_comp.at( c ) - wt_aa_comp.at( c ), 2 );
			num++;
		}
	}
	core::Real final_result = std::sqrt( dist / num );
	TR << "final RMSD AA composition " << final_result << std::endl;
	TR.flush();
	return final_result;
}

std::map< const char, core::Real >
TMsAACompFilter::find_aa_composition( core::pose::Pose const & pose,
	core::Size start,
	core::Size end ) const
{
	std::map< const char, core::Real > result;
	std::string all_aas = "ACDEFGHIKLMNPQRSTVWY";
	for ( char & c: all_aas ) result.insert( std::pair< const char, core::Size > ( c, 0 ) );

	core::Real total = 0;
	for ( core::Size i=start; i <= end; ++i ) {
		if ( is_canonical_D_aa( pose.residue( i ).aa() ) || is_canonical_L_aa_or_gly( pose.residue( i ).aa() ) ) {
			TR.Debug << "looking at " << i << " " << pose.residue( i  ).aa() << std::endl;
			result.at( static_cast< char > ( oneletter_code_from_aa( pose.residue( i ).aa() ) ) ) += 1;
			total++;
		}
	}
	std::map< const char, core::Real > normed;
	for ( char & c: all_aas ) {
		normed.insert( std::pair< const char, core::Real> ( c, result.at( c ) / total ) );
	}
	return normed;
}

std::string TMsAACompFilter::name() const {
	return class_name();
}

void TMsAACompFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "leave 0.5 for boolean behaviour" , "0.5" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the RMSD between the AA distribution in all psans in the pose copared to the WT distribution", attlist );
}

void TMsAACompFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TMsAACompFilter::provide_xml_schema( xsd );
}

std::string TMsAACompFilter::class_name() {
	return "TMsAAComp";
}

}
}
