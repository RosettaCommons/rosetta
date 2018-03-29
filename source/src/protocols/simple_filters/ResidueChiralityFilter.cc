// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueChiralityFilter.cc
/// @brief checks the chirality of a specific residues, whether it is D or L
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
#include <core/chemical/ResidueType.hh>
#include <map>
#include <numeric/random/random.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_filters/ResidueChiralityFilter.hh>
#include <protocols/simple_filters/ResidueChiralityFilterCreator.hh>
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

static basic::Tracer TR( "protocols.simple_filters.ResidueChiralityFilter" );

protocols::filters::FilterOP
ResidueChiralityFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ResidueChiralityFilter ); }

std::string
ResidueChiralityFilterCreator::keyname() const { return "ResidueChirality"; }

ResidueChiralityFilter::~ResidueChiralityFilter(){}

void
ResidueChiralityFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	required_type_ = tag->getOption< std::string >( "residue_type", "L" );
	res_num_ = tag->getOption< core::Size >( "res_num", 1 );
}

bool
ResidueChiralityFilter::apply( core::pose::Pose const & pose ) const {
	std::string type_found = compute( pose );
	bool const status = type_found == required_type_;
	return status;
}

void
ResidueChiralityFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	std::string type_found = compute( pose );
	out << "found residue " << res_num_ << " type to be " << type_found << std::endl;
}

core::Real
ResidueChiralityFilter::report_sm( core::pose::Pose const & pose ) const {
	std::string type_found = compute( pose );
	return( type_found == required_type_ );
}

std::string
ResidueChiralityFilter::compute( core::pose::Pose const & pose ) const {
	bool is_d = pose.residue(res_num_).type().is_d_aa();
	std::string residue_type = "L";
	if ( is_d ) residue_type = "D";
	return( residue_type );
}

std::string ResidueChiralityFilter::name() const {
	return class_name();
}

std::string ResidueChiralityFilter::class_name() {
	return "ResidueChirality";
}

void ResidueChiralityFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "residue_type" , xs_string , "the type you want to pass" , "L" )
		+ XMLSchemaAttribute::attribute_w_default( "res_num" , xs_integer , "which residue to test" , "1" ) ;
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "fails if the res_num residue is in the wrong chirality (D / L)", attlist );
}

void ResidueChiralityFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueChiralityFilter::provide_xml_schema( xsd );
}

}
}
