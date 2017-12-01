// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ShowConstraintsFilter.cc
/// @brief iterate over all constraints on the pose and invoke their show methods
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
#include <protocols/simple_filters/ShowConstraintsFilter.hh>
#include <protocols/simple_filters/ShowConstraintsFilterCreator.hh>
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

#include <core/scoring/constraints/ConstraintSet.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.ShowConstraintsFilter" );

protocols::filters::FilterOP
ShowConstraintsFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ShowConstraintsFilter ); }

std::string
ShowConstraintsFilterCreator::keyname() const { return "ShowConstraints"; }

ShowConstraintsFilter::~ShowConstraintsFilter(){}

void
ShowConstraintsFilter::parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
}

bool
ShowConstraintsFilter::apply( core::pose::Pose const & pose ) const {
	core::Real result = compute( pose );
	return result;
}

void
ShowConstraintsFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real result = compute( pose );
	out << "the results is " << result << std::endl;
}

core::Real
ShowConstraintsFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real result = compute( pose );
	return( result );
}

core::Real
ShowConstraintsFilter::compute( core::pose::Pose const & pose ) const {
	core::scoring::constraints::ConstraintCOPs csts = pose.constraint_set()->get_all_constraints();
	TR << "Showing all constraints" << std::endl;
	for ( auto cst : csts ) {
		TR << "==================================" << std::endl;
		cst->show( TR );
		TR << "==================================" << std::endl;
	}
	return 1;
}

std::string ShowConstraintsFilter::name() const {
	return class_name();
}

void ShowConstraintsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//attlist + XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "leave 0.5 for boolean behaviour" , "0.5" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "iterate over all constraints in the pose, and invoke their show method", attlist );
}

void ShowConstraintsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ShowConstraintsFilter::provide_xml_schema( xsd );
}

std::string ShowConstraintsFilter::class_name() {
	return "ShowConstraints";
}

}
}
