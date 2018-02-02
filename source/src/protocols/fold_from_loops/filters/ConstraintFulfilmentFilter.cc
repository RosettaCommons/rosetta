// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/ConstraintFulfilmentFilter.cc
/// @brief  Evaluate RMSD between two poses allowing to select the regions to compare in each pose through ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/filters/ConstraintFulfilmentFilter.hh>
#include <protocols/fold_from_loops/filters/ConstraintFulfilmentFilterCreator.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>


#include <algorithm>
#include <list>
#include <map>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace fold_from_loops {
namespace filters {

static basic::Tracer TR( "protocols.fold_from_loops.ConstraintFulfilmentFilter" );

ConstraintFulfilmentFilter::ConstraintFulfilmentFilter() :
	protocols::filters::Filter( class_name() ),
	distance_( default_distance() ),
	angle_( default_angle() ),
	dihedral_( default_dihedral() ),
	fulfil_threshold_( 0.33 )
{
}

ConstraintFulfilmentFilter::~ConstraintFulfilmentFilter() {}


core::Real
ConstraintFulfilmentFilter::compute( core::pose::Pose const & pose ) const
{
	core::scoring::constraints::ConstraintSetCOP constraints = pose.constraint_set();
	runtime_assert_msg( !constraints->is_empty(), "The pose needs to have constraints to be evaluated.");
	core::Real distance_total_score( 0.0 ) ;
	core::Size distance_count( 0 );
	core::Real angle_total_score( 0.0 ) ;
	core::Size angle_count( 0 );
	core::Real dihedral_total_score( 0.0 ) ;
	core::Size dihedral_count( 0 );
	for ( auto cst : constraints->get_all_constraints() ) {
		std::string csttype=cst->type();
		if ( distance_ ) {
			if ( csttype == "AtomPair" ) {
				distance_total_score += cst->show_violations( TR, pose, 500 );
				++distance_count;
			}
		}
		if ( angle_ ) {
			if ( csttype == "Angle" ) {
				angle_total_score += cst->show_violations( TR, pose, 500 );
				++angle_count;
			}
		}
		if ( dihedral_ ) {
			if ( csttype == "Dihedral" ) {
				dihedral_total_score += cst->show_violations( TR, pose, 500 );
				++dihedral_count;
			}
		}
	}
	if ( distance_ and distance_count > 0 ) {
		distance_total_score = distance_total_score / distance_count;
		TR << distance_count << " distance constraint(s) " << distance_total_score << "% are violated" << std::endl;
	}
	if ( angle_ and angle_count > 0 ) {
		angle_total_score = angle_total_score / angle_count;
		TR << angle_count << " angle constraint(s) " << angle_total_score << "% are violated" << std::endl;
	}
	if ( dihedral_ and dihedral_count > 0 ) {
		dihedral_total_score = dihedral_total_score / dihedral_count;
		TR << dihedral_count << " dihedral constraint(s) " << dihedral_total_score << "% are violated" << std::endl;
	}
	return distance_total_score + angle_total_score + dihedral_total_score;
}

bool
ConstraintFulfilmentFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ));
	return( score < fulfil_threshold_ );
}

void
ConstraintFulfilmentFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ));
	out<<"Constraint Score is: " << score<<'\n';
}

core::Real
ConstraintFulfilmentFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ));
	return( (core::Real) score );
}

void
ConstraintFulfilmentFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & , protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &  )
{
	distance( tag->getOption< bool >( "distance", default_distance() ) );
	angle( tag->getOption< bool >( "angle", default_angle() ) );
	dihedral( tag->getOption< bool >( "dihedral", default_dihedral() ) );
}

void ConstraintFulfilmentFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "distance", xsct_rosetta_bool, "Select distance (csttype=AtomPair) constraints.", std::to_string( default_distance() ) )
		+ XMLSchemaAttribute::attribute_w_default( "angle", xsct_rosetta_bool, "Select angle (csttype=Angle) constraints.", std::to_string( default_angle() ) )
		+ XMLSchemaAttribute::attribute_w_default( "dihedral", xsct_rosetta_bool, "Select dihedral (csttype=Dihedral) constraints.", std::to_string( default_dihedral() ) );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Counts how many constraints the pose violates (AtomPair/Angle/Dihedral). "
		"Returns False if more than .33%", attlist );
}

std::string ConstraintFulfilmentFilterCreator::keyname() const {
	return ConstraintFulfilmentFilter::class_name();
}

protocols::filters::FilterOP
ConstraintFulfilmentFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ConstraintFulfilmentFilter );
}

void ConstraintFulfilmentFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConstraintFulfilmentFilter::provide_xml_schema( xsd );
}

} // filters
} // fold_from_loops
} // protocols
