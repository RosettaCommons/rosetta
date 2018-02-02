// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/ScorePoseSegmentFromResidueSelectorFilter.cc
/// @brief  Evaluate RMSD between two poses allowing to select the regions to compare in each pose through ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/filters/ScorePoseSegmentFromResidueSelectorFilter.hh>
#include <protocols/fold_from_loops/filters/ScorePoseSegmentFromResidueSelectorFilterCreator.hh>
#include <protocols/fold_from_loops/movers/SplitAndMixPoseMover.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>


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

static basic::Tracer TR( "protocols.fold_from_loops.ScorePoseSegmentFromResidueSelectorFilter" );

ScorePoseSegmentFromResidueSelectorFilter::ScorePoseSegmentFromResidueSelectorFilter() :
	protocols::filters::Filter( class_name() ),
	score_id_name_( class_name() ),
	residue_select_( default_selector() ),
	scorefxn_( default_scorefxn() ),
	in_context_( false )
{
}

ScorePoseSegmentFromResidueSelectorFilter::~ScorePoseSegmentFromResidueSelectorFilter() {}


core::Real
ScorePoseSegmentFromResidueSelectorFilter::compute( core::pose::Pose const & pose ) const
{

	using namespace core::select::residue_selector;

	core::Real score = 0;

	// Get a secure, detached pose
	core::pose::Pose wpose;
	wpose.detached_copy( pose );

	if ( !in_context_ ) {
		// We need to clean constraints, in case there are constraints between stuff we keep and stuff
		// we delete (it generates an error otherwise).
		wpose.remove_constraints();

		// Remove what we don't want and evaluate.
		grafting::simple_movers::DeleteRegionMover deleter;
		ResidueSelectorCOP query_not( new NotResidueSelector( residue_select_ ) );
		deleter.set_residue_selector( query_not );
		deleter.apply( wpose );
		core::pose::Pose wpose2;
		wpose2.set_new_conformation( wpose.conformation_ptr() );
		wpose2.conformation().detect_disulfides();
		score = (*scorefxn_)( wpose2 );
	} else {
		// We directly score in this conditions.
		(*scorefxn_)( wpose );
		ResidueSubset queryres = residue_select_->apply( wpose );
		// And pick only from the residues of interest.
		for ( core::Size i=1; i<=wpose.size(); ++i ) {
			if ( queryres[i] ) {
				score += wpose.energies().residue_total_energy( i );
			}
		}
	}
	return score;
}

bool
ScorePoseSegmentFromResidueSelectorFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ));
	TR.Trace << score << std::endl;
	return( true );
}

void
ScorePoseSegmentFromResidueSelectorFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ));
	out<<"Segment Score is: " << score<<'\n';
}

core::Real
ScorePoseSegmentFromResidueSelectorFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ));
	return( (core::Real) score );
}

void
ScorePoseSegmentFromResidueSelectorFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	score_id_name( tag->getOption< std::string >( "name" ) );
	in_context( tag->getOption< bool >( "in_context", false ) );
	residue_selector( core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "residue_selector" ), data_map ) );
	if ( tag->hasOption( "scorefxn" ) ) {
		scorefxn( protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", data_map ) );
	}
}

void ScorePoseSegmentFromResidueSelectorFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "in_context", xsct_rosetta_bool, "Do not individualize the selection before scoring.", std::to_string( false ) );

	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "residue_selector", "Selector specifying the segment to be scored." );
	protocols::rosetta_scripts::attributes_for_parse_score_function_w_description( attlist, "scorefxn",
		"Score function to use for evaluation. If not specified calls Rosetta's default full atom score." );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(),
		"Applies the provided score function only to the selected region and adds the score to the pose. By default, it first individualizes the selected region and then "
		"scores it (this is ok, for example, to evaluate specific chains). By using the 'in_context' attribute, it just provides the score of the selected residues"
		"Returns always True.", attlist );
}

std::string ScorePoseSegmentFromResidueSelectorFilterCreator::keyname() const {
	return ScorePoseSegmentFromResidueSelectorFilter::class_name();
}

protocols::filters::FilterOP
ScorePoseSegmentFromResidueSelectorFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ScorePoseSegmentFromResidueSelectorFilter );
}

void ScorePoseSegmentFromResidueSelectorFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ScorePoseSegmentFromResidueSelectorFilter::provide_xml_schema( xsd );
}

} // filters
} // fold_from_loops
} // protocols
