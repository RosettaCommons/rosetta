// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/hbnet/HBNetScore.cc
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <protocols/hbnet/HBNetScore.hh>
#include <protocols/hbnet/HBNetScoreFilterCreator.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <core/scoring/hbonds/HBondSet.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>

#include <protocols/filters/filter_schemas.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/simple_filters/BuriedUnsatHbondFilter.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "protocols.hbnet.HBNetScore" );

namespace protocols {
namespace hbnet {

//Constructor
HBNetScore::HBNetScore() :
	protocols::filters::Filter( "HBNetScore" ),
	threshold_( 0 ),
	hbond_threshold_( -0.25 )
{
	scorename_ = name();
}

HBNetScore::HBNetScore( protocols::filters::Filter const & src ) :
	HBNetScore( static_cast< HBNetScore const & > (src) ){}

HBNetScore::HBNetScore( HBNetScore const & src ) :
	protocols::filters::Filter( "HBNetScore" ),
	threshold_( src.threshold_ ),
	hbond_threshold_( src.hbond_threshold_ )
{
	scorename_ = name();
}


//Destructor
HBNetScore::~HBNetScore() = default;

void HBNetScore::report( std::ostream & out, core::pose::Pose const & pose) const
{
	out << "HBNetScore: " << report_sm( pose ) << std::endl;
}

core::Real HBNetScore::get_score( core::pose::Pose & pose ) const
{
	core::select::residue_selector::ResiduePDBInfoHasLabelSelector selector( "HBNet" );
	utility::vector1< bool > residue_is_in_network;
	residue_is_in_network = selector.apply( pose );

	core::Size resid_count = 0;
	for ( core::Size resid = 1; resid <= pose.size(); ++resid ) {
		if ( residue_is_in_network[ resid ] ) {
			++resid_count;
		}
	}
	if ( resid_count < 2 ) return 9999;

	core::scoring::hbonds::HBondSetOP hbond_set( new core::scoring::hbonds::HBondSet( pose, false ) );

	core::Real reward = 0;
	for ( core::Size hbid = 1; hbid <= hbond_set->nhbonds(); ++hbid ) {
		core::scoring::hbonds::HBondCOP hbond = hbond_set->hbond_cop( hbid );
		if ( residue_is_in_network[ hbond->don_res() ] && residue_is_in_network[ hbond->acc_res() ] && hbond->energy() <= hbond_threshold_ ) {
			if ( ! hbond->don_hatm_is_backbone() || ! hbond->acc_atm_is_backbone() ) {
				reward += hbond->energy();
			}
		}
	}

	reward /= resid_count;

	protocols::simple_filters::BuriedUnsatHbondFilter buns_filter( 0 );// now updated to new unsat filter
	return reward + buns_filter.compute( pose );
}

bool HBNetScore::apply( core::pose::Pose const & pose ) const
{
	return report_sm(pose) <= threshold_;
}

void
HBNetScore::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Score must be less than or equal to this value to pass", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "hbond_threshold", xsct_real, "Disregard hbonds weaker than this", "-0.25" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, "HBNetScore", "Reads stored residue subset created by HBNet(StapleInterface) and sums all of the sc-sc and sc-bb hbonds in the network. Adds BuriedUnsatHbondFilter score for a penalty.", attlist );

}

void HBNetScore::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	threshold_ = tag->getOption<core::Real>( "threshold", 0.0 );
	hbond_threshold_ = tag->getOption<core::Real>( "hbond_threshold", 0.0 );
}

std::string HBNetScoreFilterCreator::keyname() const {
	return "HBNetScore";
}

protocols::filters::FilterOP
HBNetScoreFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new HBNetScore );
}

void HBNetScoreFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HBNetScore::provide_xml_schema( xsd );
}


} //hbnet
} //protocols
