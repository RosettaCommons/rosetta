// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/InterfaceHydrophobicResidueContactsFilter.cc
/// @brief Counts the number of hydrophobic residues on the target that have at least a certain hydrophobic score
/// @author Brian Coventry (bcov@uw.edu)

// Unit headers
#include <protocols/simple_filters/InterfaceHydrophobicResidueContactsFilter.hh>
#include <protocols/simple_filters/InterfaceHydrophobicResidueContactsFilterCreator.hh>


// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

#include <protocols/filters/filter_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>


// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.hh>


namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.InterfaceHydrophobicResidueContactsFilter" );


const core::Size InterfaceHydrophobicResidueContactsFilter::DEFAULT_THRESHOLD = 5;
const std::string InterfaceHydrophobicResidueContactsFilter::DEFAULT_APOLAR_RES = "ALA,CYS,CYD,PHE,ILE,LEU,MET,PRO,THR,VAL,TRP,TYR";
const core::Real InterfaceHydrophobicResidueContactsFilter::DEFAULT_SCORE_CUT = -0.5;
const std::set<core::scoring::ScoreType> InterfaceHydrophobicResidueContactsFilter::SCORE_TERMS = std::set<core::scoring::ScoreType> {
core::scoring::ScoreType::fa_rep,
core::scoring::ScoreType::fa_atr,
core::scoring::ScoreType::fa_sol
};


InterfaceHydrophobicResidueContactsFilter::~InterfaceHydrophobicResidueContactsFilter() = default;


core::Size InterfaceHydrophobicResidueContactsFilter::compute( core::pose::Pose const & pose ) const {
	sanity_check();
	core::pose::Pose complex = pose;
	core::pose::Pose apo = pose;

	(*scorefxn_)( complex );

	core::select::residue_selector::ResidueSubset target_subset = target_selector_->apply( apo );
	core::select::residue_selector::ResidueSubset binder_subset = binder_selector_->apply( apo );

	numeric::xyzVector<core::Real> x_unit { 1, 0, 0 };

	protocols::toolbox::pose_manipulation::rigid_body_move( x_unit, 0, x_unit, 10000, apo, binder_subset );

	(*scorefxn_)( apo );

	core::Size hydrophobic_residue_contacts = 0;

	TR << "Interface hydrophobic interactions:" << std::endl;
	for ( core::Size seqpos = 1; seqpos <= complex.size(); seqpos ++ ) {
		if ( ! target_subset[ seqpos ] ) continue;
		if ( apolar_res_set_.count( complex.residue( seqpos ).name3() ) == 0 ) continue;

		core::Real complex_score = complex.energies().residue_total_energy( seqpos );
		core::Real apo_score = apo.energies().residue_total_energy( seqpos );

		core::Real delta = complex_score - apo_score;
		bool pass = delta < score_cut_;

		if ( delta < 0 ) {
			TR << ( pass ? "* " : "  " ) << seqpos << " " << pose.residue(seqpos).name3() << " " << delta << std::endl;
		}

		if ( pass ) {
			hydrophobic_residue_contacts++;
		}
	}

	return hydrophobic_residue_contacts;
}

void
InterfaceHydrophobicResidueContactsFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
) {

	if ( tag->hasOption( "threshold" ) ) {
		hydrophobic_residue_contacts_threshold_ = tag->getOption<core::Size>( "threshold" );
	}
	if ( tag->hasOption( "target_selector" ) ) {
		target_selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data, "target_selector" );
	}
	if ( tag->hasOption( "binder_selector" ) ) {
		binder_selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data, "binder_selector" );
	}
	if ( tag->hasOption( "score_cut" ) ) {
		score_cut_ = tag->getOption<core::Real>( "score_cut" );
	}
	if ( tag->hasOption( "apolar_res" ) ) {
		set_apolar_res( tag->getOption<std::string>( "apolar_res" ) );
	}
	if ( tag->hasOption( "scorefxn" ) ) {
		set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", data ) );
	}


	TR << "Defined with threshold " << hydrophobic_residue_contacts_threshold_ << " score_cut " << score_cut_ << " apolor res ";
	bool first = true;
	for ( std::string const & res : apolar_res_set_ ) {
		if ( ! first ) TR << ",";
		first = false;
		TR << res;
	}
	TR << std::endl;

	sanity_check();
}

void
InterfaceHydrophobicResidueContactsFilter::set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ) {
	if ( ! scorefxn ) return;
	scorefxn_ = scorefxn->clone();

	utility::vector1<core::scoring::ScoreType> nonzero_scoretypes = scorefxn_->get_nonzero_weighted_scoretypes();

	for ( core::scoring::ScoreType scoretype : nonzero_scoretypes ) {
		if ( SCORE_TERMS.count( scoretype ) != 0 ) {
			if ( scorefxn_->weights()[scoretype] == 0 ) {
				TR.Warning << "Warning: " << core::scoring::name_from_score_type(scoretype) << " weight is 0!!!" << std::endl;
			}
			continue;
		}

		scorefxn_->set_weight( scoretype, 0 );
	}
}

void
InterfaceHydrophobicResidueContactsFilter::set_apolar_res( std::string const & apolar_res ) {
	apolar_res_set_.clear();
	utility::vector1<std::string> split = utility::string_split_simple( apolar_res, ',' );

	for ( std::string const & res : split ) {
		apolar_res_set_.insert( res );
	}
}

void
InterfaceHydrophobicResidueContactsFilter::sanity_check() const {
	if ( ! scorefxn_ ) {
		utility_exit_with_message( "protocols.simple_filters.InterfaceHydrophobicResidueContactsFilter: No scorefxn set!!!" );
	}
	if ( ! target_selector_ ) {
		utility_exit_with_message( "protocols.simple_filters.InterfaceHydrophobicResidueContactsFilter: No target_selector set!!!" );
	}
	if ( ! binder_selector_ ) {
		utility_exit_with_message( "protocols.simple_filters.InterfaceHydrophobicResidueContactsFilter: No binder_selector set!!!" );
	}
	if ( apolar_res_set_.size() == 0 ) {
		utility_exit_with_message( "protocols.simple_filters.InterfaceHydrophobicResidueContactsFilter: No apolar_res set!!!" );
	}
	if ( apolar_res_set_.size() == 0 ) {
		bool any_oneletter = false;
		for ( std::string const & res : apolar_res_set_ ) {
			if ( res.size() == 1 ) any_oneletter = true;
		}
		if ( any_oneletter ) TR.Warning << "Warning: apolar_res should be name3 not name1!!!" << std::endl;
	}
	if ( score_cut_ > 0 ) {
		TR.Warning << "Warning: Most people set score_cut to be less than 0 !!!" << std::endl;
	}
}

bool
InterfaceHydrophobicResidueContactsFilter::apply( core::pose::Pose const & pose ) const {
	core::Size hydrophobic_residue_contacts = compute( pose );
	return hydrophobic_residue_contacts >= hydrophobic_residue_contacts_threshold_;
}

void
InterfaceHydrophobicResidueContactsFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size const hydrophobic_residue_contacts( compute( pose ));
	out<<"hydrophobic_residue_contacts "<< hydrophobic_residue_contacts << '\n';
}

core::Real
InterfaceHydrophobicResidueContactsFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Size const hydrophobic_residue_contacts( compute( pose ));
	return( (core::Real) hydrophobic_residue_contacts );
}

std::string InterfaceHydrophobicResidueContactsFilter::name() const {
	return class_name();
}

std::string InterfaceHydrophobicResidueContactsFilter::class_name() {
	return "InterfaceHydrophobicResidueContacts";
}

void InterfaceHydrophobicResidueContactsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {

	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute(
		"name", xs_string,
		"Name of the mover");

	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "target_selector",
		"Which residues should be counted as hydrophobic target residues?" );

	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "binder_selector",
		"Which residues are part of the binder so that the interface can be defined?" );

	rosetta_scripts::attributes_for_parse_score_function_w_description_when_required(attlist, "scorefxn",
		"Which scorefunction do you want to use? Only fa_rep, fa_sol, and fa_atr will be kept.");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"threshold", xs_integer,
		"How many hydrophobic residues should have at least the score cutoff?",
		utility::to_string(DEFAULT_THRESHOLD));

	attlist + XMLSchemaAttribute::attribute_w_default(
		"score_cut", xsct_real,
		"What score cut should be used to determine whether a hydrophobic residue is contacted?",
		utility::to_string(DEFAULT_SCORE_CUT));

	attlist + XMLSchemaAttribute::attribute_w_default(
		"apolar_res", xs_string,
		"What residues should count as hydrophobic residues (comma separated name3 list)?",
		DEFAULT_APOLAR_RES);

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(),
		"Counts the number of hydrophobic residues on the target that have at least a certain hydrophobic score", attlist );
}

std::string InterfaceHydrophobicResidueContactsFilterCreator::keyname() const {
	return InterfaceHydrophobicResidueContactsFilter::class_name();
}

protocols::filters::FilterOP
InterfaceHydrophobicResidueContactsFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new InterfaceHydrophobicResidueContactsFilter );
}

void InterfaceHydrophobicResidueContactsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterfaceHydrophobicResidueContactsFilter::provide_xml_schema( xsd );
}


}
}
