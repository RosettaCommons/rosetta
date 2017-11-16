// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/HolesFilter.cc
/// @brief filter structures by will's hole value
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/simple_filters/HolesFilter.hh>
#include <protocols/simple_filters/HolesFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/database/open.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer tr( "protocols.simple_filters.HolesFilter" );

// @brief default constructor
HolesFilter::HolesFilter():
	Filter( "Holes" ),
	filtered_value_( 2.0 ),
	cmd_( "" ),
	exclude_bb_atoms_(false),
	normalize_per_atom_(false),
	normalize_per_residue_(false),
	residue_selector_( /* NULL */ )
{}

// @brief copy constructor
HolesFilter::HolesFilter( HolesFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	filtered_value_( rval.filtered_value_ ),
	cmd_( rval.cmd_ ),
	exclude_bb_atoms_( rval.exclude_bb_atoms_ ),
	normalize_per_atom_( rval.normalize_per_atom_ ),
	normalize_per_residue_( rval.normalize_per_residue_ ),
	residue_selector_( rval.residue_selector_ )
{}

// @brief set filtered value
void HolesFilter::filtered_value( Real const & value )
{
	filtered_value_ = value;
}

/// @brief
HolesFilter::Real
HolesFilter::report_sm( Pose const & pose ) const
{
	return  compute( pose );
}

/// @brief
void
HolesFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "Holes: " <<  compute( pose ) << std::endl;
}

/// @brief
HolesFilter::Real
HolesFilter::compute( Pose const & pose ) const
{

	using core::scoring::packing::HolesResult;
	using core::scoring::packing::HolesParams;
	using core::scoring::packing::compute_holes_score;

	Size MAX_RES = 5000;
	runtime_assert( pose.size() <= MAX_RES );

	if ( residue_selector_ ) { ///// sboyken 17/08/26 adding residue selector option

		// this option still uses the full pose to compute the presence/absence of holes
		//   but only sums over the atoms in the residues of the residue selector

		tr << " computing using residue selector " << std::endl;

		// according to Will Sheffler, best to use dec15 params for this
		HolesParams dec15_params;
		dec15_params.read_data_file(basic::database::full_name("scoring/rosettaholes/decoy15.params"));

		HolesResult result = compute_holes_score( pose, dec15_params );

		core::id::AtomID_Map< core::Real > per_atom_holes_scores = result.atom_scores;

		core::select::residue_selector::ResidueSubset rs( residue_selector_->apply( pose ) );

		Size total_res(0);
		Size total_atoms(0);
		core::Real total_holes_score(0.0);
		for ( Size r = 1; r <= rs.size(); ++r ) {
			if ( rs[ r ] ) {
				core::Real residue_holes_score(0.0);

				//for ( Size at_r = 1; at_r <= pose.residue( r ).natoms(); ++at_r ) {
				for ( Size at_r = 1; at_r <= per_atom_holes_scores.n_atom( r ); ++at_r ) {
					// skip backbone atoms if desired
					if ( exclude_bb_atoms_ && pose.residue(r).atom_is_backbone( at_r ) ) continue;
					core::id::AtomID const atom_id( at_r, r );
					core::Real atom_holes_score = per_atom_holes_scores[ atom_id ];
					residue_holes_score += atom_holes_score;
					total_atoms++;
				}
				total_holes_score += residue_holes_score;
				total_res++;
			}
		}
		tr << " total holes score for selection = " << total_holes_score << std::endl;
		tr << " normalize_per_atom holes score = " << total_holes_score / (core::Real)total_atoms << std::endl;
		tr << " normalize_per_residue holes score = " << total_holes_score / (core::Real)total_res << std::endl;
		if ( normalize_per_atom_ ) return total_holes_score / (core::Real)total_atoms;
		if ( normalize_per_residue_ ) return total_holes_score / (core::Real)total_res;
		return total_holes_score;

	} ///// sboyken 17/08/26 adding residue selector option

	tr << " computing using full pose " << std::endl;

	HolesResult result = compute_holes_score( pose, cmd_ );
	return result.dec15_score;
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool HolesFilter::apply( Pose const & pose ) const
{
	Real value = compute( pose );
	if ( value < filtered_value_ ) {
		tr << "Successfully filtered: " << value << std::endl;
		return true;
	} else {
		tr << "Filter failed current/threshold=" << value << "/" << filtered_value_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
HolesFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	// set filtered type
	cmd_ = tag->getOption<String>( "cmd", "" );
	if ( cmd_ == "" ) {
		tr << "cmd in xml file is emptry, so using -holes::dalphaball is expected. " << std::endl;
		tr << "if both are empty, this gonna be crash. " << std::endl;
	}

	// set threshold
	filtered_value_ = tag->getOption<Real>( "threshold", 2.0 );
	tr << "Structures which have holes value less than " << filtered_value_ << " will be filtered." << std::endl;

	exclude_bb_atoms_ = tag->getOption<bool>( "exclude_bb_atoms", false );
	normalize_per_atom_ = tag->getOption<bool>( "normalize_per_atom", false );
	normalize_per_residue_ = tag->getOption<bool>( "normalize_per_residue", false );

	// sboyken 17/08/26 adding residue selector option
	if ( tag->hasOption("residue_selector") ) {

		tr << " Holes will use residue selector " << std::endl;

		// if not user-specificied default to normalize_per_residue
		if ( ! tag->hasOption("normalize_per_atom") && ! tag->hasOption("normalize_per_residue") ) {
			normalize_per_atom_ = true;
		}

		std::string const selector_name ( tag->getOption< std::string >( "residue_selector" ) );
		if ( tr.visible() ) tr << " Set selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP selector;
		try {
			selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddCompositionConstraintMover::parse_tag()\n" + e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_message );
		}
		runtime_assert( selector );
		residue_selector_ = selector->clone();

		tr << " residue selector applied succesfully" << std::endl;

		//THIS NEEDS TO HAPPEN AT APPLY TIME, NOT PARSE TIME!!!!!
		//core_residues_ = selector->apply( pose );
	}
}

// XRW TEMP filters::FilterOP
// XRW TEMP HolesFilterCreator::create_filter() const { return filters::FilterOP( new HolesFilter ); }

// XRW TEMP std::string
// XRW TEMP HolesFilterCreator::keyname() const { return "Holes"; }

std::string HolesFilter::name() const {
	return class_name();
}

std::string HolesFilter::class_name() {
	return "Holes";
}

void HolesFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("cmd", xs_string, "expects you to use -holes::dalphaball", "XRW TO DO")
		+ XMLSchemaAttribute::attribute_w_default("threshold", xsct_real, "Threshold value for the filter", "2.0")
		+ XMLSchemaAttribute::attribute_w_default("exclude_bb_atoms", xsct_rosetta_bool, "don't include backbone (bb) atoms in residue selection (residue_selector case)", "false")
		+ XMLSchemaAttribute::attribute_w_default("normalize_per_residue", xsct_rosetta_bool, "for residue selector case, nomralize per residue", "false")
		+ XMLSchemaAttribute::attribute_w_default("normalize_per_atom", xsct_rosetta_bool,
		"for residue selector case, normalize per atom; defaults to true if residue_selector if normalize_per_atom and normalize_per_residue not explicitly defined by user", "false")
		+ XMLSchemaAttribute("residue_selector", xs_string,
		"only calculate holes score for residues in residue_selector; The holes calculation is performed on the Pose as whole (ignoring the ResidueSelector), but the score that's reported is the sum of only the atoms in the residue selector.  The Holes score is a sum of individual atoms/residues anyway (technically PoseBalls), so by only reporting a specific selection, get location-specific data.");


	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Calculates holes value(?) and filters based upon topology", attlist );

}

std::string HolesFilterCreator::keyname() const {
	return HolesFilter::class_name();
}

protocols::filters::FilterOP
HolesFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new HolesFilter );
}

void HolesFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HolesFilter::provide_xml_schema( xsd );
}



} // filters
} // protocols
