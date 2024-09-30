// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/util.cc
/// @brief  Utility functions for the residue selector classes, primarily, in constructing XML-Schema type definitions
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/select/residue_selector/util.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/OrResidueSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/SelectorLogicParser.hh>
#include <core/select/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>

// C++ headers
#include <sstream>
#include <algorithm>
#include <iterator>
#include <regex>

static basic::Tracer TR( "core.select.residue_selector.util" );

namespace core {
namespace select {
namespace residue_selector {

// Attribute::Attribute() : required( false ) {}
//
// Attribute::Attribute(
//  std::string const & name_in,
//  std::string const & type_in,
//  bool required_in
// ) :
//  name( name_in ),
//  type( type_in ),
//  required( required_in )
// {}


std::map< core::Size, core::Size >
get_residue_mapping_from_selectors(
	select::residue_selector::ResidueSelectorCOP residue_selector,
	select::residue_selector::ResidueSelectorCOP residue_selector_ref,
	core::pose::Pose const & pose,
	core::pose::Pose const & ref_pose,
	bool desymmetrize_selectors
){

	using namespace utility;

	vector1< core::Size > selector_res;
	vector1< core::Size > reference_res;

	if ( residue_selector ) {
		vector1< bool > mask = residue_selector->apply( pose );
		if ( core::pose::symmetry::is_symmetric( pose ) && desymmetrize_selectors )  {
			mask = core::select::get_master_subunit_selection(pose, mask);
		}
		selector_res = get_residues_from_subset( mask );
	}
	if ( residue_selector_ref ) {
		vector1< bool > mask = residue_selector_ref->apply( ref_pose );
		if ( core::pose::symmetry::is_symmetric( ref_pose ) && desymmetrize_selectors )  {
			mask = core::select::get_master_subunit_selection(pose, mask);
		}
		reference_res = get_residues_from_subset( mask );
	}

	//Fail if residue selector selections do not match.
	if ( residue_selector && residue_selector_ref ) {
		if ( selector_res.size() != reference_res.size() ) {
			utility_exit_with_status("Both set residue selectors must select the same number of residues in order to run RMSD calculation!");
		}
	}

	//Fail if only reference selector was set.
	if ( residue_selector_ref && (! residue_selector) ) {
		utility_exit_with_message("Cannot only set the reference residue selector! If they are both the same, please set the main selector!");
	}

	//Setup the main residue mapping we will use depending on what is set.
	std::map< core::Size, core::Size > residue_map;
	if ( residue_selector && residue_selector_ref ) {
		for ( core::Size i = 1; i <= selector_res.size(); ++i ) {
			residue_map[ selector_res[ i] ] = reference_res[ i ];
		}
	} else if ( residue_selector ) {
		for ( core::Size res : selector_res ) {
			residue_map[ res ] = res;
		}
	} else {
		for ( core::Size res = 1; res <= pose.size(); ++res ) {
			residue_map[ res ] = res;
		}
	}
	return residue_map;
}


std::vector < std::map< core::Size, core::Size > >
get_cyclic_pose_residue_mappings_from_selectors(
	select::residue_selector::ResidueSelectorCOP residue_selector,
	select::residue_selector::ResidueSelectorCOP residue_selector_ref,
	core::pose::Pose const & pose,
	core::pose::Pose const & ref_pose,
	bool desymmetrize_selectors
){
	using namespace utility;

	vector1< core::Size > selector_res;
	vector1< core::Size > reference_res;

	if ( residue_selector ) {
		vector1< bool > mask = residue_selector->apply( pose );
		if ( core::pose::symmetry::is_symmetric( pose ) && desymmetrize_selectors )  {
			mask = core::select::get_master_subunit_selection(pose, mask);
		}
		selector_res = get_residues_from_subset( mask );
	}
	if ( residue_selector_ref ) {
		vector1< bool > mask = residue_selector_ref->apply( ref_pose );
		if ( core::pose::symmetry::is_symmetric( ref_pose ) && desymmetrize_selectors )  {
			mask = core::select::get_master_subunit_selection(pose, mask);
		}
		reference_res = get_residues_from_subset( mask );
	}

	//Fail if residue selector selections do not match.
	if ( residue_selector && residue_selector_ref ) {
		if ( selector_res.size() != reference_res.size() ) {
			utility_exit_with_status("Both set residue selectors must select the same number of residues in order to run RMSD calculation!");
		}
	}

	//Fail if only reference selector was set.
	if ( residue_selector_ref && (! residue_selector) ) {
		utility_exit_with_message("Cannot only set the reference residue selector! If they are both the same, please set the main selector!");
	}

	//Setup the main residue mapping we will use depending on what is set.
	std::vector < std::map< core::Size, core::Size > > residue_maps;
	std::map< core::Size, core::Size > residue_map;
	bool matched(true);

	if ( residue_selector && residue_selector_ref ) {
		for ( core::Size i = 1; i <= selector_res.size(); ++i ) {
			residue_map[ selector_res[ i] ] = reference_res[ i ];
		}
		residue_maps.push_back(residue_map);
		for ( core::Size shift_i = 1; shift_i <= selector_res.size()-1 ; ++shift_i ) {
			matched = true;
			residue_map.clear();
			for ( core::Size i = 1; i <= selector_res.size(); ++i ) {
				core::Size ndx = (i+shift_i) % selector_res.size();
				ndx = (ndx==0) ? selector_res.size() : ndx;
				if ( TR.Debug.visible() ) {
					TR.Debug << "pose residue name: " << pose.residue(selector_res[ndx]).name()
						<< ", " << "ref pose residue name: "
						<< ref_pose.residue(reference_res[i]).name()
						<< std::endl;
				}
				if ( pose.residue(selector_res[ndx]).name() != ref_pose.residue(reference_res[i]).name() ) {
					matched = false;
					break;
				}
				residue_map[ selector_res[ ndx ] ] = reference_res[ i ];
			}
			if ( matched ) residue_maps.push_back(residue_map);
		}
	} else if ( residue_selector ) {
		for ( core::Size res : selector_res ) {
			residue_map[ res ] = res;
		}
		residue_maps.push_back(residue_map);
		for ( core::Size shift_i = 1; shift_i <= selector_res.size()-1 ; ++shift_i ) {
			matched = true;
			residue_map.clear();
			for ( core::Size i = 1; i <= selector_res.size(); ++i ) {
				core::Size ndx = (i+shift_i) % selector_res.size();
				ndx = (ndx==0) ? selector_res.size() : ndx;
				if ( TR.Debug.visible() ) {
					TR.Debug << "pose residue name: " << pose.residue(selector_res[ndx]).name() << ", " << "ref pose residue name: "
						<< ref_pose.residue(selector_res[i]).name()
						<< std::endl;
				}
				if ( pose.residue(selector_res[ndx]).name() != ref_pose.residue(selector_res[i]).name() ) {
					matched = false;
					break;
				}
				residue_map[ selector_res[ ndx ] ] = selector_res[ i ];
			}
			if ( matched ) residue_maps.push_back(residue_map);
		}
	} else {
		for ( core::Size res = 1; res <= pose.size(); ++res ) {
			residue_map[ res ] = res;
		}
		residue_maps.push_back(residue_map);
		for ( core::Size shift_i = 1; shift_i <= pose.size()-1 ; ++shift_i ) {
			matched = true;
			residue_map.clear();
			for ( core::Size i = 1; i <= pose.size(); ++i ) {
				core::Size ndx = (i+shift_i) % pose.size();
				ndx = (ndx==0) ? pose.size() : ndx;
				if ( TR.Debug.visible() ) {
					TR.Debug << "pose residue name: " << pose.residue(ndx).name()  << ", " << "ref pose residue name: "
						<< ref_pose.residue(i).name()
						<< std::endl;
				}
				if ( pose.residue(ndx).name() != ref_pose.residue(i).name() ) {
					matched = false;
					break;
				}
				residue_map[ ndx ] = i;
			}
			if ( matched ) residue_maps.push_back(residue_map);
		}
	}
	return residue_maps;
}


std::string
complex_type_name_for_residue_selector( std::string const & rs_type )
{
	return "rs_" + rs_type + "_type";
}

void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;

	std::string citation_string;
	if ( xsd.include_citation_info() ) {
		citation_string = "\n\n" + ResidueSelectorFactory::get_instance()->get_citation_humanreadable( rs_type );
	}

	ct_gen.complex_type_naming_func( & complex_type_name_for_residue_selector )
		.element_name( rs_type )
		.description( description + citation_string )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );
}

void
xsd_type_definition_w_attributes_and_optional_subselector(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelement;
	subelement.add_group_subelement( & ResidueSelectorFactory::residue_selector_xml_schema_group_name );

	std::string citation_string;
	if ( xsd.include_citation_info() ) {
		citation_string = "\n\n" + ResidueSelectorFactory::get_instance()->get_citation_humanreadable( rs_type );
	}

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_residue_selector )
		.element_name( rs_type )
		.description( description + citation_string )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.set_subelements_single_appearance_optional( subelement )
		.write_complex_type_to_schema( xsd );
}

void
xsd_type_definition_w_attributes_and_optional_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelement;
	subelement.add_group_subelement( & ResidueSelectorFactory::residue_selector_xml_schema_group_name );

	std::string citation_string;
	if ( xsd.include_citation_info() ) {
		citation_string = "\n\n" + ResidueSelectorFactory::get_instance()->get_citation_humanreadable( rs_type );
	}

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_residue_selector )
		.element_name( rs_type )
		.description( description + citation_string )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subelement )
		.write_complex_type_to_schema( xsd );
}

void
xsd_type_definition_w_attributes_and_optional_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & description,
	core::Size min_occurrence,
	core::Size max_occurrence,
	utility::tag::AttributeList const & attributes
)
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelement;
	subelement.add_group_subelement( & ResidueSelectorFactory::residue_selector_xml_schema_group_name );

	std::string citation_string;
	if ( xsd.include_citation_info() ) {
		citation_string = "\n\n" + ResidueSelectorFactory::get_instance()->get_citation_humanreadable( rs_type );
	}

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_residue_selector )
		.element_name( rs_type )
		.description( description + citation_string )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subelement, min_occurrence, max_occurrence )
		.write_complex_type_to_schema( xsd );

}

ResidueSelectorCOP
parse_residue_selector( utility::tag::TagCOP tag, basic::datacache::DataMap const & data, std::string const &option_name )
{
	std::string const selectorname = tag->getOption< std::string >( option_name, "" );
	if ( selectorname.empty() ) {
		TR.Warning << "Selector name is empty!" << std::endl;
		return ResidueSelectorCOP();
	}
	return get_residue_selector( selectorname, data );
}

/// @brief Companion function for parse_residue_selector
/// @brief This assumes the default residue selector option name ("residue_selector").
void
attributes_for_parse_residue_selector_default_option_name(
	utility::tag::AttributeList & attlist,
	std::string const & documentation_string
) {
	attributes_for_parse_residue_selector(attlist, "residue_selector", documentation_string);
}

std::string
parse_residue_selector_xsd_documentation() {
	return
		"The name of a previously declared residue selector or a logical"
		" expression of AND, NOT (!), OR, parentheses, and the names of"
		" previously declared residue selectors. Any capitalization of"
		" AND, NOT, and OR is accepted. An exclamation mark can be used"
		" instead of NOT. Boolean operators have their traditional priorities:"
		" NOT then AND then OR. For example, if selectors s1, s2, and s3 have been"
		" declared, you could write: 's1 or s2 and not s3' which would select"
		" a particular residue if that residue were selected by s1 or if it"
		" were selected by s2 but not by s3.";
}

void
attributes_for_parse_residue_selector(
	utility::tag::AttributeList & attlist,
	std::string const & option_name /* = "residue_selector" */,
	std::string const & documentation_string /* = "" */
)
{
	using namespace utility::tag;
	std::string user_docs = documentation_string;
	if ( ! utility::endswith(user_docs, ".") ) {
		user_docs+=".";
	}

	std::string const docs = user_docs +
		(user_docs == "" ? "" : " ") +
		parse_residue_selector_xsd_documentation();
	attlist + XMLSchemaAttribute( option_name, xs_string,  docs );
}

void
attributes_for_parse_residue_selector_when_required(
	utility::tag::AttributeList & attlist,
	std::string const & option_name /* = "residue_selector"*/,
	std::string const & documentation_string /* = ""*/
)
{
	using namespace utility::tag;
	std::string const docs = documentation_string +
		(documentation_string == "" ? "" : " ") +
		parse_residue_selector_xsd_documentation();
	attlist + XMLSchemaAttribute::required_attribute( option_name, xs_string, docs );
}

ResidueSelectorCOP
get_residue_selector( std::string const & selector_name, basic::datacache::DataMap const & data )
{

	core::select::residue_selector::ResidueSelectorCOP selector = nullptr;
	try {

		if ( data.has("ResidueSelector", selector_name) ) {
			selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} else {
			//JAB attempt to use logic
			TR << "Attempting to parse selector logic" << std::endl;
			SelectorLogicParser parser;
			selector =  parser.parse_string_to_residue_selector( data, selector_name );
			if ( TR.Debug.visible() ) {
				TR.Debug << "Parsed inline selector logic as:\n" << selector->debug_string() << std::endl;
			}
		}
	} catch ( utility::excn::Exception & e ) {
		std::stringstream error_msg;
		error_msg << "Failed to find ResidueSelector or parse logic of '" << selector_name << "' in the DataMap.\n";
		error_msg << e.msg();
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
	}
	if ( ! selector ) {
		std::stringstream error_msg;
		error_msg << "Failed to find ResidueSelector or parse logic of '" << selector_name << "' in the DataMap.\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
	}

	//TR << "Found residue selector " << selector_name << std::endl;
	return selector;
}

ResidueSelectorOP
get_embedded_residue_selector( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap ){
	if ( tag->size() > 1 ) {
		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		if ( tags.size() > 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  tag->getName() + " takes at most one ResidueSelector!\n" );
		}
		ResidueSelectorOP rs = ResidueSelectorFactory::get_instance()->new_residue_selector(
			tags.front()->getName(),
			tags.front(),
			datamap
		);
		return rs;
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "No selector embedded! " );
	}
}

utility::vector1< ResidueSelectorOP >
get_embedded_residue_selectors( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap ){
	if ( tag->size() > 1 ) {
		utility::vector1< ResidueSelectorOP > selectors;
		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		for ( core::Size i = 0; i <= tags.size() -1; ++i ) {
			ResidueSelectorOP rs = ResidueSelectorFactory::get_instance()->new_residue_selector(
				tags[i]->getName(),
				tags[i],
				datamap
			);
			selectors.push_back( rs );
		}
		return selectors;
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "No selector embedded! " );
	}
}

/// @brief Returns True if all the positions in the ResidueSubset are False
bool all_false_selection( ResidueSubset const & selection )
{
	return std::all_of( std::begin( selection ), std::end( selection ), []( bool i ) { return not i; } );;
}
/// @brief Returns True if all the positions in the ResidueSubset are True
bool all_true_selection( ResidueSubset const & selection )
{
	return std::all_of( std::begin( selection ), std::end( selection ), []( bool i ) { return i; } );
}
/// @brief Returns True if at least one position in the ResidueSubset is False
bool has_any_false_selection( ResidueSubset const & selection )
{
	return std::any_of( std::begin( selection ), std::end( selection ), []( bool i ) { return not i; } );;
}
/// @brief Returns True if at least one position in the ResidueSubset is True
bool has_any_true_selection( ResidueSubset const & selection )
{
	return std::any_of( std::begin( selection ), std::end( selection ), []( bool i ) { return i; } );
}
/// @brief Returns the number of selected residues in the ResidueSubset
core::Size count_selected( ResidueSubset const & selection )
{
	return std::count_if( std::begin( selection ), std::end( selection ), []( bool i ) { return i; } );
}
/// @brief Returns the Rosetta Numbering corresponding to the selected residues
utility::vector1< core::Size > selection_positions( ResidueSubset const & selection )
{
	auto it = std::find_if ( std::begin( selection ), std::end( selection ), []( bool i ) { return i; } );
	utility::vector1< core::Size > results;
	while ( it != std::end( selection ) ) {
		results.push_back( std::distance( std::begin( selection ), it ) + 1 );
		it = std::find_if( std::next( it ), std::end( selection ), []( bool i ) { return i; } );
	}
	return results;
}
/// @brief Returns the Rosetta Numbering corresponding to the unselected residues
utility::vector1< core::Size > unselection_positions( ResidueSubset const & selection )
{
	auto sele = []( bool i ) { return ! i; };
	auto it = std::find_if ( std::begin( selection ), std::end( selection ), sele );
	utility::vector1< core::Size > results;
	while ( it != std::end( selection ) ) {
		results.push_back( std::distance( std::begin( selection ), it ) + 1 );
		it = std::find_if( std::next( it ), std::end( selection ), sele );
	}
	return results;
}
/// @brief Evaluate if two ResidueSubsets are equal
bool are_selections_equal( ResidueSubset const & selection1, ResidueSubset const & selection2 )
{
	if ( selection1.size() != selection2.size() ) {
		return false;
	}
	return std::equal( selection1.begin(), selection1.end(), selection2.begin() );
}
/// @brief Returns a string representing the ResidueSubset (· for non selected, * for selected)
std::string represent_residue_selector( ResidueSubset const & selection, std::string const & is_true, std::string const & is_false )
{
	std::string output = "";
	for ( auto sele : selection ) {
		if ( sele ) {
			output.append( is_true );
		} else {
			output.append( is_false );
		}
	}
	return output;
}

ResidueSelectorOP
OR_combine( ResidueSelectorOP sele1, ResidueSelectorOP sele2 ) {
	using namespace core::select::residue_selector;
	if ( sele1 == nullptr ) { return sele2; }
	if ( sele2 == nullptr ) { return sele1; }

	OrResidueSelectorOP or_select( utility::pointer::dynamic_pointer_cast< OrResidueSelector >( sele1 ) );
	if ( or_select ) {
		// If we already have an Or selector, just add to it.
		or_select->add_residue_selector( sele2 );
		return or_select;
	} else {
		// If not, just combine with the Or selector
		return utility::pointer::make_shared< OrResidueSelector >( sele1, sele2 );
	}
}

ResidueSubset
OR_combine( ResidueSubset const & sele1, ResidueSubset const & sele2){
	debug_assert(sele1.size() == sele2.size());

	ResidueSubset new_subset(sele1.size() );
	for ( Size i = 1; i <= sele1.size(); ++i ) {
		new_subset[ i ] = sele1[ i ] || sele2[ i ];
	}
	return new_subset;
}

ResidueSelectorOP
AND_combine( ResidueSelectorOP sele1, ResidueSelectorOP sele2 ) {
	using namespace core::select::residue_selector;
	if ( sele1 == nullptr ) { return sele2; }
	if ( sele2 == nullptr ) { return sele1; }

	AndResidueSelectorOP and_select( utility::pointer::dynamic_pointer_cast< AndResidueSelector >( sele1 ) );
	if ( and_select ) {
		// If we already have an and selector, just add to it.
		and_select->add_residue_selector( sele2 );
		return and_select;
	} else {
		// If not, just combine with the Or selector
		return utility::pointer::make_shared< AndResidueSelector >( sele1, sele2 );
	}
}

ResidueSubset
AND_combine( ResidueSubset const & sele1, ResidueSubset const & sele2){

	debug_assert( sele1.size() == sele2.size() );

	ResidueSubset new_subset(sele1.size());
	for ( Size i = 1; i <= sele1.size(); ++i ) {
		new_subset[ i ] = sele1[ i ] && sele2[ i ];
	}
	return new_subset;
}

ResidueSelectorOP
parse_residue_selector_logic_string(
	basic::datacache::DataMap const & dm,
	utility::tag::TagCOP const & tag,
	std::string const & atname /*= "selector_logic" */
)
{
	if ( ! tag->hasOption( atname ) ) return nullptr;

	SelectorLogicParser parser;
	return parser.parse_string_to_residue_selector( dm, tag->getOption< std::string >( atname ) );
}

void
attributes_for_parse_residue_selector_logic_string(
	utility::tag::XMLSchemaDefinition & /*xsd*/, // not yet needed
	utility::tag::AttributeList & attributes,
	std::string const & selector_logic_attribute_name /*= "selector_logic"*/
)
{
	using namespace utility::tag;
	attributes + XMLSchemaAttribute( selector_logic_attribute_name, xs_string,
		"Logically combine already-delcared ResidueSelectors using boolean AND, OR, and ! (not)"
		" operators. As convnetional, ! (not) has the highest precedence, then AND, then OR."
		" Parentheses may be used to group operations together." );
}

// This is only usable if the regex works.
bool regex_usable() {
	std::string s("is the regex code working?");
	std::regex reg("regex");
	return std::regex_search(s,reg); // Problematic GCC versions return false for everything
}


utility::vector1<core::select::residue_selector::ResidueSubset>
identify_ss_blocks_vec( core::select::residue_selector::ResidueSubset const & selection ){
	utility::vector1<core::select::residue_selector::ResidueSubset> block_selections;
	core::select::residue_selector::ResidueSubset current_selection( selection.size(), false );

	//Note - does not take into account new chains!
	bool last_res = false;
	for ( core::Size count_resnum = 1; count_resnum <= selection.size(); ++count_resnum ) {
		if ( selection[count_resnum] ) {
			if ( !last_res ) {
				current_selection = core::select::residue_selector::ResidueSubset( selection.size(), false );
			}
			current_selection[count_resnum] = true;
		} else {
			if ( last_res ) {
				block_selections.push_back(current_selection);
			}
		}
		last_res = selection[count_resnum];
	}
	if ( last_res ) {
		block_selections.push_back(current_selection);
	}
	return block_selections;
}

bool
all_L_dssp(core::pose::Pose const & pose){
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.secstruct(i) != 'L' ) return false;
	}
	return true;
}

}
}
}
