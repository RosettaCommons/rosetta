// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/SliceResidueSelector.cc
/// @brief  The SliceResidueSelector allows slicing of the returned values of other residue selectors
/// @author Brian Coventry (bcov@uw.edu)

// Unit headers
#include <core/select/residue_selector/SliceResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.hh>

// Utility Headers
#include <utility/pointer/memory.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

#include <boost/format.hpp>

static basic::Tracer TR( "core.select.residue_selector.SliceResidueSelector" );

namespace core {
namespace select {
namespace residue_selector {


SliceResidueSelector::SliceResidueSelector() :
	selector_( nullptr ),
	slice_mode_( slice_enums::SPARSE ),
	from_( 0 ),
	to_( 0 ),
	indices_( ),
	oob_mode_( slice_enums::ERROR )
{}


SliceResidueSelector::SliceResidueSelector(
	ResidueSelectorCOP const & selector,
	slice_enums::SliceMode slice_mode,
	int from,
	int to,
	slice_enums::OutOfBoundsMode oob_mode /*= slice_enums::ERROR*/
) :
	selector_( selector ),
	slice_mode_( slice_mode ),
	from_( from ),
	to_( to ),
	indices_( ),
	oob_mode_( oob_mode )
{}

SliceResidueSelector::SliceResidueSelector(
	ResidueSelectorCOP selector,
	slice_enums::SliceMode slice_mode,
	utility::vector1< int > const & indices,
	slice_enums::OutOfBoundsMode oob_mode /*= slice_enums::ERROR*/
) :
	selector_( selector ),
	slice_mode_( slice_mode ),
	from_( 0 ),
	to_( 0 ),
	indices_( indices ),
	oob_mode_( oob_mode )
{}


/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP SliceResidueSelector::clone() const { return utility::pointer::make_shared<SliceResidueSelector>(*this); }


ResidueSubset
SliceResidueSelector::apply( core::pose::Pose const & pose ) const
{
	ResidueSubset initial_subset( pose.size(), true );

	if ( selector_ ) {
		initial_subset = selector_->apply( pose );
	}

	utility::vector1< core::Size > flat_indices;

	if ( slice_mode_ == slice_enums::SPARSE ) {
		flat_indices = get_residues_from_subset( initial_subset );
	} else if ( slice_mode_ == slice_enums::CONTIGUOUS ) {

		core::Size first = 0;
		core::Size last = 0;

		bool is_first = true;
		for ( Size i = 1; i <= pose.size(); i++ ) {
			if ( ! initial_subset[i] ) continue;

			if ( is_first ) {
				first = i;
				is_first = false;
			}
			last = i;
		}

		if ( first > 0 ) {
			flat_indices.resize( last - first + 1);
			for ( Size i = 0; i < flat_indices.size(); i++ ) {
				flat_indices[i+1] = i + first;
			}
		}
	} else {
		utility_exit();
	}

	utility::vector1< std::string > error_messages;

	utility::vector1< int > sliced_indices;
	if ( from_ > 0 || to_ > 0 ) {
		if ( indices_.size() > 0 ) {
			utility_exit_with_message("SliceResidueSelector: from-to and indices both set! Set only one of these!");
		}

		int use_from = std::max<int>( 1, wrap_index( from_, flat_indices.size(), error_messages ) );
		int use_to = std::min<int>( flat_indices.size(), wrap_index( to_, flat_indices.size(), error_messages ) );

		if ( use_from > use_to && ( oob_mode_ == slice_enums::ERROR || oob_mode_ == slice_enums::WARN ) ) {
			std::string message = "from set greater than to! from="
				+ utility::to_string(use_from) + " to=" + utility::to_string(use_to);
			error_messages.push_back( message );
		}
		for ( int i = use_from; i <= use_to; i++ ) {
			sliced_indices.push_back( flat_indices.at( i ) );
		}
	} else {
		for ( int index : indices_ ) {
			int wrapped = wrap_index( index, flat_indices.size(), error_messages );
			if ( wrapped <= 0 || wrapped > (int)flat_indices.size() ) continue;
			sliced_indices.push_back( flat_indices.at( wrapped ) );
		}
	}

	ResidueSubset final_selected( pose.size(), false );

	for ( int index : sliced_indices ) {
		final_selected[index] = true;
	}

	if ( TR.Debug.visible() ) {
		show_selection_logic( initial_subset, final_selected );
	}

	if ( error_messages.size() > 0 ) {
		error_messages.push_back("set \"-out:levels core.select.residue_selector.SliceResidueSelector:400\" to debug!");
		if ( oob_mode_ == slice_enums::ERROR ) {
			utility_exit_with_message( utility::join( error_messages, "\n" ) );
		}
		if ( oob_mode_ == slice_enums::WARN ) {
			for ( std::string const & msg : error_messages ) {
				TR.Warning << msg << std::endl;
			}
		}
	}

	return final_selected;
}

int
SliceResidueSelector::wrap_index( int index, int elements, utility::vector1< std::string > & error_messages ) const {

	int wrapped_value;
	if ( index < 0 ) {
		wrapped_value = elements + index + 1;

		if ( wrapped_value <= 0 && ( oob_mode_ == slice_enums::ERROR || oob_mode_ == slice_enums::WARN ) ) {
			std::string message = "wrapped index has gone below 1! Set="
				+ utility::to_string(index) + " wrapped=" + utility::to_string(wrapped_value);
			error_messages.push_back( message );
		}
	} else if ( index > 0 ) {

		wrapped_value = index;

		if ( wrapped_value > elements && ( oob_mode_ == slice_enums::ERROR || oob_mode_ == slice_enums::WARN ) ) {
			std::string message = "index has exceeded selection range! Set="
				+ utility::to_string(index) + " range_size=" + utility::to_string(elements);
			error_messages.push_back( message );
		}
	} else {

		wrapped_value = index;

		if ( oob_mode_ == slice_enums::ERROR || oob_mode_ == slice_enums::WARN ) {
			std::string message = "index=0 this never makes sense!";
			error_messages.push_back( message );
		}
	}

	return wrapped_value;
}

void
SliceResidueSelector::show_selection_logic(
	ResidueSubset const & initial,
	ResidueSubset const & final
) const {
	const Size width = 80;

	Size num_chars = utility::to_string( initial.size() ).length();
	// Whitespace padded number with same size as largest number
	std::string format_str = "%" + utility::to_string(num_chars) + "i";


	utility::vector1< std::string > number_bufs(num_chars, " Residue num: ");
	std::string initial_buf =                              "       Input: ";
	// std::string middle_buf =                               "  Raw select: ";
	std::string final_buf =                                "Final select: ";

	for ( Size i = 1; i <= initial.size(); i++ ) {
		std::string number = boost::str(boost::format( format_str)%i);
		debug_assert( number.size() == num_chars );

		for ( Size j = 1; j <= num_chars; j++ ) {
			number_bufs[j] += number[j-1];
		}

		initial_buf += initial[i] ? 'X' : '-';
		// middle_buf += middle[i] ? 'X' : '-';
		final_buf += final[i] ? 'X' : '-';
	}

	// The upperbound is ceil( initial_buf.size() / width )
	for ( Size i = 0; i < core::Size( (initial_buf.size() + width-1) / width ); i++ ) {

		Size start_char = width * i;
		TR << std::endl;
		for ( Size j = 1; j <= num_chars; j++ ) {
			TR << number_bufs[j].substr( start_char, width ) << std::endl;
		}

		TR.Debug << initial_buf.substr( start_char, width ) << std::endl;
		// TR.Debug << middle_buf.substr( start_char, width ) << std::endl;
		TR.Debug << final_buf.substr( start_char, width ) << std::endl;
	}
}

/// @brief Provide the citation.
/// @returns A vector of citation collections.  This allows the residue selector to provide citations for
/// itself and for any modules that it invokes.
/// @details This residue selector has no citation.  It may provide citations for the residue selector that
/// it calls, though.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::CitationCollectionCOP >
SliceResidueSelector::provide_citation_info() const {
	utility::vector1< basic::citation_manager::CitationCollectionCOP > returnvec;
	if ( selector_ != nullptr ) {
		basic::citation_manager::merge_into_citation_collection_vector( selector_->provide_citation_info(), returnvec );
	}
	return returnvec;
}

/// @brief Does this residue selector indicate that it is unpublished (and, by extension, that the author should be
/// included in publications resulting from it)?
/// @details Returns true (this residue selector is unpublished).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
SliceResidueSelector::residue_selector_is_unpublished() const {
	return true;
}

/// @brief Provide a list of authors and their e-mail addresses, as strings.
/// @returns This residue selector was created by Brian Coventry.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
SliceResidueSelector::provide_authorship_info_for_unpublished() const {
	using namespace basic::citation_manager;
	utility::vector1< UnpublishedModuleInfoCOP > returnvec{
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		"SliceResidueSelector", CitedModuleType::ResidueSelector,
		"Brian Coventry",
		"Baker Laboratory, Institute for Protein Design, Dept. of Biochemistry, University of Washington",
		"bcov@uw.edu"
		)
		};
	if ( selector_ != nullptr ) {
		merge_into_unpublished_collection_vector( selector_->provide_authorship_info_for_unpublished(), returnvec );
	}
	return returnvec;
}

void SliceResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	if ( tag->hasOption( "selector" ) ) { // fetch selector from datamap

		if ( tag->size() > 1 ) { // has subtags
			throw CREATE_EXCEPTION(utility::excn::Exception,  "SliceResidueSelector can slice ONE ResidueSelector! Either specify 'selector' option or provide subtags but not BOTH\n" );
		}
		// grab the ResidueSelector to be sliced from the selector option
		// and then grab each of the indicated residue selectors from the datamap.
		std::string selector_str;
		try {
			selector_str = tag->getOption< std::string >( "selector" );
		} catch ( utility::excn::Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access required option 'selector' from SliceResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}

		try {
			ResidueSelectorCOP selector = get_residue_selector(selector_str, datamap);
			set_residue_selector(selector);
		} catch ( utility::excn::Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selector_str << "' from the Datamap from SliceResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}
	} else if ( tag->size() > 1 ) { // attempt reading subtag
		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		if ( tags.size() > 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "SliceResidueSelector takes exactly ONE ResidueSelector! Multiple selectors were specified.\n" );
		}
		ResidueSelectorCOP rs = ResidueSelectorFactory::get_instance()->new_residue_selector(
			tags.front()->getName(),
			tags.front(),
			datamap
		);
		set_residue_selector( rs );
	}

	if ( tag->hasOption( "from" ) ) {
		set_from( tag->getOption< int >( "from" ) );
	}
	if ( tag->hasOption( "to" ) ) {
		set_to( tag->getOption< int >( "to" ) );
	}
	if ( tag->hasOption( "indices" ) ) {
		utility::vector1< int > indices = utility::string_split<int>( tag->getOption< std::string >( "indices" ), ',', int(0) );
		set_indices( indices );
	}
	if ( tag->hasOption( "slice_mode" ) ) {
		std::string str_slice = tag->getOption< std::string >( "slice_mode" );
		slice_enums::SliceMode slice_mode = slice_enums::SPARSE;
		if ( str_slice == "SPARSE" ) {
			slice_mode = slice_enums::SPARSE;
		} else if ( str_slice == "CONTIGUOUS" || str_slice == "CONTIGUOUS_OR" ) {
			slice_mode = slice_enums::CONTIGUOUS;
		} else {
			utility_exit_with_message( "SliceResidueSelector: Unknown slice_mode: " + str_slice
				+ " Available options are: SPARSE, CONTIGUOUS");
		}
		set_slice_mode( slice_mode );
	}
	if ( tag->hasOption( "oob_mode" ) ) {
		std::string str_mode = tag->getOption< std::string >( "slice_mode" );
		slice_enums::OutOfBoundsMode oob_mode = slice_enums::ERROR;
		if ( str_mode == "ERROR" ) {
			oob_mode = slice_enums::ERROR;
		} else if ( str_mode == "WARN" ) {
			oob_mode = slice_enums::WARN;
		} else if ( str_mode == "IGNORE" ) {
			oob_mode = slice_enums::IGNORE;
		} else {
			utility_exit_with_message( "SliceResidueSelector: Unknown oob_mode: " + str_mode
				+ " Available options are: ERROR, WARN, IGNORE");
		}
		set_out_of_bounds_behavior( oob_mode );
	}

}

void SliceResidueSelector::set_residue_selector( ResidueSelectorCOP selector )
{
	selector_ = selector;
}


void
SliceResidueSelector::set_slice_mode( slice_enums::SliceMode slice_mode )
{
	slice_mode_ = slice_mode;
}

void
SliceResidueSelector::set_from( int from )
{
	from_ = from;
}

void
SliceResidueSelector::set_to( int to )
{
	to_ = to;
}

void
SliceResidueSelector::set_indices( utility::vector1< int > const & indices )
{
	indices_ = indices;
}

void
SliceResidueSelector::set_out_of_bounds_behavior( slice_enums::OutOfBoundsMode oob_mode )
{
	oob_mode_ = oob_mode;
}


std::string SliceResidueSelector::get_name() const {
	return SliceResidueSelector::class_name();
}

std::string SliceResidueSelector::class_name() {
	return "Slice";
}

void
SliceResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attributes_for_parse_residue_selector(attlist, "selector" , "The selector to slice from.");
	attlist

		+ XMLSchemaAttribute::attribute_w_default( "from" , xs_integer , "Range selection: This is the first residue of the range to select." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "to" , xs_integer , "Range selection: This is the last residue of the range to select." , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "indices" , xs_string , "Index selection: Comma separated list of indices. May not use this with from-to" , "" )
		+ XMLSchemaAttribute::attribute_w_default( "slice_mode" , xs_string , "How should the previous residue selector be represented? SPARSE: all gaps are removed from the previous selection. If residues 101 and 200 are selected, only indices 1 and 2 are valid. "
		"CONTIGUOUS: gaps from previous selection are included in indices. If residues 101 and 200 are selected, all indices from 1 to 100 are valid.", "SPARSE" )
		+ XMLSchemaAttribute::attribute_w_default( "oob_mode" , xs_string , "If an index is out of bounds, how should this be handled? ERROR: Quit and display error message. WARN: Display error message. IGNORE: Ignore." , "ERROR" );

	xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(),"Residue selector that allows slicing of the selection of other residue selections. Also allows negative indexing. 1 is the first residue and -1 is the last.", attlist );
}


ResidueSelectorOP
SliceResidueSelectorCreator::create_residue_selector() const {
	return utility::pointer::make_shared<SliceResidueSelector>();
}

std::string
SliceResidueSelectorCreator::keyname() const {
	return SliceResidueSelector::class_name();
}

void
SliceResidueSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SliceResidueSelector::provide_xml_schema( xsd );
}


} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::SliceResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( selector_ ) );
	arc( CEREAL_NVP( slice_mode_ ) );
	arc( CEREAL_NVP( from_ ) );
	arc( CEREAL_NVP( to_ ) );
	arc( CEREAL_NVP( indices_ ) );
	arc( CEREAL_NVP( oob_mode_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::SliceResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector ); // ResidueSelectorCOP
	selector_ = local_selector; // copy the non-const pointer(s) into the const pointer(s)
	arc( slice_mode_ );
	arc( from_ );
	arc( to_ );
	arc( indices_ );
	arc( oob_mode_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::SliceResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::SliceResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_SliceResidueSelector )
#endif // SERIALIZATION
