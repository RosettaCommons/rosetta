// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/residue_selectors/NativeSelector.hh
/// @brief  A ResidueSelector that applies a given residue selector to the native pose.
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit headers
#include <protocols/residue_selectors/NativeSelector.hh>
#include <protocols/residue_selectors/NativeSelectorCreator.hh>

//Protocols
#include <protocols/rosetta_scripts/util.hh>

// Core Headers
#include <core/pose/selection.hh>
#include <core/pose/Pose.hh>
#include <core/pose/ref_pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh>
#include <core/types.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <utility/pointer/memory.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace residue_selectors {

static basic::Tracer TR( "protocols.residue_selectors.NativeSelector" );

NativeSelector::NativeSelector():
	selector_( 0 ),
	native_( 0 )
{}

NativeSelector::NativeSelector( NativeSelector const & src ){
	selector_ = src.selector_->clone();
	native_ = src.native_;
}

NativeSelector::NativeSelector(
	core::select::residue_selector::ResidueSelectorCOP inner_selector
):
	selector_( std::move( inner_selector ) ),
	native_( 0 )
{}

NativeSelector::~NativeSelector() = default;

core::select::residue_selector::ResidueSelectorOP
NativeSelector::clone() const{
	return utility::pointer::make_shared< NativeSelector >( * this );
}

core::select::residue_selector::ResidueSubset
NativeSelector::apply( core::pose::Pose const & pose ) const{
	using namespace basic::options;

	runtime_assert( selector_ );

	if ( ! native_ ) {
		if ( option[ OptionKeys::in::file::native ].user() ) {
			native_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ].value() );
		} else {
			utility_exit_with_message( "No Native Pose Found By NativeSelector" );
		}
	}

	core::select::residue_selector::ResidueSubset return_val = selector_->apply( * native_ );

	//resize if necessary!
	if ( return_val.size() > pose.size() ) {
		//If too large, remove elements from the end
		TR << "Warning, native pose has more residues than input pose. We are deleting the final "
			<< (return_val.size() - pose.size()) << " elements of the residue subset." << std::endl;
		return_val.resize( pose.size(), false );
	} else if ( return_val.size() < pose.size() ) {
		//If too small, append falses to the end
		auto const num_elements_to_add = pose.size() - return_val.size();
		TR << "Warning, native pose has less residues than input pose. We are appending "
			<< num_elements_to_add << " 'false' elements to the end of the residue subset." << std::endl;

#ifndef NDEBUG
		//make sure all original values are unchanged (this behavior changed with c++11 so it's good to check)
		core::select::residue_selector::ResidueSubset copy = return_val;
#endif
		return_val.resize( pose.size(), false );
#ifndef NDEBUG
		//make sure all original values are unchanged (this behavior changed with c++11 so it's good to check)
		for ( core::Size ii = 1; ii <= copy.size(); ++ii ) {
			debug_assert( copy[ ii ] == return_val[ ii ] );
		}
#endif
	}

	return return_val;
}

void
NativeSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	if ( tag->hasOption("residue_selector") ) {
		// grab the ResidueSelector from the selector option
		std::string selector_str;
		try {
			selector_str = tag->getOption< std::string >( "residue_selector" );
		} catch ( utility::excn::Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access option 'selector' from NativeSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}
		try {
			core::select::residue_selector::ResidueSelectorCOP selector =
				datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_str );
			set_residue_selector( selector );
		} catch ( utility::excn::Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selector_str << "' from the Datamap from NativeSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}
	}

	if ( datamap.has_resource( "native_pose" ) ) {
		//This overrides the -native flag
		set_native_pose( core::pose::saved_native_pose( datamap ) );
	}

}
///@details Non-static function to return class name
std::string
NativeSelector::get_name() const{
	return NativeSelector::class_name();
}

/// @details Static function to return class name; anything else returning the class name should call this function to avoid discrepancies
std::string
NativeSelector::class_name(){
	return "NativeSelector";
}

/// @details Static function allowing evaluation of XML before parsing
void
NativeSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "residue_selector", xs_string, "Name of residue selector to be applied to the native pose." );

	core::select::residue_selector::xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name() ,"A ResidueSelector that applies a given residue selector to the native pose. If the native pose is shorter than the given pose, 'false' values will be appended onto the end of the returned residue subset. If the native pose is longer than the given pose, then the subset will be shortened to be the length of the given pose.", attlist );

}


//CREATOR METHODS
core::select::residue_selector::ResidueSelectorOP
NativeSelectorCreator::create_residue_selector() const{
	return utility::pointer::make_shared< NativeSelector >();
}

std::string
NativeSelectorCreator::keyname() const{
	return NativeSelector::class_name();
}

void
NativeSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	NativeSelector::provide_xml_schema( xsd );
}

} //protocols
} //residue_selectors


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::residue_selectors::NativeSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( selector_ ) );
	arc( CEREAL_NVP( native_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::residue_selectors::NativeSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( selector_ );

	core::pose::PoseOP temp;
	arc( temp );
	native_ = temp;
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::residue_selectors::NativeSelector );
CEREAL_REGISTER_TYPE( protocols::residue_selectors::NativeSelector )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_residue_selectors_NativeSelector )
#endif // SERIALIZATION
