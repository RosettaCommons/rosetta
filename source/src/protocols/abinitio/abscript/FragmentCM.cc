// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/FragmentCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/FragmentCM.hh>
#include <protocols/abinitio/abscript/FragmentCMCreator.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <core/environment/LocalPosition.hh>
#include <protocols/environment/DofUnlock.hh>

#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/CutBiasClaim.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/select/residue_selector/ResidueSelector.hh>

#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>
#include <protocols/simple_moves/GunnCost.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>

// utility/basic headers
#include <utility/tag/Tag.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

#include <basic/datacache/DataMap.hh>

#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer tr( "protocols.environment.movers.FragmentCM", basic::t_info );

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace environment;


// creator
// XRW TEMP std::string
// XRW TEMP FragmentCMCreator::keyname() const {
// XRW TEMP  return FragmentCM::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FragmentCMCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new FragmentCM );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FragmentCM::mover_name() {
// XRW TEMP  return "FragmentCM";
// XRW TEMP }


FragmentCM::FragmentCM():
	ClientMover(),
	mover_( /* NULL */ ),
	selector_( /* NULL */ ),
	bInitialize_( true ),
	bYieldCutBias_( false )
{}

FragmentCM::FragmentCM( simple_moves::FragmentMoverOP mover,
	core::select::residue_selector::ResidueSelectorCOP selector ):
	ClientMover(),
	selector_( selector ),
	bInitialize_( true ),
	bYieldCutBias_( false )
{
	set_mover( mover );
}

FragmentCM::~FragmentCM() {}

void FragmentCM::set_selector( core::select::residue_selector::ResidueSelectorCOP selector ) {
	if ( Parent::state_check( __FUNCTION__, ( selector.get() == selector_.get() ) ) ) {
		selector_ = selector;
	}
}

void FragmentCM::set_mover( simple_moves::FragmentMoverOP mover ){

	// commented out because it doesn't play nice with subclasses making settings during
	// claiming. It would be nice to think of a way to automatically check for accidental
	// post-claming changes to the mover, but I'm not sure how to do that.

	// if( Parent::state_check( __FUNCTION__, ( mover().get() == mover_.get() ) ) ){
	mover_ = mover;
	type( get_name() );
	//}
}

core::Size determine_frag_size( std::string const& file ) {
	using namespace core::fragment;
	using namespace basic::options;

	FragmentIO frag_io( option[ OptionKeys::abinitio::number_3mer_frags ](), 1,
		option[ OptionKeys::frags::annotate ]() );

	FragSetOP fragset = frag_io.read_data( file );

	// Assumes constant length fragments. I don't think anyone uses these (or even knows its possible)
	// so this seems safe.
	return fragset->max_frag_length();
}

void FragmentCM::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const&,
	protocols::moves::Movers_map const&,
	core::pose::Pose const& ) {
	using namespace core::fragment;
	using namespace basic::options;

	std::string frag_file_tag = "fragments";

	core::Size const frag_length = determine_frag_size( tag->getOption< std::string >( frag_file_tag ) );

	// The number of frags to be randomly selected out of the total fragments. I think.
	core::Size n_frags;
	if ( frag_length == 3 ) {
		n_frags = tag->getOption< core::Size >( "nfrags", option[ OptionKeys::abinitio::number_3mer_frags ]() );
	} else if ( frag_length == 9 ) {
		n_frags = tag->getOption< core::Size >( "nfrags", option[ OptionKeys::abinitio::number_9mer_frags ]() );
	} else {
		n_frags = tag->getOption< core::Size >( "nfrags" );
	}

	FragmentIO frag_io( n_frags, 1, option[ OptionKeys::frags::annotate ]() );

	std::string const& frag_type = tag->getOption< std::string >( "frag_type", "classic" );
	if ( frag_type == "classic" ) {
		set_mover( simple_moves::FragmentMoverOP( new simple_moves::ClassicFragmentMover( frag_io.read_data( tag->getOption< std::string >( frag_file_tag ) ) ) ) );
	} else if ( frag_type == "smooth" ) {
		set_mover( simple_moves::FragmentMoverOP( new simple_moves::SmoothFragmentMover( frag_io.read_data( tag->getOption< std::string >( frag_file_tag ) ), protocols::simple_moves::FragmentCostOP( new simple_moves::GunnCost() ) ) ) );
	} else {
		std::ostringstream ss;
		ss << "The fragment type " << frag_type << " is not valid. The options "
			<< " are 'classic' and 'smooth'.";
		throw utility::excn::EXCN_RosettaScriptsOption( ss.str() );
	}

	initialize( tag->getOption< bool >( "initialize", true ) );
	yield_cut_bias(tag->getOption< bool >("yield_cut_bias", false));

	set_selector( datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", tag->getOption<std::string>( "selector" ) ) );
}

claims::EnvClaims FragmentCM::yield_claims( core::pose::Pose const& pose,
	basic::datacache::WriteableCacheableMapOP ){
	using namespace claims;
	claims::EnvClaims claim_list;

	if ( yield_cut_bias() ) {
		core::fragment::SecondaryStructureOP ss( new core::fragment::SecondaryStructure( *( mover()->fragments() ) ) );
		claim_list.push_back( protocols::environment::claims::EnvClaimOP( new environment::claims::CutBiasClaim(
			utility::pointer::static_pointer_cast< ClientMover > ( get_self_ptr() ),
			"BASE",
			*ss ) ) );
	}

	if ( selector() ) {
		utility::vector1< bool > torsion_mask;
		try {
			torsion_mask = selector()->apply( pose );
		} catch( utility::excn::EXCN_Msg_Exception& e ){
			std::ostringstream ss;
			ss << this->get_name() << " failed to apply its ResidueSelector in " << __FUNCTION__ << ".";
			e.add_msg(ss.str());
			throw;
		}
		int shift = torsion_mask.index( true )-1;

		if ( shift == -1 ) { // verify that torsion_mask isn't all-false.
			std::ostringstream ss;
			ss << "Selector used address fragment insertions of '" << this->get_name() << "' was empty.";
			throw utility::excn::EXCN_BadInput( ss.str() );
		}

		mover()->set_fragments( mover()->fragments()->clone_shifted( shift ) );
	}

	TorsionClaimOP new_claim( new TorsionClaim( utility::pointer::static_pointer_cast< ClientMover > ( get_self_ptr() ),
		"BASE",
		std::make_pair( mover()->fragments()->min_pos(),
		mover()->fragments()->max_pos() ) ) );

	if ( initialize() ) {
		new_claim->strength( CAN_CONTROL, CAN_CONTROL );
	} else {
		new_claim->strength( CAN_CONTROL, DOES_NOT_CONTROL );
	}
	claim_list.push_back( new_claim );

	return claim_list;
}

void FragmentCM::initialize( Pose& pose ){
	DofUnlock activation( pose.conformation(), passport() );
	mover()->apply_at_all_positions( pose );
}

void FragmentCM::apply( Pose& pose ) {
	if ( pose.conformation().is_protected() ) {
		DofUnlock activation( pose.conformation(), passport() );
		mover()->apply( pose );
	} else {
		mover()->apply( pose );
	}
}

void FragmentCM::initialize( bool setting ) {
	if ( Parent::state_check( __FUNCTION__, bInitialize_ == setting ) ) {
		bInitialize_ = setting;
	}
}

void FragmentCM::yield_cut_bias( bool setting ){
	if ( Parent::state_check( __FUNCTION__, bYieldCutBias_ == setting ) ) {
		bYieldCutBias_ = setting;
	}
}

void FragmentCM::passport_updated() {
	if ( mover()==NULL ) {
		tr.Warning << "error in " << get_name() << " mover is null or not cofigured " << std::endl;
		return;
	}
	debug_assert( mover() );
	if ( has_passport() ) {
		mover()->set_movemap( passport()->render_movemap() );
	} else {
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
		mm->set_bb( true );
		mover()->set_movemap( mm );
	}
}

std::string FragmentCM::get_name() const {
	return mover_name();
}

std::string FragmentCM::mover_name() {
	return "FragmentCM";
}

void FragmentCM::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"fragments", xs_string,
		"Path to fragment file to use in this mover");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"nfrags", xsct_non_negative_integer,
		"Number of fragments per position."
		"The number of 3mers or 9mers can be specified using the command line.", "25");

	XMLSchemaRestriction frag_type_enum;
	frag_type_enum.name("frag_type_name");
	frag_type_enum.base_type(xs_string);
	frag_type_enum.add_restriction(xsr_enumeration, "classic");
	frag_type_enum.add_restriction(xsr_enumeration, "smooth");
	xsd.add_top_level_element(frag_type_enum);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"frag_type", "frag_type_name", "Should classic or smooth fragment insertion be performed?", "classic");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"initialize", xsct_rosetta_bool,
		"Should the mover insert a fragment at every position after brokering is complete?",
		"true");

	attlist + XMLSchemaAttribute::required_attribute(
		"selector", xs_string,
		"Residue selector which region this claim mover will be applied to");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Fragment Claim Mover performs fragment insertion on a specified region of the pose",
		attlist );
}

std::string FragmentCMCreator::keyname() const {
	return FragmentCM::mover_name();
}

protocols::moves::MoverOP
FragmentCMCreator::create_mover() const {
	return protocols::moves::MoverOP( new FragmentCM );
}

void FragmentCMCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FragmentCM::provide_xml_schema( xsd );
}


} // abscript
} // abinitio
} // protocols
