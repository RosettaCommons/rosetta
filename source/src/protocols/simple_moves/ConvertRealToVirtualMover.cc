// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/simple_moves/ConvertRealToVirtualMover.cc
/// @brief Mover for switching a residue type to all virtual
/// @author Sebastian Raemisch (raemisch@scripps.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/simple_moves/ConvertRealToVirtualMover.hh>
#include <protocols/simple_moves/ConvertRealToVirtualMoverCreator.hh>


// Core headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.simple_moves.ConvertRealToVirtualMover" );

namespace protocols {
namespace simple_moves {

using namespace core;

// Default constructor
ConvertRealToVirtualMover::ConvertRealToVirtualMover():
	protocols::moves::Mover( ConvertRealToVirtualMover::mover_name() ),
	selector_( nullptr )
{

}

ConvertRealToVirtualMover::ConvertRealToVirtualMover(core::select::residue_selector::ResidueSelectorCOP selector):
	protocols::moves::Mover( ConvertRealToVirtualMover::mover_name() ),
	selector_(selector)
{

}

// Destructor
ConvertRealToVirtualMover::~ConvertRealToVirtualMover(){}

ConvertRealToVirtualMover::ConvertRealToVirtualMover( ConvertRealToVirtualMover const & src ):
	Mover( src )
{
	using namespace core::select::residue_selector;
	if ( src.selector_ ) selector_ = src.selector_->clone();
}

void
ConvertRealToVirtualMover::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}


// Parse my tag
void
ConvertRealToVirtualMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( tag->hasOption( "residue_selector" ) ) {
		// set the selector_ private variable
		selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data );
		if ( !selector_ ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "ResidueSelector passed to FaToVirtual mover could not be found." );
		}
	}
}

void
ConvertRealToVirtualMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	protocols::rosetta_scripts::attributes_for_parse_residue_selector(attributes);
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),

		"A mover for switching a residue type to all VIRTUAL."
		"A VIRTUAL Residue is one that is not scored or output.",

		attributes );
}


// Clone
protocols::moves::MoverOP
ConvertRealToVirtualMover::clone() const
{
	return protocols::moves::MoverOP( new ConvertRealToVirtualMover( *this ) );
}

// Fresh instance
protocols::moves::MoverOP
ConvertRealToVirtualMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ConvertRealToVirtualMover );
}

std::string
ConvertRealToVirtualMover::get_name() const
{
	return ConvertRealToVirtualMover::mover_name();
}

std::string
ConvertRealToVirtualMover::mover_name()
{
	return "ConvertRealToVirtualMover";
}

void
ConvertRealToVirtualMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, ConvertRealToVirtualMover const & mover )
{
	mover.show(os);
	return os;
}


///////////  APPLY  ////////////

void
ConvertRealToVirtualMover::apply( core::pose::Pose & pose )
{

	/// Get a ResidueSubset (list of bools) from ResidueSelector
	core::select::residue_selector::ResidueSubset subset;
	if ( selector_ ) {
		subset = selector_->apply( pose );
	} else {
		subset.resize( pose.total_residue(), true );
		TR << "Setting all residues to Virtual." << std::endl;
	}

	/// Go through the subset list and make VIRT if suset[i]=1
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( !subset[i] ) continue; //Skip residues masked by the ResidueSelector.
		TR.Debug << "Converting " << i << " to virtual " << std::endl;
		pose.real_to_virtual( i );
	}
}


/////////////// Creator ///////////////

protocols::moves::MoverOP
ConvertRealToVirtualMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new ConvertRealToVirtualMover );
}





std::string
ConvertRealToVirtualMoverCreator::keyname() const
{
	return ConvertRealToVirtualMover::mover_name();
}

void
ConvertRealToVirtualMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConvertRealToVirtualMover::provide_xml_schema( xsd );
}

} //protocols
} //simple_moves

