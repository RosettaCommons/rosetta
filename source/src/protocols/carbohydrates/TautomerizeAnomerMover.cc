// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/carbohydrates/TautomerizeAnomerMover.cc
/// @brief   Method definitions for TautomerizeAnomerMover.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <protocols/carbohydrates/TautomerizeAnomerMover.hh>
#include <protocols/carbohydrates/TautomerizeAnomerMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/exit.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Basic header
#include <basic/Tracer.hh>


// Construct tracers.
static basic::Tracer TR( "protocols.carbohydrates.TautomerizeAnomerMover" );


namespace protocols {
namespace carbohydrates {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
/// @brief  Default constructor
TautomerizeAnomerMover::TautomerizeAnomerMover() : Mover()
{
	init();
}

// Copy constructor
TautomerizeAnomerMover::TautomerizeAnomerMover( TautomerizeAnomerMover const & object_to_copy ) :
	Mover( object_to_copy )
{
	copy_data( *this, object_to_copy );
}

/// @brief  Constructor with ResidueSelector input option
TautomerizeAnomerMover::TautomerizeAnomerMover( core::select::residue_selector::ResidueSelectorCOP selector ) :
	Mover(),
	selector_( selector )
{
	init();
}

// Assignment operator
TautomerizeAnomerMover &
TautomerizeAnomerMover::operator=( TautomerizeAnomerMover const & object_to_copy )
{
	// Abort self-assignment.
	if ( this != &object_to_copy ) {
		Mover::operator=( object_to_copy );
		copy_data( *this, object_to_copy );
	}
	return *this;
}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
/// @brief  Generate string representation of TautomerizeAnomerMover for debugging purposes.
void
TautomerizeAnomerMover::show( std::ostream & output ) const
{
	Mover::show( output );  // name, type, tag
}


// Mover methods
/// @brief  Register options with the option system.
void
TautomerizeAnomerMover::register_options()
{
	// Call register_options() on all other Movers used by this class.
	Mover::register_options();  // Mover's register_options() doesn't do anything; it's just here in principle.
}

// Return the string identifier for the associated Mover (TautomerizeAnomerMover).
std::string
TautomerizeAnomerMover::get_name() const {
	return mover_name();
}

protocols::moves::MoverOP
TautomerizeAnomerMover::clone() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< TautomerizeAnomerMover >( *this ) );
}

protocols::moves::MoverOP
TautomerizeAnomerMover::fresh_instance() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< TautomerizeAnomerMover >() );
}

void
TautomerizeAnomerMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const & )
{
	// Parse the ResidueSelector tag.
	if ( tag->hasOption( "residue_selector" ) ) {
		selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data );
	}
}

void
TautomerizeAnomerMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( []( std::string const & name ) {
		return "TautomerizeAnomerMover_subelement_" + name + "Type";
	} );

	AttributeList attlist;

	rosetta_scripts::attributes_for_parse_residue_selector( attlist, "The name of a pre-defined ResidueSelector." );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"This carbohydrate-specific Mover randomly selects a free reducing end (not a glycoside) and inverts the "
		"stereochemistry, swapping alpha anomers for beta and beta for alpha.  (This could be considered an extremely "
		"limited design case; however, reducing ends readily tautomerize in solution, in contrast to other cases, in "
		"which residues do not readily mutate into others!)  It is generally not certain which form is preferred (if "
		"any) in sugar-binding proteins, and crystal structures sometimes arbitrarily assign one anomer over another "
		"when fitting density, so this Mover can assure that each anomer is sampled.",
		attlist, subelements );
}


/// @param  <input_pose>: the structure to be moved
void
TautomerizeAnomerMover::apply( Pose & input_pose )
{
	using namespace std;

	show( TR );

	TR << "Getting movable residues...." << endl;

	setup_movable_reducing_ends( input_pose );

	if ( movable_reducing_ends_.empty() ) {
		TR.Warning << "There are no movable residues available in the given pose." << endl;
		set_last_move_status( moves::FAIL_DO_NOT_RETRY );
		return;
	}

	TR << "Applying " << get_name() << " to pose...." << endl;

	core::uint const i( core::uint( numeric::random::rg().uniform() * movable_reducing_ends_.size() + 1 ) );

	TR << "Changing anomeric form of residue " << i << endl;

	core::pose::carbohydrates::tautomerize_anomer( input_pose, i );

	TR << "Move complete." << endl;
}


// Private methods ////////////////////////////////////////////////////////////
// Initialize data members from arguments.
void
TautomerizeAnomerMover::init()
{
	type( mover_name() );
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
TautomerizeAnomerMover::copy_data(
	TautomerizeAnomerMover & object_to_copy_to,
	TautomerizeAnomerMover const & object_to_copy_from )
{
	object_to_copy_to.selector_ = object_to_copy_from.selector_;
	object_to_copy_to.movable_reducing_ends_ = object_to_copy_from.movable_reducing_ends_;
}


// Setup list of movable residues from the ResidueSelector.
void
TautomerizeAnomerMover::setup_movable_reducing_ends( core::pose::Pose const & pose )
{
	using namespace core;
	using namespace select::residue_selector;
	using namespace conformation;

	movable_reducing_ends_.clear();

	conformation::Conformation const & conf( pose.conformation() );
	if ( ! conf.contains_carbohydrate_residues() ) {
		TR.Warning << "TautomerizeAnomerMover is a carbohydrate-specific Mover; ";
		TR.Warning << "this Pose contains no carbohydrate residues." << std::endl;
		return;
	}

	Size const n_res( pose.size() );

	ResidueSubset subset;
	if ( selector_ ) {
		subset = selector_->apply( pose );
	} else {
		subset = ResidueSubset( n_res, true );  // Assume that all residues can be selected.
	}

	for ( core::uint i( 1 ); i <= n_res; ++i ) {
		if ( ! subset[ i ] ) { continue; }  // Skip residues masked by the ResidueSelector.

		Residue const & res( pose.residue( i ) );
		if ( ! res.is_lower_terminus() ) { continue; }  // All reducing ends are lower termini.
		if ( ! res.is_carbohydrate() ) { continue; }
		if ( res.carbohydrate_info()->is_glycoside() ) {
			// Glycosides are chemically stable and do not tautomerize.
			continue;
		}

		// If we got this far, we are a reducing end.
		movable_reducing_ends_.push_back( i );
	}
}


// Creator methods ////////////////////////////////////////////////////////////
std::string
TautomerizeAnomerMoverCreator::keyname() const
{
	return TautomerizeAnomerMover::mover_name();
}

// Return an up-casted owning pointer (MoverOP) to the mover.
protocols::moves::MoverOP
TautomerizeAnomerMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< TautomerizeAnomerMover >() );
}

void
TautomerizeAnomerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TautomerizeAnomerMover::provide_xml_schema( xsd );
}


// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that TautomerizeAnomerMover can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, TautomerizeAnomerMover const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace carbohydrates
}  // namespace protocols
