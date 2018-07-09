// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/carbohydrates/RingPlaneFlipMover.cc
/// @brief   Method definitions for RingPlaneFlipMover.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <protocols/carbohydrates/RingPlaneFlipMover.hh>
#include <protocols/carbohydrates/RingPlaneFlipMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/movemap/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/chemical/rings/AxEqDesignation.hh>
#include <core/chemical/rings/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/pose/Pose.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Basic header
//#include <basic/options/option.hh>
#include <basic/Tracer.hh>


// Construct tracers.
static basic::Tracer TR( "protocols.carbohydrates.RingPlaneFlipMover" );


namespace protocols {
namespace carbohydrates {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
/// @brief  Default constructor
RingPlaneFlipMover::RingPlaneFlipMover() : Mover()
{
	init();
}

// Copy constructor
RingPlaneFlipMover::RingPlaneFlipMover( RingPlaneFlipMover const & object_to_copy ) : Mover( object_to_copy )
{
	copy_data( *this, object_to_copy );
}

/// @brief  Constructor with MoveMap input option
RingPlaneFlipMover::RingPlaneFlipMover( core::kinematics::MoveMapOP input_movemap ) :
	Mover(),
	movemap_( input_movemap )
{
	init();
}

/// @brief  Constructor with ResidueSelector input option
RingPlaneFlipMover::RingPlaneFlipMover( core::select::residue_selector::ResidueSelectorCOP selector ) :
	Mover(),
	selector_( selector )
{
	init();
}

// Assignment operator
RingPlaneFlipMover &
RingPlaneFlipMover::operator=( RingPlaneFlipMover const & object_to_copy )
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
/// @brief  Generate string representation of RingPlaneFlipMover for debugging purposes.
void
RingPlaneFlipMover::show( std::ostream & output ) const
{
	Mover::show( output );  // name, type, tag
}

// Mover methods
/// @brief  Register options with the option system.
void
RingPlaneFlipMover::register_options()
{
	// Call register_options() on all other Movers used by this class.
	Mover::register_options();  // Mover's register_options() doesn't do anything; it's just here in principle.
}

// Return the string identifier for the associated Mover (RingPlaneFlipMover).
std::string
RingPlaneFlipMover::get_name() const {
	return mover_name();
}

protocols::moves::MoverOP
RingPlaneFlipMover::clone() const
{
	return protocols::moves::MoverOP( new RingPlaneFlipMover( *this ) );
}

protocols::moves::MoverOP
RingPlaneFlipMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new RingPlaneFlipMover() );
}

void
RingPlaneFlipMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const & )
{
	using namespace core::select::movemap;

	// Parse the MoveMapFactory tag.
	MoveMapFactoryOP mmf_legacy( protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data ) );
	MoveMapFactoryOP mmf( parse_movemap_factory( tag, data ) );
	if ( mmf ) {
		movemap_factory( mmf );
		if ( mmf_legacy ) {
			TR.Warning << "Did you intend to pass two different MoveMapFactories for this protocol?" << std::endl;
		}
	} else if ( mmf_legacy ) {
		movemap_factory( mmf_legacy );
	}
}

void
RingPlaneFlipMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( []( std::string const & name ) {
		return "RingPlaneFlipMover_subelement_" + name + "Type";
	} );
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy( xsd, subelements );

	AttributeList attlist;

	core::select::movemap::attributes_for_parse_movemap_factory_default_attr_name(
		attlist,
		"The name of a pre-defined MoveMapFactory." );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"Based on a given ResidueSelector and limited by a MoveMap, this mover selects applicable cyclic residues "
		"and performs a 180-degree shearing move in which the anomeric bond and the main-chain bond on the opposite "
		"side of the ring are moved in opposite directions.",
		attlist, subelements );
}


/// @param  <input_pose>: the structure to be moved
void
RingPlaneFlipMover::apply( Pose & input_pose )
{
	using namespace std;
	using namespace core::id;

	show( TR );

	TR << "Getting movable torsions...." << endl;

	setup_movable_torsion_pairs( input_pose );

	if ( movable_torsion_pairs_.empty() ) {
		TR.Warning << "There are no movable torsions available in the given pose." << endl;
		set_last_move_status( moves::FAIL_DO_NOT_RETRY );
		return;
	}

	TR << "Applying " << get_name() << " to pose...." << endl;

	core::uint const i( core::uint( numeric::random::rg().uniform() * movable_torsion_pairs_.size() + 1 ) );
	pair< TorsionID, TorsionID > const & movable_torsion_pair( movable_torsion_pairs_[ i ] );

	if ( TR.visible( ) ) {
		if ( movable_torsion_pair.second.rsd() ) {
			TR << "Flipping ring plane of residue " << movable_torsion_pair.second.rsd() << endl;
		} else {
			TR << "Flipping ring plane of a terminal residue off residue ";
			TR << movable_torsion_pair.second.rsd() << endl;
		}

	}

	if ( movable_torsion_pair.first.valid() ) {  // No torsion to move if this is a lower terminus.
		input_pose.set_torsion(
			movable_torsion_pair.first, input_pose.torsion( movable_torsion_pair.first ) + 180.0 );
	}
	if ( movable_torsion_pair.second.valid() ) {  // No torsion to move if this is an upper terminus.
		input_pose.set_torsion(
			movable_torsion_pair.second, input_pose.torsion( movable_torsion_pair.second ) - 180.0 );
	}

	TR << "Move(s) complete." << endl;
}


// Accessors/Mutators
core::kinematics::MoveMapCOP
RingPlaneFlipMover::movemap( core::pose::Pose const & pose ) const
{
	using namespace core::kinematics;

	if ( ! movemap_ ) {
		if ( movemap_factory_ ) {
			movemap_ = movemap_factory_->create_movemap_from_pose( pose );
		} else {
			movemap_ = MoveMapOP( new MoveMap );
		}
	}
	return movemap_;
}


// Private methods ////////////////////////////////////////////////////////////
// Initialize data members from arguments.
void
RingPlaneFlipMover::init()
{
	type( mover_name() );
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
RingPlaneFlipMover::copy_data(
	RingPlaneFlipMover & object_to_copy_to,
	RingPlaneFlipMover const & object_to_copy_from )
{
	object_to_copy_to.movemap_ = object_to_copy_from.movemap_;
	object_to_copy_to.movemap_factory_ = object_to_copy_from.movemap_factory_;
	object_to_copy_to.selector_ = object_to_copy_from.selector_;
	object_to_copy_to.movable_torsion_pairs_ = object_to_copy_from.movable_torsion_pairs_;
}

// Setup list of movable residues from MoveMap and/or ResidueSelector.
void
RingPlaneFlipMover::setup_movable_torsion_pairs( core::pose::Pose const & pose )
{
	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace id;
	using namespace select::residue_selector;
	using namespace chemical;
	using namespace rings;
	using namespace conformation;
	using namespace conformation::carbohydrates;

	movable_torsion_pairs_.clear();

	Conformation const & conf( pose.conformation() );
	if ( ! conf.contains_carbohydrate_residues() ) {
		TR.Warning << "RingPlaneFlipMover is a carbohydrate-specific Mover; ";
		TR.Warning << "this Pose contains no carbohydrate residues." << endl;
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
		if ( ! res.is_carbohydrate() ) { continue; }

		if ( res.is_branch_point() ) { continue; }  // Only linear residues can flip.

		chemical::carbohydrates::CarbohydrateInfoCOP info( res.carbohydrate_info() );
		if ( ! info->is_pyranose() ) { continue; }  // It needs to be a 6-membered ring.

		ResidueType const & type( res.type() );
		AtomIndices const & ring_atoms( type.ring_atoms( 1 ) );  // A pyranose only has 1 ring.
		vector1< PointPosition > ring_coords( 6 );
		for ( core::uint j( 1 ); j <= 6; ++j ) {
			ring_coords[ j ] = res.xyz( ring_atoms[ j ] );
		}

		TR << "Checking if " << res.name() << ' ' << i << " can flip...." << endl;

		// The 1st potential torsion to move is the pyranose residue's phi.
		vector1< id::AtomID > const & phi_atoms( get_reference_atoms_for_phi( conf, i ) );

		if ( ! phi_atoms.empty() ) {
			// If phi_atoms IS empty, res has no parent, but we can potentially still flip it with torsion 2 below.
			// If res DOES have a parent, we check if it is equatorial to the ring, by checking if the 3rd atom of the
			// definition is.
			PointPosition const glycosidic_O( pose.residue( phi_atoms[ 3 ].rsd() ).xyz( phi_atoms[ 3 ].atomno() ) );
			PointPosition const anomeric_C( res.xyz( info->anomeric_carbon_index() ) );  // same as phi_atoms[ 2 ]
			if ( is_atom_axial_or_equatorial_to_ring( glycosidic_O, anomeric_C, ring_coords ) != EQUATORIAL ) {
				TR.Debug << "Phi bond not equatorial; flipping not possible at this residue." << endl;
				continue;
			}
		}
		// This will return a BOGUS_TORSION_ID for lower termini, but that's what we want anyway.
		TorsionID const torsion1( get_non_NU_TorsionID_from_AtomIDs( conf, phi_atoms ) );

		TorsionID torsion2;
		if ( res.is_upper_terminus() ) {
			// If the residue is an upper terminus, we can safely flip it and need no further checks.
			torsion2 = TorsionID::BOGUS_TORSION_ID();
		} else /*( ! res.is_upper_terminus() )*/ {
			// The 2nd potential torsion to move is on the opposite side of the 6-membered ring,
			// so it is found off position n + 3, where n is the anomeric position;
			// chi(n+3) gives us the information we need.
			core::uint const opposite_pos( info->anomeric_carbon() + 3 );

			// We need to ensure
			// BOTH that the main chain continues from the opposite side of the ring
			// AND that the bond is equatorial.
			core::uint const child( find_seqpos_of_saccharides_child_residue_at( res, opposite_pos ) );
			if ( ! child ) {
				// If this is false, since we already know that this is not a branch point nor an upper terminus,
				// it means that the main-chain extends from a different position.  Therefore, we cannot flip.
				TR.Debug << "Main chain not opposite anomeric carbon; flipping not possible at this residue." << endl;
				continue;
			}

			// Check if it is equatorial to the ring,
			// by checking if the 3rd atom of the definition of chi(n+3) is.
			PointPosition const opposite_O( res.xyz( type.chi_atoms( opposite_pos )[ 3 ] ) );  // OK if virtual
			PointPosition const opposite_C( res.xyz( type.chi_atoms( opposite_pos )[ 2 ] ) );
			if ( is_atom_axial_or_equatorial_to_ring( opposite_O, opposite_C, ring_coords ) != EQUATORIAL ) {
				TR.Debug << "Psi(n+1) bond not equatorial; flipping not possible at this residue." << endl;
				continue;
			}

			// Even though setting the CHI angle would work, we get the BB torsion instead,
			// because the MoveMap counts CHI separately.
			torsion2 = get_non_NU_TorsionID_from_AtomIDs( conf, get_reference_atoms_for_psi( conf, child ) );
		}

		// Check that we don't waste our time on a monosaccharide.
		if ( ( torsion1 == TorsionID::BOGUS_TORSION_ID() ) && ( torsion2 == TorsionID::BOGUS_TORSION_ID() ) ) {
			TR.Warning << "If you wish to flip monosaccharide residue " << i;
			TR.Warning << ", use a RigidBodyMover instead." << endl;
			continue;
		}

		// Finally, if the MoveMap permits it, add the torsions to the list.
		kinematics::MoveMapCOP mm( movemap( pose ) );
		if ( ( torsion1 == TorsionID::BOGUS_TORSION_ID() ) || ( mm->get( torsion1 ) ) ) {
			if ( ( torsion2 == TorsionID::BOGUS_TORSION_ID() ) || ( mm->get( torsion2 ) ) ) {
				movable_torsion_pairs_.push_back( make_pair( torsion1, torsion2 ) );
			}
		}
	}
}


// Creator methods ////////////////////////////////////////////////////////////
std::string
RingPlaneFlipMoverCreator::keyname() const
{
	return RingPlaneFlipMover::mover_name();
}

// Return an up-casted owning pointer (MoverOP) to the mover.
protocols::moves::MoverOP
RingPlaneFlipMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new RingPlaneFlipMover );
}

void
RingPlaneFlipMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RingPlaneFlipMover::provide_xml_schema( xsd );
}


// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that RingPlaneFlipMover can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, RingPlaneFlipMover const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace carbohydrates
}  // namespace protocols
