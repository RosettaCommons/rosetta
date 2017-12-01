// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/simple_moves/RingConformationMover.cc
/// @brief   Method definitions for RingConformationMover.
/// @author  Labonte  <JWLabonte@jhu.edu>

// Unit headers
#include <protocols/simple_moves/RingConformationMover.hh>
#include <protocols/simple_moves/RingConformationMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/rings/RingConformer.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/conformation/Residue.hh>

#include <protocols/rosetta_scripts/util.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Basic header
#include <basic/options/option.hh>
#include <basic/options/keys/rings.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <iostream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


// Construct tracers.
static basic::Tracer TR( "protocols.simple_moves.RingConformationMover" );


namespace protocols {
namespace simple_moves {

using namespace core;


// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
/// @details  By default, all rings within a given pose will be allowed to move.
RingConformationMover::RingConformationMover(): Mover()
{
	using namespace kinematics;

	init();
}

// Copy constructor
RingConformationMover::RingConformationMover( RingConformationMover const & object_to_copy ) : Mover( object_to_copy )
{
	copy_data( *this, object_to_copy );
}

// Constructor with MoveMap input option
/// @param    <input_movemap>: a MoveMap with desired nu torsions set to true
/// @remarks  Movable cyclic residues will generally be a subset of residues in the MoveMap whose nu
/// torsions are set to true.
RingConformationMover::RingConformationMover( core::kinematics::MoveMapOP input_movemap ):
	movemap_( input_movemap )
{
	init();
}

// Assignment operator
RingConformationMover &
RingConformationMover::operator=( RingConformationMover const & object_to_copy )
{
	// Abort self-assignment.
	if ( this == &object_to_copy ) {
		return *this;
	}

	Mover::operator=( object_to_copy );
	copy_data( *this, object_to_copy );
	return *this;
}

// Destructor
RingConformationMover::~RingConformationMover() = default;


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
RingConformationMover::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::rings::lock_rings );
	option.add_relevant( OptionKeys::rings::sample_high_energy_conformers );

	// Call register_options() on all other Movers used by this class.
	Mover::register_options();  // Mover's register_options() doesn't do anything; it's just here in principle.
}

void
RingConformationMover::show(std::ostream & output) const
{
	using namespace std;

	Mover::show( output );  // name, type, tag

	if ( locked_ ) {
		output << "This Mover was locked from the command line with the -lock_rings flag " <<
			"and will not do anything!" << endl;
	} else {
		if ( sample_all_conformers_ ) {
			output << "Sampling from all ideal ring conformations." << endl;
		} else {
			output << "Sampling from only low-energy ring conformations, if known." << endl;
		}
	}
}


// Mover methods
protocols::moves::MoverOP
RingConformationMover::clone() const
{
	return protocols::moves::MoverOP( new RingConformationMover( *this ) );
}

protocols::moves::MoverOP
RingConformationMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new RingConformationMover() );
}

void
RingConformationMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	Pose const & )
{
	using namespace core::kinematics;

	// Parse the MoveMap tag.
	movemap_factory( protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data ) );

	// Parse option specific to rings.
	if ( tag->hasOption( "sample_high_energy_conformers" ) ) {
		sample_all_conformers_ = tag->getOption< bool >( "sample_high_energy_conformers" );
	}
}


/// @details  The mover will create a list of movable residues based on the given MoveMap and select a residue from
/// the list at random.  The torsion angles of a randomly selected ring conformer will be applied to the selected
/// residue.
/// @param    <input_pose>: the structure to be moved
void
RingConformationMover::apply( Pose & input_pose )
{
	using namespace std;
	using namespace utility;
	using namespace conformation;
	using namespace chemical;

	if ( locked_ ) {
		TR.Warning << "Rings were locked from the command line with the -lock_rings flag.  " <<
			"Did you intend to call this Mover?" << endl;
		return;
	}

	show( TR );

	TR << "Getting movable residues...." << endl;

	setup_residue_list( input_pose );

	if ( residue_list_.empty() ) {
		TR.Warning << "There are no movable cyclic residues available in the given pose." << endl;
		return;
	}

	TR << "Applying " << get_name() << " to pose...." << endl;

	core::uint const i( core::uint( numeric::random::rg().uniform() * residue_list_.size() + 1 ) );
	core::uint const res_num( residue_list_[ i ] );
	Residue const & res( input_pose.residue( res_num ) );

	TR << "Selected residue " << res_num << ": " << res.name() << endl;

	Size const n_rings( res.type().n_rings() );
	for ( core::uint ring_num( 1 ); ring_num <= n_rings; ++ring_num ) {
		rings::RingConformer conformer;
		if ( sample_all_conformers_ ||
				( ! res.type().ring_conformer_set( ring_num )->low_energy_conformers_are_known() ) ) {
			conformer = res.type().ring_conformer_set( ring_num )->get_random_conformer();
		} else /* default behavior */ {
			conformer = res.type().ring_conformer_set( ring_num )->get_random_local_min_conformer();
		}

		TR << "Selected the " << conformer.specific_name <<
			" conformation to apply to ring " << ring_num << '.' << endl;

		TR << "Making move...." << endl;

		input_pose.set_ring_conformation( res_num, ring_num, conformer );
	}

	TR << "Move(s) complete." << endl;
}


// Accessors/Mutators
kinematics::MoveMapCOP
RingConformationMover::movemap( core::pose::Pose const & pose ) const
{
	if ( movemap_ ) {
		return movemap_;
	} else if ( movemap_factory_ ) {
		return movemap_factory_->create_movemap_from_pose( pose );
	} else {
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		movemap->set_nu( true );
		return movemap;
	}
}

void
RingConformationMover::movemap( kinematics::MoveMapOP new_movemap )
{
	movemap_ = new_movemap;
}

void
RingConformationMover::movemap_factory( select::movemap::MoveMapFactoryCOP new_movemap_factory )
{
	movemap_factory_ = new_movemap_factory;
}


// Private methods ////////////////////////////////////////////////////////////
// Set command-line options.  (Called by init())
void
RingConformationMover::set_commandline_options()
{
	using namespace basic::options;

	locked_ = option[ OptionKeys::rings::lock_rings ];

	sample_all_conformers_ = option[ OptionKeys::rings::sample_high_energy_conformers ];
}

// Initialize data members from arguments.
void
RingConformationMover::init()
{
	type( "RingConformationMover" );

	locked_ = false;
	sample_all_conformers_ = false;

	set_commandline_options();
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
RingConformationMover::copy_data(
	RingConformationMover & object_to_copy_to,
	RingConformationMover const & object_to_copy_from )
{
	object_to_copy_to.movemap_factory_ = object_to_copy_from.movemap_factory_;
	object_to_copy_to.movemap_ = object_to_copy_from.movemap_;
	object_to_copy_to.residue_list_ = object_to_copy_from.residue_list_;
	object_to_copy_to.locked_ = object_to_copy_from.locked_;
	object_to_copy_to.sample_all_conformers_ = object_to_copy_from.sample_all_conformers_;
}

// Setup list of movable cyclic residues from MoveMap.
void
RingConformationMover::setup_residue_list( core::pose::Pose & pose )
{
	using namespace conformation;

	residue_list_.clear();

	core::kinematics::MoveMapCOP my_movemap( movemap(pose) );
	Size const last_res_num( pose.size() );
	for ( core::uint res_num( 1 ); res_num <= last_res_num; ++res_num ) {
		Residue const & residue( pose.residue( res_num ) );
		if ( residue.type().is_cyclic() ) {
			if ( my_movemap->get_nu( res_num ) == true ) {
				residue_list_.push_back( res_num );
			}
		}
	}
}


// Creator methods ////////////////////////////////////////////////////////////
// Return an up-casted owning pointer (MoverOP) to the mover.
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RingConformationMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return moves::MoverOP( new RingConformationMover );
// XRW TEMP }

// Return the string identifier for the associated Mover (RingConformationMover).
std::string RingConformationMover::get_name() const {
	return mover_name();
}

std::string RingConformationMover::mover_name() {
	return "RingConformationMover";
}

void RingConformationMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( [] (std::string const& name) {
		return "RingConformationMover_subelement_" + name + "Type";
	});
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy(xsd, subelements);

	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"sample_high_energy_conformers", xsct_rosetta_bool,
		"XRW TO DO");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"Based on a given MoveMap, this mover selects movable cyclic residues and "
		"flips their rings to an idealized ring conformer",
		attlist, subelements );
}

std::string RingConformationMoverCreator::keyname() const {
	return RingConformationMover::mover_name();
}

protocols::moves::MoverOP
RingConformationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RingConformationMover );
}

void RingConformationMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RingConformationMover::provide_xml_schema( xsd );
}



// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that RingConformationMover can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, RingConformationMover const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace simple_moves
}  // namespace protocols
