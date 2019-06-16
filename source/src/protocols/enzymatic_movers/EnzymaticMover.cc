// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzymatic_movers/EnzymaticMover.fwd.hh
/// @brief  Method definitions for the base class EnzymaticMover
/// @author Labonte <JWLabonte@jhu.edu>


// Unit header
#include <protocols/enzymatic_movers/EnzymaticMover.hh>

// Package header
#include <core/enzymes/EnzymeManager.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <protocols/moves/mover_schemas.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Basic header
#include <basic/options/option.hh>
#include <basic/options/keys/enzymes.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


// Construct tracers.
static basic::Tracer TR( "protocols.enzymatic_movers.EnzymaticMover" );


namespace protocols {
namespace enzymatic_movers {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
EnzymaticMover::EnzymaticMover(): moves::Mover()
{
	init( "" );
}

// Constructor with enzyme family provided
EnzymaticMover::EnzymaticMover( std::string const & enzyme_family )
{
	init( enzyme_family );
}

// Copy constructor
EnzymaticMover::EnzymaticMover( EnzymaticMover const & object_to_copy ) : moves::Mover( object_to_copy )
{
	copy_data( *this, object_to_copy );
}

// Assignment operator
EnzymaticMover &
EnzymaticMover::operator=( EnzymaticMover const & object_to_copy )
{
	// Abort self-assignment.
	if ( this != &object_to_copy ) {
		moves::Mover::operator=( object_to_copy );
		copy_data( *this, object_to_copy );
	}
	return *this;
}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
EnzymaticMover::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::enzymes::species );
	option.add_relevant( OptionKeys::enzymes::enzyme );
	option.add_relevant( OptionKeys::enzymes::efficiency );

	// Mover's register_options() doesn't do anything; it's just here in principle.
	moves::Mover::register_options();
}

void
EnzymaticMover::show( std::ostream & output ) const
{
	using namespace std;

	moves::Mover::show( output );  // name, type, tag

	output << "Simulated Enzyme: " << enzyme_name_ << " from " << species_name_ << endl;

	output << "Reaction Efficiency: " << efficiency_ * 100 << '%' << endl;
	output << "Residues Excluded from Potential Reaction: ";
	core::Size const n_excluded_sites( excluded_sites_.size() );
	for ( core::uint i( 1 ); i <= n_excluded_sites; ++i ) {
		output << excluded_sites_[ i ] << ' ';
	}
	output << endl;
	output << "Residues Ensured Modification if Applicable: ";
	core::Size const n_ensured_sites( ensured_sites_.size() );
	for ( core::uint i( 1 ); i <= n_ensured_sites; ++i ) {
		output << ensured_sites_[ i ] << ' ';
	}
	output << endl;
	if ( ( get_n_co_substrates() ) && ( performs_major_reaction_only_ ) ) {
		output << "Reaction limited to " << co_substrates_[ 1 ] << endl;
	}
}


void
EnzymaticMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	Pose const & /*pose*/ )
{
	set_species( tag->getOption< std::string >( "species", "h_sapiens" ) );
	set_enzyme( tag->getOption< std::string >( "enzyme_name", "generic" ) );
	set_efficiency( tag->getOption< core::Real >( "efficiency", 1.00 ) );
	if ( tag->getOption< bool >( "perform_major_reaction_only", false ) ) { perform_major_reaction_only(); }
	if ( tag->getOption< bool >( "perform_all_reactions", true ) ) { perform_all_reactions(); }
}

utility::tag::XMLSchemaComplexTypeGeneratorOP
EnzymaticMover::xml_schema_complex_type_generator()
{
	using namespace utility::tag;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "species", xs_string, "Set the species name of this simulated enzyme." )
		+ XMLSchemaAttribute( "enzyme_name", xs_string, "Set the specific name of this simulated enzyme." )
		+ XMLSchemaAttribute( "efficiency", xsct_real, "Directly set the efficiency of this enzyme, ignoring whatever is in the database." )
		+ XMLSchemaAttribute( "perform_major_reaction_only", xsct_rosetta_bool, "Set this EnzymaticMover to perform only its major reaction." )
		+ XMLSchemaAttribute( "perform_all_reactions", xsct_rosetta_bool, "Allow this EnzymaticMover to be promiscuous, performing a random transfer from among its possible co-substrates." );

	XMLSchemaComplexTypeGeneratorOP ct_gen(
		utility::pointer::make_shared< XMLSchemaComplexTypeGenerator >() );
	ct_gen->
		add_attributes( attlist )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_optional_name_attribute();

	return ct_gen;
}


/// @details  When applied, every EnzymaticMover first obtains a list of all possible reaction sites, using enzyme
///           data.  Next, it loops through the sites and checks if they are ensured.  If not, it decides whether or
///           not to modify the site based on the enzymes efficiency, as recorded in the database.  Then it performs
///           the specific "reaction" of the particular EnzymaticMover at that site.
/// @param    <input_pose>: the structure to be post-translationally modified, i.e., "substrate 1"
void
EnzymaticMover::apply( Pose & input_pose )
{
	using namespace std;

	show( TR );

	set_pose_reactive_sites( input_pose );

	Size const n_sites( get_n_reactive_sites() );

	if ( ! n_sites ) {
		set_last_move_status( moves::FAIL_DO_NOT_RETRY );
		return;
	}

	TR << "Simulating " << get_enzyme() << " enzyme acting on the pose...." << endl;

	for ( core::uint i( 1 ); i <= n_sites; ++i ) {
		if ( ! get_ensured_sites().contains( get_reactive_site_sequence_position( i ) ) ) {
			// If the site is not ensured, randomly decide whether or not to react.
			if ( numeric::random::rg().uniform() > get_efficiency() ) {
				TR << " Modification site " << i << " skipped because of low simulated enzyme efficiency." << endl;
				continue;
			}
		}
		if ( get_n_co_substrates() ) {
			core::uint j;
			if ( performs_major_reaction_only() ) {
				j = 1;
			} else {
				j = core::uint( numeric::random::rg().uniform() * get_n_co_substrates() + 1 );
			}
			//ConsensusSequenceType const sequence_type(
			//  EnzymeManager::get_consensus_sequence_type( species_name_, enzyme_name_ ) );
			perform_reaction( input_pose, i, get_co_substrate( j ) );
		} else /*This reaction does not require a cosubstrate from Rosetta's point of view.*/ {
			perform_reaction( input_pose, i );
		}
	}

	TR << "Move(s) complete." << endl;
}


// Accessors/Mutators
// Set the family name of this simulated enzyme.
/// @details  This method is protected; it should only ever be called by a derived class.
void
EnzymaticMover::set_enzyme_family( std::string const & family_name )
{
	enzyme_family_ = family_name;
}


// Set the species name of this simulated enzyme.
/// @details  Setting the species name limits the behavior of this EnzymaticMover to reactions known to occur in the
/// given species.
/// @param    <setting>: A species name in the format "e_coli" or "h_sapiens", which must correspond to a directory
/// in the Rosetta database, e.g., "database/virtual_enzymes/glycosyltransferases/h_sapiens/"
void
EnzymaticMover::set_species( std::string const & species_name )
{
	species_name_ = species_name;
	if ( ! enzyme_name_.empty() ) {
		set_efficiency();
		set_available_co_substrates();
	}
}


// Set the specific enzyme name of this simulated enzyme.
/// @details  If set, this EnzymaticMover will use specific enzymatic details for this reaction from the database.
/// If the species name has not been set, an enzyme from "h_sapiens" is assumed.
/// @param    <setting>: An enzyme name as listed in an appropriate enzyme file in
/// "database/virtual_enzymes/glycosyltransferases/<species_name_>/"
void
EnzymaticMover::set_enzyme( std::string const & enzyme_name )
{
	enzyme_name_ = enzyme_name;
	if ( ! species_name_.empty() ) {
		set_efficiency();
		set_available_co_substrates();
	}
}


// Do not modify this site, even if it is within a consensus sequence match for this enzyme.
void
EnzymaticMover::exclude_site( core::uint seqpos )
{
	excluded_sites_.push_back( seqpos );
}

// Definitely modify this site, if it is within a consensus sequence match for this enzyme.
void
EnzymaticMover::ensure_site( core::uint seqpos )
{
	ensured_sites_.push_back( seqpos );
}


// Other Methods //////////////////////////////////////////////////////////////
// Access the EnzymeManager to determine which sites on a given Pose are able to be modified.
void
EnzymaticMover::set_pose_reactive_sites( core::pose::Pose const & pose )
{
	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace core::enzymes;

	reaction_sites_.clear();

	string const & consensus( EnzymeManager::get_consensus_sequence( enzyme_family_, species_name_, enzyme_name_ ) );
	vector1< vector1< string > > const & consensus_residues(
		EnzymeManager::get_consensus_residues( enzyme_family_, species_name_, enzyme_name_ ) );
	core::uint const site_residue_position( EnzymeManager::get_reactive_residue_consensus_sequence_position(
		enzyme_family_, species_name_, enzyme_name_ ) );
	vector1< string > const & site_residues( consensus_residues[ site_residue_position ] );  // could be more than one
	Size const n_consensus_residues( consensus_residues.size() );
	Size const n_residues_left_of_site( site_residue_position - 1 );
	Size const n_residues_right_of_site( n_consensus_residues - site_residue_position );
	string const & site_atom( EnzymeManager::get_reactive_atom( enzyme_family_, species_name_, enzyme_name_ ) );

	TR << "Searching for reactive sites within consensus sequence " << consensus << "..." << endl;

	Size const n_residues( pose.total_residue() );
	for ( core::uint i( 1 ); i <= n_residues; ++i ) {
		if ( excluded_sites_.contains( i ) ) { continue; }
		string const & res_i( pose.residue( i ).name3() );
		if ( site_residues.contains( res_i ) ) {
			TR.Debug << "Found potential site: " << res_i << i;
			TR.Debug << "  Checking if within consensus sequence..." << endl;
			if ( i < site_residue_position ) {
				TR.Trace << "Residue too close to start of sequence to fall within consensus." << endl;
				continue;
			}
			if ( i + n_residues_right_of_site > n_residues ) {
				TR.Trace << "Residue too close to end of sequence to fall within consensus." << endl;
				continue;
			}

			// Check left.
			bool left_matches( true );
			for ( core::uint j( 1 ); j <= n_residues_left_of_site; ++j ) {
				if ( ! consensus_residues[ site_residue_position - j ].contains( pose.residue( i - j ).name3() ) ) {
					TR.Trace << "Residue " << pose.residue( i - j ).name3() << i - j;
					TR.Trace << " found left of " << res_i << i;
					TR.Trace << " does not fall within consensus." << endl;
					left_matches = false;
					break;
				}
			}
			if ( ! left_matches ) { continue; }

			// Check right.
			bool right_matches( true );
			for ( core::uint j( 1 ); j <= n_residues_right_of_site; ++j ) {
				if ( ! consensus_residues[ site_residue_position + j ].contains( pose.residue( i + j ).name3() ) ) {
					TR.Trace << "Residue " << pose.residue( i + j ).name3() << i + j;
					TR.Trace << " found right of " << res_i << i;
					TR.Trace << " does not fall within consensus." << endl;
					right_matches = false;
					break;
				}
			}
			if ( ! right_matches ) { continue; }

			TR.Trace << res_i << i << " is a non-excluded match; adding..." << endl;
			reaction_sites_.push_back( make_pair( i, site_atom ) );
		}
	}

	TR << "Found " << reaction_sites_.size() << " potential reaction sites." << endl;
}


// Private methods ////////////////////////////////////////////////////////////
// Set command-line options.  (Called by init())
void
EnzymaticMover::set_commandline_options()
{
	using namespace basic::options;

	set_species( option[ OptionKeys::enzymes::species ] );
	set_enzyme( option[ OptionKeys::enzymes::enzyme ] );
	if ( option[ OptionKeys::enzymes::efficiency ].active() ) {
		set_efficiency( option[ OptionKeys::enzymes::efficiency ] );
	}
}

// Initialize data members from arguments.
void
EnzymaticMover::init( std::string const & enzyme_family )
{
	TR << "Initializing EnzymaticMover from the " << enzyme_family << " family..." << std::endl;
	set_enzyme_family( enzyme_family );
	set_commandline_options();
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
EnzymaticMover::copy_data( EnzymaticMover & object_to_copy_to, EnzymaticMover const & object_to_copy_from )
{
	object_to_copy_to.enzyme_family_ = object_to_copy_from.enzyme_family_;
	object_to_copy_to.species_name_ = object_to_copy_from.species_name_;
	object_to_copy_to.enzyme_name_ = object_to_copy_from.enzyme_name_;
	object_to_copy_to.efficiency_ = object_to_copy_from.efficiency_;
	object_to_copy_to.excluded_sites_ = object_to_copy_from.excluded_sites_;
	object_to_copy_to.ensured_sites_ = object_to_copy_from.ensured_sites_;
	object_to_copy_to.reaction_sites_ = object_to_copy_from.reaction_sites_;
	object_to_copy_to.co_substrates_ = object_to_copy_from.co_substrates_;
	object_to_copy_to.performs_major_reaction_only_ = object_to_copy_from.performs_major_reaction_only_;
}


// Access the EnzymeManager to get efficiency.
void
EnzymaticMover::set_efficiency()
{
	efficiency_ = core::enzymes::EnzymeManager::get_efficiency( enzyme_family_, species_name_, enzyme_name_ );
}

// Access the EnzymeManager to determine which co-substrates are able to take part in the reaction of this particular
// enzyme.
void
EnzymaticMover::set_available_co_substrates()
{
	co_substrates_ = core::enzymes::EnzymeManager::get_second_substrates_or_byproducts(
		enzyme_family_, species_name_, enzyme_name_ );
}

}  // namespace enzymatic_movers
}  // namespace protocols
