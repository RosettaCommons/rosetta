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
#include <protocols/rosetta_scripts/util.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Basic header
#include <basic/options/option.hh>
#include <basic/Tracer.hh>


// Construct tracers.
static THREAD_LOCAL basic::Tracer TR( "protocols.enzymatic_movers.EnzymaticMover" );


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

// Destructor
EnzymaticMover::~EnzymaticMover() {}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
EnzymaticMover::register_options()
{
	using namespace basic::options;

	//option.add_relevant( OptionKeys::foo::bar );

	moves::Mover::register_options();  // Mover's register_options() doesn't do anything; it's just here in principle.
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
	if ( performs_major_reaction_only_ ) {
		output << "Reaction limited to " << co_substrates_[ 1 ] << endl;
	}
}


void
EnzymaticMover::parse_my_tag(
	TagCOP /*tag*/,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	Pose const & /*pose*/ )
{}


/// @details  WiP
/// @param    <input_pose>: the structure to be glycosylated, i.e., "substrate 1"
void
EnzymaticMover::apply( Pose & input_pose )
{
	using namespace std;

	show( TR );

	set_pose_reactive_sites( input_pose );

	TR << "Simulating " << get_enzyme() << " acting on the pose...." << endl;

	Size const n_sites( get_n_reactive_sites() );
	for ( core::uint i( 1 ); i <= n_sites; ++i ) {
		if ( ! get_ensured_sites().contains( get_reactive_site_sequence_position( i ) ) ) {
			// If the site is not ensured, randomly decide whether or not to react.
			if ( numeric::random::rg().uniform() > get_efficiency() ) {
				continue;
			}
		}
		core::uint j;
		if ( performs_major_reaction_only() ) {
			j = 1;
		} else {
			j = core::uint( numeric::random::rg().uniform() * get_n_co_substrates() + 1 );
		}
		//ConsensusSequenceType const sequence_type(
		//  EnzymeManager::get_consensus_sequence_type( species_name_, enzyme_name_ ) );
		perform_reaction( input_pose, i, get_co_substrate( j ) );
	}

	TR << "Move(s) complete." << endl;
}


// Accessors/Mutators
// Set the family name of this simulated enzyme.
/// @details  This method is protected; it should only ever be called by a derived class.
void
EnzymaticMover::set_enzyme_family( std::string const & family_name )
{
	// TODO: Check that this family exists.  If not, throw an exit with message.
	enzyme_family_ = family_name;
}


// Set the species name of this simulated glycosyltransferase.
/// @details  Setting the species name limits the behavior of this EnzymaticMover to reactions known to occur in the
/// given species.
/// @param    <setting>: A species name in the format "e_coli" or "h_sapiens", which must correspond to a directory
/// in the Rosetta database, e.g., "database/virtual_enzymes/glycosyltransferases/h_sapiens/"
void
EnzymaticMover::set_species( std::string const & species_name )
{
	// TODO: Check that this species has the currently set enzyme.  If not, throw an exit with message.
	species_name_ = species_name;
	if ( ! enzyme_name_.empty() ) {
		set_efficiency();
		set_available_co_substrates();
	}
}


// Set the specific enzyme name of this simulated glycosyltransferase.
/// @details  If set, this EnzymaticMover will use specific enzymatic details for this reaction from the database.
/// If the species name has not been set, an enzyme from "h_sapiens" is assumed.
/// @param    <setting>: An enzyme name as listed in an appropriate enzyme file in
/// "database/virtual_enzymes/glycosyltransferases/<species_name_>/"
void
EnzymaticMover::set_enzyme( std::string const & enzyme_name )
{
	// TODO: Check that this enzyme is found in this species.  If not, throw an exit with message.
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
// Access the EnzymeManager to determine which sites on a given Pose are able to be glycosylated.
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
	Size const n_consensus_residues( site_residues.size() );
	Size const n_residues_left_of_site( site_residue_position - 1 );
	Size const n_residues_right_of_site( n_consensus_residues - site_residue_position );
	string const & site_atom( EnzymeManager::get_reactive_atom( enzyme_family_, species_name_, enzyme_name_ ) );

	TR << "Searching for reactive sites within consensus sequence " << consensus << "..." << endl;

	Size const n_residues( pose.total_residue() );
	for ( core::uint i( 1 ); i <= n_residues; ++i ) {
		if ( excluded_sites_.contains( i ) ) { continue; }
		if ( site_residues.contains( pose.residue( i ).name3() ) ) {
			TR.Debug << "Found potential site: " << pose.residue( i ).name3() << i;
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
			for ( core::uint j( 1 ); j <= n_residues_left_of_site; ++j ) {
				if ( ! consensus_residues[ site_residue_position - j ].contains( pose.residue( i - j ).name3() ) ) {
					TR.Trace << "Residue " << pose.residue( i - j ).name3() << i - j;
					TR.Trace << " found left of " << pose.residue( i ).name3() << i;
					TR.Trace << " does not fall within consensus." << endl;
					continue;
				}
			}

			// Check right.
			for ( core::uint j( 1 ); j <= n_residues_right_of_site; ++j ) {
				if ( ! consensus_residues[ site_residue_position + j ].contains( pose.residue( i + j ).name3() ) ) {
					TR.Trace << "Residue " << pose.residue( i + j ).name3() << i + j;
					TR.Trace << " found right of " << pose.residue( i ).name3() << i;
					TR.Trace << " does not fall within consensus." << endl;
					continue;
				}
			}

			TR.Trace << pose.residue( i ).name3() << i << " is a non-excluded match; adding..." << endl;
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
}

// Initialize data members from arguments.
void
EnzymaticMover::init( std::string const & enzyme_family )
{
	type( "EnzymaticMover" );
	set_enzyme_family( enzyme_family );
	set_commandline_options();

	// Set defaults.
	species_name_ = "";
	enzyme_name_ = "";
	performs_major_reaction_only_ = false;  // Allows for promiscuous enzymes by default.
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
