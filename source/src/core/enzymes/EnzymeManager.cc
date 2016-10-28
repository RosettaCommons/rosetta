// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/enzymes/EnzymeManager.hh
/// @brief   Method definitions for EnzymeManager.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <core/enzymes/EnzymeData.hh>
#include <core/enzymes/EnzymeManager.hh>
#include <core/enzymes/consensus_sequence_parsers.hh>
#include <core/enzymes/database_io.hh>

// Basic header
#include <basic/database/open.hh>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif
#endif


namespace core {
namespace enzymes {

// Public methods /////////////////////////////////////////////////////////////
// Static constant data access
// Return the consensus sequence of the requested enzyme.
std::string const &
EnzymeManager::get_consensus_sequence(
		std::string const & family,
		std::string const & species,
		std::string const & enzyme )
{
	return get_instance()->
			enzymes_.find( family )->second.find( species )->second.find( enzyme )->second.consensus_sequence;
}

// Return the consensus sequence type of the requested enzyme.
/// @return  a ConsensusSequenceType enum value
ConsensusSequenceType
EnzymeManager::get_consensus_sequence_type(
		std::string const & family,
		std::string const & species,
		std::string const & enzyme )
{
	return get_instance()->enzymes_.find( family )->second.find( species )->second.find( enzyme )->second.cs_type;
}

// Return the efficiency of the requested enzyme.
core::Real
EnzymeManager::get_efficiency(
		std::string const & family,
		std::string const & species,
		std::string const & enzyme )
{
	return get_instance()->enzymes_.find( family )->second.find( species )->second.find( enzyme )->second.efficiency;
}

// Return the second substrates or byproducts of the requested enzyme.
/// @details  Depending on whether this is a transferase or hydrolase, this will return a list of possible substrates or
/// possible byproducts, respectively.
utility::vector1< std::string > const &
EnzymeManager::get_second_substrates_or_byproducts(
		std::string const & family,
		std::string const & species,
		std::string const & enzyme )
{
	return get_instance()->enzymes_.find( family )->second.find( species )->second.find( enzyme )->
			second.second_substrates_or_byproducts;
}


// Return the identifiers (such as 3-letter codes) of the residues of the consensus sequence.
/// @return  a vector of vectors of strings, because each position may have more than one possible match.
utility::vector1< utility::vector1< std::string > > const &
EnzymeManager::get_consensus_residues(
		std::string const & family,
		std::string const & species,
		std::string const & enzyme )
{
	return get_instance()->
			enzymes_.find( family )->second.find( species )->second.find( enzyme )->second.consensus_residues;
}

// Return the position in the consensus sequence of the reactive residue.
core::uint
EnzymeManager::get_reactive_residue_consensus_sequence_position(
		std::string const & family,
		std::string const & species,
		std::string const & enzyme )
{
	return get_instance()->enzymes_.find( family )->second.find( species )->second.find( enzyme )->second.cs_resnum;
}

// Return the name of the reactive site atom.
std::string const &
EnzymeManager::get_reactive_atom(
		std::string const & family,
		std::string const & species,
		std::string const & enzyme )
{
	return get_instance()->
			enzymes_.find( family )->second.find( species )->second.find( enzyme )->second.atom_to_modify;
}


// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
EnzymeManager::EnzymeManager()
{
	// TODO: Refactor to lazily load based on the requested family-species-name combos in a thread-safe manner.
	EnzymeData enzyme_data( read_enzyme_data_from_file( basic::database::full_name(
			"virtual_enzymes/glycosyltransferases/h_sapiens/generic_N-glycosyltransferase" ) ) );  // TEMP

	parse_consensus_sequence( enzyme_data );

	enzymes_[ "glycosyltransferases" ][ "h_sapiens" ][ "generic_N-glycosyltransferase" ] = enzyme_data;  // TEMP
}


// Parse the consensus sequence within this EnzymeData and derive a list of consensus residues.
void
EnzymeManager::parse_consensus_sequence( EnzymeData & enzyme_data ) const
{
	std::string const & consensus( enzyme_data.consensus_sequence );

	switch ( enzyme_data.cs_type ) {
		case AA:
			enzyme_data.consensus_residues = get_3_letter_codes_from_peptide_consensus_sequence( consensus );
			break;
		case NA:
			enzyme_data.consensus_residues = get_codes_from_NA_consensus_sequence( consensus );
			break;
		case SACCHARIDE:
			enzyme_data.consensus_residues = get_codes_from_saccharide_consensus_sequence( consensus );
			break;
	}
}

}  // namespace enzymes
}  // namespace core
