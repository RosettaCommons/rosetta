// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/enzymes/EnzymeManager.hh
/// @brief   Declarations and simple accessor/mutator definitions for EnzymeManager.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_enzymes_EnzymeManager_HH
#define INCLUDED_core_enzymes_EnzymeManager_HH

// Unit header
#include <core/enzymes/EnzymeData.hh>
#include <core/enzymes/EnzymeManager.fwd.hh>

// Project header
#include <core/types.hh>

// Utility header
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>


namespace core {
namespace enzymes {

/// @brief  A map of enzyme family to maps of species to maps of enzyme names to enzyme data.
typedef std::map< std::string, std::map< std::string, std::map< std::string, EnzymeData > > > EnzymeDataSet;


/// @details  This class is a singleton and manages enzyme data that should only be read from the database one time as
/// needed and shared among all instances of EnzymeMovers.
class EnzymeManager : public utility::SingletonBase< EnzymeManager > {
public:  // Declare friends ///////////////////////////////////////////////////
	friend class utility::SingletonBase< EnzymeManager >;


public:  // Static constant data access ///////////////////////////////////////
	/// @brief  Return the consensus sequence of the requested enzyme.
	static std::string const & get_consensus_sequence(
			std::string const & family,
			std::string const & species,
			std::string const & enzyme );

	/// @brief  Return the consensus sequence type of the requested enzyme.
	static ConsensusSequenceType get_consensus_sequence_type(
			std::string const & family,
			std::string const & species,
			std::string const & enzyme );

	/// @brief  Return the efficiency of the requested enzyme.
	static core::Real get_efficiency(
			std::string const & family,
			std::string const & species,
			std::string const & enzyme );

	/// @brief  Return the second substrates or byproducts of the requested enzyme.
	static utility::vector1< std::string > const & get_second_substrates_or_byproducts(
			std::string const & family,
			std::string const & species,
			std::string const & enzyme );


	/// @brief  Return the identifiers (such as 3-letter codes) of the residues of the consensus sequence.
	static utility::vector1< utility::vector1< std::string > > const & get_consensus_residues(
			std::string const & family,
			std::string const & species,
			std::string const & enzyme );

	/// @brief  Return the position in the consensus sequence of the reactive residue.
	static core::uint get_reactive_residue_consensus_sequence_position(
			std::string const & family,
			std::string const & species,
			std::string const & enzyme );

	/// @brief  Return the name of the reactive site atom.
	static std::string const & get_reactive_atom(
			std::string const & family,
			std::string const & species,
			std::string const & enzyme );


private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	EnzymeManager();


	// Parse the consensus sequence within this EnzymeData and derive a list of consensus residues.
	void parse_consensus_sequence( EnzymeData & enzyme_data ) const;


private:  // Private data /////////////////////////////////////////////////////
	EnzymeDataSet enzymes_;
};  // class EnzymeManager

}  // namespace enzymes
}  // namespace core

#endif  // INCLUDED_core_enzymes_EnzymeManager_HH
