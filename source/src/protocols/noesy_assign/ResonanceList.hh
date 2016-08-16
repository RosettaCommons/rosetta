// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ResonanceList.hh
/// @brief  provides a table with atomID - chemical shift mapping
/// @detail
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_ResonanceList_hh
#define INCLUDED_protocols_noesy_assign_ResonanceList_hh


// Unit Headers
#include <protocols/noesy_assign/ResonanceList.fwd.hh>

// Package Headers
#include <protocols/noesy_assign/Resonance.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/id/NamedAtomID.hh>
#include <core/types.hh>

// Utility headers
// #include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
//#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
#include <map>

namespace protocols {
namespace noesy_assign {

/*!@detail:
the ResonanceList provides a map of chemical shifts.
each atom-chemical shift tupel has a "resonanceID" as a key. (integer)

used classes:

Resonance is an atom with chemical shift information
ResonanceIDs (typedef) is a map from ID to Resonance
ResidueMap (typedef) is a map from residue number to a vector of Resonances.


*/
class ResonanceList : public utility::pointer::ReferenceCount {

public:
	typedef utility::vector1< ResonanceOP > Resonances;

private:
	typedef std::map< core::Size, Resonances > ResidueMap;
	typedef std::map< core::Size, ResonanceOP > ResonanceIDs;

	ResonanceList( ResonanceList const& a )://private copy-c'stor to avoid that FloatingResonances have out-dated pointers.
		utility::pointer::ReferenceCount(a) //make compiler happy
	{};
	ResonanceList& operator=( ResonanceList const& ) { return *this; };

public:
	/// @brief Constructor
	ResonanceList( std::string const& sequence );

	virtual ~ResonanceList();

	/// @brief read chemical shift assignments
	void read_from_stream( std::istream& );

	/// @brief write chemical shift assignments
	void write_to_stream( std::ostream& ) const;

	/// @brief write in talos format
	void write_talos_format( std::ostream&, bool backbone_only ) const;

	/// @brief retrieve a Resonance by ResonanceID --- throws EXCN_UnknonwResonance if atom not found
	Resonance const& operator[] ( core::Size key ) const;

	/// @brief retrive a Resonance by atom --- throws EXCN_UnknonwResonance if atom not found
	Resonance const& operator[] ( core::id::NamedAtomID const& ) const;

	/// @brief all resonances of a certain residue
	/// @detail --- requires that update_residue_map() has been called (which is done by read_from_stream() )
	Resonances const& resonances_at_residue( core::Size resid ) const;

	/// @brief iterators
	typedef ResonanceIDs::const_iterator const_iterator;
	//  typedef ResonanceIDs::iterator iterator;
	const_iterator begin() const { return map_.begin(); };

	const_iterator end() const { return map_.end(); };

	/// @brief have at least one resonance for residue "resi"
	bool has_residue( core::Size resi ) const;

	/// @brief retrieve aminoacid of residue "resi"
	core::chemical::AA aa_from_resid( core::Size resi ) const;

	/// @brief retrieve the protein sequence
	std::string const& sequence() const { return sequence_; }

	/// @brief number of Resonances
	core::Size size() const { return map_.size(); };

	/// @brief first ResonanceID (given by input file)
	core::Size start_key() const { return map_.begin()->first; }

	/// @brief last ResonanceID ( given by input file )
	core::Size last_key() const { return map_.rbegin()->first; }

protected:
	/// @brief retrieve a Resonance by ResonanceID  --- no error checking
	Resonance const& operator[] ( core::Size key ) { return *map_[ key ]; };

	/// @brief sort Resonances by residue number and store in by_resid_
	void update_residue_map();

	/// @brief after this call the Proton- and LabelResonances have APs to point to their respectively connected resonance
	void update_bond_connections();

private:
	/// @brief master map...
	ResonanceIDs map_;     //Resonances are ordered by resonancesID  < resID, Resonance >

	/// @brief slave map... created by update_residue_map()
	ResidueMap by_resid_;  //Resonances are ordered by residue < res_i, vector1<Resonance> >

	/// @brief sequence of the proteion
	std::string sequence_;
};

}
}

#endif
