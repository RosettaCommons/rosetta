// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_AA_hh
#define INCLUDED_core_chemical_AA_hh


// Unit headers

// Project headers

// Utility headers

// C++ headers
// Commented by inclean daemon #include <string>

// AUTO-REMOVED #include <basic/Tracer.fwd.hh>

#include <ostream>



namespace core {
namespace chemical {

// temporary -- probably an enum?
///////////////////////////////////////////////////////////////////////////
/// @brief enumeration for amino acids and nucleotides types with the total
/// number as num_aa_types
///////////////////////////////////////////////////////////////////////////
// BUT DONT CODE TO THESE AS INTS!!!!!!
enum AA {
	// protein 1-20
	aa_ala = 1,
	aa_cys,
	aa_asp,
	aa_glu,
	aa_phe,
	aa_gly,
	aa_his,
	aa_ile,
	aa_lys,
	aa_leu,
	aa_met,
	aa_asn,
	aa_pro,
	aa_gln,
	aa_arg,
	aa_ser,
	aa_thr,
	aa_val,
	aa_trp,
	aa_tyr,
	num_canonical_aas = aa_tyr,

	// dna 21-24
	na_ade,
	first_DNA_aa = na_ade,
	na_cyt,
	na_gua,
	na_thy,
	last_DNA_aa = na_thy,

	// rna 25-28
	na_rgu,
	na_rad,
	na_rcy,
	na_ura,

	//D-amino acids 29-47.  Keep these together as a group, and ensure that aa_dal is first and aa_dty is last.
	aa_dal,
	first_D_aa = aa_dal,
	aa_dcs,
	aa_das,
	aa_dgu,
	aa_dph,
	aa_dhi,
	aa_dil,
	aa_dly,
	aa_dle,
	aa_dme,
	aa_dan,
	aa_dpr,
	aa_dgn,
	aa_dar,
	aa_dse,
	aa_dth,
	aa_dva,
	aa_dtr,
	aa_dty,
	last_D_aa = aa_dty,

	//Beta-3-amino acids 48-71.  Keep these together as a group, and ensure that last_beta3_aa is set to whatever is last.
	aa_b3a,
	first_beta3_aa = aa_b3a,
	aa_b3c,
	aa_b3d,
	aa_b3e,
	aa_b3f,
	aa_b3g,
	aa_b3h,
	aa_b3i,
	aa_b3k,
	aa_b3l,
	aa_b3m,
	aa_b3n,
	aa_b3p,
	aa_b3q,
	aa_b3r,
	aa_b3s,
	aa_b3t,
	aa_b3v,
	aa_b3w,
	aa_b3y,
	aa_b3cisACPrC,
	aa_b3cisACPC,
	aa_b3cisACHC,
	last_beta3_aa = aa_b3cisACHC,

	// h2o
	aa_h2o,

	// virtual
	aa_vrt,

  // membrane
  aa_mpr,

	// unknown polymer
	aa_unp,

	// unknown
	aa_unk,

	num_aa_types = aa_unk //keep this guy last
};

//////////////////////////////////////////////////////////
/// @brief Give an AA string name, return its enum type.
//////////////////////////////////////////////////////////
AA
aa_from_name( std::string const & name );

//////////////////////////////////////////////////////////
/// @brief Give an enum type, return true if and only if
/// it is a D-amino acid.
//////////////////////////////////////////////////////////
bool
is_D_aa( AA aa );

//////////////////////////////////////////////////////////
/// @brief Given an enum type for a D-amino acid with a
/// canonical side-chain, return the enum type for the
/// corresponding L-amino acid (or aa_unk if the
/// corresponding L-amino acid cannot be determined).
//////////////////////////////////////////////////////////
AA
get_L_equivalent( AA aa );

//////////////////////////////////////////////////////////
/// @brief Given an enum type for a L-amino acid with a
/// canonical side-chain, return the enum type for the
/// corresponding D-amino acid (or aa_unk if the
/// corresponding D-amino acid cannot be determined).
//////////////////////////////////////////////////////////
AA
get_D_equivalent( AA aa );

///////////////////////////////////////////////////////
/// @brief give a enum type and return the string name
///////////////////////////////////////////////////////
std::string
name_from_aa( AA  aa );

///////////////////////////////////////////////////////
/// @brief give a enum type and return the string name
///////////////////////////////////////////////////////
char
oneletter_code_from_aa( AA aa );

///////////////////////////////////////////////////////////
/// @brief give a 1 letter code and return the string name
///////////////////////////////////////////////////////////
AA
aa_from_oneletter_code( char onelettercode );

bool
oneletter_code_specifies_aa( char onelettercode );

/// @brief input operator for AA enum type
std::istream & operator >>( std::istream & is, AA & aa );
/// @brief output operator for AA enum type
std::ostream & operator <<( std::ostream & os, AA const & aa );


} // chemical
} // core

#endif
