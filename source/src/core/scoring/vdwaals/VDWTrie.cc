// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/vdwaals/VDWTrie.cc
/// @brief  Trie data structure for the low-resolution (centroid) repulsive energy
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/scoring/vdwaals/VDWTrie.hh>

// Package headers
#include <core/scoring/vdwaals/VDW_Energy.hh>
#include <core/scoring/trie/trie_vs_trie.hh>
#include <core/scoring/trie/trie_vs_path.hh>

#include <core/scoring/etable/etrie/CountPairData_1_1.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>

// Project headers
#include <core/conformation/Residue.hh>

// C++ headers
#include <iostream>

namespace core {
namespace scoring {
namespace vdwaals {

VDWAtom::VDWAtom() :
	xyz_( 0.0 ),
	is_hydrogen_( 0 ),
	atom_type_( 0 )
{}

VDWAtom::VDWAtom( core::conformation::Residue const & res, Size index ) :
	xyz_( res.xyz( index ) ),
	is_hydrogen_( res.atom_is_hydrogen( index ) ),
	atom_type_( res.atom( index ).type() )
{}

void VDWAtom::xyz( Vector const & setting )
{
	xyz_ = setting;
}

void
VDWAtom::atom_type( int setting ) {
	atom_type_ = setting;
}

void
VDWAtom::is_hydrogen( bool setting ) {
	is_hydrogen_ = setting ? 1 : 0;
}

/// @brief send a description of the atom to standard out
void VDWAtom::print() const
{
	print( std::cout );
}

/// @brief send a description of the atom to an output stream
void VDWAtom::print( std::ostream & os ) const
{
	os << "VDWAtom" <<  " ";
	os << "(" << xyz().x();
	os << ", " << xyz().y();
	os << ", " << xyz().z() << ") ";
	os << atom_type_ << ";";
}

VDWTrieCountPair1B::VDWTrieCountPair1B( Size res1_cp_data, Size res2_cp_data ) :
	res1_cp_data_(res1_cp_data),
	res2_cp_data_(res2_cp_data)
{}


///---------- TYPE RESOLUTION FUNCTIONS ----------///

using namespace core::scoring::etable;
using namespace core::scoring::etable::etrie;
using namespace core::scoring::trie;

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


///////////////////////// EtableEnergy -- analytic evaluation //////////////////////////////////

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

// HBONDS
void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & )
{
	utility_exit_with_message( "etable::etrie::VDWTrieCountPair1B::resolve_trie_vs_trie reached with HBondEnergy" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "etable::etrie::VDWTrieCountPair1B::resolve_trie_vs_path reached with HBondEnergy" );
}


/// Hack Elec E
void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

/////////////////////////////////// lkball /////////////////////////////////
void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}


void
VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

/////////////////////////////// MMLJEnergyInter ////////////////////////////

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_trie reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "vdwaals::VDWTrieCountPair1B::resolve_trie_vs_path reached with wrong score function" );
}

// vdw trie-vs-trie and trie-vs-path calls

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void VDWTrieCountPair1B::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


}
}
}
