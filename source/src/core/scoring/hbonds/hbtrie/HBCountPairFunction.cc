// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/HBCountPairFunction.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/hbonds/hbtrie/HBCountPairFunction.hh>

// Package Headers
#include <core/scoring/trie/trie_vs_trie.hh>
#include <core/scoring/trie/trie_vs_path.hh>
#include <core/scoring/hbonds/hbtrie/HBCPData.hh>
#include <core/scoring/hbonds/hbtrie/HBAtom.hh>
#include <core/scoring/hbonds/HBondEnergy.hh>
#include <core/scoring/lkball/lkbtrie/LKBAtom.fwd.hh>
#include <core/scoring/lkball/lkbtrie/LKBTrieEvaluator.hh>

// STL Headers
#include <iostream>

#include <utility/exit.hh>

#include <core/scoring/trie/RotamerTrie.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace hbonds {
namespace hbtrie {

using namespace trie;
using namespace etable::etrie;

HBCountPairFunction::~HBCountPairFunction() {}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with Etable Energy" );
}

/////////////////////////// EtableEnergy -- analytic evaluation //////////////////////////////////

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/

)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/

)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)

{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)

{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)

{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with CoarseEtable Energy" );
}

// HBONDS
void
HBCountPairFunction::resolve_trie_vs_trie(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData > const & trie1,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData > const & trie2,
	hbonds::HBondEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData > const & trie1,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData > const & trie2,
	hbonds::HBondEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


/// Hack Elec E
void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const &,
	elec::FA_ElecEnergy const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	elec::FA_ElecEnergy const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	elec::FA_ElecEnergy const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	elec::FA_ElecEnergy const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const &,
	elec::FA_ElecEnergy const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	elec::FA_ElecEnergy const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	elec::FA_ElecEnergy const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	elec::FA_ElecEnergy const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const &,
	elec::FA_ElecEnergy const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const &,
	elec::FA_ElecEnergy const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const &,
	elec::FA_ElecEnergy const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & ,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & ,
	elec::FA_ElecEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & ,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & ,
	elec::FA_ElecEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & ,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & ,
	elec::FA_ElecEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & ,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & ,
	elec::FA_ElecEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & ,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & ,
	elec::FA_ElecEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & ,
	trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & ,
	elec::FA_ElecEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

///////////////////////////////// lkball ////////////////////////////////////
void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const &,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const &,
	lkball::lkbtrie::LKBTrieEvaluator const &,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< core::PackerEnergy > &)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & ,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & ,
	lkball::lkbtrie::LKBTrieEvaluator const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & ,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & ,
	lkball::lkbtrie::LKBTrieEvaluator const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & ,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & ,
	lkball::lkbtrie::LKBTrieEvaluator const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & ,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & ,
	lkball::lkbtrie::LKBTrieEvaluator const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & ,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & ,
	lkball::lkbtrie::LKBTrieEvaluator const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & ,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & ,
	lkball::lkbtrie::LKBTrieEvaluator const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with FA_ElecEnergy" );
}



///////////////////////////// MMLJEnergyInter ///////////////////////////////

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}


void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

void
HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message( "hbonds::hbtrie::HBCountPairFunction type resolution reached with MMLJEnergyInter" );
}

/////////////////////////////// VDW Energy ////////////////////////////

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_trie reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}

void HBCountPairFunction::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message( "etable::etrie::HBCountPairFunction::resolve_trie_vs_path reached with VDW_Energy" );
}


void
HBCountPairFunction::print()
{
	std::cout << "HBCountPairFunction" << std::endl;
}

} // namespace hbtrie
} // namespace hbonds
} // namespace scoring
} // namespace core

