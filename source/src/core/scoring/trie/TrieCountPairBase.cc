// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/trie/TrieCountPairBase.hh
/// @brief  The functions in this file are required by the compiler, but represent paths that
///         should never be reached.  They are cases in which dynamic dispatch fails because
///         two incompatible RotamerTrie classes are used together, e.g. an HBond trie and
///         and Etable trie.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/trie/TrieCountPairBase.hh>

//#include <core/scoring/hbonds/HBondEnergy.hh>
//XRW_B_T1
//#include <core/scoring/etable/CoarseEtableEnergy.hh>
//XRW_E_T1
//#include <core/scoring/etable/EtableEnergy.hh>

#include <utility/exit.hh>

namespace core {
namespace scoring {
namespace trie {

TrieCountPairBase::~TrieCountPairBase() {}


void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}

//XRW_E_T1

// called when hbonds and etable tries get confused
void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	hbonds::HBondEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	hbonds::HBondEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


/// Hack Elec E

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}

// MMLJEnergyInter


void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}

// VDW Energy

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}

void
TrieCountPairBase::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & /*trie1*/,
	RotamerTrieBase const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit();
}


}
}
}

