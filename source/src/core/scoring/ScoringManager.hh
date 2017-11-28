// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoringManager.hh
/// @brief  Scoring manager class header.
/// @details  The ScoringManager handles the lazy loading of data for each scoretype.  Note that data load
/// must be threadsafe.  For this, the utility::thread::safely_create_load_once_object_by_OP function is used.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- added thread-safety to lazily loaded data.

#ifndef INCLUDED_core_scoring_ScoringManager_hh
#define INCLUDED_core_scoring_ScoringManager_hh

// Unit headers
#include <core/scoring/ScoringManager.fwd.hh>

// Package headers
#include <core/scoring/AtomVDW.fwd.hh>
#include <core/scoring/CenRotEnvPairPotential.fwd.hh>
#include <core/scoring/CenHBPotential.fwd.hh>
#include <core/scoring/EnvPairPotential.fwd.hh>
#include <core/scoring/GenBornPotential.fwd.hh>
#include <core/scoring/HydroxylTorsionPotential.fwd.hh>
#include <core/scoring/MultipoleElecPotential.fwd.hh>
#include <core/scoring/SASAPotential.fwd.hh>
#include <core/scoring/MembranePotential.fwd.hh>
#include <core/scoring/Membrane_FAPotential.fwd.hh> //pba
#include <core/scoring/OmegaTether.fwd.hh>
#include <core/scoring/P_AA.fwd.hh>
#include <core/scoring/P_AA_ss.fwd.hh>
#include <core/scoring/PairEPotential.fwd.hh>
#include <core/scoring/PoissonBoltzmannPotential.fwd.hh>
#include <core/scoring/ProQPotential.fwd.hh>
#include <core/scoring/RamaPrePro.fwd.hh>
#include <core/scoring/Ramachandran.fwd.hh>
#include <core/scoring/Ramachandran2B.fwd.hh>
#include <core/scoring/P_AA_ABEGO3.fwd.hh>
#include <core/scoring/SecondaryStructurePotential.fwd.hh>
#include <core/scoring/SmoothEnvPairPotential.fwd.hh>
#include <core/scoring/UnfoldedStatePotential.fwd.hh>
#include <core/scoring/VdWTinkerPotential.fwd.hh>
#include <core/scoring/WaterAdductHBondPotential.fwd.hh>

#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.fwd.hh>
#include <core/scoring/netcharge_energy/NetChargeEnergySetup.fwd.hh>

#include <core/scoring/carbohydrates/CHIEnergyFunction.fwd.hh>
#include <core/scoring/carbohydrates/OmegaPreferencesFunction.fwd.hh>

#include <core/scoring/carbon_hbonds/CarbonHBondPotential.fwd.hh>

#include <core/scoring/disulfides/FullatomDisulfidePotential.fwd.hh>
#include <core/scoring/disulfides/CentroidDisulfidePotential.fwd.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.fwd.hh>

#include <core/scoring/elec/CPRepMapType.fwd.hh>

#include <core/scoring/UnfoldedStatePotential.fwd.hh>
#include <core/scoring/PoissonBoltzmannPotential.fwd.hh>
#include <core/scoring/SplitUnfoldedTwoBodyPotential.fwd.hh>

#include <core/scoring/dna/DNA_BasePotential.fwd.hh>
#include <core/scoring/dna/DNA_DihedralPotential.fwd.hh>
#include <core/scoring/dna/DNA_EnvPairPotential.fwd.hh>
#include <core/scoring/dna/DNABFormPotential.fwd.hh>
#include <core/scoring/dna/DNATorsionPotential.fwd.hh>
#include <core/scoring/dna/DirectReadoutPotential.fwd.hh>

#include <core/scoring/facts/FACTSPotential.fwd.hh>

#include <core/scoring/geometric_solvation/DatabaseOccSolEne.fwd.hh>

#include <core/scoring/interface_/DDPlookup.fwd.hh>

#include <core/scoring/membrane/MembraneData.hh>

#include <core/scoring/mm/MMLJLibrary.fwd.hh>
#include <core/scoring/mm/MMLJEnergyTable.fwd.hh>
#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/scoring/mm/MMBondLengthLibrary.fwd.hh>

#include <core/scoring/nv/NVlookup.fwd.hh>

#include <core/scoring/orbitals/OrbitalsLookup.fwd.hh>

#include <core/scoring/rna/RNA_AtomVDW.fwd.hh>
#include <core/scoring/rna/RNA_TorsionPotential.fwd.hh>
#include <core/scoring/rna/RNA_SuitePotential.fwd.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.fwd.hh>
#include <core/scoring/rna/RNP_LowResPotential.fwd.hh>
#include <core/scoring/rna/RNP_LowResStackData.fwd.hh>
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.fwd.hh>
#include <core/scoring/rna/data/RNA_DMS_Potential.fwd.hh>
#include <core/scoring/rna/data/RNA_DMS_LowResolutionPotential.fwd.hh>

#include <core/scoring/loop_graph/evaluator/SixDTransRotPotential.fwd.hh>

#if defined(WIN32) || defined(WIN_PYROSETTA)
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh> // WIN32 INCLUDE
#include <core/scoring/etable/Etable.hh>
#endif

#include <core/scoring/ScoreType.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/memb_etable/MembEtable.fwd.hh> //pba
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>

#include <core/chemical/mainchain_potential/MainchainScoreTable.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>


// C++ headers
#include <map>
#ifdef MULTI_THREADED
#include <mutex>
#include <utility/thread/ReadWriteMutex.hh>
#include <utility/thread/threadsafe_creation.hh> // Needed for OLDER_GXX_STDLIB define
#ifdef OLDER_GXX_STDLIB
#include <atomic>
#endif
#endif

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {

//singelton class
class ScoringManager : public utility::SingletonBase< ScoringManager >
{
public:
	friend class utility::SingletonBase< ScoringManager >;

public:

	/// @brief The ScoringManager acts as an EnergyMethodFactory.  All EnergyMethods must
	/// create a helper class, an EnergyMethodCreator class, that will respond to a call to
	/// its create_energy_method by returning a new instance of that EnergyMethod its helping.
	/// This Creator class must also register itself with the ScoringManager at load time and
	/// hand an instance of itself to the singleton ScoringManager instance.
	/// @details I don't think that this function is threadsafe (VKM, 20 July 2017), but it probably
	/// doesn't matter.  Factory registration presumably happens during Rosetta initialization,
	/// before any threads are spawned.
	void factory_register( methods::EnergyMethodCreatorOP creator );

	/// @brief Get a const instance of the PairEPotential.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan.
	PairEPotential const & get_PairEPotential() const;

	/// @brief Get a const instance of the GenBornPotential.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan.
	GenBornPotential const & get_GenBornPotential() const;

	/// @brief Get a const instance of the HydroxylTorsionPotential.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan.
	HydroxylTorsionPotential const & get_HydroxylTorsionPotential() const;

	/// @brief Get an instance of the VdWTinkerPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	VdWTinkerPotential const & get_VdWTinkerPotential() const;

	/// @brief Get an instance of the MultipoleElecPotential scoring object.
	/// @details Threadsafe creation, lazily loaded.
	/// @note The MultipoleElecPotential caches pose data to the global MultipleElecPotential object during scoring.
	/// As such, it is fundamentally NOT THREADSAFE!!!  (Only the creation of this object is threadsafe.)
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	MultipoleElecPotential const & get_MultipoleElecPotential( methods::EnergyMethodOptions const & options ) const;

	/// @brief Get an instance of the SASAPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	SASAPotential const & get_SASAPotential() const;

	/// @brief Get an instance of the FACTSPotential scoring object, by const owning pointer.
	/// @details Threadsafe and lazily loaded.
	/// @note The FACTSPotential caches pose data to the global FACTSPotential object during scoring.  As such, it is
	/// fundamentally NOT THREADSAFE!!!  (Only the creation of this object is threadsafe).
	/// @author Rewritten by Vikram K. Mulligan.
	FACTSPotential const & get_FACTSPotential() const;

	/// @brief Get an instance of the DNA_BasePotential scoring object, by const owning pointer.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan.
	dna::DNA_BasePotential const & get_DNA_BasePotential() const;

	/// @brief Get an instance of the Ramachandran scoring object, by const owning pointer.
	/// @details Threadsafe and lazily loaded.
	/// @note The Ramachandran object does lazily load custom cumulative distribution functions.  However, these are finite and
	/// governed by an enum.  In the MULTI_THREADED case, these are all loaded on object creation, which should get around thread-
	/// safety issues.
	/// @author Rewritten by Vikram K. Mulligan.
	RamachandranCOP get_Ramachandran_ptr() const;

	/// @brief Get a const instance of the Ramachandran scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note The Ramachandran object does lazily load custom cumulative distribution functions.  However, these are finite and
	/// governed by an enum.  In the MULTI_THREADED case, these are all loaded on object creation, which should get around thread-
	/// safety issues.
	/// @author Rewritten by Vikram K. Mulligan.
	Ramachandran const & get_Ramachandran() const;

	/// @brief Get an instance of the Ramachandran2B scoring object, by const owning pointer.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	Ramachandran2BCOP get_Ramachandran2B_ptr() const;

	/// @brief Get an instance of the Ramachandran2B scoring object, by const instance.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	Ramachandran2B const & get_Ramachandran2B() const;

	/// @brief Get an instance of the RamaPrePro scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	RamaPrePro const & get_RamaPrePro() const;

	/// @brief Get an instance of the P_AA_ABEGO3 scoring object.
	/// @details Threadsafe and lazily loaded.  Used by AbegoEnergy.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author imv@uw.edu
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	P_AA_ABEGO3 const & get_P_AA_ABEGO3() const;

	/// @brief Get an instance of the DNABFormPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	dna::DNABFormPotential const & get_DNABFormPotential() const;

	/// @brief Get an instance of the DNATorsionPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	dna::DNATorsionPotential const & get_DNATorsionPotential() const;

	/// @brief Get an instance of the OmegaTether scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	OmegaTether const & get_OmegaTether() const;

	/// @brief Get an instance of the SmoothEnvPairPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	SmoothEnvPairPotential const & get_SmoothEnvPairPotential() const;

	/// @brief Get an instance of the CenRotEnvPairPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	CenRotEnvPairPotential const & get_CenRotEnvPairPotential() const;

	/// @brief Get an instance of the CenHBPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	CenHBPotential const & get_CenHBPotential() const;

	/// @brief Get an instance of the EnvPairPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	EnvPairPotential const & get_EnvPairPotential() const;

	/// @brief Get an instance of the DNA_EnvPairPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	dna::DNA_EnvPairPotential const & get_DNA_EnvPairPotential() const;

	/// @brief Get an instance of the DNA_DihedralPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	dna::DNA_DihedralPotential const & get_DNA_DihedralPotential() const;

	/// @brief Get an instance of the SecondaryStructurePotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	SecondaryStructurePotential const & get_SecondaryStructurePotential() const;

	/// @brief Get an instance of the AtomVDW scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Each element in the atom_vdw_ map is now threadsafe and lazily loaded (independently).  Each of these
	/// objects is also threadsafe, as far as I can tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	AtomVDW const & get_AtomVDW( std::string const & atom_type_set_name ) const;

	/// @brief Get an instance of the RNA_AtomVDW scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	rna::RNA_AtomVDW const & get_RNA_AtomVDW() const;

	/// @brief Get an instance of the DatabaseOccSolEne scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Whatever atom type set name and min occ energy are passed to this function the FIRST time determine the
	/// object that gets created.  These parameters are unused in subsequent invocations.  This aside, the targeted object
	/// is threadsafe, as far as I can tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	geometric_solvation::DatabaseOccSolEne const & get_DatabaseOccSolEne( std::string const & atom_type_set_name, Real const & min_occ_energy ) const;

	/// @brief Get an instance of the CarbonHBondPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	carbon_hbonds::CarbonHBondPotential const & get_CarbonHBondPotential() const;

	/// @brief Get an instance of the RNA_SuitePotentialCOP scoring object, by const owning pointer.
	/// @details Threadsafe and lazily loaded.
	/// @note The RNA_SuitePotential caches pose-specific scoring data in the global instance of the RNA_SuitePotential object.  As such,
	/// it is fundamentally NOT THREADSAFE.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	rna::RNA_SuitePotentialCOP get_rna_suite_potential( bool const & calculate_suiteness_bonus, std::string const & suiteness_bonus ) const;

	/// @brief Get an instance of the SixDTransRotPotential scoring object, by const owning pointer.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	loop_graph::evaluator::SixDTransRotPotentialCOP get_LoopCloseSixDPotential( std::string const & database_file ) const;

	/// @brief Get an instance of the RNA_LowResolutionPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	rna::RNA_LowResolutionPotential const & get_RNA_LowResolutionPotential() const;

	/// @brief Get an instance of the RNP_LowResPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	rna::RNP_LowResPotential const & get_RNP_LowResPotential() const;

	/// @brief Get an instance of the RNP_LowResStackData scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	rna::RNP_LowResStackData const & get_RNP_LowResStackData() const;

	/// @brief Get an instance of the RNA_ChemicalShiftPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	rna::chemical_shift::RNA_ChemicalShiftPotential const & get_RNA_ChemicalShiftPotential() const;

	/// @brief Get an instance of the RNA_DMS_Potential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note The RNA_DMS_Potential itself is fundamentally NOT THREADSAFE!!!  Note that this function
	/// returns a non-const instance (which it shouldn't).
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	rna::data::RNA_DMS_Potential & get_RNA_DMS_Potential() const;

	/// @brief Get an instance of the RNA_DMS_LowResolutionPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note The RNA_DMS_LowResolutionPotential itself is fundamentally NOT THREADSAFE!!!  Note that this function
	/// returns a non-const instance (which it shouldn't).
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	rna::data::RNA_DMS_LowResolutionPotential & get_RNA_DMS_LowResolutionPotential() const;

	/// @brief Get an instance of the DirectReadoutPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	dna::DirectReadoutPotential const & get_DirectReadoutPotential() const;

	/// @brief Get an instance of the MMLJLibrary scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	mm::MMLJLibrary const & get_MMLJLibrary() const;

	/// @brief Get an instance of the MMLJEnergyTable scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	mm::MMLJEnergyTable const & get_MMLJEnergyTable() const;

	/// @brief Get an instance of the MMTorsionLibrary scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	mm::MMTorsionLibrary const & get_MMTorsionLibrary() const;

	/// @brief Get an instance of the MMBondAngleLibrary scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	mm::MMBondAngleLibrary const & get_MMBondAngleLibrary() const;

	/// @brief Get an instance of the MMBondLengthLibrary scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	mm::MMBondLengthLibrary const & get_MMBondLengthLibrary() const;

	/// @brief Get an instance of the NVlookup scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	nv::NVlookup const & get_NVLookupTable() const;

	/// @brief Get an instance of the OrbitalsLookup scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	core::scoring::orbitals::OrbitalsLookup const & get_OrbitalsLookupTable() const;

	/// @brief Get an instance of the DDPlookup scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	interface_::DDPlookup const & get_DDPLookupTable() const;

	/// @brief Get an instance of the P_AA scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	P_AA const & get_P_AA() const;

	/// @brief Get an instance of the P_AA_ss scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	P_AA_ss const & get_P_AA_ss() const;

	/// @brief Get an instance of the UnfoldedStatePotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	UnfoldedStatePotential const & get_UnfoldedStatePotential( std::string const & type ) const;

	/// @brief Get an instance of the WaterAdductHBondPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	WaterAdductHBondPotential const & get_WaterAdductHBondPotential() const;

	/// @brief Get an instance of the MembranePotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	MembranePotential const & get_MembranePotential() const;

	/// @brief Get an instance of the MembraneData scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	membrane::MembraneData const & get_MembraneData() const;

	/// @brief Get an instance of the Membrane_FAPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note The Membrane_FAPotential object is fundamentally NOT THREADSAFE!!!
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	Membrane_FAPotential const & get_Membrane_FAPotential() const; //pba

	/// @brief Get an instance of the ProQPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	ProQPotential const & get_ProQPotential() const;

	/// @brief Get an instance of the PoissonBoltzmannPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	PoissonBoltzmannPotential const & get_PoissonBoltzmannPotential() const;

	/// @brief Get an instance of the SplitUnfoldedTwoBodyPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	SplitUnfoldedTwoBodyPotential const & get_SplitUnfoldedTwoBodyPotential(std::string const & label_type,std::string const & value_type, std::string const & score_func_type ) const;

	/// @brief Get an instance of the FullatomDisulfidePotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	disulfides::FullatomDisulfidePotential const & get_FullatomDisulfidePotential() const;

	/// @brief Get an instance of the CentroidDisulfidePotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	disulfides::CentroidDisulfidePotential const & get_CentroidDisulfidePotential() const;

	/// @brief Get an instance of the DisulfideMatchingPotential scoring object.
	/// @details Threadsafe and lazily loaded.
	/// @note Targeted object is also threadsafe, to the best of my ability to tell.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	disulfides::DisulfideMatchingPotential const & get_DisulfideMatchingPotential() const;

	/// @brief Get an instance of the CHIEnergyFunction scoring object.
	/// @details Lazily loaded, but fundamentally NOT THREADSAFE.
	carbohydrates::CHIEnergyFunction const & get_CHIEnergyFunction(
		bool const setup_sampling_data = false,
		Real const & sampling_step_size = 0.1 ) const;

	/// @brief Get an instance of the OmegaPreferencesFunction scoring object.
	/// @details Lazily loaded, but fundamentally NOT THREADSAFE.
	carbohydrates::OmegaPreferencesFunction const & get_OmegaPreferencesFunction(
		bool const setup_sampling_data = false,
		core::Real const & sampling_step_size = 0.1 ) const;

	/// @brief Test if there is an EnergyMethod class defined for a
	/// given score type.
	/// @details I THINK that this is threadsafe (VKM, 20 July 2017).
	bool has_energy_method( ScoreType t ) const;

	/// @brief When a ScoreFunction the weight for a particular ScoreType set from 0
	/// to some non-zero value, it will request an instance of the EnergyMethod class
	/// that is responsible for calculating that ScoreType.  The ScoringManager responds
	/// to that request by asking the EnergyMethodCreator that has claimed responsibility
	/// for this ScoreType for a new instance.  EnergyMethodCreators must first have
	/// registered themselves with the ScoringManager.  This should have been done at
	/// load time, using a static-variable-initialization function call.
	/// See src/core/scoring/etable/EtableEnergy.cc for an example of how the
	/// EtableEnergyCreator class registers itself with the ScoringManager.
	/// @details I THINK that this is threadsafe (VKM, 20 July 2017).
	methods::EnergyMethodOP
	energy_method( ScoreType const & t, methods::EnergyMethodOptions const & options ) const;

	/// @brief make etable for extra partially softies
	/// @details Make etable for extra partial softies, pilot app r_play_with_etables does not really work anymore
	/// the etables it added will be somehow cleared or overwriten, so I do it here now
	///
	/// table_id: i.e. FA_STANDARD_SOFT40,
	/// the number in the end is a percentage, lj_radius are given for every 5% softie
	/// from 5% to 95% in database/chemical/atom_type_sets/fa_standard/extras/extra_soft_rep_params.txt
	/// FA_STANDARD_SOFT50 would be halfway between normal softrep and normal hardrep.
	/// @note This is threadsafe but ugly -- it creates a new object every time it's invoked.
	etable::EtableOP
	make_partially_soft_etable( std::string const & name, etable::EtableOptions etable_options ) const;

	/// @brief Add a new membrane energy table to the membrane energy tables map.
	/// @details Made threadsafe on 25 July 2017.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	void
	add_memb_etable( std::string const & name, etable::MembEtableOP etable );

	/// @brief Get a membrane energy table from the membrane energy tables map.
	/// @details Made threadsafe on 25 July 2017.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu)
	etable::MembEtableCAP
	memb_etable( std::string const & table_id ) const;

	/// @brief Request an etable specifying an EnergyMethodOptions object; internally
	/// this will retrieve the EtableOptions object, and invoke the EtableOptions
	/// version of this function.
	/// @details Threadsafe, since this calls the threadsafe etable( EtableOptions) function, below.
	/// @author Hahnbeom Park (new logic for etable)
	/// @note One weak point of this method is that, change in frequently called
	/// but less relevant options such as weights, will invoke another etable
	/// construction.
	etable::EtableCAP
	etable( methods::EnergyMethodOptions const & options_in ) const;

	/// @brief Request an etable specifying an EtableOptions; internally this will
	/// query the ScoringManager's map from EtableOptions to Etables for the desired
	/// Etable, and construct a new one if needed.
	/// @details Made threadsafe on 25 July 2017.
	/// @author Thread-safety added by Vikram K. Mulligan (vmullig@uw.edu).
	etable::EtableCAP
	etable( etable::EtableOptions const & options_in ) const;

	/// @brief Create and return an etable specified only by the etable_type of the
	/// etable::EtableOptions class.  This, internally, will create an EtableOptions object,
	/// initialized from the command line, set the etable_type of this object, and then
	/// invoke the etable( EtableOptions ) method.
	/// @details Threadsafe, since this calls the threadsafe etable( EtableOptions ) function, above.
	etable::EtableCAP
	etable( std::string const & etable_id ) const;

	/// @brief Get an owning pointer to data used by the FA_ElecEnergy in beta_nov15 mode.
	/// @details If the data have not been loaded, this loads the data (lazy loading).  Lazy loading is now threadsafe.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::scoring::elec::CPRepMapTypeCOP get_cp_rep_map_byname() const;

	/// @brief Get a vector of owning pointers to data used by the AACompositionEnergy score term.
	/// @details If this vector has not yet been populated, this loads the data from disk (lazy loading).
	/// @note The lazy loading has been made threadsafe, as of the wee hours of 26 July 2017.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1< core::scoring::aa_composition_energy::AACompositionEnergySetupOP > get_cloned_aa_comp_setup_helpers( core::scoring::methods::EnergyMethodOptions const &options ) const;

	/// @brief Get a vector of owning pointers to data used by the NetChargeEnergy score term.
	/// @details If this vector has not yet been populated, this loads the data from disk (lazy loading).
	/// @note The lazy loading has been made threadsafe.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1< core::scoring::netcharge_energy::NetChargeEnergySetupOP > get_cloned_netcharge_setup_helpers( core::scoring::methods::EnergyMethodOptions const &options ) const;

	/// @brief Get a particular MainchainScoreTable for the rama_prepro score term, for a particular residue type.
	/// @details If this has not yet been populated, loads the data from disk (lazy loading)
	/// in a threadsafe manner.
	/// @note Each restype stores separate tables for general and pre-proline scoring.  The prepro_table parameter determines
	/// whether we're loading the default scoring table or the version for residues occurring before a proline.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::chemical::mainchain_potential::MainchainScoreTableCOP
	get_rama_prepro_mainchain_torsion_potential(
		core::chemical::ResidueTypeCOP restype,
		bool const use_polycubic_interpolation,
		bool const prepro_table
	) const;

private:

	//private constructor
	ScoringManager();
	~ScoringManager();


private:
	//Creation functions, needed for threadsafe lazy loading.
	//Note that these trigger database access, and are NOT to be used repeatedly.

	/// @brief Create an instance of the PairEPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static PairEPotentialOP create_pairE_potential_instance();

	/// @brief Create an instance of the GenBornPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static GenBornPotentialOP create_genborn_instance();

	/// @brief Create an instance of the VdWTinkerPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static VdWTinkerPotentialOP create_vdw_tinker_potential_instance();

	/// @brief Create an instance of the HydroxylTorsionPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static MultipoleElecPotentialOP create_multipole_elec_instance( methods::EnergyMethodOptions const & options );

	/// @brief Create an instance of the SASAPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static SASAPotentialOP create_sasa_potential_instance();

	/// @brief Create an instance of the HydroxylTorsionPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static HydroxylTorsionPotentialOP create_hxl_potential_instance();

	/// @brief Create an instance of the FactsPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static FACTSPotentialOP create_facts_potential_instance();

	/// @brief Create an instance of the FactsPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static dna::DNA_BasePotentialOP create_dnabase_potential_instance();

	/// @brief Create an instance of the Ramachandran object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static RamachandranOP create_rama_instance();

	/// @brief Create an instance of the Ramachandran2B object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static Ramachandran2BOP create_rama2b_instance();

	/// @brief Create an instance of the RamaPrePro object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static RamaPreProOP create_ramapp_instance();

	/// @brief Create an instance of the P_AA_ABEGO3 object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static P_AA_ABEGO3_OP create_p_aa_abego3_instance();

	/// @brief Create an instance of the P_AA object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static P_AAOP create_p_aa_instance();

	/// @brief Create an instance of the P_AA_ss object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static P_AA_ssOP create_p_aa_ss_instance();

	/// @brief Create an instance of the DNABFormPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static dna::DNABFormPotentialOP create_dna_bform_potential_instance();

	/// @brief Create an instance of the DNATorsionPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static dna::DNATorsionPotentialOP create_dna_torsion_potential_instance();

	/// @brief Create an instance of the OmegaTether object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static OmegaTetherOP create_omegatether_instance();

	/// @brief Create an instance of the SmoothEnvPairPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static SmoothEnvPairPotentialOP create_smoothenvpair_instance();

	/// @brief Create an instance of the CenRotEnvPairPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static CenRotEnvPairPotentialOP create_cenrotenvpair_instance();

	/// @brief Create an instance of the CenHBPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static CenHBPotentialOP create_cenhbpotential_instance();

	/// @brief Create an instance of the EnvPairPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static EnvPairPotentialOP create_envpairpotential_instance();

	/// @brief Create an instance of the DNA_EnvPairPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static dna::DNA_EnvPairPotentialOP create_dna_envpairpotential_instance();

	/// @brief Create an instance of the DNA_DihedralPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static dna::DNA_DihedralPotentialOP create_dna_dihedralpotential_instance();

	/// @brief Create an instance of the SecondaryStructurePotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static SecondaryStructurePotentialOP create_secondarystructurepotential_instance();

	/// @brief Create an instance of the AtomVDW object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static AtomVDWOP create_atomvdw_instance( std::string const & atom_type_set_name );

	/// @brief Create an instance of the RNA_AtomVDW object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static rna::RNA_AtomVDWOP create_rna_atomvdw_instance();

	/// @brief Create an instance of the DatabaseOccSolEne object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static geometric_solvation::DatabaseOccSolEneOP create_database_occsolene_instance( std::string const & atom_type_set_name, core::Real const & min_occ_energy );

	/// @brief Create an instance of the CarbonHBondPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static carbon_hbonds::CarbonHBondPotentialOP create_carbon_hbond_potential_instance();

	/// @brief Create an instance of the RNA_SuitePotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static rna::RNA_SuitePotentialOP create_rna_suitepotential_instance( bool const & calculate_suiteness_bonus, std::string const & suiteness_bonus );

	/// @brief Create an instance of the SixDTransRotPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static loop_graph::evaluator::SixDTransRotPotentialOP create_sixdtransrotpotential_instance( std::string const & database_file );

	/// @brief Create an instance of the RNA_LowResolutionPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static rna::RNA_LowResolutionPotentialOP create_rna_lowresolutionpotential_instance();

	/// @brief Create an instance of the RNP_LowResPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static rna::RNP_LowResPotentialOP create_rnp_lowrespotential_instance();

	/// @brief Create an instance of the RNP_LowResStackData object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static rna::RNP_LowResStackDataOP create_rnp_lowresstackdata_instance();

	/// @brief Create an instance of the RNA_ChemicalShiftPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static rna::chemical_shift::RNA_ChemicalShiftPotentialOP create_rna_chemshiftpotential_instance();

	/// @brief Create an instance of the RNA_DMS_Potential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static rna::data::RNA_DMS_PotentialOP create_rna_dms_potential_instance();

	/// @brief Create an instance of the RNA_DMS_LowResolutionPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static rna::data::RNA_DMS_LowResolutionPotentialOP create_rna_dms_lowrespotential_instance();

	/// @brief Create an instance of the DirectReadoutPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static dna::DirectReadoutPotentialOP create_dna_directreadoutpotential_instance();

	/// @brief Create an instance of the MMLJLibrary object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static mm::MMLJLibraryOP create_mm_lj_library_instance();

	/// @brief Create an instance of the MMLJEnergyTable object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static mm::MMLJEnergyTableOP create_mm_lj_energy_table_instance();

	/// @brief Create an instance of the MMTorsionLibrary object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static mm::MMTorsionLibraryOP create_mm_torsion_library_instance();

	/// @brief Create an instance of the MMBondAngleLibrary object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static mm::MMBondAngleLibraryOP create_mm_bondangle_library_instance();

	/// @brief Create an instance of the MMBondLengthLibrary object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static mm::MMBondLengthLibraryOP create_mm_bondlength_library_instance();

	/// @brief Create a (default) instance of the CHIEnergyFunction object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	static carbohydrates::CHIEnergyFunctionOP create_chi_energy_function_instance();

	/// @brief Create a (default) instance of the OmegaPreferencesFunction object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	static carbohydrates::OmegaPreferencesFunctionOP create_omega_preferences_function_instance();

	/// @brief Create an instance of the NVlookup object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static nv::NVlookupOP create_nvlookup_instance();

	/// @brief Create an instance of the OrbitalsLookup object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static orbitals::OrbitalsLookupOP create_orbitals_lookup_instance();

	/// @brief Create an instance of the DDPlookup object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static interface_::DDPlookupOP create_ddp_lookup_instance();

	/// @brief Create an instance of the UnfoldedStatePotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static UnfoldedStatePotentialOP create_unfolded_state_potential_instance( std::string const & type );

	/// @brief Create an instance of the WaterAdductHBondPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static WaterAdductHBondPotentialOP create_water_adduct_hbond_potential_instance();

	/// @brief Create an instance of the MembranePotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static MembranePotentialOP create_membrane_potential_instance();

	/// @brief Create an instance of the MembraneData object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static membrane::MembraneDataOP create_membrane_data_instance();

	/// @brief Create an instance of the Membrane_FAPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static Membrane_FAPotentialOP create_membrane_fa_potential_instance();

	/// @brief Create an instance of the ProQPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static ProQPotentialOP create_proq_potential_instance();

	/// @brief Create an instance of the PoissonBoltzmannPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static PoissonBoltzmannPotentialOP create_poisson_boltzmann_potential_instance();

	/// @brief Create an instance of the SplitUnfoldedTwoBodyPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static SplitUnfoldedTwoBodyPotentialOP create_split_unfolded_2body_potential_instance(std::string const & label_type,std::string const & value_type, std::string const & score_func_type);

	/// @brief Create an instance of the FullatomDisulfidePotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static disulfides::FullatomDisulfidePotentialOP create_fullatom_disulfide_potential_instance();

	/// @brief Create an instance of the CentroidDisulfidePotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static disulfides::CentroidDisulfidePotentialOP create_centroid_disulfide_potential_instance();

	/// @brief Create an instance of the DisulfideMatchingPotential object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static disulfides::DisulfideMatchingPotentialOP create_disulfide_matching_potential_instance();

	/// @brief Create an instance of a MainchainScoreTable, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.  This one is kind of silly, since it just
	/// returns a MainchainScoreTableOP that is passed in.  Still needed for threadsafe creation, though.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static core::chemical::mainchain_potential::MainchainScoreTableOP create_mainchain_scoretable_instance( core::chemical::mainchain_potential::MainchainScoreTableOP table_in );

	/// @brief Create an instance of a MembEtable, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.  This one is kind of silly, since it just
	/// returns a MembEtableOP that is passed in.  Still needed for threadsafe creation, though.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static etable::MembEtableOP create_memb_etable_instance_silly( etable::MembEtableOP table_in );

	/// @brief Create an instance of a MembEtable, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.  This is the non-silly version that builds
	/// the object based on a string.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static etable::MembEtableOP create_memb_etable_instance( std::string const & table_id );

	/// @brief Create an instance of an Etable, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.  This one is kind of silly, since it just
	/// returns an EtableOP that is passed in.  Still needed for threadsafe creation, though.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static etable::EtableOP create_etable_instance( etable::EtableOP table_in );

	/// @brief Create an instance of an AACompositionEnergySetup object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static core::scoring::aa_composition_energy::AACompositionEnergySetupOP create_aa_composition_energy_setup_instance( std::string const &filename );

	/// @brief Create an instance of an NetChargeEnergySetup object, by owning pointer.
	/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
	/// @note Not intended for use outside of ScoringManager.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	static core::scoring::netcharge_energy::NetChargeEnergySetupOP create_netcharge_energy_setup_instance( std::string const &filename );

private:

#ifdef MULTI_THREADED
	// Mutexes -- only exist in the threaded compilations.
	mutable std::mutex pairE_mutex_;
	mutable std::mutex genborn_mutex_;
	mutable std::mutex hxl_potential_mutex_;
	mutable std::mutex vdw_tinker_mutex_;
	mutable std::mutex multipole_elec_mutex_;
	mutable std::mutex sasa_potential_mutex_;
	mutable std::mutex facts_mutex_;
	mutable std::mutex dnabase_mutex_;
	mutable std::mutex rama_mutex_;
	mutable std::mutex rama2b_mutex_;
	mutable std::mutex rama_pp_mutex_;
	mutable std::mutex p_aa_abego3_mutex_;
	mutable std::mutex dnabform_mutex_;
	mutable std::mutex dnatorsion_mutex_;
	mutable std::mutex omegatether_mutex_;
	mutable std::mutex smoothenvpair_mutex_;
	mutable std::mutex cenrotenvpair_mutex_;
	mutable std::mutex cenhb_mutex_;
	mutable std::mutex envpair_mutex_;
	mutable std::mutex dnaenvpair_mutex_;
	mutable std::mutex dnadihedral_mutex_;
	mutable std::mutex secstruct_mutex_;
	mutable utility::thread::ReadWriteMutex atomvdw_mutex_;
	mutable std::mutex rna_atomvdw_mutex_;
	mutable std::mutex database_occ_sol_mutex_;
	mutable std::mutex carbonhbond_mutex_;
	mutable utility::thread::ReadWriteMutex rna_suite_mutex_;
	mutable utility::thread::ReadWriteMutex loopclose_sixdtransrot_mutex_;
	mutable std::mutex rna_lowres_mutex_;
	mutable std::mutex rnp_lowres_mutex_;
	mutable std::mutex rnp_lowresstack_mutex_;
	mutable std::mutex rna_chemshift_mutex_;
	mutable std::mutex rna_dms_mutex_;
	mutable std::mutex rna_dms_lowres_mutex_;
	mutable std::mutex dna_directreadout_mutex_;
	mutable std::mutex mm_lj_library_mutex_;
	mutable std::mutex mm_lj_energytable_mutex_;
	mutable std::mutex mm_torsionlibrary_mutex_;
	mutable std::mutex mm_bondanglelibrary_mutex_;
	mutable std::mutex mm_bondlengthlibrary_mutex_;
	mutable std::mutex nv_lookup_mutex_;
	mutable std::mutex orbitals_lookup_mutex_;
	mutable std::mutex ddp_lookup_mutex_;
	mutable std::mutex p_aa_mutex_;
	mutable std::mutex p_aa_ss_mutex_;
	mutable std::mutex unfoldedstate_mutex_;
	mutable std::mutex wateradduct_hbond_mutex_;
	mutable std::mutex membranepot_mutex_;
	mutable std::mutex membranedata_mutex_;
	mutable std::mutex membrane_fapot_mutex_;
	mutable std::mutex proq_mutex_;
	mutable std::mutex poissonboltzman_mutex_;
	mutable std::mutex splitunfolded_2body_mutex_;
	mutable std::mutex fa_disulf_potential_mutex_;
	mutable std::mutex cent_disulf_potential_mutex_;
	mutable std::mutex disulf_matching_potential_mutex_;
	mutable std::mutex carb_chienergy_mutex_;
	mutable std::mutex carb_omegapref_mutex_;
	mutable utility::thread::ReadWriteMutex memb_etable_mutex_;
	mutable utility::thread::ReadWriteMutex etable_mutex_;
	mutable std::mutex cp_rep_map_mutex_;
	mutable utility::thread::ReadWriteMutex aa_comp_mutex_;
	mutable utility::thread::ReadWriteMutex netcharge_mutex_;
#ifdef OLDER_GXX_STDLIB
	// Older versions of GCC lack the atomic_load and atomic_store functions for std::shared_ptr,
	// so we also need std::atomic_bools as a workaround
	mutable std::atomic_bool pairE_bool_;
	mutable std::atomic_bool genborn_bool_;
	mutable std::atomic_bool hxl_potential_bool_;
	mutable std::atomic_bool vdw_tinker_bool_;
	mutable std::atomic_bool multipole_elec_bool_;
	mutable std::atomic_bool sasa_potential_bool_;
	mutable std::atomic_bool facts_bool_;
	mutable std::atomic_bool dnabase_bool_;
	mutable std::atomic_bool rama_bool_;
	mutable std::atomic_bool rama2b_bool_;
	mutable std::atomic_bool rama_pp_bool_;
	mutable std::atomic_bool p_aa_abego3_bool_;
	mutable std::atomic_bool dnabform_bool_;
	mutable std::atomic_bool dnatorsion_bool_;
	mutable std::atomic_bool omegatether_bool_;
	mutable std::atomic_bool smoothenvpair_bool_;
	mutable std::atomic_bool cenrotenvpair_bool_;
	mutable std::atomic_bool cenhb_bool_;
	mutable std::atomic_bool envpair_bool_;
	mutable std::atomic_bool dnaenvpair_bool_;
	mutable std::atomic_bool dnadihedral_bool_;
	mutable std::atomic_bool secstruct_bool_;
	//mutable std::atomic_bool atomvdw_bool_;
	mutable std::atomic_bool rna_atomvdw_bool_;
	mutable std::atomic_bool database_occ_sol_bool_;
	mutable std::atomic_bool carbonhbond_bool_;
	//mutable std::atomic_bool rna_suite_bool_;
	//mutable std::atomic_bool loopclose_sixdtransrot_bool_;
	mutable std::atomic_bool rna_lowres_bool_;
	mutable std::atomic_bool rnp_lowres_bool_;
	mutable std::atomic_bool rnp_lowresstack_bool_;
	mutable std::atomic_bool rna_chemshift_bool_;
	mutable std::atomic_bool rna_dms_bool_;
	mutable std::atomic_bool rna_dms_lowres_bool_;
	mutable std::atomic_bool dna_directreadout_bool_;
	mutable std::atomic_bool mm_lj_library_bool_;
	mutable std::atomic_bool mm_lj_energytable_bool_;
	mutable std::atomic_bool mm_torsionlibrary_bool_;
	mutable std::atomic_bool mm_bondanglelibrary_bool_;
	mutable std::atomic_bool mm_bondlengthlibrary_bool_;
	mutable std::atomic_bool nv_lookup_bool_;
	mutable std::atomic_bool orbitals_lookup_bool_;
	mutable std::atomic_bool ddp_lookup_bool_;
	mutable std::atomic_bool p_aa_bool_;
	mutable std::atomic_bool p_aa_ss_bool_;
	mutable std::atomic_bool unfoldedstate_bool_;
	mutable std::atomic_bool wateradduct_hbond_bool_;
	mutable std::atomic_bool membranepot_bool_;
	mutable std::atomic_bool membranedata_bool_;
	mutable std::atomic_bool membrane_fapot_bool_;
	mutable std::atomic_bool proq_bool_;
	mutable std::atomic_bool poissonboltzman_bool_;
	mutable std::atomic_bool splitunfolded_2body_bool_;
	mutable std::atomic_bool fa_disulf_potential_bool_;
	mutable std::atomic_bool cent_disulf_potential_bool_;
	mutable std::atomic_bool disulf_matching_potential_bool_;
	mutable std::atomic_bool carb_chienergy_bool_;
	mutable std::atomic_bool carb_omegapref_bool_;
	mutable std::atomic_bool cp_rep_map_bool_;
#endif //OLDER_GXX_STDLIB
#endif //MULTI_THREADED
	// WARNING -- if you add something here don't forget to initialize to 0 in the constructor
	mutable VdWTinkerPotentialOP vdw_tinker_potential_;
	mutable PairEPotentialOP pairE_potential_;
	mutable RamachandranOP rama_;
	mutable Ramachandran2BOP rama2b_;
	mutable RamaPreProOP rama_pp_;
	mutable P_AA_ABEGO3_OP paa_abego3_;
	mutable OmegaTetherOP omega_;
	mutable EnvPairPotentialOP env_pair_potential_;
	mutable SmoothEnvPairPotentialOP smooth_env_pair_potential_;
	mutable CenRotEnvPairPotentialOP cen_rot_pair_potential_;
	mutable CenHBPotentialOP cen_hb_potential_;
	mutable SecondaryStructurePotentialOP secondary_structure_potential_;
	mutable std::map< std::string, AtomVDWOP > atom_vdw_;
	mutable rna::RNA_AtomVDWOP rna_atom_vdw_;
	mutable geometric_solvation::DatabaseOccSolEneOP occ_hbond_sol_database_;
	mutable dna::DirectReadoutPotentialOP dna_dr_potential_;
	mutable mm::MMLJLibraryOP mm_lj_library_;
	mutable mm::MMLJEnergyTableOP mm_lj_energy_table_;
	mutable mm::MMTorsionLibraryOP mm_torsion_library_;
	mutable mm::MMBondAngleLibraryOP mm_bondangle_library_;
	mutable mm::MMBondLengthLibraryOP mm_bondlength_library_;
	mutable dna::DNA_EnvPairPotentialOP dna_env_pair_potential_;
	mutable dna::DNA_DihedralPotentialOP dna_dihedral_potential_;
	mutable dna::DNABFormPotentialOP dnabform_;
	mutable dna::DNATorsionPotentialOP dna_torsion_potential_;
	mutable dna::DNA_BasePotentialOP DNA_base_potential_;
	mutable carbon_hbonds::CarbonHBondPotentialOP carbon_hbond_potential_;
	mutable rna::RNA_LowResolutionPotentialOP rna_low_resolution_potential_;
	mutable rna::RNP_LowResPotentialOP rnp_low_res_potential_;
	mutable rna::RNP_LowResStackDataOP rnp_low_res_stack_data_;
	mutable rna::chemical_shift::RNA_ChemicalShiftPotentialOP rna_chemical_shift_potential_;
	mutable rna::data::RNA_DMS_PotentialOP rna_dms_potential_;
	mutable rna::data::RNA_DMS_LowResolutionPotentialOP rna_dms_low_resolution_potential_;
	mutable std::map< std::pair< bool, std::string >, rna::RNA_SuitePotentialOP > rna_suite_potential_;
	mutable std::map< std::string, loop_graph::evaluator::SixDTransRotPotentialOP > loop_close_six_d_potential_;
	mutable P_AAOP p_aa_;
	mutable P_AA_ssOP p_aa_ss_;
	mutable WaterAdductHBondPotentialOP water_adduct_hbond_potential_;
	mutable GenBornPotentialOP gen_born_potential_;
	mutable HydroxylTorsionPotentialOP hxl_tors_potential_;
	mutable MultipoleElecPotentialOP multipole_elec_potential_;
	mutable SASAPotentialOP sasa_potential_;
	mutable FACTSPotentialOP facts_potential_;
	mutable disulfides::FullatomDisulfidePotentialOP fa_disulfide_potential_;
	mutable disulfides::CentroidDisulfidePotentialOP cen_disulfide_potential_;
	mutable disulfides::DisulfideMatchingPotentialOP disulfide_matching_potential_;
	mutable MembranePotentialOP membrane_potential_;
	mutable membrane::MembraneDataOP mp_base_potential_;
	mutable Membrane_FAPotentialOP membrane_fapotential_; //pba
	mutable ProQPotentialOP ProQ_potential_;
	mutable PoissonBoltzmannPotentialOP PB_potential_;
	mutable SplitUnfoldedTwoBodyPotentialOP sutbp_;
	//ReferenceEnergyPotential referenceEnergyPotential_;
	mutable UnfoldedStatePotentialOP unf_state_;
	mutable carbohydrates::CHIEnergyFunctionOP CHI_energy_function_;
	mutable carbohydrates::OmegaPreferencesFunctionOP carbohydrate_omega_preferences_function_;

	mutable nv::NVlookupOP NV_lookup_table_;
	mutable orbitals::OrbitalsLookupOP orbitals_lookup_table_;
	mutable interface_::DDPlookupOP DDP_lookup_table_;
	// data

	// new map for etables using EtableOptions as key
	mutable std::map< etable::EtableOptions, etable::EtableOP > etables_by_options_;

	mutable std::map< std::string, etable::MembEtableOP > memb_etables_; //pba

	/// @brief Cached data used by FA_ElecEnergy with beta_nov15.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	mutable core::scoring::elec::CPRepMapTypeOP cp_rep_map_byname_;

	/// @brief Cached data used by the AACompositionEnergy.
	/// @details Accessor function has been made threadsafe (as of the wee hours of 26 July 2017).
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	mutable std::map< std::string, core::scoring::aa_composition_energy::AACompositionEnergySetupOP > aa_composition_setup_helpers_;

	/// @brief Cached data used by the NetChargeEnergy.
	/// @details Accessor function has been made threadsafe.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	mutable std::map< std::string, core::scoring::netcharge_energy::NetChargeEnergySetupOP > netcharge_setup_helpers_;

	/// @brief Cached mainchain torsion potentials, used by rama_prepro.
	/// @details This one is for potentials for residues NOT occurring before proline.
	/// @note The lazy-loading accessor function has been made threadsafe.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	mutable std::map< std::string, core::chemical::mainchain_potential::MainchainScoreTableOP > rama_prepro_mainchain_potentials_;

	/// @brief Cached mainchain torsion potentials, used by rama_prepro.
	/// @details This one is for potentials for residues occurring before proline.
	/// @note The lazy-loading accessor function has been made threadsafe.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	mutable std::map< std::string, core::chemical::mainchain_potential::MainchainScoreTableOP > rama_prepro_mainchain_potentials_beforeproline_;

#ifdef MULTI_THREADED
	/// @brief Mutex for the map of rama_prepro potentials (for non-prepro positions).
	mutable utility::thread::ReadWriteMutex rama_prepro_mainchain_potentials_mutex_;

	/// @brief Mutex for the map of rama_prepro potentials (for prepro positions).
	mutable utility::thread::ReadWriteMutex rama_prepro_mainchain_potentials_beforeproline_mutex_;
#endif

	/// @brief The map of ( score type enum -> EnergyMethodCreatorOP ).
	/// @details Actually a simple vector (since the key is a 1-based, continuous enum).  Not threadsafe,
	/// but it doesn't really need to be, since it's initialized once by a single thread and never subsequently
	/// modified.
	utility::vector1< methods::EnergyMethodCreatorOP > method_creator_map_;


};
} // namespace core
} // namespace scoring


#endif
