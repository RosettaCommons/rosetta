// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoringManager.hh
/// @brief  Scoring manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_ScoringManager_hh
#define INCLUDED_core_scoring_ScoringManager_hh

// Unit headers
#include <core/scoring/ScoringManager.fwd.hh>

// Package headers
//#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/scoring/PairEPotential.fwd.hh>
#include <core/scoring/dna/DNA_BasePotential.fwd.hh>
#include <core/scoring/dna/DirectReadoutPotential.fwd.hh>
#include <core/scoring/Ramachandran.fwd.hh>
#include <core/scoring/Ramachandran2B.fwd.hh>
#include <core/scoring/OmegaTether.fwd.hh>
#include <core/scoring/EnvPairPotential.fwd.hh>
#include <core/scoring/SmoothEnvPairPotential.fwd.hh>
#include <core/scoring/CenHBPotential.fwd.hh>
#include <core/scoring/MembranePotential.fwd.hh>
#include <core/scoring/Membrane_FAPotential.fwd.hh> //pba
//#include <core/scoring/InterchainPotential.fwd.hh>
#include <core/scoring/SecondaryStructurePotential.fwd.hh>
#include <core/scoring/GenBornPotential.fwd.hh>
#include <core/scoring/AtomVDW.fwd.hh>
#include <core/scoring/carbon_hbonds/CarbonHBondPotential.fwd.hh>
#include <core/scoring/geometric_solvation/DatabaseOccSolEne.fwd.hh>
#include <core/scoring/rna/RNA_AtomVDW.fwd.hh>
#include <core/scoring/rna/RNA_TorsionPotential.fwd.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.fwd.hh>
#include <core/scoring/P_AA.fwd.hh>
#include <core/scoring/WaterAdductHBondPotential.fwd.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.fwd.hh>
#include <core/scoring/disulfides/CentroidDisulfidePotential.fwd.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.fwd.hh>
#include <core/scoring/UnfoldedStatePotential.fwd.hh>
#include <core/scoring/PoissonBoltzmannPotential.fwd.hh>

#include <core/scoring/mm/MMLJLibrary.fwd.hh>
#include <core/scoring/mm/MMLJEnergyTable.fwd.hh>
#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/scoring/mm/MMBondLengthLibrary.fwd.hh>

// AUTO-REMOVED #include <core/chemical/ResidueType.fwd.hh>
#include <core/scoring/nv/NVlookup.fwd.hh>
#include <core/scoring/orbitals/OrbitalsLookup.hh>
#include <core/scoring/interface/DDPlookup.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#ifdef WIN32
#include <core/scoring/methods/EnergyMethodCreator.hh> // WIN32 INCLUDE
#endif

// AUTO-REMOVED #include <core/scoring/EnergyMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/types.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/memb_etable/MembEtable.fwd.hh> //pba
//XRW_B_T1
//#include <core/coarse/CoarseEtable.fwd.hh>
//XRW_E_T1

// C++ headers
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

//singelton class
class ScoringManager
{
public:
	typedef core::scoring::mm::MMLJLibrary MMLJLibrary;
	typedef core::scoring::mm::MMLJEnergyTable MMLJEnergyTable;
	typedef core::scoring::mm::MMTorsionLibrary MMTorsionLibrary;
	typedef core::scoring::mm::MMBondAngleLibrary MMBondAngleLibrary;
	typedef core::scoring::mm::MMBondLengthLibrary MMBondLengthLibrary;

public:
	static ScoringManager * get_instance();

	void factory_register( methods::EnergyMethodCreatorOP creator );

	//P_AA const & get_P_AA() const;
	//ReferenceEnergyPotential const & get_ReferenceEnergyPotnential() const;
	PairEPotential const & get_PairEPotential() const;

	GenBornPotential const & get_GenBornPotential() const;

	dna::DNA_BasePotential const & get_DNA_BasePotential() const;

	///RotamerLibrary & get_RotamerLibrary() const;

	Ramachandran2B const & get_Ramachandran2B() const;

	Ramachandran const & get_Ramachandran() const;

	OmegaTether const & get_OmegaTether() const;

	SmoothEnvPairPotential const & get_SmoothEnvPairPotential() const;

	CenHBPotential const & get_CenHBPotential() const;

	EnvPairPotential const & get_EnvPairPotential() const;

	SecondaryStructurePotential const & get_SecondaryStructurePotential() const;

	AtomVDW const & get_AtomVDW( std::string const & atom_type_set_name ) const;

	rna::RNA_AtomVDW const & get_RNA_AtomVDW() const;

	geometric_solvation::DatabaseOccSolEne const &
	get_DatabaseOccSolEne(
		std::string const & atom_type_set_name,
		Real const & min_occ_energy
	) const;

	carbon_hbonds::CarbonHBondPotential const & get_CarbonHBondPotential() const;

	rna::RNA_LowResolutionPotential const & get_RNA_LowResolutionPotential() const;

	rna::RNA_TorsionPotential const & get_RNA_TorsionPotential() const;

	dna::DirectReadoutPotential const & get_DirectReadoutPotential() const;

	MMLJLibrary const & get_MMLJLibrary() const;

	MMLJEnergyTable const & get_MMLJEnergyTable() const;

	MMTorsionLibrary const & get_MMTorsionLibrary() const;

	MMBondAngleLibrary const & get_MMBondAngleLibrary() const;

	MMBondLengthLibrary const & get_MMBondLengthLibrary() const;

	nv::NVlookup const & get_NVLookupTable() const;
	core::scoring::orbitals::OrbitalsLookup const & get_OrbitalsLookupTable() const;



  interface::DDPlookup const & get_DDPLookupTable() const;

	P_AA const & get_P_AA() const;

	UnfoldedStatePotential const & get_UnfoldedStatePotential( std::string const & type ) const;

	WaterAdductHBondPotential const & get_WaterAdductHBondPotential() const;

	MembranePotential const & get_MembranePotential() const;

	Membrane_FAPotential const & get_Membrane_FAPotential() const; //pba

	PoissonBoltzmannPotential const & get_PoissonBoltzmannPotential() const;

	disulfides::FullatomDisulfidePotential &
	get_FullatomDisulfidePotential() const;

	disulfides::CentroidDisulfidePotential &
	get_CentroidDisulfidePotential() const;

	disulfides::DisulfideMatchingPotential &
	get_DisulfideMatchingPotential() const;

	//pack::dunbrack::SingleResidueRotamerLibraryCAP
	//get_NCAARotamerLibrary( chemical::ResidueType const & rsd_type );

	///
	methods::EnergyMethodOP
	energy_method( ScoreType const & t, methods::EnergyMethodOptions const & options ) const;

	///
	void
	add_etable( std::string const & name, etable::EtableOP etable );

	//XRW_B_T1
	/*
	///
	void
	add_coarse_etable( std::string const & name, coarse::CoarseEtableOP etable );
	*/
	//XRW_E_T1

	///pba
	void
	add_memb_etable( std::string const & name, etable::MembEtableOP etable );

	///pba
	etable::MembEtableCAP
	memb_etable( std::string const & table_id ) const;

	///
	etable::EtableCAP
	etable( std::string const & etable_id ) const;

	//XRW_B_T1
	/*
	///
	coarse::CoarseEtableCAP
	coarse_etable( std::string const & etable_id ) const;
	*/
	//XRW_E_T1

private:

	static ScoringManager * instance_;

	//private constructor
	ScoringManager();
	~ScoringManager();

	// WARNING -- if you add something here don't forget to initialize to 0 in the constructor
	mutable PairEPotential * pairE_potential_;
	//mutable RotamerLibrary * rotamer_Library_;
	mutable Ramachandran * rama_;
	mutable Ramachandran2B * rama2b_;
	mutable OmegaTether * omega_;
	mutable EnvPairPotential * env_pair_potential_;
	mutable SmoothEnvPairPotential * smooth_env_pair_potential_;
	mutable CenHBPotential * cen_hb_potential_;
	mutable SecondaryStructurePotential * secondary_structure_potential_;
	mutable std::map< std::string, AtomVDW * > atom_vdw_;
	mutable rna::RNA_AtomVDW * rna_atom_vdw_;
	mutable geometric_solvation::DatabaseOccSolEne * occ_hbond_sol_database_;
	mutable dna::DirectReadoutPotential * dna_dr_potential_;
	mutable MMLJLibrary * mm_lj_library_;
	mutable MMLJEnergyTable * mm_lj_energy_table_;
	mutable MMTorsionLibrary * mm_torsion_library_;
	mutable MMBondAngleLibrary * mm_bondangle_library_;
	mutable MMBondLengthLibrary * mm_bondlength_library_;
	//P_AA                     Paa_ppPotential_;
	mutable dna::DNA_BasePotential * DNA_base_potential_;
	mutable carbon_hbonds::CarbonHBondPotential * carbon_hbond_potential_;
	mutable rna::RNA_LowResolutionPotential * rna_low_resolution_potential_;
	mutable rna::RNA_TorsionPotential * rna_torsion_potential_;
	mutable P_AA * p_aa_;
	mutable WaterAdductHBondPotential * water_adduct_hbond_potential_;
	mutable GenBornPotential * gen_born_potential_;
	mutable disulfides::FullatomDisulfidePotential * fa_disulfide_potential_;
	mutable disulfides::CentroidDisulfidePotential * cen_disulfide_potential_;
	mutable disulfides::DisulfideMatchingPotential * disulfide_matching_potential_;
	mutable MembranePotential * membrane_potential_;
	mutable Membrane_FAPotential * membrane_fapotential_; //pba
	mutable PoissonBoltzmannPotential * PB_potential_;
	//ReferenceEnergyPotential referenceEnergyPotential_;
	mutable UnfoldedStatePotential * unf_state_;
	mutable nv::NVlookup * NV_lookup_table_;
	mutable orbitals::OrbitalsLookup * orbitals_lookup_table_;


  mutable interface::DDPlookup * DDP_lookup_table_;
	// data
	mutable std::map< std::string, etable::EtableOP > etables_;
	//XRW_B_T1
	//mutable std::map< std::string, coarse::CoarseEtableOP > coarse_etables_;
	//XRW_E_T1
  	mutable std::map< std::string, etable::MembEtableOP > memb_etables_; //pba


	// NCAA rot lib map
	/// mutable std::map< std::string, pack::dunbrack::SingleResidueRotamerLibraryCOP > ncaa_rotlibs_;

	utility::vector1< methods::EnergyMethodCreatorOP > method_creator_map_;

};
} // namespace core
} // namespace scoring


#endif
