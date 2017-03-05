// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoringManager.cc
/// @brief  Scoring manager class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/scoring/ScoringManager.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/unfolded_state.OptionKeys.gen.hh>
#include <basic/options/keys/orbitals.OptionKeys.gen.hh>

#include <core/scoring/carbon_hbonds/CarbonHBondPotential.hh>
#include <core/scoring/PairEPotential.hh>
#include <core/scoring/EnvPairPotential.hh>
#include <core/scoring/CenRotEnvPairPotential.hh>
#include <core/scoring/SmoothEnvPairPotential.hh>
#include <core/scoring/CenHBPotential.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/Membrane_FAPotential.hh> //pba
#include <core/scoring/ProQPotential.hh>
#include <core/scoring/SecondaryStructurePotential.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/Ramachandran2B.hh>
#include <core/scoring/OmegaTether.hh>
#include <core/scoring/GenBornPotential.hh>
#include <core/scoring/HydroxylTorsionPotential.hh>
#include <core/scoring/MultipoleElecPotential.hh>
#include <core/scoring/SASAPotential.hh>
#include <core/scoring/VdWTinkerPotential.hh>
#include <core/scoring/facts/FACTSPotential.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/rna/RNA_AtomVDW.hh>
#include <core/scoring/geometric_solvation/DatabaseOccSolEne.hh>
#include <core/scoring/carbohydrates/CHIEnergyFunction.hh>
#include <core/scoring/carbohydrates/OmegaPreferencesFunction.hh>
#include <core/scoring/dna/DNABFormPotential.hh>
#include <core/scoring/dna/DNATorsionPotential.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/rna/RNP_LowResPotential.hh>
#include <core/scoring/rna/RNP_LowResStackData.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/rna/RNA_SuitePotential.hh>
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.hh>
#include <core/scoring/rna/data/RNA_DMS_Potential.hh>
#include <core/scoring/rna/data/RNA_DMS_LowResolutionPotential.hh>
#include <core/scoring/loop_graph/evaluator/SixDTransRotPotential.hh>
#include <core/scoring/dna/DirectReadoutPotential.hh>
#include <core/scoring/P_AA.hh>
#include <core/scoring/P_AA_ss.hh>
#include <core/scoring/WaterAdductHBondPotential.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>
#include <core/scoring/disulfides/CentroidDisulfidePotential.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/scoring/UnfoldedStatePotential.hh>
#include <core/scoring/PoissonBoltzmannPotential.hh>
#include <core/scoring/dna/DNA_EnvPairPotential.hh>
#include <core/scoring/dna/DNA_DihedralPotential.hh>
#include <core/scoring/SplitUnfoldedTwoBodyPotential.hh>
#include <core/scoring/elec/util.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/memb_etable/MembEtable.hh>

#include <core/scoring/mm/MMTorsionLibrary.hh>
#include <core/scoring/mm/MMLJLibrary.hh>
#include <core/scoring/mm/MMLJEnergyTable.hh>
#include <core/scoring/mm/MMBondAngleLibrary.hh>
#include <core/scoring/mm/MMBondLengthLibrary.hh>

#include <core/scoring/nv/NVlookup.hh>
#include <core/scoring/orbitals/OrbitalsLookup.hh>
#include <core/scoring/interface_/DDPlookup.hh>

#include <core/scoring/types.hh>
#include <core/scoring/ScoreType.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/mainchain_potential/MainchainScoreTable.fwd.hh>
#include <core/chemical/mainchain_potential/util.hh>

#include <basic/database/open.hh>


// Utility headers
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

#include <core/scoring/methods/EnergyMethod.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>

// Basic headers
#include <basic/Tracer.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.ScoringManager" );

namespace core {
namespace scoring {

ScoringManager::~ScoringManager() = default;

///////////////////////////////////////////////////////////////////////////////
ScoringManager::ScoringManager() :
	vdw_tinker_potential_( /* 0 */ ),
	pairE_potential_( /* 0 */ ),
	rama_( /* 0 */ ),
	rama2b_( /* 0 */ ),
	rama_pp_( /* 0 */ ),
	omega_( /* 0 */ ),
	env_pair_potential_( /* 0 */ ),
	smooth_env_pair_potential_( /* 0 */ ),
	cen_rot_pair_potential_( /* 0 */ ),
	cen_hb_potential_( /* 0 */ ),
	secondary_structure_potential_( /* 0 */ ),
	atom_vdw_(),
	rna_atom_vdw_( /* 0 */ ),
	occ_hbond_sol_database_( /* 0 */ ),
	dna_dr_potential_( /* 0 */ ),
	mm_lj_library_( /* 0 */ ),
	mm_lj_energy_table_( /* 0 */ ),
	mm_torsion_library_( /* 0 */ ),
	mm_bondangle_library_( /* 0 */ ),
	mm_bondlength_library_( /* 0 */ ),
	dnabform_( /* 0 */ ),
	dna_torsion_potential_( /* 0 */ ),
	DNA_base_potential_( /* 0 */ ),
	carbon_hbond_potential_( /* 0 */ ),
	rna_low_resolution_potential_( /* 0 */ ),
	rnp_low_res_potential_( /* 0 */ ),
	rnp_low_res_stack_data_( /* 0 */ ),
	// rna_torsion_potential_( /* 0 */ ),
	rna_chemical_shift_potential_( /* 0 */ ),
	rna_dms_potential_( /* 0 */ ),
	rna_dms_low_resolution_potential_( /* 0 */ ),
	p_aa_( /* 0 */ ),
	p_aa_ss_( /* 0 */ ),
	water_adduct_hbond_potential_( /* 0 */ ),
	gen_born_potential_( /* 0 */ ),
	hxl_tors_potential_( /* 0 */ ),
	fa_disulfide_potential_( /* 0 */ ),
	cen_disulfide_potential_( /* 0 */ ),
	disulfide_matching_potential_( /* 0 */ ),
	membrane_potential_( /* 0 */ ),
	membrane_fapotential_( /* 0 */ ), //pba
	ProQ_potential_(/* 0 */),
	PB_potential_(/* 0 */),
	sutbp_( /* 0 */),
	unf_state_( /* 0 */ ),
	CHI_energy_function_( /* 0 */ ),
	NV_lookup_table_(/* 0 */),
	orbitals_lookup_table_( /* 0 */ ),
	DDP_lookup_table_(/* 0 */),
	etables_by_string_(),
	etables_by_options_(),
	memb_etables_(),
	cp_rep_map_byname_(),
	aa_composition_setup_helpers_(),
	rama_prepro_mainchain_potentials_(),
	rama_prepro_mainchain_potentials_beforeproline_(),
	method_creator_map_( n_score_types, nullptr )
{}

///////////////////////////////////////////////////////////////////////////////
VdWTinkerPotential const &
ScoringManager::get_VdWTinkerPotential() const
{
	if ( vdw_tinker_potential_ == nullptr ) {
		vdw_tinker_potential_ = VdWTinkerPotentialOP( new VdWTinkerPotential() );
	}
	return *vdw_tinker_potential_;
}

///////////////////////////////////////////////////////////////////////////////
PairEPotential const &
ScoringManager::get_PairEPotential() const
{
	if ( pairE_potential_ == nullptr ) {
		pairE_potential_ = PairEPotentialOP( new PairEPotential() );
	}
	return *pairE_potential_;
}

/// @details The ScoringManager acts as an EnergyMethodFactory.  All EnergyMethods must
/// create a helper class, an EnergyMethodCreator class, that will respond to a call to
/// its create_energy_method by returning a new instance of that EnergyMethod its helping.
/// This Creator class must also register itself with the ScoringManager at load time and
/// hand an instance of itself to the singleton ScoringManager instance.
void
ScoringManager::factory_register( methods::EnergyMethodCreatorOP creator )
{
	ScoreTypes sts = creator->score_types_for_method();
	for ( Size ii = 1; ii <= sts.size(); ++ii ) {
		///std::cout << "Registering " << (int) sts[ ii ] << " to " << creator() << std::endl;
		// make sure no two EnergyMethodCreators lay claim to the same ScoreType
		if ( method_creator_map_[ sts[ ii ] ] != nullptr ) {
			utility_exit_with_message( "Cannot register a term to two different EnergyMethodCreators. Term " + utility::to_string( sts[ ii ] ) + " has already been registered!" );
		}
		method_creator_map_[ sts[ ii ] ] = creator;
	}
}


///////////////////////////////////////////////////////////////////////////////
dna::DNA_BasePotential const &
ScoringManager::get_DNA_BasePotential() const
{
	if ( DNA_base_potential_ == nullptr ) {
		DNA_base_potential_ = dna::DNA_BasePotentialOP( new dna::DNA_BasePotential() );
	}
	return *DNA_base_potential_;
}

///////////////////////////////////////////////////////////////////////////////
dna::DNA_DihedralPotential const &
ScoringManager::get_DNA_DihedralPotential() const
{
	if ( dna_dihedral_potential_ == nullptr ) {
		dna_dihedral_potential_ = new dna::DNA_DihedralPotential();
	}
	return *dna_dihedral_potential_;
}

///////////////////////////////////////////////////////////////////////////////
dna::DNA_EnvPairPotential const &
ScoringManager::get_DNA_EnvPairPotential() const
{
	if ( dna_env_pair_potential_ == nullptr ) {
		dna_env_pair_potential_ = new dna::DNA_EnvPairPotential();
	}
	return *dna_env_pair_potential_;
}

///////////////////////////////////////////////////////////////////////////////
EnvPairPotential const &
ScoringManager::get_EnvPairPotential() const
{
	if ( env_pair_potential_ == nullptr ) {
		env_pair_potential_ = EnvPairPotentialOP( new EnvPairPotential() );
	}
	return *env_pair_potential_;
}

///////////////////////////////////////////////////////////////////////////////
SmoothEnvPairPotential const &
ScoringManager::get_SmoothEnvPairPotential() const
{
	if ( smooth_env_pair_potential_ == nullptr ) {
		smooth_env_pair_potential_ = SmoothEnvPairPotentialOP( new SmoothEnvPairPotential() );
	}
	return *smooth_env_pair_potential_;
}

///////////////////////////////////////////////////////////////////////////////
CenRotEnvPairPotential const &
ScoringManager::get_CenRotEnvPairPotential() const
{
	if ( cen_rot_pair_potential_ == nullptr ) {
		cen_rot_pair_potential_ = CenRotEnvPairPotentialOP( new CenRotEnvPairPotential() );
	}
	return *cen_rot_pair_potential_;
}

///////////////////////////////////////////////////////////////////////////////
CenHBPotential const &
ScoringManager::get_CenHBPotential() const
{
	if ( cen_hb_potential_ == nullptr ) {
		cen_hb_potential_ = CenHBPotentialOP( new CenHBPotential() );
	}
	return *cen_hb_potential_;
}

///////////////////////////////////////////////////////////////////////////////
MembranePotential const &
ScoringManager::get_MembranePotential() const
{
	if ( membrane_potential_ == nullptr ) {
		membrane_potential_ = MembranePotentialOP( new MembranePotential() );
	}
	return *membrane_potential_;
}

///////////////////////////////////////////////////////////////////////////////
membrane::MembraneData const &
ScoringManager::get_MembraneData() const
{
	if ( mp_base_potential_ == nullptr ) {
		mp_base_potential_ = membrane::MembraneDataOP( new membrane::MembraneData() );
	}
	return *mp_base_potential_;
}

///////////////////////////////////////////////////////////////////////////////
Membrane_FAPotential const &
ScoringManager::get_Membrane_FAPotential() const //pba
{
	if ( membrane_fapotential_ == nullptr ) {
		membrane_fapotential_ = Membrane_FAPotentialOP( new Membrane_FAPotential() );
	}
	return *membrane_fapotential_;
}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
ProQPotential const &
ScoringManager::get_ProQPotential() const
{
	if ( ProQ_potential_ == nullptr ) {
		ProQ_potential_ = ProQPotentialOP( new ProQPotential() );
	}
	return *ProQ_potential_;
}

///////////////////////////////////////////////////////////////////////////////


SecondaryStructurePotential const &
ScoringManager::get_SecondaryStructurePotential() const
{
	if ( secondary_structure_potential_ == nullptr ) {
		secondary_structure_potential_ = SecondaryStructurePotentialOP( new SecondaryStructurePotential() );
	}
	return *secondary_structure_potential_;
}

///////////////////////////////////////////////////////////////////////////////
GenBornPotential const &
ScoringManager::get_GenBornPotential() const
{
	if ( gen_born_potential_ == nullptr ) {
		gen_born_potential_ = GenBornPotentialOP( new GenBornPotential() );
	}
	return *gen_born_potential_;
}

HydroxylTorsionPotential const &
ScoringManager::get_HydroxylTorsionPotential() const
{
	if ( hxl_tors_potential_ == nullptr ) {
		hxl_tors_potential_ = HydroxylTorsionPotentialOP( new HydroxylTorsionPotential );
	}
	return *hxl_tors_potential_;
}


///////////////////////////////////////////////////////////////////////////////
MultipoleElecPotential const &
ScoringManager::get_MultipoleElecPotential( methods::EnergyMethodOptions const & options ) const
{
	if ( multipole_elec_potential_ == nullptr ) {
		multipole_elec_potential_ = MultipoleElecPotentialOP( new MultipoleElecPotential() );
		multipole_elec_potential_->use_polarization = options.use_polarization();
		multipole_elec_potential_->use_gen_kirkwood = options.use_gen_kirkwood();
		multipole_elec_potential_->Ep = options.protein_dielectric();
		multipole_elec_potential_->Ew = options.water_dielectric();
	}
	return *multipole_elec_potential_;
}

///////////////////////////////////////////////////////////////////////////////
SASAPotential const &
ScoringManager::get_SASAPotential() const
{
	if ( sasa_potential_ == nullptr ) {
		sasa_potential_ = SASAPotentialOP( new SASAPotential() );
	}
	return *sasa_potential_;
}

///////////////////////////////////////////////////////////////////////////////
FACTSPotential const &
ScoringManager::get_FACTSPotential() const
{
	if ( facts_potential_ == nullptr ) {
		facts_potential_ = FACTSPotentialOP( new FACTSPotential() );
	}
	return *facts_potential_;
}

///////////////////////////////////////////////////////////////////////////////
PoissonBoltzmannPotential const &
ScoringManager::get_PoissonBoltzmannPotential() const
{
	if ( PB_potential_ == nullptr ) {
		PB_potential_ = PoissonBoltzmannPotentialOP( new PoissonBoltzmannPotential );
	}
	return *PB_potential_;
}

///////////////////////////////////////////////////////////////////////////////
AtomVDW const &
ScoringManager::get_AtomVDW( std::string const & atom_type_set_name ) const
{
	if ( atom_vdw_.count( atom_type_set_name ) == 0 ) {
		atom_vdw_[ atom_type_set_name ] = AtomVDWOP( new AtomVDW( atom_type_set_name ) );
	}
	return * ( atom_vdw_[ atom_type_set_name ] );
}

///////////////////////////////////////////////////////////////////////////////
carbon_hbonds::CarbonHBondPotential const &
ScoringManager::get_CarbonHBondPotential() const
{
	if ( carbon_hbond_potential_ == nullptr ) {
		carbon_hbond_potential_ = carbon_hbonds::CarbonHBondPotentialOP( new carbon_hbonds::CarbonHBondPotential() );
	}
	return *carbon_hbond_potential_;
}

///////////////////////////////////////////////////////////////////////////////
rna::RNA_AtomVDW const &
ScoringManager::get_RNA_AtomVDW() const
{
	if ( rna_atom_vdw_ == nullptr ) {
		rna_atom_vdw_ = rna::RNA_AtomVDWOP( new rna::RNA_AtomVDW() );
	}
	return *rna_atom_vdw_;
}

///////////////////////////////////////////////////////////////////////////////
geometric_solvation::DatabaseOccSolEne const &
ScoringManager::get_DatabaseOccSolEne( std::string const & atom_type_set_name, Real const & min_occ_energy ) const
{
	if ( occ_hbond_sol_database_ == nullptr ) {
		occ_hbond_sol_database_ = geometric_solvation::DatabaseOccSolEneOP( new geometric_solvation::DatabaseOccSolEne( atom_type_set_name, min_occ_energy ) );
	}
	return *occ_hbond_sol_database_;
}

///////////////////////////////////////////////////////////////////////////////
rna::RNA_LowResolutionPotential const &
ScoringManager::get_RNA_LowResolutionPotential() const
{
	if ( rna_low_resolution_potential_ == nullptr ) {
		rna_low_resolution_potential_ = rna::RNA_LowResolutionPotentialOP( new rna::RNA_LowResolutionPotential() );
	}
	return *rna_low_resolution_potential_;
}

///////////////////////////////////////////////////////////////////////////////
rna::RNP_LowResPotential const &
ScoringManager::get_RNP_LowResPotential() const
{
	if ( rnp_low_res_potential_ == 0 ) {
		rnp_low_res_potential_ = rna::RNP_LowResPotentialOP( new rna::RNP_LowResPotential() );
	}
	return *rnp_low_res_potential_;
}
///////////////////////////////////////////////////////////////////////////////
rna::RNP_LowResStackData const &
ScoringManager::get_RNP_LowResStackData() const
{
	if ( rnp_low_res_stack_data_ == 0 ) {
		rnp_low_res_stack_data_ = rna::RNP_LowResStackDataOP( new rna::RNP_LowResStackData() );
	}
	return *rnp_low_res_stack_data_;
}
///////////////////////////////////////////////////////////////////////////////
// rna::RNA_TorsionPotential const &
// ScoringManager::get_RNA_TorsionPotential() const
// {
//  if (rna_torsion_potential_ == 0 )
//  {
//   rna_torsion_potential_ = rna::RNA_TorsionPotentialOP( new rna::RNA_TorsionPotential() );
//  }
//  return *rna_torsion_potential_;
// }

///////////////////////////////////////////////////////////////////////////////
// rna::RNA_SuitePotential const &
// ScoringManager::get_RNA_SuitePotential( bool const calculate_suiteness_bonus /* = false */ ) const
// {
//  if ( calculate_suiteness_bonus ){
//   if (rna_suite_potential_for_suiteness_bonus_ == 0 )
//    {
//     rna_suite_potential_for_suiteness_bonus_ = rna::RNA_SuitePotentialOP( new rna::RNA_SuitePotential( true ) );
//    }
//   return *rna_suite_potential_for_suiteness_bonus_;
//  }

//  if (rna_suite_potential_ == 0 )
//   {
//    rna_suite_potential_ = rna::RNA_SuitePotentialOP( new rna::RNA_SuitePotential( false ) );
//   }
//  return *rna_suite_potential_;
// }

///////////////////////////////////////////////////////////////////////////////
rna::chemical_shift::RNA_ChemicalShiftPotential const &
ScoringManager::get_RNA_ChemicalShiftPotential() const
{
	if ( rna_chemical_shift_potential_ == nullptr ) {
		rna_chemical_shift_potential_= rna::chemical_shift::RNA_ChemicalShiftPotentialOP( new rna::chemical_shift::RNA_ChemicalShiftPotential() );
	}
	return *rna_chemical_shift_potential_;
}

///////////////////////////////////////////////////////////////////////////////
rna::data::RNA_DMS_Potential &
ScoringManager::get_RNA_DMS_Potential() const
{
	if ( rna_dms_potential_ == nullptr ) {
		rna_dms_potential_ = rna::data::RNA_DMS_PotentialOP( new rna::data::RNA_DMS_Potential );
	}
	return *rna_dms_potential_;
}

///////////////////////////////////////////////////////////////////////////////
rna::data::RNA_DMS_LowResolutionPotential &
ScoringManager::get_RNA_DMS_LowResolutionPotential() const
{
	if ( rna_dms_low_resolution_potential_ == nullptr ) {
		rna_dms_low_resolution_potential_ = rna::data::RNA_DMS_LowResolutionPotentialOP( new rna::data::RNA_DMS_LowResolutionPotential );
	}
	return *rna_dms_low_resolution_potential_;
}

///////////////////////////////////////////////////////////////////////////////
loop_graph::evaluator::SixDTransRotPotentialCOP
ScoringManager::get_LoopCloseSixDPotential( std::string const & database_file ) const
{
	if ( loop_close_six_d_potential_[ database_file ] == nullptr ) {
		if ( utility::file::file_exists(  database_file ) )  {
			TR << "Reading in: " << database_file << std::endl;
			loop_close_six_d_potential_[ database_file ] = loop_graph::evaluator::SixDTransRotPotentialCOP( new loop_graph::evaluator::SixDTransRotPotential( database_file ) );
		} else {
			loop_close_six_d_potential_[ database_file ] = 0; // save information that database file does not exist.
		}
	}
	return loop_close_six_d_potential_[ database_file ];
}

///////////////////////////////////////////////////////////////////////////////
dna::DirectReadoutPotential const &
ScoringManager::get_DirectReadoutPotential() const
{
	if ( dna_dr_potential_ == nullptr ) {
		dna_dr_potential_ = dna::DirectReadoutPotentialOP( new dna::DirectReadoutPotential() );
	}
	return *dna_dr_potential_;
}

P_AA const &
ScoringManager::get_P_AA() const
{
	if ( p_aa_ == nullptr ) {
		p_aa_ = P_AAOP( new P_AA );
	}
	return *p_aa_;
}

P_AA_ss const &
ScoringManager::get_P_AA_ss() const
{
	if ( p_aa_ss_ == nullptr ) {
		p_aa_ss_ = P_AA_ssOP( new P_AA_ss );
	}
	return *p_aa_ss_;
}

WaterAdductHBondPotential const &
ScoringManager::get_WaterAdductHBondPotential() const
{
	if ( water_adduct_hbond_potential_ == nullptr ) {
		water_adduct_hbond_potential_ = WaterAdductHBondPotentialOP( new WaterAdductHBondPotential );
	}
	return *water_adduct_hbond_potential_;
}


///////////////////////////////////////////////////////////////////////////////
RamachandranCOP
ScoringManager::get_Ramachandran_ptr() const
{
	if ( rama_ == nullptr ) {
		rama_ = RamachandranOP( new Ramachandran );
	}
	return rama_;
}

/// @brief Get a non-const owning pointer to the Ramachandran object.
/// @details Warning!  Probably not threadsafe!  Needed, though, for
/// lazy loading.
RamachandranOP
ScoringManager::get_Ramachandran_ptr_nonconst()
{
	if ( rama_ == nullptr ) {
		rama_ = RamachandranOP( new Ramachandran );
	}
	return rama_;
}

Ramachandran const &
ScoringManager::get_Ramachandran() const
{
	return *get_Ramachandran_ptr();
}

/// @brief Get a non-const reference to the Ramachandran object.
/// @details Warning!  Probably not threadsafe!  Needed, though, for
/// lazy loading.
Ramachandran &
ScoringManager::get_Ramachandran_nonconst() {
	return *get_Ramachandran_ptr_nonconst();
}

///////////////////////////////////////////////////////////////////////////////
Ramachandran2BCOP
ScoringManager::get_Ramachandran2B_ptr() const
{
	if ( rama2b_ == nullptr ) {
		rama2b_ = Ramachandran2BOP( new Ramachandran2B );
	}
	return rama2b_;
}

Ramachandran2B const &
ScoringManager::get_Ramachandran2B() const
{
	return *get_Ramachandran2B_ptr();
}

///////////////////////////////////////////////////////////////////////////////
RamaPrePro const &
ScoringManager::get_RamaPrePro() const
{
	if ( rama_pp_ == nullptr ) {
		rama_pp_ = RamaPreProOP( new RamaPrePro );
	}
	return *rama_pp_;
}

///////////////////////////////////////////////////////////////////////////////
OmegaTether const &
ScoringManager::get_OmegaTether() const
{
	if ( omega_ == nullptr ) {
		omega_ = OmegaTetherOP( new OmegaTether() );
	}
	return *omega_;
}

///////////////////////////////////////////////////////////////////////////////
dna::DNABFormPotential const &
ScoringManager::get_DNABFormPotential() const
{
	if ( dnabform_ == nullptr ) {
		dnabform_ =  dna::DNABFormPotentialOP( new dna::DNABFormPotential );
	}
	return *dnabform_;
}

///////////////////////////////////////////////////////////////////////////////

dna::DNATorsionPotential const &
ScoringManager::get_DNATorsionPotential() const
{
	if ( dna_torsion_potential_ == nullptr ) {
		dna_torsion_potential_ = dna::DNATorsionPotentialOP( new dna::DNATorsionPotential() );
	}
	return *dna_torsion_potential_;
}

///////////////////////////////////////////////////////////////////////////////
core::scoring::mm::MMTorsionLibrary const &
ScoringManager::get_MMTorsionLibrary() const
{
	if ( mm_torsion_library_ == nullptr ) {
		mm_torsion_library_ = mm::MMTorsionLibraryOP( new mm::MMTorsionLibrary
			( basic::database::full_name( "chemical/mm_atom_type_sets/fa_standard/mm_torsion_params.txt" ),
			chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD ) ) );
	}
	return *mm_torsion_library_;
}
///////////////////////////////////////////////////////////////////////////////

core::scoring::mm::MMLJLibrary const &
ScoringManager::get_MMLJLibrary() const
{
	if ( mm_lj_library_ == nullptr ) {
		mm_lj_library_ = mm::MMLJLibraryOP( new mm::MMLJLibrary
			( chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD ) ) );
	}
	return *mm_lj_library_;
}

///////////////////////////////////////////////////////////////////////////////

core::scoring::mm::MMLJEnergyTable const &
ScoringManager::get_MMLJEnergyTable () const
{
	if ( mm_lj_energy_table_ == nullptr ) {
		mm_lj_energy_table_ = mm::MMLJEnergyTableOP( new mm::MMLJEnergyTable() );
	}
	return *mm_lj_energy_table_;
}

///////////////////////////////////////////////////////////////////////////////

disulfides::FullatomDisulfidePotential &
ScoringManager::get_FullatomDisulfidePotential() const
{
	if ( fa_disulfide_potential_ == nullptr ) {
		fa_disulfide_potential_ = disulfides::FullatomDisulfidePotentialOP( new disulfides::FullatomDisulfidePotential );
	}
	return *fa_disulfide_potential_;
}

disulfides::CentroidDisulfidePotential &
ScoringManager::get_CentroidDisulfidePotential() const
{
	if ( cen_disulfide_potential_ == nullptr ) {
		cen_disulfide_potential_ = disulfides::CentroidDisulfidePotentialOP( new disulfides::CentroidDisulfidePotential );
	}
	return *cen_disulfide_potential_;
}

disulfides::DisulfideMatchingPotential &
ScoringManager::get_DisulfideMatchingPotential() const
{
	if ( disulfide_matching_potential_ == nullptr ) {
		disulfide_matching_potential_ = disulfides::DisulfideMatchingPotentialOP( new disulfides::DisulfideMatchingPotential );
	}
	return *disulfide_matching_potential_;
}

///////////////////////////////////////////////////////////////////////////////


nv::NVlookup const &
ScoringManager::get_NVLookupTable() const
{
	if ( NV_lookup_table_ == nullptr ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		NV_lookup_table_ = nv::NVlookupOP( new nv::NVlookup(basic::database::full_name(option[score::NV_table]())) );
	}
	return *NV_lookup_table_;

}

///////////////////////////////////////////////////////////////////////////////
orbitals::OrbitalsLookup const &
ScoringManager::get_OrbitalsLookupTable() const
{
	if ( orbitals_lookup_table_ == nullptr ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		utility::vector1<std::string> DHO_energies;
		utility::vector1< std::string > AOH_energies;
		utility::vector1< std::string > AOO_orb_orb_energies;
		utility::vector1< std::string > DOO_orb_orb_energies;
		utility::vector1< std::string > ACO_energies;

		DHO_energies.push_back("scoring/score_functions/orbitals/BiCubic_DHO_Hpol_scOrbH.txt");//DHO_energies[1]
		DHO_energies.push_back("scoring/score_functions/orbitals/BiCubic_DHO_Hpol_bbOrbH.txt");//DHO_energies[2]
		DHO_energies.push_back("scoring/score_functions/orbitals/BiCubic_DHO_Haro_scOrbH.txt");//DHO_energies[3]

		AOH_energies.push_back("scoring/score_functions/orbitals/BiCubic_AOH_Hpol_scOrbH.txt");//AOH_energies[1]
		AOH_energies.push_back("scoring/score_functions/orbitals/BiCubic_AOH_Hpol_bbOrbH.txt");//AOH_energies[2]
		AOH_energies.push_back("scoring/score_functions/orbitals/BiCubic_AOH_Haro_scOrbH.txt");//AOH_energies[3]

		AOO_orb_orb_energies.push_back("scoring/score_functions/orbitals/BiCubic_AOD_OrbOrb.txt");
		DOO_orb_orb_energies.push_back("scoring/score_functions/orbitals/BiCubic_DOA_OrbOrb.txt");


		ACO_energies.push_back("scoring/score_functions/orbitals/BiCubic_ACO.txt");

		orbitals_lookup_table_ = orbitals::OrbitalsLookupOP( new orbitals::OrbitalsLookup(
			DHO_energies, AOH_energies, AOO_orb_orb_energies, DOO_orb_orb_energies,ACO_energies ) );

	}
	return *orbitals_lookup_table_;
}

///////////////////////////////////////////////////////////////////////////////


interface_::DDPlookup const &
ScoringManager::get_DDPLookupTable() const
{
	if ( DDP_lookup_table_ == nullptr ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		DDP_lookup_table_ = interface_::DDPlookupOP( new interface_::DDPlookup("scoring/score_functions/DDPscore/interface_ddp_score.txt") );
	}
	return *DDP_lookup_table_;
}

///////////////////////////////////////////////////////////////////////////////
core::scoring::mm::MMBondAngleLibrary const &
ScoringManager::get_MMBondAngleLibrary() const
{
	if ( mm_bondangle_library_ == nullptr ) {
		mm_bondangle_library_ = mm::MMBondAngleLibraryOP( new mm::MMBondAngleLibrary
			( basic::database::full_name( "chemical/mm_atom_type_sets/fa_standard/par_all27_prot_na.prm" ),
			chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD ) ) );
	}
	return *mm_bondangle_library_;
}

///////////////////////////////////////////////////////////////////////////////
core::scoring::mm::MMBondLengthLibrary const &
ScoringManager::get_MMBondLengthLibrary() const
{
	if ( mm_bondlength_library_ == nullptr ) {
		mm_bondlength_library_ = mm::MMBondLengthLibraryOP( new mm::MMBondLengthLibrary
			( basic::database::full_name( "chemical/mm_atom_type_sets/fa_standard/par_all27_prot_na.prm" ),
			chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD ) ) );
	}
	return *mm_bondlength_library_;
}

///////////////////////////////////////////////////////////////////////////////
UnfoldedStatePotential const &
ScoringManager::get_UnfoldedStatePotential( std::string const & type ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( unf_state_ == nullptr ) {
		if ( type == UNFOLDED_SPLIT_USER_DEFINED || option[ unfolded_state::unfolded_energies_file ].user() ) {
			TR << "Creating unfolded state potential using file: " <<  option[ unfolded_state::unfolded_energies_file ].value() << std::endl;
			unf_state_ = UnfoldedStatePotentialOP( new UnfoldedStatePotential( option[ unfolded_state::unfolded_energies_file ].value() ) );
		} else if ( type == UNFOLDED_SCORE12 ) {
			unf_state_ = UnfoldedStatePotentialOP( new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/unfolded/unfolded_state_residue_energies_score12" ) ) );
		} else if ( type == UNFOLDED_MM_STD ) {
			unf_state_ = UnfoldedStatePotentialOP( new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/unfolded/unfolded_state_residue_energies_mm_std" ) ) );
		} else if ( type == UNFOLDED_RNA ) {
			unf_state_ = UnfoldedStatePotentialOP( new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/unfolded/unfolded_state_residue_energies_rna" ) ) );  // This will later get more elaborated
		} else if ( type == UNFOLDED_SPLIT_TALARIS2013 ) {
			unf_state_ = UnfoldedStatePotentialOP( new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/split_unfolded/split_unfolded_one_body_talaris2013" ) ) ); //for the split unfolded energy one body term.
		} else if ( type == UNFOLDED_SPLIT_MM_STD ) {
			unf_state_ = UnfoldedStatePotentialOP( new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/split_unfolded/split_unfolded_one_body_mm_std" ) ) ); //for the split unfolded energy one body term.
		} else {
			utility_exit_with_message("unrecognized unfolded type: "+type );
		}
	}
	return *unf_state_;
}

///////////////////////////////////////////////////////////////////////////////
SplitUnfoldedTwoBodyPotential const &
ScoringManager::get_SplitUnfoldedTwoBodyPotential(std::string const & label_type,std::string const & value_type, std::string const & score_func_type) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( sutbp_ == nullptr ) {
		std::string database_path = "";
		std::string atom_label_type = "";

		if ( label_type == SPLIT_UNFOLDED_USER_DEFINED || value_type == SPLIT_UNFOLDED_USER_DEFINED || option[ unfolded_state::split_unfolded_energies_file ].user() ) {
			sutbp_ = SplitUnfoldedTwoBodyPotentialOP( new SplitUnfoldedTwoBodyPotential( option[ unfolded_state::split_unfolded_energies_file ].value(), option[ unfolded_state::split_unfolded_energies_atom_type ].value() ) );
			TR << "Creating split unfolded state potential using file: " <<  option[ unfolded_state::split_unfolded_energies_file ].value()
				<< " with an atom type of " << option[ unfolded_state::split_unfolded_energies_atom_type ].value() << std::endl;
			return *sutbp_;
		}

		if ( label_type == SPLIT_UNFOLDED_ELE ) {
			atom_label_type = "elemental";
			database_path = "scoring/score_functions/split_unfolded/ele_database_";
		} else if ( label_type==SPLIT_UNFOLDED_PDB ) {
			atom_label_type = "pdb";
			database_path = "scoring/score_functions/split_unfolded/pdb_database_";
		} else if ( label_type == SPLIT_UNFOLDED_ROSETTA ) {
			atom_label_type = "rosetta";
			database_path = "scoring/score_functions/split_unfolded/rosetta_database_";
		} else if ( label_type == SPLIT_UNFOLDED_MM ) {
			atom_label_type = "mm";
			database_path = "scoring/score_functions/split_unfolded/mm_database_";
		} else if ( label_type == SPLIT_UNFOLDED_UNIQUE ) {
			atom_label_type = "unique";
			database_path = "scoring/score_functions/split_unfolded/unique_database_";
		} else {
			database_path="scoring/score_functions/split_unfolded/unique_database_"; //default to unique types, since canonical design is the most common use-case(probably).
		}

		//Also need to account for the energy function being used so as to use the correct weights for the two body internal energy values, which are specified in the database files. This is specified here via the UNFOLDED_ENERGIES_TYPE weight file tag.
		//talaris2013 and mm_std are currently supported, but adding more score functions just requires replacing the weights at the top of the database files with the weights for each energy type in the new score function and defining a new option here(assuming all energies in the new score function are accounted for in the existing file).
		if ( score_func_type == UNFOLDED_SPLIT_TALARIS2013 ) {
			database_path += "talaris2013_";
		} else if ( score_func_type == UNFOLDED_SPLIT_MM_STD ) {
			database_path += "mm_std_";
		} else {
			database_path += "talaris2013_"; //default to talaris if we have no idea which to use
		}

		if ( value_type == SPLIT_UNFOLDED_MEAN ) {
			database_path += "mean";
		} else if ( value_type == SPLIT_UNFOLDED_MEDIAN ) {
			database_path += "median";
		} else if ( value_type == SPLIT_UNFOLDED_MODE ) {
			database_path += "mode";
		} else if ( value_type == SPLIT_UNFOLDED_BOLTZ ) {
			database_path += "boltz";
		} else {
			database_path += "median"; //median seems to work best, so default to that.
		}

		if ( atom_label_type != "" ) {
			TR << "Creating split unfolded state potential using file: " << basic::database::full_name( database_path ) << std::endl;
			sutbp_ = SplitUnfoldedTwoBodyPotentialOP( new SplitUnfoldedTwoBodyPotential( basic::database::full_name( database_path ), atom_label_type ) );
		} else {
			TR << "Creating split unfolded state potential using file: " << basic::database::full_name( database_path ) << std::endl;
			sutbp_ = SplitUnfoldedTwoBodyPotentialOP( new SplitUnfoldedTwoBodyPotential( basic::database::full_name( database_path ) ) );
		}
	}
	return *sutbp_;
}


///////////////////////////////////////////////////////////////////////////////
carbohydrates::CHIEnergyFunction const &
ScoringManager::get_CHIEnergyFunction( bool setup_for_sampling /* false */, Real step_size /* 0.1 */ ) const
{
	if ( CHI_energy_function_ == nullptr ) {
		TR << "Creating CHI Energy Function." << std::endl;
		CHI_energy_function_ = carbohydrates::CHIEnergyFunctionOP( new carbohydrates::CHIEnergyFunction );
	}

	if ( setup_for_sampling && ( ! CHI_energy_function_->sampling_data_setup() ) ) {
		TR << "should be setting up for sampling..." << std::endl;
		CHI_energy_function_->setup_for_sampling( step_size );
	}

	return *CHI_energy_function_;
}

///////////////////////////////////////////////////////////////////////////////
carbohydrates::OmegaPreferencesFunction const &
ScoringManager::get_OmegaPreferencesFunction(  bool setup_for_sampling /* false */, Real step_size /* 0.1 */ ) const
{
	if ( carbohydrate_omega_preferences_function_ == nullptr ) {
		TR << "Creating carbohydrate omega preferences function." << std::endl;
		carbohydrate_omega_preferences_function_ =
			carbohydrates::OmegaPreferencesFunctionOP( new carbohydrates::OmegaPreferencesFunction );

	}

	if ( setup_for_sampling && ( ! carbohydrate_omega_preferences_function_->sampling_data_setup() ) ) {
		TR << "should be setting up for sampling..." << std::endl;
		carbohydrate_omega_preferences_function_->setup_for_sampling( step_size );
	}


	return *carbohydrate_omega_preferences_function_;
}


//////////////////////////////////////////////////////////////////////////////
/// @details Make etable for extra partial softies, pilot app r_play_with_etables does not really work anymore
/// the etables it added will be somehow cleared or overwriten, so I do it here now
///
/// table_id: i.e. FA_STANDARD_SOFT40,
/// the number in the end is a percentage, lj_radius are given for every 5% softie
/// from 5% to 95% in database/chemical/atom_type_sets/fa_standard/extras/extra_soft_rep_params.txt
/// FA_STANDARD_SOFT50 would be halfway between normal softrep and normal hardrep.
etable::EtableOP
ScoringManager::make_partially_soft_etable( std::string const & table_id, etable::EtableOptions etable_options ) const
{
	using namespace etable;
	using namespace utility;

	debug_assert( utility::startswith( table_id, "FA_STANDARD_SOFT" ) );
	std::string table_name( string_split( table_id, '_').back() );
	std::string table_value( trim( table_id, "FA_STANDARD_SOFT" ) );
	if ( string2float( table_value ) != -1 ) {
		core::Real weight( string2Real( table_value ) * 0.01 );
		if ( weight < 1 && weight > 0 && ( string2int( table_value ) % 5 ==0 ) ) {
			etable_options.lj_switch_dis2sigma = 0.6*(1.0-weight) + weight*0.91;
			return etable::EtableOP(new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
				etable_options, table_name ));
		} else {
			std::string msg = "lj_radius is not given for partially soft " + trim(table_id, "FA_STANDARD_SOFT") + "%";
			utility_exit_with_message( msg );
		}
	} else {
		std::string msg = "input " + table_id + " format error, expecting a integer between 5 and 95 after FA_STANDARD_SOFT";
		utility_exit_with_message( msg );
	}
}


///////////////////////////////////////////////////p////////////////////////////
void
ScoringManager::add_memb_etable( std::string const & name, etable::MembEtableOP etable ) //pba
{
	debug_assert( memb_etables_.count(name) == 0 );
	memb_etables_[ name ] = etable;
}

///////////////////////////////////////////////////////////////////////////////
etable::MembEtableCAP
ScoringManager::memb_etable( std::string const & table_id ) const //pba
{
	using namespace etable;

	if ( memb_etables_.find( table_id ) == memb_etables_.end() ) {
		// try to build if possible
		if ( table_id == FA_STANDARD_DEFAULT ) {
			MembEtableOP etable_ptr( new MembEtable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
				EtableOptions() ) );
			memb_etables_[ table_id ] = etable_ptr;
		} else {
			std::string msg = "unrecognized etable: "+table_id;
			utility_exit_with_message( msg );
		}
	}
	return (memb_etables_.find( table_id )->second);
}

///////////////////////////////////////////////////////////////////////////////
// New logic for etable - 01/07/14 Hahnbeom Park
// one weak point of this method is that, change in frequently called but less relevant options
// such as weights, will invoke another etable construction.
// we can pro
etable::EtableCAP
ScoringManager::etable( methods::EnergyMethodOptions const &options_in ) const
{
	return etable( options_in.etable_options() );
}

etable::EtableCAP
ScoringManager::etable( etable::EtableOptions const & options_in ) const
{
	using namespace etable;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::map< EtableOptions, EtableOP >::const_iterator it = etables_by_options_.find( options_in );

	// add if no matching EtableOption is found
	if ( it == etables_by_options_.end() ) {
		EtableOP etable_ptr;
		std::string const & table_id = options_in.etable_type;
		if ( table_id == FA_STANDARD_SOFT ) {
			// soft rep etable: modified radii and also change to lj_switch_dis2sigma
			EtableOptions options_local( options_in );
			options_local.lj_switch_dis2sigma = 0.91;
			if ( option[corrections::beta_nov15 ]() || option[ corrections::beta_nov15_cart ]() ) {
				// hacky route for beta energy function
				etable_ptr = EtableOP( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
					options_local, "SOFTBETANOV15" ) );
			} else { // default
				etable_ptr = EtableOP( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
					options_local, "SOFT" ) );
			}
		} else if ( table_id == FA_STANDARD_MULTIPOLE ) {
			// multipole etable: change to lj_switch_dis2sigma to make harder repulsion.
			// Necessary to stop oppositely charged atoms from approaching each other.
			EtableOptions options_local( options_in );
			options_local.lj_switch_dis2sigma = 0.1;
			etable_ptr = EtableOP( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
				options_in ) );
		} else if ( utility::startswith( table_id, "FA_STANDARD_SOFT" ) ) {
			// add more softies, radii and lj_switch_dis2sigma are linear-interpolated from standard to SOFT
			// radii are given in database/chemical/atom_type_sets/fa_standard/extras/extra_soft_rep_params.txt
			EtableOptions options_loc( options_in );
			etable_ptr = make_partially_soft_etable( table_id, options_loc );
		} else if ( table_id.substr(0, FA_STANDARD_DEFAULT.size() + 1 ) == FA_STANDARD_DEFAULT+"_" ) {
			// original comments: note we check for soft rep 1st since that would match this as well -- confusing??
			// apply a modification of the radii/wdepths
			std::string const alternate_parameters( table_id.substr( FA_STANDARD_DEFAULT.size() + 1 ) );
			etable_ptr = EtableOP( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
				options_in, alternate_parameters ) );
		} else { // General way of adding
			etable_ptr = EtableOP( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
				options_in ) );
		}

		etables_by_options_[ options_in ] = etable_ptr;
		it = etables_by_options_.find( options_in );
	}

	return it->second;
}

etable::EtableCAP
ScoringManager::etable( std::string const & etable_id ) const
{
	etable::EtableOptions default_options;
	default_options.etable_type = etable_id;
	return etable( default_options );
}

/// @brief Get an owning pointer to data used by the FA_ElecEnergy in beta_nov15 mode.
/// @details If the data have not been loaded, this loads the data (lazy loading).  NOT THREADSAFE.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
core::scoring::elec::CPRepMapTypeCOP
ScoringManager::get_cp_rep_map_byname() const {
	if ( !cp_rep_map_byname_ ) {
		cp_rep_map_byname_ = core::scoring::elec::read_cp_tables_from_db( "scoring/score_functions/elec_cp_reps.dat" );
	}
	debug_assert( cp_rep_map_byname_ );
	return cp_rep_map_byname_;
}

/// @brief Get a vector of owning pointers to data used by the AACompositionEnergy score term.
/// @details If this vector has not yet been populated, this loads the data from disk (lazy loading).  NOT THREADSAFE.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
utility::vector1< core::scoring::aa_composition_energy::AACompositionEnergySetupOP >
ScoringManager::get_cloned_aa_comp_setup_helpers(
	core::scoring::methods::EnergyMethodOptions const &options
) const {
	core::Size const n_setup_helpers( options.aa_composition_setup_file_count() );

	//Load the data (once) if necessary:
	if ( aa_composition_setup_helpers_.size() != n_setup_helpers ) {
		aa_composition_setup_helpers_.clear();
		aa_composition_setup_helpers_.reserve( n_setup_helpers );
		for ( core::Size i=1; i<=n_setup_helpers; ++i ) {
			aa_composition_setup_helpers_.push_back( core::scoring::aa_composition_energy::AACompositionEnergySetupOP( new core::scoring::aa_composition_energy::AACompositionEnergySetup ) );
			aa_composition_setup_helpers_[i]->initialize_from_file( options.aa_composition_setup_file(i) );
		}
	}

	//Return a clone of the data:
	utility::vector1< core::scoring::aa_composition_energy::AACompositionEnergySetupOP > return_vect;
	return_vect.reserve( n_setup_helpers );
	for ( core::Size i=1; i<=n_setup_helpers; ++i ) {
		return_vect.push_back( aa_composition_setup_helpers_[i]->clone() );
	}
	return return_vect;
}

/// @brief Get a particular MainchainScoreTable for the rama_prepro score term.
/// @details If this has not yet been populated, loads the data from disk (lazy loading).  NOT THREADSAFE.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
core::chemical::mainchain_potential::MainchainScoreTableCOP
ScoringManager::get_rama_prepro_mainchain_torsion_potential(
	core::chemical::ResidueTypeCOP restype,
	bool const use_polycubic_interpolation,
	bool const prepro_table
) const {
	using namespace core::chemical::mainchain_potential;

	std::string const & mapname( restype->get_rama_prepro_mainchain_torsion_potential_name(prepro_table) );

	// First, check to see whether we need to create the potential.  THIS NEEDS TO BE MADE THREADSAFE.
	if (
			(!prepro_table && !rama_prepro_mainchain_potentials_.count( mapname )) ||
			(prepro_table && !rama_prepro_mainchain_potentials_beforeproline_.count( mapname ))
			) {
		std::string const mapfile( restype->get_rama_prepro_map_file_name(prepro_table) );
		if ( mapfile.empty() ) return MainchainScoreTableCOP(); //Return null pointer if there's no mapfile.
		utility::vector1< std::pair< std::string, MainchainScoreTableOP > > newtables;
		core::chemical::mainchain_potential::read_rama_map_file_shapovalov( mapfile, use_polycubic_interpolation, newtables ); //Read the file, and put results for ALL types in "newtables".
		bool mytype_found(false);
		for ( core::Size i=1, imax=newtables.size(); i<=imax; ++i ) {
			if ( newtables[i].first.compare( mapname ) == 0 ) mytype_found = true; //Found the type we're trying to load.
			if (
					( prepro_table && rama_prepro_mainchain_potentials_beforeproline_.count( newtables[i].first ) ) ||
					( ! prepro_table && rama_prepro_mainchain_potentials_.count( newtables[i].first ) )
					) {
				TR.Warning << "When trying to load RamaPrePro tables for residue type " << mapname << ", re-encountered tables for " << newtables[i].first << " (which were already loaded).  Skipping." << std::endl;
			} else {
				if ( prepro_table ) {
					rama_prepro_mainchain_potentials_beforeproline_[ newtables[i].first ] = newtables[i].second;
				} else {
					rama_prepro_mainchain_potentials_[ newtables[i].first ] = newtables[i].second;
				}
			}
		}
		runtime_assert_string_msg( mytype_found, "Error in core::scoring::ScoringManager::get_rama_prepro_mainchain_torsion_potential().  Could not load RamaPrePro scoring table for " + mapname + " from file " + mapfile + "." );
	}

	// Next, find and return the potenial.
	if ( prepro_table ) return rama_prepro_mainchain_potentials_beforeproline_.at( mapname ); //If we're fetching a preproline table.
	return rama_prepro_mainchain_potentials_.at( mapname ); //Otherwise.
}


// ///////////////////////////////////////////////////////////////////////////////
// // Support for previous version;
// // This function may not be called through usual score function calls,
// // but still necessary for whom direct call is necessary
// etable::EtableCAP
// ScoringManager::etable( std::string const & table_id ) const
// {
//  using namespace etable;
//
//  if ( etables_by_string_.find( table_id ) == etables_by_string_.end() ) {
//   // try to build if possible
//   if ( table_id == FA_STANDARD_DEFAULT ) {
//    EtableOP etable_ptr( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
//            EtableOptions() ) );
//    etables_by_string_[ table_id ] = etable_ptr;
//   } else if ( table_id == FA_STANDARD_SOFT ) {
//    // soft rep etable: modified radii and also change to lj_switch_dis2sigma
//    EtableOptions options;
//    options.lj_switch_dis2sigma = 0.91;
//    EtableOP etable_ptr( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
//            options, "SOFT" ) );
//    etables_by_string_[ table_id ] = etable_ptr;
//
//   } else if ( utility::startswith( table_id, "FA_STANDARD_SOFT" ) ) {
//    EtableOptions options;
//    EtableOP etable_ptr( make_partially_soft_etable( table_id, options ));
//    etables_by_string_[ table_id ] = etable_ptr;
//   } else if ( table_id.substr(0, FA_STANDARD_DEFAULT.size() + 1 ) == FA_STANDARD_DEFAULT+"_" ) {
//    // note we check for soft rep 1st since that would match this as well -- confusing??
//    // apply a modification of the radii/wdepths
//    std::string const alternate_parameters( table_id.substr( FA_STANDARD_DEFAULT.size() + 1 ) );
//    EtableOP etable_ptr( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
//            EtableOptions(), alternate_parameters ) );
//    etables_by_string_[ table_id ] = etable_ptr;
//
//   } else {
//    utility_exit_with_message("unrecognized etable: "+table_id );
//   }
//  }
//  return (etables_by_string_.find( table_id )->second);
// }

///////////////////////////////////////////////////////////////////////////////
/// @details Test if there is an EnergyMethod class defined for a
/// given score type.
bool
ScoringManager::has_energy_method(
	ScoreType score_type
) const {
	if ( score_type > n_score_types ) {
		return false;
	}

	if ( score_type == python ) return false;

	if ( method_creator_map_[score_type] == nullptr ) {
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

/// @details When a ScoreFunction the weight for a particular ScoreType set from 0
/// to some non-zero value, it will request an instance of the EnergyMethod class
/// that is responsible for calculating that ScoreType.  The ScoringManager responds
/// to that request by asking the EnergyMethodCreator that has claimed responsibility
/// for this ScoreType for a new instance.  EnergyMethodCreators must first have
/// registered themselves with the ScoringManager.  This should have been done at
/// load time, using a static-variable-initialization function call.
/// See src/core/scoring/etable/EtableEnergy.cc for an example of how the
/// EtableEnergyCreator class registers itself with the ScoringManager.
methods::EnergyMethodOP
ScoringManager::energy_method(
	ScoreType const & score_type,
	methods::EnergyMethodOptions const & options
) const
{
	if ( score_type > n_score_types ) {
		/// Inactive score type requested.  The program must be recompiled such that the desired score type
		/// appears before the n_score_types element in the ScoreType enumeration.  The program must now quit
		/// or it will later produce a segmentation fault when the EnergyMethod responsible for this term
		/// attempts to write into an EnergyMap object at a position that the EnergyMap has not allocated space
		/// for.
		std::cerr << "Critical error in ScoringManager::energy_method().\nRequested an inactive score_type '" << score_type;
		std::cerr << "' defined at position " << (int) score_type << " in the ScoreType enumeration.\n";
		std::cerr << "Active score types must appear before the n_score_types element ";
		std::cerr << "(at position " << (int) n_score_types << ") as this element marks the end of the active score types.\n";
		std::cerr << "Rosetta must be recompiled after src/core/scoring/ScoreType.hh is modified to include " << score_type;
		std::cerr << " as an active score type." << std::endl;
		utility_exit_with_message( "ERROR: Attempted to use an inactive score type" );
	}

	if ( score_type == python ) return nullptr; /// python special case; this could now be changed...

	if ( method_creator_map_[ score_type ] == nullptr ) {
		throw utility::excn::EXCN_Msg_Exception( "Requested ScoreType '" + utility::to_string( score_type ) + "' does not have a registered EnergyMethodCreator." );
	}

	return method_creator_map_[ score_type ]->create_energy_method( options );
}

/// global etable_id
std::string const FA_STANDARD_DEFAULT( "FA_STANDARD_DEFAULT" ); // keep this string the same as the etable_type in the default ctor in EtableOptions.cc
std::string const FA_STANDARD_SOFT   ( "FA_STANDARD_SOFT" );
std::string const FA_STANDARD_MULTIPOLE   ( "FA_STANDARD_MULTIPOLE" );

std::string const UNFOLDED_SCORE12( "UNFOLDED_SCORE12" );
std::string const UNFOLDED_MM_STD( "UNFOLDED_MM_STD" );
std::string const UNFOLDED_RNA( "UNFOLDED_RNA" ); // This will later get more elaborated

std::string const UNFOLDED_SPLIT_TALARIS2013( "UNFOLDED_SPLIT_TALARIS2013" ); //to use the split unfolded energy one body term.
std::string const UNFOLDED_SPLIT_MM_STD( "UNFOLDED_SPLIT_MM_STD" ); //to use the split unfolded energy one body term.

std::string const UNFOLDED_SPLIT_USER_DEFINED( "UNFOLDED_SPLIT_USER_DEFINED" );

std::string const SPLIT_UNFOLDED_ELE("SPLIT_UNFOLDED_ELE");
std::string const SPLIT_UNFOLDED_PDB("SPLIT_UNFOLDED_PDB");
std::string const SPLIT_UNFOLDED_ROSETTA("SPLIT_UNFOLDED_ROSETTA");
std::string const SPLIT_UNFOLDED_MM("SPLIT_UNFOLDED_MM");
std::string const SPLIT_UNFOLDED_UNIQUE("SPLIT_UNFOLDED_UNIQUE");

std::string const SPLIT_UNFOLDED_MEAN("SPLIT_UNFOLDED_MEAN");
std::string const SPLIT_UNFOLDED_MEDIAN("SPLIT_UNFOLDED_MEDIAN");
std::string const SPLIT_UNFOLDED_MODE("SPLIT_UNFOLDED_MODE");
std::string const SPLIT_UNFOLDED_BOLTZ("SPLIT_UNFOLDED_BOLTZ");

std::string const SPLIT_UNFOLDED_USER_DEFINED( "SPLIT_UNFOLDED_USER_DEFINED" );

} // namespace core
} // namespace scoring
