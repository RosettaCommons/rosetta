// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoringManager.cc
/// @brief  Scoring manager class.
/// @details  The ScoringManager handles the lazy loading of data for each scoretype.  Note that data load
/// must be threadsafe.  For this, the utility::thread::safely_create_load_once_object_by_OP function is used.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- added thread-safety to lazily loaded data.

// Unit headers
#include <core/scoring/ScoringManager.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
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
#include <core/scoring/P_AA_ABEGO3.hh>
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
	// Mutexes -- only exist in the threaded compilations.
#ifdef MULTI_THREADED
	pairE_mutex_(),
	genborn_mutex_(),
	hxl_potential_mutex_(),
	vdw_tinker_mutex_(),
	multipole_elec_mutex_(),
	sasa_potential_mutex_(),
	facts_mutex_(),
	dnabase_mutex_(),
	rama_mutex_(),
	rama2b_mutex_(),
	rama_pp_mutex_(),
	p_aa_abego3_mutex_(),
	dnabform_mutex_(),
	dnatorsion_mutex_(),
	omegatether_mutex_(),
	smoothenvpair_mutex_(),
	cenrotenvpair_mutex_(),
	cenhb_mutex_(),
	envpair_mutex_(),
	dnaenvpair_mutex_(),
	dnadihedral_mutex_(),
	secstruct_mutex_(),
	atomvdw_mutex_(),
	rna_atomvdw_mutex_(),
	database_occ_sol_mutex_(),
	carbonhbond_mutex_(),
	rna_suite_mutex_(),
	loopclose_sixdtransrot_mutex_(),
	rna_lowres_mutex_(),
	rnp_lowres_mutex_(),
	rnp_lowresstack_mutex_(),
	rna_chemshift_mutex_(),
	rna_dms_mutex_(),
	rna_dms_lowres_mutex_(),
	dna_directreadout_mutex_(),
	mm_lj_library_mutex_(),
	mm_lj_energytable_mutex_(),
	mm_torsionlibrary_mutex_(),
	mm_bondanglelibrary_mutex_(),
	mm_bondlengthlibrary_mutex_(),
	nv_lookup_mutex_(),
	orbitals_lookup_mutex_(),
	ddp_lookup_mutex_(),
	p_aa_mutex_(),
	p_aa_ss_mutex_(),
	unfoldedstate_mutex_(),
	wateradduct_hbond_mutex_(),
	membranepot_mutex_(),
	membranedata_mutex_(),
	membrane_fapot_mutex_(),
	proq_mutex_(),
	poissonboltzman_mutex_(),
	splitunfolded_2body_mutex_(),
	fa_disulf_potential_mutex_(),
	cent_disulf_potential_mutex_(),
	disulf_matching_potential_mutex_(),
	carb_chienergy_mutex_(),
	carb_omegapref_mutex_(),
	memb_etable_mutex_(),
	etable_mutex_(),
	cp_rep_map_mutex_(),
	aa_comp_mutex_(),
#ifdef OLDER_GCC
	pairE_bool_(false),
	genborn_bool_(false),
	hxl_potential_bool_(false),
	vdw_tinker_bool_(false),
	multipole_elec_bool_(false),
	sasa_potential_bool_(false),
	facts_bool_(false),
	dnabase_bool_(false),
	rama_bool_(false),
	rama2b_bool_(false),
	rama_pp_bool_(false),
	p_aa_abego3_bool_(false),
	dnabform_bool_(false),
	dnatorsion_bool_(false),
	omegatether_bool_(false),
	smoothenvpair_bool_(false),
	cenrotenvpair_bool_(false),
	cenhb_bool_(false),
	envpair_bool_(false),
	dnaenvpair_bool_(false),
	dnadihedral_bool_(false),
	secstruct_bool_(false),
	//atomvdw_bool_(false),
	rna_atomvdw_bool_(false),
	database_occ_sol_bool_(false),
	carbonhbond_bool_(false),
	//rna_suite_bool_(false),
	//loopclose_sixdtransrot_bool_(false),
	rna_lowres_bool_(false),
	rnp_lowres_bool_(false),
	rnp_lowresstack_bool_(false),
	rna_chemshift_bool_(false),
	rna_dms_bool_(false),
	rna_dms_lowres_bool_(false),
	dna_directreadout_bool_(false),
	mm_lj_library_bool_(false),
	mm_lj_energytable_bool_(false),
	mm_torsionlibrary_bool_(false),
	mm_bondanglelibrary_bool_(false),
	mm_bondlengthlibrary_bool_(false),
	nv_lookup_bool_(false),
	orbitals_lookup_bool_(false),
	ddp_lookup_bool_(false),
	p_aa_bool_(false),
	p_aa_ss_bool_(false),
	unfoldedstate_bool_(false),
	wateradduct_hbond_bool_(false),
	membranepot_bool_(false),
	membranedata_bool_(false),
	membrane_fapot_bool_(false),
	proq_bool_(false),
	poissonboltzman_bool_(false),
	splitunfolded_2body_bool_(false),
	fa_disulf_potential_bool_(false),
	cent_disulf_potential_bool_(false),
	disulf_matching_potential_bool_(false),
	carb_chienergy_bool_(false),
	carb_omegapref_bool_(false),
	cp_rep_map_bool_(false),
#endif //OLDER_GCC
#endif //MULTI_THREADED
	vdw_tinker_potential_( /* 0 */ ),
	pairE_potential_( /* 0 */ ),
	rama_( /* 0 */ ),
	rama2b_( /* 0 */ ),
	rama_pp_( /* 0 */ ),
	paa_abego3_( /* 0 */ ),
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
	etables_by_options_(),
	memb_etables_(),
	cp_rep_map_byname_(),
	aa_composition_setup_helpers_(),
	rama_prepro_mainchain_potentials_(),
	rama_prepro_mainchain_potentials_beforeproline_(),
#ifdef MULTI_THREADED
	rama_prepro_mainchain_potentials_mutex_(),
	rama_prepro_mainchain_potentials_beforeproline_mutex_(),
#endif
	method_creator_map_( n_score_types, nullptr )
{}

/// @brief Get an instance of the VdWTinkerPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
VdWTinkerPotential const &
ScoringManager::get_VdWTinkerPotential() const
{
	boost::function< VdWTinkerPotentialOP () > creator( boost::bind( &ScoringManager::create_vdw_tinker_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, vdw_tinker_potential_, SAFELY_PASS_MUTEX(vdw_tinker_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(vdw_tinker_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *vdw_tinker_potential_;
}

/// @brief Get a const instance of the PairEPotential.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan.PairEPotential const &
PairEPotential const &
ScoringManager::get_PairEPotential() const
{
	boost::function< PairEPotentialOP () > creator( boost::bind( &ScoringManager::create_pairE_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, pairE_potential_, SAFELY_PASS_MUTEX(pairE_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(pairE_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *pairE_potential_;
}

/// @brief The ScoringManager acts as an EnergyMethodFactory.  All EnergyMethods must
/// create a helper class, an EnergyMethodCreator class, that will respond to a call to
/// its create_energy_method by returning a new instance of that EnergyMethod its helping.
/// This Creator class must also register itself with the ScoringManager at load time and
/// hand an instance of itself to the singleton ScoringManager instance.
/// @details I don't think that this function is threadsafe (VKM, 20 July 2017), but it probably
/// doesn't matter.  Factory registration presumably happens during Rosetta initialization,
/// before any threads are spawned.
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

/// @brief Get an instance of the DNA_BasePotential scoring object, by const owning pointer.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan.
dna::DNA_BasePotential const &
ScoringManager::get_DNA_BasePotential() const
{
	boost::function< dna::DNA_BasePotentialOP () > creator( boost::bind( &ScoringManager::create_dnabase_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, DNA_base_potential_, SAFELY_PASS_MUTEX(dnabase_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(dnabase_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *DNA_base_potential_;
}

/// @brief Get a const instance of the GenBornPotential.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan.
GenBornPotential const &
ScoringManager::get_GenBornPotential() const
{
	boost::function< GenBornPotentialOP () > creator( boost::bind( &ScoringManager::create_genborn_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, gen_born_potential_, SAFELY_PASS_MUTEX(genborn_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(genborn_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *gen_born_potential_;
}

/// @brief Get a const instance of the HydroxylTorsionPotential.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan.
HydroxylTorsionPotential const &
ScoringManager::get_HydroxylTorsionPotential() const
{
	boost::function< HydroxylTorsionPotentialOP () > creator( boost::bind( &ScoringManager::create_hxl_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, hxl_tors_potential_, SAFELY_PASS_MUTEX(hxl_potential_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(hxl_potential_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *hxl_tors_potential_;
}


/// @brief Get an instance of the MultipoleElecPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
MultipoleElecPotential const &
ScoringManager::get_MultipoleElecPotential( methods::EnergyMethodOptions const & options ) const
{
	boost::function< MultipoleElecPotentialOP () > creator( boost::bind( &ScoringManager::create_multipole_elec_instance, options ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, multipole_elec_potential_, SAFELY_PASS_MUTEX(multipole_elec_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(multipole_elec_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	//if ( multipole_elec_potential_ == nullptr ) {
	// multipole_elec_potential_ = MultipoleElecPotentialOP( new MultipoleElecPotential() );
	// multipole_elec_potential_->use_polarization = options.use_polarization();
	// multipole_elec_potential_->use_gen_kirkwood = options.use_gen_kirkwood();
	// multipole_elec_potential_->Ep = options.protein_dielectric();
	// multipole_elec_potential_->Ew = options.water_dielectric();
	//}
	return *multipole_elec_potential_;
}

/// @brief Get an instance of the SASAPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
SASAPotential const &
ScoringManager::get_SASAPotential() const
{
	boost::function< SASAPotentialOP () > creator( boost::bind( &ScoringManager::create_sasa_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, sasa_potential_, SAFELY_PASS_MUTEX(sasa_potential_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(sasa_potential_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *sasa_potential_;
}

/// @brief Get an instance of the FACTSPotential scoring object, by const owning pointer.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan.
FACTSPotential const &
ScoringManager::get_FACTSPotential() const
{
	boost::function< FACTSPotentialOP () > creator( boost::bind( &ScoringManager::create_facts_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, facts_potential_, SAFELY_PASS_MUTEX(facts_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(facts_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *facts_potential_;
}

/// @brief Get an instance of the P_AA scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
P_AA const &
ScoringManager::get_P_AA() const
{
	boost::function< P_AAOP () > creator( boost::bind( &ScoringManager::create_p_aa_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, p_aa_, SAFELY_PASS_MUTEX(p_aa_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(p_aa_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *p_aa_;
}

/// @brief Get an instance of the P_AA_ss scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
P_AA_ss const &
ScoringManager::get_P_AA_ss() const
{
	boost::function< P_AA_ssOP () > creator( boost::bind( &ScoringManager::create_p_aa_ss_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, p_aa_ss_, SAFELY_PASS_MUTEX(p_aa_ss_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(p_aa_ss_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *p_aa_ss_;
}

/// @brief Get an instance of the Ramachandran scoring object, by const owning pointer.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan.
RamachandranCOP
ScoringManager::get_Ramachandran_ptr() const
{
	boost::function< RamachandranOP () > creator( boost::bind( &ScoringManager::create_rama_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, rama_, SAFELY_PASS_MUTEX(rama_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(rama_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return rama_;
}

/// @brief Get a const instance of the Ramachandran scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan.
Ramachandran const &
ScoringManager::get_Ramachandran() const
{
	return *get_Ramachandran_ptr();
}

/// @brief Get an instance of the Ramachandran2B scoring object, by const owning pointer.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
Ramachandran2BCOP
ScoringManager::get_Ramachandran2B_ptr() const
{
	boost::function< Ramachandran2BOP () > creator( boost::bind( &ScoringManager::create_rama2b_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, rama2b_, SAFELY_PASS_MUTEX(rama2b_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(rama2b_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return rama2b_;
}

/// @brief Get an instance of the Ramachandran2B scoring object, by const instance.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
Ramachandran2B const &
ScoringManager::get_Ramachandran2B() const
{
	return *get_Ramachandran2B_ptr();
}

///////////////////////////////////////////////////////////////////////////////
RamaPrePro const &
ScoringManager::get_RamaPrePro() const
{
	boost::function< RamaPreProOP () > creator( boost::bind( &ScoringManager::create_ramapp_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, rama_pp_, SAFELY_PASS_MUTEX(rama_pp_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(rama_pp_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *rama_pp_;
}


/// @brief Get an instance of the P_AA_ABEGO3 scoring object.
/// @details Threadsafe and lazily loaded.  Used by AbegoEnergy.
/// @author imv@uw.edu
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
P_AA_ABEGO3 const &
ScoringManager::get_P_AA_ABEGO3() const
{
	boost::function< P_AA_ABEGO3_OP () > creator( boost::bind( &ScoringManager::create_p_aa_abego3_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, paa_abego3_, SAFELY_PASS_MUTEX(p_aa_abego3_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(p_aa_abego3_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *paa_abego3_;
}

/// @brief Get an instance of the DNABFormPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
dna::DNABFormPotential const &
ScoringManager::get_DNABFormPotential() const
{
	boost::function< dna::DNABFormPotentialOP () > creator( boost::bind( &ScoringManager::create_dna_bform_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, dnabform_, SAFELY_PASS_MUTEX(dnabform_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(dnabform_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *dnabform_;
}

/// @brief Get an instance of the DNATorsionPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
dna::DNATorsionPotential const &
ScoringManager::get_DNATorsionPotential() const
{
	boost::function< dna::DNATorsionPotentialOP () > creator( boost::bind( &ScoringManager::create_dna_torsion_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, dna_torsion_potential_, SAFELY_PASS_MUTEX(dnatorsion_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(dnatorsion_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *dna_torsion_potential_;
}

/// @brief Get an instance of the OmegaTether scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
OmegaTether const &
ScoringManager::get_OmegaTether() const
{
	boost::function< OmegaTetherOP () > creator( boost::bind( &ScoringManager::create_omegatether_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, omega_, SAFELY_PASS_MUTEX(omegatether_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(omegatether_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *omega_;
}

/// @brief Get an instance of the SmoothEnvPairPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
SmoothEnvPairPotential const &
ScoringManager::get_SmoothEnvPairPotential() const
{
	boost::function< SmoothEnvPairPotentialOP () > creator( boost::bind( &ScoringManager::create_smoothenvpair_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, smooth_env_pair_potential_, SAFELY_PASS_MUTEX(smoothenvpair_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(smoothenvpair_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *smooth_env_pair_potential_;
}

/// @brief Get an instance of the CenRotEnvPairPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
CenRotEnvPairPotential const &
ScoringManager::get_CenRotEnvPairPotential() const
{
	boost::function< CenRotEnvPairPotentialOP () > creator( boost::bind( &ScoringManager::create_cenrotenvpair_instance) );
	utility::thread::safely_create_load_once_object_by_OP( creator, cen_rot_pair_potential_, SAFELY_PASS_MUTEX(cenrotenvpair_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(cenrotenvpair_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *cen_rot_pair_potential_;
}

/// @brief Get an instance of the CenHBPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
CenHBPotential const &
ScoringManager::get_CenHBPotential() const
{
	boost::function< CenHBPotentialOP () > creator( boost::bind( &ScoringManager::create_cenhbpotential_instance) );
	utility::thread::safely_create_load_once_object_by_OP( creator, cen_hb_potential_, SAFELY_PASS_MUTEX(cenhb_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(cenhb_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *cen_hb_potential_;
}

/// @brief Get an instance of the EnvPairPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
EnvPairPotential const &
ScoringManager::get_EnvPairPotential() const
{
	boost::function< EnvPairPotentialOP () > creator( boost::bind( &ScoringManager::create_envpairpotential_instance) );
	utility::thread::safely_create_load_once_object_by_OP( creator, env_pair_potential_, SAFELY_PASS_MUTEX(envpair_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(envpair_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *env_pair_potential_;
}

/// @brief Get an instance of the DNA_EnvPairPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
dna::DNA_EnvPairPotential const &
ScoringManager::get_DNA_EnvPairPotential() const
{
	boost::function< dna::DNA_EnvPairPotentialOP () > creator( boost::bind( &ScoringManager::create_dna_envpairpotential_instance) );
	utility::thread::safely_create_load_once_object_by_OP( creator, dna_env_pair_potential_, SAFELY_PASS_MUTEX( dnaenvpair_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( dnaenvpair_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *dna_env_pair_potential_;
}

/// @brief Get an instance of the DNA_DihedralPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
dna::DNA_DihedralPotential const &
ScoringManager::get_DNA_DihedralPotential() const
{
	boost::function< dna::DNA_DihedralPotentialOP () > creator( boost::bind( &ScoringManager::create_dna_dihedralpotential_instance) );
	utility::thread::safely_create_load_once_object_by_OP( creator, dna_dihedral_potential_, SAFELY_PASS_MUTEX( dnadihedral_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( dnadihedral_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *dna_dihedral_potential_;
}

/// @brief Get an instance of the SecondaryStructurePotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
SecondaryStructurePotential const &
ScoringManager::get_SecondaryStructurePotential() const
{
	boost::function< SecondaryStructurePotentialOP () > creator( boost::bind( &ScoringManager::create_secondarystructurepotential_instance) );
	utility::thread::safely_create_load_once_object_by_OP( creator, secondary_structure_potential_, SAFELY_PASS_MUTEX( secstruct_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( secstruct_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *secondary_structure_potential_;
}

/// @brief Get an instance of the AtomVDW scoring object.
/// @details Threadsafe and lazily loaded.
/// @note Each element in the atom_vdw_ map is now threadsafe and lazily loaded (independently).
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
AtomVDW const &
ScoringManager::get_AtomVDW( std::string const & atom_type_set_name ) const
{
	boost::function< AtomVDWOP () > creator( boost::bind( &ScoringManager::create_atomvdw_instance, atom_type_set_name ) );
	return *( utility::thread::safely_check_map_for_key_and_insert_if_absent( creator, SAFELY_PASS_MUTEX( atomvdw_mutex_ ), atom_type_set_name, atom_vdw_ ) );
}

/// @brief Get an instance of the RNA_AtomVDW scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
rna::RNA_AtomVDW const &
ScoringManager::get_RNA_AtomVDW() const
{
	boost::function< rna::RNA_AtomVDWOP () > creator( boost::bind( &ScoringManager::create_rna_atomvdw_instance) );
	utility::thread::safely_create_load_once_object_by_OP( creator, rna_atom_vdw_, SAFELY_PASS_MUTEX( rna_atomvdw_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( rna_atomvdw_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *rna_atom_vdw_;
}

/// @brief Get an instance of the DatabaseOccSolEne scoring object.
/// @details Threadsafe and lazily loaded.
/// @note Whatever atom type set name and min occ energy are passed to this function the FIRST time determine the
/// object that gets created.  These parameters are unused in subsequent invocations.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
geometric_solvation::DatabaseOccSolEne const &
ScoringManager::get_DatabaseOccSolEne( std::string const & atom_type_set_name, Real const & min_occ_energy ) const
{
	boost::function< geometric_solvation::DatabaseOccSolEneOP () > creator( boost::bind( &ScoringManager::create_database_occsolene_instance, atom_type_set_name, min_occ_energy ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, occ_hbond_sol_database_, SAFELY_PASS_MUTEX( database_occ_sol_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( database_occ_sol_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *occ_hbond_sol_database_;
}

/// @brief Get an instance of the CarbonHBondPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
carbon_hbonds::CarbonHBondPotential const &
ScoringManager::get_CarbonHBondPotential() const
{
	boost::function< carbon_hbonds::CarbonHBondPotentialOP () > creator( boost::bind( &ScoringManager::create_carbon_hbond_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, carbon_hbond_potential_, SAFELY_PASS_MUTEX( carbonhbond_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( carbonhbond_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *carbon_hbond_potential_;
}

/// @brief Get an instance of the RNA_SuitePotential scoring object, by const owning pointer.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
rna::RNA_SuitePotentialCOP
ScoringManager::get_rna_suite_potential( bool const & calculate_suiteness_bonus, std::string const & suiteness_bonus ) const
{
	std::pair< bool, std::string > const key = std::make_pair( calculate_suiteness_bonus, suiteness_bonus );
	boost::function< rna::RNA_SuitePotentialOP () > creator( boost::bind( &ScoringManager::create_rna_suitepotential_instance, calculate_suiteness_bonus, suiteness_bonus ) );
	return rna::RNA_SuitePotentialCOP(
		utility::thread::safely_check_map_for_key_and_insert_if_absent( creator, SAFELY_PASS_MUTEX( rna_suite_mutex_ ), key, rna_suite_potential_ )
	);
}

/// @brief Get an instance of the SixDTransRotPotential scoring object, by const owning pointer.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
loop_graph::evaluator::SixDTransRotPotentialCOP
ScoringManager::get_LoopCloseSixDPotential( std::string const & database_file ) const
{
	boost::function< loop_graph::evaluator::SixDTransRotPotentialOP () > creator( boost::bind( &ScoringManager::create_sixdtransrotpotential_instance, database_file ) );
	return loop_graph::evaluator::SixDTransRotPotentialCOP(
		utility::thread::safely_check_map_for_key_and_insert_if_absent( creator, SAFELY_PASS_MUTEX( loopclose_sixdtransrot_mutex_ ), database_file, loop_close_six_d_potential_ )
	);
}

/// @brief Get an instance of the RNA_LowResolutionPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
rna::RNA_LowResolutionPotential const &
ScoringManager::get_RNA_LowResolutionPotential() const
{
	boost::function< rna::RNA_LowResolutionPotentialOP () > creator( boost::bind( &ScoringManager::create_rna_lowresolutionpotential_instance) );
	utility::thread::safely_create_load_once_object_by_OP( creator, rna_low_resolution_potential_, SAFELY_PASS_MUTEX( rna_lowres_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( rna_lowres_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *rna_low_resolution_potential_;
}

/// @brief Get an instance of the RNP_LowResPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
rna::RNP_LowResPotential const &
ScoringManager::get_RNP_LowResPotential() const
{
	boost::function< rna::RNP_LowResPotentialOP () > creator( boost::bind( &ScoringManager::create_rnp_lowrespotential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, rnp_low_res_potential_, SAFELY_PASS_MUTEX( rnp_lowres_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( rnp_lowres_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *rnp_low_res_potential_;
}

/// @brief Get an instance of the RNP_LowResStackData scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).rna::RNP_LowResStackData const &
rna::RNP_LowResStackData const &
ScoringManager::get_RNP_LowResStackData() const
{
	boost::function< rna::RNP_LowResStackDataOP () > creator( boost::bind( &ScoringManager::create_rnp_lowresstackdata_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, rnp_low_res_stack_data_, SAFELY_PASS_MUTEX( rnp_lowresstack_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( rnp_lowresstack_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *rnp_low_res_stack_data_;
}

/// @brief Get an instance of the RNA_ChemicalShiftPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
rna::chemical_shift::RNA_ChemicalShiftPotential const &
ScoringManager::get_RNA_ChemicalShiftPotential() const
{
	boost::function< rna::chemical_shift::RNA_ChemicalShiftPotentialOP () > creator( boost::bind( &ScoringManager::create_rna_chemshiftpotential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, rna_chemical_shift_potential_, SAFELY_PASS_MUTEX( rna_chemshift_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( rna_chemshift_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *rna_chemical_shift_potential_;
}

/// @brief Get an instance of the RNA_DMS_Potential scoring object.
/// @details Threadsafe and lazily loaded.
/// @note The RNA_DMS_Potential itself is fundamentally NOT THREADSAFE!!!
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
rna::data::RNA_DMS_Potential &
ScoringManager::get_RNA_DMS_Potential() const
{
	boost::function< rna::data::RNA_DMS_PotentialOP () > creator( boost::bind( &ScoringManager::create_rna_dms_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, rna_dms_potential_, SAFELY_PASS_MUTEX( rna_dms_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( rna_dms_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *rna_dms_potential_;
}

/// @brief Get an instance of the RNA_DMS_LowResolutionPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @note The RNA_DMS_LowResolutionPotential itself is fundamentally NOT THREADSAFE!!!
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
rna::data::RNA_DMS_LowResolutionPotential &
ScoringManager::get_RNA_DMS_LowResolutionPotential() const
{
	boost::function< rna::data::RNA_DMS_LowResolutionPotentialOP () > creator( boost::bind( &ScoringManager::create_rna_dms_lowrespotential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, rna_dms_low_resolution_potential_, SAFELY_PASS_MUTEX( rna_dms_lowres_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( rna_dms_lowres_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *rna_dms_low_resolution_potential_;
}

/// @brief Get an instance of the DirectReadoutPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
dna::DirectReadoutPotential const &
ScoringManager::get_DirectReadoutPotential() const
{
	boost::function< dna::DirectReadoutPotentialOP () > creator( boost::bind( &ScoringManager::create_dna_directreadoutpotential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, dna_dr_potential_, SAFELY_PASS_MUTEX( dna_directreadout_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( dna_directreadout_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *dna_dr_potential_;
}

/// @brief Get an instance of the MMLJLibrary scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
core::scoring::mm::MMLJLibrary const &
ScoringManager::get_MMLJLibrary() const
{
	boost::function< core::scoring::mm::MMLJLibraryOP () > creator( boost::bind( &ScoringManager::create_mm_lj_library_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, mm_lj_library_, SAFELY_PASS_MUTEX( mm_lj_library_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( mm_lj_library_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *mm_lj_library_;
}

/// @brief Get an instance of the MMLJEnergyTable scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
core::scoring::mm::MMLJEnergyTable const &
ScoringManager::get_MMLJEnergyTable () const
{
	boost::function< core::scoring::mm::MMLJEnergyTableOP () > creator( boost::bind( &ScoringManager::create_mm_lj_energy_table_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, mm_lj_energy_table_, SAFELY_PASS_MUTEX( mm_lj_energytable_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( mm_lj_energytable_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *mm_lj_energy_table_;
}

/// @brief Get an instance of the MMTorsionLibrary scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
core::scoring::mm::MMTorsionLibrary const &
ScoringManager::get_MMTorsionLibrary() const
{
	boost::function< mm::MMTorsionLibraryOP () > creator( boost::bind( &ScoringManager::create_mm_torsion_library_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, mm_torsion_library_, SAFELY_PASS_MUTEX( mm_torsionlibrary_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( mm_torsionlibrary_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *mm_torsion_library_;
}

/// @brief Get an instance of the MMBondAngleLibrary scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
core::scoring::mm::MMBondAngleLibrary const &
ScoringManager::get_MMBondAngleLibrary() const
{
	boost::function< mm::MMBondAngleLibraryOP () > creator( boost::bind( &ScoringManager::create_mm_bondangle_library_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, mm_bondangle_library_, SAFELY_PASS_MUTEX( mm_bondanglelibrary_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( mm_bondanglelibrary_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *mm_bondangle_library_;
}

/// @brief Get an instance of the MMBondLengthLibrary scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
core::scoring::mm::MMBondLengthLibrary const &
ScoringManager::get_MMBondLengthLibrary() const
{
	boost::function< mm::MMBondLengthLibraryOP () > creator( boost::bind( &ScoringManager::create_mm_bondlength_library_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, mm_bondlength_library_, SAFELY_PASS_MUTEX( mm_bondlengthlibrary_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( mm_bondlengthlibrary_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *mm_bondlength_library_;
}

/// @brief Get an instance of the NVlookup scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
nv::NVlookup const &
ScoringManager::get_NVLookupTable() const
{
	boost::function< nv::NVlookupOP () > creator( boost::bind( &ScoringManager::create_nvlookup_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator,  NV_lookup_table_, SAFELY_PASS_MUTEX( nv_lookup_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( nv_lookup_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *NV_lookup_table_;
}

/// @brief Get an instance of the OrbitalsLookup scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
orbitals::OrbitalsLookup const &
ScoringManager::get_OrbitalsLookupTable() const
{
	boost::function< orbitals::OrbitalsLookupOP () > creator( boost::bind( &ScoringManager::create_orbitals_lookup_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator,  orbitals_lookup_table_, SAFELY_PASS_MUTEX( orbitals_lookup_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( orbitals_lookup_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *orbitals_lookup_table_;
}

/// @brief Get an instance of the DDPlookup scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
interface_::DDPlookup const &
ScoringManager::get_DDPLookupTable() const
{
	boost::function< interface_::DDPlookupOP () > creator( boost::bind( &ScoringManager::create_ddp_lookup_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator,  DDP_lookup_table_, SAFELY_PASS_MUTEX( ddp_lookup_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( ddp_lookup_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *DDP_lookup_table_;
}

/// @brief Get an instance of the UnfoldedStatePotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
UnfoldedStatePotential const &
ScoringManager::get_UnfoldedStatePotential( std::string const & type ) const
{
	boost::function< UnfoldedStatePotentialOP () > creator( boost::bind( &ScoringManager::create_unfolded_state_potential_instance, type ) );
	utility::thread::safely_create_load_once_object_by_OP( creator,  unf_state_, SAFELY_PASS_MUTEX( unfoldedstate_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( unfoldedstate_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *unf_state_;
}

/// @brief Get an instance of the WaterAdductHBondPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
WaterAdductHBondPotential const &
ScoringManager::get_WaterAdductHBondPotential() const
{
	boost::function< WaterAdductHBondPotentialOP () > creator( boost::bind( &ScoringManager::create_water_adduct_hbond_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator,  water_adduct_hbond_potential_, SAFELY_PASS_MUTEX( wateradduct_hbond_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( wateradduct_hbond_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *water_adduct_hbond_potential_;
}

/// @brief Get an instance of the MembranePotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
MembranePotential const &
ScoringManager::get_MembranePotential() const
{
	boost::function< MembranePotentialOP () > creator( boost::bind( &ScoringManager::create_membrane_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator,  membrane_potential_, SAFELY_PASS_MUTEX( membranepot_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( membranepot_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *membrane_potential_;
}

/// @brief Get an instance of the MembraneData scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
membrane::MembraneData const &
ScoringManager::get_MembraneData() const
{
	boost::function< membrane::MembraneDataOP () > creator( boost::bind( &ScoringManager::create_membrane_data_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator,  mp_base_potential_, SAFELY_PASS_MUTEX( membranedata_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( membranedata_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *mp_base_potential_;
}

/// @brief Get an instance of the Membrane_FAPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @note The Membrane_FAPotential object is fundamentally NOT THREADSAFE!!!
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
Membrane_FAPotential const &
ScoringManager::get_Membrane_FAPotential() const //pba
{
	boost::function< Membrane_FAPotentialOP () > creator( boost::bind( &ScoringManager::create_membrane_fa_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, membrane_fapotential_, SAFELY_PASS_MUTEX( membrane_fapot_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( membrane_fapot_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *membrane_fapotential_;
}

/// @brief Get an instance of the ProQPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
ProQPotential const &
ScoringManager::get_ProQPotential() const
{
	boost::function< ProQPotentialOP () > creator( boost::bind( &ScoringManager::create_proq_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, ProQ_potential_, SAFELY_PASS_MUTEX( proq_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( proq_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *ProQ_potential_;
}

/// @brief Get an instance of the PoissonBoltzmannPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
PoissonBoltzmannPotential const &
ScoringManager::get_PoissonBoltzmannPotential() const
{
	boost::function< PoissonBoltzmannPotentialOP () > creator( boost::bind( &ScoringManager::create_poisson_boltzmann_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, PB_potential_, SAFELY_PASS_MUTEX( poissonboltzman_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( poissonboltzman_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *PB_potential_;
}

/// @brief Get an instance of the SplitUnfoldedTwoBodyPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
SplitUnfoldedTwoBodyPotential const &
ScoringManager::get_SplitUnfoldedTwoBodyPotential(std::string const & label_type,std::string const & value_type, std::string const & score_func_type) const
{
	boost::function< SplitUnfoldedTwoBodyPotentialOP () > creator( boost::bind( &ScoringManager::create_split_unfolded_2body_potential_instance, label_type, value_type, score_func_type ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, sutbp_, SAFELY_PASS_MUTEX( splitunfolded_2body_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( splitunfolded_2body_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *sutbp_;
}

/// @brief Get an instance of the FullatomDisulfidePotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
disulfides::FullatomDisulfidePotential const &
ScoringManager::get_FullatomDisulfidePotential() const
{
	boost::function< disulfides::FullatomDisulfidePotentialOP () > creator( boost::bind( &ScoringManager::create_fullatom_disulfide_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, fa_disulfide_potential_, SAFELY_PASS_MUTEX( fa_disulf_potential_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( fa_disulf_potential_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *fa_disulfide_potential_;
}

/// @brief Get an instance of the CentroidDisulfidePotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
disulfides::CentroidDisulfidePotential const &
ScoringManager::get_CentroidDisulfidePotential() const
{
	boost::function< disulfides::CentroidDisulfidePotentialOP () > creator( boost::bind( &ScoringManager::create_centroid_disulfide_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, cen_disulfide_potential_, SAFELY_PASS_MUTEX( cent_disulf_potential_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( cent_disulf_potential_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *cen_disulfide_potential_;
}

/// @brief Get an instance of the DisulfideMatchingPotential scoring object.
/// @details Threadsafe and lazily loaded.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
disulfides::DisulfideMatchingPotential const &
ScoringManager::get_DisulfideMatchingPotential() const
{
	boost::function< disulfides::DisulfideMatchingPotentialOP () > creator( boost::bind( &ScoringManager::create_disulfide_matching_potential_instance ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, disulfide_matching_potential_, SAFELY_PASS_MUTEX( disulf_matching_potential_mutex_ ), SAFELY_PASS_THREADSAFETY_BOOL( disulf_matching_potential_bool_ ) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return *disulfide_matching_potential_;
}

/// @brief Get an instance of the CHIEnergyFunction scoring object.
/// @details Lazily loaded, but fundamentally NOT THREADSAFE.
carbohydrates::CHIEnergyFunction const &
ScoringManager::get_CHIEnergyFunction( bool const setup_for_sampling /* false */, core::Real const & step_size /* 0.1 */ ) const
{
#ifdef MULTI_THREADED
	utility_exit_with_message( "Error in ScoringManager: the carbohydrate CHIEnergyFunction is fundamentally not threadsafe, and cannot be used in a multithreaded environment.  Please contact Jason Labonte (JWLabonte@jhu.edu) to complain about this." );
#endif

	if ( CHI_energy_function_ == nullptr ) {
		TR << "Creating CHI Energy Function." << std::endl;
		CHI_energy_function_ = carbohydrates::CHIEnergyFunctionOP( new carbohydrates::CHIEnergyFunction );
	}

	// VKM 20 July 2017: The following is fundamentally not threadsafe, since different threads could simultaneously
	// be trying to configure the global CHI energy function for sampling or not for sampling.
	if ( setup_for_sampling && ( ! CHI_energy_function_->sampling_data_setup() ) ) {
		TR << "should be setting up for sampling..." << std::endl;
		CHI_energy_function_->setup_for_sampling( step_size );
	}

	return *CHI_energy_function_;
}

/// @brief Get an instance of the OmegaPreferencesFunction scoring object.
/// @details Lazily loaded, but fundamentally NOT THREADSAFE.
carbohydrates::OmegaPreferencesFunction const &
ScoringManager::get_OmegaPreferencesFunction(  bool const setup_for_sampling /* false */, core::Real const & step_size /* 0.1 */ ) const
{
#ifdef MULTI_THREADED
	utility_exit_with_message( "Error in ScoringManager: the carbohydrate OmegaPreferencesFunction is fundamentally not threadsafe, and cannot be used in a multithreaded environment.  Please contact Jason Labonte (JWLabonte@jhu.edu) to complain about this." );
#endif

	if ( carbohydrate_omega_preferences_function_ == nullptr ) {
		TR << "Creating carbohydrate omega preferences function." << std::endl;
		carbohydrate_omega_preferences_function_ =
			carbohydrates::OmegaPreferencesFunctionOP( new carbohydrates::OmegaPreferencesFunction );

	}

	// VKM 20 July 2017: The following is fundamentally not threadsafe, since different threads could simultaneously
	// be trying to configure the global OmegaPreferences energy function for sampling or not for sampling.
	if ( setup_for_sampling && ( ! carbohydrate_omega_preferences_function_->sampling_data_setup() ) ) {
		TR << "should be setting up for sampling..." << std::endl;
		carbohydrate_omega_preferences_function_->setup_for_sampling( step_size );
	}

	return *carbohydrate_omega_preferences_function_;
}

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
	return etable::EtableOP(nullptr); //To keep all compilers happy.
}


/// @brief Add a new membrane energy table to the membrane energy tables map.
/// @details Made threadsafe on 25 July 2017.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
void
ScoringManager::add_memb_etable( std::string const & name, etable::MembEtableOP etable ) //pba
{
#ifndef MULTI_THREADED
	debug_assert( memb_etables_.count(name) == 0 );
#endif
	boost::function< etable::MembEtableOP () > builder( boost::bind( &ScoringManager::create_memb_etable_instance_silly, etable ) ); //Note that this is a silly "builder" function that just returns the object passed to it.
	utility::thread::safely_check_map_for_key_and_insert_if_absent( builder, SAFELY_PASS_MUTEX( memb_etable_mutex_ ), name, memb_etables_ );
}

/// @brief Get a membrane energy table from the membrane energy tables map.
/// @details Made threadsafe on 25 July 2017.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu)
etable::MembEtableCAP
ScoringManager::memb_etable( std::string const & table_id ) const //pba
{
	boost::function< etable::MembEtableOP () > builder( boost::bind( &ScoringManager::create_memb_etable_instance, table_id ) ); //Note that this calls the non-silly version of the function.
	return ( utility::thread::safely_check_map_for_key_and_insert_if_absent(builder, SAFELY_PASS_MUTEX(memb_etable_mutex_), table_id, memb_etables_ ) );
}

/// @brief Request an etable specifying an EnergyMethodOptions object; internally
/// this will retrieve the EtableOptions object, and invoke the EtableOptions
/// version of this function.
/// @details Threadsafe, since this calls the threadsafe etable( EtableOptions) function, below.
/// @author Hahnbeom Park (new logic for etable)
/// @note One weak point of this method is that, change in frequently called
/// but less relevant options such as weights, will invoke another etable
/// construction.
etable::EtableCAP
ScoringManager::etable( methods::EnergyMethodOptions const &options_in ) const
{
	return etable( options_in.etable_options() );
}

/// @brief Request an etable specifying an EtableOptions; internally this will
/// query the ScoringManager's map from EtableOptions to Etables for the desired
/// Etable, and construct a new one if needed.
/// @details Made threadsafe on 25 July 2017.
/// @author Thread-safety added by Vikram K. Mulligan (vmullig@uw.edu).
etable::EtableCAP
ScoringManager::etable( etable::EtableOptions const & options_in ) const
{
	using namespace etable;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::map< EtableOptions, EtableOP >::const_iterator it, it_end;

	{ //Scope for the read lock
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard lock( etable_mutex_ );
#endif
		it = etables_by_options_.find( options_in );
		it_end = etables_by_options_.end();
	} //End scope for read lock

	// add if no matching EtableOption is found
	if ( it == it_end ) {
		EtableOP etable_ptr;
		std::string const & table_id = options_in.etable_type;
		if ( table_id == FA_STANDARD_SOFT ) {
			// soft rep etable: modified radii and also change to lj_switch_dis2sigma
			EtableOptions options_local( options_in );
			options_local.lj_switch_dis2sigma = 0.91;
			/*
			if ( option[corrections::beta_nov15 ]() || option[ corrections::beta_nov15_cart ]() ) {
			// hacky route for beta energy function
			etable_ptr = EtableOP( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
			options_local, "SOFTBETANOV15" ) );
			*/
			if ( option[ mistakes::restore_pre_talaris_2013_behavior ]() || option[ corrections::restore_talaris_behavior ]() ) {
				etable_ptr = EtableOP( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
					options_local, "SOFT" ) );
			} else { // default
				etable_ptr = EtableOP( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
					options_local, "SOFTBETANOV15" ) );
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

#ifdef MULTI_THREADED
		boost::function< EtableOP () > builder( boost::bind( &ScoringManager::create_etable_instance, etable_ptr ) ); //Note that this is a silly "builder" function that just returns the object passed to it.
		it = utility::thread::create_and_insert( builder, SAFELY_PASS_MUTEX( etable_mutex_ ), options_in, etables_by_options_ );
#else
		etables_by_options_[ options_in ] = etable_ptr;
		it = etables_by_options_.find( options_in );
#endif

	}

	return it->second;
}

/// @brief Create and return an etable specified only by the etable_type of the
/// etable::EtableOptions class.  This, internally, will create an EtableOptions object,
/// initialized from the command line, set the etable_type of this object, and then
/// invoke the etable( EtableOptions ) method.
/// @details Threadsafe, since this calls the threadsafe etable( EtableOptions) function, above.
etable::EtableCAP
ScoringManager::etable( std::string const & etable_id ) const
{
	etable::EtableOptions default_options;
	default_options.etable_type = etable_id;
	return etable( default_options );
}

/// @brief Get an owning pointer to data used by the FA_ElecEnergy in beta_nov15 mode.
/// @details If the data have not been loaded, this loads the data (lazy loading).  Lazy loading is now threadsafe.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
core::scoring::elec::CPRepMapTypeCOP
ScoringManager::get_cp_rep_map_byname() const {
	boost::function< core::scoring::elec::CPRepMapTypeOP () > creator( boost::bind( &core::scoring::elec::read_cp_tables_from_db, std::string("scoring/score_functions/elec_cp_reps.dat") ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, cp_rep_map_byname_, SAFELY_PASS_MUTEX(cp_rep_map_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(cp_rep_map_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return cp_rep_map_byname_;
}

/// @brief Get a vector of owning pointers to data used by the AACompositionEnergy score term.
/// @details If this vector has not yet been populated, this loads the data from disk (lazy loading).
/// @note The lazy loading has been made threadsafe, as of the wee hours of 26 July 2017.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
utility::vector1< core::scoring::aa_composition_energy::AACompositionEnergySetupOP >
ScoringManager::get_cloned_aa_comp_setup_helpers(
	core::scoring::methods::EnergyMethodOptions const &options
) const {
	core::Size const n_setup_helpers( options.aa_composition_setup_file_count() );

	// Vector in which to return a clone of the data:
	utility::vector1< core::scoring::aa_composition_energy::AACompositionEnergySetupOP > return_vect;
	return_vect.reserve(n_setup_helpers);

	// Load the data if necessary (once, in a threadsafe manner), and populate the return vector:
	for ( core::Size i(1); i<=n_setup_helpers; ++i ) {
		boost::function< core::scoring::aa_composition_energy::AACompositionEnergySetupOP () > creator( boost::bind( &ScoringManager::create_aa_composition_energy_setup_instance, options.aa_composition_setup_file(i) ) );
		return_vect.push_back( ( utility::thread::safely_check_map_for_key_and_insert_if_absent(creator, SAFELY_PASS_MUTEX(aa_comp_mutex_), options.aa_composition_setup_file(i), aa_composition_setup_helpers_) )->clone() );
	}

	return return_vect;
}

/// @brief Get a particular MainchainScoreTable for the rama_prepro score term.
/// @details If this has not yet been populated, loads the data from disk (lazy loading)
/// in a threadsafe manner.
/// @note Each restype stores separate tables for general and pre-proline scoring.  The prepro_table parameter determines
/// whether we're loading the default scoring table or the version for residues occurring before a proline.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
core::chemical::mainchain_potential::MainchainScoreTableCOP
ScoringManager::get_rama_prepro_mainchain_torsion_potential(
	core::chemical::ResidueTypeCOP restype,
	bool const use_polycubic_interpolation,
	bool const prepro_table
) const {
	using namespace core::chemical::mainchain_potential;

	std::string const & mapname( restype->get_rama_prepro_mainchain_torsion_potential_name(prepro_table) );

	// First, check to see whether we need to create the potential.  This is now threadsafe:
#ifdef MULTI_THREADED
	utility::thread::ReadWriteMutex & mut( prepro_table ? rama_prepro_mainchain_potentials_mutex_ : rama_prepro_mainchain_potentials_beforeproline_mutex_ );
#endif
	std::map< std::string, MainchainScoreTableOP >::const_iterator iter, iter_end;
	{ //Scope for the read lock
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard lock( mut );
#endif
		if ( !prepro_table ) {
			iter = rama_prepro_mainchain_potentials_.find(mapname);
			iter_end = rama_prepro_mainchain_potentials_.end();
		} else {
			iter = rama_prepro_mainchain_potentials_beforeproline_.find(mapname);
			iter_end = rama_prepro_mainchain_potentials_beforeproline_.end();
		}
	} //Release read lock

	MainchainScoreTableCOP temptable, table_to_return;

	if ( iter == iter_end /*The map key was not found, meaning we need to load the data and add them to the map.*/ ) {
		std::string const mapfile( restype->get_rama_prepro_map_file_name(prepro_table) );
		if ( mapfile.empty() ) return MainchainScoreTableCOP(); //Return null pointer if there's no mapfile.
		utility::vector1< std::pair< std::string, MainchainScoreTableOP > > newtables;
		core::chemical::mainchain_potential::read_rama_map_file_shapovalov( mapfile, use_polycubic_interpolation, newtables ); //Read the file, and put results for ALL types in "newtables".
		bool mytype_found(false);
		for ( core::Size i=1, imax=newtables.size(); i<=imax; ++i ) {
			// Add the newly-created mainchain scoretable to the map.  This does this in a special way:
			// 1.  It checks to see if it's already in the map (with read guarded by a ReadWriteMutex).  It's possible that between the check above and the check here, another
			// thread already created the MainchainScoreTable in question and added it to the map.  In that case, we discard this one and
			// just return that one.
			// 2.  If it's not in the map, it gets a write lock for the ReadWriteMutex and adds it to the map (checking one more time, once the write lock is on, before doing
			// the addition).
			boost::function< MainchainScoreTableOP () > builder( boost::bind( &ScoringManager::create_mainchain_scoretable_instance, newtables[i].second ) ); //Note that this is a silly "builder" function that just returns the object passed to it.
			temptable = utility::thread::safely_check_map_for_key_and_insert_if_absent( builder, SAFELY_PASS_MUTEX( mut ), newtables[i].first, prepro_table ? rama_prepro_mainchain_potentials_beforeproline_ : rama_prepro_mainchain_potentials_ );
			if ( newtables[i].first.compare( mapname ) == 0 ) {
				mytype_found = true; //Found the type we're trying to load.
				table_to_return = temptable; //So let's store it so that we can return it.
			}
		}
		runtime_assert_string_msg( mytype_found, "Error in core::scoring::ScoringManager::get_rama_prepro_mainchain_torsion_potential().  Could not load RamaPrePro scoring table for " + mapname + " from file " + mapfile + "." );
	} else {
		// If the map key WAS found, then we just need to return it.
		table_to_return = iter->second;
	}

	return table_to_return;
}

/// @brief Test if there is an EnergyMethod class defined for a
/// given score type.
/// @details I THINK that this is threadsafe (VKM, 20 July 2017).
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

/// @brief Create an instance of the PairEPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
PairEPotentialOP
ScoringManager::create_pairE_potential_instance() {
	return PairEPotentialOP( new PairEPotential );
}

/// @brief Create an instance of the GenBornPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
GenBornPotentialOP
ScoringManager::create_genborn_instance() {
	return GenBornPotentialOP( new GenBornPotential );
}

/// @brief Create an instance of the HydroxylTorsionPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
HydroxylTorsionPotentialOP
ScoringManager::create_hxl_potential_instance() {
	return HydroxylTorsionPotentialOP( new HydroxylTorsionPotential );
}

/// @brief Create an instance of the VdWTinkerPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
VdWTinkerPotentialOP
ScoringManager::create_vdw_tinker_potential_instance() {
	return VdWTinkerPotentialOP( new VdWTinkerPotential );
}

/// @brief Create an instance of the HydroxylTorsionPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
MultipoleElecPotentialOP
ScoringManager::create_multipole_elec_instance(
	methods::EnergyMethodOptions const & options
) {
	return MultipoleElecPotentialOP( new MultipoleElecPotential(options) );
}

/// @brief Create an instance of the SASAPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
SASAPotentialOP
ScoringManager::create_sasa_potential_instance() {
	return SASAPotentialOP( new SASAPotential );
}


/// @brief Create an instance of the FactsPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
FACTSPotentialOP
ScoringManager::create_facts_potential_instance() {
	return FACTSPotentialOP( new FACTSPotential );
}

/// @brief Create an instance of the FactsPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
dna::DNA_BasePotentialOP
ScoringManager::create_dnabase_potential_instance() {
	return dna::DNA_BasePotentialOP( new dna::DNA_BasePotential );
}

/// @brief Create an instance of the Ramachandran object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
RamachandranOP
ScoringManager::create_rama_instance() {
	return RamachandranOP( new Ramachandran );
}

/// @brief Create an instance of the Ramachandran2B object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @author Vikram K. Mulligan (vmullig@uw.edu)
Ramachandran2BOP
ScoringManager::create_rama2b_instance() {
	return Ramachandran2BOP( new Ramachandran2B );
}

/// @brief Create an instance of the RamaPrePro object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @author Vikram K. Mulligan (vmullig@uw.edu)
RamaPreProOP
ScoringManager::create_ramapp_instance() {
	return RamaPreProOP( new RamaPrePro );
}

/// @brief Create an instance of the P_AA_ABEGO3 object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
P_AA_ABEGO3_OP
ScoringManager::create_p_aa_abego3_instance() {
	return P_AA_ABEGO3_OP( new P_AA_ABEGO3 );
}

/// @brief Create an instance of the P_AA object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
P_AAOP
ScoringManager::create_p_aa_instance() {
	return P_AAOP( new P_AA );
}

/// @brief Create an instance of the P_AA_ss object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
P_AA_ssOP
ScoringManager::create_p_aa_ss_instance() {
	return P_AA_ssOP( new P_AA_ss );
}

/// @brief Create an instance of the DNABFormPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
dna::DNABFormPotentialOP
ScoringManager::create_dna_bform_potential_instance() {
	return dna::DNABFormPotentialOP( new dna::DNABFormPotential );
}

/// @brief Create an instance of the DNATorsionPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
dna::DNATorsionPotentialOP
ScoringManager::create_dna_torsion_potential_instance() {
	return dna::DNATorsionPotentialOP( new dna::DNATorsionPotential );
}

/// @brief Create an instance of the OmegaTether object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
OmegaTetherOP
ScoringManager::create_omegatether_instance() {
	return OmegaTetherOP( new OmegaTether );
}

/// @brief Create an instance of the SmoothEnvPairPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
SmoothEnvPairPotentialOP
ScoringManager::create_smoothenvpair_instance() {
	return SmoothEnvPairPotentialOP( new SmoothEnvPairPotential );
}

/// @brief Create an instance of the CenRotEnvPairPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
CenRotEnvPairPotentialOP
ScoringManager::create_cenrotenvpair_instance() {
	return CenRotEnvPairPotentialOP( new CenRotEnvPairPotential );
}

/// @brief Create an instance of the CenHBPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
CenHBPotentialOP
ScoringManager::create_cenhbpotential_instance() {
	return CenHBPotentialOP( new CenHBPotential );
}

/// @brief Create an instance of the EnvPairPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
EnvPairPotentialOP
ScoringManager::create_envpairpotential_instance() {
	return EnvPairPotentialOP( new EnvPairPotential );
}

/// @brief Create an instance of the DNA_EnvPairPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
dna::DNA_EnvPairPotentialOP
ScoringManager::create_dna_envpairpotential_instance() {
	return dna::DNA_EnvPairPotentialOP( new dna::DNA_EnvPairPotential );
}

/// @brief Create an instance of the DNA_DihedralPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
dna::DNA_DihedralPotentialOP
ScoringManager::create_dna_dihedralpotential_instance() {
	return dna::DNA_DihedralPotentialOP( new dna::DNA_DihedralPotential );
}

/// @brief Create an instance of the SecondaryStructurePotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
SecondaryStructurePotentialOP
ScoringManager::create_secondarystructurepotential_instance() {
	return SecondaryStructurePotentialOP( new SecondaryStructurePotential );
}


/// @brief Create an instance of the AtomVDW object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
AtomVDWOP
ScoringManager::create_atomvdw_instance(
	std::string const & atom_type_set_name
) {
	return AtomVDWOP( new AtomVDW( atom_type_set_name ) );
}

/// @brief Create an instance of the RNA_AtomVDW object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
rna::RNA_AtomVDWOP
ScoringManager::create_rna_atomvdw_instance() {
	return rna::RNA_AtomVDWOP( new rna::RNA_AtomVDW );
}

/// @brief Create an instance of the DatabaseOccSolEne object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
geometric_solvation::DatabaseOccSolEneOP
ScoringManager::create_database_occsolene_instance(
	std::string const & atom_type_set_name,
	core::Real const & min_occ_energy
) {
	return geometric_solvation::DatabaseOccSolEneOP( new geometric_solvation::DatabaseOccSolEne( atom_type_set_name, min_occ_energy ) );
}

/// @brief Create an instance of the CarbonHBondPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
carbon_hbonds::CarbonHBondPotentialOP
ScoringManager::create_carbon_hbond_potential_instance() {
	return carbon_hbonds::CarbonHBondPotentialOP( new carbon_hbonds::CarbonHBondPotential );
}

/// @brief Create an instance of the RNA_SuitePotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
rna::RNA_SuitePotentialOP
ScoringManager::create_rna_suitepotential_instance(
	bool const & calculate_suiteness_bonus,
	std::string const & suiteness_bonus
) {
	return rna::RNA_SuitePotentialOP( new rna::RNA_SuitePotential( calculate_suiteness_bonus, suiteness_bonus ) );
}

/// @brief Create an instance of the SixDTransRotPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
loop_graph::evaluator::SixDTransRotPotentialOP
ScoringManager::create_sixdtransrotpotential_instance(
	std::string const & database_file
) {
	if ( utility::file::file_exists(  database_file ) )  {
		TR << "Reading in: " << database_file << std::endl;
		return loop_graph::evaluator::SixDTransRotPotentialOP( new loop_graph::evaluator::SixDTransRotPotential( database_file ) );
	}
	TR.Warning << "File " << database_file << " does not exist!" << std::endl;
	return 0; // save information that database file does not exist. (?)
}

/// @brief Create an instance of the RNA_LowResolutionPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
rna::RNA_LowResolutionPotentialOP
ScoringManager::create_rna_lowresolutionpotential_instance() {
	return rna::RNA_LowResolutionPotentialOP( new rna::RNA_LowResolutionPotential );
}

/// @brief Create an instance of the RNP_LowResPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
rna::RNP_LowResPotentialOP
ScoringManager::create_rnp_lowrespotential_instance() {
	return rna::RNP_LowResPotentialOP( new rna::RNP_LowResPotential );
}

/// @brief Create an instance of the RNP_LowResStackData object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
rna::RNP_LowResStackDataOP
ScoringManager::create_rnp_lowresstackdata_instance() {
	return rna::RNP_LowResStackDataOP( new rna::RNP_LowResStackData );
}

/// @brief Create an instance of the RNA_ChemicalShiftPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
rna::chemical_shift::RNA_ChemicalShiftPotentialOP
ScoringManager::create_rna_chemshiftpotential_instance() {
	return rna::chemical_shift::RNA_ChemicalShiftPotentialOP( new rna::chemical_shift::RNA_ChemicalShiftPotential );
}

/// @brief Create an instance of the RNA_DMS_Potential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
rna::data::RNA_DMS_PotentialOP
ScoringManager::create_rna_dms_potential_instance() {
	return rna::data::RNA_DMS_PotentialOP( new rna::data::RNA_DMS_Potential );
}

/// @brief Create an instance of the RNA_DMS_LowResolutionPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
rna::data::RNA_DMS_LowResolutionPotentialOP
ScoringManager::create_rna_dms_lowrespotential_instance() {
	return rna::data::RNA_DMS_LowResolutionPotentialOP( new rna::data::RNA_DMS_LowResolutionPotential );
}

/// @brief Create an instance of the DirectReadoutPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
dna::DirectReadoutPotentialOP
ScoringManager::create_dna_directreadoutpotential_instance() {
	return dna::DirectReadoutPotentialOP( new dna::DirectReadoutPotential );
}

/// @brief Create an instance of the MMLJLibrary object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
mm::MMLJLibraryOP
ScoringManager::create_mm_lj_library_instance() {
	return mm::MMLJLibraryOP( new mm::MMLJLibrary( chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD ) ) );
}

/// @brief Create an instance of the MMLJEnergyTable object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
mm::MMLJEnergyTableOP
ScoringManager::create_mm_lj_energy_table_instance() {
	return mm::MMLJEnergyTableOP( new mm::MMLJEnergyTable );
}

/// @brief Create an instance of the MMTorsionLibrary object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
mm::MMTorsionLibraryOP
ScoringManager::create_mm_torsion_library_instance() {
	return mm::MMTorsionLibraryOP(
		new mm::MMTorsionLibrary(
		basic::database::full_name( "chemical/mm_atom_type_sets/fa_standard/mm_torsion_params.txt" ),
		chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD )
		)
	);
}

/// @brief Create an instance of the MMBondAngleLibrary object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
mm::MMBondAngleLibraryOP
ScoringManager::create_mm_bondangle_library_instance() {
	return mm::MMBondAngleLibraryOP(
		new mm::MMBondAngleLibrary(
		basic::database::full_name( "chemical/mm_atom_type_sets/fa_standard/par_all27_prot_na.prm" ),
		chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD )
		)
	);
}

/// @brief Create an instance of the MMBondLengthLibrary object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
mm::MMBondLengthLibraryOP
ScoringManager::create_mm_bondlength_library_instance() {
	return mm::MMBondLengthLibraryOP(
		new mm::MMBondLengthLibrary(
		basic::database::full_name( "chemical/mm_atom_type_sets/fa_standard/par_all27_prot_na.prm" ),
		chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD )
		)
	);
}

/// @brief Create an instance of the NVlookup object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
nv::NVlookupOP
ScoringManager::create_nvlookup_instance() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	return nv::NVlookupOP( new nv::NVlookup(basic::database::full_name(option[score::NV_table]())) );
}

/// @brief Create an instance of the OrbitalsLookup object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
orbitals::OrbitalsLookupOP
ScoringManager::create_orbitals_lookup_instance() {
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

	return orbitals::OrbitalsLookupOP( new orbitals::OrbitalsLookup( DHO_energies, AOH_energies, AOO_orb_orb_energies, DOO_orb_orb_energies,ACO_energies ) );
}

/// @brief Create an instance of the DDPlookup object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
interface_::DDPlookupOP
ScoringManager::create_ddp_lookup_instance() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	return interface_::DDPlookupOP( new interface_::DDPlookup("scoring/score_functions/DDPscore/interface_ddp_score.txt") );
}

/// @brief Create an instance of the UnfoldedStatePotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
UnfoldedStatePotentialOP
ScoringManager::create_unfolded_state_potential_instance(std::string const &type) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( type == UNFOLDED_SPLIT_USER_DEFINED || option[ unfolded_state::unfolded_energies_file ].user() ) {
		TR << "Creating unfolded state potential using file: " <<  option[ unfolded_state::unfolded_energies_file ].value() << std::endl;
		return UnfoldedStatePotentialOP( new UnfoldedStatePotential( option[ unfolded_state::unfolded_energies_file ].value() ) );
	} else if ( type == UNFOLDED_SCORE12 ) {
		return UnfoldedStatePotentialOP( new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/unfolded/unfolded_state_residue_energies_score12" ) ) );
	} else if ( type == UNFOLDED_MM_STD ) {
		return UnfoldedStatePotentialOP( new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/unfolded/unfolded_state_residue_energies_mm_std" ) ) );
	} else if ( type == UNFOLDED_RNA ) {
		return UnfoldedStatePotentialOP( new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/unfolded/unfolded_state_residue_energies_rna" ) ) );  // This will later get more elaborated
	} else if ( type == UNFOLDED_SPLIT_TALARIS2013 ) {
		return UnfoldedStatePotentialOP( new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/split_unfolded/split_unfolded_one_body_talaris2013" ) ) ); //for the split unfolded energy one body term.
	} else if ( type == UNFOLDED_SPLIT_MM_STD ) {
		return UnfoldedStatePotentialOP( new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/split_unfolded/split_unfolded_one_body_mm_std" ) ) ); //for the split unfolded energy one body term.
	} else {
		utility_exit_with_message("unrecognized unfolded type: "+type );
	}
	return UnfoldedStatePotentialOP(nullptr); //To keep compiler happy
}

/// @brief Create an instance of the WaterAdductHBondPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
WaterAdductHBondPotentialOP
ScoringManager::create_water_adduct_hbond_potential_instance() {
	return WaterAdductHBondPotentialOP( new WaterAdductHBondPotential );
}

/// @brief Create an instance of the MembranePotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
MembranePotentialOP
ScoringManager::create_membrane_potential_instance() {
	return MembranePotentialOP( new MembranePotential );
}

/// @brief Create an instance of the MembraneData object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
membrane::MembraneDataOP
ScoringManager::create_membrane_data_instance() {
	return membrane::MembraneDataOP( new membrane::MembraneData );
}

/// @brief Create an instance of the Membrane_FAPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
Membrane_FAPotentialOP
ScoringManager::create_membrane_fa_potential_instance() {
	return Membrane_FAPotentialOP( new Membrane_FAPotential );
}

/// @brief Create an instance of the ProQPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
ProQPotentialOP
ScoringManager::create_proq_potential_instance() {
	return ProQPotentialOP( new ProQPotential );
}

/// @brief Create an instance of the PoissonBoltzmannPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
PoissonBoltzmannPotentialOP
ScoringManager::create_poisson_boltzmann_potential_instance() {
	return PoissonBoltzmannPotentialOP( new PoissonBoltzmannPotential );
}

/// @brief Create an instance of the SplitUnfoldedTwoBodyPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
SplitUnfoldedTwoBodyPotentialOP
ScoringManager::create_split_unfolded_2body_potential_instance(
	std::string const & label_type,
	std::string const & value_type,
	std::string const & score_func_type
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string database_path = "";
	std::string atom_label_type = "";

	if ( label_type == SPLIT_UNFOLDED_USER_DEFINED || value_type == SPLIT_UNFOLDED_USER_DEFINED || option[ unfolded_state::split_unfolded_energies_file ].user() ) {
		TR << "Creating split unfolded state potential using file: " <<  option[ unfolded_state::split_unfolded_energies_file ].value()
			<< " with an atom type of " << option[ unfolded_state::split_unfolded_energies_atom_type ].value() << std::endl;
		return SplitUnfoldedTwoBodyPotentialOP( new SplitUnfoldedTwoBodyPotential( option[ unfolded_state::split_unfolded_energies_file ].value(), option[ unfolded_state::split_unfolded_energies_atom_type ].value() ) );
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
		return SplitUnfoldedTwoBodyPotentialOP( new SplitUnfoldedTwoBodyPotential( basic::database::full_name( database_path ), atom_label_type ) );
	}
	TR << "Creating split unfolded state potential using file: " << basic::database::full_name( database_path ) << std::endl;
	return SplitUnfoldedTwoBodyPotentialOP( new SplitUnfoldedTwoBodyPotential( basic::database::full_name( database_path ) ) );
}

/// @brief Create an instance of the FullatomDisulfidePotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
disulfides::FullatomDisulfidePotentialOP
ScoringManager::create_fullatom_disulfide_potential_instance() {
	return disulfides::FullatomDisulfidePotentialOP( new disulfides::FullatomDisulfidePotential );
}

/// @brief Create an instance of the CentroidDisulfidePotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
disulfides::CentroidDisulfidePotentialOP
ScoringManager::create_centroid_disulfide_potential_instance() {
	return disulfides::CentroidDisulfidePotentialOP( new disulfides::CentroidDisulfidePotential );
}

/// @brief Create an instance of the DisulfideMatchingPotential object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
disulfides::DisulfideMatchingPotentialOP
ScoringManager::create_disulfide_matching_potential_instance() {
	return disulfides::DisulfideMatchingPotentialOP( new disulfides::DisulfideMatchingPotential );
}

/// @brief Create an instance of a MainchainScoreTable, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.  This one is kind of silly, since it just
/// returns a MainchainScoreTableOP that is passed in.  Still needed for threadsafe creation, though.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::chemical::mainchain_potential::MainchainScoreTableOP
ScoringManager::create_mainchain_scoretable_instance(
	core::chemical::mainchain_potential::MainchainScoreTableOP table_in
) {
	return table_in; //Silly, yes, but still needed.
}

/// @brief Create an instance of a MembEtable, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.  This one is kind of silly, since it just
/// returns a MembEtableOP that is passed in.  Still needed for threadsafe creation, though.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
etable::MembEtableOP ScoringManager::create_memb_etable_instance_silly( etable::MembEtableOP table_in ) { return table_in; }

/// @brief Create an instance of a MembEtable, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.  This is the non-silly version that builds
/// the object based on a string.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
etable::MembEtableOP
ScoringManager::create_memb_etable_instance(
	std::string const & table_id
) {
	if ( table_id != FA_STANDARD_DEFAULT ) {
		std::string msg = "unrecognized etable: "+table_id;
		utility_exit_with_message( msg );
	}
	return etable::MembEtableOP ( new etable::MembEtable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ), etable::EtableOptions() ) );
}

/// @brief Create an instance of an Etable, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.  This one is kind of silly, since it just
/// returns an EtableOP that is passed in.  Still needed for threadsafe creation, though.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
etable::EtableOP ScoringManager::create_etable_instance( etable::EtableOP table_in ) { return table_in; }

/// @brief Create an instance of an AACompositionEnergySetup object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::scoring::aa_composition_energy::AACompositionEnergySetupOP
ScoringManager::create_aa_composition_energy_setup_instance(
	std::string const &filename
) {
	using namespace core::scoring::aa_composition_energy;
	AACompositionEnergySetupOP op_to_return( new AACompositionEnergySetup );
	op_to_return->initialize_from_file( filename );
	return op_to_return;
}

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
