// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoringManager.cc
/// @brief  Scoring manager class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/scoring/ScoringManager.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/orbitals.OptionKeys.gen.hh>

#include <core/scoring/carbon_hbonds/CarbonHBondPotential.hh>
#include <core/scoring/PairEPotential.hh>
#include <core/scoring/EnvPairPotential.hh>
#include <core/scoring/CenRotEnvPairPotential.hh>
#include <core/scoring/SmoothEnvPairPotential.hh>
#include <core/scoring/CenHBPotential.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/Membrane_FAPotential.hh> //pba
#include <core/scoring/ProQPotential.hh>
#include <core/scoring/SecondaryStructurePotential.hh>
//#include <core/pack/dunbrack/RotamerLibrary.hh>
//#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
//#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/scoring/Rama2BOffset.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/Ramachandran2B.hh>
#include <core/scoring/OmegaTether.hh>
#include <core/scoring/GenBornPotential.hh>
#include <core/scoring/facts/FACTSPotential.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/rna/RNA_AtomVDW.hh>
#include <core/scoring/geometric_solvation/DatabaseOccSolEne.hh>
#include <core/scoring/dna/DNABFormPotential.hh>
#include <core/scoring/dna/DNATorsionPotential.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.hh>
#include <core/scoring/dna/DirectReadoutPotential.hh>
#include <core/scoring/P_AA.hh>
#include <core/scoring/WaterAdductHBondPotential.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>
#include <core/scoring/disulfides/CentroidDisulfidePotential.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/scoring/UnfoldedStatePotential.hh>
#include <core/scoring/PoissonBoltzmannPotential.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/memb_etable/MembEtable.hh>

#include <core/scoring/mm/MMTorsionLibrary.hh>
#include <core/scoring/mm/MMLJLibrary.hh>
#include <core/scoring/mm/MMLJEnergyTable.hh>
#include <core/scoring/mm/MMBondAngleLibrary.hh>
#include <core/scoring/mm/MMBondLengthLibrary.hh>

#include <core/scoring/nv/NVlookup.hh>
#include <core/scoring/orbitals/OrbitalsLookup.hh>
#include <core/scoring/interface/DDPlookup.hh>

#include <core/scoring/types.hh>
#include <core/scoring/ScoreType.hh>

#include <core/chemical/ChemicalManager.hh>

#include <basic/database/open.hh>


// Utility headers
#include <utility/string_util.hh>

#include <core/scoring/methods/EnergyMethod.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

namespace core {
namespace scoring {


///////////////////////////////////////////////////////////////////////////////
ScoringManager* ScoringManager::instance_( 0 );

///////////////////////////////////////////////////////////////////////////////
//singleton class
/// SAFE singleton initialization; static from within a function ensures
/// proper load-time behavior
ScoringManager * ScoringManager::get_instance()
{
	if ( instance_ == 0 ) {
		 instance_ = new ScoringManager();
	}
	return instance_;
}


ScoringManager::~ScoringManager() {}

///////////////////////////////////////////////////////////////////////////////
ScoringManager::ScoringManager() :
	pairE_potential_( 0 ),
	rama_( 0 ),
	rama2b_( 0 ),
	rama2bo_( 0 ),
	omega_( 0 ),
	env_pair_potential_( 0 ),
	smooth_env_pair_potential_( 0 ),
	cen_rot_pair_potential_( 0 ),
	cen_hb_potential_( 0 ),
	secondary_structure_potential_( 0 ),
	atom_vdw_(),
	rna_atom_vdw_( 0 ),
	occ_hbond_sol_database_( 0 ),
	dna_dr_potential_( 0 ),
	mm_lj_library_( 0 ),
	mm_lj_energy_table_( 0 ),
	mm_torsion_library_( 0 ),
	mm_bondangle_library_( 0 ),
	dnabform_( 0 ),
	dna_torsion_potential_( 0 ),
	DNA_base_potential_( 0 ),
	carbon_hbond_potential_( 0 ),
	rna_low_resolution_potential_( 0 ),
	rna_torsion_potential_( 0 ),
    rna_chemical_shift_potential_( 0 ),
	p_aa_( 0 ),
	water_adduct_hbond_potential_( 0 ),
	gen_born_potential_( 0 ),
	fa_disulfide_potential_( 0 ),
	cen_disulfide_potential_( 0 ),
	disulfide_matching_potential_( 0 ),
	membrane_potential_( 0 ),
  membrane_fapotential_( 0 ), //pba
	PB_potential_(0),
	unf_state_( 0 ),
	NV_lookup_table_(0),
	orbitals_lookup_table_( 0 ),
  DDP_lookup_table_(0),
	method_creator_map_( n_score_types, 0 )
{}

///////////////////////////////////////////////////////////////////////////////
PairEPotential const &
ScoringManager::get_PairEPotential() const
{
	if (pairE_potential_ == 0 )
	{
		pairE_potential_ = new PairEPotential();
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
		if ( method_creator_map_[ sts[ ii ] ] != 0 ) {
			utility_exit_with_message( "Cannot register a term to two different EnergyMethodCreators. Term " + utility::to_string( sts[ ii ] ) + " has already been registered!" );
		}
		method_creator_map_[ sts[ ii ] ] = creator;
	}
}


///////////////////////////////////////////////////////////////////////////////
dna::DNA_BasePotential const &
ScoringManager::get_DNA_BasePotential() const
{
	if (DNA_base_potential_ == 0 )
	{
		DNA_base_potential_ = new dna::DNA_BasePotential();
	}
	return *DNA_base_potential_;
}

///////////////////////////////////////////////////////////////////////////////
EnvPairPotential const &
ScoringManager::get_EnvPairPotential() const
{
	if (env_pair_potential_ == 0 )
	{
		env_pair_potential_ = new EnvPairPotential();
	}
	return *env_pair_potential_;
}

///////////////////////////////////////////////////////////////////////////////
SmoothEnvPairPotential const &
ScoringManager::get_SmoothEnvPairPotential() const
{
	if (smooth_env_pair_potential_ == 0 )
	{
		smooth_env_pair_potential_ = new SmoothEnvPairPotential();
	}
	return *smooth_env_pair_potential_;
}

///////////////////////////////////////////////////////////////////////////////
CenRotEnvPairPotential const &
ScoringManager::get_CenRotEnvPairPotential() const
{
	if (cen_rot_pair_potential_ == 0 )
	{
		cen_rot_pair_potential_ = new CenRotEnvPairPotential();
	}
	return *cen_rot_pair_potential_;
}

///////////////////////////////////////////////////////////////////////////////
CenHBPotential const &
ScoringManager::get_CenHBPotential() const
{
	if (cen_hb_potential_ == 0 )
	{
		cen_hb_potential_ = new CenHBPotential();
	}
	return *cen_hb_potential_;
}

///////////////////////////////////////////////////////////////////////////////
MembranePotential const &
ScoringManager::get_MembranePotential() const
{
	if (membrane_potential_ == 0 )
	{
		membrane_potential_ = new MembranePotential();
	}
	return *membrane_potential_;
}

///////////////////////////////////////////////////////////////////////////////
Membrane_FAPotential const &
ScoringManager::get_Membrane_FAPotential() const //pba
{
  if (membrane_fapotential_ == 0 )
  {
    membrane_fapotential_ = new Membrane_FAPotential();
  }
  return *membrane_fapotential_;
}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
ProQPotential const &
ScoringManager::get_ProQPotential() const
{
	if (ProQ_potential_ == 0 )
	{
		ProQ_potential_ = new ProQPotential();
	}
	return *ProQ_potential_;
}

///////////////////////////////////////////////////////////////////////////////


SecondaryStructurePotential const &
ScoringManager::get_SecondaryStructurePotential() const
{
	if (secondary_structure_potential_ == 0 )
	{
		secondary_structure_potential_ = new SecondaryStructurePotential();
	}
	return *secondary_structure_potential_;
}

///////////////////////////////////////////////////////////////////////////////
GenBornPotential const &
ScoringManager::get_GenBornPotential() const
{
	if (gen_born_potential_ == 0 )
	{
		gen_born_potential_ = new GenBornPotential();
	}
	return *gen_born_potential_;
}

///////////////////////////////////////////////////////////////////////////////
FACTSPotential const &
ScoringManager::get_FACTSPotential() const
{
	if (facts_potential_ == 0 )
	{
		facts_potential_ = new FACTSPotential();
	}
	return *facts_potential_;
}

///////////////////////////////////////////////////////////////////////////////
PoissonBoltzmannPotential const &
ScoringManager::get_PoissonBoltzmannPotential() const
{
	if ( PB_potential_ == 0 )
	{
		PB_potential_ =  new PoissonBoltzmannPotential;
	}
	return *PB_potential_;
}

///////////////////////////////////////////////////////////////////////////////
AtomVDW const &
ScoringManager::get_AtomVDW( std::string const & atom_type_set_name ) const
{
	if ( atom_vdw_.count( atom_type_set_name ) == 0 ) {
		atom_vdw_[ atom_type_set_name ] = new AtomVDW( atom_type_set_name );
	}
	return * ( atom_vdw_[ atom_type_set_name ] );
}

///////////////////////////////////////////////////////////////////////////////
carbon_hbonds::CarbonHBondPotential const &
ScoringManager::get_CarbonHBondPotential() const
{
	if (carbon_hbond_potential_ == 0 )
	{
		carbon_hbond_potential_ = new carbon_hbonds::CarbonHBondPotential();
	}
	return *carbon_hbond_potential_;
}

///////////////////////////////////////////////////////////////////////////////
rna::RNA_AtomVDW const &
ScoringManager::get_RNA_AtomVDW() const
{
	if (rna_atom_vdw_ == 0 )
	{
		rna_atom_vdw_ = new rna::RNA_AtomVDW();
	}
	return *rna_atom_vdw_;
}

///////////////////////////////////////////////////////////////////////////////
geometric_solvation::DatabaseOccSolEne const &
ScoringManager::get_DatabaseOccSolEne( std::string const & atom_type_set_name, Real const & min_occ_energy ) const
{
	if (occ_hbond_sol_database_ == 0 )
	{
		occ_hbond_sol_database_ = new geometric_solvation::DatabaseOccSolEne( atom_type_set_name, min_occ_energy );
	}
	return *occ_hbond_sol_database_;
}

///////////////////////////////////////////////////////////////////////////////
rna::RNA_LowResolutionPotential const &
ScoringManager::get_RNA_LowResolutionPotential() const
{
	if (rna_low_resolution_potential_ == 0 )
	{
		rna_low_resolution_potential_ = new rna::RNA_LowResolutionPotential();
	}
	return *rna_low_resolution_potential_;
}

///////////////////////////////////////////////////////////////////////////////
rna::RNA_TorsionPotential const &
ScoringManager::get_RNA_TorsionPotential() const
{
	if (rna_torsion_potential_ == 0 )
	{
		rna_torsion_potential_ = new rna::RNA_TorsionPotential();
	}
	return *rna_torsion_potential_;
}

///////////////////////////////////////////////////////////////////////////////
rna::chemical_shift::RNA_ChemicalShiftPotential const &
ScoringManager::get_RNA_ChemicalShiftPotential() const
{
	if (rna_chemical_shift_potential_ == 0 )
	{
		rna_chemical_shift_potential_= new rna::chemical_shift::RNA_ChemicalShiftPotential();
	}
	return *rna_chemical_shift_potential_;
}

///////////////////////////////////////////////////////////////////////////////
dna::DirectReadoutPotential const &
ScoringManager::get_DirectReadoutPotential() const
{
	if (dna_dr_potential_ == 0 )
	{
		dna_dr_potential_ = new dna::DirectReadoutPotential();
	}
	return *dna_dr_potential_;
}

P_AA const &
ScoringManager::get_P_AA() const
{
	if ( p_aa_ == 0 ) {
		p_aa_ = new P_AA;
	}
	return *p_aa_;
}

WaterAdductHBondPotential const &
ScoringManager::get_WaterAdductHBondPotential() const
{
	if ( water_adduct_hbond_potential_ == 0 ) {
		water_adduct_hbond_potential_ = new WaterAdductHBondPotential;
	}
	return *water_adduct_hbond_potential_;
}


/////////////////////////////////////
//////////////////////////////
////////////
//ScoringManager::RotamerLibrary &
//ScoringManager::get_RotamerLibrary() const
//{
//	if (rotamer_Library_ == 0 )
//	{
//		rotamer_Library_ = new RotamerLibrary();
//		read_dunbrack_library( *rotamer_Library_ );
//	}
//	return *rotamer_Library_;
//}

///////////////////////////////////////////////////////////////////////////////
Ramachandran const &
ScoringManager::get_Ramachandran() const
{
	if ( rama_ == 0 )
	{
		rama_ =  new Ramachandran;
	}
	return *rama_;
}

///////////////////////////////////////////////////////////////////////////////
Ramachandran2B const &
ScoringManager::get_Ramachandran2B() const
{
	if ( rama2b_ == 0 )
	{
		rama2b_ =  new Ramachandran2B;
	}
	return *rama2b_;
}

///////////////////////////////////////////////////////////////////////////////
Rama2BOffset const &
ScoringManager::get_Rama2BOffset() const
{
	if ( rama2bo_ == 0 )
	{
		rama2bo_ =  new Rama2BOffset;
	}
	return *rama2bo_;
}

///////////////////////////////////////////////////////////////////////////////
OmegaTether const &
ScoringManager::get_OmegaTether() const
{
	if ( omega_ == 0 )
	{
		omega_ =  new OmegaTether();
	}
	return *omega_;
}

///////////////////////////////////////////////////////////////////////////////
dna::DNABFormPotential const &
ScoringManager::get_DNABFormPotential() const
{
	if( dnabform_ == 0 )
	{
		dnabform_ =  new dna::DNABFormPotential;
	}
	return *dnabform_;
}

///////////////////////////////////////////////////////////////////////////////

dna::DNATorsionPotential const &
ScoringManager::get_DNATorsionPotential() const
{
	if (dna_torsion_potential_ == 0 )
	{
		dna_torsion_potential_ = new dna::DNATorsionPotential();
	}
	return *dna_torsion_potential_;
}

///////////////////////////////////////////////////////////////////////////////
core::scoring::mm::MMTorsionLibrary const &
ScoringManager::get_MMTorsionLibrary() const
{
	if ( mm_torsion_library_ == 0 )
		{
			mm_torsion_library_ = new mm::MMTorsionLibrary
				( basic::database::full_name( "chemical/mm_atom_type_sets/fa_standard/mm_torsion_params.txt" ),
					chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD ) );
		}
	return *mm_torsion_library_;
}
///////////////////////////////////////////////////////////////////////////////

core::scoring::mm::MMLJLibrary const &
ScoringManager::get_MMLJLibrary() const
{
	if ( mm_lj_library_ == 0 )
		{
			mm_lj_library_ = new mm::MMLJLibrary
				( chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD ) );
		}
	return *mm_lj_library_;
}

///////////////////////////////////////////////////////////////////////////////

core::scoring::mm::MMLJEnergyTable const &
ScoringManager::get_MMLJEnergyTable () const
{
	if ( mm_lj_energy_table_ == 0 )
		{
			mm_lj_energy_table_ = new mm::MMLJEnergyTable();
		}
	return *mm_lj_energy_table_;
}

///////////////////////////////////////////////////////////////////////////////

disulfides::FullatomDisulfidePotential &
ScoringManager::get_FullatomDisulfidePotential() const
{
	if ( fa_disulfide_potential_ == 0 ) {
		fa_disulfide_potential_ = new disulfides::FullatomDisulfidePotential;
	}
	return *fa_disulfide_potential_;
}

disulfides::CentroidDisulfidePotential &
ScoringManager::get_CentroidDisulfidePotential() const
{
	if ( cen_disulfide_potential_ == 0 ) {
		cen_disulfide_potential_ = new disulfides::CentroidDisulfidePotential;
	}
	return *cen_disulfide_potential_;
}

disulfides::DisulfideMatchingPotential &
ScoringManager::get_DisulfideMatchingPotential() const
{
	if ( disulfide_matching_potential_ == 0 ) {
		disulfide_matching_potential_ = new disulfides::DisulfideMatchingPotential;
	}
	return *disulfide_matching_potential_;
}

///////////////////////////////////////////////////////////////////////////////


nv::NVlookup const &
ScoringManager::get_NVLookupTable() const
{
	if(NV_lookup_table_ == 0)
	{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		NV_lookup_table_ = new nv::NVlookup(basic::database::full_name(option[score::NV_table]()));
	}
	return *NV_lookup_table_;

}

///////////////////////////////////////////////////////////////////////////////
	orbitals::OrbitalsLookup const &
	ScoringManager::get_OrbitalsLookupTable() const
	{
		if(orbitals_lookup_table_ == 0){
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

			orbitals_lookup_table_ = new orbitals::OrbitalsLookup(
					DHO_energies, AOH_energies, AOO_orb_orb_energies, DOO_orb_orb_energies,ACO_energies );

		}
		return *orbitals_lookup_table_;
	}




///////////////////////////////////////////////////////////////////////////////


interface::DDPlookup const &
ScoringManager::get_DDPLookupTable() const
{
	if(DDP_lookup_table_ == 0)
	{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		DDP_lookup_table_ = new interface::DDPlookup("scoring/score_functions/DDPscore/interface_ddp_score.txt");
	}
	return *DDP_lookup_table_;
}

///////////////////////////////////////////////////////////////////////////////
core::scoring::mm::MMBondAngleLibrary const &
ScoringManager::get_MMBondAngleLibrary() const
{
	if ( mm_bondangle_library_ == 0 )
		{
			mm_bondangle_library_ = new mm::MMBondAngleLibrary
				( basic::database::full_name( "chemical/mm_atom_type_sets/fa_standard/par_all27_prot_na.prm" ),
					chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD ) );
		}
	return *mm_bondangle_library_;
}

///////////////////////////////////////////////////////////////////////////////
core::scoring::mm::MMBondLengthLibrary const &
ScoringManager::get_MMBondLengthLibrary() const
{
	if ( mm_bondlength_library_ == 0 )
		{
			mm_bondlength_library_ = new mm::MMBondLengthLibrary
				( basic::database::full_name( "chemical/mm_atom_type_sets/fa_standard/par_all27_prot_na.prm" ),
					chemical::ChemicalManager::get_instance()->mm_atom_type_set( chemical::FA_STANDARD ) );
		}
	return *mm_bondlength_library_;
}

///////////////////////////////////////////////////////////////////////////////
UnfoldedStatePotential const &
ScoringManager::get_UnfoldedStatePotential( std::string const & type ) const
{
	if ( unf_state_ == 0 ) {
		if ( type == UNFOLDED_SCORE12 ) {
			unf_state_ = new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/unfolded/unfolded_state_residue_energies_score12" ) );
		} else if ( type == UNFOLDED_MM_STD ) {
			unf_state_ = new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/unfolded/unfolded_state_residue_energies_mm_std" ) );
		} else if ( type == UNFOLDED_RNA ) {
			unf_state_ = new UnfoldedStatePotential( basic::database::full_name( "scoring/score_functions/unfolded/unfolded_state_residue_energies_rna" ) );  // This will later get more elaborated
		} else {
			utility_exit_with_message("unrecognized unfolded type: "+type );
		}
	}
	return *unf_state_;
}

///////////////////////////////////////////////////////////////////////////////
void
ScoringManager::add_etable( std::string const & name, etable::EtableOP etable )
{
	assert( etables_.count(name) == 0 );
	etables_[ name ] = etable;
}

//XRW_B_T1
/*
///////////////////////////////////////////////////////////////////////////////
void
ScoringManager::add_coarse_etable( std::string const &name, coarse::CoarseEtableOP etable )
{
	assert( coarse_etables_.count(name) == 0);
	coarse_etables_ [name] = etable;
}
*/
 //XRW_E_T1

///////////////////////////////////////////////////p////////////////////////////
void
ScoringManager::add_memb_etable( std::string const & name, etable::MembEtableOP etable ) //pba
{
  assert( memb_etables_.count(name) == 0 );
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
      MembEtableOP etable_ptr
        ( new MembEtable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
                      EtableOptions() ) );
      memb_etables_[ table_id ] = etable_ptr;
    } else {
      std::string msg = "unrecognized etable: "+table_id;
      utility_exit_with_message( msg );
    }
  }
  return (memb_etables_.find( table_id )->second)();
}

///////////////////////////////////////////////////////////////////////////////
etable::EtableCAP
ScoringManager::etable( std::string const & table_id ) const
{
	using namespace etable;

	if ( etables_.find( table_id ) == etables_.end() ) {
		// try to build if possible
		if ( table_id == FA_STANDARD_DEFAULT ) {
			EtableOP etable_ptr
				( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
											EtableOptions() ) );
			etables_[ table_id ] = etable_ptr;
		} else if ( table_id == FA_STANDARD_SOFT ) {
			// soft rep etable: modified radii and also change to lj_switch_dis2sigma
			EtableOptions options;
			options.lj_switch_dis2sigma = 0.91;
			EtableOP etable_ptr
				( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
											options, "SOFT" ) );
			etables_[ table_id ] = etable_ptr;

		} else if ( table_id.substr(0, FA_STANDARD_DEFAULT.size() + 1 ) == FA_STANDARD_DEFAULT+"_" ) {
			// note we check for soft rep 1st since that would match this as well -- confusing??
			// apply a modification of the radii/wdepths
			std::string const alternate_parameters( table_id.substr( FA_STANDARD_DEFAULT.size() + 1 ) );
			EtableOP etable_ptr
				( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
											EtableOptions(), alternate_parameters ) );
			etables_[ table_id ] = etable_ptr;

		} else {
			std::cout << "1:" << table_id.substr(0, FA_STANDARD_DEFAULT.size() + 1 ) << '\n' <<
				"2:" << FA_STANDARD_DEFAULT+"_" << std::endl;
			utility_exit_with_message("unrecognized etable: "+table_id );
		}
	}
	return (etables_.find( table_id )->second)();
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//XRW_B_T1
/*
coarse::CoarseEtableCAP
ScoringManager::coarse_etable( std::string const & table_id ) const
{
	using namespace etable;
	if ( coarse_etables_.find( table_id ) == coarse_etables_.end() ) {
		// try to build if possible
		if ( table_id == chemical::COARSE_TWO_BEAD ) {
			coarse::CoarseEtableOP etable_ptr
				( new coarse::CoarseEtable(
							chemical::ChemicalManager::get_instance()->
							atom_type_set( chemical::COARSE_TWO_BEAD ), chemical::COARSE_TWO_BEAD )
				);
			coarse_etables_[ table_id ] = etable_ptr;
		} else {
			std::string msg = "unrecognized etable: "+table_id;
			utility_exit_with_message( msg );
		}
	}
	return (coarse_etables_.find( table_id )->second)();
}
*/
 //XRW_E_T1

///////////////////////////////////////////////////////////////////////////////
/// alot of this was pulled from RotamerLibrary.cc
/*pack::dunbrack::SingleResidueRotamerLibraryCAP
ScoringManager::get_NCAARotamerLibrary( chemical::ResidueType const & rsd_type )
{
	using namespace pack::dunbrack;

	// get some info about amino acid type
	std::string aa_name3( rsd_type.name3() );
	Size n_rotlib_chi( rsd_type.nchi() - rsd_type.n_proton_chi() );
	chemical::AA aan( rsd_type.aa() );
	bool dun02( true );

	if ( ncaa_rotlibs_.find( aa_name3 ) == ncaa_rotlibs_.end() ) {

		// create izstream from path
		std::string dir_name = basic::database::full_name( "/rotamer/ncaa_rotlibs/" );
		std::string file_name = rsd_type.get_ncaa_rotlib_path();
		utility::io::izstream rotlib_in( dir_name + file_name );
		std::cout << "Reading in rot lib " << dir_name + file_name << "...";

		// get an instance of RotamericSingleResidueDunbrackLibrary, but need a RotamerLibrary to do it
		// this means that when ever you read in the NCAA libraries you will also read in the Dunbrack libraries
		// this may need to be a pointer to the full type and not just a SRRLOP
		std::string empty_string("");

		// this comes almost directally from RotmerLibrary.cc::create_rotameric_dunlib()
		SingleResidueRotamerLibraryOP ncaa_rotlib;
		RotamericSingleResidueDunbrackLibrary< ONE > * r1;
 		RotamericSingleResidueDunbrackLibrary< TWO > * r2;
 		RotamericSingleResidueDunbrackLibrary< THREE > * r3;
 		RotamericSingleResidueDunbrackLibrary< FOUR > * r4;

		switch ( n_rotlib_chi ) {
		case 1:
			r1 = new RotamericSingleResidueDunbrackLibrary< ONE >( aan, dun02 );
			r1->set_n_chi_bins( rsd_type.get_ncaa_rotlib_n_bin_per_rot() );
			r1->read_from_file( rotlib_in, false );
			ncaa_rotlib = r1; break;
		case 2:
			r2 = new RotamericSingleResidueDunbrackLibrary< TWO >( aan, dun02 );
			r2->set_n_chi_bins( rsd_type.get_ncaa_rotlib_n_bin_per_rot() );
			r2->read_from_file( rotlib_in, false );
			ncaa_rotlib = r2; break;
		case 3:
			r3 = new RotamericSingleResidueDunbrackLibrary< THREE >( aan, dun02 );
			r3->set_n_chi_bins( rsd_type.get_ncaa_rotlib_n_bin_per_rot() );
			r3->read_from_file( rotlib_in, false );
			ncaa_rotlib = r3; break;
		case 4:
			r4 = new RotamericSingleResidueDunbrackLibrary< FOUR >( aan, dun02 );
			r4->set_n_chi_bins( rsd_type.get_ncaa_rotlib_n_bin_per_rot() );
			r4->read_from_file( rotlib_in, false );
			ncaa_rotlib = r4; break;
		default:
			utility_exit_with_message( "ERROR: too many chi angles desired for ncaa library: " + n_rotlib_chi );
			break;
		}

		// add new rotamer library to map
		ncaa_rotlibs_[ aa_name3 ] = ncaa_rotlib;
		std::cout << "done!" << std::endl;
	}
	return ( ncaa_rotlibs_.find( aa_name3 )->second)();
}*/

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

	if ( method_creator_map_[score_type] == 0 ) {
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

	if ( score_type == python ) return 0; /// python special case; this could now be changed...

	if ( method_creator_map_[ score_type ] == 0 ) {
		throw utility::excn::EXCN_Msg_Exception( "Requested ScoreType '" + utility::to_string( score_type ) + "' does not have a registered EnergyMethodCreator." );
	}

	return method_creator_map_[ score_type ]->create_energy_method( options );

	/* OLD WAY
	using namespace methods;
	switch( score_type ) {
	case fa_atr:
	case fa_rep:
	case fa_sol:
	case fa_intra_atr:
	case fa_intra_rep:
	case fa_intra_sol:
		return new etable::EtableEnergy( *etable( options.etable_type() ), options );
	case lk_hack:
		return new methods::LK_hack( *etable( options.etable_type() ));
	case lk_costheta:
	case lk_polar:
	case lk_nonpolar:
		return new methods::LK_CosThetaEnergy( *etable( options.etable_type() ));
  case fa_mbenv:  //pba
    return new methods::Fa_MbenvEnergy( *memb_etable( options.etable_type() ));
  case fa_mbsolv: //pba
    return new methods::Fa_MbsolvEnergy( *etable( options.etable_type() ), *memb_etable( options.etable_type() ));
	case coarse_fa_atr:
	case coarse_fa_rep:
	case coarse_fa_sol:
	case coarse_beadlj:
		return new etable::CoarseEtableEnergy( *etable( options.etable_type() ), options,
			coarse_etable( chemical::COARSE_TWO_BEAD ) );
	case hbond_lr_bb:
	case hbond_sr_bb:
	case hbond_bb_sc:
	case hbond_sr_bb_sc:
	case hbond_lr_bb_sc:
	case hbond_sc:
		return new hbonds::HBondEnergy( options.hbond_options() );
	case ch_bond:
	case ch_bond_sc_sc:
	case ch_bond_bb_sc:
	case ch_bond_bb_bb:
		return new hbonds::CarbonHBondEnergy();
	case geom_sol:
		return new hbonds::GeometricSolEnergy( options );
	case occ_sol_fitted:
		return new geometric_solvation::OccludedHbondSolEnergy( options );
	case occ_sol_fitted_onebody:
		return new geometric_solvation::OccludedHbondSolEnergy_onebody( options );
	case occ_sol_exact:
		return new geometric_solvation::ExactOccludedHbondSolEnergy(
			basic::options::option[basic::options::OptionKeys::score::exact_occ_pairwise],
			basic::options::option[basic::options::OptionKeys::score::exact_occ_split_between_res],
			! basic::options::option[ basic::options::OptionKeys::score::exact_occ_self_res_no_occ ],
			basic::options::option[basic::options::OptionKeys::score::exact_occ_radius_scaling] );
	case dslf_ss_dst:
	case dslf_cs_ang:
	case dslf_ss_dih:
	case dslf_ca_dih:
	case dslf_cbs_ds:
		return new disulfides::FullatomDisulfideEnergy( get_FullatomDisulfidePotential() );
	case dslfc_cen_dst:
	case dslfc_cb_dst:
	case dslfc_ang:
	case dslfc_cb_dih:
	case dslfc_bb_dih:
		return new disulfides::CentroidDisulfideEnergy( get_CentroidDisulfidePotential() );
	case neigh_count:
	case neigh_vect:
	case neigh_vect_raw:
		return new nv::NVscore();
	case symE_bonus:
		return new symDesign::symE();
	case gb_elec:
		return new methods::GenBornEnergy( options );
	case fa_pair:
	case fa_pair_aro_aro:
	case fa_pair_aro_pol:
	case fa_pair_pol_pol:
		return new PairEnergy;
	case fa_dun:
		return new DunbrackEnergy;
	case dna_chi:
		return new dna::DNAChiEnergy;
	case p_aa_pp:
		return new P_AA_pp_Energy;
	case pro_close:
		return new ProClosureEnergy;
	case h2o_intra:
		return new WaterAdductIntraEnergy;
	case h2o_hbond:
		return new WaterAdductHBondEnergy;
	case ref:
		if ( options.has_method_weights( ref ) ) {
			return new ReferenceEnergy( options.method_weights( ref ) );
		} else {
			return new ReferenceEnergy;
		}
	case rama:
		return new RamachandranEnergy;
	case rama2b:
		return new RamachandranEnergy2B;
	case omega:
		return new OmegaTetherEnergy;
	case chainbreak:
		return new ChainbreakEnergy;
	case linear_chainbreak:
	case overlap_chainbreak:
		return new LinearChainbreakEnergy;
	case distance_chainbreak:
		return new DistanceChainbreakEnergy;
	case envsmooth:
		return new EnvSmoothEnergy;
	case e_pH:
		return new pHEnergy;
	case env:
	case cbeta:
		return new EnvEnergy;
	case pair:
	case cenpack:
		return new CenPairEnergy;
	case Menv:
			return new MembraneEnvEnergy;
	case Menv_non_helix:
	case Menv_termini:
	case Menv_tm_proj:
		return new MembraneEnvPenalties;
	case Mcbeta:
			return new MembraneCbetaEnergy;
	case Mpair:
		return new MembraneCenPairEnergy;
	case Mlipo:
		return new MembraneLipo;
	case interchain_env:
	case interchain_contact:
		return new InterchainEnvEnergy;
	case interchain_pair:
	case interchain_vdw:
		return new InterchainPairEnergy;
	case hs_pair:
	case ss_pair:
	case rsigma:
	case sheet:
		return new SecondaryStructureEnergy;
	case rdc:
		return new ResidualDipolarCouplingEnergy;
	case rdc_rohl:
		return new ResidualDipolarCouplingEnergy_Rohl;
	case holes:
	case holes_decoy:
	case holes_resl:
	case holes_min:
		return new packing::HolesEnergy;
	case dab_sasa:
	case dab_sev:
		return new packing::SurfVolEnergy;
	case vdw:
		return new VDW_Energy( options );
	case hybrid_vdw:
		return new HybridVDW_Energy();
	case rna_rg:
		return new rna::RG_Energy_RNA;
	case rna_vdw:
		return new rna::RNA_VDW_Energy;
	case rna_base_pair:
	case rna_base_axis:
	case rna_base_stagger:
	case rna_base_stack:
	case rna_base_stack_axis:
	case rna_base_pair_pairwise:
	case rna_base_axis_pairwise:
	case rna_base_stagger_pairwise:
	case rna_base_stack_pairwise:
	case rna_base_stack_axis_pairwise:
	case rna_base_backbone:
	case rna_backbone_backbone:
	case rna_repulsive:
		return new rna::RNA_PairwiseLowResolutionEnergy;
	case fa_stack:
		return new rna::RNA_FullAtomStackingEnergy;
	case rna_torsion:
	case rna_sugar_close:
		return new rna::RNA_TorsionEnergy;
	case rna_fa_atr_base:
	case rna_fa_rep_base:
		return new rna::RNA_LJ_BaseEnergy( *etable( options.etable_type() ) );
  case dna_bb_torsion:
  case dna_sugar_close:
  case dna_base_distance:
    return new dna::DNATorsionEnergy;
	case rg:
		return new RG_Energy_Fast;
	case co:
		return new ContactOrderEnergy;
	case rms:
		return new RMS_Energy;
	case mm_twist:
		return new MMTorsionEnergy;
	case mm_bend:
		return new MMBondAngleEnergy( options );
//	case csd_torsion:
//		return new CSD_TorsionEnergy();
	case fa_cust_pair_dist:
		return new custom_pair_distance::FullatomCustomPairDistanceEnergy;
	case python:
		return NULL;
	case fa_elec:
	case fa_elec_bb_bb:
	case fa_elec_bb_sc:
	case fa_elec_sc_sc:
		return new elec::FA_ElecEnergy( options );
	case fa_elec_rna_phos_phos:
	case fa_elec_rna_phos_sugr:
	case fa_elec_rna_phos_base:
	case fa_elec_rna_sugr_sugr:
	case fa_elec_rna_sugr_base:
	case fa_elec_rna_base_base:
		return new elec::RNA_FA_ElecEnergy( options );
	case dna_bp:
	case dna_bs:
		return new DNA_BaseEnergy;
	case dna_dr:
		return new DirectReadoutEnergy;
	case atom_pair_constraint:
	case angle_constraint:
	case dihedral_constraint:
	case constant_constraint:
	case coordinate_constraint:
	case dof_constraint:
	case res_type_constraint:
	case backbone_stub_constraint:
  case backbone_stub_linear_constraint:
	case big_bin_constraint:
	case rna_bond_geometry:
		return new constraints::ConstraintsEnergy;
	case suck:
		return new SuckerEnergy;
	case gauss:
		return new GaussianOverlapEnergy;
	case pack_stat:
		return new PackStatEnergy;
	case surface:
		return new SurfaceEnergy;
	case p_aa:
		return new P_AA_Energy;
	case unfolded:
		if ( options.has_method_weights( unfolded ) ) {
			return new UnfoldedStateEnergy( options.method_weights( unfolded ) );
		} else {
			return new UnfoldedStateEnergy;
		}
	case elec_dens_window:
		return new ElecDensEnergy;
	case elec_dens_whole_structure_ca:
		return new ElecDensCenEnergy;
	case elec_dens_whole_structure_allatom:
		return new ElecDensAllAtomCenEnergy;
	case peptide_bond:
		return new PeptideBondEnergy;
	case pcs:
		return new pcs::PCS_Energy;
	default:
		utility_exit_with_message(
			"no energymethod for this score_type " + name_from_score_type( score_type )
		);
	}
	return 0;
	*/
}

/// global etable_id
std::string const FA_STANDARD_DEFAULT( "FA_STANDARD_DEFAULT" );
std::string const FA_STANDARD_SOFT   ( "FA_STANDARD_SOFT" );

std::string const UNFOLDED_SCORE12( "UNFOLDED_SCORE12" );
std::string const UNFOLDED_MM_STD( "UNFOLDED_MM_STD" );
std::string const UNFOLDED_RNA( "UNFOLDED_RNA" ); // This will later get more elaborated

} // namespace core
} // namespace scoring

