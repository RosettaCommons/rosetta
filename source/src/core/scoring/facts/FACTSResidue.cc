// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// @file:   core/scoring/facts/FACTSPotential.cc
// @brief:  The definitions of 3 classes of the FACTS algorithm resides here (see devel/khorvash/FACTSPotential.hh
// @author: Hahnbeom Park

// Unit headers
#include <core/scoring/facts/FACTSPotential.fwd.hh>
#include <core/scoring/facts/FACTSResidue.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/residue_io.hh>

#include <basic/prof.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <math.h>
#include <stdio.h>
#include <utility/assert.hh>
#include <utility/assert.hh>

static thread_local basic::Tracer TR( "core.scoring.FACTSPotential" );

# define Math_PI 3.14159265358979323846

using namespace std;

namespace core {
namespace scoring {

/**************************************************************************************************/
/*                                                                                                */
/*    @brief: The FACTSRsdTypeInfo class provides all the constants and parameters for given      */
/*            residue type                                                                        */
/*                                                                                                */
/**************************************************************************************************/
void FACTSRsdTypeInfo::create_info( chemical::ResidueType const & rsd )

{
  initialize_parameters( rsd );
  initialize_intrascale( rsd );
}

// This function initializes native parameters that are used for empirical function calculations
void FACTSRsdTypeInfo::initialize_parameters( chemical::ResidueType const & rsd ){

  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  // First define natoms
  natoms_ = rsd.natoms();

  //grab data from "database/chemical/atom_type_sets/fa_standard/extras/facts_params.txt
  Size const FACTS_RADIUS_INDEX( rsd.atom_type_set().extra_parameter_index( "FACTS_RADIUS" ) );
  Size const FACTS_CUT_INDEX   ( rsd.atom_type_set().extra_parameter_index( "FACTS_CUT"    ) );
  Size const FACTS_B1_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_B1"     ) );
  Size const FACTS_B2_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_B2"     ) );
  Size const FACTS_D1_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_D1"     ) );
  Size const FACTS_D2_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_D2"     ) );
  Size const FACTS_A0_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_A0"     ) );
  Size const FACTS_A1_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_A1"     ) );
  Size const FACTS_A2_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_A2"     ) );
  Size const FACTS_A3_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_A3"     ) );
  Size const FACTS_C0_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_C0"     ) );
  Size const FACTS_C1_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_C1"     ) );
  Size const FACTS_C2_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_C2"     ) );
  Size const FACTS_C3_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_C3"     ) );

  Size FACTS_ALPHA_INDEX;
  Size const asp_patch( option[ score::facts_asp_patch ]() );
  if ( asp_patch == 2 ){
    FACTS_ALPHA_INDEX = rsd.atom_type_set().extra_parameter_index( "FACTS_ALPHA2"  );
  } else if ( asp_patch == 3 ){
    FACTS_ALPHA_INDEX = rsd.atom_type_set().extra_parameter_index( "FACTS_ALPHA3"  );
  } else {
    FACTS_ALPHA_INDEX = rsd.atom_type_set().extra_parameter_index( "FACTS_ALPHA"  );
  }

  // Initialize array sizes & put in parameters
  q_.resize( natoms() );
  COradius2_.resize( natoms() );
  volume_.resize( natoms() );
  alpha_.resize( natoms() );
  a0_.resize( natoms() );
  a1_.resize( natoms() );
  a2_.resize( natoms() );
  a3_.resize( natoms() );
  b1_.resize( natoms() );
  b2_.resize( natoms() );
  c0_.resize( natoms() );
  c1_.resize( natoms() );
  c2_.resize( natoms() );
  c3_.resize( natoms() );
  d1_.resize( natoms() );
  d2_.resize( natoms() );
  not_using_.resize( natoms() );
  charged_.resize( natoms(), true );
  is_chargedH_.resize( natoms() );
	is_freedof_.resize( natoms(), false );

	// Read new charge set specified by "score::facts_charge_dir"
	std::string filename;
	if( option[ score::facts_eq_type ]().compare("apprx") == 0 ){ // Call neutral aliphatic charge set
		filename = option[ in::path::database ](1).name() + "/"
			+ option[ score::facts_eff_charge_dir ]() + "/"
			+ rsd.name() + ".params";
	} else {
		filename = option[ in::path::database ](1).name() + "/"
			+ option[ score::facts_charge_dir ]() + "/"
			+ rsd.name() + ".params";
	}

 	// Option for binding affinity calculation
 	bool const binding_affinity( option[ score::facts_binding_affinity ]() );

	chemical::ResidueTypeOP	rsd_for_charge;

	if( utility::file::file_exists( filename ) ) {
		chemical::AtomTypeSetCOP atom_types = chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
		chemical::ElementSetCOP elements = chemical::ChemicalManager::get_instance()->element_set("default");
		chemical::MMAtomTypeSetCOP mm_atom_types = chemical::ChemicalManager::get_instance()->mm_atom_type_set("fa_standard");
		chemical::orbitals::OrbitalTypeSetCOP orbital_types = chemical::ChemicalManager::get_instance()->orbital_type_set("fa_standard");
		chemical::ResidueTypeSetCOP rsd_type_set = chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		rsd_for_charge = chemical::read_topology_file( filename, atom_types, elements, mm_atom_types, orbital_types, rsd_type_set );

	} else {
		rsd_for_charge = rsd.clone();
	}


	// Assign parameters
  for(Size i = 1; i <= natoms(); ++i){
    core::chemical::AtomType const &type = rsd.atom_type(i);

    // Partial charge
    q_[i] = rsd_for_charge->atom(i).charge();

    if( std::abs(q_[i]) < 1.0e-3 ) charged_[i] = false;

    // Residue polarity
    string atmname( type.atom_type_name() );
    string rsdatomname( rsd.atom_name( i ) );
    bool is_chargedH( atmname.compare( 0, 4, "Hpol" ) == 0 &&
		      (	rsd.aa() == core::chemical::aa_arg ||
			rsd.aa() == core::chemical::aa_lys ||
			rsd.aa() == core::chemical::aa_his ) );

    is_chargedH_[i] = is_chargedH;
    Real vdw_radius = 0.0;

    //Corrections for atomic parameters
    // 1. HIS aromatic carbons to be consistent with CHARMM definition
    if( rsd.aa() == core::chemical::aa_his && atmname == "aroC" ){
      vdw_radius = 1.8;
      alpha_[i] = 0.02;
      COradius2_[i] = 8.37057*8.37057;
      b1_[i] =  0.853399e+1; b2_[i] = -0.109161e+1;
      d1_[i] = -0.227929e+5; d2_[i] = -0.825630e+1;
      a0_[i] = -0.123827e+3; a1_[i] =  0.123827e+3; a2_[i] = 0.185613e-2; a3_[i] =  0.347902e+3;
      c0_[i] =  0.315402e+2; c1_[i] = -0.313502e+2; c2_[i] = 0.567384e-3; c3_[i] = -0.142830e+5;

      // 2. ARG aromatic carbons to be consistent with CHARMM definition
      // copied from backbone carbon
    } else if( rsd.aa() == core::chemical::aa_arg && atmname == "aroC" ){
      vdw_radius = 2.0;
      alpha_[i] = 0.02;
      COradius2_[i] = 8.81363*8.81363;
      b1_[i] = 0.147227e+3; b2_[i] = -0.811304e+0;
      d1_[i] = -0.811741e+4; d2_[i] = -0.217625e+1;
      a0_[i] = -0.858924e+2; a1_[i] = 0.858904e+2; a2_[i] = 0.196363e-2; a3_[i] = 0.900140e+3;
      c0_[i] = 0.168481e+3; c1_[i] = -0.168287e+3; c2_[i] = 0.113765e-2; c3_[i] = -0.672543e+4;

 			// 3. weaken OOC solvation energy - this is useful for binding affinity calculation
 			// if there is uncertainty of interacting explicit water
		} else if( binding_affinity && atmname == "OOC" ){
 			a0_[i] = -0.095500e+3; a1_[i] = 0.095000e+3; a2_[i] = 0.180000e-2; a3_[i] = 0.850000e+3;

      // Otherwise just use default
    } else {
      vdw_radius = type.extra_parameter( FACTS_RADIUS_INDEX );
      if ( vdw_radius <= 1.0e-6 ){
				not_using_[i] = true;
      } else {
				not_using_[i] = false;
      }

      alpha_[i] = type.extra_parameter( FACTS_ALPHA_INDEX );
      COradius2_[i] = type.extra_parameter( FACTS_CUT_INDEX )*type.extra_parameter( FACTS_CUT_INDEX );
      if( COradius2_[i] <= 1.0e-3 && !not_using_[i] ){ // Don't do this for Virtual atoms
	TR << "Unrealistic cutoff radii: set to new cut " << 100.0 << std::endl;
	COradius2_[i] = 10.0*10.0;
      }

      b1_[i] = type.extra_parameter( FACTS_B1_INDEX );
      b2_[i] = type.extra_parameter( FACTS_B2_INDEX );
      d1_[i] = type.extra_parameter( FACTS_D1_INDEX );
      d2_[i] = type.extra_parameter( FACTS_D2_INDEX );
      a0_[i] = type.extra_parameter( FACTS_A0_INDEX );
      a1_[i] = type.extra_parameter( FACTS_A1_INDEX );
      a2_[i] = type.extra_parameter( FACTS_A2_INDEX );
      a3_[i] = type.extra_parameter( FACTS_A3_INDEX );
      c0_[i] = type.extra_parameter( FACTS_C0_INDEX );
      c1_[i] = type.extra_parameter( FACTS_C1_INDEX );
      c2_[i] = type.extra_parameter( FACTS_C2_INDEX );
      c3_[i] = type.extra_parameter( FACTS_C3_INDEX );
    }
    volume_[i] = (4.0/3.0) * Math_PI * vdw_radius * vdw_radius * vdw_radius;

		// Free DOF atoms
		if( rsd.is_lower_terminus() &&
				(	rsdatomname.compare( "1H  " ) == 0 ||
					rsdatomname.compare( "2H  " ) == 0 ||
					rsdatomname.compare( "3H  " ) == 0 )
				){
			is_freedof_[i] = true;
		} else if( rsd.is_upper_terminus() &&
							 (	rsdatomname.compare( " O  " ) == 0 ||
									rsdatomname.compare( " OXT" ) == 0 )
							 ){
			is_freedof_[i] = true;
		}
		else if( rsd.aa() == chemical::aa_ser &&
							 rsdatomname.compare( " HG " ) == 0 ){
			is_freedof_[i] = true;
		} else if( rsd.aa() == chemical::aa_thr &&
							 rsdatomname.compare( " HG1" ) == 0 ){
			is_freedof_[i] = true;
		}
	}

  // Reduce Born radii for charged polarH: sidechain Hpol of ARG/LYS/HIS
  /*
    for(Size i = 1; i <= natoms(); ++i){
    if( option[ score::facts_saltbridge_correction ].user() &&is_chargedH ){
    a0_[i] *= option[ score::facts_saltbridge_correction ]();
    a1_[i] *= option[ score::facts_saltbridge_correction ]();
    TR.Debug << rsd.seqpos() << " " << atmname << " " << a0_[i] << std::endl;
    }
    }
  */

}

void FACTSRsdTypeInfo::initialize_intrascale( chemical::ResidueType const & rsd )
{

	// initialize
  intra_solv_scale_.resize( rsd.natoms() );
  intra_elec_scale_.resize( rsd.natoms() );
  for( Size atm1 = 1; atm1 <= rsd.natoms(); ++atm1 ){
    intra_solv_scale_[atm1].resize( rsd.natoms(), 0.0 );
		intra_elec_scale_[atm1].resize( rsd.natoms(), 0.0 );
	}

	// If is not amino acid, turn on full scale
	if( !rsd.is_protein() ){
		for( Size atm1 = 1; atm1 <= rsd.natoms(); ++atm1 ){
			for( Size atm2 = 2; atm2 <= rsd.natoms(); ++atm2 )
				intra_solv_scale_[atm1][atm2] = 1.0;
		}
		return;
	}

	// Start here if amino acid
  bool const plane_to_self( basic::options::option[ basic::options::OptionKeys::score::facts_plane_to_self ]() );

	utility::vector1< Real > const intbb_solv_scale =
		basic::options::option[ basic::options::OptionKeys::score::facts_intbb_solv_scale ]();
	utility::vector1< Real > const intbs_solv_scale =
		basic::options::option[ basic::options::OptionKeys::score::facts_intbs_solv_scale ]();
	utility::vector1< Real > const intsc_solv_scale =
		basic::options::option[ basic::options::OptionKeys::score::facts_intsc_solv_scale ]();
	utility::vector1< Real > const intbb_elec_scale =
		basic::options::option[ basic::options::OptionKeys::score::facts_intbb_elec_scale ]();
	utility::vector1< Real > const intbs_elec_scale =
		basic::options::option[ basic::options::OptionKeys::score::facts_intbs_elec_scale ]();
	utility::vector1< Real > const intsc_elec_scale =
		basic::options::option[ basic::options::OptionKeys::score::facts_intsc_elec_scale ]();

	utility::vector1< std::string > plane_aa;
	if( basic::options::option[ basic::options::OptionKeys::score::facts_plane_aa ].user() ){
		plane_aa = basic::options::option[ basic::options::OptionKeys::score::facts_plane_aa ]();
	}

debug_assert( intbb_elec_scale.size() == 3 );
debug_assert( intbb_solv_scale.size() == 3 );
debug_assert( intsc_elec_scale.size() == 3 );
debug_assert( intsc_solv_scale.size() == 3 );
debug_assert( intbs_elec_scale.size() == 5 );
debug_assert( intbs_solv_scale.size() == 5 );

	bool const intrascale_by_level =
		basic::options::option[ basic::options::OptionKeys::score::facts_intrascale_by_level ]();

  // 1-4 is important for solvation free energy of ASN/GLN/SER/THR
  // backbone rule is important to give reasonable "distinct" bb solvation for ALA/VAL/ILE/LEU
  // >= 1-6 is turned off; this is important to get rid of artifacts from ASP OD-O & GLU OE-O at exposed

  for( Size atm1 = 1; atm1 <= rsd.natoms(); ++atm1 ){
    for( Size atm2 = 1; atm2 <= rsd.natoms(); ++atm2 ){

			Size path_dist( rsd.path_distance(atm1,atm2) );

			bool const is_atm1_bb( rsd.atom_is_backbone( atm1 ) );
			bool const is_atm2_bb( rsd.atom_is_backbone( atm2 ) );

			// Turn full strength for <= 1-3
			if( path_dist <= 2 ){
				intra_solv_scale_[atm1][atm2] = 1.0;
				intra_elec_scale_[atm1][atm2] = 0.0;

			} else if( is_atm1_bb && is_atm2_bb ){ //bb-bb
				if( path_dist == 3 ){
					intra_solv_scale_[atm1][atm2] = intbb_solv_scale[1];
					intra_elec_scale_[atm1][atm2] = intbb_elec_scale[1];
				} else if( path_dist == 4 ){
					intra_solv_scale_[atm1][atm2] = intbb_solv_scale[2];
					intra_elec_scale_[atm1][atm2] = intbb_elec_scale[2];
				} else {
					intra_solv_scale_[atm1][atm2] = intbb_solv_scale[3];
					intra_elec_scale_[atm1][atm2] = intbb_elec_scale[3];
				}

      } else if( is_atm1_bb || is_atm2_bb ){ // bb-sc
				// Terminus correction: turn full strength since fa_dun won't work here
				/*
				if( rsd.is_terminus() ){
					intra_solv_scale_[atm1][atm2] = 1.0;
					intra_elec_scale_[atm1][atm2] = 1.0;

				// aa_specific rules
				} else
				*/
				if(( rsd.aa() == core::chemical::aa_ser || rsd.aa() == core::chemical::aa_thr )
					 && path_dist == 3 ){
					intra_solv_scale_[atm1][atm2] = 1.0;

				} else {
					// Override definition for path_dist if intrascale_by_level
					if( intrascale_by_level && rsd.has( "CA" ) ){
						Size const i_CA = rsd.atom_index("CA");
						if( is_atm1_bb ){
							path_dist = rsd.path_distance( i_CA, atm2 );
						} else if ( is_atm2_bb ){
							path_dist = rsd.path_distance( i_CA, atm1 );
						}
						path_dist ++; // To convert indexing start from Gamma
					}

					if( path_dist <= 5 ){
						intra_solv_scale_[atm1][atm2] = intbs_solv_scale[ path_dist - 2 ];
						intra_elec_scale_[atm1][atm2] = intbs_elec_scale[ path_dist - 2 ];
					} else {
						intra_solv_scale_[atm1][atm2] = intbs_solv_scale[4];
						intra_elec_scale_[atm1][atm2] = intbs_elec_scale[4];
					}
				}

      } else { // sc-sc
				if( path_dist == 3 ){
					intra_solv_scale_[atm1][atm2] = intsc_solv_scale[1];
					intra_elec_scale_[atm1][atm2] = intsc_elec_scale[2];
				} else if( path_dist == 4 ){
					intra_solv_scale_[atm1][atm2] = intsc_solv_scale[2];
					intra_elec_scale_[atm1][atm2] = intsc_elec_scale[2];
				} else {
					intra_solv_scale_[atm1][atm2] = intsc_solv_scale[3];
					intra_elec_scale_[atm1][atm2] = intsc_elec_scale[3];
				}
			}

			// Override special rule for freeDOF atoms - removed
			//if( path_dist >= 4 && ( is_freedof( atm1 ) || is_freedof( atm2 ) )){
			//	intra_elec_scale_[atm1][atm2] = intra_solv_scale_[atm1][atm2];
			//}
		} // atm2
	} // atm1


  // Plane rule overrides
  for( Size atm1 = 1; atm1 <= rsd.natoms(); ++atm1 ){
    for( Size atm2 = 1; atm2 <= rsd.natoms(); ++atm2 ){
      if( plane_to_self && rsd.aa() == core::chemical::aa_arg ) {
				// Add xHHx - NHx / xHHx - xHHx / xHHx - NE / xHHx - HE / NHx - HE / xHHx - CD / NE or NH - CD
				// or in other words, add all pairs in H or E position
				if( ( rsd.atom_name(atm1).compare( 2, 1, "H" ) == 0 || rsd.atom_name(atm1).compare( 2, 1, "E" ) == 0
							|| atm1 == rsd.atom_index("CD") ) &&
						(	rsd.atom_name(atm2).compare( 2, 1, "H" ) == 0 || rsd.atom_name(atm2).compare( 2, 1, "E" ) == 0
							|| atm2 == rsd.atom_index("CD") )
						)
					intra_solv_scale_[atm1][atm2] = 1.0;
      }

      if( plane_to_self && rsd.aa() == core::chemical::aa_his && plane_aa.has_value( "his" ) ) {
				// Add all pairs in G or D or E position
				if(
					 (( rsd.atom_name(atm1).compare( 2, 1, "G" ) == 0 || rsd.atom_name(atm1).compare( 2, 1, "D" ) == 0)&&
						rsd.atom_name(atm1).compare( 2, 1, "E" ) == 0 ) &&
					 (( rsd.atom_name(atm2).compare( 2, 1, "G" ) == 0 || rsd.atom_name(atm2).compare( 2, 1, "D" ) == 0)&&
						rsd.atom_name(atm2).compare( 2, 1, "E" ) == 0 )
						){
					/*
					std::cout << rsd.name() << " " << rsd.atom_name(atm1) << " " << rsd.atom_name(atm2);
					std::cout << " " << intra_solv_scale_[atm1][atm2];
					std::cout << std::endl;
					*/
					intra_solv_scale_[atm1][atm2] = 1.0;
				}
      }

      // TYR OH-HEX
      if( plane_to_self && rsd.aa() == core::chemical::aa_tyr && plane_aa.has_value( "tyr" ) ) {
				// Add OH-aromatic ring pairs
				if( ( atm1 == rsd.atom_index("OH") && rsd.atom_type(atm2).atom_type_name().compare( "aroC" ) == 0 ) ||
						( atm2 == rsd.atom_index("OH") && rsd.atom_type(atm1).atom_type_name().compare( "aroC" ) == 0 ) ||
						( atm1 == rsd.atom_index("OH") && rsd.atom_type(atm2).atom_type_name().compare( "Haro" ) == 0 ) ||
						( atm2 == rsd.atom_index("OH") && rsd.atom_type(atm1).atom_type_name().compare( "Haro" ) == 0 )
						){
					/*
					std::cout << rsd.name() << " " << rsd.atom_name(atm1) << " " << rsd.atom_name(atm2);
					std::cout << " " << intra_solv_scale_[atm1][atm2];
					std::cout << std::endl;
					*/
					intra_solv_scale_[atm1][atm2] = 1.0;
				}
      }

      // TRP dipole - ring
      if( plane_to_self && rsd.aa() == core::chemical::aa_trp && plane_aa.has_value( "trp" ) ) {
				// Add OH-aromatic ring pairs
				if( ( atm1 == rsd.atom_index("NE1") || atm1 == rsd.atom_index("HE1") ) &&
						( rsd.atom_type(atm2).atom_type_name().compare( "aroC" ) == 0 || rsd.atom_type(atm2).atom_type_name().compare( "Haro" ) == 0 ) )
					{
						/*
						std::cout << rsd.name() << " " << rsd.atom_name(atm1) << " " << rsd.atom_name(atm2);
						std::cout << " " << intra_solv_scale_[atm1][atm2];
						std::cout << std::endl;
						*/
						intra_solv_scale_[atm1][atm2] = 1.0;
						intra_solv_scale_[atm2][atm1] = 1.0;
					}
      }

			if( plane_to_self && rsd.aa() == core::chemical::aa_gln && plane_aa.has_value( "gln" ) ) {
				// Add all pairs in E position
				if( rsd.atom_name(atm1).compare( 2, 1, "E" ) == 0 && rsd.atom_name(atm2).compare( 2, 1, "E" ) == 0 )
					intra_solv_scale_[atm1][atm2] = 1.0;

			} else if( plane_to_self && rsd.aa() == core::chemical::aa_asn && plane_aa.has_value( "asn" ) ) {
				// Add all pairs in D position
				if( rsd.atom_name(atm1).compare( 2, 1, "D" ) == 0 && rsd.atom_name(atm2).compare( 2, 1, "D" ) == 0 )
					intra_solv_scale_[atm1][atm2] = 1.0;
			}

    } //atm2
  } //atm1

} // END void initialize_intrascale

/// FACTSRsdTypeInfo

/**************************************************************************************************/
/*                                                                                                */
/*    @brief: The FACTSResidueInfo class provides all the functions, constants and parameters     */
/*            for different atoms, which are required to calculate the solvation free energy of   */
/*                      of a molecule embedded in water using FACTS method                        */
/*                                                                                                */
/**************************************************************************************************/

//This function initializes all the values for FACTS original parameters, atomic volume, Ai, Bi, esolvE, sasa...
void FACTSResidueInfo::initialize(
		conformation::Residue const & rsd,
		FACTSRsdTypeInfoCOP restypeinfo,
		bool const is_rotamer )
{

  natoms_ = rsd.natoms();

  // This variable is used in res_res_burial and evaluate_polar(nonpolar)_energy
  // 06/27/13: Just fully turn on for fast evaluation...
  flag_for_calculation_ = utility::vector1< bool >( natoms(), true );

  // Initialize Arrays
  Vector i( 0.0 );
  nmtr_ = utility::vector1< Vector >( natoms(), i ); // dnmtr of Ai (equation 4 on page 704 of FACTS paper)
  dnmtr_ = utility::vector1< Real >( natoms(), 1.0 ); // dnmtr of Bi (equation 4 on page 704 of FACTS paper)
  Ai_ = utility::vector1< Real >( natoms(), 0.0 );// Ai (equation 3 on page 704 of FACTS paper)
  Bi_ = utility::vector1< Real >( natoms(), 0.0 );// Bi (equation 4 on page 704 of FACTS paper)
  Ci_ = utility::vector1< Real >( natoms(), 0.0 );// Ci (equation 6 on page 704 of FACTS paper)
  Di_ = utility::vector1< Real >( natoms(), 0.0 );// Di (equation 10 on page 704 of FACTS paper)

  esolvE_ = utility::vector1< Real >( natoms(), 0.0 ); // DeltaGi (equation 7 on page 704 of FACTS paper)
  sasa_ = utility::vector1< Real >( natoms(), 0.0 ); // atomic SASA (equation 11 on page 706 of FACTS paper)
  BR_ = utility::vector1< Real >( natoms(), 0.0 ); // BornRadius

  // auxiliary arrays for calculating derivatives
  if( !is_rotamer ){
    dG_dCi_ = utility::vector1< Real >( natoms() );
    dSA_dDi_ = utility::vector1< Real >( natoms() );
    dsolv_dBR_ = utility::vector1< Real >( natoms() );
    dB_dBnmtr_ = utility::vector1< Real >( natoms() );
    dB_dBdnmtr_ = utility::vector1< Real >( natoms() );
    dBR_dG_ = utility::vector1< Real >( natoms() );
    elecF2_ = utility::vector1< Vector >( natoms(), i );
    solvF2d_ = utility::vector1< Vector >( natoms(), i );
    solvF2BR_ = utility::vector1< Vector >( natoms(), i );
    sasaF2_ = utility::vector1< Vector >( natoms(), i );
  }
  restypeinfo_ = restypeinfo;
}

void FACTSResidueInfo::refresh_energy_cache( Size const nres ){
  E_elec_ = utility::vector1< Real >( nres, 0.0 );
  E_solv_ = utility::vector1< Real >( nres, 0.0 );
  E_solv_pair_ = utility::vector1< Real >( nres, 0.0 );
  E_solv_self_ = utility::vector1< Real >( nres, 0.0 );
}

void FACTSResidueInfo::store_xyz( Residue const &rsd ){
  xyz_.resize( 0 );
  for( Size iatm = 1; iatm <= rsd.natoms(); ++iatm ){
    xyz_.push_back( rsd.xyz(iatm) );
  }
}

/**************************************************************************************************/
/*                                                                                                */
/*    @brief: The  class    FACTSRotamerSetInfo                                                   */
/*                                                                                                */
/**************************************************************************************************/

void FACTSRotamerSetInfo::initialize( RotamerSet const & rotamer_set, FACTSRsdTypeMap &rsdtypemap )
{
  Size const nrot( rotamer_set.num_rotamers() );
  residue_info_.resize( nrot );
  for ( Size i=1; i<= nrot; ++i ) {
    core::chemical::ResidueType const &rsdtype = rotamer_set.rotamer(i)->type();
    FACTSRsdTypeMap::const_iterator it = rsdtypemap.find( &rsdtype );
    if ( it == rsdtypemap.end() ) {
      TR << "Adding new FACTS residue type info: " << rsdtype.name() << std::endl;
      FACTSRsdTypeInfoOP rsdtypeinfo( new FACTSRsdTypeInfo );
      rsdtypeinfo->create_info( rsdtype );
      rsdtypemap[ &rsdtype ] = rsdtypeinfo;
      it = rsdtypemap.find( &rsdtype );
    }

    residue_info_[i] = FACTSResidueInfoOP( new FACTSResidueInfo( *rotamer_set.rotamer(i), it->second, true ) );
  }
} // FACTSRotamerSetInfo

} // namespace scoring
} // namespace core
