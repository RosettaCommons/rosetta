// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/RNA_TorsionPotential.cc
/// @brief  RNA_TorsionPotential potential class implementation
/// @author Rhiju Das

// Unit Headers
#include <core/scoring/rna/RNA_TorsionPotential.hh>


// Package Headers
#include <core/scoring/rna/RNA_Util.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AA.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/file/file_sys_util.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/CircularGeneral1D_Func.hh>
#include <core/scoring/constraints/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/FadeFunc.hh>
#include <core/scoring/constraints/SumFunc.hh>
#include <core/scoring/EnergyMap.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/deriv/dihedral_deriv.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>
//Auto Headers
#include <platform/types.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <utility/Bound.hh>
#include <utility/vector0_bool.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/options/BooleanOption.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end

///////////////////////////


static basic::Tracer tr( "core.scoring.rna.RNA_TorsionPotential" );

using namespace ObjexxFCL::fmt;
using numeric::conversions::radians;

namespace core {
namespace scoring {
namespace rna {

// @brief Auto-generated virtual destructor
RNA_TorsionPotential::~RNA_TorsionPotential() {}

	RNA_TorsionPotential::RNA_TorsionPotential():
		rna_tight_torsions_( true ),
		delta_fade_( 10.0 ),
		alpha_fade_( 10.0 ),
		skip_chainbreak_torsions_( basic::options::option[ basic::options::OptionKeys::score::rna_torsion_skip_chainbreak ]() ),
		verbose_( false ),
		use_new_potential_( false ),
		use_2prime_OH_potential_( basic::options::option[ basic::options::OptionKeys::score::use_2prime_OH_potential ]() )
	{
		if( basic::options::option[ basic::options::OptionKeys::score::rna_torsion_potential ].user() ){

			path_to_torsion_files_="scoring/rna/torsion_potentials/" + basic::options::option[ basic::options::OptionKeys::score::rna_torsion_potential ]();

			std::cout << "-----------------------------------------------------------------------------------" << std::endl;
			std::cout << "USER INPUTTED path_to_torsion_files_=" <<  basic::database::full_name(path_to_torsion_files_) << std::endl;


			//Turn on the new torsional potential if the folder name ends wirh "new"
			if ( path_to_torsion_files_ .compare(path_to_torsion_files_.size() - 3, 3, "new") == 0 ) {
				use_new_potential_ = true;
				std::cout << "Path name ends with 'new'... Turn on the new torsional potential" << std::endl;
			}
			std::cout << "-----------------------------------------------------------------------------------" << std::endl;
		}else{
			//path_to_torsion_files_( "scoring/rna/torsion_potentials/rd2008/" ),
			//path_to_torsion_files_( "scoring/rna/torsion_potentials/FINAL_Mar_24_2010_new_delta_zeta_chi/" ),
			//path_to_torsion_files_( "scoring/rna/torsion_potentials/FINAL_April_28_OLD_syn_chi/" ),

			path_to_torsion_files_= "scoring/rna/torsion_potentials/ps_04282011/";

			std::cout << "-----------------------------------------------------------------------------------" << std::endl;
			std::cout << "DEFAULT path_to_torsion_files_=" <<  basic::database::full_name(path_to_torsion_files_) << std::endl;
			std::cout << "-----------------------------------------------------------------------------------" << std::endl;
		}

		//		init_potentials_from_gaussian_parameters();
		init_potentials_from_rna_torsion_database_files();

		init_fade_functions();

	}

	//////////////////////////////////////////////////////////////////////////
	Real
	RNA_TorsionPotential::eval_intrares_energy(core::conformation::Residue const & rsd, pose::Pose const & pose) const
	{
		using namespace core::id;

		if (!rsd.is_RNA() ) return 0.0;

		if(verbose_){
			std::cout << std::endl;
			std::cout << "Intra_res: " << " rsd.seqpos()= " << rsd.seqpos() << std::endl;
			std::cout << std::endl;
		}
		Real score( 0.0 );

		Real const beta= numeric::principal_angle_degrees( rsd.mainchain_torsion( BETA ) );
		Real const gamma= numeric::principal_angle_degrees(  rsd.mainchain_torsion( GAMMA ) );
		Real const delta= numeric::principal_angle_degrees(  rsd.mainchain_torsion( DELTA ) );
		Real const chi= numeric::principal_angle_degrees(  rsd.chi( CHI - NUM_RNA_MAINCHAIN_TORSIONS ) );
		Real const nu2= numeric::principal_angle_degrees(  rsd.chi( NU2 - NUM_RNA_MAINCHAIN_TORSIONS ) );
		Real const nu1= numeric::principal_angle_degrees(  rsd.chi( NU1 - NUM_RNA_MAINCHAIN_TORSIONS ) );
		Real const o2h = numeric::principal_angle_degrees(  rsd.chi( O2H - NUM_RNA_MAINCHAIN_TORSIONS ) );

		if(verbose_) std::cout << "Beta torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::BB, BETA ) ) ) {
			score += beta_potential_->func( beta ); //beta
		}

		if(verbose_) std::cout << "Gamma torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::BB, GAMMA ) ) ) {
			score += gamma_potential_->func( gamma ); //gamma
		}

		if(verbose_) std::cout << "Delta torsion" << std::endl;


		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::BB, DELTA ) ) ) {
			score += ( fade_delta_north_->func( delta ) * delta_north_potential_->func( delta ) +
								 fade_delta_south_->func( delta ) * delta_south_potential_->func( delta ) ); //delta
		}

		if(verbose_)  std::cout << "Chi torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::CHI, CHI - NUM_RNA_MAINCHAIN_TORSIONS ) ) ) {
			if (use_new_potential_) {
				if ( rsd.aa() == core::chemical::na_rgu || rsd.aa() == core::chemical::na_rad  ) {
					score += ( fade_delta_north_->func( delta ) * chi_purine_north_potential_->func( chi ) +
									 fade_delta_south_->func( delta ) * chi_purine_south_potential_->func( chi ) ); //chi
				}else{
					score += ( fade_delta_north_->func( delta ) * chi_pyrimidine_north_potential_->func( chi ) +
									 fade_delta_south_->func( delta ) * chi_pyrimidine_south_potential_->func( chi ) ); //chi
				}
			} else {
				if ( rsd.aa() == core::chemical::na_rgu ) {
					score += ( fade_delta_north_->func( delta ) * chi_north_potential_guanosine_->func( chi ) +
									 fade_delta_south_->func( delta ) * chi_south_potential_guanosine_->func( chi ) ); //chi
				}else{
					score += ( fade_delta_north_->func( delta ) * chi_north_potential_others_->func( chi ) +
									 fade_delta_south_->func( delta ) * chi_south_potential_others_->func( chi ) ); //chi
				}
			}
		}

		if(verbose_)  std::cout << "nu2 torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::CHI, NU2 - NUM_RNA_MAINCHAIN_TORSIONS ) ) ) {
			score += ( fade_delta_north_->func( delta ) * nu2_north_potential_->func( nu2 ) +
								 fade_delta_south_->func( delta ) * nu2_south_potential_->func( nu2 ) ); //nu2
		}

		if(verbose_)  std::cout << "nu1 torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::CHI, NU1 - NUM_RNA_MAINCHAIN_TORSIONS ) ) ) {
			score += ( fade_delta_north_->func( delta ) * nu1_north_potential_->func( nu1 ) +
							 fade_delta_south_->func( delta ) * nu1_south_potential_->func( nu1 ) ); //nu1
		}

		if(use_2prime_OH_potential_ && use_new_potential_ && Should_score_torsion(pose, TorsionID( rsd.seqpos(), id::CHI, O2H - NUM_RNA_MAINCHAIN_TORSIONS ) ) ) {
			score += ( fade_delta_north_->func( delta ) * o2h_north_potential_->func( o2h ) +
							 fade_delta_south_->func( delta ) * o2h_south_potential_->func( o2h ) ); //o2h
		}

		return score;

	}

	//////////////////////////////////////////////////////////////////////////
	Real
	RNA_TorsionPotential::residue_pair_energy(core::conformation::Residue const & rsd1, core::conformation::Residue const & rsd2, pose::Pose const & pose) const
	{

		using namespace core::id;

		// CAN WE ASSUME rsd1 < rsd2 ???

//THESE ASSERTION FAIL!! TRY RERUN ~/minirosetta/test/Mar_27_torsion_potential_fix/1zih/new_RNA_TorssionPotential/trail_2/output.txt
//		if((rsd1.seqpos() < rsd2.seqpos())==false){
//			std::cout << "rsd1.seqpos()= " << rsd1.seqpos() << " rsd2.seqpos()= " << rsd2.seqpos() << std::endl;
//			utility_exit_with_message( "(rsd1.seqpos() < rsd2.seqpos())==false" );
//		}

//		assert( rsd1.seqpos() < rsd2.seqpos() );

		if ( rsd1.seqpos() != (rsd2.seqpos() - 1) ) return 0.0;
		if (!rsd1.is_RNA() ) return 0.0;
		if (!rsd2.is_RNA() ) return 0.0;

		if(verbose_)  {
			std::cout << std::endl;
			std::cout << "Between_res= " << " rsd1.seqpos()= " << rsd1.seqpos() << " rsd2.seqpos()= " << rsd2.seqpos() << std::endl;
			std::cout << std::endl;
		}

		Real score( 0.0 );

		Real const delta= numeric::principal_angle_degrees(  rsd1.mainchain_torsion( DELTA ) );
		Real const epsilon= numeric::principal_angle_degrees(  rsd1.mainchain_torsion( EPSILON ) );
		Real const zeta= numeric::principal_angle_degrees(  rsd1.mainchain_torsion( ZETA ) );
		Real const next_alpha= numeric::principal_angle_degrees(  rsd2.mainchain_torsion( ALPHA ) );

		if(verbose_)  std::cout << "epsilon torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd1.seqpos(), id::BB, EPSILON ) ) ) {
			score += ( fade_delta_north_->func( delta ) * epsilon_north_potential_->func( epsilon ) +
								 fade_delta_south_->func( delta ) * epsilon_south_potential_->func( epsilon ) );
		}

		if(verbose_)  std::cout << "zeta torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd1.seqpos(), id::BB, ZETA ) ) ) {
			score += ( fade_alpha_sc_minus_->func( next_alpha ) * zeta_alpha_sc_minus_potential_->func( zeta ) +
							 fade_alpha_sc_plus_->func( next_alpha )  * zeta_alpha_sc_plus_potential_->func( zeta ) +
							 fade_alpha_ap_->func( next_alpha )       * zeta_alpha_ap_potential_->func( zeta ) );
		}

		if(verbose_)  std::cout << "next alpha torsion" << std::endl;

		if(Should_score_torsion(pose, TorsionID( rsd2.seqpos(), id::BB, ALPHA ) ) ) {
			score += alpha_potential_->func( next_alpha );
		}

		return score;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_TorsionPotential::get_f1_f2( core::id::TorsionID const & torsion_id, pose::Pose const & pose, core::id::AtomID const & id, Vector & f1, Vector & f2 ) const
	{

			conformation::Conformation const & conformation( pose.conformation() );

			if ( !pose.residue_type( torsion_id.rsd() ).is_RNA() ) return false;

			// Check that torsion is intraresidue.
			id::AtomID id1,id2,id3,id4;
			if  ( conformation.get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 ) ) return false;

			//Kinda hacky, but most succinct this way
			if(verbose_) {
				tr.Info << "In get_f1_f2_function: ";
				tr.Info << " atom_id: " << id << std::endl;
			}
		  if(Should_score_torsion(pose, torsion_id)==false) return false;

			Real theta( 0.0 );
			// copied from DihedralConstraint. Better work, damnit.
			if ( id == id1 ) {
				numeric::deriv::dihedral_p1_cosine_deriv( conformation.xyz( id1 ),conformation.xyz( id2 ),
																									conformation.xyz( id3 ),conformation.xyz( id4 ), theta, f1, f2 );
			} else if ( id == id2 ) {
				numeric::deriv::dihedral_p2_cosine_deriv( conformation.xyz( id1 ),conformation.xyz( id2 ),
																									conformation.xyz( id3 ),conformation.xyz( id4 ), theta, f1, f2 );
			} else if ( id == id3 ) {
				numeric::deriv::dihedral_p2_cosine_deriv( conformation.xyz( id4 ),conformation.xyz( id3 ),
																									conformation.xyz( id2 ),conformation.xyz( id1 ), theta, f1, f2 );
			} else if ( id == id4 ) {
				numeric::deriv::dihedral_p1_cosine_deriv( conformation.xyz( id4 ),conformation.xyz( id3 ),
																									conformation.xyz( id2 ),conformation.xyz( id1 ), theta, f1, f2 );
			} else {
				return false;
			}
			return true;
	}

	//////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::eval_atom_derivative(
																						 id::AtomID const & id,
																						 pose::Pose const & pose,
																						 EnergyMap const & weights,
																						 Vector & F1,
																						 Vector & F2
																						 ) const{

		Real const radians2degrees = 1.0 / radians( 1.0 );

		Size const current_seqpos( id.rsd() );
		Size const nres( pose.total_residue() );

		// Look for torsions in current, previous, and next residues that might match up with this atom_id.
		// This isn't totally efficient, but the extra checks should not be rate limiting...
		for ( int offset = -1; offset <= +1; offset++ ) {

			Size const seqpos = current_seqpos + offset;

			if (seqpos < 1 ) continue;
			if (seqpos > nres ) continue;
			if ( !pose.residue_type( seqpos ).is_RNA() ) continue;

			conformation::Residue const & rsd( pose.residue( seqpos ) );

			Real const alpha= numeric::principal_angle_degrees(   rsd.mainchain_torsion( ALPHA ) );
			Real const beta = numeric::principal_angle_degrees(   rsd.mainchain_torsion( BETA ) );
			Real const gamma= numeric::principal_angle_degrees(   rsd.mainchain_torsion( GAMMA ) );
			Real const delta= numeric::principal_angle_degrees(   rsd.mainchain_torsion( DELTA ) );
			Real const epsilon= numeric::principal_angle_degrees( rsd.mainchain_torsion( EPSILON ) );
			Real const zeta= numeric::principal_angle_degrees(    rsd.mainchain_torsion( ZETA ) );
			Real const chi  = numeric::principal_angle_degrees(   rsd.chi( CHI - NUM_RNA_MAINCHAIN_TORSIONS ) );
			Real const nu2  = numeric::principal_angle_degrees(   rsd.chi( NU2 - NUM_RNA_MAINCHAIN_TORSIONS ) );
			Real const nu1  = numeric::principal_angle_degrees(   rsd.chi( NU1 - NUM_RNA_MAINCHAIN_TORSIONS ) );
			Real const o2h  = numeric::principal_angle_degrees(   rsd.chi( O2H - NUM_RNA_MAINCHAIN_TORSIONS ) );

			Vector f1( 0.0 ), f2( 0.0 );

			///////////////////////////////ALPHA//////////////////////////////
			if ( seqpos > 1 && pose.residue( seqpos - 1 ).is_RNA() && get_f1_f2( id::TorsionID( seqpos, id::BB, ALPHA ),	pose, id, f1, f2 ) ){

				Real dE_dtorsion = alpha_potential_->dfunc( alpha );

				Real const previous_zeta = numeric::principal_angle_degrees( pose.residue( seqpos - 1 ).mainchain_torsion( ZETA ) );  ///NEED TO CHANGE THIS AS WELL
				dE_dtorsion += ( fade_alpha_sc_minus_->dfunc( alpha ) * zeta_alpha_sc_minus_potential_->func( previous_zeta ) +
												 fade_alpha_sc_plus_->dfunc( alpha )  * zeta_alpha_sc_plus_potential_->func( previous_zeta ) +
												 fade_alpha_ap_->dfunc( alpha )       * zeta_alpha_ap_potential_->func( previous_zeta ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////Beta/////////////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::BB, BETA ), pose, id, f1, f2 ) ){
				Real const dE_dtorsion = beta_potential_->dfunc( beta );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////Gamma/////////////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::BB, GAMMA ),	pose, id, f1, f2 ) ){
				Real const dE_dtorsion = gamma_potential_->dfunc( gamma );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////Delta////////////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::BB, DELTA ),	pose, id, f1, f2 ) ){
				Real dE_dtorsion = ( fade_delta_north_->dfunc( delta ) * delta_north_potential_->func( delta ) +
														 fade_delta_north_->func( delta )  * delta_north_potential_->dfunc( delta ) +
														 fade_delta_south_->dfunc( delta ) * delta_south_potential_->func( delta ) +
														 fade_delta_south_->func( delta )  * delta_south_potential_->dfunc( delta ) );

				if (use_new_potential_) {
					if (  rsd.aa() == core::chemical::na_rgu || rsd.aa() == core::chemical::na_rad ) {
						dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * chi_purine_north_potential_->func( chi ) +
														 fade_delta_south_->dfunc( delta ) * chi_purine_south_potential_->func( chi ) );
					}else{
						dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * chi_pyrimidine_north_potential_->func( chi ) +
														 fade_delta_south_->dfunc( delta ) * chi_pyrimidine_south_potential_->func( chi ) );
					}
				} else {
					if ( rsd.aa() == core::chemical::na_rgu ) {
						dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * chi_north_potential_guanosine_->func( chi ) +
														 fade_delta_south_->dfunc( delta ) * chi_south_potential_guanosine_->func( chi ) );
					}else{
						dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * chi_north_potential_others_->func( chi ) +
														 fade_delta_south_->dfunc( delta ) * chi_south_potential_others_->func( chi ) );
					}
				}

				dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * nu2_north_potential_->func( nu2 ) +
												 fade_delta_south_->dfunc( delta ) * nu2_south_potential_->func( nu2 ) );

				dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * nu1_north_potential_->func( nu1 ) +
												 fade_delta_south_->dfunc( delta ) * nu1_south_potential_->func( nu1 ) );

				if ( seqpos < nres && pose.residue( seqpos+1 ).is_RNA() ){
					dE_dtorsion += ( fade_delta_north_->dfunc( delta ) * epsilon_north_potential_->func( epsilon ) +
																	 fade_delta_south_->dfunc( delta ) * epsilon_south_potential_->func( epsilon ) );
				}

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			////////////////////////////////////EPSILON////////////////////////////////////////////////
			if ( seqpos < nres && pose.residue( seqpos+1 ).is_RNA() && get_f1_f2( id::TorsionID( seqpos, id::BB, EPSILON ),	pose, id, f1, f2 ) ){
				Real const dE_dtorsion = ( fade_delta_north_->func( delta ) * epsilon_north_potential_->dfunc( epsilon ) +
																	 fade_delta_south_->func( delta ) * epsilon_south_potential_->dfunc( epsilon ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			////////////////////////////////////ZETA////////////////////////////////////////////////
			if ( seqpos < nres && pose.residue( seqpos+1 ).is_RNA() && get_f1_f2( id::TorsionID( seqpos, id::BB, ZETA ),	pose, id, f1, f2 ) ){

				Real const next_alpha = numeric::principal_angle_degrees(   pose.residue( seqpos+1 ).mainchain_torsion( ALPHA ) );  ///NEED TO CHANGE THIS AS WELL
				Real const dE_dtorsion = ( fade_alpha_sc_minus_->func( next_alpha ) * zeta_alpha_sc_minus_potential_->dfunc( zeta ) +
																	 fade_alpha_sc_plus_->func( next_alpha )  * zeta_alpha_sc_plus_potential_->dfunc( zeta ) +
																	 fade_alpha_ap_->func( next_alpha )       * zeta_alpha_ap_potential_->dfunc( zeta ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////////CHI///////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::CHI, CHI - NUM_RNA_MAINCHAIN_TORSIONS ),	pose, id, f1, f2 ) ){

				Real dE_dtorsion;

				if (use_new_potential_) {
					if ( rsd.aa() == core::chemical::na_rgu || rsd.aa() == core::chemical::na_rad ) {
						dE_dtorsion = ( fade_delta_north_->func( delta ) * chi_purine_north_potential_->dfunc( chi ) +
														fade_delta_south_->func( delta ) * chi_purine_south_potential_->dfunc( chi ) );
					} else {
						dE_dtorsion = ( fade_delta_north_->func( delta ) * chi_pyrimidine_north_potential_->dfunc( chi ) +
														fade_delta_south_->func( delta ) * chi_pyrimidine_south_potential_->dfunc( chi ) );
					}
				} else {
					if ( rsd.aa() == core::chemical::na_rgu ) {
						dE_dtorsion = ( fade_delta_north_->func( delta ) * chi_north_potential_guanosine_->dfunc( chi ) +
														fade_delta_south_->func( delta ) * chi_south_potential_guanosine_->dfunc( chi ) );

					} else {
						dE_dtorsion = ( fade_delta_north_->func( delta ) * chi_north_potential_others_->dfunc( chi ) +
														fade_delta_south_->func( delta ) * chi_south_potential_others_->dfunc( chi ) );
					}
				}

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////////NU2//////////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::CHI, NU2 - NUM_RNA_MAINCHAIN_TORSIONS ),	pose, id, f1, f2 ) ){
				Real const dE_dtorsion = ( fade_delta_north_->func( delta ) * nu2_north_potential_->dfunc( nu2 ) +
																	 fade_delta_south_->func( delta ) * nu2_south_potential_->dfunc( nu2 ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////////NU1//////////////////////////////////////////////////////
			if ( get_f1_f2( id::TorsionID( seqpos, id::CHI, NU1 - NUM_RNA_MAINCHAIN_TORSIONS ),	pose, id, f1, f2 ) ){
				Real const dE_dtorsion = ( fade_delta_north_->func( delta ) * nu1_north_potential_->dfunc( nu1 ) +
																	 fade_delta_south_->func( delta ) * nu1_south_potential_->dfunc( nu1 ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}

			/////////////////////////////////O2H/////////////////////////////////////////////////////////
			if ( use_2prime_OH_potential_ && use_new_potential_ &&
			     get_f1_f2( id::TorsionID( seqpos, id::CHI, O2H - NUM_RNA_MAINCHAIN_TORSIONS ),	pose, id, f1, f2 ) ){
				Real const dE_dtorsion = ( fade_delta_north_->func( delta ) * o2h_north_potential_->dfunc( o2h ) +
																	 fade_delta_south_->func( delta ) * o2h_south_potential_->dfunc( o2h ) );

				F1 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f1;
				F2 += radians2degrees * dE_dtorsion * weights[ rna_torsion ] * f2;
			}
		}
	}


	////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_TorsionPotential::check_intra_residue( id::TorsionID const & torsion_id, pose::Pose const & pose, Size const seqpos ) const{

		// Check that torsion is intraresidue.
		id::AtomID id1,id2,id3,id4;
		bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
		if (fail) return false;

		if ( id1.rsd() != seqpos ) return false;
		if ( id2.rsd() != seqpos ) return false;
		if ( id3.rsd() != seqpos ) return false;
		if ( id4.rsd() != seqpos ) return false;

		return true;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::init_potentials_from_rna_torsion_database_files() {

		// Initialize potential functions.
		initialize_potential_from_file( alpha_potential_, "alpha_potential.txt" );
		initialize_potential_from_file( beta_potential_, "beta_potential.txt" );
		initialize_potential_from_file( gamma_potential_, "gamma_potential.txt" );
		initialize_potential_from_file( delta_north_potential_, "delta_north_potential.txt" );
		initialize_potential_from_file( delta_south_potential_, "delta_south_potential.txt" );
		initialize_potential_from_file( epsilon_north_potential_, "epsilon_north_potential.txt" );
		initialize_potential_from_file( epsilon_south_potential_, "epsilon_south_potential.txt" );
		initialize_potential_from_file( zeta_alpha_sc_minus_potential_, "zeta_alpha_sc_minus_potential.txt" );
		initialize_potential_from_file( zeta_alpha_sc_plus_potential_, "zeta_alpha_sc_plus_potential.txt" );
		initialize_potential_from_file( zeta_alpha_ap_potential_, "zeta_alpha_ap_potential.txt" );
		initialize_potential_from_file( nu2_north_potential_, "nu2_north_potential.txt" );
		initialize_potential_from_file( nu2_south_potential_, "nu2_south_potential.txt" );
		initialize_potential_from_file( nu1_north_potential_, "nu1_north_potential.txt" );
		initialize_potential_from_file( nu1_south_potential_, "nu1_south_potential.txt" );

		if (use_new_potential_) {
			initialize_potential_from_file( chi_purine_north_potential_, "chi_purine_north_potential.txt" );
			initialize_potential_from_file( chi_purine_south_potential_, "chi_purine_south_potential.txt" );
			initialize_potential_from_file( chi_pyrimidine_north_potential_, "chi_pyrimidine_north_potential.txt" );
			initialize_potential_from_file( chi_pyrimidine_south_potential_, "chi_pyrimidine_south_potential.txt" );
			initialize_potential_from_file( o2h_north_potential_, "o2h_north_potential.txt" );
			initialize_potential_from_file( o2h_south_potential_, "o2h_south_potential.txt" );
		} else {
			initialize_potential_from_file( chi_north_potential_guanosine_, "chi_north_potential_guanosine.txt" );
			initialize_potential_from_file( chi_south_potential_guanosine_, "chi_south_potential_guanosine.txt" );
			initialize_potential_from_file( chi_north_potential_others_, "chi_north_potential_others.txt" );
			initialize_potential_from_file( chi_south_potential_others_, "chi_south_potential_others.txt" );
		}

	}

	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::initialize_potential_from_file( core::scoring::constraints::FuncOP & func,
																	std::string const & filename ) {

		std::string const full_filename = basic::database::full_name( path_to_torsion_files_ + "/"+filename  );

		//std::cout << "full_torsional_potential_filename= " << full_filename << std::endl;
		if(utility::file::file_exists( full_filename )==false){
			utility_exit_with_message( "full_torsional_potential_filename " + full_filename + " doesn't exist!" );
		}

		func = new core::scoring::constraints::CircularGeneral1D_Func( full_filename );

	}


	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::init_fade_functions()
	{

		using namespace scoring::constraints;

		RNA_FittedTorsionInfo rna_torsion_fitted_info;
		Real const DELTA_CUTOFF_ = rna_torsion_fitted_info.delta_cutoff();

		// FadeFunc initialized with min, max, fade-width, and function value.
		fade_delta_north_ = new FadeFunc(
																		 -180.0 -delta_fade_,
																		 DELTA_CUTOFF_ + 0.5*delta_fade_ ,
																		 delta_fade_ ,
																		 1.0  );
		fade_delta_south_ = new FadeFunc(
																		 DELTA_CUTOFF_ - 0.5*delta_fade_,
																		 180.0 +delta_fade_ ,
																		 delta_fade_ ,
																		 1.0  );


		// FadeFunc initialized with min, max, fade-width, and function value.
		fade_alpha_sc_minus_ = new FadeFunc(
																				-120.0 - 0.5 * alpha_fade_ ,
																				0.0 + 0.5 * alpha_fade_ ,
																				alpha_fade_ ,
																				1.0  );
		fade_alpha_sc_plus_ = new FadeFunc(
																			 0.0 - 0.5 * alpha_fade_ ,
																			 100.0 + 0.5 * alpha_fade_ ,
																			 alpha_fade_ ,
																			 1.0  );

		fade_alpha_ap_ = new SumFunc();
		fade_alpha_ap_->add_func(
														 new FadeFunc(
																					-180.0 - alpha_fade_ ,
																					-120.0 + 0.5 * alpha_fade_ ,
																					alpha_fade_ ,
																					1.0  ) );

		fade_alpha_ap_->add_func(
														 new FadeFunc(
																					 100.0 - 0.5 * alpha_fade_ ,
																					 180.0 + alpha_fade_ ,
																					 alpha_fade_ ,
																					 1.0  ) );
	}

	//HACKY	/////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::Output_boolean(std::string const & tag, bool boolean) const {

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;
		std::cout << tag;

		if(boolean==true){
			std::cout << A(4,"T");
		} else {
			std::cout << A(4,"F");
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_TorsionPotential::Is_cutpoint_closed_atom(core::conformation::Residue const & rsd, id::AtomID const & id) const{
		std::string const & atom_name=rsd.type().atom_name(id.atomno());

		if(atom_name=="OVU1" || atom_name=="OVL1" || atom_name=="OVL2"){
			return true;
		}else{
			return false;
		}

	}
	/////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_TorsionPotential::Is_cutpoint_closed_torsion(pose::Pose const & pose, id::TorsionID const & torsion_id) const
	{
		using namespace ObjexxFCL;

		Size torsion_seq_num=torsion_id.rsd();
		Size lower_seq_num=0;
		Size upper_seq_num=0;


		if ( torsion_id.type() != id::BB) return false;

		if ( torsion_id.torsion() == ALPHA){ //COULD BE A UPPER RESIDUE OF A CHAIN_BREAK_CLOSE

			lower_seq_num=torsion_seq_num-1;
			upper_seq_num=torsion_seq_num;

		}else if(torsion_id.torsion() == EPSILON || torsion_id.torsion() == ZETA){
			lower_seq_num=torsion_seq_num;
			upper_seq_num=torsion_seq_num+1;
		}else{
			if( torsion_id.torsion()!=DELTA && torsion_id.torsion() != BETA && torsion_id.torsion() != GAMMA){
				utility_exit_with_message("The torsion should be DELTA(lower), BETA(upper) or GAMMA(upper) !!" );
			}
			return false;
		}

		if(upper_seq_num==1) return false;

		if(lower_seq_num==pose.total_residue()) return false;

		if(pose.residue( lower_seq_num ).has_variant_type( chemical::CUTPOINT_LOWER ) ){
			if( pose.residue( upper_seq_num ).has_variant_type( chemical::CUTPOINT_UPPER )==false ){
				utility_exit_with_message("seq_num " + string_of(lower_seq_num) + " is a CUTPOINT_LOWER but seq_num " + string_of(upper_seq_num) + " is not a cutpoint CUTPOINT_UPPER??" );
			}
			return true;
		}

		return false;
	}


	///////////////////////////////////////////////Dec 26, 2010//////////////////////////////////////////////////////////
	void
	RNA_TorsionPotential::print_torsion_info(pose::Pose const & pose, id::TorsionID const & torsion_id) const
	{

		id::AtomID id1,id2,id3,id4;
		pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );

		std::cout << "torsion_id: " << torsion_id << std::endl;

		bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
		if (fail){
			std::cout << "fail to get torsion!, perhap this torsion is located at a chain_break " << std::endl;
			return;
		}

		conformation::Residue const & rsd_1=pose.residue(id1.rsd());
		conformation::Residue const & rsd_2=pose.residue(id2.rsd());
		conformation::Residue const & rsd_3=pose.residue(id3.rsd());
		conformation::Residue const & rsd_4=pose.residue(id4.rsd());

		tr.Info << " Torsion containing one or more virtual atom(s)" << std::endl;
		tr.Info << "  torsion_id: " << torsion_id;
		tr.Info << "  atom_id: " << id1 << " " << id2 << " " << id3 << " " << id4 << std::endl;
		tr.Info << "  name: " << rsd_1.type().atom_name(id1.atomno()) << " " << rsd_2.type().atom_name(id2.atomno()) << " " << rsd_3.type().atom_name(id3.atomno()) << " " << rsd_4.type().atom_name(id4.atomno()) << std::endl;
		tr.Info << "  type: " << rsd_1.atom_type(id1.atomno()).name() << " " << rsd_2.atom_type(id2.atomno()).name() << " " << rsd_3.atom_type(id3.atomno()).name() << " " << rsd_4.atom_type(id4.atomno()).name() << std::endl;
		tr.Info << "		atom_type_index: " << rsd_1.atom_type_index( id1.atomno()) << " " << rsd_2.atom_type_index( id2.atomno()) << " " << rsd_3.atom_type_index( id3.atomno())  << " " << rsd_4.atom_type_index( id4.atomno()) << std::endl;
		tr.Info << "		atomic_charge: " << rsd_1.atomic_charge( id1.atomno())	<< " " << rsd_2.atomic_charge( id2.atomno())	<< " " << rsd_3.atomic_charge( id3.atomno())	<< " " << rsd_4.atomic_charge( id4.atomno()) << std::endl;


	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool
	RNA_TorsionPotential::Should_score_torsion(pose::Pose const & pose, id::TorsionID const & torsion_id) const
	{

		id::AtomID id1,id2,id3,id4;
		pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );

		if(verbose_) std::cout << "torsion_id: " << torsion_id << std::endl;

		bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
		if (fail){
			if(verbose_) std::cout << "fail to get torsion!, perhap this torsion is located at a chain_break " << std::endl;
			return false;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		conformation::Residue const & rsd_1=pose.residue(id1.rsd());
		conformation::Residue const & rsd_2=pose.residue(id2.rsd());
		conformation::Residue const & rsd_3=pose.residue(id3.rsd());
		conformation::Residue const & rsd_4=pose.residue(id4.rsd());

		if ( !rsd_1.is_RNA() || !rsd_2.is_RNA() || !rsd_3.is_RNA() || !rsd_4.is_RNA() ) return false;  //Rhiju's added this between r47018 and r46616
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		bool Is_virtual_torsion=( rsd_1.is_virtual(id1.atomno()) || rsd_2.is_virtual(id2.atomno()) || rsd_3.is_virtual(id3.atomno()) || rsd_4.is_virtual(id4.atomno()) );

		if(Is_virtual_torsion && verbose_) print_torsion_info(pose, torsion_id);


		/////////////////////Check for cutpoint_closed (Since these torsions will contain virtual atom(s), but want to score these torsions!///////////////////////////////
		//Method 1:
		bool const METHOD_ONE_Is_cutpoint_closed_torsion = Is_cutpoint_closed_torsion(pose, torsion_id);

    //Method 2:
		bool const METHOD_TWO_Is_cutpoint_closed_torsion = ( Is_cutpoint_closed_atom(rsd_1, id1) || Is_cutpoint_closed_atom(rsd_2, id2) || Is_cutpoint_closed_atom(rsd_3, id3) || Is_cutpoint_closed_atom(rsd_4, id4) );


		if( METHOD_ONE_Is_cutpoint_closed_torsion != METHOD_TWO_Is_cutpoint_closed_torsion){
			Output_boolean(" METHOD_ONE_Is_cutpoint_closed_torsion= ", METHOD_ONE_Is_cutpoint_closed_torsion);
			Output_boolean(" METHOD_TWO_Is_cutpoint_closed_torsion= ", METHOD_TWO_Is_cutpoint_closed_torsion);
			Output_boolean(" Is_virtual_torsion= ", Is_virtual_torsion); std::cout << std::endl;

			print_torsion_info(pose, torsion_id);

			utility_exit_with_message( "METHOD_ONE_Is_cutpoint_closed_torsion != METHOD_TWO_Is_cutpoint_closed_torsion !!" );
		}

		if( METHOD_ONE_Is_cutpoint_closed_torsion==true && Is_virtual_torsion==false ){

			print_torsion_info(pose, torsion_id);

			utility_exit_with_message( "METHOD_ONE_Is_cutpoint_closed_torsion==true && Is_virtual_torsion==false !!" );
		}


		////////////////////////////////////////////////Jan 19, 2012: New code//////////////////////////////////////////
		if(rsd_1.seqpos()>rsd_2.seqpos()) utility_exit_with_message("rsd_1.seqpos()>rsd_2.seqpos()");
		if(rsd_2.seqpos()>rsd_3.seqpos()) utility_exit_with_message("rsd_1.seqpos()>rsd_2.seqpos()");
		if(rsd_3.seqpos()>rsd_4.seqpos()) utility_exit_with_message("rsd_1.seqpos()>rsd_2.seqpos()");

		if( (rsd_1.seqpos()!=rsd_4.seqpos()) && (rsd_1.seqpos()!=(rsd_4.seqpos()-1)) ){
			utility_exit_with_message("(rsd_1.seqpos()!=rsd_4.seqpos()) && (rsd_1.seqpos()!=(rsd_4.seqpos()-1))");
		}

		bool const inter_residue_torsion= (rsd_1.seqpos()!=rsd_4.seqpos());

		bool Is_chain_break_torsion=false;

		if(inter_residue_torsion){

			//Note that chain_break_torsion does not neccessarily have to be located at a cutpoint_open. For example, in RNA might contain multiple strands, but the user not have specified them as cutpoint_open
			//This happen frequently, for example when modeling single-stranded RNA loop (PNAS December 20, 2011 vol. 108 no. 51 20573-20578).
			//Actually if chain_break is cutpoint_open, pose.conformation().get_torsion_angle_atom_ids() should fail, which leads to the EARLY RETURN FALSE statement at the beginning of this function.

			bool const violate_max_O3_prime_to_P_bond_dist= is_rna_chainbreak(pose, rsd_1.seqpos());

			//Note that cutpoint_closed_torsions are NOT considered as chain_break_torsion since we want to score them EVEN when skip_chainbreak_torsions_=true!
			//Necessary since for cutpoint_closed_torsions, the max O3_prime_to_P_bond_dist might be violated during stages of the Fragment Assembly and Stepwise Assembly where the chain is not yet closed.

			Is_chain_break_torsion= (violate_max_O3_prime_to_P_bond_dist==true && METHOD_ONE_Is_cutpoint_closed_torsion==false);

		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bool should_score_this_torsion= true; //Warning before Jan 20, 2012, this used to be "Size should_score_this_torsion= true;"

		if( Is_virtual_torsion==true && METHOD_ONE_Is_cutpoint_closed_torsion==false ) should_score_this_torsion=false;

		if( skip_chainbreak_torsions_ && Is_chain_break_torsion ) should_score_this_torsion=false;

		if(verbose_){
			Output_boolean(" should_score_torsion= ", should_score_this_torsion);
			Output_boolean(" | Is_cutpoint_closed_torsion= ", METHOD_ONE_Is_cutpoint_closed_torsion);
			Output_boolean(" | Is_virtual_torsion= ", Is_virtual_torsion);
			Output_boolean(" | skip_chainbreak_torsions_= ", skip_chainbreak_torsions_);
			Output_boolean(" | Is_chain_break_torsion   = ", Is_chain_break_torsion) ;
			std::cout << std::endl;
		}

		return should_score_this_torsion;


	}



}
}
}
