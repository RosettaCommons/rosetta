// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBonds.cxxtest.hh
/// @brief  Test evaluation of hbond potentials for each HBEvalType across parameter ranges.
/// @author Matthew O'Meara

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Package headers
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondTypeManager.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>

// Project headers
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility headers
#include <utility/vector1.hh>

//Auto Headers
#include <core/scoring/hbonds/HBondOptions.hh>
#include <numeric/constants.hh>

using namespace core;
  using namespace conformation;
  using namespace scoring;
    using namespace hbonds;

using pose::Pose;
using conformation::Residue;

static basic::Tracer TR("core.scoring.hbonds.HBonds.cxxtest");

class HBondsTest : public CxxTest::TestSuite {

public:
  void setUp() {
		core_init();

	}

	void tearDown(){}

	// Classic hbond types

//	void dont_test_hbond_compute_energy()
//  {
//		test::UTracer UT("core/scoring/hbonds/hbond_compute_energy.u");
//
//		Real energy, dE_dr, dE_dxD, dE_dxH;
//
//		UT << "# Computed energies and derivatives given HBEvalType and geometric parameters\n";
//		UT << "# HBEvalType\tAHdis\txD\txH\tenergy\n";
//		utility::vector1<HBEvalType> hbe_types;
//
//		hbe_types.push_back(hbe_BB);
//		hbe_types.push_back(hbe_BBTURN);
//		hbe_types.push_back(hbe_BBHELIX);
//		hbe_types.push_back(hbe_BBOTHER);
//		hbe_types.push_back(hbe_SP2B);
//		hbe_types.push_back(hbe_SP3B);
//		hbe_types.push_back(hbe_RINGB);
//		hbe_types.push_back(hbe_BSC);
//		hbe_types.push_back(hbe_SP2SC);
//		hbe_types.push_back(hbe_SP3SC);
//		hbe_types.push_back(hbe_RINGSC);
//
//    for( Size hbe = 1; hbe <= hbe_types.size(); hbe++){
//			HBEvalType hbe_type = hbe_types[hbe];
//			for (Real AHdis = MIN_R; AHdis <MAX_R; AHdis += AHdis_step){
//				for (Real xD = MIN_xD; xD < MAX_xD; xD += xD_step){
//					for (Real xH = MIN_xH; xH < MAX_xH; xH += xH_step){
//						hbond_compute_energy(hbe_type, AHdis, xD, xH, energy, dE_dr, dE_dxD, dE_dxH);
//						if (energy < 0){
//							UT << hbe_type << "\t" << AHdis << "\t" << xD << "\t" << xH << "\t";
//							UT << energy << "\t" << dE_dr << "\t" << dE_dxD << "\t" << dE_dxH << "\n";
//						}
//					}
//				}
//			}
//		}
//	}


	void test_hbond_compute_energy2()
	{
		//To debug at full precision set:
		//UT.precision(16);
		HBondDatabaseCOP database(HBondDatabase::get_database());
		HBondOptions hboptions;

		test::UTracer UT("core/scoring/hbonds/hbond_compute_energy.u");
		bool dummy_chipenalty( false );
		HBGeoDimType AHD_geometric_dimension;
		Real energy, dE_dr, dE_dxD, dE_dxH;
		Size i(0);
		for ( Size ii = 1; ii <= hbdon_MAX; ++ii ) {
			for ( Size jj = 1; jj <= hbacc_MAX; ++jj ) {
				for ( Size kk = 1; kk <= seq_sep_MAX; ++kk ) {
					HBEvalTuple hbt = HBEvalTuple( HBDonChemType(ii), HBAccChemType(jj), HBSeqSep(kk) );
					if ( hbt.eval_type() == 0 ) continue;
					for (Real AHdis = MIN_R; AHdis <MAX_R; AHdis += .8){
						for (Real xD = MIN_xD; xD <= MAX_xD; xD += .5){
							for (Real xH = MIN_xH; xH <= MAX_xH; xH += .5){
								for (Real chi = 0; chi < numeric::constants::d::pi_2; chi += numeric::constants::d::pi_over_2){
									hbond_compute_energy( *database, hboptions, hbt,
										AHdis, xD, xH, chi, energy,
										dummy_chipenalty, AHD_geometric_dimension, dE_dr, dE_dxD, dE_dxH);
									if( energy < 0){
										UT << i << "\t" << ii << "\t" << jj << "\t" << kk << "\t" << hbt.eval_type() << "\t" << AHdis << "\t" << xD << "\t" << xH << "\t" << chi << "\t";
										UT << energy << "\t" << dE_dr << "\t" << dE_dxD << "\t" << dE_dxH << std::endl;
									}
									i++;
								}
							}
						}
					}
				}
			}
		}
	}

	void test_hbond_compute_energy_sp2()
	{
		//To debug at full precision set:
		//std::cout.precision(16);
		HBondDatabaseCOP database(HBondDatabase::get_database("sp2_params"));


		HBondOptions hboptions;
		hboptions.use_sp2_chi_penalty(true);

		test::UTracer UT("core/scoring/hbonds/hbond_compute_energy_sp2.u");
		bool chipenalty;
		HBGeoDimType AHD_geometric_dimension;
		Real energy, dE_dr, dE_dxD, dE_dxH, dE_dBAH, dE_dchi;
		Size i(0);
		HBEvalTuple hbt = HBEvalTuple( hbdon_PBA, hbacc_PBA, seq_sep_other );
		Real xD(.9);
		for (Real AHdis = 1.7; AHdis <= 1.9; AHdis += .2){
			for (Real xH = MIN_xH; xH < MAX_xH; xH += .4){
				for (Real chi = 0.0; chi < numeric::constants::d::pi_2; chi += numeric::constants::d::pi_over_3){
					hbond_compute_energy( *database, hboptions, hbt,
						AHdis, xD, xH, chi, energy,
						chipenalty, AHD_geometric_dimension, dE_dr, dE_dxD, dE_dxH, dE_dBAH, dE_dchi);
					if( energy < 0){
						UT << i << "\t" << hbt.eval_type() << "\t" << AHdis << "\t" << xD << "\t" << xH << "\t" << chi << "\t";
						UT << energy << "\t" << dE_dr << "\t" << dE_dxD << "\t" << dE_dxH << "\t" << dE_dBAH << "\t" << dE_dchi << std::endl;
					}
					i++;
				}
			}
		}
	}

	void test_specific_instance() {

		HBondDatabaseCOP database(HBondDatabase::get_database("sp2_params"));


		HBondOptions hboptions;
		hboptions.use_sp2_chi_penalty(true);

		bool dummy_chipenalty;
		HBGeoDimType AHD_geometric_dimension;
		Real energy, dE_dr, dE_dxD, dE_dxH, dE_dBAH, dE_dchi;
		HBEvalTuple hbt = HBEvalTuple( hbdon_PBA, hbacc_PBA, seq_sep_other );

		test::UTracer UT("core/scoring/hbonds/hbond_compute_energy_sp2_specific_cases.u");
		{
			Real AHdis(1.7), xD(1.0), xH(0.36), chi(0);

			hbond_compute_energy( *database, hboptions, hbt,
				AHdis, xD, xH, chi, energy,
				dummy_chipenalty, AHD_geometric_dimension, dE_dr, dE_dxD, dE_dxH, dE_dBAH, dE_dchi);
			UT << "\t" << hbt.eval_type() << "\t" << AHdis << "\t" << xD << "\t" << xH << "\t" << chi << "\t";
			UT << energy << "\t" << dE_dr << "\t" << dE_dxD << "\t" << dE_dxH << "\t" << dE_dBAH << "\t" << dE_dchi << std::endl;
		}

		{
			Real AHdis(1.7), xD(.9), xH(-0.07), chi(0.03);
			hbond_compute_energy( *database, hboptions, hbt,
				AHdis, xD, xH, chi, energy,
				dummy_chipenalty, AHD_geometric_dimension, dE_dr, dE_dxD, dE_dxH, dE_dBAH, dE_dchi);
			UT << "\t" << hbt.eval_type() << "\t" << AHdis << "\t" << xD << "\t" << xH << "\t" << chi << "\t";
			UT << energy << "\t" << dE_dr << "\t" << dE_dxD << "\t" << dE_dxH << "\t" << dE_dBAH << "\t" << dE_dchi << std::endl;
		}
		{
			Real AHdis(1.7), xD(.9), xH(-0.07), chi(0.04);
			hbond_compute_energy( *database, hboptions, hbt,
				AHdis, xD, xH, chi, energy,
				dummy_chipenalty, AHD_geometric_dimension, dE_dr, dE_dxD, dE_dxH, dE_dBAH, dE_dchi);
			UT << "\t" << hbt.eval_type() << "\t" << AHdis << "\t" << xD << "\t" << xH << "\t" << chi << "\t";
			UT << energy << "\t" << dE_dr << "\t" << dE_dxD << "\t" << dE_dxH << "\t" << dE_dBAH << "\t" << dE_dchi << std::endl;
		}
	}


	void test_hbond_set(){

		HBondDatabaseCOP database(HBondDatabase::get_database());
		HBondOptions hboptions;
		Pose pose = create_trpcage_ideal_pose();
		ScoreFunctionOP score_function( get_score_function() );
		score_function->score( pose );

		HBondSet hbond_set( pose.n_residue() );

		fill_hbond_set(pose,
			true, /* calculate derivative */
			hbond_set);
		
		
		for (Size i = 1; i <= pose.total_residue(); ++i){
		    if (hbond_set.nhbonds(i, false) > 0){
			utility::vector1< HBondCOP > const residue_hbonds = hbond_set.residue_hbonds(i, false /* include all*/);
			TS_ASSERT(residue_hbonds.size() >= 1);
			for (Size k = 1; i <= pose.residue(i).natoms(); ++i){
			    id::AtomID atom = id::AtomID(k, i);
			    if (hbond_set.nhbonds(atom, false) >= 1){
				utility::vector1< HBondCOP > const atom_hbonds = hbond_set.atom_hbonds(atom, false /* include only allowed*/);
				TS_ASSERT(atom_hbonds.size() >= 1);
				break;
			    }
			}
			
			break;
		    }
		}
		
		for( Size i=1; i <= hbond_set.nhbonds(); ++i ){
			const HBond hbond( hbond_set.hbond(i) );
			const Residue acc_res( pose.residue(hbond.acc_res()));
			const Vector Axyz(acc_res.atom(hbond.acc_atm()).xyz());
			const Vector Bxyz(acc_res.atom(acc_res.atom_base(hbond.acc_atm())).xyz());
			const Vector B2xyz(acc_res.atom(acc_res.abase2(hbond.acc_atm())).xyz());
			const Residue don_res( pose.residue(hbond.don_res()));
			const Vector Hxyz(don_res.atom(hbond.don_hatm()).xyz());
			const Vector Dxyz(don_res.atom(don_res.atom_base(hbond.don_hatm())).xyz());

			Real energy;
			HBondDerivs deriv;

			hb_energy_deriv(*database, hboptions,
				HBEvalTuple( don_res.atom_base( hbond.don_hatm()), don_res, hbond.acc_atm(), acc_res ),
				Dxyz, Hxyz, Axyz, Bxyz, B2xyz,
				energy,
				hbderiv_ABE_GO, deriv);

			TS_ASSERT( energy == hbond.energy() );
			TS_ASSERT( deriv.h_deriv.f1() == hbond.derivs().h_deriv.f1() );
			TS_ASSERT( deriv.h_deriv.f2() == hbond.derivs().h_deriv.f2() );
			TS_ASSERT( deriv.acc_deriv.f1() == hbond.derivs().acc_deriv.f1() );
			TS_ASSERT( deriv.acc_deriv.f2() == hbond.derivs().acc_deriv.f2() );
		}
	}



	// A tests that can be run if the hbond types might change.  Go
	// through all chemical possible hydrogen bond donor and acceptors
	// and document non-geometric information about them.

	// Before use:
	// create blank file 'test/core/scoring/hbonds/hbond_non_geometric.u' before building the unit tests.
	// Otherwise the test won't have anything to compare the output against and it will complain

	// This test might exceed the time limit and get killed...

//	void do_for_each_type(
//		HBSeqSep const seq_sep,
//		Residue const & don_rsd,
//		Size const & datm,
//		Residue const & acc_rsd,
//		Size const & aatm,
//		test::UTracer & UT // write output here
//	){
//		HBDonChemType don_chem_type( get_hb_don_chem_type(datm, don_rsd));
//		HBAccChemType acc_chem_type( get_hb_acc_chem_type(aatm, acc_rsd));
//		HBEvalType hbe(HBEval_lookup( don_chem_type, acc_chem_type, seq_sep));
//
//		UT << don_rsd.atom_name(datm) << "\t";
//		UT << don_rsd.name() << "\t";
//		UT << don_chem_type << "\t";
//		UT << acc_rsd.atom_name(aatm) << "\t";
//		UT << acc_rsd.name() << "\t";
//		UT << acc_chem_type << "\t";
//		UT << seq_sep << "\t";
//		UT << hbe << "\t";
//		UT << hbe_is_BB_type( hbe ) << "\t";
//		UT << get_hbond_weight_type( hbe ) << "\t";
//		UT << get_hbe_acc_hybrid( hbe ) << "\n";
//	}
//
//
//
//
//	void dont_test_hbond_non_geometric(){
//		using namespace chemical;
//		test::UTracer UT("core/scoring/hbonds/hbond_non_geometric.u");
//
//		UT << "# Compute non-geometric parameters for hydrogen bonds\n";
//		UT << "# Since the parameters are non-geometric can identify them by non-coordinate information in  residues\n";
//
//		utility::vector1<HBSeqSep> sc_seq_seps;
//		sc_seq_seps.push_back(seq_sep_other);
//
//		utility::vector1<HBSeqSep> bsc_seq_seps;
//		bsc_seq_seps.push_back(seq_sep_other);
//		bsc_seq_seps.push_back(seq_sep_PM1);
//
//		utility::vector1<HBSeqSep> bb_seq_seps;
//		bb_seq_seps.push_back(seq_sep_other);
//		bb_seq_seps.push_back( seq_sep_M4);
//		bb_seq_seps.push_back( seq_sep_M3);
//		bb_seq_seps.push_back( seq_sep_M2);
//		bb_seq_seps.push_back( seq_sep_PM1);
//		bb_seq_seps.push_back( seq_sep_P2);
//		bb_seq_seps.push_back( seq_sep_P3);
//		bb_seq_seps.push_back( seq_sep_P4);
//
//		utility::vector1<ResidueTypeSetCAP> residue_type_sets;
//		residue_type_sets.push_back( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD));
//		residue_type_sets.push_back( ChemicalManager::get_instance()->residue_type_set( "rna" ));
//
//		for( Size i_rset = 1; i_rset != residue_type_sets.size(); ++i_rset){
//			ResidueTypeSetCAP rset(residue_type_sets[i_rset]);
//
//			UT << "\n";
//			UT << "ResidueTypeSet:" << residue_type_sets[i_rset] << "\n";
//			UT << "datm\tdon_rsd\tHBDonChemType\taatm\tacc_rsd\tHBAccChemType\tHBSeqSep\tHBEvalType\tis_BB_hb\tHBondWeighType\tHybridization\n";
//
//
//			// iterate through all donor possible donor and acceptor pairs
//			for( ResidueTypeSet::const_residue_iterator don_rsdtype_it=rset->all_residues_begin(), edon_rsdtype_it=rset->all_residues_end();
//					 don_rsdtype_it!=edon_rsdtype_it; ++don_rsdtype_it ) {
//
//				conformation::Residue don_rsd( *(don_rsdtype_it->second), false /*dummy arg*/);
//				for ( AtomIndices::const_iterator hnum  = don_rsd.Hpos_polar().begin(),	hnume = don_rsd.Hpos_polar().end();
//							hnum != hnume; ++hnum ) {
//					Size const datm(don_rsd.atom_base(*hnum));
//					HBDonChemType don_chem_type( get_hb_don_chem_type(datm, don_rsd));
//
//					for( ResidueTypeSet::const_residue_iterator acc_rsdtype_it=rset->all_residues_begin(), eacc_rsdtype_it=rset->all_residues_end();
//							 acc_rsdtype_it!=eacc_rsdtype_it; ++acc_rsdtype_it ) {
//						conformation::Residue acc_rsd( *(acc_rsdtype_it->second), false /*dummy arg*/);
//						for ( AtomIndices::const_iterator	anum  = acc_rsd.accpt_pos().begin(), anume = acc_rsd.accpt_pos().end();
//									anum != anume; ++anum ) {
//							Size const aatm( *anum );
//							HBAccChemType acc_chem_type( get_hb_acc_chem_type(aatm, acc_rsd));
//
//							if (don_chem_type == hbdon_PBA && acc_chem_type == hbacc_PBA){
//								for( Size seq_sep = 1; seq_sep <= bb_seq_seps.size(); ++seq_sep){
//									do_for_each_type(bb_seq_seps[seq_sep], don_rsd, datm, acc_rsd, aatm, UT);
//								}
//							} else if (don_chem_type == hbdon_PBA || acc_chem_type == hbacc_PBA){
//								for( Size seq_sep = 1; seq_sep <= bsc_seq_seps.size(); ++seq_sep){
//									do_for_each_type(bsc_seq_seps[seq_sep], don_rsd, datm, acc_rsd, aatm, UT);
//								}
//							}else{
//								do_for_each_type(seq_sep_other, don_rsd, datm, acc_rsd, aatm, UT);
//							}
// 						}
// 					}
// 				}
//			}
//		}
//	}

	void test_hbond_type_manager(){

		test::UTracer UT("core/scoring/hbonds/hbond_type_manager.u");
		UT << "Hydrogen Bond Types"<< std::endl << std::endl;

		UT << "HBondWeightType:" << std::endl;
		for (Size hbw=1; hbw <= hbw_MAX; ++hbw){
			std::string hbw_name(HBondTypeManager::name_from_weight_type(HBondWeightType(hbw)));
			UT << "\t" << hbw_name << std::endl;
			TS_ASSERT(HBondTypeManager::weight_type_from_name(hbw_name) == HBondWeightType(hbw));
			TS_ASSERT(HBondTypeManager::is_weight_type(hbw_name));
		}
		UT << std::endl;

		UT << "HBDerivType:" << std::endl;
		for (Size hbderiv=1; hbderiv <= hbderiv_MAX; ++hbderiv){
			std::string hbderiv_name(HBondTypeManager::name_from_deriv_type(HBDerivType(hbderiv)));
			UT << "\t" << hbderiv_name << std::endl;
			TS_ASSERT(HBondTypeManager::deriv_type_from_name(hbderiv_name) == HBDerivType(hbderiv));
			TS_ASSERT(HBondTypeManager::is_deriv_type(hbderiv_name));
		}
		UT << std::endl;

		UT << "HBDonChemType:" << std::endl;
		for (Size hbdon=1; hbdon <= hbdon_MAX; ++hbdon){
			std::string hbdon_name(HBondTypeManager::name_from_don_chem_type(HBDonChemType(hbdon)));
			UT << "\t" << hbdon_name << std::endl;
			TS_ASSERT(HBondTypeManager::don_chem_type_from_name(hbdon_name) == HBDonChemType(hbdon));
			TS_ASSERT(HBondTypeManager::is_don_chem_type(hbdon_name));
		}
		UT << std::endl;

		UT << "HBAccChemType:" << std::endl;
		for (Size hbacc=1; hbacc <= hbacc_MAX; ++hbacc){
			std::string hbacc_name(HBondTypeManager::name_from_acc_chem_type(HBAccChemType(hbacc)));
			UT << "\t" << hbacc_name << std::endl;
			TS_ASSERT(HBondTypeManager::acc_chem_type_from_name(hbacc_name) == HBAccChemType(hbacc));
			TS_ASSERT(HBondTypeManager::is_acc_chem_type(hbacc_name));
		}
		UT << std::endl;

		UT << "HBSeqSep:" << std::endl;
		for (Size hbseq_sep=1; hbseq_sep <= seq_sep_MAX; ++hbseq_sep){
			std::string hbseq_sep_name(HBondTypeManager::name_from_seq_sep_type(HBSeqSep(hbseq_sep)));
			UT << "\t" << hbseq_sep_name << std::endl;
			TS_ASSERT(HBondTypeManager::seq_sep_type_from_name(hbseq_sep_name) == HBSeqSep(hbseq_sep));
			TS_ASSERT(HBondTypeManager::is_seq_sep_type(hbseq_sep_name));
		}
		UT << std::endl;

		UT << "HBGeoDimType:" << std::endl;
		for (Size hbgd=1; hbgd <= hbgd_MAX; ++hbgd){
			std::string hbgd_name(HBondTypeManager::name_from_geo_dim_type(HBGeoDimType(hbgd)));
			UT << "\t" << hbgd_name << std::endl;
			TS_ASSERT(HBondTypeManager::geo_dim_type_from_name(hbgd_name) == HBGeoDimType(hbgd));
			TS_ASSERT(HBondTypeManager::is_geo_dim_type(hbgd_name));
		}
		UT << std::endl;

	}


private:

	//Real TOLERATED_ERROR;
};
