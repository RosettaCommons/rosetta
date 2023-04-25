// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// C/C++ headers
#include <map>
#include <string>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/lddt.hh>
#include <basic/Tracer.hh>


static basic::Tracer TR("test.core.scoring.lddt.cxxtest");

class lDDTTest : public CxxTest::TestSuite {
public:

	core::pose::PoseCOP ref_;
	core::pose::PoseCOP mod_;

	void setUp() {
		core_init_with_additional_options(
			"-extra_res_fa core/scoring/1pqc.params"
		);

		ref_ = core::import_pose::pose_from_file("core/scoring/4dO8B.pdb", core::import_pose::PDB_file);
		mod_ = core::import_pose::pose_from_file("core/scoring/model-5.pdb", core::import_pose::PDB_file);
	}

	void test_global_lddt() {
		TR << "test_global_lddt()" << std::endl;

		core::scoring::lDDT_Calculator calc;
		core::Real global_lddt = calc( *ref_, *mod_ );

		// Value from the https://swissmodel.expasy.org/lddt/ site
		TS_ASSERT_DELTA(0.8171, global_lddt, 0.0001);
	}

	void test_residue_lddt() {
		TR << "test_residue_lddt()" << std::endl;

		// Values from the https://swissmodel.expasy.org/lddt/ site
		utility::vector1< core::Real > ground_truth = {
			0.891233,
			0.886029,
			0.900916,
			0.875492,
			0.740529,
			0.546237,
			0.392538,
			0.377485,
			0.49878,
			0.33283 ,
			0.510046,
			0.528785,
			0.778052,
			0.847431,
			0.73371 ,
			0.898645,
			0.922372,
			0.929828,
			0.870083,
			0.941716,
			0.919282,
			0.863006,
			0.924087,
			0.920779,
			0.892754,
			0.889142,
			0.742158,
			0.843124,
			0.629973,
			0.908732,
			0.883245,
			0.912879,
			0.887409,
			0.826923,
			0.86844 ,
			0.859364,
			0.825912,
			0.862894,
			0.84984 ,
			0.863383,
			0.891406,
			0.898077,
			0.921585,
			0.923265,
			0.92232 ,
			0.919221,
			0.893962,
			0.723097,
			0.924697,
			0.803685,
			0.822586,
			0.913801,
			0.831568,
			0.815603,
			0.909329,
			0.901797,
			0.595484,
			0.927887,
			0.836212,
			0.771011,
			0.804084,
			0.799305,
			0.779029,
			0.806148,
			0.867251,
			0.719833
			};

		core::scoring::lDDT_Calculator calc;
		utility::vector1< core::Real > residue_lddt = calc.residue_lDDT(*ref_, *mod_);

		for ( core::Size ii(1); ii <= ref_->size(); ++ii ) {
			TSM_ASSERT_DELTA("Residue "+std::to_string(ii), ground_truth[ii], residue_lddt[ii], 0.0001);
		}
	}


	void check_lddt_for_subsection(
		core::pose::Pose const & ref,
		core::pose::Pose const & mod,
		core::Size start, // Start and end are inclusive
		core::Size end,
		utility::vector1< core::Real > const & ground_truth
	) {
		std::map< core::Size, core::Size > res_map;
		for ( core::Size ii(start); ii <= end; ++ii ) {
			res_map[ii] = ii;
		}

		core::scoring::lDDT_Calculator calc;

		std::map< core::Size, core::Real > residue_lddt;
		calc.residue_lDDT( ref, mod, res_map, residue_lddt );

		for ( core::Size ii(start); ii <= end; ++ii ) {
			TSM_ASSERT_DELTA("Residue "+std::to_string(ii), ground_truth[ii-(start-1)], residue_lddt[ii], 0.0001);
		}
	}

	void no_test_residue_lddt_noalt_subsection() { // Disabled for time -- global and DE subsection will test
		TR << "test_lddt_noalt_subsection()" << std::endl;
		/// Test a subsection of the protein without any chemical alternatives.

		// Values from the https://swissmodel.expasy.org/lddt/ site, calculated on a trimmed model & reference.
		utility::vector1< core::Real > ground_truth = {
			0.843023, // (725/860)
			0.712644, // (744/1044)
			0.915899, // (795/868)
			0.817972, // (710/868)
			0.777311 // (740/952)
			};

		check_lddt_for_subsection(*ref_, *mod_, 47, 51, ground_truth);
	}

	void no_test_residue_lddt_VL_subsection() { // Disabled for time -- global and DE subsection will test
		TR << "test_lddt_VL_subsection()" << std::endl;
		/// Test a subsection of the protein with V/L

		// Values from the https://swissmodel.expasy.org/lddt/ site, calculated on a trimmed model & reference.
		utility::vector1< core::Real > ground_truth = {
			0.931915, // (876/940)
			0.833333, // (910/1092)
			0.83125 , // (931/1120)
			0.72252 , // (1078/1492)
			0.703378, // (1041/1480)
			0.673437, // (862/1280)
			0.573155, // (901/1572)
			0.660417, // (634/960)
			0.494469  // (447/904)
			};

		check_lddt_for_subsection(*ref_, *mod_, 1, 9, ground_truth);
	}

	void test_residue_lddt_DE_subsection() {
		TR << "test_lddt_DE_subsection()" << std::endl;
		/// Test a subsection of the protein with D/E

		// Values from the https://swissmodel.expasy.org/lddt/ site, calculated on a trimmed model & reference.
		utility::vector1< core::Real > ground_truth = {
			0.87936 , // (1210/1376)
			0.831301, // (1227/1476)
			0.7881  , // (1510/1916)
			0.781103, // (1331/1704)
			0.813158, // (1545/1900)
			0.85198 , // (1721/2020)
			0.84543 , // (1258/1488)
			0.877349, // (1681/1916)
			0.868449  // (1657/1908)
			};

		check_lddt_for_subsection(*ref_, *mod_, 58, 66, ground_truth);
	}

	void no_test_calpha() { // Disabled for time -- test_backbone() will test predicates
		TR << "test_calpha()" << std::endl;

		core::scoring::lDDT_Calculator calc( utility::pointer::make_shared< core::scoring::IsProteinCAPredicate >() );

		// Don't have a gold standard, so we're just checking that with & without alt states are the same.
		core::Real with_alt = calc( *ref_, *mod_ );

		calc.consider_alt_states( false );
		core::Real without_alt = calc( *ref_, *mod_ );

		TS_ASSERT_EQUALS( with_alt, without_alt ); // Exactly equal, despite real -- this should be the same interger/integer division.
	}

	void test_backbone() {
		TR << "test_backbone()" << std::endl;

		core::scoring::lDDT_Calculator calc( utility::pointer::make_shared< core::scoring::IsProteinBackboneIncludingOPredicate >() );
		calc.consider_alt_states( false ); // For speed

		core::Real global = calc( *ref_, *mod_ );

		// Gold Standard from the https://swissmodel.expasy.org/lddt/ site, for an all-GLY version.
		TS_ASSERT_DELTA( global, 0.9130, 0.0001);

		TR << "test_backbone() -- residue" << std::endl;

		utility::vector1< core::Real > ground_truth = {
			0.958879, // (2052/2140)
			0.944676, // (2271/2404)
			0.939985, // (2553/2716)
			0.927878, // (2676/2884)
			0.829359, // (2226/2684)
			0.582737, // (1303/2236)
			0.51129 , // (951/1860)
			0.392857, // (550/1400)
			0.510514, // (437/856)
			0.471354, // (362/768)
			0.511194, // (685/1340)
			0.664678, // (1114/1676)
			0.859345, // (1784/2076)
			0.892665, // (1996/2236)
			0.953846, // (2232/2340)
			0.980354, // (1996/2036)
			0.97998 , // (2007/2048)
			0.987875, // (1711/1732)
			0.958333, // (1242/1296)
			0.990183, // (1513/1528)
			0.952495, // (1985/2084)
			0.982759, // (2166/2204)
			0.985791, // (2567/2604)
			0.978112, // (2860/2924)
			0.958696, // (3087/3220)
			0.956755, // (3031/3168)
			0.945442, // (2738/2896)
			0.953333, // (2288/2400)
			0.955567, // (1785/1868)
			0.993802, // (1443/1452)
			0.996552, // (1156/1160)
			1       , // (992/992)
			0.952273, // (838/880)
			0.953196, // (835/876)
			0.985965, // (1124/1140)
			0.939759, // (1560/1660)
			0.90392 , // (1891/2092)
			0.888681, // (2371/2668)
			0.907672, // (2792/3076)
			0.915129, // (2976/3252)
			0.941307, // (2967/3152)
			0.954993, // (2716/2844)
			0.975   , // (2535/2600)
			0.972378, // (2077/2136)
			0.973611, // (2103/2160)
			0.965485, // (2070/2144)
			0.966216, // (2145/2220)
			0.989201, // (1832/1852)
			0.994885, // (1556/1564)
			0.946013, // (1139/1204)
			0.973032, // (1335/1372)
			0.980205, // (1337/1364)
			0.990554, // (1573/1588)
			0.983004, // (1793/1824)
			0.966964, // (2166/2240)
			0.968905, // (2337/2412)
			0.967344, // (2666/2756)
			0.972358, // (2392/2460)
			0.969298, // (2431/2508)
			0.906643, // (2020/2228)
			0.915636, // (2073/2264)
			0.889815, // (1922/2160)
			0.877295, // (2102/2396)
			0.915748, // (2326/2540)
			0.855609, // (2471/2888)
			0.799073, // (2068/2588)
			};

		utility::vector1< core::Real > residue_lddt = calc.residue_lDDT(*ref_, *mod_);

		for ( core::Size ii(1); ii <= ref_->size(); ++ii ) {
			TSM_ASSERT_DELTA("Residue "+std::to_string(ii), ground_truth[ii], residue_lddt[ii], 0.0001);
		}

		TR << "test_backbone() -- alt states" << std::endl;

		// Check that the alt state calcuation also comes up with the proper value
		calc.consider_alt_states( true );
		TS_ASSERT_DELTA( global, 0.9130, 0.0001);
	}

	void swap_coords(core::pose::Pose & pose, core::Size resi, std::string const & atom1, std::string const & atom2 ) {
		core::id::AtomID const id1( pose.residue(resi).atom_index( atom1 ), resi );
		core::id::AtomID const id2( pose.residue(resi).atom_index( atom2 ), resi );

		core::Vector tmp = pose.xyz( id1 );
		pose.set_xyz( id1, pose.xyz( id2 ) );
		pose.set_xyz( id2, tmp );
	}

	void test_ligand() {
		TR << "test_ligand()" << std::endl;

		core::pose::PoseOP ref = core::import_pose::pose_from_file("core/scoring/1pqc_0001.pdb", core::import_pose::PDB_file);
		core::pose::PoseOP mod = core::import_pose::pose_from_file("core/scoring/1pqc_0003.pdb", core::import_pose::PDB_file);

		core::scoring::lDDT_Calculator calc;
		calc.seqsep( 0 ); // This is a single-residue pose, so we need to include the self to get any entries.
		calc.consider_alt_states( true );

		core::Real const with_alt = calc( *ref, *mod );

		TS_ASSERT_DIFFERS( with_alt, 0.0 );

		calc.consider_alt_states( false );

		core::Real const base = calc( *ref, *mod );
		TS_ASSERT_DIFFERS( base, with_alt ); // with_alt should be lower

		// Now we swap around coordinates
		swap_coords( *mod, 1, "F3", "F2" );
		swap_coords( *mod, 1, "C6", "C7" );
		swap_coords( *mod, 1, "C4", "C5" );
		swap_coords( *mod, 1, "C16", "C17" );
		swap_coords( *mod, 1, "F9", "F4" ); // F9 -> F4
		swap_coords( *mod, 1, "F7", "F6" ); // F7 -> F6
		swap_coords( *mod, 1, "F8", "F5" ); // F8 -> F5
		swap_coords( *mod, 1, "F7", "F9" ); // F6 -> F7 -> F9
		// F4 -> F9 -> F7 ( Should be F8 )
		// F5 -> F8 ( Should be F7 )
		swap_coords( *mod, 1, "F7", "F8" );

		core::Real const swapped = calc( *ref, *mod );
		TS_ASSERT_DELTA( with_alt, swapped, 0.0001 );
	}

	void test_large_protein() {
		TR << "test_large_protein()" << std::endl;

		// This is chosen mainly for being a 3156 residue protein.
		// This test is mainly to confirm that we're not horribly slow with dealing with large proteins.
		core::pose::PoseOP mod = core::import_pose::pose_from_file("core/pose/2WDQ__tr.pdb", core::import_pose::PDB_file);

		core::scoring::lDDT_Calculator calc;

		calc.consider_alt_states( true );

		//utility::vector1< core::Real > residue_lddt = calc.residue_lDDT( *mod, *mod );

		std::map< core::Size, core::Real > residue_lddt;
		std::map< core::id::AtomID, core::Real > atom_lddt;

		TR << "STARTING CALC FOR LARGE PROTEIN" << std::endl;

		core::Real const with_alt = calc.all_lDDT( *mod, *mod, residue_lddt, atom_lddt );

		TR << " NUMBER OF RESIDUES " << residue_lddt.size() << std::endl;
		TR << " NUMBER OF ATOMS " << atom_lddt.size() << std::endl;

		TS_ASSERT_EQUALS( with_alt, 1.0 );
	}

};
