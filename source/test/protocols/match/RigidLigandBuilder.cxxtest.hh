// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/RigidLigandBuilderTests.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/match/downstream/RigidLigandBuilder.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.hh>
#include <protocols/match/downstream/ClassicMatchAlgorithm.hh>
#include <protocols/match/upstream/ProteinSCSampler.hh>
#include <protocols/match/upstream/ProteinUpstreamBuilder.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>
#include <protocols/match/output/WriteUpstreamCoordinateKineamge.hh>

#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>

// C++ headers
#include <string>
#include <iostream>
#include <sstream>

//Auto Headers
#include <core/id/AtomID.hh>
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>
#include <boost/unordered/unordered_map_fwd.hpp>


using namespace protocols::match;
using namespace protocols::match::downstream;
using namespace protocols::match::output;
using namespace protocols::match::upstream;
using namespace protocols::toolbox::match_enzdes_util;


// --------------- Test Class --------------- //

class RigidLigandBuilderTests : public CxxTest::TestSuite {

public:

	typedef core::Real Real;
	typedef core::Size Size;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	utility::vector1< Real > d_1, ang_U2D1_1, tor_U3D1_1, ang_U1D2_1, tor_U1D3_1, tor_U2D2_1;
	utility::vector1< Real > d_2, ang_U2D1_2, tor_U3D1_2, ang_U1D2_2, tor_U1D3_2, tor_U2D2_2;
	utility::vector1< Real > d_3, ang_U2D1_3, tor_U3D1_3, ang_U1D2_3, tor_U1D3_3, tor_U2D2_3;

	// Shared initialization goes here.
	void setUp() {
		using namespace core::chemical;

		core_init();

		ResidueTypeSet & restype_set(
			ChemicalManager::get_instance()->nonconst_residue_type_set( FA_STANDARD ));
		if ( ! restype_set.has_name( "CARBARYL_LG1" ) ) {
			utility::vector1< std::string > carbaryl_list;
			carbaryl_list.push_back( "protocols/match/carbaryl_LG1.params" );
			restype_set.read_files_for_custom_residue_types(carbaryl_list);
		}

		if ( ! restype_set.has_name( "MBH_LG1" ) ) {
			utility::vector1< std::string > carbaryl_list;
			carbaryl_list.push_back( "protocols/match/MBH_LG.params" );
			restype_set.read_files_for_custom_residue_types(carbaryl_list);
		}


		d_1.clear(); ang_U2D1_1.clear(); tor_U3D1_1.clear(); ang_U1D2_1.clear(); tor_U1D3_1.clear(); tor_U2D2_1.clear();
		d_2.clear(); ang_U2D1_2.clear(); tor_U3D1_2.clear(); ang_U1D2_2.clear(); tor_U1D3_2.clear(); tor_U2D2_2.clear();
		d_3.clear(); ang_U2D1_3.clear(); tor_U3D1_3.clear(); ang_U1D2_3.clear(); tor_U1D3_3.clear(); tor_U2D2_3.clear();

	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_rigid_ligand_ds_builder_recover_input_coords() {
		using namespace core;
		using namespace core::chemical;
		using namespace core::io::pdb;
		using namespace core::pose;

		Pose carbaryl_pose;
		core::import_pose::pose_from_file( carbaryl_pose,
			*(ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )),
			"protocols/match/E1cb_carbaryl_1his_oxy_1bb_10_2.pdb", core::import_pose::PDB_file );

		std::string const  at4( "C12" ),  at5( "C9" ),   at6( "C13" );
		std::string const oat1( "C17" ), oat2( "C19" ), oat3( "C22" ); // orientation atoms, oats

		Size const at4id( carbaryl_pose.residue(1).atom_index( at4 ) );
		Size const at5id( carbaryl_pose.residue(1).atom_index( at5 ) );
		Size const at6id( carbaryl_pose.residue(1).atom_index( at6 ) );

		Size const oat1id( carbaryl_pose.residue(1).atom_index( oat1 ) );
		Size const oat2id( carbaryl_pose.residue(1).atom_index( oat2 ) );
		Size const oat3id( carbaryl_pose.residue(1).atom_index( oat3 ) );

		RigidLigandBuilderOP rigid_builder( new RigidLigandBuilder );
		rigid_builder->initialize_from_residue(
			at4id, at5id, at6id, oat1id, oat2id, oat3id, carbaryl_pose.residue(1) );

		numeric::HomogeneousTransform< Real > oat3frame(
			carbaryl_pose.residue(1).xyz( oat1id ),
			carbaryl_pose.residue(1).xyz( oat2id ),
			carbaryl_pose.residue(1).xyz( oat3id ) );

		Vector euler = oat3frame.euler_angles_deg();

		Hit hit;
		hit.first()[ 1 ] = 1;
		hit.first()[ 2 ] = 1;
		hit.first()[ 3 ] = 1;
		hit.second()[ 1 ] = carbaryl_pose.residue(1).xyz( oat3id ).x();
		hit.second()[ 2 ] = carbaryl_pose.residue(1).xyz( oat3id ).y();
		hit.second()[ 3 ] = carbaryl_pose.residue(1).xyz( oat3id ).z();
		hit.second()[ 4 ] = euler( 1 );
		hit.second()[ 5 ] = euler( 2 );
		hit.second()[ 6 ] = euler( 3 );

		utility::vector1< core::id::AtomID > atomids( carbaryl_pose.residue(1).natoms() );
		for ( Size ii = 1; ii <= carbaryl_pose.residue(1).natoms(); ++ii ) atomids[ ii ] = core::id::AtomID( ii, 1 );
		utility::vector1< Vector > coords( carbaryl_pose.residue(1).natoms()  );

		rigid_builder->coordinates_from_hit( hit, atomids, coords );
		for ( Size ii = 1; ii <= carbaryl_pose.residue(1).natoms(); ++ii ) {
			TS_ASSERT( coords[ ii ].distance( carbaryl_pose.residue(1).xyz( ii ) ) < 1e-6 );
			//std::cout << coords[ ii ].x() << " " << coords[ ii ].y() << " " << coords[ ii ].z() << " vs ";
			//std::cout << carbaryl_pose.residue(1).atom_name( ii ) << " " << carbaryl_pose.residue(1).xyz( ii ).x() << " " << carbaryl_pose.residue(1).xyz( ii ).y() << " " << carbaryl_pose.residue(1).xyz( ii ).z() << std::endl;
		}

		//SingleDownstreamResidueWriter writer;
		//writer.set_restype( & (carbaryl_pose.residue(1).type()) );
		//writer.set_downstream_builder( rigid_builder );
		//writer.write_downstream_coordinates( hit, std::cout );

		utility::vector1< Vector > rotated_coords( carbaryl_pose.residue(1).natoms() );
		numeric::xyzMatrix< double > zrot = numeric::z_rotation_matrix_degrees( 30 );
		for ( Size ii = 1; ii <= carbaryl_pose.residue(1).natoms(); ++ii ) rotated_coords[ ii ] = zrot * carbaryl_pose.residue(1).xyz( ii );

		numeric::HomogeneousTransform< Real > oat3frame2(
			rotated_coords[ oat1id ],
			rotated_coords[ oat2id ],
			rotated_coords[ oat3id ] );

		Vector euler2 = oat3frame2.euler_angles_deg();

		Hit hit2;
		hit2.first()[ 1 ] = 1;
		hit2.first()[ 2 ] = 1;
		hit2.first()[ 3 ] = 1;
		hit2.second()[ 1 ] = rotated_coords[ oat3id ].x();
		hit2.second()[ 2 ] = rotated_coords[ oat3id ].y();
		hit2.second()[ 3 ] = rotated_coords[ oat3id ].z();
		hit2.second()[ 4 ] = euler2( 1 );
		hit2.second()[ 5 ] = euler2( 2 );
		hit2.second()[ 6 ] = euler2( 3 );

		rigid_builder->coordinates_from_hit( hit2, atomids, coords );
		for ( Size ii = 1; ii <= carbaryl_pose.residue(1).natoms(); ++ii ) {
			TS_ASSERT( coords[ ii ].distance( rotated_coords[ ii ] ) < 1e-6 );
			//std::cout << coords[ ii ].x() << " " << coords[ ii ].y() << " " << coords[ ii ].z() << " vs ";
			//std::cout << carbaryl_pose.residue(1).atom_name( ii ) << " " << rotated_coords[ ii ].x() << " " << rotated_coords[ ii ].y() << " " << rotated_coords[ ii ].z() << std::endl;
		}

	}

	void test_lig_downstream_geom() {

		using namespace core;
		using namespace core::chemical;
		using namespace core::io::pdb;
		using namespace core::pose;

		// TEMP TEMP TEMP
		// SHORT CIRCUIT
		return;

		d_1.push_back( 3.2 );
		d_2.push_back( 2.6 );
		d_3.push_back( 2.75 );

		ang_U2D1_1.push_back(  90.0 );
		ang_U2D1_2.push_back( 110.0 ); ang_U2D1_2.push_back( 120.0 ); ang_U2D1_2.push_back( 130.0 );
		ang_U2D1_3.push_back( 180.0 ); ang_U2D1_3.push_back( 160.0 ); ang_U2D1_3.push_back( 140.0 ); ang_U2D1_3.push_back( 120.0 );

		tor_U3D1_1.push_back( 90.0  );
		tor_U3D1_2.push_back( 180.0 ); tor_U3D1_2.push_back( 170.0 ); tor_U3D1_2.push_back( 190.0 );
		tor_U3D1_3.push_back( 0.0   ); tor_U3D1_3.push_back( 180.0 );

		ang_U1D2_1.push_back( 100.0 ); ang_U1D2_1.push_back( 110.0 ); ang_U1D2_1.push_back( 90.0 );
		ang_U1D2_2.push_back( 120.0 );
		ang_U1D2_3.push_back( 120.0 );

		tor_U1D3_1.push_back( 0.0 ); tor_U1D3_1.push_back( 30.0 ); tor_U1D3_1.push_back( 60.0 ); tor_U1D3_1.push_back( 90.0 ); tor_U1D3_1.push_back( 120.0 ); tor_U1D3_1.push_back( 150.0 );
		tor_U1D3_1.push_back(180.0); tor_U1D3_1.push_back(-30.0 ); tor_U1D3_1.push_back(-60.0 ); tor_U1D3_1.push_back(-90.0 ); tor_U1D3_1.push_back(-120.0 ); tor_U1D3_1.push_back(-150.0 );
		tor_U1D3_2.push_back( 90.0 ); tor_U1D3_2.push_back( 270.0 );
		tor_U1D3_3.push_back( 0.0 ); tor_U1D3_3.push_back( 30.0 ); tor_U1D3_3.push_back( 60.0 ); tor_U1D3_3.push_back( 90.0 ); tor_U1D3_3.push_back( 120.0 ); tor_U1D3_3.push_back( 150.0 );
		tor_U1D3_3.push_back(180.0); tor_U1D3_3.push_back(-30.0 ); tor_U1D3_3.push_back(-60.0 ); tor_U1D3_3.push_back(-90.0 ); tor_U1D3_3.push_back(-120.0 ); tor_U1D3_3.push_back(-150.0 );
		tor_U1D3_3.push_back( 0.0  + 15.0); tor_U1D3_3.push_back( 30.0  + 15.0); tor_U1D3_3.push_back( 60.0  + 15.0); tor_U1D3_3.push_back( 90.0  + 15.0); tor_U1D3_3.push_back( 120.0  + 15.0); tor_U1D3_3.push_back( 150.0  + 15.0);
		tor_U1D3_3.push_back(180.0 + 15.0); tor_U1D3_3.push_back(-30.0  + 15.0); tor_U1D3_3.push_back(-60.0  + 15.0); tor_U1D3_3.push_back(-90.0  + 15.0); tor_U1D3_3.push_back(-120.0  + 15.0); tor_U1D3_3.push_back(-150.0  + 15.0);


		tor_U2D2_1.push_back( 90.0 );
		tor_U2D2_2.push_back( 180.0 );
		tor_U2D2_3.push_back( 180.0 ); tor_U2D2_3.push_back( 200.0 ); tor_U2D2_3.push_back( 160.0 );

		ResidueTypeSet & restype_set(
			ChemicalManager::get_instance()->nonconst_residue_type_set( FA_STANDARD ));

		ResidueTypeCOP cys_restype( ResidueTypeOP( new ResidueType( restype_set.name_map( "CYS" ))));

		Pose trpcage = create_trpcage_ideal_pose();
		Pose carbaryl_pose;
		core::import_pose::pose_from_file( carbaryl_pose,
			*(ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )),
			"protocols/match/E1cb_carbaryl_1his_oxy_1bb_10_2.pdb", core::import_pose::PDB_file );

		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		BuildSet build_set;
		build_set.set_residue_type( cys_restype );

		SampleStrategyData strat; strat.set_strategy( rotameric_chi_mimic_EX_flags );
		strat.set_sample_level( core::pack::task::EX_ONE_STDDEV );

		build_set.set_sample_strategy_for_chi( 1, strat );
		build_set.set_sample_strategy_for_chi( 2, strat );

		//std::ostringstream sout;
		std::ofstream fout;
		fout.open("hits.kin");

		ExternalGeomSampler exsampler;
		exsampler.set_tor_U3D1_samples(     tor_U3D1_1    );
		exsampler.set_dis_U1D1_samples(      d_1     );
		exsampler.set_ang_U2D1_samples( ang_U2D1_1 );
		exsampler.set_ang_U1D2_samples( ang_U1D2_1 );
		exsampler.set_tor_U2D2_samples(    tor_U2D2_1   );
		exsampler.set_tor_U1D3_samples(     tor_U1D3_1    );


		std::string const  at4( "C12" ),  at5( "C9" ),   at6( "C13" );
		std::string const oat1( "C17" ), oat2( "C19" ), oat3( "C22" ); // orientation atoms, oats

		Size const at4id( carbaryl_pose.residue(1).atom_index( at4 ) );
		Size const at5id( carbaryl_pose.residue(1).atom_index( at5 ) );
		Size const at6id( carbaryl_pose.residue(1).atom_index( at6 ) );

		Size const oat1id( carbaryl_pose.residue(1).atom_index( oat1 ) );
		Size const oat2id( carbaryl_pose.residue(1).atom_index( oat2 ) );
		Size const oat3id( carbaryl_pose.residue(1).atom_index( oat3 ) );

		RigidLigandBuilderOP rigid_builder( new RigidLigandBuilder );
		rigid_builder->initialize_from_residue(
			at4id, at5id, at6id, oat1id, oat2id, oat3id, carbaryl_pose.residue(1) );
		rigid_builder->initialize_upstream_residue( cys_restype ); /// do not initialize bonded data.

		BumpGridOP bb_grid = bump_grid_to_enclose_pose( trpcage );
		for ( Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			BumpGridOP resbgop = bump_grid_to_enclose_residue_backbone( trpcage.residue( ii ), *bb_grid );
			fill_grid_with_backbone_heavyatom_spheres( trpcage.residue( ii ), *resbgop );
			bb_grid->or_with( *resbgop );
		}

		ClassicMatchAlgorithmOP match_algorithm( new ClassicMatchAlgorithm( 1 ) );
		match_algorithm->set_residue_type( cys_restype );
		match_algorithm->add_external_geom_sampler(
			exsampler, 1, "CA", "CB", "SG", rigid_builder );
		WriteUpstreamCoordinateKinemageOP dsalgorithm( new WriteUpstreamCoordinateKinemage( fout ) );
		dsalgorithm->set_match_algorithm( match_algorithm );
		dsalgorithm->set_n_downstream_to_output( 5 );
		match_algorithm->set_bb_grid( bb_grid );
		rigid_builder->set_bb_grid( bb_grid );

		SingleDownstreamResidueWriterOP downstream_writer( new SingleDownstreamResidueWriter );
		downstream_writer->set_restype( carbaryl_pose.residue(1).type().get_self_ptr() );
		downstream_writer->set_downstream_builder( rigid_builder );
		downstream_writer->set_downstream_master( "carbaryl" );

		dsalgorithm->set_downstream_writer( downstream_writer );

		build_set.set_downstream_algorithm( dsalgorithm );

		ProteinUpstreamBuilder scbuilder;
		scbuilder.add_build_set( build_set );
		scbuilder.set_sampler( ProteinSCSamplerCOP( new DunbrackSCSampler ) );

		scbuilder.build( *res2bp );

		fout.close();

	}

	void test_downstream_placement_from_sinisas_param_file() {
		using namespace core;
		using namespace core::chemical;
		using namespace core::io::pdb;
		using namespace core::pose;

		//d          2.25  0.25  5.00     1
		//angleA    112.0 30.00  0.50     2
		//torsionA  60.00 20.00  0.01     2
		//angleB    103.0 30.00  0.01     2
		//torsionAB 60.00 120.0  0.01     1
		//torsionB   60.0 120.0  0.01     1

		d_1.push_back( 2.0 ); d_1.push_back( 2.25 ); d_1.push_back( 2.5 );
		ang_U2D1_1.push_back(  52 ); ang_U2D1_1.push_back(  82 ); ang_U2D1_1.push_back( 112 ); ang_U2D1_1.push_back(  142 ); ang_U2D1_1.push_back( 172 );
		tor_U3D1_1.push_back( 20.0  ); tor_U3D1_1.push_back( 40.0  );tor_U3D1_1.push_back( 60.0  );tor_U3D1_1.push_back( 80.0  );tor_U3D1_1.push_back( 100.0  );
		ang_U1D2_1.push_back( 43.0 ); ang_U1D2_1.push_back( 73.0 ); ang_U1D2_1.push_back( 103.0 ); ang_U1D2_1.push_back( 133.0 ); ang_U1D2_1.push_back( 163.0 );
		tor_U2D2_1.push_back( -60.0 ); tor_U2D2_1.push_back( 60.0 ); tor_U2D2_1.push_back( 180.0 );
		tor_U1D3_1.push_back( -60.0 ); tor_U1D3_1.push_back( 60.0 ); tor_U1D3_1.push_back( 180.0 );

		ResidueTypeSet & restype_set(
			ChemicalManager::get_instance()->nonconst_residue_type_set( FA_STANDARD ));

		ResidueTypeCOP cys_restype( ResidueTypeOP( new ResidueType( restype_set.name_map( "CYS" ))));

		Pose trpcage = create_trpcage_ideal_pose();
		Pose mbh_pose;
		core::import_pose::pose_from_file( mbh_pose,
			*(ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )),
			"protocols/match/mbh_lg_0001.pdb", core::import_pose::PDB_file );

		BuildSet build_set;
		build_set.set_residue_type( cys_restype );

		SampleStrategyData strat; strat.set_strategy( rotameric_chi_mimic_EX_flags );
		strat.set_sample_level( core::pack::task::EX_ONE_STDDEV );

		build_set.set_sample_strategy_for_chi( 1, strat );
		build_set.set_sample_strategy_for_chi( 2, strat );

		//std::ostringstream sout;
		std::ofstream fout;
		fout.open("hits_mbh.kin");

		ExternalGeomSampler exsampler;
		exsampler.set_tor_U3D1_samples( tor_U3D1_1 );
		exsampler.set_dis_U1D1_samples(      d_1   );
		exsampler.set_ang_U2D1_samples( ang_U2D1_1 );
		exsampler.set_ang_U1D2_samples( ang_U1D2_1 );
		exsampler.set_tor_U2D2_samples( tor_U2D2_1 );
		exsampler.set_tor_U1D3_samples( tor_U1D3_1 );


		std::string const  at4( "C1" ),  at5( "C5" ),   at6( "C6" );
		std::string const oat1( "C20" ), oat2( "C19" ), oat3( "C18" ); // orientation atoms, oats

		Size const at4id( mbh_pose.residue(1).atom_index( at4 ) );
		Size const at5id( mbh_pose.residue(1).atom_index( at5 ) );
		Size const at6id( mbh_pose.residue(1).atom_index( at6 ) );

		Size const oat1id( mbh_pose.residue(1).atom_index( oat1 ) );
		Size const oat2id( mbh_pose.residue(1).atom_index( oat2 ) );
		Size const oat3id( mbh_pose.residue(1).atom_index( oat3 ) );

		RigidLigandBuilderOP rigid_builder( new RigidLigandBuilder );
		rigid_builder->initialize_from_residue(
			at4id, at5id, at6id, oat1id, oat2id, oat3id, mbh_pose.residue(1) );
		rigid_builder->initialize_upstream_residue( cys_restype ); /// do not initialize bonded data.

		BumpGridOP bb_grid = bump_grid_to_enclose_pose( trpcage );
		for ( Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			BumpGridOP resbgop = bump_grid_to_enclose_residue_backbone( trpcage.residue( ii ), *bb_grid );
			fill_grid_with_backbone_heavyatom_spheres( trpcage.residue( ii ), *resbgop );
			bb_grid->or_with( *resbgop );
		}

		ClassicMatchAlgorithmOP match_algorithm( new ClassicMatchAlgorithm( 1 ) );
		match_algorithm->set_residue_type( cys_restype );
		match_algorithm->add_external_geom_sampler(
			exsampler, 1, "CA", "CB", "SG", rigid_builder );
		WriteUpstreamCoordinateKinemageOP dsalgorithm( new WriteUpstreamCoordinateKinemage( fout ) );
		dsalgorithm->set_match_algorithm( match_algorithm );
		dsalgorithm->set_n_downstream_to_output( 15 );
		match_algorithm->set_bb_grid( bb_grid );
		rigid_builder->set_bb_grid( bb_grid );

		SingleDownstreamResidueWriterOP downstream_writer( new SingleDownstreamResidueWriter );
		downstream_writer->set_restype( mbh_pose.residue(1).type().get_self_ptr() );
		downstream_writer->set_downstream_builder( rigid_builder );
		downstream_writer->set_downstream_master( "MBH" );

		dsalgorithm->set_downstream_writer( downstream_writer );

		build_set.set_downstream_algorithm( dsalgorithm );

		ProteinUpstreamBuilder scbuilder;
		scbuilder.add_build_set( build_set );
		scbuilder.set_sampler( ProteinSCSamplerCOP( new DunbrackSCSampler ) );

		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );
		scbuilder.build( *res2bp );

		fout.close();


	}
};

