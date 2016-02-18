// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/match/upstream/ProteinSCSampler.hh>
#include <protocols/match/upstream/ProteinUpstreamBuilder.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>
#include <protocols/match/upstream/UpstreamResTypeGeometry.hh>
#include <protocols/match/output/WriteUpstreamCoordinateKineamge.hh>

#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

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
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>
#include <boost/unordered/unordered_map_fwd.hpp>


using core::Size;
using namespace protocols::match;
using namespace protocols::match::upstream;
using namespace protocols::match::output;

// --------------- Test Class --------------- //

std::string
trpcage_res2phe_kinemage();

utility::vector1< std::string >
trpcage_res2phe_kinemages();


class ProteinUpstreamBuilderTests : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior -override_rsd_type_limit -extra_res_fa protocols/match/PBF.params" ); // kinemages below made with dun02 library
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	void test_sc_upstream_geom_map_atoms_to_controlling_chi() {
		using namespace core;
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		UpstreamResTypeGeometry leu_geom( trpcage.residue_type( 2 ));

		/*
		std::cout << leu_geom.natoms() << std::endl;
		std::cout << leu_geom.nchi() << std::endl;
		for ( Size ii = 1; ii <= leu_geom.natoms(); ++ii ) {
		//std::cout << "Atom: " << ii << " " << trpcage.residue_type( 2 ).atom_name( ii );
		//std::cout << " " << leu_geom.atom_controlled_by_any_chi( ii );
		//std::cout << " " << leu_geom.atom_is_chitip( ii );
		//std::cout << " " << leu_geom.controlling_chi_for_atom()[ ii ];
		//std::cout << " " << leu_geom.which_point_for_atom()[ ii ] << std::endl;
		std::cout << "TS_ASSERT( leu_geom.atom_controlled_by_any_chi( " << ii << " ) = " << leu_geom.atom_controlled_by_any_chi( ii ) << "); " << std::endl;
		std::cout << "TS_ASSERT( leu_geom.atom_is_chitip( " << ii << " ) = " << leu_geom.atom_is_chitip( ii ) << "); " << std::endl;
		std::cout << "TS_ASSERT( leu_geom.controlling_chi_for_atom()[ " << ii << " ] = " << leu_geom.controlling_chi_for_atom()[ ii ] << "); " << std::endl;
		std::cout << "TS_ASSERT( leu_geom.which_point_for_atom()[ " << ii << " ] = " << leu_geom.which_point_for_atom()[ ii ] << "); " << std::endl;

		}

		for ( Size ii = 1; ii <= leu_geom.nchi(); ++ii ) {
		//std::cout << "Chi: " << ii;
		//std::cout << " " << leu_geom.chitip_atoms()[ ii ];
		//std::cout << " " << leu_geom.nonchitip_atoms()[ ii ].size() << std::endl;
		//for ( Size jj = 1; jj <= leu_geom.nonchitip_atoms()[ ii ].size(); ++jj ) {
		// std::cout << "   " << jj << " " << leu_geom.nonchitip_atoms()[ ii ][ jj ] << std::endl;
		//}
		std::cout << "TS_ASSERT( leu_geom.chitip_atoms()[ " << ii << " ] = " << leu_geom.chitip_atoms()[ ii ] << "); " << std::endl;
		std::cout << "TS_ASSERT( leu_geom.nonchitip_atoms()[ " << ii << " ].size() = " << leu_geom.nonchitip_atoms()[ ii ].size() << "); " << std::endl;
		for ( Size jj = 1; jj <= leu_geom.nonchitip_atoms()[ ii ].size(); ++jj ) {
		std::cout << "TS_ASSERT( leu_geom.nonchitip_atoms()[ " << ii << " ][ " << jj << "] = " << leu_geom.nonchitip_atoms()[ ii ][ jj ] << "); " << std::endl;
		}

		}*/

		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 1 ) == 0);
		TS_ASSERT( leu_geom.atom_is_chitip( 1 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 1 ] == 0);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 1 ] == 0);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 2 ) == 0);
		TS_ASSERT( leu_geom.atom_is_chitip( 2 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 2 ] == 0);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 2 ] == 0);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 3 ) == 0);
		TS_ASSERT( leu_geom.atom_is_chitip( 3 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 3 ] == 0);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 3 ] == 0);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 4 ) == 0);
		TS_ASSERT( leu_geom.atom_is_chitip( 4 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 4 ] == 0);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 4 ] == 0);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 5 ) == 0);
		TS_ASSERT( leu_geom.atom_is_chitip( 5 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 5 ] == 0);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 5 ] == 0);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 6 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 6 ) == 1);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 6 ] == 1);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 6 ] == 0);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 7 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 7 ) == 1);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 7 ] == 2);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 7 ] == 0);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 8 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 8 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 8 ] == 2);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 8 ] == 1);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 9 ) == 0);
		TS_ASSERT( leu_geom.atom_is_chitip( 9 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 9 ] == 0);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 9 ] == 0);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 10 ) == 0);
		TS_ASSERT( leu_geom.atom_is_chitip( 10 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 10 ] == 0);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 10 ] == 0);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 11 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 11 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 11 ] == 1);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 11 ] == 1);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 12 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 12 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 12 ] == 1);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 12 ] == 2);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 13 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 13 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 13 ] == 2);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 13 ] == 2);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 14 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 14 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 14 ] == 2);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 14 ] == 3);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 15 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 15 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 15 ] == 2);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 15 ] == 4);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 16 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 16 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 16 ] == 2);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 16 ] == 5);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 17 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 17 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 17 ] == 2);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 17 ] == 6);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 18 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 18 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 18 ] == 2);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 18 ] == 7);
		TS_ASSERT( leu_geom.atom_controlled_by_any_chi( 19 ) == 1);
		TS_ASSERT( leu_geom.atom_is_chitip( 19 ) == 0);
		TS_ASSERT( leu_geom.controlling_chi_for_atom()[ 19 ] == 2);
		TS_ASSERT( leu_geom.which_point_for_atom()[ 19 ] == 8);
		TS_ASSERT( leu_geom.chitip_atoms()[ 1 ] == 6);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 1 ].size() == 2);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 1 ][ 1] == 11);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 1 ][ 2] == 12);
		TS_ASSERT( leu_geom.chitip_atoms()[ 2 ] == 7);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 2 ].size() == 8);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 2 ][ 1] == 8);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 2 ][ 2] == 13);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 2 ][ 3] == 14);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 2 ][ 4] == 15);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 2 ][ 5] == 16);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 2 ][ 6] == 17);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 2 ][ 7] == 18);
		TS_ASSERT( leu_geom.nonchitip_atoms()[ 2 ][ 8] == 19);
	}

	void test_sc_upstream_geom_coords() {
		using namespace core;
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		UpstreamResTypeGeometry leu_geom( trpcage.residue_type( 2 ));

		utility::vector1< UpstreamResTypeGeometry::HTReal > chitip_start_frames, chitip_frames, restype_icoor_chi;
		chitip_start_frames.reserve( leu_geom.nchi() );
		restype_icoor_chi.reserve( leu_geom.nchi() );

		core::chemical::ResidueType const & res( trpcage.residue_type( 2 ) );

		for ( Size ii = 1; ii <= leu_geom.nchi(); ++ii ) {
			//std::cout << "Chi: " << ii;
			//std::cout << " " << leu_geom.chitip_atoms()[ ii ];
			//std::cout << " " << leu_geom.nonchitip_atoms()[ ii ].size() << std::endl;
			//for ( Size jj = 1; jj <= leu_geom.nonchitip_atoms()[ ii ].size(); ++jj ) {
			// std::cout << "   " << jj << " " << leu_geom.nonchitip_atoms()[ ii ][ jj ] << std::endl;
			//}
			UpstreamResTypeGeometry::HTReal ht = leu_geom.ht_for_chitip_atoms()[ ii ];
			//std::cout << "ht:" << ii << std::endl;
			//std::cout << "  " << ht.xx() << " " << ht.yx() << " " << ht.zx() << " " << ht.px() << std::endl;
			//std::cout << "  " << ht.xy() << " " << ht.yy() << " " << ht.zy() << " " << ht.py() << std::endl;
			//std::cout << "  " << ht.xz() << " " << ht.yz() << " " << ht.zz() << " " << ht.pz() << std::endl;

			Size const
				chiat1( res.chi_atoms( ii )[ 1 ] ),
				chiat2( res.chi_atoms( ii )[ 2 ] ),
				chiat3( res.chi_atoms( ii )[ 3 ] ),
				chiat4( res.chi_atoms( ii )[ 4 ] );

			chitip_start_frames.push_back( UpstreamResTypeGeometry::HTReal(
				res.atom( chiat1 ).ideal_xyz(),
				res.atom( chiat2 ).ideal_xyz(),
				res.atom( chiat3 ).ideal_xyz() ));

			chitip_frames.push_back( UpstreamResTypeGeometry::HTReal(
				res.atom( chiat2 ).ideal_xyz(),
				res.atom( chiat3 ).ideal_xyz(),
				res.atom( chiat4 ).ideal_xyz() ));

			//ht = chitip_start_frames[ ii ];
			//std::cout << "chitip_start_frames:" << ii << std::endl;
			//std::cout << "  " << ht.xx() << " " << ht.yx() << " " << ht.zx() << " " << ht.px() << std::endl;
			//std::cout << "  " << ht.xy() << " " << ht.yy() << " " << ht.zy() << " " << ht.py() << std::endl;
			//std::cout << "  " << ht.xz() << " " << ht.yz() << " " << ht.zz() << " " << ht.pz() << std::endl;


			Real const dihedral = numeric::dihedral_degrees(
				res.atom( chiat1 ).ideal_xyz(),
				res.atom( chiat2 ).ideal_xyz(),
				res.atom( chiat3 ).ideal_xyz(),
				res.atom( chiat4 ).ideal_xyz()
			);
			//std::cout << "Ideal dihedral " << ii << " " << dihedral << std::endl;

			//std::cout << "Ideal improper bond angle: "<< ii << " " << res.icoor( chiat4 ).theta() << std::endl;

			UpstreamResTypeGeometry::HTReal chi_transform;
			chi_transform.set_zaxis_rotation_deg( dihedral );
			restype_icoor_chi.push_back( chi_transform );

			//ht = restype_icoor_chi[ ii ];
			//std::cout << "restype_icoor_chi:" << ii << std::endl;
			//std::cout << "  " << ht.xx() << " " << ht.yx() << " " << ht.zx() << " " << ht.px() << std::endl;
			//std::cout << "  " << ht.xy() << " " << ht.yy() << " " << ht.zy() << " " << ht.py() << std::endl;
			//std::cout << "  " << ht.xz() << " " << ht.yz() << " " << ht.zz() << " " << ht.pz() << std::endl;

		}

		for ( Size ii = 1; ii <= leu_geom.natoms(); ++ii ) {
			//std::cout << "Atom " << ii << " Icoor: ";
			//std::cout << trpcage.residue_type( 2 ).xyz( ii ).x() << " " <<trpcage.residue_type( 2 ).xyz( ii ).y() << " " <<trpcage.residue_type( 2 ).xyz( ii ).z() << std::endl;
			Vector ideal = trpcage.residue_type( 2 ).atom( ii ).ideal_xyz();

			Size controlling_chi = leu_geom.controlling_chi_for_atom()[ ii ];
			if ( controlling_chi == 0 ) continue;
			Size which_vector    = leu_geom.which_point_for_atom()[ ii ];
			if ( which_vector == 0 ) continue;

			UpstreamResTypeGeometry::HTReal ht = leu_geom.ht_for_chitip_atoms()[ controlling_chi ];
			UpstreamResTypeGeometry::Vector v  = leu_geom.points_for_nonchitip_atoms()[ controlling_chi ][ which_vector ];
			//std::cout << "Atom " << ii << " v: ";
			//std::cout << v.x() << " " <<v.y() << " " <<v.z() << std::endl;

			//Vector launch_point = ( chitip_start_frames[ controlling_chi ] * restype_icoor_chi[ controlling_chi ] * ht ).point();
			//std::cout << "Atom " << ii << " launch_point: ";
			//std::cout << launch_point.x() << " " <<launch_point.y() << " " <<launch_point.z() << std::endl;

			Vector recreated = chitip_start_frames[ controlling_chi ] * restype_icoor_chi[ controlling_chi ] * ht * v;

			//std::cout << "Atom " << ii << " recreated: ";
			//std::cout << recreated.x() << " " <<recreated.y() << " " <<recreated.z() << std::endl;


			Vector recreated2 = chitip_frames[ controlling_chi ] * v;
			//std::cout << "Atom " << ii << " recreated2: ";
			//std::cout << recreated2.x() << " " <<recreated2.y() << " " <<recreated2.z() << std::endl;
			TS_ASSERT( recreated.distance_squared( ideal ) < 1e-6 );
			TS_ASSERT( recreated.distance_squared( recreated2 ) < 1e-6 );
		}
	}

	void test_build_set_ctor() {
		using namespace core::chemical;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		ResidueTypeSetCOP res2_set( trpcage.residue( 2 ).residue_type_set() );
		//ResidueTypeCOPs const & aas( ResidueTypeFinder( *res2_set ).aa( aa_phe ).get_all_possible_residue_types() );

		BuildSet build_set;
		build_set.set_residue_type( ResidueTypeFinder( *res2_set ).aa( aa_phe ).variants( trpcage.residue( 2 ).type().variant_types() ).get_representative_type() );

		/// Find the matching phe residue type for residue 2.
		/*for ( ResidueTypeCOPs::const_iterator
		aas_iter = aas.begin(),
		aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( variants_match( trpcage.residue( 2 ).type(), **aas_iter ) ) {
		build_set.set_residue_type( *aas_iter );
		break;
		}
		}*/

		TS_ASSERT( build_set.restype_geometry().N_atom_id()  == 1 );
		TS_ASSERT( build_set.restype_geometry().CA_atom_id() == 2 );
		TS_ASSERT( build_set.restype_geometry().C_atom_id()  == 3 );
		TS_ASSERT( build_set.restype_geometry().O_atom_id()  == 4 );
		TS_ASSERT( build_set.restype_geometry().CB_atom_id() == 5 );

		TS_ASSERT( build_set.restype_geometry().nchi() == 2 );
		TS_ASSERT( build_set.restype_geometry().natoms() == 20 );
	}

	void test_full_sample_set_ctor()
	{
		using namespace core;
		using namespace core::chemical;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		ResidueTypeSetCOP res2_set( trpcage.residue( 2 ).residue_type_set() );
		//ResidueTypeCOPs const & aas( ResidueTypeFinder( *res2_set ).aa( aa_phe ).get_all_possible_residue_types() );

		BuildSet build_set;
		build_set.set_residue_type( ResidueTypeFinder( *res2_set ).aa( aa_phe ).variants( trpcage.residue( 2 ).type().variant_types() ).get_representative_type() );

		/// Find the matching phe residue type for residue 2.
		/*for ( ResidueTypeCOPs::const_iterator
		aas_iter = aas.begin(),
		aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( variants_match( trpcage.residue( 2 ).type(), **aas_iter ) ) {
		build_set.set_residue_type( *aas_iter );
		break;
		}
		}*/

		SampleStrategyData strat; strat.set_strategy( rotameric_chi_mimic_EX_flags );
		strat.set_sample_level( core::pack::task::EX_ONE_STDDEV );

		build_set.set_sample_strategy_for_chi( 1, strat );
		build_set.set_sample_strategy_for_chi( 2, strat );

		DunbrackSCSampler dunsampler;
		DunbrackSCSampler::DunbrackRotamerSampleDataVector dun_samps(
			dunsampler.samples( *res2bp, build_set.restype() ));
		FullChiSampleSet expanded_samples( build_set, dun_samps[ 1 ], false );


		/*
		for ( Size ii = 1; ii <= 2; ++ii ) {
		//std::cout << "Chi " << ii << " with " << expanded_samples.n_samples_per_chi()[ ii ] << ":" << std::endl;
		std::cout << "TS_ASSERT( expanded_samples.n_samples_per_chi()[ " << ii << " ] = " << expanded_samples.n_samples_per_chi()[ ii ] << " );" << std::endl;
		for ( Size jj = 1; jj <= expanded_samples.n_samples_per_chi()[ ii ]; ++jj ) {
		//std::cout << jj << ": " << expanded_samples.chi_sample( ii, jj ) << std::endl;
		std::cout << "TS_ASSERT_DELTA( expanded_samples.chi_sample( " << ii << ", " << jj << " ), " << expanded_samples.chi_sample( ii, jj ) << ", 1e-3 );" << std::endl;
		}
		}*/

		TS_ASSERT( expanded_samples.n_samples_per_chi()[ 1 ] = 3 );
		TS_ASSERT_DELTA( expanded_samples.chi_sample( 1, 1 ), 168.735, 3e-3 );
		TS_ASSERT_DELTA( expanded_samples.chi_sample( 1, 2 ), 178.019, 3e-3 );
		TS_ASSERT_DELTA( expanded_samples.chi_sample( 1, 3 ), 187.303, 3e-3 );

		TS_ASSERT( expanded_samples.n_samples_per_chi()[ 2 ] = 3 );
		TS_ASSERT_DELTA( expanded_samples.chi_sample( 2, 1 ), 68.3786, 3e-3 );
		TS_ASSERT_DELTA( expanded_samples.chi_sample( 2, 2 ), 80.3798, 3e-3 );
		TS_ASSERT_DELTA( expanded_samples.chi_sample( 2, 3 ), 92.381,  3e-3 );
	}

	void test_ProteinUpstreamBuilder_build()
	{
		using namespace core;
		using namespace core::chemical;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		ResidueTypeSetCOP res2_set( trpcage.residue( 2 ).residue_type_set() );
		//ResidueTypeCOPs const & aas( ResidueTypeFinder( *res2_set ).aa( aa_phe ).get_all_possible_residue_types() );

		BuildSet build_set;
		build_set.set_residue_type( ResidueTypeFinder( *res2_set ).aa( aa_phe ).variants( trpcage.residue( 2 ).type().variant_types() ).get_representative_type() );

		/// Find the matching phe residue type for residue 2.
		/*for ( ResidueTypeCOPs::const_iterator
		aas_iter = aas.begin(),
		aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( variants_match( trpcage.residue( 2 ).type(), **aas_iter ) ) {
		build_set.set_residue_type( *aas_iter );
		break;
		}
		}*/
		SampleStrategyData strat; strat.set_strategy( rotameric_chi_mimic_EX_flags );
		strat.set_sample_level( core::pack::task::EX_ONE_STDDEV );

		build_set.set_sample_strategy_for_chi( 1, strat );
		build_set.set_sample_strategy_for_chi( 2, strat );

		std::ostringstream sout;

		WriteUpstreamCoordinateKinemageOP dsalgorithm( new WriteUpstreamCoordinateKinemage( sout ) );

		build_set.set_downstream_algorithm( dsalgorithm );

		ProteinUpstreamBuilder scbuilder;
		scbuilder.add_build_set( build_set );
		scbuilder.set_sampler( ProteinSCSamplerCOP( new DunbrackSCSampler ) );

		BumpGridOP bbgrid( new BumpGrid );
		scbuilder.set_bb_grid( bbgrid );

		scbuilder.build( *res2bp );

		std::string correct_kinemage = trpcage_res2phe_kinemage();
		TS_ASSERT( sout.str() == correct_kinemage );
		if ( sout.str() != correct_kinemage ) {
			std::cout << sout.str() << std::endl;
		}

	}

	void test_ProteinUpstreamBuilder_build_w_pre_chitip_transform()
	{
		// use a residue type (a non canonical amino acid) that has a rigid segment
		// between two of the chis, so that the pre_chitip_transform is needed to
		// correctly create the HT before the chi-tip atom can be placed.
		using namespace core;
		using namespace core::chemical;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		ResidueTypeSetCOP res2_set( trpcage.residue( 2 ).residue_type_set() );
		ResidueTypeCOP restype( res2_set->name_mapOP( "PBF" ));

		BuildSet build_set;
		build_set.set_residue_type( restype );

		//// Find the matching phe residue type for residue 2.
		//for ( ResidueTypeCOPs::const_iterator
		//  aas_iter = aas.begin(),
		//  aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		// if ( variants_match( trpcage.residue( 2 ).type(), **aas_iter ) ) {
		//  build_set.set_residue_type( *aas_iter );
		//  break;
		// }
		//}
		SampleStrategyData strat;
		//strat.set_strategy( rotameric_chi_mimic_EX_flags );
		//strat.set_sample_level( core::pack::task::EX_ONE_STDDEV );

		build_set.set_sample_strategy_for_chi( 1, strat );
		build_set.set_sample_strategy_for_chi( 2, strat );
		build_set.set_sample_strategy_for_chi( 3, strat );
		build_set.set_sample_strategy_for_chi( 4, strat );

		std::ostringstream sout;

		WriteUpstreamCoordinateKinemageOP dsalgorithm( new WriteUpstreamCoordinateKinemage( sout ) );

		build_set.set_downstream_algorithm( dsalgorithm );

		ProteinUpstreamBuilder scbuilder;
		scbuilder.add_build_set( build_set );
		scbuilder.set_sampler( ProteinSCSamplerCOP( new DunbrackSCSampler ) );

		BumpGridOP bbgrid( new BumpGrid );
		scbuilder.set_bb_grid( bbgrid );

		scbuilder.build( *res2bp );

		std::ifstream ifs( "protocols/match/pbf_on_trpcage_res2_gold.kin" );
		std::string gold_kin;
		utility::slurp( ifs, gold_kin );

		TS_ASSERT_EQUALS( gold_kin, sout.str() );


		//std::cout << "PBF: " <<  sout.str() << std::endl;

	}

	void test_ProteinUpstreamBuilder_recover_hit()
	{
		using namespace core;
		using namespace core::chemical;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		ResidueTypeSetCOP res2_set( trpcage.residue( 2 ).residue_type_set() );
		//ResidueTypeCOPs const & aas( ResidueTypeFinder( *res2_set ).aa( aa_phe ).get_all_possible_residue_types() );

		BuildSet build_set;
		build_set.set_residue_type( ResidueTypeFinder( *res2_set ).aa( aa_phe ).variants( trpcage.residue( 2 ).type().variant_types() ).get_representative_type() );

		/// Find the matching phe residue type for residue 2.
		/*for ( ResidueTypeCOPs::const_iterator
		aas_iter = aas.begin(),
		aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( variants_match( trpcage.residue( 2 ).type(), **aas_iter ) ) {
		build_set.set_residue_type( *aas_iter );
		break;
		}
		}*/
		SampleStrategyData strat; strat.set_strategy( rotameric_chi_mimic_EX_flags );
		strat.set_sample_level( core::pack::task::EX_ONE_STDDEV );

		build_set.set_sample_strategy_for_chi( 1, strat );
		build_set.set_sample_strategy_for_chi( 2, strat );

		std::ostringstream sout, sout2;

		WriteUpstreamCoordinateKinemageOP dsalgorithm( new WriteUpstreamCoordinateKinemage( sout ) );
		dsalgorithm->return_pseudo_hits( true );

		WriteUpstreamHitKinemageOP kin_processor( new WriteUpstreamHitKinemage( sout2 ) );
		kin_processor->set_master( "rotamers" );
		kin_processor->animate( true );
		kin_processor->dominant( true );
		kin_processor->group( true );

		build_set.set_downstream_algorithm( dsalgorithm );

		ProteinUpstreamBuilder scbuilder;
		scbuilder.add_build_set( build_set );
		scbuilder.set_sampler( ProteinSCSamplerCOP( new DunbrackSCSampler ) );

		BumpGridOP bbgrid( new BumpGrid );
		scbuilder.set_bb_grid( bbgrid );

		std::list< Hit > hitlist = scbuilder.build( *res2bp );
		Size counter( 0 );
		//std::cout << "hits: " << std::endl;

		utility::vector1< std::string > correct_rotamer_kins = trpcage_res2phe_kinemages();

		if ( hitlist.size() != 36 ) {
			TS_ASSERT( false );
			return;
		}

		for ( std::list< Hit >::const_iterator iter = hitlist.begin(),
				iter_end = hitlist.end();
				iter != iter_end; ++iter ) {
			/*std::cout << "Hit " << ++counter << " ";
			for ( Size ii = 1; ii <= 4; ++ii ) {
			std::cout << iter->first()[ ii ] << " ";
			}
			for ( Size ii = 1; ii <= 6; ++ii ) {
			std::cout << iter->second()[ ii ]<< " ";
			}
			std::cout << std::endl;*/

			scbuilder.recover_hit( *iter, *res2bp, *kin_processor );

			TS_ASSERT( correct_rotamer_kins[ ++counter ] == sout2.str() );
			if ( correct_rotamer_kins[ counter ] != sout2.str() ) {
				std::cout << "Correct: " << std::endl;
				std::cout << correct_rotamer_kins[ counter ] << std::endl << "Incorrect:" << std::endl;
				std::cout << sout2.str() << std::endl;
			}
			//sout2.clear();
			sout2.str("");
		}

	}


	void test_ProteinUpstreamBuilder_recover_hits()
	{
		using namespace core;
		using namespace core::chemical;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		ResidueTypeSetCOP res2_set( trpcage.residue( 2 ).residue_type_set() );
		//ResidueTypeCOPs const & aas( ResidueTypeFinder( *res2_set ).aa( aa_phe ).get_all_possible_residue_types() );

		BuildSet build_set;
		build_set.set_residue_type( ResidueTypeFinder( *res2_set ).aa( aa_phe ).variants( trpcage.residue( 2 ).type().variant_types() ).get_representative_type() );

		/// Find the matching phe residue type for residue 2.
		/*for ( ResidueTypeCOPs::const_iterator
		aas_iter = aas.begin(),
		aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( variants_match( trpcage.residue( 2 ).type(), **aas_iter ) ) {
		build_set.set_residue_type( *aas_iter );
		break;
		}
		}*/
		SampleStrategyData strat; strat.set_strategy( rotameric_chi_mimic_EX_flags );
		strat.set_sample_level( core::pack::task::EX_ONE_STDDEV );

		build_set.set_sample_strategy_for_chi( 1, strat );
		build_set.set_sample_strategy_for_chi( 2, strat );

		std::ostringstream sout, sout2;

		WriteUpstreamCoordinateKinemageOP dsalgorithm( new WriteUpstreamCoordinateKinemage( sout ) );
		dsalgorithm->return_pseudo_hits( true );

		WriteUpstreamHitKinemageOP kin_processor( new WriteUpstreamHitKinemage( sout2 ) );
		kin_processor->set_master( "rotamers" );
		kin_processor->animate( true );
		kin_processor->dominant( true );
		kin_processor->group( true );

		build_set.set_downstream_algorithm( dsalgorithm );

		ProteinUpstreamBuilder scbuilder;
		scbuilder.add_build_set( build_set );
		scbuilder.set_sampler( ProteinSCSamplerCOP( new DunbrackSCSampler ) );

		BumpGridOP bbgrid( new BumpGrid );
		scbuilder.set_bb_grid( bbgrid );

		std::list< Hit > hitlist = scbuilder.build( *res2bp );
		//std::cout << "hits: " << std::endl;

		utility::vector1< std::string > correct_rotamer_kins = trpcage_res2phe_kinemages();

		if ( hitlist.size() != 36 ) {
			TS_ASSERT( false );
			return;
		}


		std::list< Hit > hit_subset;
		std::string hit_string;
		Size counter( 0 );
		for ( std::list< Hit >::const_iterator iter = hitlist.begin(),
				iter_end = hitlist.end();
				iter != iter_end; ++iter ) {
			hit_subset.push_back( *iter );
			hit_string += correct_rotamer_kins[ ++counter ];
			++iter;
			if ( iter == iter_end ) break;
			++counter;
		}

		scbuilder.recover_hits( hit_subset.begin(), hit_subset.end(), *res2bp, *kin_processor );

		TS_ASSERT( hit_string == sout2.str() );

	}

	void test_FullChiSampleSet_exflags_dryrun_vs_regular()
	{
		using namespace core;
		using namespace core::chemical;


		typedef utility::pointer::shared_ptr< FullChiSampleSet > FullChiSampleSetOP;

		FullChiSampleSetOP sampset_dry_run, sampset_regular;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		ResidueTypeSetCOP res2_set( trpcage.residue( 2 ).residue_type_set() );
		//ResidueTypeCOPs const & aas( ResidueTypeFinder( *res2_set ).aa( aa_phe ).get_all_possible_residue_types() );

		BuildSet build_set;
		build_set.set_residue_type( ResidueTypeFinder( *res2_set ).aa( aa_phe ).variants( trpcage.residue( 2 ).type().variant_types() ).get_representative_type() );

		/// Find the matching phe residue type for residue 2.
		/*for ( ResidueTypeCOPs::const_iterator
		aas_iter = aas.begin(),
		aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( variants_match( trpcage.residue( 2 ).type(), **aas_iter ) ) {
		build_set.set_residue_type( *aas_iter );
		break;
		}
		}*/

		ProteinSCSamplerOP sampler_( new DunbrackSCSampler );
		ProteinSCSampler::DunbrackRotamerSampleDataVector
			rotamer_samples = sampler_->samples( *res2bp, build_set.restype() );

		for ( Size ii = 1; ii <= rotamer_samples.size(); ++ii ) {

			SampleStrategyData strat; strat.set_strategy( rotameric_chi_mimic_EX_flags );

			for ( Size jj = 0; jj < core::pack::task::ExtraRotSampleCardinality; ++jj ) {
				strat.set_sample_level( core::pack::task::ExtraRotSample( jj ) );
				build_set.set_sample_strategy_for_chi( 1, strat );
				build_set.set_sample_strategy_for_chi( 2, strat );

				sampset_dry_run = FullChiSampleSetOP( new FullChiSampleSet( build_set, rotamer_samples[ ii ], true ) );
				sampset_regular = FullChiSampleSetOP( new FullChiSampleSet( build_set, rotamer_samples[ ii ], false ) );

				TS_ASSERT( sampset_dry_run->num_chi_samples_total() == sampset_regular->num_chi_samples_total() );
			}

		}

	}

	void test_FullChiSampleSet_step_wi_sdrange_dryrun_vs_regular()
	{
		using namespace core;
		using namespace core::chemical;


		typedef utility::pointer::shared_ptr< FullChiSampleSet > FullChiSampleSetOP;

		FullChiSampleSetOP sampset_dry_run, sampset_regular;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		ResidueTypeSetCOP res2_set( trpcage.residue( 2 ).residue_type_set() );
		//ResidueTypeCOPs const & aas( ResidueTypeFinder( *res2_set ).aa( aa_phe ).get_all_possible_residue_types() );

		BuildSet build_set;
		build_set.set_residue_type( ResidueTypeFinder( *res2_set ).aa( aa_phe ).variants( trpcage.residue( 2 ).type().variant_types() ).get_representative_type() );

		/// Find the matching phe residue type for residue 2.
		/*for ( ResidueTypeCOPs::const_iterator
		aas_iter = aas.begin(),
		aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( variants_match( trpcage.residue( 2 ).type(), **aas_iter ) ) {
		build_set.set_residue_type( *aas_iter );
		break;
		}
		}*/

		ProteinSCSamplerOP sampler_( new DunbrackSCSampler );
		ProteinSCSampler::DunbrackRotamerSampleDataVector
			rotamer_samples = sampler_->samples( *res2bp, build_set.restype() );

		SampleStrategyData default_strat; default_strat.set_strategy( rotameric_chi_mimic_EX_flags );
		default_strat.set_sample_level( core::pack::task::EX_ONE_STDDEV );

		for ( Size ii = 1; ii <= rotamer_samples.size(); ++ii ) {

			SampleStrategyData strat1; strat1.set_strategy( rotameric_chi_step_wi_sd_range );
			strat1.set_sd_range( 2.0 );
			strat1.set_step_size( 3.5 );

			build_set.set_sample_strategy_for_chi( 1, strat1 );
			build_set.set_sample_strategy_for_chi( 2, default_strat );

			sampset_dry_run = FullChiSampleSetOP( new FullChiSampleSet( build_set, rotamer_samples[ ii ], true ) );
			sampset_regular = FullChiSampleSetOP( new FullChiSampleSet( build_set, rotamer_samples[ ii ], false ) );

			TS_ASSERT( sampset_dry_run->num_chi_samples_total() == sampset_regular->num_chi_samples_total() );

		}

	}


	void test_FullChiSampleSet_semirotameric_dryrun_vs_regular()
	{
		using namespace core;
		using namespace core::chemical;


		typedef utility::pointer::shared_ptr< FullChiSampleSet > FullChiSampleSetOP;
		typedef core::pack::dunbrack::DunbrackRotamerSampleData DunbrackRotamerSampleData;

		FullChiSampleSetOP sampset_dry_run, sampset_regular;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		OriginalBackboneBuildPointOP res2bp( new OriginalBackboneBuildPoint( trpcage.residue( 2 ), 1 ) );

		ResidueTypeSetCOP res2_set( trpcage.residue( 2 ).residue_type_set() );
		//ResidueTypeCOPs const & aas( ResidueTypeFinder( *res2_set ).aa( aa_phe ).get_all_possible_residue_types() );

		BuildSet build_set;
		build_set.set_residue_type( ResidueTypeFinder( *res2_set ).aa( aa_phe ).variants( trpcage.residue( 2 ).type().variant_types() ).get_representative_type() );
		/// Find the matching phe residue type for residue 2.
		/*
		for ( ResidueTypeCOPs::const_iterator
		aas_iter = aas.begin(),
		aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( variants_match( trpcage.residue( 2 ).type(), **aas_iter ) ) {
		build_set.set_residue_type( *aas_iter );
		break;
		}
		}*/


		//ProteinSCSamplerOP sampler_ = new DunbrackSCSampler;
		//ProteinSCSampler::DunbrackRotamerSampleDataVector
		// rotamer_samples = sampler_->samples( *res2bp, build_set.restype() );


		SampleStrategyData strat1, strat2;
		strat1.set_strategy( rotameric_chi_mimic_EX_flags );

		DunbrackRotamerSampleData nrchi_sample;

		/// This sample isn't meant to reflect any phi/psi in particular; it's just for testing
		nrchi_sample.set_nrchi_sample( true );
		nrchi_sample.set_nchi( 2 );
		nrchi_sample.set_rotwell(  1, 2 );
		nrchi_sample.set_rotwell(  2, 1 );
		nrchi_sample.set_chi_mean( 1, 178 );
		nrchi_sample.set_chi_mean( 2, 95 );
		nrchi_sample.set_chi_sd(   1, 8.5 );
		nrchi_sample.set_chi_sd(   2, 20.5 );
		nrchi_sample.set_prob( 0.4 );

		nrchi_sample.set_nrchi_lower_boundary( 75 );
		nrchi_sample.set_nrchi_upper_boundary( 105 );
		nrchi_sample.set_nrchi_probability( 0.3 );

		// 1. Regular ex flag behavior
		strat2.set_strategy( rotameric_chi_mimic_EX_flags );
		strat2.set_sample_level( core::pack::task::EX_ONE_STDDEV );
		build_set.set_sample_strategy_for_chi( 2, strat2 );

		for ( Size ii = 0; ii < core::pack::task::ExtraRotSampleCardinality; ++ii ) {
			strat1.set_sample_level( core::pack::task::ExtraRotSample( ii ) );
			build_set.set_sample_strategy_for_chi( 1, strat1 );

			sampset_dry_run = FullChiSampleSetOP( new FullChiSampleSet( build_set, nrchi_sample, true ) );
			sampset_regular = FullChiSampleSetOP( new FullChiSampleSet( build_set, nrchi_sample, false ) );

			TS_ASSERT( sampset_dry_run->num_chi_samples_total() == sampset_regular->num_chi_samples_total() );
		}

		// 2. Expanded samples in the nrchi bin
		strat2.set_strategy( nonrotameric_chi_sample_wi_nrchi_bin );
		strat2.set_n_samples_per_side_of_nrchi_bin( 2 );
		build_set.set_sample_strategy_for_chi( 2, strat2 );

		for ( Size ii = 0; ii < core::pack::task::ExtraRotSampleCardinality; ++ii ) {
			strat1.set_sample_level( core::pack::task::ExtraRotSample( ii ) );
			build_set.set_sample_strategy_for_chi( 1, strat1 );

			sampset_dry_run = FullChiSampleSetOP( new FullChiSampleSet( build_set, nrchi_sample, true ) );
			sampset_regular = FullChiSampleSetOP( new FullChiSampleSet( build_set, nrchi_sample, false ) );

			TS_ASSERT( sampset_dry_run->num_chi_samples_total() == sampset_regular->num_chi_samples_total() );
		}

	}


};

std::string trpcage_res2phe_kinemage() {
	utility::vector1< std::string > kins = trpcage_res2phe_kinemages();
	std::string concatenated( "@kinemage {1}\n"
		"@title { matcher }\n"
		"@1center -4.3903 5.36323 -2.52079\n"
		"@1span 25\n" );
	for ( Size ii = 1; ii <= kins.size(); ++ii ) {
		concatenated += kins[ ii ];
	}
	return concatenated;
}


utility::vector1< std::string >
trpcage_res2phe_kinemages()
{
	utility::vector1< std::string > kins;
	kins.reserve( 36 );
	kins.push_back(
		"@group { rot1 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.88935 5.37612 -2.58147\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.74274 5.87776 -3.41379\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.67091 5.93644 -1.63817\n"
		"{PHE 2  CG  (aroC)}P -2.88935 5.37612 -2.58147\n"
		"{PHE 2  CD1 (aroC)} -2.22223 4.90542 -3.70266\n"
		"{PHE 2  CG  (aroC)}P -2.88935 5.37612 -2.58147\n"
		"{PHE 2  CD2 (aroC)} -2.1415 5.85796 -1.51789\n"
		"{PHE 2  CD1 (aroC)}P -2.22223 4.90542 -3.70266\n"
		"{PHE 2  CE1 (aroC)} -0.841269 4.91666 -3.75952\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.22223 4.90542 -3.70266\n"
		"{PHE 2  HD1 (Haro)} 'h' -2.79998 4.52436 -4.54516\n"
		"{PHE 2  CD2 (aroC)}P -2.1415 5.85796 -1.51789\n"
		"{PHE 2  CE2 (aroC)} -0.761363 5.87156 -1.57235\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.1415 5.85796 -1.51789\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.65554 6.23023 -0.630985\n"
		"{PHE 2  CE1 (aroC)}P -0.841269 4.91666 -3.75952\n"
		"{PHE 2  CZ  (aroC)} -0.110832 5.39981 -2.69485\n"
		"{PHE 2  CE1 (aroC)}P 'h' -0.841269 4.91666 -3.75952\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.329952 4.5436 -4.64657\n"
		"{PHE 2  CE2 (aroC)}P -0.761363 5.87156 -1.57235\n"
		"{PHE 2  CZ  (aroC)} -0.110832 5.39981 -2.69485\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.761363 5.87156 -1.57235\n"
		"{PHE 2  HE2 (Haro)} 'h' -0.185372 6.25357 -0.729709\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.110832 5.39981 -2.69485\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.97735 5.40812 -2.73835\n" );
	kins.push_back(
		"@group { rot2 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.88935 5.37612 -2.58147\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.74274 5.87776 -3.41379\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.67091 5.93644 -1.63817\n"
		"{PHE 2  CG  (aroC)}P -2.88935 5.37612 -2.58147\n"
		"{PHE 2  CD1 (aroC)} -2.22731 5.14317 -3.77776\n"
		"{PHE 2  CG  (aroC)}P -2.88935 5.37612 -2.58147\n"
		"{PHE 2  CD2 (aroC)} -2.13643 5.62032 -1.44299\n"
		"{PHE 2  CD1 (aroC)}P -2.22731 5.14317 -3.77776\n"
		"{PHE 2  CE1 (aroC)} -0.846352 5.15464 -3.83472\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.22731 5.14317 -3.77776\n"
		"{PHE 2  HD1 (Haro)} 'h' -2.80907 4.95022 -4.6795\n"
		"{PHE 2  CD2 (aroC)}P -2.13643 5.62032 -1.44299\n"
		"{PHE 2  CE2 (aroC)} -0.756278 5.6336 -1.49711\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.13643 5.62032 -1.44299\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.64645 5.80441 -0.496731\n"
		"{PHE 2  CE1 (aroC)}P -0.846352 5.15464 -3.83472\n"
		"{PHE 2  CZ  (aroC)} -0.110835 5.40004 -2.69486\n"
		"{PHE 2  CE1 (aroC)}P 'h' -0.846352 5.15464 -3.83472\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.339058 4.96983 -4.78128\n"
		"{PHE 2  CE2 (aroC)}P -0.756278 5.6336 -1.49711\n"
		"{PHE 2  CZ  (aroC)} -0.110835 5.40004 -2.69486\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.756278 5.6336 -1.49711\n"
		"{PHE 2  HE2 (Haro)} 'h' -0.176269 5.82746 -0.595038\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.110835 5.40004 -2.69486\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.977339 5.40826 -2.73859\n" );
	kins.push_back(
		"@group { rot3 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.88935 5.37612 -2.58147\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.74274 5.87776 -3.41379\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.67091 5.93644 -1.63817\n"
		"{PHE 2  CG  (aroC)}P -2.88935 5.37612 -2.58147\n"
		"{PHE 2  CD1 (aroC)} -2.23041 5.39137 -3.80182\n"
		"{PHE 2  CG  (aroC)}P -2.88935 5.37612 -2.58147\n"
		"{PHE 2  CD2 (aroC)} -2.13334 5.37227 -1.41911\n"
		"{PHE 2  CD1 (aroC)}P -2.23041 5.39137 -3.80182\n"
		"{PHE 2  CE1 (aroC)} -0.84946 5.40308 -3.85883\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.23041 5.39137 -3.80182\n"
		"{PHE 2  HD1 (Haro)} 'h' -2.81462 5.39476 -4.72241\n"
		"{PHE 2  CD2 (aroC)}P -2.13334 5.37227 -1.41911\n"
		"{PHE 2  CE2 (aroC)} -0.753169 5.38517 -1.47297\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.13334 5.37227 -1.41911\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.6409 5.35992 -0.453895\n"
		"{PHE 2  CE1 (aroC)}P -0.84946 5.40308 -3.85883\n"
		"{PHE 2  CZ  (aroC)} -0.110835 5.40025 -2.69483\n"
		"{PHE 2  CE1 (aroC)}P 'h' -0.84946 5.40308 -3.85883\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.344625 5.4148 -4.82449\n"
		"{PHE 2  CE2 (aroC)}P -0.753169 5.38517 -1.47297\n"
		"{PHE 2  CZ  (aroC)} -0.110835 5.40025 -2.69483\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.753169 5.38517 -1.47297\n"
		"{PHE 2  HE2 (Haro)} 'h' -0.170703 5.38261 -0.551853\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.110835 5.40025 -2.69483\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.977329 5.40845 -2.73879\n" );
	kins.push_back(
		"@group { rot4 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.91477 5.37573 -2.8025\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.86453 5.93496 -3.31743\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.53026 5.87946 -1.57198\n"
		"{PHE 2  CG  (aroC)}P -2.91477 5.37573 -2.8025\n"
		"{PHE 2  CD1 (aroC)} -2.42659 4.97601 -4.03764\n"
		"{PHE 2  CG  (aroC)}P -2.91477 5.37573 -2.8025\n"
		"{PHE 2  CD2 (aroC)} -2.01202 5.78626 -1.83339\n"
		"{PHE 2  CD1 (aroC)}P -2.42659 4.97601 -4.03764\n"
		"{PHE 2  CE1 (aroC)} -1.06918 4.98696 -4.29788\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.42659 4.97601 -4.03764\n"
		"{PHE 2  HD1 (Haro)} 'h' -3.1272 4.6514 -4.80744\n"
		"{PHE 2  CD2 (aroC)}P -2.01202 5.78626 -1.83339\n"
		"{PHE 2  CE2 (aroC)} -0.655041 5.79941 -2.09099\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.01202 5.78626 -1.83339\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.38466 6.10234 -0.858383\n"
		"{PHE 2  CE1 (aroC)}P -1.06918 4.98696 -4.29788\n"
		"{PHE 2  CZ  (aroC)} -0.183472 5.39877 -3.32504\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.06918 4.98696 -4.29788\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.699255 4.6701 -5.27268\n"
		"{PHE 2  CE2 (aroC)}P -0.655041 5.79941 -2.09099\n"
		"{PHE 2  CZ  (aroC)} -0.183472 5.39877 -3.32504\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.655041 5.79941 -2.09099\n"
		"{PHE 2  HE2 (Haro)} 'h' 0.0438573 6.12496 -1.32074\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.183472 5.39877 -3.32504\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.88634 5.40677 -3.52886\n" );
	kins.push_back(
		"@group { rot5 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.91477 5.37573 -2.8025\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.86453 5.93496 -3.31743\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.53026 5.87946 -1.57198\n"
		"{PHE 2  CG  (aroC)}P -2.91477 5.37573 -2.8025\n"
		"{PHE 2  CD1 (aroC)} -2.43979 5.2181 -4.09603\n"
		"{PHE 2  CG  (aroC)}P -2.91477 5.37573 -2.8025\n"
		"{PHE 2  CD2 (aroC)} -1.99886 5.54429 -1.77519\n"
		"{PHE 2  CD1 (aroC)}P -2.43979 5.2181 -4.09603\n"
		"{PHE 2  CE1 (aroC)} -1.08239 5.22928 -4.35636\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.43979 5.2181 -4.09603\n"
		"{PHE 2  HD1 (Haro)} 'h' -3.1508 5.08502 -4.91184\n"
		"{PHE 2  CD2 (aroC)}P -1.99886 5.54429 -1.77519\n"
		"{PHE 2  CE2 (aroC)} -0.641819 5.5571 -2.03249\n"
		"{PHE 2  CD2 (aroC)}P 'h' -1.99886 5.54429 -1.77519\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.36107 5.66876 -0.754057\n"
		"{PHE 2  CE1 (aroC)}P -1.08239 5.22928 -4.35636\n"
		"{PHE 2  CZ  (aroC)} -0.183473 5.39899 -3.32504\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.08239 5.22928 -4.35636\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.722931 5.10411 -5.37744\n"
		"{PHE 2  CE2 (aroC)}P -0.641819 5.5571 -2.03249\n"
		"{PHE 2  CZ  (aroC)} -0.183473 5.39899 -3.32504\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.641819 5.5571 -2.03249\n"
		"{PHE 2  HE2 (Haro)} 'h' 0.0675251 5.69107 -1.21603\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.183473 5.39899 -3.32504\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.886296 5.40692 -3.52908\n" );
	kins.push_back(
		"@group { rot6 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.91477 5.37573 -2.8025\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.86453 5.93496 -3.31743\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.53026 5.87946 -1.57198\n"
		"{PHE 2  CG  (aroC)}P -2.91477 5.37573 -2.8025\n"
		"{PHE 2  CD1 (aroC)} -2.44336 5.46734 -4.10367\n"
		"{PHE 2  CG  (aroC)}P -2.91477 5.37573 -2.8025\n"
		"{PHE 2  CD2 (aroC)} -1.99532 5.29521 -1.76771\n"
		"{PHE 2  CD1 (aroC)}P -2.44336 5.46734 -4.10367\n"
		"{PHE 2  CE1 (aroC)} -1.08597 5.47876 -4.36404\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.44336 5.46734 -4.10367\n"
		"{PHE 2  HD1 (Haro)} 'h' -3.15717 5.53141 -4.92537\n"
		"{PHE 2  CD2 (aroC)}P -1.99532 5.29521 -1.76771\n"
		"{PHE 2  CE2 (aroC)} -0.638233 5.30763 -2.02477\n"
		"{PHE 2  CD2 (aroC)}P 'h' -1.99532 5.29521 -1.76771\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.35472 5.22243 -0.740604\n"
		"{PHE 2  CE1 (aroC)}P -1.08597 5.47876 -4.36404\n"
		"{PHE 2  CZ  (aroC)} -0.183465 5.39921 -3.32499\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.08597 5.47876 -4.36404\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.729347 5.55095 -5.39122\n"
		"{PHE 2  CE2 (aroC)}P -0.638233 5.30763 -2.02477\n"
		"{PHE 2  CZ  (aroC)} -0.183465 5.39921 -3.32499\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.638233 5.30763 -2.02477\n"
		"{PHE 2  HE2 (Haro)} 'h' 0.0739382 5.24435 -1.20225\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.183465 5.39921 -3.32499\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.886259 5.40712 -3.52927\n" );
	kins.push_back(
		"@group { rot7 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.97259 5.38959 -3.0169\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.97002 5.9862 -3.20066\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.38208 5.81795 -1.5311\n"
		"{PHE 2  CG  (aroC)}P -2.97259 5.38959 -3.0169\n"
		"{PHE 2  CD1 (aroC)} -2.67734 5.06851 -4.33348\n"
		"{PHE 2  CG  (aroC)}P -2.97259 5.38959 -3.0169\n"
		"{PHE 2  CD2 (aroC)} -1.93151 5.73457 -2.16849\n"
		"{PHE 2  CD1 (aroC)}P -2.67734 5.06851 -4.33348\n"
		"{PHE 2  CE1 (aroC)} -1.37329 5.09227 -4.79097\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.67734 5.06851 -4.33348\n"
		"{PHE 2  HD1 (Haro)} 'h' -3.488 4.79591 -5.00975\n"
		"{PHE 2  CD2 (aroC)}P -1.93151 5.73457 -2.16849\n"
		"{PHE 2  CE2 (aroC)} -0.627463 5.76037 -2.62318\n"
		"{PHE 2  CD2 (aroC)}P 'h' -1.93151 5.73457 -2.16849\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.15198 5.98854 -1.13104\n"
		"{PHE 2  CE1 (aroC)}P -1.37329 5.09227 -4.79097\n"
		"{PHE 2  CZ  (aroC)} -0.34849 5.43835 -3.93628\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.37329 5.09227 -4.79097\n"
		"{PHE 2  HE1 (Haro)} 'h' -1.15548 4.83753 -5.82787\n"
		"{PHE 2  CE2 (aroC)}P -0.627463 5.76037 -2.62318\n"
		"{PHE 2  CZ  (aroC)} -0.34849 5.43835 -3.93628\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.627463 5.76037 -2.62318\n"
		"{PHE 2  HE2 (Haro)} 'h' 0.181586 6.03387 -1.94616\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.34849 5.43835 -3.93628\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.679449 5.45636 -4.29561\n" );
	kins.push_back(
		"@group { rot8 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.97259 5.38959 -3.0169\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.97002 5.9862 -3.20066\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.38208 5.81795 -1.5311\n"
		"{PHE 2  CG  (aroC)}P -2.97259 5.38959 -3.0169\n"
		"{PHE 2  CD1 (aroC)} -2.69605 5.31388 -4.3739\n"
		"{PHE 2  CG  (aroC)}P -2.97259 5.38959 -3.0169\n"
		"{PHE 2  CD2 (aroC)} -1.91287 5.48933 -2.12825\n"
		"{PHE 2  CD1 (aroC)}P -2.69605 5.31388 -4.3739\n"
		"{PHE 2  CE1 (aroC)} -1.39202 5.33788 -4.83146\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.69605 5.31388 -4.3739\n"
		"{PHE 2  HD1 (Haro)} 'h' -3.52145 5.2354 -5.08197\n"
		"{PHE 2  CD2 (aroC)}P -1.91287 5.48933 -2.12825\n"
		"{PHE 2  CE2 (aroC)} -0.60872 5.51477 -2.58267\n"
		"{PHE 2  CD2 (aroC)}P 'h' -1.91287 5.48933 -2.12825\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.11857 5.54909 -1.05889\n"
		"{PHE 2  CE1 (aroC)}P -1.39202 5.33788 -4.83146\n"
		"{PHE 2  CZ  (aroC)} -0.348488 5.43857 -3.93626\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.39202 5.33788 -4.83146\n"
		"{PHE 2  HE1 (Haro)} 'h' -1.18904 5.27743 -5.9004\n"
		"{PHE 2  CE2 (aroC)}P -0.60872 5.51477 -2.58267\n"
		"{PHE 2  CZ  (aroC)} -0.348488 5.43857 -3.93626\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.60872 5.51477 -2.58267\n"
		"{PHE 2  HE2 (Haro)} 'h' 0.215133 5.59409 -1.87365\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.348488 5.43857 -3.93626\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.679374 5.45653 -4.29582\n" );
	kins.push_back(
		"@group { rot9 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.97259 5.38959 -3.0169\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.97002 5.9862 -3.20066\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.38208 5.81795 -1.5311\n"
		"{PHE 2  CG  (aroC)}P -2.97259 5.38959 -3.0169\n"
		"{PHE 2  CD1 (aroC)} -2.69764 5.5631 -4.36521\n"
		"{PHE 2  CG  (aroC)}P -2.97259 5.38959 -3.0169\n"
		"{PHE 2  CD2 (aroC)} -1.91133 5.24027 -2.13708\n"
		"{PHE 2  CD1 (aroC)}P -2.69764 5.5631 -4.36521\n"
		"{PHE 2  CE1 (aroC)} -1.39363 5.58736 -4.8228\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.69764 5.5631 -4.36521\n"
		"{PHE 2  HD1 (Haro)} 'h' -3.52425 5.68176 -5.06625\n"
		"{PHE 2  CD2 (aroC)}P -1.91133 5.24027 -2.13708\n"
		"{PHE 2  CE2 (aroC)} -0.607104 5.2653 -2.5913\n"
		"{PHE 2  CD2 (aroC)}P 'h' -1.91133 5.24027 -2.13708\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.11579 5.1028 -1.07467\n"
		"{PHE 2  CE1 (aroC)}P -1.39363 5.58736 -4.8228\n"
		"{PHE 2  CZ  (aroC)} -0.34847 5.43878 -3.9362\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.39363 5.58736 -4.8228\n"
		"{PHE 2  HE1 (Haro)} 'h' -1.19192 5.72425 -5.8849\n"
		"{PHE 2  CE2 (aroC)}P -0.607104 5.2653 -2.5913\n"
		"{PHE 2  CZ  (aroC)} -0.34847 5.43878 -3.9362\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.607104 5.2653 -2.5913\n"
		"{PHE 2  HE2 (Haro)} 'h' 0.218013 5.1474 -1.88915\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.34847 5.43878 -3.9362\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.679312 5.45675 -4.29598\n" );
	kins.push_back(
		"@group { rot10 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.44241 5.97572 -3.8915\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.93001 6.03771 -1.85733\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.34037 5.33498 -2.23237\n"
		"{PHE 2  CG  (aroC)}P -4.44241 5.97572 -3.8915\n"
		"{PHE 2  CD1 (aroC)} -5.60721 6.56569 -4.35931\n"
		"{PHE 2  CG  (aroC)}P -4.44241 5.97572 -3.8915\n"
		"{PHE 2  CD2 (aroC)} -3.32769 5.96269 -4.71606\n"
		"{PHE 2  CD1 (aroC)}P -5.60721 6.56569 -4.35931\n"
		"{PHE 2  CE1 (aroC)} -5.65633 7.1294 -5.62035\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.60721 6.56569 -4.35931\n"
		"{PHE 2  HD1 (Haro)} 'h' -6.48995 6.58199 -3.71953\n"
		"{PHE 2  CD2 (aroC)}P -3.32769 5.96269 -4.71606\n"
		"{PHE 2  CE2 (aroC)} -3.37345 6.52632 -5.97628\n"
		"{PHE 2  CD2 (aroC)}P 'h' -3.32769 5.96269 -4.71606\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.40647 5.50112 -4.35869\n"
		"{PHE 2  CE1 (aroC)}P -5.65633 7.1294 -5.62035\n"
		"{PHE 2  CZ  (aroC)} -4.53978 7.11006 -6.42869\n"
		"{PHE 2  CE1 (aroC)}P 'h' -5.65633 7.1294 -5.62035\n"
		"{PHE 2  HE1 (Haro)} 'h' -6.57816 7.58947 -5.97539\n"
		"{PHE 2  CE2 (aroC)}P -3.37345 6.52632 -5.97628\n"
		"{PHE 2  CZ  (aroC)} -4.53978 7.11006 -6.42869\n"
		"{PHE 2  CE2 (aroC)}P 'h' -3.37345 6.52632 -5.97628\n"
		"{PHE 2  HE2 (Haro)} 'h' -2.48998 6.50994 -6.6142\n"
		"{PHE 2  CZ  (aroC)}P 'h' -4.53978 7.11006 -6.42869\n"
		"{PHE 2  HZ  (Haro)} 'h' -4.57786 7.55309 -7.42287\n" );
	kins.push_back(
		"@group { rot11 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.44241 5.97572 -3.8915\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.93001 6.03771 -1.85733\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.34037 5.33498 -2.23237\n"
		"{PHE 2  CG  (aroC)}P -4.44241 5.97572 -3.8915\n"
		"{PHE 2  CD1 (aroC)} -5.42106 6.90363 -4.21538\n"
		"{PHE 2  CG  (aroC)}P -4.44241 5.97572 -3.8915\n"
		"{PHE 2  CD2 (aroC)} -3.51395 5.62507 -4.85984\n"
		"{PHE 2  CD1 (aroC)}P -5.42106 6.90363 -4.21538\n"
		"{PHE 2  CE1 (aroC)} -5.47005 7.46769 -5.47627\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.42106 6.90363 -4.21538\n"
		"{PHE 2  HD1 (Haro)} 'h' -6.1563 7.18713 -3.46181\n"
		"{PHE 2  CD2 (aroC)}P -3.51395 5.62507 -4.85984\n"
		"{PHE 2  CE2 (aroC)} -3.55969 6.18801 -6.12037\n"
		"{PHE 2  CD2 (aroC)}P 'h' -3.51395 5.62507 -4.85984\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.74017 4.89611 -4.61635\n"
		"{PHE 2  CE1 (aroC)}P -5.47005 7.46769 -5.47627\n"
		"{PHE 2  CZ  (aroC)} -4.53954 7.11032 -6.42859\n"
		"{PHE 2  CE1 (aroC)}P 'h' -5.47005 7.46769 -5.47627\n"
		"{PHE 2  HE1 (Haro)} 'h' -6.24454 8.19539 -5.71732\n"
		"{PHE 2  CE2 (aroC)}P -3.55969 6.18801 -6.12037\n"
		"{PHE 2  CZ  (aroC)} -4.53954 7.11032 -6.42859\n"
		"{PHE 2  CE2 (aroC)}P 'h' -3.55969 6.18801 -6.12037\n"
		"{PHE 2  HE2 (Haro)} 'h' -2.82352 5.90419 -6.87219\n"
		"{PHE 2  CZ  (aroC)}P 'h' -4.53954 7.11032 -6.42859\n"
		"{PHE 2  HZ  (Haro)} 'h' -4.57797 7.55349 -7.42268\n" );
	kins.push_back(
		"@group { rot12 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.44241 5.97572 -3.8915\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.93001 6.03771 -1.85733\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.34037 5.33498 -2.23237\n"
		"{PHE 2  CG  (aroC)}P -4.44241 5.97572 -3.8915\n"
		"{PHE 2  CD1 (aroC)} -5.12119 7.16537 -4.10982\n"
		"{PHE 2  CG  (aroC)}P -4.44241 5.97572 -3.8915\n"
		"{PHE 2  CD2 (aroC)} -3.81381 5.36368 -4.96524\n"
		"{PHE 2  CD1 (aroC)}P -5.12119 7.16537 -4.10982\n"
		"{PHE 2  CE1 (aroC)} -5.16992 7.72973 -5.37059\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.12119 7.16537 -4.10982\n"
		"{PHE 2  HD1 (Haro)} 'h' -5.61904 7.65571 -3.27285\n"
		"{PHE 2  CD2 (aroC)}P -3.81381 5.36368 -4.96524\n"
		"{PHE 2  CE2 (aroC)} -3.85978 5.92595 -6.22606\n"
		"{PHE 2  CD2 (aroC)}P 'h' -3.81381 5.36368 -4.96524\n"
		"{PHE 2  HD2 (Haro)} 'h' -3.27743 4.42767 -4.80524\n"
		"{PHE 2  CE1 (aroC)}P -5.16992 7.72973 -5.37059\n"
		"{PHE 2  CZ  (aroC)} -4.53922 7.11048 -6.42853\n"
		"{PHE 2  CE1 (aroC)}P 'h' -5.16992 7.72973 -5.37059\n"
		"{PHE 2  HE1 (Haro)} 'h' -5.707 8.66473 -5.52803\n"
		"{PHE 2  CE2 (aroC)}P -3.85978 5.92595 -6.22606\n"
		"{PHE 2  CZ  (aroC)} -4.53922 7.11048 -6.42853\n"
		"{PHE 2  CE2 (aroC)}P 'h' -3.85978 5.92595 -6.22606\n"
		"{PHE 2  HE2 (Haro)} 'h' -3.36091 5.43499 -7.06142\n"
		"{PHE 2  CZ  (aroC)}P 'h' -4.53922 7.11048 -6.42853\n"
		"{PHE 2  HZ  (Haro)} 'h' -4.57793 7.55391 -7.4225\n" );
	kins.push_back(
		"@group { rot13 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.68824 6.08572 -3.80373\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.78408 5.98232 -1.71577\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.30554 5.30948 -2.4387\n"
		"{PHE 2  CG  (aroC)}P -4.68824 6.08572 -3.80373\n"
		"{PHE 2  CD1 (aroC)} -5.90978 6.71419 -3.9949\n"
		"{PHE 2  CG  (aroC)}P -4.68824 6.08572 -3.80373\n"
		"{PHE 2  CD2 (aroC)} -3.74857 6.13791 -4.82203\n"
		"{PHE 2  CD1 (aroC)}P -5.90978 6.71419 -3.9949\n"
		"{PHE 2  CE1 (aroC)} -6.18505 7.3791 -5.17496\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.90978 6.71419 -3.9949\n"
		"{PHE 2  HD1 (Haro)} 'h' -6.65628 6.6799 -3.20094\n"
		"{PHE 2  CD2 (aroC)}P -3.74857 6.13791 -4.82203\n"
		"{PHE 2  CE2 (aroC)} -4.02038 6.80267 -6.00193\n"
		"{PHE 2  CD2 (aroC)}P 'h' -3.74857 6.13791 -4.82203\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.78475 5.64691 -4.68279\n"
		"{PHE 2  CE1 (aroC)}P -6.18505 7.3791 -5.17496\n"
		"{PHE 2  CZ  (aroC)} -5.24063 7.42368 -6.17824\n"
		"{PHE 2  CE1 (aroC)}P 'h' -6.18505 7.3791 -5.17496\n"
		"{PHE 2  HE1 (Haro)} 'h' -7.14904 7.86842 -5.31194\n"
		"{PHE 2  CE2 (aroC)}P -4.02038 6.80267 -6.00193\n"
		"{PHE 2  CZ  (aroC)} -5.24063 7.42368 -6.17824\n"
		"{PHE 2  CE2 (aroC)}P 'h' -4.02038 6.80267 -6.00193\n"
		"{PHE 2  HE2 (Haro)} 'h' -3.27282 6.83673 -6.79421\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.24063 7.42368 -6.17824\n"
		"{PHE 2  HZ  (Haro)} 'h' -5.45703 7.94649 -7.1088\n" );
	kins.push_back(
		"@group { rot14 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.68824 6.08572 -3.80373\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.78408 5.98232 -1.71577\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.30554 5.30948 -2.4387\n"
		"{PHE 2  CG  (aroC)}P -4.68824 6.08572 -3.80373\n"
		"{PHE 2  CD1 (aroC)} -5.6941 7.03888 -3.86214\n"
		"{PHE 2  CG  (aroC)}P -4.68824 6.08572 -3.80373\n"
		"{PHE 2  CD2 (aroC)} -3.96432 5.81353 -4.9546\n"
		"{PHE 2  CD1 (aroC)}P -5.6941 7.03888 -3.86214\n"
		"{PHE 2  CE1 (aroC)} -5.9692 7.70413 -5.04205\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.6941 7.03888 -3.86214\n"
		"{PHE 2  HD1 (Haro)} 'h' -6.26976 7.26131 -2.96328\n"
		"{PHE 2  CD2 (aroC)}P -3.96432 5.81353 -4.9546\n"
		"{PHE 2  CE2 (aroC)} -4.23618 6.47762 -6.13486\n"
		"{PHE 2  CD2 (aroC)}P 'h' -3.96432 5.81353 -4.9546\n"
		"{PHE 2  HD2 (Haro)} 'h' -3.17131 5.06562 -4.92037\n"
		"{PHE 2  CE1 (aroC)}P -5.9692 7.70413 -5.04205\n"
		"{PHE 2  CZ  (aroC)} -5.24037 7.42392 -6.17817\n"
		"{PHE 2  CE1 (aroC)}P 'h' -5.9692 7.70413 -5.04205\n"
		"{PHE 2  HE1 (Haro)} 'h' -6.76247 8.45059 -5.07386\n"
		"{PHE 2  CE2 (aroC)}P -4.23618 6.47762 -6.13486\n"
		"{PHE 2  CZ  (aroC)} -5.24037 7.42392 -6.17817\n"
		"{PHE 2  CE2 (aroC)}P 'h' -4.23618 6.47762 -6.13486\n"
		"{PHE 2  HE2 (Haro)} 'h' -3.65929 6.25472 -7.03222\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.24037 7.42392 -6.17817\n"
		"{PHE 2  HZ  (Haro)} 'h' -5.4571 7.94688 -7.10857\n" );
	kins.push_back(
		"@group { rot15 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.68824 6.08572 -3.80373\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.78408 5.98232 -1.71577\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.30554 5.30948 -2.4387\n"
		"{PHE 2  CG  (aroC)}P -4.68824 6.08572 -3.80373\n"
		"{PHE 2  CD1 (aroC)} -5.37526 7.29054 -3.79446\n"
		"{PHE 2  CG  (aroC)}P -4.68824 6.08572 -3.80373\n"
		"{PHE 2  CD2 (aroC)} -4.28311 5.5622 -5.02211\n"
		"{PHE 2  CD1 (aroC)}P -5.37526 7.29054 -3.79446\n"
		"{PHE 2  CE1 (aroC)} -5.65008 7.95608 -4.97427\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.37526 7.29054 -3.79446\n"
		"{PHE 2  HD1 (Haro)} 'h' -5.69855 7.71186 -2.84221\n"
		"{PHE 2  CD2 (aroC)}P -4.28311 5.5622 -5.02211\n"
		"{PHE 2  CE2 (aroC)} -4.55526 6.22565 -6.20266\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.28311 5.5622 -5.02211\n"
		"{PHE 2  HD2 (Haro)} 'h' -3.74251 4.61521 -5.04137\n"
		"{PHE 2  CE1 (aroC)}P -5.65008 7.95608 -4.97427\n"
		"{PHE 2  CZ  (aroC)} -5.24004 7.42408 -6.17816\n"
		"{PHE 2  CE1 (aroC)}P 'h' -5.65008 7.95608 -4.97427\n"
		"{PHE 2  HE1 (Haro)} 'h' -6.19092 8.90186 -4.95246\n"
		"{PHE 2  CE2 (aroC)}P -4.55526 6.22565 -6.20266\n"
		"{PHE 2  CZ  (aroC)} -5.24004 7.42408 -6.17816\n"
		"{PHE 2  CE2 (aroC)}P 'h' -4.55526 6.22565 -6.20266\n"
		"{PHE 2  HE2 (Haro)} 'h' -4.23068 5.80359 -7.15358\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.24004 7.42408 -6.17816\n"
		"{PHE 2  HZ  (Haro)} 'h' -5.45702 7.94728 -7.10836\n" );
	kins.push_back(
		"@group { rot16 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.91129 6.18867 -3.66269\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.61516 5.91526 -1.60915\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.31051 5.30085 -2.64925\n"
		"{PHE 2  CG  (aroC)}P -4.91129 6.18867 -3.66269\n"
		"{PHE 2  CD1 (aroC)} -6.13594 6.83302 -3.56932\n"
		"{PHE 2  CG  (aroC)}P -4.91129 6.18867 -3.66269\n"
		"{PHE 2  CD2 (aroC)} -4.17876 6.32202 -4.83242\n"
		"{PHE 2  CD1 (aroC)}P -6.13594 6.83302 -3.56932\n"
		"{PHE 2  CE1 (aroC)} -6.61633 7.59261 -4.61939\n"
		"{PHE 2  CD1 (aroC)}P 'h' -6.13594 6.83302 -3.56932\n"
		"{PHE 2  HD1 (Haro)} 'h' -6.72056 6.73544 -2.65416\n"
		"{PHE 2  CD2 (aroC)}P -4.17876 6.32202 -4.83242\n"
		"{PHE 2  CE2 (aroC)} -4.65574 7.08145 -5.88299\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.17876 6.32202 -4.83242\n"
		"{PHE 2  HD2 (Haro)} 'h' -3.21458 5.81942 -4.91705\n"
		"{PHE 2  CE1 (aroC)}P -6.61633 7.59261 -4.61939\n"
		"{PHE 2  CZ  (aroC)} -5.87646 7.71715 -5.77592\n"
		"{PHE 2  CE1 (aroC)}P 'h' -6.61633 7.59261 -4.61939\n"
		"{PHE 2  HE1 (Haro)} 'h' -7.58029 8.09335 -4.53263\n"
		"{PHE 2  CE2 (aroC)}P -4.65574 7.08145 -5.88299\n"
		"{PHE 2  CZ  (aroC)} -5.87646 7.71715 -5.77592\n"
		"{PHE 2  CE2 (aroC)}P 'h' -4.65574 7.08145 -5.88299\n"
		"{PHE 2  HE2 (Haro)} 'h' -4.06975 7.17865 -6.79673\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.87646 7.71715 -5.77592\n"
		"{PHE 2  HZ  (Haro)} 'h' -6.25466 8.31464 -6.60422\n" );
	kins.push_back(
		"@group { rot17 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.91129 6.18867 -3.66269\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.61516 5.91526 -1.60915\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.31051 5.30085 -2.64925\n"
		"{PHE 2  CG  (aroC)}P -4.91129 6.18867 -3.66269\n"
		"{PHE 2  CD1 (aroC)} -5.89358 7.14536 -3.45411\n"
		"{PHE 2  CG  (aroC)}P -4.91129 6.18867 -3.66269\n"
		"{PHE 2  CD2 (aroC)} -4.42115 6.00997 -4.9474\n"
		"{PHE 2  CD1 (aroC)}P -5.89358 7.14536 -3.45411\n"
		"{PHE 2  CE1 (aroC)} -6.37378 7.90529 -4.50403\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.89358 7.14536 -3.45411\n"
		"{PHE 2  HD1 (Haro)} 'h' -6.28629 7.29475 -2.44799\n"
		"{PHE 2  CD2 (aroC)}P -4.42115 6.00997 -4.9474\n"
		"{PHE 2  CE2 (aroC)} -4.89825 6.76876 -5.99838\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.42115 6.00997 -4.9474\n"
		"{PHE 2  HD2 (Haro)} 'h' -3.64887 5.26022 -5.12313\n"
		"{PHE 2  CE1 (aroC)}P -6.37378 7.90529 -4.50403\n"
		"{PHE 2  CZ  (aroC)} -5.87619 7.71738 -5.77587\n"
		"{PHE 2  CE1 (aroC)}P 'h' -6.37378 7.90529 -4.50403\n"
		"{PHE 2  HE1 (Haro)} 'h' -7.14587 8.65339 -4.32601\n"
		"{PHE 2  CE2 (aroC)}P -4.89825 6.76876 -5.99838\n"
		"{PHE 2  CZ  (aroC)} -5.87619 7.71738 -5.77587\n"
		"{PHE 2  CE2 (aroC)}P 'h' -4.89825 6.76876 -5.99838\n"
		"{PHE 2  HE2 (Haro)} 'h' -4.50405 6.61878 -7.00329\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.87619 7.71738 -5.77587\n"
		"{PHE 2  HZ  (Haro)} 'h' -6.25468 8.315 -6.60394\n" );
	kins.push_back(
		"@group { rot18 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.91129 6.18867 -3.66269\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.61516 5.91526 -1.60915\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.31051 5.30085 -2.64925\n"
		"{PHE 2  CG  (aroC)}P -4.91129 6.18867 -3.66269\n"
		"{PHE 2  CD1 (aroC)} -5.5634 7.39003 -3.42789\n"
		"{PHE 2  CG  (aroC)}P -4.91129 6.18867 -3.66269\n"
		"{PHE 2  CD2 (aroC)} -4.75124 5.76561 -4.97344\n"
		"{PHE 2  CD1 (aroC)}P -5.5634 7.39003 -3.42789\n"
		"{PHE 2  CE1 (aroC)} -6.0433 8.15023 -4.47775\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.5634 7.39003 -3.42789\n"
		"{PHE 2  HD1 (Haro)} 'h' -5.69481 7.73279 -2.40121\n"
		"{PHE 2  CD2 (aroC)}P -4.75124 5.76561 -4.97344\n"
		"{PHE 2  CE2 (aroC)} -5.2287 6.52379 -6.02469\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.75124 5.76561 -4.97344\n"
		"{PHE 2  HD2 (Haro)} 'h' -4.24033 4.82231 -5.16983\n"
		"{PHE 2  CE1 (aroC)}P -6.0433 8.15023 -4.47775\n"
		"{PHE 2  CZ  (aroC)} -5.87586 7.71754 -5.77591\n"
		"{PHE 2  CE1 (aroC)}P 'h' -6.0433 8.15023 -4.47775\n"
		"{PHE 2  HE1 (Haro)} 'h' -6.55398 9.09212 -4.27891\n"
		"{PHE 2  CE2 (aroC)}P -5.2287 6.52379 -6.02469\n"
		"{PHE 2  CZ  (aroC)} -5.87586 7.71754 -5.77591\n"
		"{PHE 2  CE2 (aroC)}P 'h' -5.2287 6.52379 -6.02469\n"
		"{PHE 2  HE2 (Haro)} 'h' -5.09579 6.18017 -7.05036\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.87586 7.71754 -5.77591\n"
		"{PHE 2  HZ  (Haro)} 'h' -6.25455 8.31538 -6.60373\n" );
	kins.push_back(
		"@group { rot19 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.62778 6.05837 -3.83117\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.82345 5.99754 -1.74851\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.31062 5.31437 -2.38572\n"
		"{PHE 2  CG  (aroC)}P -4.62778 6.05837 -3.83117\n"
		"{PHE 2  CD1 (aroC)} -4.68607 5.33733 -5.01454\n"
		"{PHE 2  CG  (aroC)}P -4.62778 6.05837 -3.83117\n"
		"{PHE 2  CD2 (aroC)} -4.79258 7.43414 -3.88357\n"
		"{PHE 2  CD1 (aroC)}P -4.68607 5.33733 -5.01454\n"
		"{PHE 2  CE1 (aroC)} -4.90436 5.97588 -6.22077\n"
		"{PHE 2  CD1 (aroC)}P 'h' -4.68607 5.33733 -5.01454\n"
		"{PHE 2  HD1 (Haro)} 'h' -4.55809 4.25493 -4.98542\n"
		"{PHE 2  CD2 (aroC)}P -4.79258 7.43414 -3.88357\n"
		"{PHE 2  CE2 (aroC)} -5.01208 8.07488 -5.0874\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.79258 7.43414 -3.88357\n"
		"{PHE 2  HD2 (Haro)} 'h' -4.74826 8.01127 -2.95924\n"
		"{PHE 2  CE1 (aroC)}P -4.90436 5.97588 -6.22077\n"
		"{PHE 2  CZ  (aroC)} -5.06762 7.34427 -6.25739\n"
		"{PHE 2  CE1 (aroC)}P 'h' -4.90436 5.97588 -6.22077\n"
		"{PHE 2  HE1 (Haro)} 'h' -4.94763 5.39709 -7.14306\n"
		"{PHE 2  CE2 (aroC)}P -5.01208 8.07488 -5.0874\n"
		"{PHE 2  CZ  (aroC)} -5.06762 7.34427 -6.25739\n"
		"{PHE 2  CE2 (aroC)}P 'h' -5.01208 8.07488 -5.0874\n"
		"{PHE 2  HE2 (Haro)} 'h' -5.14064 9.15676 -5.11467\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.06762 7.34427 -6.25739\n"
		"{PHE 2  HZ  (Haro)} 'h' -5.23866 7.8484 -7.2075\n" );
	kins.push_back(
		"@group { rot20 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.62778 6.05837 -3.83117\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.82345 5.99754 -1.74851\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.31062 5.31437 -2.38572\n"
		"{PHE 2  CG  (aroC)}P -4.62778 6.05837 -3.83117\n"
		"{PHE 2  CD1 (aroC)} -5.08343 5.3523 -4.93458\n"
		"{PHE 2  CG  (aroC)}P -4.62778 6.05837 -3.83117\n"
		"{PHE 2  CD2 (aroC)} -4.39551 7.41895 -3.96359\n"
		"{PHE 2  CD1 (aroC)}P -5.08343 5.3523 -4.93458\n"
		"{PHE 2  CE1 (aroC)} -5.30212 5.99083 -6.14076\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.08343 5.3523 -4.93458\n"
		"{PHE 2  HD1 (Haro)} 'h' -5.26973 4.282 -4.84209\n"
		"{PHE 2  CD2 (aroC)}P -4.39551 7.41895 -3.96359\n"
		"{PHE 2  CE2 (aroC)} -4.61432 8.05998 -5.16739\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.39551 7.41895 -3.96359\n"
		"{PHE 2  HD2 (Haro)} 'h' -4.03673 7.9841 -3.10261\n"
		"{PHE 2  CE1 (aroC)}P -5.30212 5.99083 -6.14076\n"
		"{PHE 2  CZ  (aroC)} -5.06795 7.34437 -6.25728\n"
		"{PHE 2  CE1 (aroC)}P 'h' -5.30212 5.99083 -6.14076\n"
		"{PHE 2  HE1 (Haro)} 'h' -5.66006 5.42383 -6.99976\n"
		"{PHE 2  CE2 (aroC)}P -4.61432 8.05998 -5.16739\n"
		"{PHE 2  CZ  (aroC)} -5.06795 7.34437 -6.25728\n"
		"{PHE 2  CE2 (aroC)}P 'h' -4.61432 8.05998 -5.16739\n"
		"{PHE 2  HE2 (Haro)} 'h' -4.42842 9.13002 -5.25794\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.06795 7.34437 -6.25728\n"
		"{PHE 2  HZ  (Haro)} 'h' -5.23903 7.84815 -7.20757\n" );
	kins.push_back(
		"@group { rot21 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.62778 6.05837 -3.83117\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.82345 5.99754 -1.74851\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.31062 5.31437 -2.38572\n"
		"{PHE 2  CG  (aroC)}P -4.62778 6.05837 -3.83117\n"
		"{PHE 2  CD1 (aroC)} -5.44104 5.4868 -4.79842\n"
		"{PHE 2  CG  (aroC)}P -4.62778 6.05837 -3.83117\n"
		"{PHE 2  CD2 (aroC)} -4.03824 7.28433 -4.09975\n"
		"{PHE 2  CD1 (aroC)}P -5.44104 5.4868 -4.79842\n"
		"{PHE 2  CE1 (aroC)} -5.66011 6.12541 -6.00448\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.44104 5.4868 -4.79842\n"
		"{PHE 2  HD1 (Haro)} 'h' -5.91008 4.52311 -4.59813\n"
		"{PHE 2  CD2 (aroC)}P -4.03824 7.28433 -4.09975\n"
		"{PHE 2  CE2 (aroC)} -4.25632 7.92544 -5.30364\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.03824 7.28433 -4.09975\n"
		"{PHE 2  HD2 (Haro)} 'h' -3.39651 7.74293 -3.34657\n"
		"{PHE 2  CE1 (aroC)}P -5.66011 6.12541 -6.00448\n"
		"{PHE 2  CZ  (aroC)} -5.06822 7.34456 -6.25713\n"
		"{PHE 2  CE1 (aroC)}P 'h' -5.66011 6.12541 -6.00448\n"
		"{PHE 2  HE1 (Haro)} 'h' -6.30126 5.66487 -6.75569\n"
		"{PHE 2  CE2 (aroC)}P -4.25632 7.92544 -5.30364\n"
		"{PHE 2  CZ  (aroC)} -5.06822 7.34456 -6.25713\n"
		"{PHE 2  CE2 (aroC)}P 'h' -4.25632 7.92544 -5.30364\n"
		"{PHE 2  HE2 (Haro)} 'h' -3.78741 8.88905 -5.50194\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.06822 7.34456 -6.25713\n"
		"{PHE 2  HZ  (Haro)} 'h' -5.23945 7.84801 -7.20756\n" );
	kins.push_back(
		"@group { rot22 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.82701 6.14935 -3.72412\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.68383 5.94283 -1.64655\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.30378 5.30219 -2.5666\n"
		"{PHE 2  CG  (aroC)}P -4.82701 6.14935 -3.72412\n"
		"{PHE 2  CD1 (aroC)} -5.08614 5.51341 -4.92915\n"
		"{PHE 2  CG  (aroC)}P -4.82701 6.14935 -3.72412\n"
		"{PHE 2  CD2 (aroC)} -4.97887 7.52583 -3.65411\n"
		"{PHE 2  CD1 (aroC)}P -5.08614 5.51341 -4.92915\n"
		"{PHE 2  CE1 (aroC)} -5.48786 6.23572 -6.03698\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.08614 5.51341 -4.92915\n"
		"{PHE 2  HD1 (Haro)} 'h' -4.97028 4.43131 -4.99587\n"
		"{PHE 2  CD2 (aroC)}P -4.97887 7.52583 -3.65411\n"
		"{PHE 2  CE2 (aroC)} -5.38135 8.25015 -4.75923\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.97887 7.52583 -3.65411\n"
		"{PHE 2  HD2 (Haro)} 'h' -4.77746 8.03646 -2.71171\n"
		"{PHE 2  CE1 (aroC)}P -5.48786 6.23572 -6.03698\n"
		"{PHE 2  CZ  (aroC)} -5.63575 7.60372 -5.9522\n"
		"{PHE 2  CE1 (aroC)}P 'h' -5.48786 6.23572 -6.03698\n"
		"{PHE 2  HE1 (Haro)} 'h' -5.68793 5.72329 -6.97766\n"
		"{PHE 2  CE2 (aroC)}P -5.38135 8.25015 -4.75923\n"
		"{PHE 2  CZ  (aroC)} -5.63575 7.60372 -5.9522\n"
		"{PHE 2  CE2 (aroC)}P 'h' -5.38135 8.25015 -4.75923\n"
		"{PHE 2  HE2 (Haro)} 'h' -5.4975 9.3316 -4.69063\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.63575 7.60372 -5.9522\n"
		"{PHE 2  HZ  (Haro)} 'h' -5.95125 8.17381 -6.82487\n" );
	kins.push_back(
		"@group { rot23 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.82701 6.14935 -3.72412\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.68383 5.94283 -1.64655\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.30378 5.30219 -2.5666\n"
		"{PHE 2  CG  (aroC)}P -4.82701 6.14935 -3.72412\n"
		"{PHE 2  CD1 (aroC)} -5.46507 5.52376 -4.78486\n"
		"{PHE 2  CG  (aroC)}P -4.82701 6.14935 -3.72412\n"
		"{PHE 2  CD2 (aroC)} -4.60024 7.51526 -3.79842\n"
		"{PHE 2  CD1 (aroC)}P -5.46507 5.52376 -4.78486\n"
		"{PHE 2  CE1 (aroC)} -5.86717 6.24604 -5.89258\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.46507 5.52376 -4.78486\n"
		"{PHE 2  HD1 (Haro)} 'h' -5.64889 4.45009 -4.73732\n"
		"{PHE 2  CD2 (aroC)}P -4.60024 7.51526 -3.79842\n"
		"{PHE 2  CE2 (aroC)} -5.00203 8.23988 -4.90361\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.60024 7.51526 -3.79842\n"
		"{PHE 2  HD2 (Haro)} 'h' -4.09897 8.01758 -2.97028\n"
		"{PHE 2  CE1 (aroC)}P -5.86717 6.24604 -5.89258\n"
		"{PHE 2  CZ  (aroC)} -5.63605 7.60381 -5.95203\n"
		"{PHE 2  CE1 (aroC)}P 'h' -5.86717 6.24604 -5.89258\n"
		"{PHE 2  HE1 (Haro)} 'h' -6.36732 5.74175 -6.71903\n"
		"{PHE 2  CE2 (aroC)}P -5.00203 8.23988 -4.90361\n"
		"{PHE 2  CZ  (aroC)} -5.63605 7.60381 -5.95203\n"
		"{PHE 2  CE2 (aroC)}P 'h' -5.00203 8.23988 -4.90361\n"
		"{PHE 2  HE2 (Haro)} 'h' -4.8183 9.31314 -4.94919\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.63605 7.60381 -5.95203\n"
		"{PHE 2  HZ  (Haro)} 'h' -5.95163 8.17356 -6.82489\n" );
	kins.push_back(
		"@group { rot24 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -4.82701 6.14935 -3.72412\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.68383 5.94283 -1.64655\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.30378 5.30219 -2.5666\n"
		"{PHE 2  CG  (aroC)}P -4.82701 6.14935 -3.72412\n"
		"{PHE 2  CD1 (aroC)} -5.79402 5.64929 -4.58348\n"
		"{PHE 2  CG  (aroC)}P -4.82701 6.14935 -3.72412\n"
		"{PHE 2  CD2 (aroC)} -4.27164 7.38961 -3.99976\n"
		"{PHE 2  CD1 (aroC)}P -5.79402 5.64929 -4.58348\n"
		"{PHE 2  CE1 (aroC)} -6.19647 6.37165 -5.69101\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.79402 5.64929 -4.58348\n"
		"{PHE 2  HD1 (Haro)} 'h' -6.23788 4.67513 -4.37655\n"
		"{PHE 2  CD2 (aroC)}P -4.27164 7.38961 -3.99976\n"
		"{PHE 2  CE2 (aroC)} -4.67272 8.11431 -5.10515\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.27164 7.38961 -3.99976\n"
		"{PHE 2  HD2 (Haro)} 'h' -3.51011 7.79248 -3.33104\n"
		"{PHE 2  CE1 (aroC)}P -6.19647 6.37165 -5.69101\n"
		"{PHE 2  CZ  (aroC)} -5.63629 7.60399 -5.95182\n"
		"{PHE 2  CE1 (aroC)}P 'h' -6.19647 6.37165 -5.69101\n"
		"{PHE 2  HE1 (Haro)} 'h' -6.95712 5.9667 -6.35803\n"
		"{PHE 2  CE2 (aroC)}P -4.67272 8.11431 -5.10515\n"
		"{PHE 2  CZ  (aroC)} -5.63629 7.60399 -5.95182\n"
		"{PHE 2  CE2 (aroC)}P 'h' -4.67272 8.11431 -5.10515\n"
		"{PHE 2  HE2 (Haro)} 'h' -4.22867 9.08825 -5.3101\n"
		"{PHE 2  CZ  (aroC)}P 'h' -5.63629 7.60399 -5.95182\n"
		"{PHE 2  HZ  (Haro)} 'h' -5.95205 8.17343 -6.82482\n" );
	kins.push_back(
		"@group { rot25 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -5.00499 6.23311 -3.58012\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.53032 5.88072 -1.57264\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.32649 5.30273 -2.74658\n"
		"{PHE 2  CG  (aroC)}P -5.00499 6.23311 -3.58012\n"
		"{PHE 2  CD1 (aroC)} -5.46462 5.68409 -4.76797\n"
		"{PHE 2  CG  (aroC)}P -5.00499 6.23311 -3.58012\n"
		"{PHE 2  CD2 (aroC)} -5.12426 7.60165 -3.39167\n"
		"{PHE 2  CD1 (aroC)}P -5.46462 5.68409 -4.76797\n"
		"{PHE 2  CE1 (aroC)} -6.03022 6.48352 -5.74338\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.46462 5.68409 -4.76797\n"
		"{PHE 2  HD1 (Haro)} 'h' -5.37622 4.60906 -4.92704\n"
		"{PHE 2  CD2 (aroC)}P -5.12426 7.60165 -3.39167\n"
		"{PHE 2  CE2 (aroC)} -5.69015 8.40289 -4.36415\n"
		"{PHE 2  CD2 (aroC)}P 'h' -5.12426 7.60165 -3.39167\n"
		"{PHE 2  HD2 (Haro)} 'h' -4.76585 8.04426 -2.46159\n"
		"{PHE 2  CE1 (aroC)}P -6.03022 6.48352 -5.74338\n"
		"{PHE 2  CZ  (aroC)} -6.1433 7.84256 -5.54162\n"
		"{PHE 2  CE1 (aroC)}P 'h' -6.03022 6.48352 -5.74338\n"
		"{PHE 2  HE1 (Haro)} 'h' -6.38705 6.03898 -6.67211\n"
		"{PHE 2  CE2 (aroC)}P -5.69015 8.40289 -4.36415\n"
		"{PHE 2  CZ  (aroC)} -6.1433 7.84256 -5.54162\n"
		"{PHE 2  CE2 (aroC)}P 'h' -5.69015 8.40289 -4.36415\n"
		"{PHE 2  HE2 (Haro)} 'h' -5.77854 9.47715 -4.20323\n"
		"{PHE 2  CZ  (aroC)}P 'h' -6.1433 7.84256 -5.54162\n"
		"{PHE 2  HZ  (Haro)} 'h' -6.5879 8.47338 -6.31006\n" );
	kins.push_back(
		"@group { rot26 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -5.00499 6.23311 -3.58012\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.53032 5.88072 -1.57264\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.32649 5.30273 -2.74658\n"
		"{PHE 2  CG  (aroC)}P -5.00499 6.23311 -3.58012\n"
		"{PHE 2  CD1 (aroC)} -5.81502 5.68527 -4.56368\n"
		"{PHE 2  CG  (aroC)}P -5.00499 6.23311 -3.58012\n"
		"{PHE 2  CD2 (aroC)} -4.77417 7.60027 -3.59596\n"
		"{PHE 2  CD1 (aroC)}P -5.81502 5.68527 -4.56368\n"
		"{PHE 2  CE1 (aroC)} -6.38098 6.48465 -5.53892\n"
		"{PHE 2  CD1 (aroC)}P 'h' -5.81502 5.68527 -4.56368\n"
		"{PHE 2  HD1 (Haro)} 'h' -6.0037 4.61139 -4.56102\n"
		"{PHE 2  CD2 (aroC)}P -4.77417 7.60027 -3.59596\n"
		"{PHE 2  CE2 (aroC)} -5.33938 8.4018 -4.56859\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.77417 7.60027 -3.59596\n"
		"{PHE 2  HD2 (Haro)} 'h' -4.13849 8.04184 -2.82762\n"
		"{PHE 2  CE1 (aroC)}P -6.38098 6.48465 -5.53892\n"
		"{PHE 2  CZ  (aroC)} -6.14357 7.84263 -5.5414\n"
		"{PHE 2  CE1 (aroC)}P 'h' -6.38098 6.48465 -5.53892\n"
		"{PHE 2  HE1 (Haro)} 'h' -7.0153 6.04098 -6.30591\n"
		"{PHE 2  CE2 (aroC)}P -5.33938 8.4018 -4.56859\n"
		"{PHE 2  CZ  (aroC)} -6.14357 7.84263 -5.5414\n"
		"{PHE 2  CE2 (aroC)}P 'h' -5.33938 8.4018 -4.56859\n"
		"{PHE 2  HE2 (Haro)} 'h' -5.15047 9.47513 -4.56933\n"
		"{PHE 2  CZ  (aroC)}P 'h' -6.14357 7.84263 -5.5414\n"
		"{PHE 2  HZ  (Haro)} 'h' -6.58827 8.47314 -6.31004\n" );
	kins.push_back(
		"@group { rot27 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -5.00499 6.23311 -3.58012\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.53032 5.88072 -1.57264\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -3.32649 5.30273 -2.74658\n"
		"{PHE 2  CG  (aroC)}P -5.00499 6.23311 -3.58012\n"
		"{PHE 2  CD1 (aroC)} -6.10524 5.79719 -4.30337\n"
		"{PHE 2  CG  (aroC)}P -5.00499 6.23311 -3.58012\n"
		"{PHE 2  CD2 (aroC)} -4.48428 7.48822 -3.85617\n"
		"{PHE 2  CD1 (aroC)}P -6.10524 5.79719 -4.30337\n"
		"{PHE 2  CE1 (aroC)} -6.67152 6.59664 -5.27837\n"
		"{PHE 2  CD1 (aroC)}P 'h' -6.10524 5.79719 -4.30337\n"
		"{PHE 2  HD1 (Haro)} 'h' -6.52333 4.81206 -4.09472\n"
		"{PHE 2  CD2 (aroC)}P -4.48428 7.48822 -3.85617\n"
		"{PHE 2  CE2 (aroC)} -5.04882 8.28985 -4.82912\n"
		"{PHE 2  CD2 (aroC)}P 'h' -4.48428 7.48822 -3.85617\n"
		"{PHE 2  HD2 (Haro)} 'h' -3.619 7.84111 -3.29389\n"
		"{PHE 2  CE1 (aroC)}P -6.67152 6.59664 -5.27837\n"
		"{PHE 2  CZ  (aroC)} -6.14377 7.8428 -5.54115\n"
		"{PHE 2  CE1 (aroC)}P 'h' -6.67152 6.59664 -5.27837\n"
		"{PHE 2  HE1 (Haro)} 'h' -7.53569 6.24154 -5.83925\n"
		"{PHE 2  CE2 (aroC)}P -5.04882 8.28985 -4.82912\n"
		"{PHE 2  CZ  (aroC)} -6.14377 7.8428 -5.54115\n"
		"{PHE 2  CE2 (aroC)}P 'h' -5.04882 8.28985 -4.82912\n"
		"{PHE 2  HE2 (Haro)} 'h' -4.63023 9.27462 -5.03586\n"
		"{PHE 2  CZ  (aroC)}P 'h' -6.14377 7.8428 -5.54115\n"
		"{PHE 2  HZ  (Haro)} 'h' -6.58868 8.473 -6.30991\n" );
	kins.push_back(
		"@group { rot28 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.90708 5.37467 -2.75876\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.84144 5.92398 -3.33837\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.55904 5.89123 -1.5832\n"
		"{PHE 2  CG  (aroC)}P -2.90708 5.37467 -2.75876\n"
		"{PHE 2  CD1 (aroC)} -2.18565 4.19044 -2.78724\n"
		"{PHE 2  CG  (aroC)}P -2.90708 5.37467 -2.75876\n"
		"{PHE 2  CD2 (aroC)} -2.23012 6.56896 -2.95381\n"
		"{PHE 2  CD1 (aroC)}P -2.18565 4.19044 -2.78724\n"
		"{PHE 2  CE1 (aroC)} -0.820916 4.19979 -3.00595\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.18565 4.19044 -2.78724\n"
		"{PHE 2  HD1 (Haro)} 'h' -2.70749 3.24524 -2.63532\n"
		"{PHE 2  CD2 (aroC)}P -2.23012 6.56896 -2.95381\n"
		"{PHE 2  CE2 (aroC)} -0.866526 6.58132 -3.17381\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.23012 6.56896 -2.95381\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.78725 7.50629 -2.93309\n"
		"{PHE 2  CE1 (aroC)}P -0.820916 4.19979 -3.00595\n"
		"{PHE 2  CZ  (aroC)} -0.161367 5.39479 -3.19947\n"
		"{PHE 2  CE1 (aroC)}P 'h' -0.820916 4.19979 -3.00595\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.266363 3.26194 -3.0256\n"
		"{PHE 2  CE2 (aroC)}P -0.866526 6.58132 -3.17381\n"
		"{PHE 2  CZ  (aroC)} -0.161367 5.39479 -3.19947\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.866526 6.58132 -3.17381\n"
		"{PHE 2  HE2 (Haro)} 'h' -0.34656 7.52689 -3.32632\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.161367 5.39479 -3.19947\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.914113 5.4032 -3.37087\n" );
	kins.push_back(
		"@group { rot29 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.90708 5.37467 -2.75876\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.84144 5.92398 -3.33837\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.55904 5.89123 -1.5832\n"
		"{PHE 2  CG  (aroC)}P -2.90708 5.37467 -2.75876\n"
		"{PHE 2  CD1 (aroC)} -2.26432 4.25899 -3.27432\n"
		"{PHE 2  CG  (aroC)}P -2.90708 5.37467 -2.75876\n"
		"{PHE 2  CD2 (aroC)} -2.1515 6.50014 -2.46709\n"
		"{PHE 2  CD1 (aroC)}P -2.26432 4.25899 -3.27432\n"
		"{PHE 2  CE1 (aroC)} -0.899669 4.26834 -3.49352\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.26432 4.25899 -3.27432\n"
		"{PHE 2  HD1 (Haro)} 'h' -2.84839 3.36835 -3.50763\n"
		"{PHE 2  CD2 (aroC)}P -2.1515 6.50014 -2.46709\n"
		"{PHE 2  CE2 (aroC)} -0.787773 6.51283 -2.68624\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.1515 6.50014 -2.46709\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.64637 7.38306 -2.06092\n"
		"{PHE 2  CE1 (aroC)}P -0.899669 4.26834 -3.49352\n"
		"{PHE 2  CZ  (aroC)} -0.161433 5.39496 -3.19988\n"
		"{PHE 2  CE1 (aroC)}P 'h' -0.899669 4.26834 -3.49352\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.407417 3.38469 -3.89888\n"
		"{PHE 2  CE2 (aroC)}P -0.787773 6.51283 -2.68624\n"
		"{PHE 2  CZ  (aroC)} -0.161433 5.39496 -3.19988\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.787773 6.51283 -2.68624\n"
		"{PHE 2  HE2 (Haro)} 'h' -0.205547 7.40416 -2.45329\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.161433 5.39496 -3.19988\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.914042 5.40289 -3.37132\n" );
	kins.push_back(
		"@group { rot30 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.90708 5.37467 -2.75876\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.84144 5.92398 -3.33837\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.55904 5.89123 -1.5832\n"
		"{PHE 2  CG  (aroC)}P -2.90708 5.37467 -2.75876\n"
		"{PHE 2  CD1 (aroC)} -2.33321 4.52307 -3.69103\n"
		"{PHE 2  CG  (aroC)}P -2.90708 5.37467 -2.75876\n"
		"{PHE 2  CD2 (aroC)} -2.08268 6.23596 -2.05081\n"
		"{PHE 2  CD1 (aroC)}P -2.33321 4.52307 -3.69103\n"
		"{PHE 2  CE1 (aroC)} -0.968636 4.53263 -3.91068\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.33321 4.52307 -3.69103\n"
		"{PHE 2  HD1 (Haro)} 'h' -2.97175 3.84161 -4.25378\n"
		"{PHE 2  CD2 (aroC)}P -2.08268 6.23596 -2.05081\n"
		"{PHE 2  CE2 (aroC)} -0.718804 6.2486 -2.26905\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.08268 6.23596 -2.05081\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.52303 6.90975 -1.31494\n"
		"{PHE 2  CE1 (aroC)}P -0.968636 4.53263 -3.91068\n"
		"{PHE 2  CZ  (aroC)} -0.161483 5.39528 -3.20018\n"
		"{PHE 2  CE1 (aroC)}P 'h' -0.968636 4.53263 -3.91068\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.530943 3.85801 -4.64607\n"
		"{PHE 2  CE2 (aroC)}P -0.718804 6.2486 -2.26905\n"
		"{PHE 2  CZ  (aroC)} -0.161483 5.39528 -3.20018\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.718804 6.2486 -2.26905\n"
		"{PHE 2  HE2 (Haro)} 'h' -0.0820568 6.93096 -1.70632\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.161483 5.39528 -3.20018\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.913956 5.4028 -3.37186\n" );
	kins.push_back(
		"@group { rot31 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.95924 5.38589 -2.97714\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.95146 5.97704 -3.22416\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.41049 5.82985 -1.53682\n"
		"{PHE 2  CG  (aroC)}P -2.95924 5.38589 -2.97714\n"
		"{PHE 2  CD1 (aroC)} -2.26499 4.20397 -3.18872\n"
		"{PHE 2  CG  (aroC)}P -2.95924 5.38589 -2.97714\n"
		"{PHE 2  CD2 (aroC)} -2.30415 6.5884 -3.19504\n"
		"{PHE 2  CD1 (aroC)}P -2.26499 4.20397 -3.18872\n"
		"{PHE 2  CE1 (aroC)} -0.948234 4.22363 -3.60844\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.26499 4.20397 -3.18872\n"
		"{PHE 2  HD1 (Haro)} 'h' -2.77005 3.25238 -3.02077\n"
		"{PHE 2  CD2 (aroC)}P -2.30415 6.5884 -3.19504\n"
		"{PHE 2  CE2 (aroC)} -0.988667 6.61115 -3.61566\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.30415 6.5884 -3.19504\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.84018 7.52398 -3.03145\n"
		"{PHE 2  CE1 (aroC)}P -0.948234 4.22363 -3.60844\n"
		"{PHE 2  CZ  (aroC)} -0.310123 5.4268 -3.82217\n"
		"{PHE 2  CE1 (aroC)}P 'h' -0.948234 4.22363 -3.60844\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.414592 3.28747 -3.77063\n"
		"{PHE 2  CE2 (aroC)}P -0.988667 6.61115 -3.61566\n"
		"{PHE 2  CZ  (aroC)} -0.310123 5.4268 -3.82217\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.988667 6.61115 -3.61566\n"
		"{PHE 2  HE2 (Haro)} 'h' -0.48555 7.56315 -3.78389\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.310123 5.4268 -3.82217\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.727711 5.44327 -4.15192\n" );
	kins.push_back(
		"@group { rot32 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.95924 5.38589 -2.97714\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.95146 5.97704 -3.22416\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.41049 5.82985 -1.53682\n"
		"{PHE 2  CG  (aroC)}P -2.95924 5.38589 -2.97714\n"
		"{PHE 2  CD1 (aroC)} -2.41468 4.30416 -3.65314\n"
		"{PHE 2  CG  (aroC)}P -2.95924 5.38589 -2.97714\n"
		"{PHE 2  CD2 (aroC)} -2.15458 6.48797 -2.73097\n"
		"{PHE 2  CD1 (aroC)}P -2.41468 4.30416 -3.65314\n"
		"{PHE 2  CE1 (aroC)} -1.09807 4.32385 -4.07334\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.41468 4.30416 -3.65314\n"
		"{PHE 2  HD1 (Haro)} 'h' -3.03812 3.43215 -3.85249\n"
		"{PHE 2  CD2 (aroC)}P -2.15458 6.48797 -2.73097\n"
		"{PHE 2  CE2 (aroC)} -0.838831 6.51099 -3.15076\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.15458 6.48797 -2.73097\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.57216 7.3441 -2.19987\n"
		"{PHE 2  CE1 (aroC)}P -1.09807 4.32385 -4.07334\n"
		"{PHE 2  CZ  (aroC)} -0.310246 5.427 -3.82255\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.09807 4.32385 -4.07334\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.682961 3.46694 -4.6033\n"
		"{PHE 2  CE2 (aroC)}P -0.838831 6.51099 -3.15076\n"
		"{PHE 2  CZ  (aroC)} -0.310246 5.427 -3.82255\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.838831 6.51099 -3.15076\n"
		"{PHE 2  HE2 (Haro)} 'h' -0.217258 7.38372 -2.95145\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.310246 5.427 -3.82255\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.727569 5.44299 -4.15238\n" );
	kins.push_back(
		"@group { rot33 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -2.95924 5.38589 -2.97714\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -4.95146 5.97704 -3.22416\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.41049 5.82985 -1.53682\n"
		"{PHE 2  CG  (aroC)}P -2.95924 5.38589 -2.97714\n"
		"{PHE 2  CD1 (aroC)} -2.54174 4.59487 -4.03717\n"
		"{PHE 2  CG  (aroC)}P -2.95924 5.38589 -2.97714\n"
		"{PHE 2  CD2 (aroC)} -2.02765 6.19718 -2.34737\n"
		"{PHE 2  CD1 (aroC)}P -2.54174 4.59487 -4.03717\n"
		"{PHE 2  CE1 (aroC)} -1.22527 4.6148 -4.45778\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.54174 4.59487 -4.03717\n"
		"{PHE 2  HD1 (Haro)} 'h' -3.26563 3.9531 -4.54008\n"
		"{PHE 2  CD2 (aroC)}P -2.02765 6.19718 -2.34737\n"
		"{PHE 2  CE2 (aroC)} -0.711625 6.2201 -2.76629\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.02765 6.19718 -2.34737\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.3447 6.82311 -1.51244\n"
		"{PHE 2  CE1 (aroC)}P -1.22527 4.6148 -4.45778\n"
		"{PHE 2  CZ  (aroC)} -0.310337 5.42734 -3.82282\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.22527 4.6148 -4.45778\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.910789 3.98802 -5.29188\n"
		"{PHE 2  CE2 (aroC)}P -0.711625 6.2201 -2.76629\n"
		"{PHE 2  CZ  (aroC)} -0.310337 5.42734 -3.82282\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.711625 6.2201 -2.76629\n"
		"{PHE 2  HE2 (Haro)} 'h' 0.0105032 6.86277 -2.26307\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.310337 5.42734 -3.82282\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.727403 5.44294 -4.1529\n" );
	kins.push_back(
		"@group { rot34 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -3.04329 5.41138 -3.18407\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -5.04252 6.0229 -3.09162\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.25747 5.76517 -1.51723\n"
		"{PHE 2  CG  (aroC)}P -3.04329 5.41138 -3.18407\n"
		"{PHE 2  CD1 (aroC)} -2.40325 4.24371 -3.57209\n"
		"{PHE 2  CG  (aroC)}P -3.04329 5.41138 -3.18407\n"
		"{PHE 2  CD2 (aroC)} -2.41313 6.62364 -3.42071\n"
		"{PHE 2  CD1 (aroC)}P -2.40325 4.24371 -3.57209\n"
		"{PHE 2  CE1 (aroC)} -1.16381 4.28682 -4.18229\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.40325 4.24371 -3.57209\n"
		"{PHE 2  HD1 (Haro)} 'h' -2.88939 3.28466 -3.39125\n"
		"{PHE 2  CD2 (aroC)}P -2.41313 6.62364 -3.42071\n"
		"{PHE 2  CE2 (aroC)} -1.17505 6.66989 -4.03141\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.41313 6.62364 -3.42071\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.90699 7.54816 -3.11943\n"
		"{PHE 2  CE1 (aroC)}P -1.16381 4.28682 -4.18229\n"
		"{PHE 2  CZ  (aroC)} -0.549811 5.49951 -4.41223\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.16381 4.28682 -4.18229\n"
		"{PHE 2  HE1 (Haro)} 'h' -0.672111 3.36163 -4.48187\n"
		"{PHE 2  CE2 (aroC)}P -1.17505 6.66989 -4.03141\n"
		"{PHE 2  CZ  (aroC)} -0.549811 5.49951 -4.41223\n"
		"{PHE 2  CE2 (aroC)}P 'h' -1.17505 6.66989 -4.03141\n"
		"{PHE 2  HE2 (Haro)} 'h' -0.690864 7.62937 -4.21221\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.549811 5.49951 -4.41223\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.427249 5.53439 -4.89207\n" );
	kins.push_back(
		"@group { rot35 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -3.04329 5.41138 -3.18407\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -5.04252 6.0229 -3.09162\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.25747 5.76517 -1.51723\n"
		"{PHE 2  CG  (aroC)}P -3.04329 5.41138 -3.18407\n"
		"{PHE 2  CD1 (aroC)} -2.61936 4.3742 -4.00151\n"
		"{PHE 2  CG  (aroC)}P -3.04329 5.41138 -3.18407\n"
		"{PHE 2  CD2 (aroC)} -2.19718 6.49293 -2.99164\n"
		"{PHE 2  CD1 (aroC)}P -2.61936 4.3742 -4.00151\n"
		"{PHE 2  CE1 (aroC)} -1.38014 4.41737 -4.61215\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.61936 4.3742 -4.00151\n"
		"{PHE 2  HD1 (Haro)} 'h' -3.27642 3.5187 -4.16026\n"
		"{PHE 2  CD2 (aroC)}P -2.19718 6.49293 -2.99164\n"
		"{PHE 2  CE2 (aroC)} -0.95872 6.5394 -3.60155\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.19718 6.49293 -2.99164\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.52003 7.31402 -2.35056\n"
		"{PHE 2  CE1 (aroC)}P -1.38014 4.41737 -4.61215\n"
		"{PHE 2  CZ  (aroC)} -0.549987 5.49973 -4.41258\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.38014 4.41737 -4.61215\n"
		"{PHE 2  HE1 (Haro)} 'h' -1.05958 3.59543 -5.25179\n"
		"{PHE 2  CE2 (aroC)}P -0.95872 6.5394 -3.60155\n"
		"{PHE 2  CZ  (aroC)} -0.549987 5.49973 -4.41258\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.95872 6.5394 -3.60155\n"
		"{PHE 2  HE2 (Haro)} 'h' -0.303502 7.39563 -3.44251\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.549987 5.49973 -4.41258\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.427038 5.53414 -4.89251\n" );
	kins.push_back(
		"@group { rot36 } animate dominant\n"
		"@vectorlist {PHE 1} color= bluetint master= {rotamers}\n"
		"{PHE 2  N   (Nbb)}P -6.414 4.03 -2.127\n"
		"{PHE 2  CA  (CAbb)} -4.993 3.959 -2.449\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  C   (CObb)} -4.238 3.117 -1.429\n"
		"{PHE 2  CA  (CAbb)}P -4.993 3.959 -2.449\n"
		"{PHE 2  CB  (CH2)} -4.3903 5.36323 -2.52079\n"
		"{PHE 2  C   (CObb)}P -4.238 3.117 -1.429\n"
		"{PHE 2  O   (OCbb)} -3.492 2.207 -1.791\n"
		"{PHE 2  CB  (CH2)}P -4.3903 5.36323 -2.52079\n"
		"{PHE 2  CG  (aroC)} -3.04329 5.41138 -3.18407\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 1HB  (Hapo)} 'h' -5.04252 6.0229 -3.09162\n"
		"{PHE 2  CB  (CH2)}P 'h' -4.3903 5.36323 -2.52079\n"
		"{PHE 2 2HB  (Hapo)} 'h' -4.25747 5.76517 -1.51723\n"
		"{PHE 2  CG  (aroC)}P -3.04329 5.41138 -3.18407\n"
		"{PHE 2  CD1 (aroC)} -2.79872 4.68954 -4.34286\n"
		"{PHE 2  CG  (aroC)}P -3.04329 5.41138 -3.18407\n"
		"{PHE 2  CD2 (aroC)} -2.01802 6.17754 -2.65069\n"
		"{PHE 2  CD1 (aroC)}P -2.79872 4.68954 -4.34286\n"
		"{PHE 2  CE1 (aroC)} -1.5597 4.73298 -4.95389\n"
		"{PHE 2  CD1 (aroC)}P 'h' -2.79872 4.68954 -4.34286\n"
		"{PHE 2  HD1 (Haro)} 'h' -3.59756 4.08374 -4.77143\n"
		"{PHE 2  CD2 (aroC)}P -2.01802 6.17754 -2.65069\n"
		"{PHE 2  CE2 (aroC)} -0.779155 6.22385 -3.25978\n"
		"{PHE 2  CD2 (aroC)}P 'h' -2.01802 6.17754 -2.65069\n"
		"{PHE 2  HD2 (Haro)} 'h' -2.19897 6.74895 -1.73955\n"
		"{PHE 2  CE1 (aroC)}P -1.5597 4.73298 -4.95389\n"
		"{PHE 2  CZ  (aroC)} -0.550112 5.50008 -4.4128\n"
		"{PHE 2  CE1 (aroC)}P 'h' -1.5597 4.73298 -4.95389\n"
		"{PHE 2  HE1 (Haro)} 'h' -1.38119 4.16067 -5.86388\n"
		"{PHE 2  CE2 (aroC)}P -0.779155 6.22385 -3.25978\n"
		"{PHE 2  CZ  (aroC)} -0.550112 5.50008 -4.4128\n"
		"{PHE 2  CE2 (aroC)}P 'h' -0.779155 6.22385 -3.25978\n"
		"{PHE 2  HE2 (Haro)} 'h' 0.018005 6.83053 -2.8306\n"
		"{PHE 2  CZ  (aroC)}P 'h' -0.550112 5.50008 -4.4128\n"
		"{PHE 2  HZ  (Haro)} 'h' 0.426795 5.53412 -4.89301\n" );

	return kins;
}
