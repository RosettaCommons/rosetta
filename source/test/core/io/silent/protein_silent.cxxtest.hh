// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/silent/protein_silent.cxxtest.hh
/// @brief  test suite for protein silent-file format
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/rms_util.hh>

#include <core/scoring/ScoreType.hh>

#include <utility/file/file_sys_util.hh>

#include <string>

//Auto Headers
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/symmetry/SymDof.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/io/silent/EnergyNames.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/util.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


static basic::Tracer TR("test.core.io.silent.protein_silent");

using namespace core;

class ProteinSilentTests : public CxxTest::TestSuite {
public:
	ProteinSilentTests() {}

	// shared data
	pose::PoseOP start_pose;
	core::chemical::ResidueTypeSetCOP rsd_set;

	pose::PoseOP centroid_pose_;
	core::chemical::ResidueTypeSetCOP cen_rsd_set_;

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-mute core.io.pdb -mute core.conformation");

		rsd_set =
			core::chemical::ChemicalManager::get_instance()->residue_type_set(
			"fa_standard"
		);
		cen_rsd_set_ =
			core::chemical::ChemicalManager::get_instance()->residue_type_set(
			"centroid"
		);
		start_pose = core::import_pose::pose_from_file("core/io/test_in_idealized.pdb", core::import_pose::PDB_file);
		centroid_pose_ = core::import_pose::pose_from_file(*cen_rsd_set_, "core/io/test_in_idealized.pdb", core::import_pose::PDB_file);

	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_save_and_restore() {
		using namespace core::io::silent;

		// configuration information for tests

		double const RMS_ERROR( 1e-3 );
		double const CHI_ERROR( 1e-2 );
		double const BB_ERROR ( 1e-2 );
		std::string const silent_outfile( "test.silent.out" );
		pose::Pose restored_pose;

		utility::file::file_delete( silent_outfile );

		core::io::silent::SilentFileData sfd;

		core::io::silent::ProteinSilentStruct pss( *start_pose, "tag", true );
		sfd.write_silent_struct( pss, silent_outfile );

		sfd.read_file( silent_outfile );
		core::io::silent::SilentFileData::iterator iter = sfd.begin();
		TS_ASSERT( iter->decoy_tag() == "tag" );
		iter->fill_pose( restored_pose, *rsd_set );

		TS_ASSERT( start_pose->total_residue() == restored_pose.total_residue() );
		for ( Size seqpos = 1; seqpos <= restored_pose.total_residue(); ++seqpos ) {
			TS_ASSERT_DELTA(
				start_pose->phi( seqpos ), restored_pose.phi( seqpos ),
				BB_ERROR
			);
			TS_ASSERT_DELTA(
				start_pose->psi( seqpos ), restored_pose.psi( seqpos ),
				BB_ERROR
			);
			TS_ASSERT_DELTA(
				start_pose->omega( seqpos ), restored_pose.omega( seqpos ),
				BB_ERROR
			);
			for ( Size chi_idx = 1; chi_idx <= start_pose->residue_type(seqpos).nchi();
					++chi_idx
					) {
				TS_ASSERT_DELTA(
					start_pose->chi( chi_idx, seqpos ),
					restored_pose.chi( chi_idx, seqpos ), CHI_ERROR
				);
			}
		}

		core::io::silent::ProteinSilentStruct original( *start_pose, "start", false );
		original.fill_pose( *start_pose, *rsd_set );

		Real rms_to_restored = scoring::CA_rmsd( *start_pose, restored_pose );
		TR << "RMS error from save/restore: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < RMS_ERROR );

		utility::file::file_delete( silent_outfile );
	} // test_save_and_restore

	void test_strict_column_mode() {
		// test strict_column_mode for the SilentFileData class, which I added as an
		// expression of affection for Oliver Lange.
		using namespace core::io::silent;

		std::string const silent_outfile( "test.silent.out" );
		utility::file::file_delete( silent_outfile );

		SilentStructOP ss1( new ProteinSilentStruct( *start_pose, "ss1", false ) );
		SilentStructOP ss2( new ProteinSilentStruct( *start_pose, "ss2", false ) );
		ss1->add_energy( "energy1",     1.0 );
		ss1->add_energy( "energy2",     2.0 );
		ss2->add_energy( "irrelevant", 10.0 );
		ss2->add_energy( "energy2",     5.0 );

		SilentFileData sfd;
		sfd.strict_column_mode( true );

		sfd.write_silent_struct( *ss1, silent_outfile );
		sfd.write_silent_struct( *ss2, silent_outfile );

		sfd.clear();
		sfd.read_file( silent_outfile );

		SilentFileData::iterator restore1 = sfd.begin();
		SilentFileData::iterator restore2 = sfd.begin();
		++restore2;

		TS_ASSERT(  restore1->get_energy("energy1") == 1.0 );
		TS_ASSERT(  restore2->get_energy("energy1") == 0.0 );
		TS_ASSERT(  restore2->get_energy("energy2") == 5.0 );
		TS_ASSERT( !restore2->has_energy("irrelevant")     );
		utility::file::file_delete( silent_outfile );

	}

	void test_arbitrary_energies() {
		using namespace core::io::silent;

		std::string const silent_outfile( "test.silent.out" );
		utility::file::file_delete( silent_outfile );

		core::pose::Pose restored_pose;
		SilentStructOP ss1( new ProteinSilentStruct( *start_pose, "ss1", false ) );
		SilentStructOP ss2( new ProteinSilentStruct( *start_pose, "ss2", false ) );

		// test initialization of arbitrary energies into the pose
		ss1->clear_energies();
		ss2->clear_energies();

		ss1->add_energy( "vdw"      , 1.0 );
		ss1->add_energy( "arbitrary", 5.0 );
		ss1->fill_pose( restored_pose );

		core::scoring::EnergyMap & emap( restored_pose.energies().total_energies() );
		TS_ASSERT( emap[ core::scoring::vdw ] == 1.0 );

		ss1->fill_struct( restored_pose );
		TS_ASSERT( ss1->has_energy("vdw")       );
		TS_ASSERT( ss1->has_energy("arbitrary") );
		TS_ASSERT( ss1->get_energy("vdw")       == 1.0 );
		TS_ASSERT( ss1->get_energy("arbitrary") == 5.0 );

		utility::file::file_delete( silent_outfile );
	}

	void test_redundant_tags() {
		// test inserting structures into the SFD class that have redundant tags. Do
		// this by writing the same SilentStruct twice, reading it in, and checking
		// that the two SilentStruct objects are equal except for their tags.
		using namespace core::io::silent;

		std::string const silent_outfile( "test.silent.out" );
		utility::file::file_delete( silent_outfile );

		core::io::silent::SilentFileData sfd;
		core::pose::Pose restored_pose;
		SilentStructOP ss1( new ProteinSilentStruct( *start_pose, "ss1", false ) );
		SilentStructOP ss2( new ProteinSilentStruct( *start_pose, "ss2", false ) );

		utility::file::file_delete( silent_outfile );
		sfd.write_silent_struct( *ss1, silent_outfile );
		sfd.write_silent_struct( *ss1, silent_outfile );
		sfd.clear();

		sfd.read_file( silent_outfile );
		TS_ASSERT( sfd.size() == 2 );
		SilentFileData::iterator it1 = sfd.begin();
		SilentFileData::iterator it2 = sfd.begin();
		++it2;

		TS_ASSERT( it1->decoy_tag() != it2->decoy_tag() );
		utility::file::file_delete( silent_outfile );
	} // test_redundant_tags

	void test_score_cut() {
		using core::Size;
		using core::Real;
		using namespace core::io::silent;

		Real const filter ( 0.1 );
		Size const nstruct( 20 );
		std::string const silent_outfile( "test.silent.out" );
		utility::file::file_delete( silent_outfile );
		SilentFileData sfd;

		for ( Size ii = 1; ii <= nstruct; ++ii ) {
			std::string const tag( "S" + string_of(ii) );
			SilentStructOP ss( new ProteinSilentStruct( *start_pose, tag, false ) );
			ss->add_energy( "score", ii );
			sfd.add_structure( ss );
		}
		sfd.write_all( silent_outfile );
		sfd.clear();

		sfd.read_file( silent_outfile );
		TS_ASSERT( sfd.size() == nstruct );
		sfd.score_filter( filter );
		TS_ASSERT( sfd.size() == 2 );
		SilentFileData::iterator it = sfd.begin();
		TS_ASSERT( it->get_energy( "score" ) == 1.0 );
		++it;
		TS_ASSERT( it->get_energy( "score" ) == 2.0 );
		utility::file::file_delete( silent_outfile );
	} // test_score_cuts

	void test_reverse_score_cut() {
		using core::Size;
		using core::Real;
		using namespace core::io::silent;

		Real const filter ( -0.1 );
		Size const nstruct( 20 );
		std::string const silent_outfile( "test.silent.out" );
		utility::file::file_delete( silent_outfile );
		SilentFileData sfd;

		for ( Size ii = 1; ii <= nstruct; ++ii ) {
			std::string const tag( "S" + string_of(ii) );
			SilentStructOP ss( new ProteinSilentStruct( *start_pose, tag, false ) );
			ss->add_energy( "score", ii );
			sfd.add_structure( ss );
		}
		sfd.write_all( silent_outfile );
		sfd.clear();

		sfd.read_file( silent_outfile );
		TS_ASSERT( sfd.size() == nstruct );
		sfd.reverse_score_filter( filter );
		TS_ASSERT( sfd.size() == 2 );
		SilentFileData::iterator it = sfd.begin();
		TS_ASSERT( it->get_energy( "score" ) == 19.0 );
		++it;
		TS_ASSERT( it->get_energy( "score" ) == 20.0 );
		utility::file::file_delete( silent_outfile );
	} // test_reverse_score_cut

	void test_string_into_pose() {
		using namespace core::pose;
		using namespace core::io::silent;
		core::pose::Pose restored_pose;
		SilentStructOP ss( new ProteinSilentStruct );

		ss->add_string_value( "short_key", "tag" );
		TS_ASSERT( ss->get_comment( "short_key" ) == "" );
		TS_ASSERT( ss->has_energy( "short_key" ) );
		TS_ASSERT_DELTA( ss->get_energy( "short_key" ), 0.0, 1e-5 );
		// add checks for filling Pose here.
		ss = SilentStructOP( new ProteinSilentStruct( *start_pose, "tag", false ) );

		ss->add_comment( "comment", "tag" );
		TS_ASSERT_EQUALS( ss->get_comment( "comment" ), "tag" );
		ss->fill_pose( restored_pose, *rsd_set );
	}

	void test_string_from_pose() {
		using namespace core::pose;
		using namespace core::io::silent;
		Pose pose(*start_pose);
		add_comment( pose, "comment", "tag" );
		SilentStructOP ss( new ProteinSilentStruct );
		ss->fill_struct(pose);

		std::string retval;
		get_comment(pose,"comment",retval);
		TS_ASSERT_EQUALS( ss->get_comment("comment"), retval );
	}

	void test_save_and_restore_centroid() {
		using namespace core::io::silent;

		// configuration information for tests
		double const RMS_ERROR( 1e-3 );
		//double const CHI_ERROR( 1e-2 );
		double const BB_ERROR ( 1e-2 );
		std::string const silent_outfile( "test.silent_centroid.out" );
		pose::Pose restored_pose;

		utility::file::file_delete( silent_outfile );

		core::io::silent::SilentFileData sfd;

		core::io::silent::ProteinSilentStruct pss( *centroid_pose_, "tag", true );
		sfd.write_silent_struct( pss, silent_outfile );

		TS_ASSERT( !(centroid_pose_->is_fullatom()) );
		sfd.read_file( silent_outfile );
		core::io::silent::SilentFileData::iterator iter = sfd.begin();
		TS_ASSERT( iter->decoy_tag() == "tag" );
		iter->fill_pose( restored_pose );

		TS_ASSERT( !restored_pose.is_fullatom() )
			TS_ASSERT( centroid_pose_->total_residue() == restored_pose.total_residue() );
		for ( Size seqpos = 1; seqpos <= restored_pose.total_residue(); ++seqpos ) {
			TS_ASSERT_DELTA(
				centroid_pose_->phi( seqpos ), restored_pose.phi( seqpos ),
				BB_ERROR
			);
			TS_ASSERT_DELTA(
				centroid_pose_->psi( seqpos ), restored_pose.psi( seqpos ),
				BB_ERROR
			);
			TS_ASSERT_DELTA(
				centroid_pose_->omega( seqpos ), restored_pose.omega( seqpos ),
				BB_ERROR
			);
		}

		core::io::silent::ProteinSilentStruct original( *centroid_pose_, "start", false );
		original.fill_pose( *centroid_pose_ );

		Real rms_to_restored = scoring::CA_rmsd( *centroid_pose_, restored_pose );
		TR << "RMS error from save/restore: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < RMS_ERROR );

		utility::file::file_delete( silent_outfile );
	} // test_save_and_restore


}; // ProteinSilentTests
