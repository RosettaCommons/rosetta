// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Oliver Lange

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <protocols/jumping/JumpSetup.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <protocols/jumping/SheetBuilder.hh>
#include <protocols/jumping/SameStrand.hh>
#include <protocols/jumping/PairingLibrary.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <protocols/simple_filters/JumpEvaluator.hh>

// project headers
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>


#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragID_Iterator.hh>

#include <core/scoring/rms_util.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentFileData.hh>

// utility headers
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <utility/file/file_sys_util.hh>

//#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

// #include <core/fragment/FragData.hh>
// #include <core/fragment/FragCache.hh>
// #include <core/fragment/BBTorsionSRFD.hh>
// #include <protocols/abinitio/util.hh>


// #include <core/pack/rotamer_trials.hh>
// #include <core/pack/task/PackerTask.hh>
// #include <core/pack/task/TaskFactory.hh>


// #include <core/scoring/ScoreFunction.hh>

// #include <core/types.hh>

#include <basic/Tracer.hh>


#include <ObjexxFCL/format.hh>
#include <fstream>

//Auto Headers
#include <core/fragment/FragData.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/io/silent/EnergyNames.fwd.hh>
#include <utility/fix_boinc_read.hh>
#include <utility/vector1.hh>

// #include <numeric/random/random.hh>
// #include <numeric/numeric.functions.hh>

static basic::Tracer tr("protocol.jumping.cxxtest");
//MY_TRACERS("core.fragment.ConstantLengthFragments.cxxtest")

using namespace core;

using namespace ObjexxFCL;

class JumpingTest : public CxxTest::TestSuite {
public:
	JumpingTest() {}
	void setUp();
	void test_rt_library();
	void test_strand_fraction();
	void test_SheetBuilder();
	void test_jump_geometry();
	void test_save_and_restore_silentio_with_jumps();
private:
	fragment::ConstantLengthFragSet fragset3mer_;
	pose::Pose pose_;
};

using namespace protocols::jumping;
void
get_distance( pose::Pose &pose, JumpSetup jump_def, core::Real &d1, core::Real &d2 ) {
	JumpSetup::const_iterator it = jump_def.begin();
	Size const res1( it->jump_.start_ );
	Size const res2( it->jump_.end_ );
	chemical::ResidueType const& rt1 ( pose.residue_type ( res1 ) );
	chemical::ResidueType const& rt2 ( pose.residue_type ( res2 ) );
	PointPosition pO1;
	PointPosition pHN1;
	PointPosition pO2;
	PointPosition pHN2;

	// read CAs of 3 consecutive residues
	id::AtomID O1( rt1.atom_index ("O") , res1 );
	id::AtomID O2( rt2.atom_index ("O") , res2 );
	id::AtomID HN1( rt1.atom_index ("H") , res1 );
	id::AtomID HN2( rt2.atom_index ("H") , res2 );

	pO1 = pose.xyz( O1 );
	pO2 = pose.xyz( O2 );
	pHN1 = pose.xyz( HN1 );
	pHN2 = pose.xyz( HN2 );

	Vector r1 = pO1-pHN2;
	Vector r2 = pO2-pHN1;
	d1 = r1.length();
	d2 = r2.length();
}

void
apply_ss_jumps( pose::Pose &pose, JumpSetup jump_def, std::string tag ) {
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;
	using namespace fragment;
	using namespace protocols;
	using namespace jumping;

	JumpSample jump( jump_def );
	tr.Info << jump << std::endl;

	jump.set_fold_tree_in_pose( pose );
	kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_jump( true );

	using namespace format;
	FrameList jump_geometries;

	jump.steal_orientation_and_pleating( pose );
	jump.generate_jump_frags( *(jumping::StandardPairingLibrary::get_instance()), mm, true /* bWithTorsion */, jump_geometries );
	tr.Debug << *jump_geometries.front() << std::endl;
	tr.Info << "found " << jump_geometries.size() << " frames with "
		<<  jump_geometries[ 1 ]->nr_frags() << " frags in first frame " << std::endl;
	int ct = 1;
	std::string out_fn( "distance_"+tag+".dat" );
	std::string in_fn ("protocols/jumping/gold_distance_"+tag+".dat" );
	std::ofstream dist_out_file( out_fn.c_str() );
	std::ifstream dist_in_file( in_fn.c_str() );
	if ( !dist_in_file.good() ) {
		tr.Fatal << "can't find file " << in_fn.c_str() << std::endl;
	}
	bool success = true;
	for ( FragID_Iterator it=jump_geometries.begin(), eit=jump_geometries.end();
			it!=eit; ++it, ++ct ) {
		it->apply( mm, pose );
		std::ostringstream fn;
		tr.Trace << "apply frag_nr " << ct << std::endl;
		if ( ct<20 ) {
			fn << "sspair_" << tag << "_" << ct << ".pdb";
			if ( tr.Trace.visible() ) pose.dump_pdb( fn.str() );
		}
		Real d1,d2;
		get_distance( pose, jump_def, d1, d2 );
		tr.Trace << "distance between HN-O " << F(10,3,d1) << F(10,3,d2) <<std::endl;
		dist_out_file << F(10,3,d1) << F(10,3,d2) << std::endl;
		Real gold1, gold2;
		dist_in_file >> gold1 >> gold2;
		if ( success ) { //suppress output if failed, but run to the end to produce new 'gold' file
			TS_ASSERT_DELTA( d1, gold1, 0.002);
			TS_ASSERT_DELTA( d2, gold2, 0.002);
			if ( std::abs(d1-gold1) > 0.002 ) success = false;
			if ( std::abs(d2-gold2) > 0.002 ) success = false;
		}
	}
} //apply_ss_jumps


/// @detail apply jumps from SS-pair library to a small protein fragment ( two strands forming a short sheet )
/// the distances between O-HN for the jump residues are measured and compared to canned data
/// if the library changes the canned data has to be updated accordingly.
/// To this end copy files distance_xxx_dat (produced each time test runs) to gold_distance_xxx_dat
void
JumpingTest::test_rt_library()
{
	using namespace core;
	using namespace pose;
	using namespace protocols;
	using namespace jumping;


	Size nres = pose_.size();

	JumpSetup jump1_def( nres );
	jump1_def.add_jump( Interval( 4, 15), Interval( 8, 9) );
	apply_ss_jumps( pose_, jump1_def, "jump1" );

	JumpSetup jump2_def( nres );
	jump2_def.add_jump( Interval( 5, 14), Interval( 8, 8) );
	apply_ss_jumps( pose_, jump2_def, "jump2" );

}

void JumpingTest::test_save_and_restore_silentio_with_jumps()
{
	using namespace io::silent;
	using namespace core::chemical;
	//ResidueTypeSetCAP rsd_set
	// = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	//should live in core/io/silent but I'd like to use the JumpSetup stuff
	pose::Pose native_pose;
	//pose::Pose native_pose( create_test_in_pdb_pose() ); // jumping/test_in.pdb is different from all the other test_in.pdbs
	core::import_pose::pose_from_file( native_pose, "protocols/jumping/test_in.pdb" , core::import_pose::PDB_file); //has to be idealized or a decoy

	JumpSetup jump_def( native_pose.size() );
	jump_def.read_file( "protocols/jumping/jumps.def" );
	JumpSample jumps ( jump_def );
	jumps.set_fold_tree_in_pose( native_pose );

	//cleanup
	utility::file::file_delete( "test.out" );

	SilentFileOptions opts;
	SilentFileData sfd_out(opts);
	ProteinSilentStruct pss(opts);
	pss.fill_struct( native_pose, "native_structure" );
	sfd_out.write_silent_struct( pss, "test.out" );

	// if ( !utility::file::file_exists( "test_backward.out" ) )  {
	//  sfd_out.write_silent_struct( pss, "test_backward.out" );
	//  }

	// test internal compatibility
	{
		SilentFileData sfd(opts);
		sfd.read_file( "test.out" );

		for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
			pose::Pose pose;
			std::string tag = it->decoy_tag();
			it->fill_pose( pose );
			tr.Info << "RMSD between structures is " << scoring::CA_rmsd( pose, native_pose ) << std::endl;
			pose.dump_pdb( "silent_reread.pdb");
		}
	}

	// test backward compatibility
	{
		//SilentFileData sfd;
		//sfd.read_file( "test_backward.out" );
		/*
		for ( SilentFileData::const_iterator it=sfd.begin_const(), eit=sfd.end_const(); it!=eit; ++it ) {
		pose::Pose pose;
		std::string tag = it->decoy_tag();
		it->fill_pose( pose );
		tr.Info << "RMSD between structures is " << scoring::CA_rmsd( pose, native_pose ) << std::endl;
		pose.dump_pdb( "backward_silent_reread.pdb");
		}
		*/
	}
}


void JumpingTest::test_strand_fraction() {
	using namespace protocols::jumping;
	using namespace core::fragment;
	using namespace format;

	SecondaryStructure ss( fragset3mer_ );
	for ( Size i = 1; i<= ss.total_residue(); i++ ) {
		tr.Debug << "pos " << i
			<< " E " << RJ( 3, ss.strand_fraction( i ) )
			<< " L " << RJ( 3, ss.loop_fraction( i ) )
			<< " H " << RJ( 3, ss.helix_fraction( i ) )
			<< std::endl;
	}
	if ( tr.Debug.visible() ) ss.show( tr );
	// std::ofstream os( "gb3_secondary_structure.dat" );
	// ss.show( os );
	SecondaryStructure ss_gold;
	ss_gold.read_from_file( "protocols/jumping/gb3_secondary_structure.dat");
	// if ( tr.Debug.visible() )
	ss_gold.show( tr.Debug );
	TS_ASSERT_EQUALS( ss_gold.total_residue(), ss.total_residue() );
	for ( Size i = 1; i<= ss.total_residue(); i++ ) {
		if ( i > ss_gold.total_residue() ) break; // avoid RT errors in unit test
		TS_ASSERT_DELTA( ss_gold.helix_fraction( i ), ss.helix_fraction( i ), 0.001);
		TS_ASSERT_DELTA( ss_gold.loop_fraction( i ), ss.loop_fraction( i ), 0.001 );
		TS_ASSERT_DELTA( ss_gold.strand_fraction( i ), ss.strand_fraction( i ), 0.001 );
	}

	{
		SecondaryStructureOP ss( new SecondaryStructure( fragset3mer_ ) );
		SameStrand same_strand( ss );
		if ( tr.Debug.visible() ) {
			same_strand.show( tr.Debug );
		}
	}
}


void JumpingTest::test_SheetBuilder() {
	using namespace protocols::jumping;
	using namespace core::fragment;
	using namespace format;

	SecondaryStructureOP ss( new SecondaryStructure( fragset3mer_ ) );

	core::scoring::dssp::PairingsList pairings;
	read_pairing_list( "protocols/jumping/pairings.dat", pairings );

	{ //scoping
		SheetBuilder::SheetTopology sheets;
		sheets.push_back( 2 ); // make one sheet with 3 strands

		SheetBuilder sheet_jumps( ss, pairings, sheets );

		for ( int i = 1; i <= 10; i++ ) {
			JumpSample jumps = sheet_jumps.create_jump_sample();
			tr.Info << "SheetBuilder created jumps: 3 strand "
				<< jumps << std::endl;
		}
	}

	{ //scoping
		SheetBuilder::SheetTopology sheets;
		sheets.push_back( 3 ); // make one sheet with 3 strands

		SheetBuilder sheet_jumps( ss, pairings, sheets );

		for ( int i = 1; i <= 10; i++ ) {
			JumpSample jumps = sheet_jumps.create_jump_sample();
			tr.Info << "SheetBuilder created jumps: 4 strand "
				<< jumps << std::endl;
		}
	}

	{ //scoping
		SheetBuilder::SheetTopology sheets;
		sheets.push_back( 1 ); // make two sheets with 2 strands each
		sheets.push_back( 1 );

		SheetBuilder sheet_jumps( ss, pairings, sheets );

		for ( int i = 1; i <= 10; i++ ) {
			JumpSample jumps = sheet_jumps.create_jump_sample();
			tr.Info << "SheetBuilder created jumps: 2x2 strand "
				<< jumps << std::endl;
		}
	}
}

Real check_jump( pose::Pose& pose, pose::Pose const& native, Size jump_nr ) {
	protocols::simple_filters::JumpEvaluator eval( native, jump_nr );
	return eval.apply( pose );
}


void JumpingTest::test_jump_geometry() {
	using namespace core::fragment;
	using namespace protocols::jumping;

	pose::Pose native_pose;
	//pose::Pose native_pose( create_test_in_pdb_pose() ); // jumping/test_in.pdb is different from all the other test_in.pdbs
	core::import_pose::pose_from_file( native_pose, "protocols/jumping/test_in.pdb" , core::import_pose::PDB_file); //has to be idealized or a decoy

	JumpSetup jump_def( native_pose.size() );
	jump_def.read_file( "protocols/jumping/jumps.def" );
	JumpSample jumps ( jump_def );
	jumps.set_fold_tree_in_pose( native_pose );


	kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_jump( true );

	//get SS-fragments
	FrameList jump_frags;
	jumps.steal_orientation_and_pleating( native_pose );
	jumps.generate_jump_frags( *StandardPairingLibrary::get_instance(), mm, false /* bWithTorsion */, jump_frags );

	tr.Info << native_pose.fold_tree() << std::endl;
	//apply an SS-fragment for each jump-position ( = Frame )
	Size good_jumps[6] = { 3029 , 45 , 514, 678, 336, 263 };
	// to find these numbers
	// call application r_frag_quality
	// as r_frag_quality.macosgccdebug -in:file:native test_in.pdb -jumpss jumps.def
	//             -database ~/minirosetta_database/ | grep jump | sort -n -k 11 | head -n 20

	pose::Pose ref_pose( native_pose );
	Size jump_nr = 1;
	for ( FrameList::const_iterator it = jump_frags.begin(), eit = jump_frags.end(); it!=eit; ++it, ++jump_nr ) {
		tr.Debug  << "apply " << **it << std::endl;
		(*it )->fragment( good_jumps[ jump_nr - 1 ] ).apply( ref_pose, **it );
	}

	// --- if jump-geometries are not completely off numbers should be small say < 2
	for ( Size ii = 1; ii<=jumps.size();  ii++ ) {
		tr.Info << "Fragment quality " <<  " jump: " << ii <<" ";
		Real jump_RT_rms = check_jump( ref_pose, native_pose, ii ); //);
		tr.Info << jump_RT_rms << std::endl;
		TS_ASSERT_LESS_THAN( jump_RT_rms, 1.0 );
	}


	//copy pose
	pose::Pose pose( ref_pose );

	//trivial test --- should be fine and numbers should be 0
	for ( Size ii = 1; ii<=jumps.size();  ii++ ) {
		tr.Info << "START " <<  " jump: " << ii <<" ";
		Real jump_RT_rms = check_jump( pose, ref_pose, ii ); //);
		tr.Info << jump_RT_rms << std::endl;
		TS_ASSERT_DELTA( jump_RT_rms, 0.0, 0.00001 );
	}

	//big whack on structure
	for ( Size pos = 1; pos <= pose.size(); pos++ ) {
		///    if ( pos == 21 || pos == 20 || pos == 22 ) continue;
		pose.set_phi( 128, -45 );
		pose.set_psi( pos, -45 );
		pose.set_omega( pos, 180 );
	}

	//test jump geometry is still the same ---  numbers should be 0
	for ( Size ii = 1; ii<=jumps.size();  ii++ ) {
		tr.Info << "MOVED " <<  " jump: " << ii <<" ";
		Real jump_RT_rms = check_jump( pose, ref_pose, ii ); //);
		tr.Info << jump_RT_rms << std::endl;
		TS_ASSERT_DELTA( jump_RT_rms, 0.0, 0.00001 );
	}
	tr.Info << "FoldTree " << pose.fold_tree() << std::endl;
}

void JumpingTest::setUp() {
	core_init();
	fragset3mer_.read_fragment_file( "protocols/abinitio/mfr_aa2GB3_03_05.200_v1_3" );

	core::import_pose::pose_from_file( pose_, "protocols/jumping/jump_test.pdb" , core::import_pose::PDB_file);

}
