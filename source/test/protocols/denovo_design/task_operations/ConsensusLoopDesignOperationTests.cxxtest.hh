// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/denovo_design/ConsensusLoopDesignOperation.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::components::ConsensusLoopDesignOperation
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/task_operations/ConsensusLoopDesignOperation.hh>

// Protocol headers
#include <protocols/denovo_design/util.hh>

// Core headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/sequence/ABEGOManager.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Boost
#include <boost/assign.hpp>

// C++ headers

static THREAD_LOCAL basic::Tracer TR("protocols.denovo_design.ConsensusLoopDesignOperation.cxxtest");

// --------------- Test Class --------------- //
class ConsensusLoopDesignOperationTests : public CxxTest::TestSuite {
public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();

		// set preserve header always
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_loop_info_and_task_manipulation() {
		using namespace protocols::denovo_design::task_operations;
		core::pose::Pose input_pose;
		// loop is 17-18
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/disulf_test.pdb" );

		// create simple database
		ConsensusLoopDatabase & db = *ConsensusLoopDatabase::get_instance();

		ConsensusLoopDatabase::Abego const loop_abego = "BGBA";
		SurroundingSS ss_around( 'E', 'H' );

		AAFrequencies const aa_freq1 = db.frequencies( ss_around, loop_abego, 1 );

		/* Database for EH ss, BGBA abego
EH GB 1 100 A 0.04 0.0627617553498 -0.450466442605
EH GB 1 100 C 0.04 0.0127242399251 1.1453706243
EH GB 1 100 D 0.13 0.073877482684 0.565126368742
EH GB 1 100 E 0.02 0.0619159530468 -1.13004559582
EH GB 1 100 F 0.03 0.0368457278923 -0.205542297885
EH GB 1 100 G 0.06 0.105467330576 -0.564056679952
EH GB 1 100 H 0.03 0.0252153261362 0.173745392005
EH GB 1 100 I 0.05 0.0463777462221 0.0752032685265
EH GB 1 100 K 0.04 0.0611681249109 -0.424746764933
EH GB 1 100 L 0.06 0.0737757879789 -0.206686039389
EH GB 1 100 M 0.14 0.0170561863245 2.1051294504
EH GB 1 100 N 0.03 0.0539676437213 -0.587187294969
EH GB 1 100 P 0.14 0.060285117228 0.842557161479
EH GB 1 100 Q 0.02 0.0333856275617 -0.512393221295
EH GB 1 100 R 0.01 0.0488196593227 -1.58554799371
EH GB 1 100 S 0.03 0.0660457502155 -0.789150305306
EH GB 1 100 T 0.0 0.0585451456281 -10
EH GB 1 100 V 0.1 0.0582797472515 0.539915541468
EH GB 1 100 W 0.01 0.0112633086744 -0.118965329754
EH GB 1 100 Y 0.02 0.0322223393503 -0.476927707092
EH GB 2 100 A 0.04 0.0627617553498 -0.450466442605
EH GB 2 100 C 0.0 0.0127242399251 -10
EH GB 2 100 D 0.09 0.073877482684 0.197401588617
EH GB 2 100 E 0.01 0.0619159530468 -1.82319277638
EH GB 2 100 F 0.01 0.0368457278923 -1.30415458655
EH GB 2 100 G 0.4 0.105467330576 1.33306330493
EH GB 2 100 H 0.0 0.0252153261362 -10
EH GB 2 100 I 0.01 0.0463777462221 -1.53423464391
EH GB 2 100 K 0.05 0.0611681249109 -0.201603213619
EH GB 2 100 L 0.02 0.0737757879789 -1.30529832806
EH GB 2 100 M 0.0 0.0170561863245 -10
EH GB 2 100 N 0.14 0.0539676437213 0.953257745978
EH GB 2 100 P 0.0 0.060285117228 -10
EH GB 2 100 Q 0.01 0.0333856275617 -1.20554040186
EH GB 2 100 R 0.02 0.0488196593227 -0.892400813155
EH GB 2 100 S 0.04 0.0660457502155 -0.501468232854
EH GB 2 100 T 0.0 0.0585451456281 -10
EH GB 2 100 V 0.14 0.0582797472515 0.876387778089
EH GB 2 100 W 0.0 0.0112633086744 -10
EH GB 2 100 Y 0.02 0.0322223393503 -0.476927707092
EH GB 3 100 A 0.03 0.0627617553498 -0.738148515057
EH GB 3 100 C 0.02 0.0127242399251 0.452223443739
EH GB 3 100 D 0.08 0.073877482684 0.0796185529608
EH GB 3 100 E 0.01 0.0619159530468 -1.82319277638
EH GB 3 100 F 0.01 0.0368457278923 -1.30415458655
EH GB 3 100 G 0.06 0.105467330576 -0.564056679952
EH GB 3 100 H 0.06 0.0252153261362 0.866892572565
EH GB 3 100 I 0.01 0.0463777462221 -1.53423464391
EH GB 3 100 K 0.03 0.0611681249109 -0.712428837385
EH GB 3 100 L 0.14 0.0737757879789 0.640611820998
EH GB 3 100 M 0.01 0.0170561863245 -0.533927879211
EH GB 3 100 N 0.24 0.0539676437213 1.49225424671
EH GB 3 100 P 0.04 0.060285117228 -0.410205807016
EH GB 3 100 Q 0.02 0.0333856275617 -0.512393221295
EH GB 3 100 R 0.03 0.0488196593227 -0.486935705047
EH GB 3 100 S 0.05 0.0660457502155 -0.27832468154
EH GB 3 100 T 0.12 0.0585451456281 0.717693566029
EH GB 3 100 V 0.01 0.0582797472515 -1.76266955153
EH GB 3 100 W 0.0 0.0112633086744 -10
EH GB 3 100 Y 0.03 0.0322223393503 -0.0714625989839
EH GB 4 100 A 0.07 0.0627617553498 0.109149345331
EH GB 4 100 C 0.0 0.0127242399251 -10
EH GB 4 100 D 0.03 0.073877482684 -0.901210700051
EH GB 4 100 E 0.22 0.0619159530468 1.26784967698
EH GB 4 100 F 0.03 0.0368457278923 -0.205542297885
EH GB 4 100 G 0.01 0.105467330576 -2.35581614918
EH GB 4 100 H 0.01 0.0252153261362 -0.924866896663
EH GB 4 100 I 0.16 0.0463777462221 1.23835407833
EH GB 4 100 K 0.02 0.0611681249109 -1.11789394549
EH GB 4 100 L 0.02 0.0737757879789 -1.30529832806
EH GB 4 100 M 0.0 0.0170561863245 -10
EH GB 4 100 N 0.01 0.0539676437213 -1.68579958364
EH GB 4 100 P 0.25 0.060285117228 1.42237565673
EH GB 4 100 Q 0.03 0.0333856275617 -0.106928113187
EH GB 4 100 R 0.0 0.0488196593227 -10
EH GB 4 100 S 0.04 0.0660457502155 -0.501468232854
EH GB 4 100 T 0.02 0.0585451456281 -1.0740659032
EH GB 4 100 V 0.05 0.0582797472515 -0.153231639092
EH GB 4 100 W 0.02 0.0112633086744 0.574181850806
EH GB 4 100 Y 0.01 0.0322223393503 -1.17007488765
*/

		ConsensusLoopDesignOperation loopdesign;
		loopdesign.set_enrichment_threshold( 0.0 );

		utility::vector1< ConsensusLoopDesignOperation::AAs > const allowed_aas = boost::assign::list_of
			("CDHIMPV")("DGNV")("CDHLNT")("AEIPW");

		for ( core::Size loop_res=1; loop_res<=loop_abego.size(); ++loop_res ) {
			TR << "Loop res " << loop_res << std::endl;
			AAFrequencies const & aa_freqs = db.frequencies( ss_around, loop_abego, loop_res );
			ConsensusLoopDesignOperation::AAs const res_aas = loopdesign.forbidden_aas( aa_freqs );
			TS_ASSERT_EQUALS( 20 - res_aas.size(), allowed_aas[loop_res].size() );
			for ( std::string::const_iterator aa=allowed_aas[loop_res].begin(); aa!=allowed_aas[loop_res].end(); ++aa ) {
				TS_ASSERT_EQUALS( std::find( res_aas.begin(), res_aas.end(), *aa ), res_aas.end() );
			}
		}


		core::pack::task::TaskFactory tf;
		core::pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( input_pose );
		LoopInfo info;
		info.startres = 17;
		info.ss_around = ss_around;
		info.abego = loop_abego;
		loopdesign.disallow_aas( *task, info );

		TS_ASSERT_EQUALS( task->total_residue(), input_pose.total_residue() );
		for ( core::Size i=1; i<=input_pose.total_residue(); ++i ) {
			// all aas should be allowed except for residues 17, 18
			TS_ASSERT( task->being_packed( i ) );
			if ( i == 17 ) {
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), allowed_aas[2].size() );
				for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter t=task->residue_task( i ).allowed_residue_types_begin(), endt = task->residue_task( i ).allowed_residue_types_end(); t != endt; ++t ) {
					TS_ASSERT_DIFFERS( std::find( allowed_aas[2].begin(), allowed_aas[2].end(), (*t)->name1() ), allowed_aas[2].end() );
				}
			} else if ( i == 18 ) {
				// the 1 + is because HIS and HIS_D are two different types
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), 1 + allowed_aas[3].size() );
				for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter t=task->residue_task( i ).allowed_residue_types_begin(), endt = task->residue_task( i ).allowed_residue_types_end(); t != endt; ++t ) {
					TS_ASSERT_DIFFERS( std::find( allowed_aas[3].begin(), allowed_aas[3].end(), (*t)->name1() ), allowed_aas[3].end() );
				}
			} else if ( input_pose.residue( i ).name() == "CYS:disulfide" ) {
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), 1 );
			} else {
				// the +1 is for HIS_D
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), core::Size( core::chemical::num_canonical_aas )+1 );
			}
		}

		task = tf.create_task_and_apply_taskoperations( input_pose );
		loopdesign.set_include_adjacent_residues( true );
		loopdesign.disallow_aas( *task, info );
		TS_ASSERT_EQUALS( task->total_residue(), input_pose.total_residue() );
		for ( core::Size i=1; i<=input_pose.total_residue(); ++i ) {
			// all aas should be allowed except for residues 17, 18
			TS_ASSERT( task->being_packed( i ) );
			if ( ( i >= 16 ) && ( i <= 19 ) ) {
				core::Size const lres = i - 15;
				core::Size const n_his = std::count( allowed_aas[lres].begin(), allowed_aas[lres].end(), 'H' );
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), n_his + allowed_aas[lres].size() );
				for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter t=task->residue_task( i ).allowed_residue_types_begin(); t!=task->residue_task( i ).allowed_residue_types_end(); ++t ) {
					TS_ASSERT_DIFFERS( std::find( allowed_aas[lres].begin(), allowed_aas[lres].end(), (*t)->name1() ), allowed_aas[lres].end() );
				}
			} else if ( input_pose.residue( i ).name() == "CYS:disulfide" ) {
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), 1 );
			} else {
				// the +1 is for HIS_D
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), core::Size( core::chemical::num_canonical_aas )+1 );
			}
		}

		// test getting loop info
		LoopInfoVec loopinfo = loopdesign.get_loop_info( input_pose );
		TS_ASSERT_EQUALS( loopinfo.size(), 4 );

		// get ss/abego for comparision
		core::scoring::dssp::Dssp dssp( input_pose );
		std::string const ss = dssp.get_dssp_secstruct();
		utility::vector1< std::string > const abego = core::sequence::get_abego( input_pose, 1 );

		for ( LoopInfoVec::const_iterator l=loopinfo.begin(); l!=loopinfo.end(); ++l ) {
			TR << "Start: " << l->startres << " abego: " << l->abego << " ss: " << ss << std::endl;
			TS_ASSERT_EQUALS( ss[ l->startres - 2 ], l->ss_around.before );
			TS_ASSERT_EQUALS( ss[ l->startres + l->abego.size() - 1 ], l->ss_around.after );
			for ( core::Size res=l->startres; res<=l->startres+l->abego.size()-3; ++res ) {
				core::Size const abegoidx = res - l->startres + 1;
				TR << "Res: " << res << " abegoidx: " << abegoidx << std::endl;
				TS_ASSERT_EQUALS( ss[ res - 1 ], 'L' );
				TS_ASSERT_EQUALS( abego[ res ].size(), 1 );
				TS_ASSERT_EQUALS( abego[ res ][ 0 ], l->abego[ abegoidx ] );
			}
		}

	}

	void test_read_db() {
		using namespace protocols::denovo_design::task_operations;

		ConsensusLoopDatabase & db = *ConsensusLoopDatabase::get_instance();

		// just look at alanine frequency data
		utility::vector1< core::Real > const ala_freqs = boost::assign::list_of
			(0.152334152334)(0.012285012285)(0.044226044226)(0.0405405405405);

		utility::vector1< core::Real > const ala_enrichment = boost::assign::list_of
			(0.886730581866)(-1.63096589075)(-0.350032045283)(-0.437043422273);

		// get stats from a known motif -- GB connecting two helices
		for ( core::Size lres=1; lres<=4; ++lres ) {
			AAFrequencies const aa_freqs = db.frequencies( SurroundingSS( 'H', 'H' ), "AGBA", lres );
			TS_ASSERT_DELTA( aa_freqs.frequency('A').frequency(), ala_freqs[lres], 1e-6 );
			TS_ASSERT_DELTA( aa_freqs.frequency('A').enrichment(), ala_enrichment[lres], 1e-6 );
		}
	}


};
