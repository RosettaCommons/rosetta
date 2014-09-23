// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file full_length_model.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <core/chemical/util.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/PDBSilentStruct.hh>

#include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>

#include <protocols/relax/MiniRelax.cc>

#include <protocols/loops/loops_main.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/ThreadingMover.hh>
#include <protocols/comparative_modeling/StealSideChainsMover.hh>
#include <protocols/comparative_modeling/RecoverSideChainsMover.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>


void restore_hack( core::pose::Pose & pose ) {
	core::import_pose::PDBSilentStruct ss;

	ss.fill_struct( pose );
	ss.fill_pose( pose, *(core::chemical::rsd_set_from_cmd_line().lock()) );
}

//class FullLengthPose {
//public:
//
//void apply() {
//		// make an alignment from the input pose to the full-length fasta
//		ScoringSchemeOP ss( new SimpleScoringScheme( 6, 1, -4, -1 ) );
//		SequenceOP full_length( new Sequence(
//			sequence, "full_length", 1
//		) );
//		SequenceOP pdb_seq( new Sequence(
//			start_pose.sequence(), "pdb", 1
//		) );
//		SequenceAlignment aln = sw_align.align( full_length, pdb_seq, ss );
//		tr.Debug << "rebuilding pose with alignment: " << std::endl;
//		tr.Debug << aln << std::endl;
//		tr.flush();
//
//		// build a model using threading (close loops)
//		ThreadingMover threader( aln, start_pose );
//		threader.build_loops( true );
//		if ( option[ cm::min_loop_size ]() ) {
//			threader.min_loop_size( option[ cm::min_loop_size ]() );
//		}
//		core::pose::Pose full_length_pose;
//		core::pose::make_pose_from_sequence( full_length_pose, sequence, *rsd_set );
//		threader.frag_libs( frag_libs );
//
//		// repack rebuilt sidechains
//		RecoverSideChainsMover mover( threader );
//		mover.apply( full_length_pose );
//		vector1< bool > residues_to_repack(
//			full_length_pose.total_residue(), false
//		);
//		MoveMapOP mm( new MoveMap() );
//		mm->set_bb (false);
//		mm->set_chi(false);
//		SequenceMapping map = aln.sequence_mapping(1,2);
//		for ( core::Size ii = 1; ii <= full_length->length(); ++ii ) {
//			if ( map[ii] == 0 ) {
//				residues_to_repack[ii] = true;
//				mm->set_bb ( ii, true );
//				mm->set_chi( ii, true );
//			}
//		}
//		StealSideChainsMover sc_mover( start_pose, map );
//		sc_mover.apply( full_length_pose );
//
//		task::PackerTaskOP task
//			= task::TaskFactory::create_packer_task( full_length_pose );
//		task->initialize_from_command_line();
//		task->restrict_to_repacking();
//		task->restrict_to_residues(residues_to_repack);
//		ScoreFunctionOP scorefxn(
//			get_score_function()
//		);
//		restore_hack( full_length_pose );
//
//};

int
main( int argc, char* argv [] ) {
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::sequence;
	using namespace core::fragment;
	using namespace core::import_pose::pose_stream;
	using namespace protocols::comparative_modeling;
	using namespace core::scoring;
	using namespace core::pack;
	using namespace core::kinematics;
	using core::Size;
	using utility::vector1;

	basic::Tracer tr( "full_length_model" );

	// options, random initialization
	devel::init( argc, argv );

	ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
	MetaPoseInputStream input = streams_from_cmd_line();
	std::string sequence = core::sequence::read_fasta_file(
		option[ in::file::fasta ]()[1]
	)[1]->sequence();
	NWAligner sw_align;

	// initialize fragments
	utility::vector1< core::fragment::FragSetOP > frag_libs;
	protocols::loops::read_loop_fragments( frag_libs );

	while ( input.has_another_pose() ) {
		core::pose::Pose start_pose;
		input.fill_pose( start_pose, *rsd_set );

		// make an alignment from the input pose to the full-length fasta
		ScoringSchemeOP ss( new SimpleScoringScheme( 6, 1, -4, -1 ) );
		SequenceOP full_length( new Sequence(
			sequence, "full_length", 1
		) );
		SequenceOP pdb_seq( new Sequence(
			start_pose.sequence(), "pdb", 1
		) );
		SequenceAlignment aln = sw_align.align( full_length, pdb_seq, ss );
		tr.Debug << "rebuilding pose with alignment: " << std::endl;
		tr.Debug << aln << std::endl;
		tr.flush();

		// build a model using threading (close loops)
		ThreadingMover threader( aln, start_pose );
		threader.build_loops( true );
		if ( option[ cm::min_loop_size ]() ) {
			threader.min_loop_size( option[ cm::min_loop_size ]() );
		}
		core::pose::Pose full_length_pose;
		core::pose::make_pose_from_sequence( full_length_pose, sequence, *rsd_set );
		threader.frag_libs( frag_libs );

		// repack rebuilt sidechains
		RecoverSideChainsMover mover( threader );
		mover.apply( full_length_pose );
		vector1< bool > residues_to_repack(
			full_length_pose.total_residue(), false
		);
		MoveMapOP mm( new MoveMap() );
		mm->set_bb (false);
		mm->set_chi(false);
		core::id::SequenceMapping map = aln.sequence_mapping(1,2);
		for ( core::Size ii = 1; ii <= full_length->length(); ++ii ) {
			if ( map[ii] == 0 ) {
				residues_to_repack[ii] = true;
				mm->set_bb ( ii, true );
				mm->set_chi( ii, true );
			}
		}
		StealSideChainsMover sc_mover( start_pose, map );
		sc_mover.apply( full_length_pose );

		task::PackerTaskOP task
			= task::TaskFactory::create_packer_task( full_length_pose );
		task->initialize_from_command_line();
		task->restrict_to_repacking();
		task->restrict_to_residues(residues_to_repack);
		ScoreFunctionOP scorefxn(
			get_score_function()
		);
		restore_hack( full_length_pose );
		//core::import_pose::pose_from_pdb( pose, *rsd_set, "debug.pdb" );

		//(*scorefxn)(full_length_pose);
		//pack_rotamers( full_length_pose, *scorefxn, task );
		//full_length_pose.dump_pdb("debug.pdb");
		alignment_into_pose( aln, full_length_pose );

		if ( option[ relax::mini ]() ) {
			protocols::relax::MiniRelax relax( scorefxn );
			mm->show( std::cout, full_length_pose.total_residue() );
			//relax.set_movemap(mm);
			relax.apply(full_length_pose);
		}

		std::string output_prefix( core::pose::tag_from_pose(start_pose) );
		// output PDB
		utility::io::ozstream output( output_prefix + "_full_length.pdb" );
		output << "REMARK REBUILT_RESIDUES";
		for ( Size ii = 1; ii <= residues_to_repack.size(); ++ii ) {
			if ( residues_to_repack[ii] ) {
				output << " " << ii;
			}
		}
		output << std::endl;
		output << "REMARK query_aln    " << aln.sequence(1)->to_string() << std::endl;
		output << "REMARK template_aln " << aln.sequence(2)->to_string()<< std::endl;
		output << std::endl;
		core::io::pdb::dump_pdb( full_length_pose, output );
		output.close();
		//full_length_pose.dump_pdb( output_prefix + "_full_length.pdb" );
	} // has_another_pose()

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

} // int main( int argc, char * argv [] )
