// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic_thread.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueType.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/comparative_modeling/ThreadingMover.hh>
#include <protocols/comparative_modeling/StealSideChainsMover.hh>
#include <protocols/comparative_modeling/RecoverSideChainsMover.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <core/pack/optimizeH.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char* argv [] ) {
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::sequence;
	using namespace core::fragment;
	using namespace protocols::comparative_modeling;
	using namespace core::scoring;
	using namespace core::pack;
	using namespace core::io::pdb;
	using utility::vector1;
	using std::string;

	basic::Tracer tr( "basic_thread" );

	// options, random initialization
	devel::init( argc, argv );

	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	vector1< SequenceOP > seqs = core::sequence::read_fasta_file(
		option[ in::file::fasta ]()[1]
	);
	SWAligner aligner;

	// read in template pose
	core::pose::Pose template_pose;
	string const templ_fn( option[ in::file::template_pdb ]()[1] );
	string const templ_prefix( utility::file_basename( templ_fn ) );
	core::import_pose::pose_from_pdb( template_pose, *rsd_set, templ_fn );

	typedef vector1< SequenceOP >::const_iterator iter;
	for ( iter sequence = seqs.begin(), end = seqs.end(); sequence != end; ++sequence ) {
		// make an alignment from the template pose to the full-length fasta
		ScoringSchemeOP ss( new SimpleScoringScheme( 6, 1, -4, -1 ) );
		SequenceOP pdb_seq( new Sequence(
			template_pose.sequence(), "pdb", 1
		) );
		tr.Error << "sequence = " << *(*sequence) << std::endl;
		SequenceAlignment aln;
		if ( option[ james::thread_unaligned ]() ) {
			aln.add_sequence( *sequence );
			aln.add_sequence( pdb_seq );
		} else {
			aln = aligner.align( *sequence, pdb_seq, ss );
		}

		tr.Error << "rebuilding pose with alignment: " << std::endl;
		tr.Error << aln << std::endl;
		tr.flush();

		// build a model using threading (don't close loops)
		ThreadingMover threader( aln, template_pose );
		threader.build_loops( false );
		threader.repack_query( false );

		core::pose::Pose modeled_pose;
		core::pose::make_pose_from_sequence(
			modeled_pose, (*sequence)->sequence(), *rsd_set
		);
		tr.Error << "sequence of modeled pose " << modeled_pose.sequence() << std::endl;

		// steal aligned sidechains
		threader.apply( modeled_pose );
		SequenceMapping map = aln.sequence_mapping(1,2);
		StealSideChainsMover sc_mover( template_pose, map );
		sc_mover.apply( modeled_pose );

		// repack missing sidechains
		vector1< bool > residues_to_repack(
			modeled_pose.total_residue(), false
		);
		for ( core::Size ii = 1; ii <= (*sequence)->length(); ++ii ) {
			string const name3( modeled_pose.residue_type(ii).name3() );
			string const name3_src( template_pose.residue_type(map[ii]).name3());

			if ( map[ii] == 0 ) residues_to_repack[ii] = true;
			if ( name3 != name3_src ) residues_to_repack[ii] = true;
		}
		task::PackerTaskOP task
			= task::TaskFactory::create_packer_task( modeled_pose );
		task->initialize_from_command_line();
		task->restrict_to_repacking();
		task->restrict_to_residues(residues_to_repack);
		ScoreFunctionOP scorefxn(
			core::scoring::getScoreFunction()
		);
		pack_rotamers( modeled_pose, *scorefxn, task );

		// optimize hydrogens
		tr.Error << "setting up ideal hydrogen geometry on all residues."
			<< std::endl;
		//typedef core::conformation::ResidueOPs::iterator iter;
		//for ( iter it = modeled_pose.res_begin(), end = modeled_pose.res_end();
		//	it != end; ++it
		//) {
		//	core::conformation::idealize_hydrogens( *(*it), modeled_pose.conformation() );
		//}
		for ( core::Size ii = 1; ii <= modeled_pose.total_residue(); ++ii ) {
			core::conformation::ResidueOP iires = modeled_pose.residue(ii).clone();
			core::conformation::idealize_hydrogens( *iires, modeled_pose.conformation() );
			modeled_pose.replace_residue( ii, *iires, false );
		}

		core::id::AtomID_Mask missing( true );
		core::pose::initialize_atomid_map( missing, modeled_pose );
		core::pack::optimize_H_and_notify( modeled_pose, missing );

		// output PDB
		tr.Error << "sequence of modeled pose " << modeled_pose.sequence() << std::endl;
		string prefix( "thread" );
		string output_name( (*sequence)->id() + ".thread." + templ_prefix );
		if ( output_name.substr( output_name.length() - 4, 4 ) != ".pdb" ) {
			output_name += ".pdb";
		}
		tr.Error << "dumping pose to " << output_name << std::endl;
		modeled_pose.dump_pdb( output_name );
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
