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
/// @author TJ Brunette


#include <core/types.hh>
#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <core/chemical/util.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/SWAligner.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>


#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <map>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <utility/io/ozstream.hh>

std::map< std::string, core::pose::Pose >
poses_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
	using std::map;
	using std::string;
	using core::pose::Pose;
	using utility::vector1;
	using utility::file::file_exists;
	using core::import_pose::pose_from_file;

	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
	map< string, Pose > poses;

	typedef vector1< string >::const_iterator iter;
	for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
		if ( file_exists(*it) ) {
			Pose pose;
			core::import_pose::pose_from_file( pose, *rsd_set, *it , core::import_pose::PDB_file);
			string name = utility::file_basename( *it );
			name = name.substr( 0, 5 );
			poses[name] = pose;
		}
	}

	return poses;
}

int
main( int argc, char * argv [] ) {
	try {

		devel::init( argc, argv );

		using std::map;
		using std::string;
		using core::Real;
		using core::Size;
		using core::pose::Pose;
		using utility::vector1;
		using core::sequence::SequenceAlignment;
		using core::sequence::SequenceProfile;
		using core::id::SequenceMapping;
		using core::import_pose::pose_from_file;
		using namespace core::chemical;
		using namespace core::sequence;
		using namespace core::io::silent;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		basic::Tracer tr( "fix_alignment_to_match_pdb" );

		// read pdbs into poses
		map< string, Pose > poses;
		poses = poses_from_cmd_line(
			option[ in::file::template_pdb ]()
		);

		// Set up alignment options
		std::string output_fn = option[ out::file::alignment ]();
		utility::io::ozstream output( output_fn );
		SWAligner sw_align;
		ScoringSchemeOP ss( new SimpleScoringScheme( 120, 0, -140, 0 ) );
		Size query_index_ = 1 ;
		Size template_index_ = 2;

		// read in alignments, fix and output each one
		vector1< std::string > align_fns = option[ in::file::alignment ]();
		typedef vector1< string >::const_iterator aln_iter;
		for ( aln_iter aln_fn = align_fns.begin(), aln_end = align_fns.end();
				aln_fn != aln_end; ++aln_fn
				) {
			vector1< SequenceAlignment > alns = core::sequence::read_aln(
				option[ cm::aln_format ](), *aln_fn
			);

			for ( vector1< SequenceAlignment >::iterator it = alns.begin(),
					end = alns.end();
					it != end; ++it
					) {
				string const template_id( it->sequence(2)->id().substr(0,5) );
				string const template_id_full( it->sequence(2)->id());
				tr << *it << std::endl;
				tr << "template " << it->sequence(2)->id() << " => " << template_id << std::endl;
				string target_ungapped(it->sequence(1)->ungapped_sequence());
				string template_ungapped( it->sequence(2)->ungapped_sequence() );
				map< string, Pose >::iterator pose_it = poses.find( template_id );
				if ( pose_it == poses.end() ) {
					string msg( "Error: can't find pose (id = " + template_id + ")");
					//utility_exit_with_message(msg);
					tr.Error << msg << std::endl;
					continue;
				} else {
					Pose template_pose;
					template_pose = pose_it->second;
					string pdbTemplate_ungapped(template_pose.sequence());

					SequenceOP t_target_seq( new Sequence(
						target_ungapped,
						it->sequence(1)->id(),
						it->sequence(1)->start()
						) );

					SequenceOP t_align_seq( new Sequence(
						template_ungapped,
						template_id + "_align_seq",
						it->sequence(2)->start()
						) );

					SequenceOP t_pdb_seq( new Sequence (
						pdbTemplate_ungapped,
						template_id + "_pdb_seq",
						1
						) );
					SequenceOP t_pdb_seq_out( new Sequence (
						pdbTemplate_ungapped,
						template_id_full,
						1
						) );
					//SequenceAlignment intermediate = sw_align.align( t_align_seq, t_pdb_seq, ss );
					SequenceAlignment intermediate = sw_align.align( t_align_seq, t_pdb_seq, ss );
					if ( intermediate.identities() != intermediate.length() ) {
						tr.Warning << "Mismatch between sequence from alignment! Picked up error ";
						tr.Warning << "alignment: " << std::endl << intermediate << std::endl;
					}

					SequenceMapping query_to_fullseq = it->sequence_mapping(1,2);
					tr.Debug << "Query:    " << it->sequence( query_index_ ) << std::endl;
					tr.Debug << "Template: " << it->sequence( template_index_ ) << std::endl;
					tr.Debug << "Original Mapping:" << query_index_ << "-->" << template_index_
						<< std::endl;

					tr.Debug << "query to fullseq" << std::endl;
					query_to_fullseq.show( tr.Debug );

					SequenceMapping intermed_map = intermediate.sequence_mapping(1,2);

					tr.Debug << "intermediate map" << std::endl;
					intermed_map.show(tr.Debug);
					SequenceMapping query_to_pdbseq = core::sequence::transitive_map(
						query_to_fullseq, intermed_map
					);

					tr.Debug << "Transitive Map" << std::endl;
					query_to_pdbseq.show( tr.Debug );

					SequenceAlignment new_aln = mapping_to_alignment(query_to_pdbseq,t_target_seq,t_pdb_seq_out);

					new_aln.scores( it->scores() );
					std::map< std::string, core::Real > scores( new_aln.scores() );
					new_aln.printGrishinFormat(output);
				} // if found a template pdb
			} // for alns
		} // for aln_fn
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
} // main
