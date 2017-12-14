// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file score_aln2.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <core/chemical/util.hh>

#include <basic/options/option.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceFactory.hh>
#include <core/sequence/CompositeSequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/CompositeScoringScheme.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>

#include <protocols/comparative_modeling/Align_RmsdEvaluator.hh>
#include <protocols/comparative_modeling/PartialThreadingMover.hh>

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
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>


std::map< std::string, core::pose::Pose >
poses_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
	using std::map;
	using std::string;
	using core::pose::Pose;
	using utility::file::file_exists;
	using utility::vector1;
	using core::import_pose::pose_from_file;
	using namespace core::chemical;

	ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
	map< string, Pose > poses;

	for ( auto const & it : fn_list ) {
		if ( file_exists(it) ) {
			Pose pose;
			core::import_pose::pose_from_file( pose, *rsd_set, it , core::import_pose::PDB_file);
			string name = utility::file_basename( it );
			name = name.substr( 0, 5 );
			poses[name] = pose;
		}
	}

	return poses;
}

void print_seq_map(
	std::ostream & out,
	std::map< std::string, core::sequence::CompositeSequenceOP > const & seqs
) {
	using std::map;
	using std::string;
	using core::sequence::CompositeSequenceOP;
	for ( auto const & seq : seqs ) {
		out << seq.first << " => " << *seq.second << std::endl;
	}
}

std::map< std::string, core::sequence::CompositeSequenceOP >
composite_sequences_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
	using core::Size;
	using std::map;
	using std::string;
	using utility::file::file_exists;
	using namespace core::sequence;

	typedef map< string, CompositeSequenceOP > container;
	map< string, CompositeSequenceOP > seqs;
	runtime_assert( (fn_list.size() % 2) == 0 );

	for ( Size ii = 1; ii <= fn_list.size() - 1; ii += 2 ) {
		//std::cout << "making sequence from " << fn_list[ii] << " and " << fn_list[ii+1]
		// << std::endl;
		SequenceOP seq = SequenceFactory::get_instance()->seq_from_file(
			fn_list[ii], fn_list[ii+1]
		);
		//std::cout << "read sequence " << seq->to_string() << std::endl;
		string name = utility::file_basename(fn_list[ii]).substr(0,5);
		container::const_iterator iter = seqs.find(name);
		if ( iter == seqs.end() ) {
			//std::cout << "making new composite sequence for " << name << std::endl;
			seqs[name] = core::sequence::CompositeSequenceOP( new CompositeSequence );
		}
		seqs[name]->add_sequence(seq);
	}
	//print_seq_map(std::cout,seqs);

	return seqs;
}

void save_per_residue_scores(
	std::string const & fn,
	core::sequence::SequenceAlignment const & aln,
	core::sequence::ScoringSchemeOP ss,
	std::string const & aln_id
) {
	using core::Size;
	using core::Real;
	using utility::vector1;

	vector1< Real > scores( aln.calculate_per_position_scores(ss) );
	utility::io::ozstream output(fn);
	output << "resi score aln_id" << std::endl;
	for ( Size ii = 1; ii <= aln.length(); ++ii ) {
		if ( !aln.sequence(1)->is_gap(ii) ) {
			Size const resi ( aln.sequence(1)->resnum(ii) );
			Real const score( scores[ii] );
			output << resi << ' ' << score << ' ' << aln_id << std::endl;
		}
	}
	output.close();
}

int
main( int argc, char* argv [] ) {
	try {

		// options, random initialization
		devel::init( argc, argv );

		using std::map;
		using std::string;
		using core::Real;
		using core::Size;
		using core::pose::Pose;
		using core::pose::PoseOP;
		using utility::vector1;
		using core::import_pose::pose_from_file;
		using protocols::comparative_modeling::Align_RmsdEvaluator;
		using core::pose::make_pose_from_sequence;
		using protocols::comparative_modeling::PartialThreadingMover;

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using namespace core::io::silent;
		using namespace core::sequence;

		basic::Tracer tr( "score_aln" );
		//std::cout << "making sequence from " << option[in::file::seq]()[1] << " and "
		// << option[ in::file::seq ]()[2] << std::endl;

		SequenceOP qseq = SequenceFactory::get_instance()->seq_from_file(
			option[ in::file::seq ]()[1], option[ in::file::seq ]()[2]
		);
		CompositeSequenceOP query_seq( new CompositeSequence );
		query_seq->add_sequence(qseq);
		map< string, CompositeSequenceOP > seqs = composite_sequences_from_cmd_line(
			option[ in::file::seq ]()
		);

		vector1< std::string > align_fns = option[ in::file::alignment ]();

		PoseOP native_pose( new Pose );
		bool have_native( false );
		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_file(
				*native_pose, *(rsd_set_from_cmd_line().lock()), option[ in::file::native ](), core::import_pose::PDB_file
			);
			have_native = true;
		}

		map< string, Pose > poses;
		if ( have_native ) {
			poses = poses_from_cmd_line(
				option[ in::file::template_pdb ]()
			);
		}

		// scoring scheme for aligning profiles
		//std::string const scoring_scheme_type( option[ cm::seq_score ]()[1] );

		CompositeScoringSchemeOP ss( new CompositeScoringScheme );
		ScoringSchemeFactory ssf;
		vector1< string > ss_names( option[ cm::seq_score ]() );
		for ( Size ii = 1; ii <= ss_names.size(); ++ii ) {
			ss->add_scoring_scheme( ssf.get_scoring_scheme( ss_names[ii] ) );
		}

		// for the ProfSim scoring scheme, the optimal opening and extension
		// penalties were 2 and 0.2, with a scoring shift of -0.45 applied to
		// all ungapped aligned pairs for database searches.
		Real const gap_open  ( option[ cm::min_gap_open ]() );
		Real const gap_extend( option[ cm::min_gap_extend ]() );
		ss->gap_open  ( gap_open   );
		ss->gap_extend( gap_extend );

		SilentFileOptions opts; // initialized from the command line
		SilentFileData sfd(opts);

		using aln_iter = vector1<string>::const_iterator;
		for ( aln_iter aln_fn = align_fns.begin(), aln_end = align_fns.end();
				aln_fn != aln_end; ++aln_fn
				) {
			vector1< SequenceAlignment > alns = core::sequence::read_aln(
				option[ cm::aln_format ](), *aln_fn
			);

			for ( auto & aln : alns ) {
				string const template_id( aln.sequence(2)->id().substr(0,5) );
				tr << aln << std::endl;
				tr << "id " << aln.sequence(2)->id() << " => " << template_id
					<< std::endl;
				string const ungapped_query( aln.sequence(1)->ungapped_sequence() );

				SilentStructOP ss_out( new ScoreFileSilentStruct(opts) );
				map< string, CompositeSequenceOP >::iterator seq_it = seqs.find( template_id );
				if ( seq_it == seqs.end() ) {
					//print_seq_map( std::cerr, seqs );
					string msg( "Error: can't find seq (id = " + template_id + ")" );
					//utility_exit_with_message(msg);
					tr.Error << msg << std::endl;
				} else {
					SequenceOP template_seq = seq_it->second;

					utility::vector1< SequenceOP > my_seqs;
					my_seqs.push_back( query_seq->clone() );
					my_seqs.push_back( template_seq->clone() );
					//std::cout << "query sequence " << std::endl;
					//std::cout << query_seq->to_string() << std::endl;
					//std::cout << std::endl;
					//std::cout << "template sequence " << std::endl;
					//std::cout << template_seq->to_string() << std::endl;
					//std::cout << "template sequence length " << template_seq->sequence().length() << std::endl;
					//std::cout << std::endl;

					SequenceAlignment rescore_aln = steal_alignment( aln, my_seqs );
					std::cout << "stolen alignment:" << std::endl;
					std::cout << rescore_aln << std::endl;
					rescore_aln.score( rescore_aln.calculate_score_sum_of_pairs( ss ) );
					//string const ungapped_templ( it->sequence(2)->ungapped_sequence() );
					if ( option[ run::debug ]() ) {
						string const id_out( aln.sequence(2)->id() );
						save_per_residue_scores( id_out + ".dat", rescore_aln, ss, id_out );
					}

					ss_out->add_energy( "aln_score", rescore_aln.score() );
					//ss_out->add_energy( "n_ali_templ", ungapped_templ.length() );
					tr.Debug << "score(" << ss_out->decoy_tag() << ") = " << rescore_aln.score()
						<< std::endl;
				}

				if ( have_native ) {
					// calc rmsd/gdt stats
					map< string, Pose >::iterator pose_it = poses.find( template_id );
					if ( pose_it == poses.end() ) {
						string msg( "Error: can't find pose (id = "
							+ template_id + ")"
						);
						//utility_exit_with_message(msg);
						tr.Error << msg << std::endl;
					} else {
						Pose query_pose, template_pose;
						make_pose_from_sequence(
							query_pose,
							ungapped_query,
							*(rsd_set_from_cmd_line().lock())
						);
						template_pose = pose_it->second;
						PartialThreadingMover mover(aln,template_pose);
						Align_RmsdEvaluator eval( native_pose, "ali", true );
						eval.apply(query_pose,"ali",*ss_out);
						if ( option[ run::debug ]() ) {
							string const id_out( aln.sequence(2)->id() );
							query_pose.dump_pdb(id_out + ".pdb");
						}
					} // template pdb check
				} // have native
				ss_out->scoreline_prefix( "" );
				ss_out->decoy_tag( aln.sequence(2)->id() );
				ss_out->add_energy( "n_ali_query", ungapped_query.length() );
				ss_out->add_string_value( "template", template_id );
				sfd.write_silent_struct( *ss_out, option[ out::file::silent ]() );
			} // alns
		} // for ( it in aligns )

		tr.Debug << "finished rescoring alignments." << std::endl;
		tr.flush_all_channels();
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
} // int main( int argc, char * argv [] )
