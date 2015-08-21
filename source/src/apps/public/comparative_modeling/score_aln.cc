// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file score_aln.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/rms_util.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>

#include <protocols/comparative_modeling/PartialThreadingMover.hh>
#include <protocols/comparative_modeling/coord_util.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/model_quality/rms.hh>

#include <map>
#include <iostream>
#include <string>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/io/ozstream.hh>
#include <ObjexxFCL/format.hh>

///////////////////////////////////////////////////////////////////////////////

using std::map;
using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using namespace basic;
using namespace core::sequence;
using namespace basic::options;
using namespace basic::options::OptionKeys;

std::map< std::string, core::pose::Pose >
poses_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
	using std::map;
	using std::string;
	using core::pose::Pose;
	using utility::file::file_exists;
	using core::import_pose::pose_from_pdb;
	using namespace core::chemical;

	ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
	map< string, Pose > poses;

	typedef vector1< string >::const_iterator iter;
	for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
		if ( file_exists(*it) ) {
			Pose pose;
			core::import_pose::pose_from_pdb( pose, *rsd_set, *it );
			string name = utility::file_basename( *it );
			name = name.substr( 0, 5 );
			poses[name] = pose;
		}
	}

	return poses;
}

void print_seq_map(
	std::ostream & out,
	std::map< std::string, core::sequence::SequenceOP > const & seqs
) {
	using std::map;
	using std::string;
	using core::sequence::SequenceOP;
	typedef map< string, SequenceOP >::const_iterator iter;
	for ( iter it = seqs.begin(), end = seqs.end(); it != end; ++it ) {
		out << it->first << " => " << *it->second << std::endl;
	}
}

std::map< std::string, core::sequence::SequenceOP >
sequences_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
	using std::map;
	using std::string;
	using utility::file::file_exists;
	using namespace core::sequence;

	map< string, SequenceOP > seqs;

	typedef vector1< string >::const_iterator iter;
	for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
		SequenceProfileOP prof( new SequenceProfile );
		if ( file_exists(*it) ) {
			prof->read_from_file( *it );
			prof->convert_profile_to_probs( 1.0 ); // was previously implicit in read_from_file()
			string name = utility::file_basename( *it );
			name = name.substr( 0, 5 );
			seqs[name] = prof;
		}
	}

	return seqs;
}

void save_per_residue_scores(
	std::string const & fn,
	core::sequence::SequenceAlignment const & aln,
	core::sequence::ScoringSchemeOP ss,
	std::string const & aln_id
) {
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
		using utility::vector1;
		using core::sequence::SequenceAlignment;
		using core::sequence::SequenceProfile;
		using core::import_pose::pose_from_pdb;
		using namespace core::chemical;
		using namespace core::io::silent;

		basic::Tracer tr( "score_aln" );

		SequenceProfileOP query_prof( new SequenceProfile );
		map< string, SequenceOP > seqs;
		if ( option[ in::file::pssm ].user() ) {
			query_prof->read_from_file( option[ in::file::pssm ]()[1] );
			query_prof->convert_profile_to_probs( 1.0 ); // was previously implicit in read_from_file()
			seqs = sequences_from_cmd_line(
				option[ in::file::pssm ]()
			);
		}

		vector1< std::string > align_fns = option[ in::file::alignment ]();
		//std::string fasta_seq = core::sequence::read_fasta_file_return_str( option[ in::file::fasta ](1) );
		std::string fasta_seq = read_fasta_file_str(option[ in::file::fasta ](1))[1] ;
		Pose native_pose;
		bool have_native( false );
		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_pdb(
				native_pose,
				*(rsd_set_from_cmd_line().lock()),
				option[ in::file::native ]()
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
		std::string const scoring_scheme_type( option[ cm::seq_score ]()[1] );
		ScoringSchemeFactory ssf;
		ScoringSchemeOP ss( ssf.get_scoring_scheme( scoring_scheme_type ) );

		// for the ProfSim scoring scheme, the optimal opening and extension
		// penalties were 2 and 0.2, with a scoring shift of -0.45 applied to
		// all ungapped aligned pairs for database searches.
		Real const gap_open  ( option[ cm::min_gap_open ]() );
		Real const gap_extend( option[ cm::min_gap_extend ]() );
		ss->gap_open  ( gap_open   );
		ss->gap_extend( gap_extend );

		SilentFileData sfd;

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
				tr << *it << std::endl;
				tr << "id " << it->sequence(2)->id() << " => " << template_id
					<< std::endl;
				string const ungapped_query( it->sequence(1)->ungapped_sequence() );

				SilentStructOP ss_out( new ScoreFileSilentStruct );
				//bool encountered_error(false);
				map< string, SequenceOP >::iterator seq_it = seqs.find( template_id );
				if ( seq_it == seqs.end() ) {
					//print_seq_map( std::cerr, seqs );
					string msg( "Error: can't find seq (id = " + template_id + ")" );
					//utility_exit_with_message(msg);
					tr.Error << msg << std::endl;
					//continue;
				} else {
					SequenceOP template_prof = seq_it->second;

					utility::vector1< SequenceOP > my_seqs;
					my_seqs.push_back( query_prof->clone() );
					my_seqs.push_back( template_prof->clone() );
					SequenceAlignment rescore_aln = steal_alignment( *it, my_seqs );
					Real const aln_score( rescore_aln.calculate_score_sum_of_pairs( ss ) );
					rescore_aln.score( aln_score );
					//string const ungapped_templ( it->sequence(2)->ungapped_sequence() );
					save_per_residue_scores(
						it->sequence(2)->id() + ".dat", rescore_aln, ss,
						it->sequence(2)->id()
					);

					ss_out->add_energy( "aln_score", aln_score );
					//ss_out->add_energy( "n_ali_templ", ungapped_templ.length() );
					tr.Debug << "score(" << ss_out->decoy_tag() << ") = " << aln_score
						<< std::endl;
				}

				if ( have_native ) {
					// calc rmsd/gdt stats
					Pose query_pose, template_pose;
					core::pose::make_pose_from_sequence(query_pose, fasta_seq, *(rsd_set_from_cmd_line().lock()));

					map< string, Pose >::iterator pose_it = poses.find( template_id );
					if ( pose_it == poses.end() ) {
						string msg( "Error: can't find pose (id = "
							+ template_id + ")"
						);
						//utility_exit_with_message(msg);
						tr.Error << msg << std::endl;
						continue;
					} else {
						template_pose = pose_it->second;

						using utility::vector1;
						using numeric::model_quality::calc_rms;
						using numeric::model_quality::rms_wrapper;
						using core::scoring::xyz_gdtmm;

						int natoms = 0;
						ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
						ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );

						protocols::comparative_modeling::PartialThreadingMover pt(*it,template_pose);
						pt.apply(query_pose);
						SequenceAlignment aln = core::sequence::align_poses_naive(query_pose,native_pose);

						protocols::comparative_modeling::gather_coords(
							query_pose, native_pose, aln,
							natoms, p1a, p2a
						);

						Real const rmsd_ali  = rms_wrapper( natoms, p1a, p2a );
						Real const gdtmm_ali = xyz_gdtmm( p1a, p2a );
						Real const coverage  = (Real) natoms / (Real) native_pose.total_residue();
						Real const gdtmm_overall( gdtmm_ali * coverage );

						ss_out->add_energy( "coverage", coverage );
						ss_out->add_energy( "rmsd_ali", rmsd_ali );
						ss_out->add_energy( "gdtmm_ali", gdtmm_ali );
						ss_out->add_energy( "gdtmm_overall", gdtmm_overall );
						ss_out->add_energy( "n_ali_query", natoms );
						ss_out->add_energy( "nres_query", native_pose.total_residue() );
					} // template pdb check
				} // have native
				ss_out->scoreline_prefix( "" );
				ss_out->decoy_tag( it->sequence(2)->id() );
				ss_out->add_string_value( "template", template_id );
				sfd.write_silent_struct( *ss_out, option[ out::file::silent ]() );
			} // alns
		} // for ( it in aligns )

		tr.Debug << "finished rescoring alignments." << std::endl;
		tr.flush_all_channels();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
} // int main( int argc, char * argv [] )
