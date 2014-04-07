// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file thread_fast.cc
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
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>

#include <protocols/comparative_modeling/ThreadingMover.hh>
#include <protocols/jobdist/not_universal_main.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/model_quality/rms.hh>

// C++ headers
#include <map>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <protocols/jobdist/Jobs.hh>
#include <utility/io/mpistream.hh>
#include <ObjexxFCL/format.hh>

#include <utility/excn/Exceptions.hh>



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

	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
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

	basic::Tracer tr( "thread_fast" );

	SequenceProfileOP query_prof( new SequenceProfile );
	query_prof->read_from_file( option[ in::file::pssm ]()[1] );
	query_prof->convert_profile_to_probs( 1.0 ); // was previously implicit in read_from_file()
	vector1< std::string > align_fns = option[ in::file::alignment ]();

	Pose native_pose;
	bool have_native( false );
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb(
			native_pose,
			*(rsd_set_from_cmd_line()),
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

	map< string, SequenceOP > seqs = sequences_from_cmd_line(
		option[ in::file::pssm ]()
	);

	// scoring scheme for aligning profiles
	std::string const scoring_scheme_type( option[ cm::seq_score ]() );
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

			map< string, SequenceOP >::iterator seq_it = seqs.find( template_id );
			if ( seq_it == seqs.end() ) {
				print_seq_map( std::cerr, seqs );
				string msg( "Error: can't find seq (id = " + template_id + ")" );
				//utility_exit_with_message(msg);
				std::cerr << msg << std::endl;
				continue;
			}
			SequenceOP template_prof = seq_it->second;

			utility::vector1< SequenceOP > my_seqs;
			my_seqs.push_back( query_prof->clone() );
			my_seqs.push_back( template_prof->clone() );
			SequenceAlignment rescore_aln = steal_alignment( *it, my_seqs );
			Real const aln_score( rescore_aln.calculate_score_sum_of_pairs( ss ) );
			rescore_aln.score( aln_score );
			string const ungapped_query( it->sequence(1)->ungapped_sequence() );
			//string const ungapped_templ( it->sequence(2)->ungapped_sequence() );

			SilentStructOP ss_out( new ScoreFileSilentStruct );
			ss_out->decoy_tag( it->sequence(2)->id() );
			ss_out->add_energy( "aln_score", aln_score );
			ss_out->add_energy( "n_ali_query", ungapped_query.length() );
			//ss_out->add_energy( "n_ali_templ", ungapped_templ.length() );
			tr.Debug << "score(" << ss_out->decoy_tag() << ") = " << aln_score
				<< std::endl;

			using core::pose::setPoseExtraScores;
			if ( have_native ) {
				// calc rmsd/gdt stats

				Pose query_pose, template_pose;
				core::pose::make_pose_from_sequence(
					query_pose,
					ungapped_query,
					*(rsd_set_from_cmd_line())
				);

				map< string, Pose >::iterator pose_it = poses.find( template_id );
				if ( pose_it == poses.end() ) {
					string msg( "Error: can't find pose (id = "
						+ template_id + ")"
					);
					//utility_exit_with_message(msg);
					std::cerr << msg << std::endl;
					continue;
				}
				template_pose = pose_it->second;

				static string const atm( "CA" );
				SequenceMapping mapping = it->sequence_mapping(1,2);
				using utility::vector1;
				using numeric::xyzVector;
				using numeric::model_quality::calc_rms;
				using numeric::model_quality::rms_wrapper;
				using core::scoring::xyz_gdtmm;
				vector1< xyzVector< Real > > native_coords, template_coords;

				int natoms(0);
				for ( Size ii = 1; ii <= native_pose.total_residue(); ++ii ) {
					Size const templ_ii( mapping[ii] );
					if ( templ_ii != 0 ) ++natoms;
				}

				ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
				ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
				Size n_gap(0);
				for ( Size ii = 1; ii <= native_pose.total_residue(); ++ii ) {
					Size const templ_ii( mapping[ii] );
					if ( templ_ii == 0 ) {
						n_gap++;
						continue;
					}

					core::Vector native_xyz  ( native_pose.residue(ii).xyz(atm) );
					core::Vector template_xyz( template_pose.residue(templ_ii).xyz(atm) );
					native_coords.push_back(native_xyz);
					template_coords.push_back(template_xyz);
					for ( Size jj = 1; jj <= 3; ++jj ) {
						p1a(jj,ii - n_gap) = native_xyz  [jj-1];
						p2a(jj,ii - n_gap) = template_xyz[jj-1];
					}
				}
				runtime_assert( native_coords.size() == template_coords.size() );

				Real const debug_rmsd_ali( calc_rms( native_coords, template_coords ) );
				Real const rmsd_ali ( rms_wrapper( natoms, p1a, p2a ) );
				Real const gdtmm_ali( xyz_gdtmm( p1a, p2a ) );
				Real const coverage(
					(Real) ungapped_query.length() / (Real) native_pose.total_residue()
				);

				//Real const gdtmm_ali( 1.0 );
				ss_out->add_energy( "coverage", coverage );
				ss_out->add_energy( "rmsd_ali", rmsd_ali );
				ss_out->add_energy( "gdtmm_ali", gdtmm_ali );
				ss_out->add_energy( "debug_rmsd_ali", debug_rmsd_ali );
			}
			ss_out->add_string_value( "template", template_id );
			ss_out->scoreline_prefix( "" );
			sfd.write_silent_struct( *ss_out, option[ out::file::silent ]() );
		} // alns
	} // for ( it in aligns )

	tr.Debug << "finished rescoring alignments." << std::endl;
	tr.flush();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
