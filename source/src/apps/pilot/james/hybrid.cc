// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file hybrid.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/SequenceMapping.hh>

#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/rms_util.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>

#include <core/io/pdb/pdb_writer.hh>
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
#include <utility/io/mpistream.hh>
#include <utility/io/ozstream.hh>
#include <ObjexxFCL/format.hh>

#include <utility/excn/Exceptions.hh>

///////////////////////////////////////////////////////////////////////////////

using std::map;
using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::string_of;
using namespace basic;
using namespace core::sequence;
using namespace basic::options;
using namespace basic::options::OptionKeys;

void superimpose_via_alignment(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	core::sequence::SequenceAlignment aln,
	std::string const & atom_name = "CA"
) {
	using core::id::AtomID;
	using core::id::AtomID_Map;
	AtomID_Map< AtomID > atom_map;
  core::pose::initialize_atomid_map( atom_map, pose1, core::id::BOGUS_ATOM_ID );

	vector1< Size > residues;
	for ( Size ii = 1; ii <= pose1.size(); ++ii )
		residues.push_back(ii);

	core::id::SequenceMapping mapping( aln.sequence_mapping(1,2) );

	typedef vector1< Size >::const_iterator iter;
	for ( iter it = residues.begin(), end = residues.end(); it != end; ++it ) {
		Size const templ_ii( mapping[*it] );
		if ( templ_ii == 0 ) {
			continue;
		}
		if ( ! pose1.residue(*it).has(atom_name) ) continue;
		if ( ! pose2.residue(templ_ii).has(atom_name) ) continue;
		//std::cout << *it << " => " << templ_ii << std::endl;
		AtomID const id1( pose1.residue(*it).atom_index(atom_name), *it );
		AtomID const id2( pose2.residue(templ_ii).atom_index(atom_name), templ_ii );
		atom_map.set( id1, id2 );
	}

	using core::scoring::superimpose_pose;
	superimpose_pose( pose1, pose2, atom_map );
}

std::map< std::string, core::pose::Pose >
poses_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
	using std::map;
	using std::string;
	using core::pose::Pose;
	using utility::file::file_exists;
	using core::import_pose::pose_from_file;
	using namespace core::chemical;

	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
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
	using core::import_pose::pose_from_file;
	using namespace core::chemical;
	using namespace core::io::silent;

	basic::Tracer tr( "hybrid" );
	vector1< std::string > align_fns = option[ in::file::alignment ]();

	Pose native_pose;
	core::import_pose::pose_from_file(
		native_pose, *(rsd_set_from_cmd_line()), option[ in::file::native ]()
	);

	map< string, Pose > poses;
	poses = poses_from_cmd_line(
		option[ in::file::template_pdb ]()
	);

	Real const upper_dist_limit(12.0);
	string const atom_name("CA");
	string const query_seq( read_fasta_file( option[ in::file::fasta ]()[1] )[1]->sequence() );
	SilentFileData sfd;

	map< string, vector1< Real > > distances;

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
			string const aln_id( it->sequence(2)->id() );
			string const template_id( it->sequence(2)->id().substr(0,5) );
			tr << *it << std::endl;
			tr << "id " << it->sequence(2)->id() << " => " << template_id
				<< std::endl;
			string const ungapped_query( it->sequence(1)->ungapped_sequence() );

			Pose query_pose, template_pose;
			core::pose::make_pose_from_sequence(
				query_pose, query_seq, *(rsd_set_from_cmd_line())
			);

			map< string, Pose >::iterator pose_it = poses.find( template_id );
			if ( pose_it == poses.end() ) {
				string msg( "Error: can't find pose (id = " + template_id + ")" );
				//utility_exit_with_message(msg);
				tr.Error << msg << std::endl;
			} else {
				template_pose = pose_it->second;

				// build a partial threading model
				protocols::comparative_modeling::PartialThreadingMover mover(*it,template_pose);
				using core::sequence::align_poses_naive;
				mover.apply(query_pose);
				SequenceAlignment aln( align_poses_naive(query_pose, native_pose) );
				superimpose_via_alignment( query_pose, native_pose, aln, atom_name );

				core::id::SequenceMapping mapping( aln.sequence_mapping(2,1) );
				distances[aln_id] = vector1< Real >();
				for ( Size native_res = 1; native_res <= native_pose.size(); ++native_res ) {
					Size const query_res( mapping[native_res] );
					Real distance( upper_dist_limit );
					if ( mapping[native_res] ) {
						distance = native_pose.residue(native_res).xyz(atom_name).distance(
							query_pose.residue(query_res).xyz(atom_name)
						);
					}

					distances[aln_id].push_back( std::min( distance, upper_dist_limit ) );
				}
			} // template pdb check
		} // for alns
	} // for ( it in aligns )

	typedef map< string, vector1< Real > >::const_iterator dev_iter;
	vector1< string > aln_ids;
	for ( dev_iter it = distances.begin(), end = distances.end(); it != end; ++it ) {
		string const & aln_id( it->first );
		aln_ids.push_back(aln_id);
	}

	for ( Size ii = 1; ii <= native_pose.size(); ++ii ) {
		SilentStructOP ss_out( new ScoreFileSilentStruct );
		ss_out->scoreline_prefix( "" );
		ss_out->decoy_tag( string_of(ii) );

		Real min_dist( upper_dist_limit );
		for ( Size jj = 1; jj <= aln_ids.size(); ++jj ) {
			string const & aln_id( aln_ids[jj] );
			Real const & dist( distances[aln_id][ii] );
			ss_out->add_energy( "aln_" + aln_id, dist );

			min_dist = std::min( dist, min_dist );
		}
		ss_out->add_energy( "min_dist", min_dist );

		sfd.write_silent_struct( *ss_out, option[ out::file::silent ]() );
	}

	tr.Debug << "finished rescoring alignments." << std::endl;
	tr.flush_all_channels();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
