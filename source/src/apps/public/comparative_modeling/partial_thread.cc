// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file partial_thread.cc
/// @brief
/// @author James Thompson
/// @author Ray Wang

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <core/chemical/util.hh>

#include <core/id/SequenceMapping.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <basic/options/option.hh>

#include <utility/io/ozstream.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <utility/excn/Exceptions.hh>


#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>

#include <protocols/comparative_modeling/PartialThreadingMover.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <map>
#include <iostream>
#include <string>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/partial_thread.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>

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

int
main( int argc, char* argv [] ) {
	try{
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
		using core::pose::make_pose_from_sequence;
		using protocols::comparative_modeling::PartialThreadingMover;

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using namespace core::sequence;

		basic::Tracer tr( "partial_thread" );

		SequenceOP fasta_seq = core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1];
		map< string, Pose > poses = poses_from_cmd_line( option[ in::file::template_pdb ]() );
		bool skip_repack = option[ partial_thread::skip_repack ]();

		vector1< string > align_fns = option[ in::file::alignment ]();
		using aln_iter = vector1<string>::const_iterator;
		for ( aln_iter aln_fn = align_fns.begin(), aln_end = align_fns.end(); aln_fn != aln_end; ++aln_fn ) {
			vector1< SequenceAlignment > alns = core::sequence::read_aln(option[ cm::aln_format ](), *aln_fn);

			//fd  One thing kind of dumb about this is that:
			//fd     "input id" is it->sequence(2)->id().substr(0,5)
			//fd     "output id" is it->sequence(2)->id()
			//fd  So input and outputs are "linked" and it is really easy to overwrite your inputs
			for ( auto & aln : alns ) {
				string const template_id( aln.sequence(2)->id().substr(0,5) );

				string const ungapped_query( aln.sequence(1)->ungapped_sequence() );

				map< string, Pose >::iterator pose_it = poses.find( template_id );
				if ( pose_it == poses.end() ) {
					tr.Error << "can't find pose (id = " << template_id << ")" << std::endl;
				} else {
					tr << aln << std::endl;
					tr << "id " << aln.sequence(2)->id() << " => " << template_id
						<< std::endl;

					Pose query_pose, template_pose;
					make_pose_from_sequence( query_pose, fasta_seq->sequence(), *(rsd_set_from_cmd_line().lock()) );
					template_pose = pose_it->second;
					PartialThreadingMover mover(aln,template_pose,skip_repack);
					mover.apply(query_pose);
					string const id_out( aln.sequence(2)->id() );

					// print out query-anchored alignment
					utility::io::ozstream output( id_out+".pdb" );

					core::id::SequenceMapping map( aln.sequence_mapping(1,2) );
					output << "REMARK query_anchored_aln ";
					for ( core::Size ii = 1; ii <= fasta_seq->sequence().size(); ++ii ) {
						if ( map[ii] ) output << fasta_seq->at(ii);
						else           output << "-";
					}
					output << std::endl;
					core::io::pdb::dump_pdb( query_pose, output );
					output.close();

				} // template pdb check
			} // alns
		} // for ( it in aligns )

		tr.Debug << "finished building partial models." << std::endl;
		tr.flush_all_channels();
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
} // int main( int argc, char * argv [] )
