// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson


#include <core/types.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/MCAligner.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/SimpleScoringScheme.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/ExtendedPoseInputStream.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/comparative_modeling/ThreadingMover.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


///////////////////////////////////////////////////////////////////////////////

using core::Size;
using utility::file::FileName;
using namespace core::sequence;
using namespace basic::options;
using namespace core::chemical;
using namespace basic::options::OptionKeys;

int
main( int argc, char* argv [] )
{
	try {

	// options, random initialization
	devel::init( argc, argv );

	utility::vector1< utility::file::FileName >
		align_files( option[ in::file::alignment ]() );


	typedef utility::vector1< SequenceAlignment > alignlist;
	alignlist aligns
		= read_general_aln_file( option[ in::file::alignment ]()[1] );

	core::pose::Pose template_pose;

	using namespace protocols::comparative_modeling;
	using namespace protocols::jobdist;

	SequenceAlignment true_aln;
	if ( align_files.size() > 1 ) {
		std::string align_fn  = option[ in::file::alignment ]()[2];
		true_aln.read_from_file( align_fn );
		std::cout << true_aln << std::endl;
		true_aln.comment( "mammoth_alignment" );

		protocols::comparative_modeling::ThreadingMover m(
			true_aln,
			template_pose
		);
		protocols::jobdist::not_universal_main( m );
	}

	std::cout << "aligns.size() = " << aligns.size() << std::endl;
	for ( alignlist::iterator it = aligns.begin(), end = aligns.end();
			it != end; ++it
	) {
			std::string template_id = it->sequence(2).id();
			core::import_pose::pose_from_file(
				template_pose,
				template_id
			);
			std::cout << "template_id = " << template_id << std::endl;
			std::cout << "template_sequence = " << template_pose.sequence() << std::endl;
			it->comment( template_id );
			protocols::comparative_modeling::ThreadingMover mover(
				*it,
				template_pose
			);

			protocols::jobdist::not_universal_main( mover );
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
