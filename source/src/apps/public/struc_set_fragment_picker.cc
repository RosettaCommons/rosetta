// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @details
///

#include <protocols/jobdist/Jobs.hh>
#include <core/types.hh>
#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <protocols/viewer/viewers.hh>
#include <numeric/random/random.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/BBTorsionSRFD.hh>

THREAD_LOCAL basic::Tracer tr( "struc_set_fragment_picker" );

using namespace core;
using namespace protocols;
using namespace fragment;
using namespace basic::options;

// add options
OPT_1GRP_KEY( Integer, struc_set_fragment_picker, frag_length )
OPT_1GRP_KEY( String, struc_set_fragment_picker, frag_name )
OPT_1GRP_KEY( Integer, struc_set_fragment_picker, sequence_length )

////////////////////////////////////////////////////////////////////////////////////////////
void
register_options()
{

	NEW_OPT( struc_set_fragment_picker::frag_length, "fragment length", 3 );
	NEW_OPT( struc_set_fragment_picker::frag_name, "fragment name", "struc_set_frag" );
	NEW_OPT( struc_set_fragment_picker::sequence_length, "sequence length", 30 );

}

void
run()
{

	using namespace basic::options;

	Size frag_length( option[ OptionKeys::struc_set_fragment_picker::frag_length ]() );
	Size sequence_length( option[ OptionKeys::struc_set_fragment_picker::sequence_length ]() );
	utility::vector1<utility::file::FileName> files( option[ OptionKeys::in::file::l ]() );
	Size frag_library_size( option[ OptionKeys::frags::n_frags ]() );
	utility::vector1<core::pose::Pose> poses;
	FragSetOP fragset(new ConstantLengthFragSet);

	for ( size_t i = 1; i <= files.size(); ++i ) {
		utility::io::izstream list( files[i] );
		std::string fname;
		while ( list >> fname ) {
			tr << "reading pose from " << fname << std::endl;
			core::pose::Pose tmp;
			core::import_pose::pose_from_file(tmp,fname, core::import_pose::PDB_file);
			core::scoring::dssp::Dssp dssp(tmp);
			dssp.insert_ss_into_pose(tmp);
			poses.push_back( tmp );
		}
	}

	FrameOP frame(new Frame( 1, frag_length ));
	for ( Size i=1; i <= poses.size(); ++i ) {
		for ( Size pos=1; pos <= poses[i].total_residue()-frag_length; ++pos ) {
			FragDataOP frag(new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), frag_length ));
			bool steal_success = frag->steal( poses[i], pos, pos + frag_length -1 );
			if ( steal_success && frag->is_valid() ) {
				bool success = frame->add_fragment( frag );
				assert( success );
			}
		}
	}
	for ( Size pos=1; pos <= sequence_length-frag_length+1; ++pos ) {
		FrameOP new_frame(new Frame( pos, frag_length ));
		for ( Size i=1; i<=frag_library_size; ++i  ) {
			Size random_frag_nbr = numeric::random::random_range( 1, frame->nr_frags() );
			FragDataCOP frag = (frame->fragment_ptr(random_frag_nbr) );
			new_frame->add_fragment( frag );
		}
		if ( new_frame->is_valid() ) {
			fragset->add( new_frame );
		}
	}

	utility::io::ozstream dump_frag( option[ OptionKeys::struc_set_fragment_picker::frag_name ]() );
	Size ct( 1 );
	for ( ConstFrameIterator it=fragset->begin(), eit=fragset->end(); it!=eit; ++it, ++ct ) {
		(*it)->show( dump_frag );
	}
}


void *
main_local( void* ) {

	run();
	return 0;

}


////////////////////////////////////////////////////////////////////////////////////////////


int
main( int argc, char * argv [] )
{
	try {
		register_options();
		devel::init( argc, argv );

		main_local(NULL);
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception caught: " << std::endl;
		excn.show( std::cerr );
		std::cerr.flush();
		return -1;
	}

	return 0;
}


