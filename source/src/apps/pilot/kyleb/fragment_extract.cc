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


// Example command line:
//  fragment_extract.macosgccrelease  -database ~/minirosetta_database -in:file:silent  region_10_18_sample.cluster.out  region_11_19_sample.cluster.out -frag_res 10 11 -o test.frags
//

// libRosetta headers
#include <core/scoring/Energies.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>

#include <protocols/viewer/viewers.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>
#include <core/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/io/pose_stream/ExtendedPoseInputStream.hh>
#include <core/io/pose_stream/PoseInputStream.fwd.hh>
#include <core/io/pose_stream/PDBPoseInputStream.hh>
#include <core/io/pose_stream/SilentFilePoseInputStream.hh>

#include <core/util/basic.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;
using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;
using utility::vector1;
using ObjexxFCL::format::F;
using io::pdb::dump_pdb;


OPT_KEY( IntegerVector, frag_res )

///////////////////////////////////////////////////////////////////////
// .Emacs, rosetta.el --> biox
// type "gtags" in mini/src/
// in emacs, F1, F2 lets you navigate.
// F8 = goto line

void
frag_extract_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::io::silent;
	using namespace core::chemical;
	using namespace core::io::pose_stream;


	PoseInputStreamOP input;
	utility::vector1< std::string > const & silent_files = option[ in::file::silent ]();

	// will get mysterious error if frag_res is both variable name and option name.
	utility::vector1< Size > const & which_res = option[ frag_res ]();

	if ( which_res.size()  != silent_files.size() ) utility_exit_with_message( " -in::file::silent and -frag_res must be specified, and have same number of elements!" );

	core::pose::Pose pose;

	core::chemical::ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set(
		option[ out::file::residue_type_set ]()
	);

	std::string outfile  = option[ out::file::o ]();
	utility::io::ozstream out( outfile );

	for ( Size n = 1; n <= silent_files.size(); n++ ) {

		// Could also hack this to accept normal PDBs as a PDBPoseInputStream ( see extract_pdbs.cc)
		input = new SilentFilePoseInputStream( silent_files[ n ] );

		out << "Position " << which_res[ n ] << std::endl;

		while ( input->has_another_pose() ) {

			input->fill_pose( pose, *rsd_set );
			if ( n == 1 ) protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

			for ( Size i = 1; i <= pose.total_residue(); i++ ) {
				out << ' ' << F( 8, 3, pose.phi( i ) );
				out << ' ' << F( 8, 3, pose.psi( i ) );
				out << ' ' << F( 8, 3, pose.omega( i ) );
				out << std::endl;
			}
			out << std::endl;

		}


	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	frag_extract_test();

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
    try {
	using namespace core::options;

	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;

	//Uh, options?
	NEW_OPT( frag_res, "residues for file", blank_size_vector ); //I am here.

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
	//exit( 0 );

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////
}
