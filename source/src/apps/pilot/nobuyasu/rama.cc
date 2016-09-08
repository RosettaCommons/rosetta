// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/nobuyasu/make_blueprint.cc
/// @brief makes blueprint file
/// @author Nobuyasu Koga ( 02/07/2010 )

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/ABEGOManager.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/moves/DsspMover.hh>
#include <utility/excn/Exceptions.hh>

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>

#include <fstream>

static THREAD_LOCAL basic::Tracer TR( "rama" );

typedef core::Size Size;
typedef std::string String;
using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

OPT_KEY( Integer, abego )
OPT_KEY( File, output )

void ThisApplication::register_options() {
	OPT( in::file::s );
	NEW_OPT( abego, "abego output level [ 1-3 ] ", 1 );
	NEW_OPT( output, "output filename", "rama.dat" );
}


int
main( int argc, char * argv [] )
{
	try{
		using core::chemical::oneletter_code_from_aa;

		ThisApplication::register_options();
		devel::init(argc, argv);

		// blueprint output file
		std::ofstream out;
		std::ostringstream filename;
		filename <<  basic::options::option[ output ]();
		out.open( filename.str().c_str() ,std::ios::out );

		// calc secondary structure info
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, option[ in::file::s ].value().at( 1 ) , core::import_pose::PDB_file);
		protocols::moves::DsspMover dsm;
		dsm.apply( pose );

		utility::vector1< std::string > abego_ = core::sequence::get_abego( pose, option[ abego ]() );

		using namespace ObjexxFCL::format;
		out << "# resn aa ss abego phi psi omega" << std::endl;
		for ( core::Size ii=1; ii<=pose.size(); ii++ ) {
			out << I( 5, pose.pdb_info()->number( ii ) ) << " " << oneletter_code_from_aa( pose.aa( ii ) ) << " " << pose.secstruct( ii ) << " " << abego_[ ii ] << " "
				<< F( 8, 2, pose.phi( ii ) ) << F( 8, 2, pose.psi( ii ) ) << F( 8, 2, pose.omega( ii ) )
				<< std::endl;
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

