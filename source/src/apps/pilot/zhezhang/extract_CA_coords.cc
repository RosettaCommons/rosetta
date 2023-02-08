// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange


#include <core/pose/Pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <devel/init.hh>
#include <core/types.hh>







#include <core/chemical/ChemicalManager.hh>

// Auto-header: duplicate removed #include <protocols/loops/Loops.hh>

// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// ObjexxFCL includes
#include <ObjexxFCL/format.hh>
#include <ostream>
#include <string>

//Auto Headers
#include <utility/excn/Exceptions.hh>

#include <core/io/silent/SilentStruct.hh> // AUTO IWYU For SilentStruct
#include <ObjexxFCL/FArray2D.hh> // AUTO IWYU For FArray2D, FArray2D<>::size_type, FArray2D::IR


static basic::Tracer tr( "main" );

using namespace core;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;


OPT_KEY( IntegerVector, reslist )
OPT_KEY( File, points )
OPT_KEY( String, bfac )

void register_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::silent );
	OPT(out::file::residue_type_set);
	NEW_OPT( reslist, "extract CA coords of atoms in reslist", 1 );
	NEW_OPT( points, "output file name", "test.pdb" );
	NEW_OPT( bfac, "which energy column to put in bfactor","");
}


void write_for_resnum(int resnum, char /*chainID*/){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace io::silent;
	// string chainID = "Z";
	SilentFileOptions opts; // initialized from the command line
	SilentFileData sfd( opts );
	sfd.read_file( option[ in::file::silent ]()[1] );
	int ct = 1;
	utility::io::ozstream out( option[points ]() );
	for ( SilentFileData::iterator it = sfd.begin(); it!=sfd.end(); ++it, ++ct ) {
		ObjexxFCL::FArray2D< Real > coords( it->get_CA_xyz() );
		Real x,y,z;
		x = coords( 1, resnum );
		y = coords( 2, resnum );
		z = coords( 3, resnum );
		char outbuf[200];

		snprintf(outbuf, sizeof(outbuf), "ATOM  %5d %4s %s %s%4d    %8.3f%8.3f%8.3f  1.00  1.00\n", ct, "CA", "TRJ","Z",ct, x, y, z);

		out << outbuf ;
	}
}

void run(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace io::silent;
	SilentFileOptions opts; // initialized from the command line
	SilentFileData sfd( opts );
	core::pose::Pose pose;
	core::chemical::ResidueTypeSetCOP rsd_set;
	int ct = 1;
	if ( !option[ reslist ]().empty() ) {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( option[ out::file::residue_type_set ]() );
		sfd.read_file( option[ in::file::silent ]()[1] );

		sfd.begin()->fill_pose( pose, *rsd_set );
		pose.dump_pdb("./firstPose.pdb");   // make it easier to later compare to the right coordinates system.

		utility::io::ozstream out( option[points ]() );
		for ( SilentFileData::iterator it = sfd.begin(); it!=sfd.end(); ++it, ++ct ) {
			ObjexxFCL::FArray2D< Real > coords( it->get_CA_xyz() );
			Real x,y,z;
			Real energy;
			char outbuf[200];
			char chainID ='Z';
			for ( Size i=1; i <= option[ reslist ]().size(); i++ ) {
				int res = option[ reslist ]()[i];
				if ( !option[ bfac ]().empty() && it->has_energy( option[ bfac ]() ) ) {
					energy = it->get_energy( option[ bfac ]() );
				} else {
					energy = 1.00;
				}
				x = coords( 1, res );
				y = coords( 2, res );
				z = coords( 3, res );
				chainID = chainID - i + 1;
				//    snprintf(outbuf, sizeof(outbuf), "ATOM  %5d %4s %s %s%4d    %8.3f%8.3f%8.3f  1.00  1.00\n", ct, "CA", "TRJ","Z",ct, x, y, z);
				//     snprintf(outbuf, sizeof(outbuf), "ATOM  %5d %4s %s %c%4d    %8.3f%8.3f%8.3f  1.00  1.00\n", 1, "CA", "TRJ",chainID,1, x, y, z);
				snprintf(outbuf, sizeof(outbuf), "ATOM  %5d %4s %s %c%4d    %8.3f%8.3f%8.3f  1.00  %3.2f\n", 1, "CA", "TRJ",chainID,1, x, y, z,energy);
				out << outbuf;
			}
		}
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		register_options();
		devel::init( argc, argv );
		tr.Trace << "test in main" << std::endl;

		run();
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}


