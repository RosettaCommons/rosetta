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

#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/types.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobInputter.hh>

#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <protocols/toolbox/DecoySetEvaluation.impl.hh>
#include <protocols/toolbox/InteratomicVarianceMatrix.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/toolbox/Cluster.hh>
#include <protocols/toolbox/Cluster.impl.hh>

#include <protocols/loops/Loops.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>


#include <core/id/AtomID.hh>

#include <core/chemical/ChemicalManager.hh>

// Auto-header: duplicate removed #include <protocols/loops/Loops.hh>

// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

// ObjexxFCL includes
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <vector>
#include <ostream>
#include <algorithm>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <utility/excn/Exceptions.hh>


static basic::Tracer tr( "main" );

using namespace core;
using namespace protocols;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace toolbox;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;


// OPT_KEY( IntegerVector, reslist )
// OPT_KEY( File, points )
// OPT_KEY( String, bfac )

// void register_options(){
//  using namespace basic::options;
//  using namespace basic::options::OptionKeys;
//  OPT( in::file::silent );
//  OPT(out::file::residue_type_set);
//  NEW_OPT( reslist, "extract CA coords of atoms in reslist", 1 );
//   NEW_OPT( points, "which atom to extract", "test.pdb" );
//  NEW_OPT( bfac, "which energy column to put in bfactor","");
// }


// void write_for_resnum(int resnum, char chainID){
//  using namespace basic::options;
//  using namespace basic::options::OptionKeys;
//  using namespace io::silent;
//  // string chainID = "Z";
//  SilentFileData sfd;
//  sfd.read_file( option[ in::file::silent ]()[1] );
//  int ct = 1;
//  utility::io::ozstream out( option[points ]() );
//  for ( SilentFileData::iterator it = sfd.begin(); it!=sfd.end(); ++it, ++ct ) {
//   ObjexxFCL::FArray2D< Real > coords( it->get_CA_xyz() );
//   Real x,y,z;
//    x = coords( 1, resnum );
//    y = coords( 2, resnum );
//    z = coords( 3, resnum );
//    char outbuf[200];

//    sprintf(outbuf, "ATOM  %5d %4s %s %s%4d    %8.3f%8.3f%8.3f  1.00  1.00\n", ct, "CA", "TRJ","Z",ct, x, y, z);

//    out << outbuf ;
//  }
// }

// void run(){
//  using namespace basic::options;
//  using namespace basic::options::OptionKeys;
//  using namespace io::silent;
//  SilentFileData sfd;
//  core::pose::Pose pose;
//  core::chemical::ResidueTypeSetCAP rsd_set;
//  int ct = 1;
//  if (!option[ reslist ]().empty() ){
//   rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( option[ out::file::residue_type_set ]() );
//   sfd.read_file( option[ in::file::silent ]()[1] );

//   sfd.begin()->fill_pose( pose, *rsd_set );
//    pose.dump_pdb("./firstPose.pdb");   // make it easier to later compare to the right coordinates system.

//   utility::io::ozstream out( option[points ]() );
//   for ( SilentFileData::iterator it = sfd.begin(); it!=sfd.end(); ++it, ++ct ){
//    ObjexxFCL::FArray2D< Real > coords( it->get_CA_xyz() );
//    Real x,y,z;
//    Real energy;
//    char outbuf[200];
//    char chainID ='Z';
//    for ( int i=1; i <= option[ reslist ]().size(); i++ ){
//     int res = option[ reslist ]()[i];
//     if ( !option[ bfac ]().empty() && it->has_energy( option[ bfac ]() )){
//       energy = it->get_energy( option[ bfac ]() );
//     }
//     else {
//      energy = 1.00;
//     }
//     x = coords( 1, res );
//     y = coords( 2, res );
//     z = coords( 3, res );
//     chainID = chainID - i + 1;
//     //    sprintf(outbuf, "ATOM  %5d %4s %s %s%4d    %8.3f%8.3f%8.3f  1.00  1.00\n", ct, "CA", "TRJ","Z",ct, x, y, z);
//     //     sprintf(outbuf, "ATOM  %5d %4s %s %c%4d    %8.3f%8.3f%8.3f  1.00  1.00\n", 1, "CA", "TRJ",chainID,1, x, y, z);
//      sprintf(outbuf, "ATOM  %5d %4s %s %c%4d    %8.3f%8.3f%8.3f  1.00  %3.2f\n", 1, "CA", "TRJ",chainID,1, x, y, z,energy);
//     out << outbuf;
//    }
//   }
//  }
// }


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		// register_options();
		devel::init( argc, argv );
		tr.Trace << "test in main" << std::endl;

		// string sequence = "CASEDELVAEFLQDQN";
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose,"CASEDELVAEFLQDQN", core::chemical::FA_STANDARD, true );
		pose.dump_pdb();

		//  try{
		//   run();
		//  } catch (utility::excn::Exception& excn ) {
		//   excn.show( std::cerr );
		//  }
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


