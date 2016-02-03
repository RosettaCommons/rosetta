// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/types.hh>

#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>

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
#include <iomanip>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <utility/excn/EXCN_Base.hh>


static THREAD_LOCAL basic::Tracer tr( "main" );

using namespace core;
using namespace protocols;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace toolbox;
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
 	NEW_OPT( points, "which atom to extract", "test.pdb" );
	NEW_OPT( bfac, "which energy column to put in bfactor","");
}

core::io::silent::SilentStructOP
extract_replica(std::string filename, std::string jobname, core::Real temp_level)
{
        using namespace io::silent;
        SilentFileData sfd;
        SilentFileData::iterator final;
        Real energy;

        sfd.read_file( filename );
        for ( SilentFileData::iterator it=sfd.begin(); it!=sfd.end(); ++it)
        {
                if ( it->has_energy( "temp_level" ) )
                {
                        energy = it->get_energy( "temp_level" );
                        if ( energy == temp_level ) final = it;
                }

        }

        final->print_score_header( std::cout );
        return *final;
}

core::io::silent::SilentStructOP
extract_replica(std::string filename, std::string jobname, int replica_id)
{
        using namespace io::silent;
				using namespace std;
        SilentFileData sfd;
        SilentFileData::iterator final;
				ostringstream replica_id_str;
				replica_id_str << setw(3) <<setfill('0') << replica_id;
				std::string tag = jobname+"_"+replica_id_str.str();
				std::string full_tag;
				std::string::size_type pos;
        sfd.read_file( filename );

        for ( SilentFileData::iterator it=sfd.begin(); it!=sfd.end(); ++it)
        {
					full_tag = it->decoy_tag();
					pos = full_tag.find( tag );
					if ( pos!=std::string::npos ){
						final = it;
					}
        }

        final->print_score_header( std::cout );
        return *final;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

// 	register_options();
 	devel::init( argc, argv );
	tr.Trace << "test in main" << std::endl;
	core::io::silent::SilentStructOP ss = extract_replica( "for_test.out" , "test" , 1.000);
	ss->print_scores( std::cout );

	core::io::silent::SilentStructOP st = extract_replica( "for_test.out" , "P_0002" , 1);
	st->print_scores( std::cout );
	 } catch ( utility::excn::EXCN_Base const & e ) {
		 std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}


