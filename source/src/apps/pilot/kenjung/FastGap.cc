// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ken Jung
/// @brief  lhalign
/// Given a partial threaded model pdb and a fasta, outputs the bin radius required to find a fragment
/// to fit inside each gap. Does not check for clashes between gaps for now.

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/pose/annotated_sequence.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>

#include <protocols/loophash/FastGapMover.hh>

#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

// AUTO-REMOVED #include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/conformation/util.hh>
#include <core/chemical/util.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/util/SwitchResidueTypeSet.hh>


// C++ headers
#include <map>
#include <devel/init.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


using namespace core;
using namespace io;
using namespace pdb;
using namespace chemical;
using namespace conformation;
using namespace protocols::loophash;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace utility;

using core::pose::PoseOP;
using core::pose::PoseAP;
using core::pose::Pose;
using utility::vector1;
using core::Size;
using std::string;

static thread_local basic::Tracer TR( "main" );

std::map< std::string, core::pose::Pose > poses_from_cmd_line(utility::vector1< std::string > const & fn_list){

  using utility::file::file_exists;
  using core::import_pose::pose_from_pdb;

  ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
  std::map< std::string, core::pose::Pose > poses;
  typedef utility::vector1< std::string >::const_iterator iter;
  for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
    if ( file_exists(*it) ) {
                        core::pose::Pose pose;
      core::import_pose::pose_from_pdb( pose, *rsd_set, *it );
                        std::string name = utility::file_basename( *it );
      name = name.substr( 0, 5 );
      poses[name] = pose;
    }
  }
  return poses;
}

int
main( int argc, char * argv [] )
{
    try {
    	// initialize core
    	devel::init(argc, argv);

      // read in poses
    	std::map< string, Pose > input_poses = poses_from_cmd_line(option[ in::file::s ]());

    	typedef std::map< string, Pose >::iterator iter;
    	for( iter it = input_poses.begin(), end = input_poses.end(); it != end; ++it ) {
    			FastGapMoverOP fastgap( new FastGapMover() );
    			fastgap->apply(it->second);
    			// dump resulting pdb
    			string m = it->first + ".closed.pdb";
    		 (it->second).dump_pdb(m);
    	}

    } catch ( utility::excn::EXCN_Base const & e ) {
			std::cerr << "caught exception " << e.msg() << std::endl;
			return -1;
    }
    return 0;
}


