// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @brief
/// @author simon kelow - simon.kelow@gmail.com

// Project Headers
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <fstream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/file_data.hh>

using core::Size;
using core::Real;
OPT_KEY( String, mov_pose )

static basic::Tracer TR( "apps.pilot.simon_superimpose.main" );

Real
calpha_pdb_superimpose_pose(
														core::pose::Pose & mod_pose,
														core::pose::Pose const & ref_pose
														)
{
 core::id::AtomID_Map< core::id::AtomID > atom_map;
 core::pose::initialize_atomid_map( atom_map, mod_pose, core::id::BOGUS_ATOM_ID );
  for ( Size ii = 1; ii <= mod_pose.total_residue(); ++ii ) {
    if ( ! mod_pose.residue(ii).has("CA") ) continue;
    if ( ! mod_pose.residue(ii).is_protein() ) continue;
    for ( Size jj = 1; jj <= ref_pose.total_residue(); ++jj ) {
      if ( ! ref_pose.residue(jj).has("CA") ) continue;
      if ( ! ref_pose.residue(jj).is_protein() ) continue;
      if ( mod_pose.pdb_info()->chain(ii) != ref_pose.pdb_info()->chain(jj)) continue;
      if ( mod_pose.pdb_info()->number(ii) != ref_pose.pdb_info()->number(jj)) continue;
      core::id::AtomID const id1( mod_pose.residue(ii).atom_index("CA"), ii );
      core::id::AtomID const id2( ref_pose.residue(jj).atom_index("CA"), jj );
      atom_map.set( id1, id2 );
      break;
    }

  }
  return core::scoring::superimpose_pose( mod_pose, ref_pose, atom_map );
}


int main( int argc, char * argv [] ){
	try{
		NEW_OPT( mov_pose, "the movable structure", "" );

		devel::init(argc, argv);

		// create reference pose from pdb
		core::pose::Pose ref_pose;
		std::string const ref_pdb_name( basic::options::start_file() );
        	core::import_pose::pose_from_pdb( ref_pose, ref_pdb_name );

		// create movable pose from pdb
		core::pose::Pose mov_pose;
		std::string const mov_pdb_name( basic::options::option[basic::options::OptionKeys::mov_pose] );
		core::import_pose::pose_from_pdb( mov_pose, mov_pdb_name );
		
		// superimpose mov_pose to ref_pose
		calpha_pdb_superimpose_pose( mov_pose, ref_pose );

		// print superimposed pose to output
		std::ofstream super_os( "superposed.pdb", std::ios_base::out ); 
		core::io::pdb::FileData::dump_pdb( mov_pose, super_os );
	}
	catch ( utility::excn::EXCN_Base const & e ){
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	return 0; 

}
