//Get rmsd at interface and for all residue alpha carbons
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University
/// @brief
/// @author jk + dj

#include <iostream>
#include <iomanip>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/rms_util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//Need -interface_list -first_structure -second_structure
OPT_KEY( String, interface_list )
OPT_KEY( String, first_structure )
OPT_KEY( String, second_structure )

static thread_local basic::Tracer TR( "apps.pilot.karen_compare_pocket_rmsd.main" );

//set to store pdb info keys
std::set <std::string> interface;

//stores resid of the ligand residue
core::Size lig_res_num;

bool is_interface_heavyatom(
		       core::pose::Pose const & pose,
		       core::pose::Pose const & ,//pose2,
		       core::Size resno,
		       core::Size atomno
		       )
{
  // ws get residue "key" for set
  std::ostringstream residuestream;
  residuestream << pose.pdb_info()->chain(resno) << pose.pdb_info()->number(resno);
  std::string res_id = residuestream.str();

  core::conformation::Residue const & rsd = pose.residue(resno);

  if ( interface.count( res_id ) > 0 ) return rsd.is_protein() && !rsd.atom_is_hydrogen(atomno);

  return false;
}

/// General testing code
int
main( int argc, char * argv [] ){

	try{

  NEW_OPT( interface_list, "interface residues", "interface residues" );
  NEW_OPT( first_structure, "first structure", "-1" );
  NEW_OPT( second_structure, "comparison structure", "-1" );

  devel::init(argc, argv);

  std::string const ifilename = option[ interface_list ];
  if ( ifilename != "" ){
    std::ifstream ifs(ifilename.c_str(), std::ifstream::in);
    if (!ifs.is_open()){
      std::cout<< "Error opening contact list file "<<ifilename<<std::endl;
      return -100;
    }
    std::string intres;
    while (ifs.good()){
      ifs >> intres;
      interface.insert(intres);
    }
  }

  std::string const structure1 = option[ first_structure ] ;
  std::string const structure2 = option[ second_structure ] ;

  TR << "Starting recomputing scores and rmsds" << std::endl;

  // create pose from pdb
  pose::Pose pose1;
  core::import_pose::pose_from_pdb( pose1, structure1 );
  pose::Pose pose2;
  core::import_pose::pose_from_pdb( pose2, structure2 );


  TR << "Defined interface" << std::endl;

core::Real CA_rms = rmsd_with_super( pose1, pose2, is_protein_CA );
core::Real heavyatom_rms = rmsd_with_super( pose1, pose2, is_interface_heavyatom );

TR << "All residue rmsd: " << CA_rms << " Interface rmsd: " << heavyatom_rms <<std::endl;



TR << "Done computing rmsds" << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }

return 0;

}
