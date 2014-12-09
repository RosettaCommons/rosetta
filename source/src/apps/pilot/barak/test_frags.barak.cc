// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   FlexPepDock.cc
//
/// @brief Main file for applying flexible peptide docking protocol
/// @author Barak Raveh
/// @date August 05, 2008

//#define GL_GRAPHICS

// Mini-Rosetta headers

#include <numeric/random/random.hh>

#include <devel/init.hh>
#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/io/pdb/pose_io.hh>

#include <basic/options/util.hh>//option.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>

//#include <protocols/jobdist/standard_mains.hh>
//#include <protocols/moves/Mover.hh>
#include <protocols/frags/VallData.hh>
#include <protocols/frags/TorsionFragment.hh>

#include <ObjexxFCL/format.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>

using basic::T;
using basic::Error;
using basic::Warning;




///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

  using namespace core;
  using namespace protocols::frags;

  devel::init(argc, argv);

  // test reading of file
  VallData vall( "/vol/klee/ravehb/tmp/flexpeps/frags/aa2fm8C03_05.200_v1_3" );
  std::cout << "Part 1: " << std::endl;
 {
    SingleResidueTorsionFragmentLibrary lib;
    vall.get_frags(
      500/*nstruct*/, "LLL"/*target seq*/, "---" /*target_ss*/,
      1.0/*seq_weight*/, 0.0/*ss_weight*/,
      true/*XGly*/, true/*XPro*/, true/*XCys*/,
      lib );
  }


  // test matching to protein of 3-mers
  std::cout << "Part 2: " << std::endl;
  {
    pose::Pose pose;
    core::import_pose::pose_from_pdb( pose, basic::options::start_file() );
    std::string const sequence( pose.sequence() );

    TorsionFragmentLibrary lib;
    Size const nres( pose.total_residue() );
    Size const frag_size(3);
    lib.resize( nres - frag_size + 1 );
    for ( Size i=1; i<= nres-frag_size+1; ++i )
      {
	std::string const frag_seq( sequence.substr(i-1,3) );
	vall.get_frags( 200, frag_seq, "---", 1.0, 0.0, false, false, true, lib[i] );
      }
  }


  { // try building an ideal peptide:
    using namespace pose;
    using namespace chemical;
    using namespace conformation;

    Pose pose;
    // read centroid residue set
    chemical::ResidueTypeSetCAP rsd_set( chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" ) );
    for ( Size i=1; i<= 20; ++i ) {
      ResidueTypeCOPs const & rsd_list( rsd_set->aa_map( static_cast<AA>(i) ) /*BAD*/ );
      for ( ResidueTypeCOPs::const_iterator iter=rsd_list.begin(), iter_end= rsd_list.end(); iter!= iter_end; ++iter ) {
	ResidueType const & rsd_type( **iter );
	if ( ( rsd_type.is_lower_terminus() == ( i == 1 ) ) &&
	  ( rsd_type.is_upper_terminus() == ( i == 20 ) ) ) {
	  ResidueOP new_rsd( ResidueFactory::create_residue( rsd_type ) );
	  pose.append_residue_by_bond( *new_rsd, true );
	}
      }
    }
    for ( Size i=1; i<= 20; ++i ) {
      pose.set_omega(i,180.0);
    }
    io::pdb::dump_pdb( pose, "test_ideal.pdb" );
  }


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
