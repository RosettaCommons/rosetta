// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
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
#include <algorithm>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <basic/MetricValue.hh>
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

#include <core/pose/metrics/CalculatorFactory.hh>
#include <numeric/random/random.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>


// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );


OPT_KEY( String, fa_file )

//set to store pdb info keys
std::vector<std::string> surface;
std::set <std::string> interface;


/// General testing code
int
main( int argc, char * argv [] )
{
	try {

std::vector<std::string> surface;
std::set <std::string> interface;
  NEW_OPT ( fa_file, "File name for fasta data","");

  using namespace core;
  using namespace core::scoring;

  devel::init(argc, argv);
  pose::Pose pose;
  std::string const input_pdb_name ( basic::options::start_file() );
  core::import_pose::pose_from_pdb( pose, input_pdb_name );

  std::string const ffilename = option[ fa_file ];
  std::string seq;
  if ( ffilename != "" ){
    std::ifstream ifs(ffilename.c_str(), std::ifstream::in);
    if (!ifs.is_open()){
      std::cout<< "Error opening fasta file "<<ffilename<<std::endl;
      return -100;
    }
    while (ifs.good()){
      std::string intres;
      ifs >> intres;
      if (intres.length() >0){
        if (intres[0]!= '>'){
          seq.append(intres);
        }
      }
    }
  }

  if (seq.length()==0){
      std::cout<< "Fasta file contains no sequence: "<<ffilename<<std::endl;
      return -101;

  }

  int last = -100;
  core::chemical::ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );

  for( Size i = 1; i <= pose.total_residue(); ++i){
    if (last == -100){
      last = pose.pdb_info()->number(i);
      std::cout<<i<<" "<< pose.pdb_info()->number(i)<<" "<<pose.residue(i).name1()<<"\n" ;
      continue;
    }
    if (last > pose.pdb_info()->number(i)+1){
      for (int j = last-1; j > pose.pdb_info()->number(i); j--){
        std::cout<<seq[j-1]<<"\n";
        core::chemical::ResidueTypeCOP new_rsd_type( core::chemical::ResidueSelector().set_name1( seq[j-1] ).exclude_variants().select( rsd_set )[1] );
        core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );
        pose.append_polymer_residue_after_seqpos(*new_rsd, i, false );
        pose.pdb_info()->number(i+2, pose.pdb_info()->number(i)-1);
        pose.pdb_info()->chain(i+2, pose.pdb_info()->chain(i));
      }

    }else if (last < pose.pdb_info()->number(i)-1){
      for (int j = last+1; j < pose.pdb_info()->number(i); j++){
        std::cout<<seq[j-1]<<"\n";
        core::chemical::ResidueTypeCOP new_rsd_type( core::chemical::ResidueSelector().set_name1( seq[j-1] ).exclude_variants().select( rsd_set )[1] );
        core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );
        pose.append_polymer_residue_after_seqpos(*new_rsd, i-1, false );
        pose.pdb_info()->number(i, pose.pdb_info()->number(i-1)+1);
        pose.pdb_info()->chain(i, pose.pdb_info()->chain(i-1));
        i++;
      }


    }
    std::cout<<i<<" "<< pose.pdb_info()->number(i)<<" "<<pose.residue(i).name1()<<"\n" ;
    last = pose.pdb_info()->number(i);
  }

  for( Size i = 1; i <= pose.total_residue(); ++i){
    std::cout<<pose.pdb_info()->number(i)<<"\n";
  }
  for( int i = 1; i <pose.pdb_info()->number(1); ++i){
    std::cout<<" ";
  }
  for( Size i = 1; i <= pose.total_residue(); ++i){
    std::cout<<pose.residue(i).name1();
	}
  std::cout<<"\n"<<seq<<"\n";

  pose.pdb_info()->obsolete(false);
  pose.dump_pdb("test.pdb");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

  return 0;

}



