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
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
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
#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>


// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );


OPT_KEY( Integer, max_residues )
OPT_KEY( Real, min_sasa )
OPT_KEY( String, contact_list )

//set to store pdb info keys
std::vector<std::string> surface;
std::set <std::string> interface;

//stores resid of the ligand residue
void register_metrics() {

  core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new core::pose::metrics::simple_calculators::SasaCalculator;
  core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );

}

/// General testing code
int
main( int argc, char * argv [] )
{
std::vector<std::string> surface;
std::set <std::string> interface;
  NEW_OPT( max_residues, "Maximum number of residues to return", 1000);
  NEW_OPT( min_sasa, "Minimum SASA to classify a residue as surface", 10);
  NEW_OPT ( contact_list, "File name for optional list of contact residues to check","");

  using namespace core;
  using namespace core::scoring;

  devel::init(argc, argv);
  register_metrics();
  pose::Pose pose;
  std::string const input_pdb_name ( basic::options::start_file() );
  core::import_pose::pose_from_pdb( pose, input_pdb_name );

  std::string const cfilename = option[ contact_list ];
  if ( cfilename != "" ){
    std::ifstream ifs(cfilename.c_str(), std::ifstream::in);
    if (!ifs.is_open()){
      std::cout<< "Error opening contact list file "<<cfilename<<std::endl;
      return -100;
    }
    //ifb.open (cfilename,std::ios::in);
    //std::ostream ios(&ifb);
    std::string intres;
    while (ifs.good()){
      ifs >> intres;
      interface.insert(intres);
    }
  }

  core::Size sasa_threshold = option [ min_sasa ];
  core::Size max_resi = option [ max_residues ];
  basic::MetricValue< utility::vector1< Real > > resisasa;
  pose.metric( "sasa", "residue_sasa", resisasa );

  for( Size i = 1; i <= pose.total_residue(); ++i){
    if (resisasa.value()[i]>= sasa_threshold){
      std::ostringstream residuestream;
      residuestream << pose.pdb_info()->chain(i) << pose.pdb_info()->number(i);
      std::string res_id = residuestream.str();
      if (interface.find(res_id) == interface.end()){
        //std::cout<<res_id<<std::endl;
        surface.push_back(res_id);
      }
    }
  }

  std::filebuf fb;
  std::stringstream filename;
  filename<<option[ OptionKeys::out::output_tag ]()<<".randomsurface";
  fb.open (filename.str().c_str(),std::ios::out);
  std::ostream os(&fb);


  if ( max_resi < surface.size()){
    std::set<int> indeces;
    for (core::Size i = 0; i < max_resi; i++){
      int r=(int) (numeric::random::uniform() * surface.size());
      if (indeces.find(r) == indeces.end()){
        indeces.insert(r);
        /*if (indeces.size() == 0) indeces.insert(i);
        else{
          for (std::set<int>::iterator it=indeces.begin(); true; it++){
            if (it == indeces.end()) {
              indeces.insert(i);
              break;
            }
            if (r < *it){
              indeces.insert(it, r);
              break;
            }
          }
        }
      */
      }else{
        i--;
      }
    }
    //sort (indeces.begin(), indeces.end());
    for (std::set<int>::iterator it=indeces.begin(); it != indeces.end(); it++){
      os<<surface[*it]<<std::endl;
    }
  }else{
    for (std::vector<std::string>::iterator it=surface.begin(); it != surface.end(); it++){
      os<<*it<<std::endl;
    }
  }

	return 0;

}



