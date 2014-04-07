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
#include <core/scoring/EnergyGraph.hh>

#include <core/scoring/rms_util.hh>

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

OPT_KEY( String, apo_chain )

static basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );

//set to store pdb info keys
std::set <std::string> interface;
//stores resid of the ligand residue

//mjo commenting out 'input_pose' because it is unused and causes a warning
void define_interface( core::pose::Pose & /*input_pose*/ ) {
}

/// General testing code
int
main( int argc, char * argv [] )
{
	try {

  NEW_OPT ( apo_chain, "Chain to find contacts on","A");
  devel::init(argc, argv);
  pose::Pose input_pose;
  std::string chain = option [ apo_chain ];

  //read in pdb file from command line
  std::string const input_pdb_name ( basic::options::start_file() );
  core::import_pose::pose_from_pdb( input_pose, input_pdb_name );

  std::filebuf fb;
  std::stringstream filename;
  filename<<option[ OptionKeys::out::output_tag ]()<<".contacts";
  fb.open (filename.str().c_str(),std::ios::out);
  std::ostream os(&fb);

  scoring::ScoreFunctionOP scorefxn( getScoreFunction() );
  (*scorefxn)(input_pose);

  EnergyGraph & energy_graph(input_pose.energies().energy_graph());

  for ( int i = 1, resnum = input_pose.total_residue(); i <= resnum; ++i ) {
  std::ostringstream tmp;
  tmp<< input_pose.pdb_info()->chain(i);
  std::string rch = tmp.str();
  //if (rch != chain) continue;
	 //std::cout<<std::endl<<input_pose.pdb_info()->number(i)<<":";
    for ( graph::Graph::EdgeListIter
           iru  = energy_graph.get_node( i )->edge_list_begin(),
           irue = energy_graph.get_node( i )->edge_list_end();
           iru != irue; ++iru ) {
	EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
	Size const j( edge->get_other_ind( i ) );
	//Size const j( edge->get_first_node_ind() );

	// the pair energies cached in the link
	EnergyMap const & emap( edge->fill_energy_map());
	Real const attr( emap[ fa_atr ] );
	//TR<<"\n"<<j<<": "<<attr<<"\n";
	if (attr < -.2){
	 //std::cout<<input_pose.pdb_info()->number(j)<<input_pose.pdb_info()->chain(j)<<":";
	  std::ostringstream tmp2;
	  tmp2<< input_pose.pdb_info()->chain(j);
	  std::string bch = tmp2.str();
	  // create string id to store in set
	  if (bch != chain && rch == chain){
	    std::ostringstream residuestream;
      //std::cout << rch<<" "<<bch<<" "<<chain<<std::endl;
	    TR << "resi "<< input_pose.pdb_info()->number(i)<<" or ";

	    residuestream << input_pose.pdb_info()->chain(i) << input_pose.pdb_info()->number(i);
	    std::string res_id = residuestream.str();
	    interface.insert(res_id);
	    os<<res_id<<std::endl;
      break;
 	  }
	}
    }
  }


	TR << std::endl;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}



