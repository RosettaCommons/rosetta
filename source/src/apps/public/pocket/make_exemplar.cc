// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>

#include <protocols/jd2/util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>



// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <utility/io/mpistream.hh>
#include <core/kinematics/MoveMap.hh>

//Protocol Headers
#include <protocols/pockets/PocketGrid.hh>
//#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;
using namespace conformation;
using namespace core::pose::datacache;
using namespace core::id;
using namespace protocols::rigid;
using namespace protocols::simple_moves;
using namespace protocols;
using namespace protocols::moves;


OPT_KEY( String, central_relax_pdb_num )

static thread_local basic::Tracer TR( "apps.public.make_exemplar.main", basic::t_debug );

//This Mover creates an exemplar
class ExemplarMover : public moves::Mover {

public:
  ExemplarMover();

  ~ExemplarMover();

  virtual MoverOP clone() const;
  virtual MoverOP fresh_instance() const;

  virtual void apply( Pose & pose );
  virtual std::string get_name() const;
  virtual void test_move( Pose & pose )
  {
    apply(pose);
  }
};

ExemplarMover::ExemplarMover() :
  Mover( "benchmark" )
{

}

ExemplarMover::~ExemplarMover() {}

MoverOP ExemplarMover::clone() const {
  return MoverOP( new ExemplarMover( *this ) );
}

MoverOP ExemplarMover::fresh_instance() const {
  return MoverOP( new ExemplarMover );
}

void ExemplarMover::apply( Pose & pose ) {
  using namespace pose;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core;

  std::string const resid_c = option[central_relax_pdb_num];

  //validate that the exemplar_target_pdb_num is properly formatted
  if (resid_c.length()){
    std::string resid_tag = resid_c;
    while ( true) {
      std::size_t fpos( resid_tag.find(','));
      if ( fpos == std::string::npos ) break;
        resid_tag[fpos]='-';
    }
    std::string out_exfname = get_current_tag().substr(0,get_current_tag().find_last_of("_")) +".pdb."+ resid_tag + ".exemplar.pdb";

    std::vector< conformation::ResidueOP > residues = protocols::pockets::PocketGrid::getRelaxResidues(pose, resid_c);

    protocols::pockets::PocketGrid comparison_pg( residues );
    comparison_pg.zeroAngle();
    comparison_pg.autoexpanding_pocket_eval( residues, pose ) ;
    comparison_pg.dumpExemplarToFile( out_exfname.c_str() );
  }

}

std::string ExemplarMover::get_name() const {
  return "ExemplarMover";
}


/// General testing code
int main( int argc, char * argv [] ) {

	try{
  using namespace core;
  using namespace protocols::moves;
  using namespace scoring;
  using namespace basic::options;
  using namespace protocols::jobdist;
  using namespace basic::options::OptionKeys;
  using protocols::moves::MoverOP;

	NEW_OPT( central_relax_pdb_num, "target residue(s)", "-1");

	// APL NOTE: Tracers cannot be written to before devel::init gets called. TR << "Calling init" << std::endl;
	//initializes Rosetta functions
  jd2::register_options();
  option.add_relevant( OptionKeys::in::file::fullatom );
  option.add_relevant( OptionKeys::in::file::movemap );
	devel::init(argc, argv);

  MoverOP protocol( new ExemplarMover() );
  protocols::jd2::JobDistributor::get_instance()->go( protocol,  jd2::JobOutputterOP( new jd2::NoOutputJobOutputter ) );

  }
	catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
    return -1;
  }
	return 0;
}

