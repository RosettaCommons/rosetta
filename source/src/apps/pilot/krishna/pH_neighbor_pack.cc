// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/krishna/pH_protocol.cc
/// @Calculates the pKa value for a residue

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/pHEnergy.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <devel/init.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>

#include <basic/Tracer.hh>
#include <string>
#include <iostream>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

//old JD and relevant headers
//#include <protocols/jobdist/standard_mains.hh>
//#include <core/io/raw_data/ScoreMap.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.krishna.PhProtocol" );

class PhProtocol : public protocols::moves::Mover {

public:
  PhProtocol()
  {
    // curr_pH_ = 0.1; //start value for pH
  }

  virtual ~PhProtocol(){};

  virtual
  void
  apply( core::pose::Pose & pose ){

    using namespace core;
    using namespace basic::options;
    using namespace chemical;
    using namespace conformation;

    Real shift_pH = 1.0, pka_value = 1.0, ipka = 1.0;

    //setting up residue details
    std::string pdb_file_name = option[ OptionKeys::in::file::s ]()[1];
    int pdb_res_no = option[ OptionKeys::pH::calc_pka::pka_for_resno](); //int deliberately used for the following condition
    Real repack_rad = option[ OptionKeys::pH::calc_pka::pka_rad]();

    if ( pdb_res_no < 0 ){
      TR << "PLEASE ENTER RESIDUE NO FOR WHICH PKA IS TO BE CALCULATED AND TRY AGAIN" << std::endl;
      return;
    }

    std::string pdb_chain_no = option[ OptionKeys::pH::calc_pka::pka_for_chainno]();
    Size res_no = pose.pdb_info()->pdb2pose( pdb_chain_no[0], pdb_res_no );
    std::string res_name = pose.residue( res_no ).name();
    std::string res_name3 = pose.residue( res_no ).name3();

    //ipKa details
    switch ( pose.residue( res_no ).type().aa() ){
      case aa_asp: ipka = 4.0; break;
      case aa_glu: ipka = 4.4; break;
      case aa_his: ipka = 6.3; break;
      case aa_lys: ipka = 10.4; break;
      case aa_tyr: ipka = 10.0; break;
      default: ipka = 0.0; break;
    }

    TR << " RES NO BEING PACKED IS " << pdb_res_no << " FROM CHAIN " << pdb_chain_no << std::endl;

    pose::Pose old_pose = pose;  //initialize the pose

    for ( Real curr_pH = 0.0; curr_pH <= 13.0; curr_pH += 1.0 ){

      pose::Pose curr_pose = pose;

      TR << " CURRENTLY SIMULATING AT pH " << curr_pH << std::endl;

      pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( curr_pose ));

      conformation::Residue const & res1 = curr_pose.residue( res_no );
    
      for ( Size ii = 1; ii <= curr_pose.size() ; ++ii) {

        conformation::Residue const & res2 = pose.residue( ii );

        if ( ii == res_no ) {
          task->nonconst_residue_task( ii ).restrict_to_repacking();
        }

//        else if ( curr_pose.residue(res_no).xyz("CA").distance(curr_pose.residue(ii).xyz("CA")) < repack_rad ){

        else if ( res1.xyz( res1.nbr_atom() ).distance( res2.xyz( res2.nbr_atom()) ) < repack_rad ) {
          task->nonconst_residue_task( ii ).restrict_to_repacking();
        }

        else {
          task->nonconst_residue_task( ii ).prevent_repacking();
        }

      }

      std::cout << *task << std::endl;

      scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );
      scoring::methods::pHEnergy::set_pH ( curr_pH );

      protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover( score_fxn, task ) );

      pack_mover->apply(curr_pose);

      (*score_fxn)(curr_pose);
//      (*score_fxn).show(curr_pose);

      if ( ( old_pose.residue( res_no ).name() != curr_pose.residue( res_no ).name() ) && ( curr_pH != 0.0 ) )
        break;

      old_pose = curr_pose;
      shift_pH = curr_pH;

    }

    for ( Real curr_pH = shift_pH; curr_pH < shift_pH + 1.0; curr_pH += 0.1 ){

      pose::Pose curr_pose = old_pose; //look at smaller intervals of pH taking the previous pose

      pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( curr_pose ));

      conformation::Residue const & res1 = curr_pose.residue( res_no );

      for ( Size ii = 1; ii <= curr_pose.size() ; ++ii) {

        conformation::Residue const & res2 = pose.residue( ii );

        if ( ii == res_no ) {
          task->nonconst_residue_task( ii ).restrict_to_repacking();
        }

//        else if ( curr_pose.residue(res_no).xyz("CA").distance(curr_pose.residue(ii).xyz("CA")) < repack_rad ){

        else if ( res1.xyz( res1.nbr_atom() ).distance( res2.xyz( res2.nbr_atom()) ) < repack_rad ) {
          task->nonconst_residue_task( ii ).restrict_to_repacking();
        }

        else {
          task->nonconst_residue_task( ii ).prevent_repacking();
        }

      }

/*
      for ( Size ii = 1; ii <= curr_pose.size() ; ++ii) {

        if ( nbr_atom.distance( curr_pose.residue(ii).nbr_atom_xyz() ) < 5.0 || ii == res_no ) {
          task->nonconst_residue_task( ii ).restrict_to_repacking();
        }
        else {
          task->nonconst_residue_task( ii ).prevent_repacking();
        }

      }
*/
      scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );
      scoring::methods::pHEnergy::set_pH ( curr_pH );

      protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover( score_fxn, task ) );

      TR << " CURRENTLY SIMULATING AT pH " << curr_pH << std::endl;

      pack_mover->apply(curr_pose);

      (*score_fxn)(curr_pose);

//      (*score_fxn).show(curr_pose);

      if ( old_pose.residue( res_no ).name() != curr_pose.residue( res_no ).name() ){
//	score_map_["pKa"] = curr_pH;
//	core::io::raw_data::ScoreMap::nonzero_energies( score_map_, score_fxn, curr_pose );
        break;
      }

      old_pose = curr_pose;
      pka_value = curr_pH;

    }

    if ( pka_value > 13.9 ){
      TR << " ENTERED THE LOOP " << std::endl;
      pka_value = 14.0;   //If it doesn't titrate return 14
    }

    TR << "PKA FOR\t" << res_name3 << "\t" << pdb_chain_no << "\t" << pdb_res_no << "\t" << pdb_file_name << "\t" << ipka << "\t" << pka_value << std::endl;

  }

      std::string get_name() const { return "PhProtocol"; }

  virtual
  protocols::moves::MoverOP
  fresh_instance() const {
    return new PhProtocol;
  }

  virtual
  bool
  reinitialize_for_each_job() const { return false; }

  virtual
  bool
  reinitialize_for_new_input() const { return false; }

private:

  //core::Real curr_pH_;
//  std::map < std::string, core::Real > score_map_;

};

typedef utility::pointer::owning_ptr< PhProtocol > PhProtocolOP;


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
//  using namespace protocols;
//  using namespace protocols::jobdist;
//  using namespace protocols::moves;

  devel::init(argc, argv);

  PhProtocolOP pH_test(new PhProtocol);

//  protocols::jobdist::main_plain_mover( *pH_test );

  protocols::jd2::JobDistributor::get_instance()->go(pH_test);

}


