// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <protocols/moves/mc_convergence_checks/HierarchicalLevel.hh>
#include <protocols/moves/mc_convergence_checks/HierarchicalLevel.fwd.hh>

#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

#include <iostream>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::toolbox;
using namespace core::io::silent;
using namespace protocols::moves::mc_convergence_checks;

core::Size
first_zero_pos( Address& addr ) {
  core::Size first_zero_pos;
  for( first_zero_pos = 1; first_zero_pos <= addr.size(); first_zero_pos++ ) {
    if( addr[ first_zero_pos ] == 0 ) break;
  }
  return first_zero_pos;
}

void
assign_tag( Address& addr, std::string & newtag, HierarchicalLevelOP hlevel ) {
  std::ostringstream q;
  Address tmp_addr = addr;
  for( core::Size ii = 1; ii <= tmp_addr.size(); ii++ ) {
    if( tmp_addr[ ii ] == 0 ) {
      q << hlevel->pool_size( tmp_addr, ii - 1 ) + 1 << ".";
      tmp_addr[ ii ] = hlevel->pool_size( tmp_addr, ii - 1 ) + 1;
      //TR.Debug << " at level " << ii << " assigning " << ( hlevel_.level_size( ii ) + 1 ) << std::endl;
    } else {
      q << tmp_addr[ ii ] << ".";
    }
  }
  if( addr[ addr.size() ] == 0 ) {
    q << 1;
  }else {
    q << hlevel->pool_size( addr, hlevel->nlevels() ) + 1;
  }
  newtag = "c." + q.str();
  //std::cout << " newtag: " << newtag << " address: ";
  //for( core::Size ii = 1; ii <= tmp_addr.size(); ii++ ) {
  //  std::cout << tmp_addr[ ii ] << " ";
  //}
  //std::cout << std::endl;

}

int main(int argc, char *argv[])
{
    try {
  devel::init(argc, argv);

  std::string silentin = option[ in::file::silent ]()[1];
  SilentFileData sfd;
  sfd.read_file( silentin );


  core::Size K_level = option[ cluster::K_level ]();
  utility::vector1< core::Real > level_radii = option[ cluster::K_radius ]();

  if( level_radii.size() != K_level ) {
    utility_exit_with_message("you cannot specify a different number of levels and radii! ");
  }

  utility::vector1< std::string > tags = sfd.tags();

  HierarchicalLevelOP hlevel = new HierarchicalLevel( K_level, option[ in::file::silent ]()[ 1 ] );
  double progress = 0;
  time_t t_start, t_finish;
  (void) time(&t_start);
  for( core::Size ii = 1; ii <= tags.size(); ii++ ) {
    if( (double)((double)ii/(double)tags.size()) > (progress + 0.01) ) {
      (void) time(&t_finish);
      std::cout << (double)((double)ii/(double)tags.size()) << " done... finished " << ii << " of " << tags.size() << " time: " << (t_finish-t_start) << std::endl;
      progress += 0.01;
    } else {
      //std::cout << ii <<  " of " << tags.size() << std::endl;
    }
    SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_in();
    ss = sfd[ tags[ ii ] ];
    Address addr( hlevel->nlevels(), 0 );
    utility::vector1< core::Real > best_rmsd;
    std::string best_decoy;
    hlevel->evaluate( *ss, best_decoy, best_rmsd, addr );
    core::Size new_level = first_zero_pos(addr);
    std::string new_tag;
    assign_tag( addr, new_tag, hlevel );
    //std::cout << "tag: " << ss->decoy_tag();
    ss->decoy_tag( new_tag );
    //std::cout << "new_tag: " << ss->decoy_tag() << " address: ";
    //for( core::Size ii = 1; ii <= addr.size(); ii++ ) {
    //  std::cout << addr[ ii ] << " ";
    //}
    //std::cout << std::endl;
    hlevel->add_new( *ss, best_decoy, addr, true, new_level );
  }
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}

