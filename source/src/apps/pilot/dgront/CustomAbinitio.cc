// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/*
 * CustomAbinitio.cc
 *
 *  Created on: Jan 12, 2009
 *      Author: dgront
 */

// Project Headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragmentIO.hh>

#include "CustomAbinitio.hh"
#include "ApplicationContext.hh"
#include "ApplicationContext.cc"

OPT_1GRP_KEY( IntegerVector, dgront, declare_domain )

void protocols::abinitio::CustomAbinitio::register_options() {
  NEW_OPT( dgront::declare_domain, "declare a protein fragment (possibly a domain) that will be treated as a rigid body",0);
}

namespace protocols {
  namespace abinitio {

    CustomAbinitio::CustomAbinitio(ApplicationContext &context,
        core::fragment::FragSetCOP fragset_small,
        core::fragment::FragSetCOP fragset_large,
        core::kinematics::MoveMapCOP movemap) :
      ClassicAbinitio(fragset_small, fragset_large, movemap) {
      context_ = context;
      setDefaults();
    }

    void CustomAbinitio::setDefaults() {

      init(*(context_.startingPose_));

      bSkipStage1_ = false;
      bSkipStage2_ = false;
      bSkipStage3_ = true;
      bSkipStage4_ = true;
    }

    void CustomAbinitio::run() {

      for (int i = 1; i <= context_.nStruct_; i++) {
        apply(*(context_.startingPose_));
        io::pdb::old_dump_pdb(*(context_.startingPose_), "resulting_pose.pdb");

        io::silent::SilentFileData outsfd;
        outsfd.strict_column_mode(true);

        io::silent::SilentStructOP pss =
            io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
        if (context_.silentScoreFile_) {
          context_.silentScoreFile_ -> write_silent_struct(*pss,
              context_.silentScoreFile_->filename(), true);
        }
      }
    }
  }
}

int main(int argc, char * argv[]) {
try {
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using std::string;
  using utility::vector1;
  using namespace protocols::abinitio;

  ApplicationContext::register_options();
  ClassicAbinitio::register_options();
  CustomAbinitio::register_options();
  devel::init(argc, argv);

  ApplicationContext* context = new ApplicationContext();

  kinematics::MoveMapOP movemap = new kinematics::MoveMap;
  movemap->set_bb(true);
  //-------------- SET UP DOMAINS -------------
  if (option[dgront::declare_domain].user()) {
    tr.Info << "dgront::declare_domain noticed at command line" << std::endl;
    tr.Info << "The fixed residues are:" << std::endl;
    utility::vector1<int> const& domainBoundaries(
        option[dgront::declare_domain]());
    for (Size i = 1; i <= domainBoundaries.size() / 2; i++) {
      for (Size iRes = domainBoundaries[i * 2 - 1]; iRes <= domainBoundaries[i
          * 2]; iRes++) {
        movemap->set_bb(iRes, false);
        tr.Info << context->startingPose_->residue(iRes).seqpos()<<" "<< context->startingPose_->residue(iRes).name()<< std::endl;
      }
    }
  }

  //-------------- GO! ---------------------
  tr.Info << "Creating the sampler" << std::endl;
  protocols::abinitio::CustomAbinitio abinitio(*context,
      context->fragsetSmall_, context->fragsetLarge_, movemap);

  tr.Info << "starting the calculations" << std::endl;
  abinitio.run();
} catch ( utility::excn::EXCN_Base const & e ) {
                          std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                              }
  return 0;
}
