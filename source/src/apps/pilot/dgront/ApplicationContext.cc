// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/*
 * ApplicationContext.cc
 *
 *  Created on: Jan 13, 2009
 *      Author: dgront
 */
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.hh>


#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option_macros.hh>

#include "ApplicationContext.hh"

//Auto Headers
#include <core/pose/annotated_sequence.hh>


static THREAD_LOCAL basic::Tracer tr( "CustomAbinitio" );

OPT_KEY( Boolean, start_native )
OPT_KEY( File, sf )

void protocols::abinitio::ApplicationContext::register_options() {

  OPT( in::file::native );
  OPT( in::file::s );
  OPT( in::file::frag3 );
  OPT( in::file::frag9 );
  OPT( in::file::fasta );
  OPT( out::file::silent );
  OPT( out::nstruct );
  NEW_OPT( start_native, "start from extended structure (instead of native)", false );
  NEW_OPT( sf, "filename for score output", "score.fsc" );
}

namespace protocols {
  namespace abinitio {

    using namespace core;
    using namespace fragment;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    using std::string;
    using utility::vector1;

    ApplicationContext::ApplicationContext() {

      tr.Info << "staring ApplicationContext" << std::endl;
      //-------------- SOME INPUT
      if (option[out::nstruct].user())
        nStruct_ = std::max(1, option[out::nstruct]());
      else
        nStruct_ = 5;

      //-------------- I/O SETUP   ---------------------
      if (option[sf].user()) {
        tr.Info << "sf noticed at command line" << std::endl;
        silentScoreFile_ = new io::silent::SilentFileData;
        silentScoreFile_-> set_filename(std::string(option[sf]()));
      }
      //-------------- READ INITIAL POSE
      // native pose
      if (option[in::file::native].user()) {
        nativePose_ = new core::pose::Pose;
        tr.Info << "in::file::native noticed at command line" << std::endl;
        core::import_pose::pose_from_pdb(*nativePose_, option[in::file::native]());
        pose::set_ss_from_phipsi(*nativePose_);
        sequence_ = nativePose_->sequence();
      }
      // starting pose
      startingPose_ = new core::pose::Pose;
      generateExtendedPose(*startingPose_, sequence_);

      if (option[start_native]()) {
        copyStructure(*startingPose_, *nativePose_);
      } else if (option[in::file::s].user()) {
        core::pose::PoseOP tmp_pose(new core::pose::Pose);
        std::string fn = option[in::file::s](1);
        core::import_pose::pose_from_pdb(*tmp_pose, fn);
        copyStructure(*startingPose_, *tmp_pose);
      }

      //      if (option[in::file::s].user()) {
      //        tr.Info << "in::file::s noticed at command line" << std::endl;
      //        std::string fn = option[in::file::s](1);
      //        core::chemical::ResidueTypeSetCAP rsd_set =
      //            core::chemical::ChemicalManager::get_instance()->residue_type_set(
      //                core::chemical::CENTROID);
      //        core::import_pose::pose_from_pdb(*startingPose_, *rsd_set, fn);
      //      }

      if (option[start_native]()) {
        tr.Info << "start_native noticed at command line" << std::endl;
        ApplicationContext::copyStructure(*nativePose_, *startingPose_);
      }

      //-------------- READ FRAGMENTS ---------------------
      //declare fragments file
      std::string frag_large_file, frag_small_file;

      if (option[in::file::fragA].user())
        frag_large_file = option[in::file::fragA]();
      else
        frag_large_file = option[in::file::frag9]();

      if (option[in::file::fragB].user())
        frag_small_file = option[in::file::fragB]();
      else
        frag_small_file = option[in::file::frag3]();

      // declare 9mer fragments
      fragsetLarge_ = FragmentIO().read(frag_large_file);

      // declare 3mer fragments
      fragsetSmall_ = FragmentIO().read(frag_small_file);
    }

    void ApplicationContext::generateExtendedPose(
        core::pose::Pose &extended_pose, std::string  &sequence)  {

      core::pose::make_pose_from_sequence(extended_pose, sequence,
          *(chemical::ChemicalManager::get_instance()->residue_type_set(
              chemical::CENTROID)));

      // make extended chain
      for (Size pos = 1; pos <= extended_pose.total_residue(); pos++) {
        if (!extended_pose.residue(pos).is_protein())
          continue;
        extended_pose.set_phi(pos, -150);
        extended_pose.set_psi(pos, 150);
        extended_pose.set_omega(pos, 180);
      }
    }

    void ApplicationContext::copyStructure(core::pose::Pose &source,
        core::pose::Pose &destination) {

      Size seg_len = std::min(destination.total_residue(),
          source.total_residue());
      Size protein_len = 0;
      for (Size i = 1; i <= seg_len; ++i) {
        if (destination.residue(i).is_protein()
            && source.residue(i).is_protein()) {
          protein_len++;
        }
      }
      seg_len = protein_len;
      fragment::Frame long_frame(1, seg_len);
      FragData frag(new BBTorsionSRFD, seg_len);
      frag.steal(source, long_frame);
      frag.apply(destination, long_frame);
    }
  }
}
