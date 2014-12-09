// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef APPLICATIONCONTEXT_HH_
#define APPLICATIONCONTEXT_HH_

namespace protocols {
  namespace abinitio {

    class ApplicationContext {
      public:
        int nStruct_;
        std::string sequence_;
        core::io::silent::SilentFileDataOP silentScoreFile_;
        core::pose::PoseOP nativePose_;
        core::pose::PoseOP startingPose_;
        core::fragment::FragSetOP fragsetLarge_;
        core::fragment::FragSetOP fragsetSmall_;
        ApplicationContext();

        static void generateExtendedPose(core::pose::Pose &extended_pose,
            std::string &sequence);
        static void copyStructure(core::pose::Pose & source,
            core::pose::Pose & destination);
        static void register_options();
    };
  }
}

#endif /* APPLICATIONCONTEXT_HH_ */
