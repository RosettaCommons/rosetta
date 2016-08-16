// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef CUSTOMABINITIO_HH_
#define CUSTOMABINITIO_HH_

// Unit Headers
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include "ApplicationContext.hh"

namespace protocols {
  namespace abinitio {

    class CustomAbinitio: public ClassicAbinitio {
      public:
        ApplicationContext context_;
        CustomAbinitio(ApplicationContext &context,core::fragment::FragSetCOP fragset_small,
            core::fragment::FragSetCOP fragset_large,
            core::kinematics::MoveMapCOP movemap);
        void  setDefaults();
        void run();

        static void register_options();
    };
  }
}
#endif /* CUSTOMABINITIO_HH_ */
