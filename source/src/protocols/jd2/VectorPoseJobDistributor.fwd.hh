// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/VectorPoseJobDistributor.fwd.hh
/// @brief  Job distributor for running RECON MSD under MPI
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_protocols_jd2_VectorPoseJobDistributor_fwd_hh
#define INCLUDED_protocols_jd2_VectorPoseJobDistributor_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd2 {

class VectorPoseJobDistributor;
typedef utility::pointer::shared_ptr< VectorPoseJobDistributor > VectorPoseJobDistributorOP;
typedef utility::pointer::shared_ptr< VectorPoseJobDistributor const > VectorPoseJobDistributorCOP;

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_VectorPoseJobDistributor_fwd_hh
