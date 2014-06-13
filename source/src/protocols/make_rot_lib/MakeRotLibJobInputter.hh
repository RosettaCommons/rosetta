// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/make_rot_lib/MakeRotLibJobInputter.hh
/// @brief  Header file for MakeRotLibJobInputter class
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )


#ifndef INCLUDED_protocols_make_rot_lib_MakeRotLibJobInputter_hh
#define INCLUDED_protocols_make_rot_lib_MakeRotLibJobInputter_hh

//unit headers
#include <protocols/make_rot_lib/MakeRotLibJobInputter.fwd.hh>
#include <protocols/make_rot_lib/MakeRotLibOptionsData.fwd.hh>

#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

//utility headers

namespace protocols {
namespace make_rot_lib {

///@details JobInputter that creats jobs based on MakeRotLib option files.
class MakeRotLibJobInputter : public protocols::jd2::JobInputter
{
public:

  MakeRotLibJobInputter();

  virtual ~MakeRotLibJobInputter();

  ///@brief this function is responsible for filling the pose reference with the pose indicated by the job
  virtual void pose_from_job( core::pose::Pose & pose, jd2::JobOP job );

  ///@brief Determines what jobs exist from the make_rot_lib options file
  virtual void fill_jobs( jd2::Jobs & jobs );

  /// @brief Return the type of input source that the MakeRotLibJobInputter is currently using.
  virtual jd2::JobInputterInputSource::Enum input_source() const;

private:

	MakeRotLibOptionsDataOP mrlod_;

}; // MakeRotLibJobInputter

} // namespace make_rot_lib
} // namespace protocols

#endif //INCLUDED_protocols_make_rot_lib_MakeRotLibJobInputter_HH
