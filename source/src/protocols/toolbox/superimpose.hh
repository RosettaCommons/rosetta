// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/Job.hh
/// @brief  header file for ThreadingJob classes, part of August 2008 job distributor as planned at RosettaCon08.  This file is responsible for three ideas: "inner" jobs, "outer" jobs (with which the job distributor works) and job container (currently just typdefed in the .fwd.hh)
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_toolbox_superimpose_hh
#define INCLUDED_protocols_toolbox_superimpose_hh

// AUTO-REMOVED #include <ObjexxFCL/FArray3D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <numeric/xyzMatrix.hh>
// AUTO-REMOVED #include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>


namespace protocols {
namespace toolbox {

typedef  numeric::xyzMatrix< core::Real > Matrix;
typedef  numeric::xyzVector< core::Real > Vector;

void CA_superimpose( ObjexxFCL::FArray1_double const& weights, core::pose::Pose const& ref_pose, core::pose::Pose& fit_pose );
void CA_superimpose( core::pose::Pose const& ref_pose, core::pose::Pose& fit_pose );
/*
void superimpose(
     core::Size natoms,
     ObjexxFCL::FArray1_double& weights,
     ObjexxFCL::FArray2_double& ref_coords,
     ObjexxFCL::FArray2_double& coords,
     Matrix &R //returns rotation matrix
);

void superimpose(
     core::Size natoms,
     ObjexxFCL::FArray1_double& weights,
     ObjexxFCL::FArray2_double& ref_coords,
     ObjexxFCL::FArray2_double& coords
); */

void fit_centered_coords(
     core::Size natoms,
     ObjexxFCL::FArray1_double const& weights,
     ObjexxFCL::FArray2_double const& ref_coords,
     ObjexxFCL::FArray2_double& coords,
     Matrix &R
);

void reset_x(
     core::Size n,
     ObjexxFCL::FArray2_double& x,
     ObjexxFCL::FArray1_double const& wts,
     ObjexxFCL::FArray1_double& transvec
);

///@brief write a CA ALA pdb
void dump_as_pdb(
		 std::string filename,
		 core::Size n,
		 ObjexxFCL::FArray2_double& coords,
		 ObjexxFCL::FArray1D_double transvec = ObjexxFCL::FArray1D_double( 3, 0.0 )
);

void fill_CA_coords( core::pose::Pose const& pose, core::Size natoms, ObjexxFCL::FArray2_double& coords );
void fill_CA_coords( core::pose::Pose const& pose, ObjexxFCL::FArray2_double& coords );

}
}

#endif
