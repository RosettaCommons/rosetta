// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cryst/spacegroup_symmops4.hh
/// @brief Contains helper functions for the symmetry operations.  Split from spacegroups.cc to prevent the
/// 32-bit compilation from running out of memory during compilation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_protocols_cryst_spacegroup_symmops4_hh
#define INCLUDED_protocols_cryst_spacegroup_symmops4_hh

#include <protocols/cryst/CheshireCell.hh>

// Rosetta core includes:
#include <core/types.hh>
#include <core/kinematics/RT.fwd.hh>

// Rosetta numeric includes:
#include <numeric/xyzVector.hh>

namespace protocols {
namespace cryst {

void get_symmops_P6422( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P6322( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P6mm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P6cc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P63cm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P63mc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pminus6m2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pminus6c2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pminus62m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pminus62c( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P6slashmmm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P6slashmcc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P63slashmcm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P63slashmmc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P23( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_F23( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_I23( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P213( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_I213( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pmminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pnminus3__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Fmminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Fdminus3__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Imminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Paminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Iaminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P432( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P4232( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_F432( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_F4132( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_I432( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P4332( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P4132( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_I4132( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pminus43m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Fminus43m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Iminus43m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pminus43n( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Fminus43c( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Iminus43d( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pmminus3m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pnminus3n__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pmminus3n( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pnminus3m__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Fmminus3m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Fmminus3c( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Fdminus3m__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Fdminus3c__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Imminus3m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Iaminus3d( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_B11m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );


}
}

#endif
