// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cryst/spacegroup_symmops1.hh
/// @brief Contains helper functions for the symmetry operations.  Split from spacegroups.cc to prevent the
/// 32-bit compilation from running out of memory during compilation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_protocols_cryst_spacegroup_symmops1_hh
#define INCLUDED_protocols_cryst_spacegroup_symmops1_hh

#include <protocols/cryst/CheshireCell.hh>

// Rosetta core includes:
#include <core/types.hh>
#include <core/kinematics/RT.fwd.hh>

// Rosetta numeric includes:
#include <numeric/xyzVector.hh>

namespace protocols {
namespace cryst {

void get_symmops_P1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pminus1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P121( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P1211( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_C121( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P1m1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P1c1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_C1m1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_C1c1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P12slashm1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P121slashm1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_C12slashm1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P12slashc1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P121slashc1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_C12slashc1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P222( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P2221( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P21212( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_P212121( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_C2221( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_C222( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_F222( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_I222( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_I212121( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pmm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pmc21( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pcc2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pma2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pca21( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pnc2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pmn21( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pba2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pna21( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pnn2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Cmm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Cmc21( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Ccc2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Amm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Abm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Ama2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Aba2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Fmm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Fdd2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Imm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Iba2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Ima2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pmmm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pnnn__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pccm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pban__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pmma( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pnna( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pmna( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pcca( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pbam( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pccn( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pbcm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pnnm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pmmn__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );
void get_symmops_Pbcn( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc );

}
}

#endif
