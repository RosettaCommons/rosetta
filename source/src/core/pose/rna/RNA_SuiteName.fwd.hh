// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/rna/RNA_SuiteName.fwd.hh
/// @brief RNA suite assignment ported from suitename program
/// @author Fang-Chieh Chou


#ifndef INCLUDED_core_pose_rna_RNA_SuiteName_fwd_HH
#define INCLUDED_core_pose_rna_RNA_SuiteName_fwd_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pose {
namespace rna {

class RNA_SuiteAssignment;
class RNA_SuiteName;
class RNA_SuiteInfo;
typedef utility::pointer::owning_ptr< RNA_SuiteInfo > RNA_SuiteInfoOP;

}
}
}

#endif
