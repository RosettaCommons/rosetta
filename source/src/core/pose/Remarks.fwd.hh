// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data.hh
///
/// @brief
/// @author Sergey Lyskov

#ifndef INCLUDED_core_pose_Remarks_fwd_hh
#define INCLUDED_core_pose_Remarks_fwd_hh


#include <utility/pointer/owning_ptr.hh>

// C++ headers

namespace core {
namespace pose {

class Remarks;

typedef utility::pointer::shared_ptr< Remarks > RemarksOP;

} // namespace pose
} // namespace core
#ifdef USEBOOSTSERIALIZE
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#endif


#endif // INCLUDED_core_pdb_remarks_FWD_HH

