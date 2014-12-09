// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/ShortestPathInFoldTree.fwd.hh
/// @brief  kinematics::ShortestPathInFoldTree forward declarations header
/// @author Oliver Lange


#ifndef INCLUDED_protocols_jd2_archive_ArchiveBase_fwd_hh
#define INCLUDED_protocols_jd2_archive_ArchiveBase_fwd_hh


// Utility headers
// AUTO-REMOVED #include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace jd2 {
namespace archive {
// Forward
class AbstractArchiveBase;
// Types
typedef  utility::pointer::shared_ptr< AbstractArchiveBase >  AbstractArchiveBaseOP;
typedef  utility::pointer::shared_ptr< AbstractArchiveBase const >  AbstractArchiveBaseCOP;

class ArchiveBase;
// Types
typedef  utility::pointer::shared_ptr< ArchiveBase >  ArchiveBaseOP;
typedef  utility::pointer::shared_ptr< ArchiveBase const >  ArchiveBaseCOP;

}
} // namespace kinematics
} // namespace core

#endif
