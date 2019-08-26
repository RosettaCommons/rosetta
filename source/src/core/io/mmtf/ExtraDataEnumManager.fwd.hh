// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmtf/ExtraDataEnumManager.fwd.hh
/// @brief Enum string/enum functions for pose extra data we will be storing/retrieving.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_io_mmtf_ExtraDataEnumManager_fwd_hh
#define INCLUDED_core_io_mmtf_ExtraDataEnumManager_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace io {
namespace mmtf {

class ExtraDataEnumManager;

typedef utility::pointer::shared_ptr< ExtraDataEnumManager > ExtraDataEnumManagerOP;
typedef utility::pointer::shared_ptr< ExtraDataEnumManager const > ExtraDataEnumManagerCOP;

} //core
} //io
} //mmtf

#endif //INCLUDED_core_io_mmtf_ExtraDataEnumManager_fwd_hh
