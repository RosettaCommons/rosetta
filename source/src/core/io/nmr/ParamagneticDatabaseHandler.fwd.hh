// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/ParamagneticDatabaseHandler.fwd.hh
/// @brief   Forward declaration for ParamagneticDatabaseHandler
/// @details last Modified: 05/11/16
/// @author  Georg Kuenze georg.kuenze@vanderbilt.edu

#ifndef INCLUDED_core_io_nmr_ParamagneticDatabaseHandler_FWD_HH
#define INCLUDED_core_io_nmr_ParamagneticDatabaseHandler_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace io {
namespace nmr {

/// @brief  A singleton class for handling static paramagnetic ion properties.
class ParamagneticDatabaseHandler;

typedef utility::pointer::shared_ptr< ParamagneticDatabaseHandler > ParamagneticDatabaseHandlerOP;
typedef utility::pointer::shared_ptr< ParamagneticDatabaseHandler const > ParamagneticDatabaseHandlerCOP;

} // namespace nmr
} // namespace io
} // namespace core


#endif // INCLUDED_core_io_nmr_ParamagneticDatabaseHandler_FWD_HH
