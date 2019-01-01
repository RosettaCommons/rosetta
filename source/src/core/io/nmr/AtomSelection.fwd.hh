// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/AtomSelection.fwd.hh
/// @brief   forward declaration for AromSelection
/// @details last Modified: 06/08/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_io_nmr_AtomSelection_FWD_HH
#define INCLUDED_core_io_nmr_AtomSelection_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace io {
namespace nmr {

class AtomSelection;

typedef utility::pointer::shared_ptr< AtomSelection > AtomSelectionOP;
typedef utility::pointer::shared_ptr< AtomSelection const > AtomSelectionCOP;
typedef utility::pointer::weak_ptr< AtomSelection > AtomSelectionAP;
typedef utility::pointer::weak_ptr< AtomSelection const > AtomSelectionCAP;


} // namespace nmr
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_nmr_AtomSelection_FWD_HH
