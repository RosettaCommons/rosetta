// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/ParaIon.fwd.hh
/// @brief   forward declaration for ParaIon
/// @details last Modified: 05/11/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_io_nmr_ParaIon_FWD_HH
#define INCLUDED_core_io_nmr_ParaIon_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace io {
namespace nmr {

class ParaIon;

typedef utility::pointer::shared_ptr< ParaIon > ParaIonOP;
typedef utility::pointer::shared_ptr< ParaIon const > ParaIonCOP;
typedef utility::pointer::weak_ptr< ParaIon > ParaIonAP;
typedef utility::pointer::weak_ptr< ParaIon const > ParaIonCAP;


} // namespace nmr
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_nmr_ParaIon_FWD_HH
