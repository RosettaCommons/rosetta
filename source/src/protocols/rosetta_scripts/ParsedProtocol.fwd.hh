// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Forward declarations for ParsedProtocol
/// @author Sarel Fleishman


#ifndef INCLUDED_protocols_rosetta_scripts_ParsedProtocol_FWD_HH
#define INCLUDED_protocols_rosetta_scripts_ParsedProtocol_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace rosetta_scripts {

class ParsedProtocol;

typedef utility::pointer::shared_ptr< ParsedProtocol > ParsedProtocolOP;
typedef utility::pointer::shared_ptr< ParsedProtocol const > ParsedProtocolCOP;
typedef utility::pointer::weak_ptr< ParsedProtocol > ParsedProtocolAP;
typedef utility::pointer::weak_ptr< ParsedProtocol const > ParsedProtocolCAP;

} // namespace protocols
} // namespace rosetta_scripts

#endif //INCLUDED_protocols_rosetta_scripts_ParsedProtocol_FWD_HH
