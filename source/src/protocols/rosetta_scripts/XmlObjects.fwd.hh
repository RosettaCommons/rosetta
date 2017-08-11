// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rosetta_scripts/XmlObjects.fwd.hh
/// @brief Class to load objects from xml
/// @details Construct this class with a xml string and a pose, then use the
///   get_* methods to extract objects by name.
/// @author Brian Coventry


#ifndef INCLUDED_protocols_rosetta_scripts_XmlObjects_fwd_hh
#define INCLUDED_protocols_rosetta_scripts_XmlObjects_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace rosetta_scripts {

// Forward
class XmlObjects;

// Types
typedef  utility::pointer::shared_ptr< XmlObjects >  XmlObjectsOP;
typedef  utility::pointer::shared_ptr< XmlObjects const >  XmlObjectsCOP;

typedef  utility::pointer::weak_ptr< XmlObjects >  XmlObjectsAP;
typedef  utility::pointer::weak_ptr< XmlObjects const >  XmlObjectsCAP;


} // rosetta_scripts
} // protocols

#endif
