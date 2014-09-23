// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/ContactMap.fwd.hh
/// @brief  ContactMap forward declarations header
/// @author Joerg Schaarschmidt

#ifndef INCLUDED_protocols_contact_map_ContactMap_fwd_hh
#define INCLUDED_protocols_contact_map_ContactMap_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols{
namespace contact_map{

//Forwards and OP/COP/AP typedefs
class ContactMap;
class Contact;
class ContactPartner;

typedef utility::pointer::shared_ptr< ContactMap > ContactMapOP;
typedef utility::pointer::shared_ptr< ContactMap const > ContactMapCOP;
typedef utility::pointer::weak_ptr< ContactPartner > ContactPartnerAP;

}//moves
}//protocols

#endif //INCLUDED_protocols_moves_ContactMap_FWD_HH
