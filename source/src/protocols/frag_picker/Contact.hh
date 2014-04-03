// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/Contact.hh
/// @brief  A contact.
/// @author David E. Kim (dekim@u.washington.edu)

#ifndef INCLUDED_protocols_frag_picker_contact_hh
#define INCLUDED_protocols_frag_picker_contact_hh

// unit headers
#include <protocols/frag_picker/Contact.fwd.hh>

// project headers
#include <protocols/frag_picker/ContactTypes.hh>

// type headers
#include <core/types.hh>

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>
#include <cmath>

namespace protocols {
namespace frag_picker {

class Contact: public utility::pointer::ReferenceCount {
public:

	Contact( core::Size i, core::Size j, core::Real dist_squared, ContactType type ) {
		i_ = i;
		j_ = j;
		dist_squared_ = dist_squared;
		type_ = type;
	}

	~Contact(){};

  core::Size & i() {
    return i_;
  }

	core::Size & j() {
		return j_;
	}

	core::Real & dist_squared() {
		return dist_squared_;
	}

	core::Real dist() {
		return std::sqrt(dist_squared_);
	}

	ContactType & type() {
		return type_;
	}

	std::string type_name() {
		return contact_name(type_);
	}

private:
	core::Size i_;
	core::Size j_;
	core::Real dist_squared_;
	ContactType type_;
};


} // namespace frag_picker
} // namespace protocols

#endif // INCLUDED_protocols_frag_picker_contact_HH
