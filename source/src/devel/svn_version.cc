// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/svn_version.cc
///
/// @brief
/// @author Ian W. Davis
/// @author Andrew Leaver-Fay
/// @author Sergey Lyskov

#include <core/svn_version.hh>

#include <utility/version.hh>


namespace devel {

std::string rosetta_svn_version() { return utility::Version::commit_id() + ':' + utility::Version::commit() + ' ' + utility::Version::date(); }
std::string rosetta_svn_url() { return utility::Version::url(); }

class VersionRegistrator
{
public:
	VersionRegistrator() {
		core::set_svn_version_and_url( rosetta_svn_version(), rosetta_svn_url() );
	}
};

// There should only ever be one instance of this class
// so that core::set_svn_version_and_url is called only once
VersionRegistrator vr;

void
register_version_with_core() {
	// oh -- there's nothing in this function.  But
	// forcing devel::init to call this function ensures
	// that the vr variable in this file gets instantiated
}

} // namespace devel
