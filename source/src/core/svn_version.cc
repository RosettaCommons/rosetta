// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/svn_version.cc
///
/// @brief
/// @author Ian W. Davis
/// @author Andrew Leaver-Fay

#include <core/svn_version.hh>
#include <utility/exit.hh>
#include <string>

namespace core {

class SVNVersion {
public:
	// Singelton accessor pattern
	static
	SVNVersion *
	get_instance() {
		if ( ! instance_ ) {
			instance_ = new SVNVersion;
		}
		return instance_;
	}

	/// Initialize singelton data.  Should be called at most once.
	void set_svn_version_and_url(
		std::string const & version_in,
		std::string const & url_in
	)
	{
		if ( version_ != "" || url_ != "" ) {
			utility_exit_with_message( "SVNVersion::set_svn_version_and_url called multiple times.  Should only be called once.  Was devel::init() called twice?" );
		}
		if ( version_in == "" ) {
			utility_exit_with_message( "SVNVersion::set_svn_version_and_url called with a blank version" );
		}
		if ( url_in == "" ) {
			utility_exit_with_message( "SVNVersion::set_svn_version_and_url called with a blank url" );
		}
		version_ = version_in;
		url_ = url_in;
	}

	std::string const & version() const { return version_; }
	std::string const & url() const { return url_; }

private :

	static SVNVersion * instance_;

	std::string version_;
	std::string url_;

	SVNVersion() :
		version_( "" ),
		url_( "" )
	{}


};

SVNVersion * SVNVersion::instance_( 0 );

void
set_svn_version_and_url(
	std::string const & version,
	std::string const & url
)
{
	SVNVersion::get_instance()->set_svn_version_and_url( version, url );
}

std::string minirosetta_svn_version() { return SVNVersion::get_instance()->version(); }
std::string minirosetta_svn_url() { return SVNVersion::get_instance()->url(); }

} // namespace core
