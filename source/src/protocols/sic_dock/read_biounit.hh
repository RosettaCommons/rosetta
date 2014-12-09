// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_sic_dock_read_biounit_hh
#define INCLUDED_protocols_sic_dock_read_biounit_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/sic_dock/types.hh>
#include <map>

namespace protocols {
namespace sic_dock {

bool
read_biounit(
	std::string const & fname,
	core::pose::Pose & pose,
	int const & max_res = 9999999,
	bool const & DEBUG = false
);

bool
read_biounit(
	std::string const & fname,
	core::pose::Pose & pose,
	utility::vector1<core::Real> & bfactors,
	utility::vector1<core::Real> & occupancy,
	int const & max_res = 9999999,
	bool const & DEBUG = false
);

bool
read_biounit(
	std::string const & fname,
	core::pose::Pose & pose,
	utility::vector1<core::Real> & bfactors,
	utility::vector1<core::Real> & occupancy,
	utility::vector1<int>  & pdbres,
	std::map<int,char> & pdbchain,
	int& nresmodel1,
	int const & max_res = 9999999,
	bool const & DEBUG = false
);

}
}

#endif
