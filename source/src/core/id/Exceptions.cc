// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/Job.hh
/// @brief  header file for ThreadingJob classes, part of August 2008 job distributor as planned at RosettaCon08.  This file is responsible for three ideas: "inner" jobs, "outer" jobs (with which the job distributor works) and job container (currently just typdefed in the .fwd.hh)
/// @author Steven Lewis smlewi@gmail.com


#include <core/id/Exceptions.hh>
#include <utility/excn/Exceptions.hh>
#include <core/id/NamedAtomID.hh>

namespace core {
namespace id {

EXCN_AtomNotFound::EXCN_AtomNotFound(char const *file, int line, NamedAtomID const& id ) :
	Exception(file, line, "Atom '"+id.to_string()+"' not found"), id_( id )
{
}

}
}
