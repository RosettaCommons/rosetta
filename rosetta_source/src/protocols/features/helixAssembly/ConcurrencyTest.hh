// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ConcurrencyTest.hh
///
/// @brief
/// @author tim

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <protocols/features/helixAssembly/ConcurrencyTest.fwd.hh>

//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

//External
#include <boost/uuid/uuid.hpp>

//Devel
#include <protocols/features/helixAssembly/HelixBundleFeatures.hh>
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility and basic
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//C++
#include <string>
#include <math.h>

//External Headers
#include <cppdb/frontend.h>

//Basic
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/helixAssembly.OptionKeys.gen.hh>

#ifndef INCLUDED_ConcurrencyTest_hh
#define INCLUDED_ConcurrencyTest_hh

namespace protocols {
namespace features {
namespace helixAssembly {

class ConcurrencyTest : public protocols::features::FeaturesReporter
{

public:
	
	ConcurrencyTest(){}
	
	virtual
	std::string
	type_name() const  {
		return "HelixBundleFeatures";
	}
	
	///@brief return sql statements that sets up the appropriate tables
	///to contain the features.
	virtual
	std::string
	schema() const;
	
	///@brief collect all the feature data for the pose
	virtual
	core::Size
	report_features(
					core::pose::Pose const & pose,
					utility::vector1<bool> const & relevant_residues,
					boost::uuids::uuid struct_id,
					utility::sql_database::sessionOP db_session
					);
};
	
}
}
}
#endif
