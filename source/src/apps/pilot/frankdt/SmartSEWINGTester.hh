// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/frankdt/SmartSEWINGTester.hh
/// @brief testing app for SEWING file IO
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_apps_pilot_frankdt_SmartSEWINGTester_hh
#define INCLUDED_apps_pilot_frankdt_SmartSEWINGTester_hh

#include <apps/pilot/frankdt/SmartSEWINGTester.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace apps {
namespace pilot {
namespace frankdt {

///@brief testing app for SEWING file IO
class SmartSEWINGTester : public utility::pointer::ReferenceCount {

public:

	SmartSEWINGTester();
	SmartSEWINGTester(SmartSEWINGTester const & src);

	virtual ~SmartSEWINGTester();

	SmartSEWINGTesterOP
	clone() const;

private:

};


} //apps
} //pilot
} //frankdt



#endif //INCLUDED_apps_pilot_frankdt_SmartSEWINGTester_hh





