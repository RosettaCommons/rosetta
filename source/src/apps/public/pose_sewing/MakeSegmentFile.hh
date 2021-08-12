// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/frankdt/MakeSegmentFile.hh
/// @brief an app to make segment files
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_apps_pilot_frankdt_MakeSegmentFile_hh
#define INCLUDED_apps_pilot_frankdt_MakeSegmentFile_hh

#include <apps/public/pose_sewing/MakeSegmentFile.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace apps {
namespace pilot {
namespace frankdt {

/// @brief an app to make segment files
class MakeSegmentFile : public utility::VirtualBase {

public:

	MakeSegmentFile();
	MakeSegmentFile(MakeSegmentFile const & src);

	~MakeSegmentFile() override;

	MakeSegmentFileOP
	clone() const;

private:

};


} //apps
} //pilot
} //frankdt



#endif //INCLUDED_apps_pilot_frankdt_MakeSegmentFile_hh





