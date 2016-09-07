// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/PosePropertyReporterCreator.hh
/// @brief  Base class for PosePropertyReporters for the load-time factory registration scheme
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_moves_PosePropertyReporterCreator_hh
#define INCLUDED_protocols_moves_PosePropertyReporterCreator_hh

// Unit Headers
#include <protocols/rosetta_scripts/PosePropertyReporter.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <string>

namespace protocols {
namespace rosetta_scripts {

/// @brief Abstract base class for a PosePropertyReporter factory; the Creator class is responsible for
/// creating a particular PosePropertyReporter class.
class PosePropertyReporterCreator : public utility::pointer::ReferenceCount
{
public:
	PosePropertyReporterCreator();
	~PosePropertyReporterCreator() override;

	virtual PosePropertyReporterOP create_reporter() const = 0;
	virtual std::string keyname() const = 0;
};

typedef utility::pointer::shared_ptr< PosePropertyReporterCreator > PosePropertyReporterCreatorOP;
typedef utility::pointer::shared_ptr< PosePropertyReporterCreator const > PosePropertyReporterCreatorCOP;

} //namespace rosetta_scripts
} //namespace protocols

#endif
