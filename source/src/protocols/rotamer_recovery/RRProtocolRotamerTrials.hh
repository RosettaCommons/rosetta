// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRProtocolRotamerTrials.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_rotamer_recovery_RRProtocolRotamerTrials_hh
#define INCLUDED_protocols_rotamer_recovery_RRProtocolRotamerTrials_hh

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocolRotamerTrials.fwd.hh>
#include <protocols/rotamer_recovery/RRProtocol.hh>

// Project Headers
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>
#include <protocols/rotamer_recovery/RRReporter.fwd.hh>

// Platform Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

//Auto Headers
#include <utility/vector1.hh>

namespace protocols {
namespace rotamer_recovery {

class RRProtocolRotamerTrials : public RRProtocol {

public: // constructors destructors

	RRProtocolRotamerTrials();

	~RRProtocolRotamerTrials() override;

	RRProtocolRotamerTrials( RRProtocolRotamerTrials const & );

public: // public interface

	
	std::string
	get_name() const override;

	
	std::string
	get_parameters() const override;

	void
	run(
		RRComparerOP comparer,
		RRReporterOP reporter,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & score_function,
		core::pack::task::PackerTask const & packer_task) override;

};

} // namespace
} // namespace

#endif // include guard
