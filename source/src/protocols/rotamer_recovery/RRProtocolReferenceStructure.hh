// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRProtocolReferenceStructure.hh
/// @brief  Preform the rotamer recovery against a reference structure
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_rotamer_recovery_RRProtocolReferenceStructure_hh
#define INCLUDED_protocols_rotamer_recovery_RRProtocolReferenceStructure_hh

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocolReferenceStructure.fwd.hh>
#include <protocols/rotamer_recovery/RRProtocol.hh>

// Project Headers
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>
#include <protocols/rotamer_recovery/RRReporter.fwd.hh>

// Platform Headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

//Auto Headers
#include <core/types.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace rotamer_recovery {

class RRProtocolReferenceStructure : public RRProtocol {

public: // constructors destructors

	RRProtocolReferenceStructure();

	RRProtocolReferenceStructure(
		core::pose::PoseCOP reference_pose);

	~RRProtocolReferenceStructure() override;

	RRProtocolReferenceStructure(
		RRProtocolReferenceStructure const & src);

public: // public interface

	std::string
	get_name() const override;

	std::string
	get_parameters() const override;

	void
	reference_structure(
		core::pose::PoseCOP reference_pose);

	void
	run(
		RRComparerOP comparer,
		RRReporterOP reporter,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::pack::task::PackerTask const & packer_task) override;

private: // member data
	core::pose::PoseCOP reference_pose_;
};

} // namespace
} // namespace

#endif // include guard
