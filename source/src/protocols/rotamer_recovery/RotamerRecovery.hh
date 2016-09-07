// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RotamerRecovery.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_rotamer_recovery_RotamerRecovery_hh
#define INCLUDED_protocols_rotamer_recovery_RotamerRecovery_hh

// Unit Headers
#include <protocols/rotamer_recovery/RotamerRecovery.fwd.hh>
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>
#include <protocols/rotamer_recovery/RRProtocol.fwd.hh>
#include <protocols/rotamer_recovery/RRReporter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#ifdef PYROSETTA
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRProtocol.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>
#endif


namespace protocols {
namespace rotamer_recovery {

class RotamerRecovery : public utility::pointer::ReferenceCount {

public: // constructors destructors

	/// @brief default constructor
	RotamerRecovery();

	/// @brief specify comparer and reporter
	RotamerRecovery(
		RRProtocolOP protocol,
		RRComparerOP comparer,
		RRReporterOP reporter);

	/// @brief destructor
	~RotamerRecovery() override;

	/// @brief copy constructor
	RotamerRecovery( RotamerRecovery const & src);

public: // public interface

	virtual
	void
	reset_recovery();

	virtual
	void
	register_options() const;

	core::Real
	run(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & score_function,
		core::pack::task::PackerTask const & packer_task
	);

	void
	init_rotamer_recovery_with_options(
		RotamerRecovery & rotamer_recovery
	);

	void
	init_with_options();

	void
	set_ignore_unrecognized_res(
		bool const ignore_unrecognized_res
	);

	bool
	get_ignore_unrecognized_res();

	virtual
	void
	show(std::ostream & out ) const;

	virtual
	void
	show();

	virtual
	core::Real
	recovery_rate() const;

private: // data

	RRProtocolOP protocol_;
	RRComparerOP comparer_;
	RRReporterOP reporter_;

	bool ignore_unrecognized_res_;
};

} // namespace rotamer_recovery
} // namespace protocols

#endif
