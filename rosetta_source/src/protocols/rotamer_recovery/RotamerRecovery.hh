// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RotamerRecovery.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_rotamer_recovery_RotamerRecovery_hh
#define INCLUDED_protocols_rotamer_recovery_RotamerRecovery_hh

// Unit Headers
#include <protocols/rotamer_recovery/RotamerRecovery.fwd.hh>
#include <protocols/rotamer_recovery/RRReporter.fwd.hh>
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/Mover.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
// AUTO-REMOVED #include <ostream>

#include <protocols/moves/Mover.fwd.hh>


namespace protocols {
namespace rotamer_recovery {

class RotamerRecovery : public utility::pointer::ReferenceCount {

public: // constructors destructors

	///@brief default constructor
	RotamerRecovery();

	///@brief specify comparer and reporter
	RotamerRecovery(
		RRReporterOP reporter,
		RRComparerOP comparer);

	//@brief convience constructor that selects reporter and comparer based on the names given.
	RotamerRecovery(
		std::string const & reporter,
		std::string const & output_fname,
		std::string const & comparer);


	///@brief destructor
	virtual ~RotamerRecovery();

	///@brief copy constructor
	RotamerRecovery( RotamerRecovery const & src);

public: // public interface

	virtual
	void
	reset_recovery();

	virtual
	void
	register_options() const;

	virtual
	bool
	measure_rotamer_recovery(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2
	);

	core::Real
	rtmin_rotamer_recovery(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & score_function,
		core::pack::task::PackerTask const & packer_task
	);


	virtual
	void
	compare_rotamers(
		core::pose::Pose const & pose,
		moves::Mover & mover,
		utility::vector1< core::Size > const & res_ids
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

	RRReporterOP reporter_;
	RRComparerOP comparer_;

	bool ignore_unrecognized_res_;
};

} // namespace rotamer_recovery
} // namespace protocols

#endif
