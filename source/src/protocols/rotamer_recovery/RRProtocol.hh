// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RRProtocol.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_rotamer_recovery_RRProtocol_hh
#define INCLUDED_protocols_rotamer_recovery_RRProtocol_hh

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocol.fwd.hh>

// Project Headers
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>
#include <protocols/rotamer_recovery/RRReporter.fwd.hh>

// Platform Headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
#include <utility/vector1.hh>
namespace protocols {
namespace rotamer_recovery {

/// @brief The protocol to run to compute the rotamer recovery the rotamer recovery test
///
/// Besides implementing the interface given in the base class
/// RRProtocol each RRProtocol should have an entry in the convenience
/// RotamerRecovery constructor so its use can be indicated by name.
class RRProtocol : public utility::pointer::ReferenceCount {

public: // constructors destructors
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~RRProtocol();

//	RRProtocol();
//
//	~RRProtocol();
//
//	RRProtocol( RRProtocol const & src );

protected:
	bool
	measure_rotamer_recovery(
		RRComparerOP comparer,
		RRReporterOP reporter,
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2);

public: // public interface

	virtual
	std::string
	get_name() const = 0;

	virtual
	std::string
	get_parameters() const = 0;

	virtual
	void
	run(
		RRComparerOP comparer,
		RRReporterOP reporter,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & score_function,
		core::pack::task::PackerTask const & packer_task) = 0;

};

} // namespace
} // namespace

#endif // include guard
