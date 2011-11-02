// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RRReporter.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_rotamer_recovery_RRReporter_hh
#define INCLUDED_protocols_rotamer_recovery_RRReporter_hh

// Unit Headers
#include <protocols/rotamer_recovery/RRReporter.fwd.hh>

// Project Headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
// AUTO-REMOVED #include <ostream>

#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace rotamer_recovery {

///@brief The reporting functionality for the rotamer recovery test
///
/// Becides implementing the interface given in the base class
/// RRReporter each RRReporter should have an entry in the conevience
/// RotamerRecovery constructor so its use can be indicated by name.
class RRReporter : public utility::pointer::ReferenceCount {

public: // constructors destructors


public: // public interface

	virtual
	void
	set_comparer_info(
		std::string const & comparer_name,
		std::string const & comparer_params) = 0;

	virtual
	void
	reset_recovery() = 0;

	virtual
	void
	report_rotamer_recovery(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real score,
		bool recovered
	) = 0;

	virtual
	core::Real
	recovery_rate() const = 0;

	virtual
	void
	show(std::ostream & out ) const = 0;

	virtual
	void
	show() const = 0;

};

class RRReporterSimple : public RRReporter {

public: // constructors destructors

	RRReporterSimple();

	RRReporterSimple( RRReporterSimple const & );

	~RRReporterSimple();


public: // public interface

	void
	set_comparer_info(
		std::string const & /*comparer_name*/,
		std::string const & /*comparer_params*/){}

	virtual
	void
	reset_recovery();

	virtual
	void
	report_rotamer_recovery(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real score,
		bool recovered
	);

	virtual
	core::Real
	recovery_rate() const;

	virtual
	void
	show(std::ostream & out ) const;

	virtual
	void
	show() const;

private: // data members

	core::Size residues_considered_;
	core::Size rotamers_recovered_;

};

} // namespace rotamer_recovery
} // namespace protocols


#endif // include guard
