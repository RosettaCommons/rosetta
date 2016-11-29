// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRComparer.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_rotamer_recovery_RRComparer_hh
#define INCLUDED_protocols_rotamer_recovery_RRComparer_hh

// Unit Headers
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers

#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace rotamer_recovery {

/// @brief The comparison functionality for the rotamer recovery test
///
/// Besides implementing the interface given in the base class
/// RRComparer each RRComparer should have an entry in the conevience
/// RotamerRecovery constructor so its use can be indicated by name.
class RRComparer : public utility::pointer::ReferenceCount {

public: // constructors destructors
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~RRComparer() override;

	// RRComparer();
	//
	// ~RRComparer();
	//
	// RRComparer( RRComparer const & src );

public: // public interface

	virtual
	std::string
	get_name() const = 0;

	virtual
	std::string
	get_parameters() const = 0;

	virtual
	void
	set_recovery_threshold(
		core::Real const recovery_threshold) = 0;

	virtual
	void
	set_absolute_threshold(
		core::Real const absolute_threshold) = 0;

	virtual
	bool
	measure_rotamer_recovery(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & score,
		bool & recovered ) = 0;

};

class RRComparerRotBins : public RRComparer {

public: // constructors destructors

	RRComparerRotBins();

	~RRComparerRotBins() override;

	RRComparerRotBins( RRComparerRotBins const & );

public: // public interface


	std::string
	get_name() const override;


	std::string
	get_parameters() const override;


	virtual
	void
	set_recovery_threshold(
		core::Real const recovery_threshold) override;
	
	virtual
	void
	set_absolute_threshold(
		core::Real const absolute_threshold) override;

	bool
	measure_rotamer_recovery(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & score,
		bool & recovered) override;

private: // data members

	core::Real recovery_threshold_;
	core::Real absolute_threshold_;

};

class RRComparerChiDiff : public RRComparer {

public: // constructors destructors

	RRComparerChiDiff();

	~RRComparerChiDiff() override;

	RRComparerChiDiff( RRComparerChiDiff const & );

public: // public interface

	void set_recovery_threshold( core::Real const setting ) override;
	void set_absolute_threshold( core::Real const setting ) override;

	virtual
	void
	set_max_chi_considered( core::Size const max_chi );


	std::string
	get_name() const override;



	std::string
	get_parameters() const override;


	bool
	measure_rotamer_recovery(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & score,
		bool & recovered) override;

private: // data members
	core::Real tolerance_;
	core::Real absolute_threshold_;
	bool limit_chi_angles_;
	core::Size max_chi_considered_; // only relevant if limit_chi_angles_ is set

};

} // rotamer_recovery
} // protocols

#endif // include guard
