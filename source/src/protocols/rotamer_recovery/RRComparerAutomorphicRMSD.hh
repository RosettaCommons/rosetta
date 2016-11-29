// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRComparerAutomorphicRMSD.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_rotamer_recovery_RRComparerAutomorphicRMSD_hh
#define INCLUDED_protocols_rotamer_recovery_RRComparerAutomorphicRMSD_hh

// Unit Headers
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRComparerAutomorphicRMSD.fwd.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers

#include <utility/vector1.hh>


namespace protocols {
namespace rotamer_recovery {


class RRComparerAutomorphicRMSD : public RRComparer {

public: // constructors destructors

	RRComparerAutomorphicRMSD();

	~RRComparerAutomorphicRMSD() override;

	RRComparerAutomorphicRMSD( RRComparerAutomorphicRMSD const & );

public: // public interface


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

	virtual
	void
	set_include_backbone_atoms(
		bool const include_backbone_atoms
	);

	virtual
	bool
	get_include_backbone_atoms() const;


	void
	set_recovery_threshold(
		core::Real const recovery_threshold
	) override;

	virtual
	core::Real
	get_recovery_threshold() const;

	virtual
	void
	set_absolute_threshold(
		core::Real const absolute_threshold
	) override;

	virtual
	core::Real
	get_absolute_threshold() const;


private: // data members

	bool include_backbone_atoms_;
	core::Real recovery_threshold_;
	core::Real absolute_threshold_;

};

} // rotamer_recovery
} // protocols

#endif // include guard
