// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/movers/RotateSegmentMover.hh
/// @brief Rotates a segment in the pose
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_movers_RotateSegmentMover_hh
#define INCLUDED_protocols_denovo_design_movers_RotateSegmentMover_hh

// Unit headers
#include <protocols/denovo_design/movers/RotateSegmentMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <numeric/xyzVector.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace movers {

class RotationRange {
public:
	RotationRange( core::Real start_val, core::Real const stop_val );
	core::Real start;
	core::Real stop;

private:
	RotationRange();
};

///@brief Rotates a segment in the pose
class RotateSegmentMover : public protocols::moves::Mover {
public:
	typedef utility::vector1< std::string > Residues;
	typedef utility::vector1< std::string > Atoms;
	typedef utility::vector1< core::id::AtomID > AtomIDs;

public:
	RotateSegmentMover();

	virtual
	~RotateSegmentMover();

	virtual void
	apply( core::pose::Pose & pose );

	virtual std::string
	get_name() const;

	static std::string
	class_name();

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	virtual protocols::moves::MoverOP
	fresh_instance() const;

	virtual protocols::moves::MoverOP
	clone() const;

public:
	void
	apply( core::pose::Pose & pose, core::Real const random ) const;

	void
	set_rotation_range( core::Real const start, core::Real const stop );

private:
	core::Real
	select_rotation( core::Real const random ) const;

	AtomIDs
	target_atoms(
		core::pose::Pose const & pose,
		components::StructureData const & sd ) const;

private:
	RotationRange range_;
	Residues residues_;
	Atoms atoms_;
	SegmentNames frozen_segments_;
};

// calculates a dihedral from four coordinates in degrees
core::Real
calc_dihedral(
	core::Vector const & p1,
	core::Vector const & p2,
	core::Vector const & p3,
	core::Vector const & p4 );

} //protocols
} //denovo_design
} //movers


#endif //protocols/denovo_design/movers_RotateSegmentMover_hh

