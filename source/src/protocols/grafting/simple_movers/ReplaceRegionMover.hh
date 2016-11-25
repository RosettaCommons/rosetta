// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/grafting/ReplaceRegionMover.hh
/// @brief   Base class for graftmovers
/// @author  Jared Adolf-Bryfogle

#ifndef INCLUDED_protocols_grafting_simple_movers_ReplaceRegionMover_HH
#define INCLUDED_protocols_grafting_simple_movers_ReplaceRegionMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/grafting/simple_movers/ReplaceRegionMover.fwd.hh>

#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace grafting {
namespace simple_movers {

/// @brief Replace a region of residues with another from a pose piece.  Wrapper to grafting utility function
/// @details Specify the start residue of the source pose, the residue where replacement will begin,
///  In addition to the number of residues (span) .
///
class ReplaceRegionMover : public  protocols::moves::Mover {

public:
	ReplaceRegionMover(bool copy_pdbinfo = false);
	ReplaceRegionMover(core::pose::Pose const & src_pose,
		core::Size const src_pose_start,
		core::Size const target_pose_start,
		core::Size const span,
		bool copy_pdbinfo = false);

	ReplaceRegionMover(ReplaceRegionMover const & src);

	virtual ~ReplaceRegionMover();

	void
	apply(core::pose::Pose & pose) override;

public:

	core::Size src_pose_start() const;
	void           src_pose_start(core::Size start);

	core::Size target_pose_start() const;
	void           target_pose_start(core::Size start);

	// Undefined, commenting out to fix PyRosetta build  core::Size span() const;
	// Undefined, commenting out to fix PyRosetta build  void       span(core::Size span);

	void src_pose(core::pose::Pose const & from_pose);

public:

	// XRW TEMP  virtual std::string
	// XRW TEMP  get_name() const;

	protocols::moves::MoverOP
	clone() const override;

	protocols::moves::MoverOP
	fresh_instance() const override;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		moves::Movers_map const & movers,
		Pose const & pose
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Size src_pose_start_;
	core::Size target_pose_start_;
	core::Size span_;
	bool copy_pdbinfo_;
	core::pose::PoseOP src_pose_;
	TagCOP tag_; //This is so pdb_num can be parsed at apply time instead of construction time.
};


}
}
}


#endif  // INCLUDED_protocols_grafting_ReplaceRegionMover_HH
