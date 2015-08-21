// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/DeleteRegionMover.hh
/// @brief   Base class for graftmovers
/// @author  Jared Adolf-Bryfogle

#ifndef INCLUDED_protocols_grafting_DeleteRegionMover_HH
#define INCLUDED_protocols_grafting_DeleteRegionMover_HH

#include <protocols/grafting/simple_movers/DeleteRegionMover.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>


namespace protocols {
namespace grafting {
namespace simple_movers {

/// @brief Delete a region of a pose. Mover Wrapper to grafting utility function.
class DeleteRegionMover : public  protocols::moves::Mover {

public:

	DeleteRegionMover();
	DeleteRegionMover(core::Size res_start, core::Size res_end);

	DeleteRegionMover(DeleteRegionMover const & src);

	virtual ~DeleteRegionMover();

	virtual void
	apply(core::pose::Pose & pose);


public:

	/// @brief Set the region of the pose where deletion will occur
	void
	region(core::Size res_start, core::Size res_end);

	std::pair<core::Size, core::Size>
	region() const;

	/// @brief Set the first residue that will be deleted.
	void
	start(core::Size res_start);

	core::Size
	start() const;

	/// @brief Set the last residue that will be deleted.
	void
	end(core::Size res_end);

	core::Size
	end() const;

public:

	virtual std::string
	get_name() const;

	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;

	virtual void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		moves::Movers_map const & movers,
		Pose const & pose);

private:
	core::Size start_;
	core::Size end_;
	core::Size nter_overhang_;
	core::Size cter_overhang_;
	TagCOP tag_; //This is so pdb_num can be parsed at apply time instead of construction time.
};


}
}
}

#endif  // INCLUDED_protocols_grafting_DeleteRegionMover_HH
