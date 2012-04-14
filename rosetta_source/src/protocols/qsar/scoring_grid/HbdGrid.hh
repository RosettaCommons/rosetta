// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/HbaGrid.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_HbdGrid_hh
#define INCLUDED_protocols_qsar_scoring_grid_HbdGrid_hh

#include <protocols/qsar/scoring_grid/HbdGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>

#include <utility/vector1.hh>



namespace protocols {
namespace qsar {
namespace scoring_grid {

class HbdGrid : public SingleGrid
{
public:
	HbdGrid();
	HbdGrid(core::Real weight);
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude);
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center);
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> ligand_chain_ids_to_exclude);

	void parse_my_tag(utility::tag::TagPtr const tag);

private:

	core::Real radius_;
	core::Real width_;
	core::Real magnitude_;

};

}
}
}

#endif /* HBAGRID_FWD_HH_ */
