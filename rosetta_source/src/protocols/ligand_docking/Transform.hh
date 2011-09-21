// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /grid/src/protocols/ligand_docking/Transform.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_ligand_docking_Transform_hh
#define INCLUDED_protocols_ligand_docking_Transform_hh

#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/ligand_docking/Transform.fwd.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>

namespace protocols {
namespace ligand_docking {

struct Transform_info{ // including default values

public:
	std::string chain;
	core::Size chain_id;
	core::Size jump_id;
	//Distribution distribution;
	core::Real move_distance;
	core::Real box_size;
	core::Real angle;
	core::Size cycles;
	core::Real temperature;
	Transform_info(): chain(""), move_distance(0),box_size(0), angle(0), cycles(0){};
};

class Transform: public protocols::moves::Mover
{
public:
	Transform();
	virtual ~Transform();
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual std::string get_name() const;

	void parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	void apply(core::pose::Pose & pose);

	//core::Size get_chain_id(core::pose::Pose const & pose);

	/*
	void add_excluded_chains(
		std::set<core::Size>::const_iterator begin,
		std::set<core::Size>::const_iterator end
	);
	*/

private:
	void transform_ligand(core::pose::Pose & pose);


	void change_conformer(core::pose::Pose & pose, core::Size const & seqpos);

private:
	//qsar::scoring_grid::GridManagerOP grid_manager_;
	Transform_info transform_info_;
	utility::vector1< core::conformation::ResidueOP >  ligand_conformers_;


};

}
}

#endif /* TRANSFORM_HH_ */
