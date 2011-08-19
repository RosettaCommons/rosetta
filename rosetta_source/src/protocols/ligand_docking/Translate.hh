// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_Translate_hh
#define INCLUDED_protocols_ligand_docking_Translate_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/ligand_docking/Translate.fwd.hh>
#include <protocols/ligand_docking/DistributionMap.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

//// Scoring grid headers
#include <protocols/qsar/scoring_grid/GridManager.fwd.hh>

//// Project Headers
#include <core/pose/Pose.hh>
#include <core/grid/CartGrid.hh>
#include <utility/vector1.hh>

#include <set>

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

struct Translate_info{ // including default values

public:
	std::string chain;
	Distribution distribution;
	core::Real angstroms;
	core::Size cycles;
	Translate_info(): chain(""), distribution(Uniform), angstroms(0), cycles(0){};
};

class Translate : public protocols::moves::Mover
{
public:
	Translate();
	virtual ~Translate();
	Translate(Translate const & that);

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

	core::Size get_chain_id(core::pose::Pose const & pose);
	void add_excluded_chains(
			std::set<core::Size>::const_iterator begin,
			std::set<core::Size>::const_iterator end
	);

private:
	Translate_info translate_info_;
	utility::pointer::owning_ptr<core::grid::CartGrid<int> > grid_;
	utility::vector1<core::Size> chain_ids_to_exclude_; // these are invisible the translation grid, so ligand can land on top.
	qsar::scoring_grid::GridManagerOP grid_manager_;

	void translate_ligand(
			utility::pointer::owning_ptr<core::grid::CartGrid<int> >  const & grid,
			core::Size const jump_id,
			core::pose::Pose & pose
	);

	void translate_ligand(core::Size const jump_id,core::pose::Pose & pose, core::Size const & residue_id);


};

} //namespace ligand_docking
} //namespace protocols

#endif
