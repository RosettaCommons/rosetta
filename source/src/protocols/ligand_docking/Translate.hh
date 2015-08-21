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
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

//// Scoring grid headers

//// Project Headers
#include <utility/vector1.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

struct Translate_info{ // including default values

public:
	core::Size chain_id;
	core::Size jump_id;
	Distribution distribution;
	core::Real angstroms;
	core::Size cycles;
	bool force;
	Translate_info(): chain_id(0), jump_id(0), distribution(Uniform), angstroms(0), cycles(0), force(false){};
};

class Translate : public protocols::moves::Mover
{
public:
	Translate();
	Translate(Translate_info translate_info);
	virtual ~Translate();
	Translate(Translate const & that);

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual std::string get_name() const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map,
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
	utility::pointer::shared_ptr<core::grid::CartGrid<int> > grid_;
	utility::vector1<core::Size> chain_ids_to_exclude_; // these are invisible the translation grid, so ligand can land on top.
	std::vector<core::Size> tag_along_jumps_; // these guys tag along, such as waters and metals

	void translate_ligand(
		utility::pointer::shared_ptr<core::grid::CartGrid<int> >  const & grid,
		core::Size const jump_id,
		core::pose::Pose & pose
	);

	void translate_ligand(core::Size const jump_id,core::pose::Pose & pose, core::Size const & residue_id);

	void uniform_translate_ligand(
		const utility::pointer::shared_ptr<core::grid::CartGrid<int> >  & grid,
		const core::Size jump_id,
		core::pose::Pose & pose
	);

	void gaussian_translate_ligand(
		const utility::pointer::shared_ptr<core::grid::CartGrid<int> >  & grid,
		const core::Size jump_id,
		core::pose::Pose & pose
	);

};

} //namespace ligand_docking
} //namespace protocols

#endif
