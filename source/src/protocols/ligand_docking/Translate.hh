// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_Translate_hh
#define INCLUDED_protocols_ligand_docking_Translate_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/ligand_docking/Translate.fwd.hh>
#include <protocols/ligand_docking/DistributionMap.hh>
#include <protocols/qsar/scoring_grid/GridSet.fwd.hh>

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

	Distribution distribution = Uniform;
	core::Real angstroms = 0;
	core::Size cycles = 0;
	bool force = false;

	Translate_info() = default;

	core::Size chain_id( core::pose::Pose const & pose ) const;
	char chain_letter( core::pose::Pose const & pose ) const;
	core::Size jump_id( core::pose::Pose const & pose ) const;

	void set_chain_id( core::Size id );
	void set_chain_letter( std::string const & str);

private:
	bool by_string_ = true; // Is the chain represented by a chain letter or a chain ID?
	std::string chain_string_ = "X";
	core::Size chain_number_ = 0;
};

class Translate : public protocols::moves::Mover
{
public:
	Translate();
	Translate(Translate_info translate_info);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map
	) override;

	void apply(core::pose::Pose & pose) override;

	core::Size get_chain_id(core::pose::Pose const & pose);
	void add_excluded_chains(
		std::set<core::Size>::const_iterator begin,
		std::set<core::Size>::const_iterator end
	);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	void translate_ligand(
		utility::pointer::shared_ptr<core::grid::CartGrid<int> >  const & grid,
		core::Size const jump_id,
		core::pose::Pose & pose
	);

	void translate_ligand(
		qsar::scoring_grid::GridSetCOP grid_set,
		core::Size const jump_id,
		core::pose::Pose & pose,
		core::Size const & residue_id);

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

private:
	protocols::qsar::scoring_grid::GridSetCOP grid_set_prototype_;
	Translate_info translate_info_;
	utility::pointer::shared_ptr<core::grid::CartGrid<int> > grid_;
	utility::vector1<core::Size> chain_ids_to_exclude_; // these are invisible the translation grid, so ligand can land on top.
	utility::vector1<std::string> tag_along_chains_; // these chains tag along, such as waters and metals (will also be excluded)

};

} //namespace ligand_docking
} //namespace protocols

#endif
