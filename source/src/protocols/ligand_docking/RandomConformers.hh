// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  this class chooses a random conformer for every residue in a specified chain
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_RandomConformers_hh
#define INCLUDED_protocols_ligand_docking_RandomConformers_hh

#include <protocols/moves/Mover.hh>
#include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

class RandomConformers : public protocols::moves::Mover
{
public:
	RandomConformers();
	~RandomConformers() override;
	RandomConformers(RandomConformers const & that);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	std::string get_name() const override;

	//void set_chain(std::string chain);
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	void apply(core::pose::Pose & pose) override;

	// Undefined, commenting out to make PyRosetta compile
	//RandomConformers(
	// core::pose::Pose & pose,
	// std::set< protocols::ligand_docking::ResidueTorsionRestraintsOP > const & restraints
	//);

private:
	std::string chain_;
	void apply_residue( core::Size const residue_id, core::pose::Pose & pose );
};

} //namespace ligand_docking
} //namespace protocols

#endif
