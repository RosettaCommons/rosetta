// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/LigandDesign.hh
///
/// @brief
/// @author Gordon Lemmon


#ifndef INCLUDED_protocols_ligand_docking_LigandDesign_hh
#define INCLUDED_protocols_ligand_docking_LigandDesign_hh

#include <protocols/ligand_docking/LigandDesign.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

/// @brief
///
/// @details
///
class LigandDesign : public protocols::moves::Mover{

public:
	LigandDesign();
	virtual ~LigandDesign();
	LigandDesign(LigandDesign const & that);
	virtual void apply( core::pose::Pose & pose );

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual std::string get_name() const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);


private:
	std::string option_file_;

	void
	add_scores_to_job(
		core::pose::Pose & pose
	);

	void set_fragments();
	void fragments_to_string() const;
	core::conformation::ResidueCOPs fragments_;

}; // class LigandDesign

bool grow(core::pose::Pose, core::Size start, core::Size end);
bool has_incomplete_connections(core::pose::Pose pose, core::Size start, core::Size const end);
bool passes_filters(core::pose::Pose const & pose, core::Size start, core::Size const end);
core::Size random_connection(core::conformation::ResidueCOP residue);
utility::vector1<core::Size> get_incomplete_connections(core::conformation::ResidueCOP residue);
utility::vector1<core::Size> find_unconnected_residues( core::pose::Pose const & pose, core::Size start, core::Size const end );

} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_LigandDesign_HH
