// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief PeptideCyclizeMover a pose, by default with a peptide bond to N and C terminal and constraints.
/// @author Parisa Hosseinzadeh (parisah@uw.edu) & Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_PeptideCyclizeMover_hh
#define INCLUDED_protocols_cyclic_peptide_PeptideCyclizeMover_hh
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <core/pose/Pose.hh>

//Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

#include <protocols/cyclic_peptide/PeptideCyclizeMover.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace cyclic_peptide {

class PeptideCyclizeMover : public moves::Mover {
public:
	PeptideCyclizeMover();
	~PeptideCyclizeMover() override;

	void apply( Pose & ) override;
	std::string get_name() const override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & ) override;

	virtual void set_bond(core::Size, std::string, core::Size, std::string, bool, bool);

	virtual void set_distance(core::Size, std::string, core::Size, std::string, std::string);

	virtual void set_angle(core::Size, std::string, core::Size, std::string, core::Size, std::string, std::string);

	virtual void set_torsion(core::Size, std::string, core::Size, std::string, core::Size, std::string, core::Size, std::string,std::string);

	virtual void set_selector(core::select::residue_selector::ResidueSelectorCOP);

	virtual void set_default
	();

private:

	core::Size counter_;

	utility::vector1<core::Size> res1_;
	utility::vector1<std::string> atom1_;
	utility::vector1<core::Size> res2_;
	utility::vector1<std::string> atom2_;
	utility::vector1<bool> add_termini_;
	utility::vector1<bool> rebuild_fold_tree_;
	bool bond_assigned_;

	utility::vector1<core::Size> res1_dist_;
	utility::vector1<std::string> atom1_dist_;
	utility::vector1<core::Size> res2_dist_;
	utility::vector1<std::string> atom2_dist_;
	utility::vector1<std::string> cst_func_dist_;
	bool distance_assigned_;

	utility::vector1<core::Size> res1_angle_;
	utility::vector1<std::string> atom1_angle_;
	utility::vector1<std::string> atom1_angle2_;
	utility::vector1<core::Size> res2_angle_;
	utility::vector1<std::string> atom2_angle_;
	utility::vector1<std::string> atom2_angle2_;
	utility::vector1<Size> res_center_;
	utility::vector1<Size> res_center2_;
	utility::vector1<std::string> atom_center_;
	utility::vector1<std::string> atom_center2_;
	utility::vector1<std::string> cst_func_angle_;
	utility::vector1<std::string> cst_func_angle2_;
	bool angle_assigned_;

	utility::vector1<core::Size> res1_torsion_;
	utility::vector1<std::string> atom1_torsion_;
	utility::vector1<core::Size> res2_torsion_;
	utility::vector1<std::string> atom2_torsion_;
	utility::vector1<core::Size> res3_torsion_;
	utility::vector1<std::string> atom3_torsion_;
	utility::vector1<core::Size> res4_torsion_;
	utility::vector1<std::string> atom4_torsion_;
	utility::vector1<std::string> cst_func_torsion_;
	bool torsion_assigned_;

	core::Real angle1_;
	core::Real angle2_;
	core::Real distance_;

	std::string distance_func_;
	std::string angle1_func_;
	std::string angle2_func_;

	core::select::residue_selector::ResidueSelectorCOP selector_;

	void get_all(core::select::residue_selector::ResidueSubset, core::pose::Pose const & );
	void get_values();

};

} // moves
} // protocols

#endif
