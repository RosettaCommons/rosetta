// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#ifndef INCLUDED_protocols_cyclic_peptide_DeclareBond_hh
#define INCLUDED_protocols_cyclic_peptide_DeclareBond_hh

#include <protocols/cyclic_peptide/DeclareBond.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>


namespace protocols {
namespace cyclic_peptide {

class DeclareBond : public moves::Mover {
public:
	DeclareBond();
	~DeclareBond() override;

	void apply( Pose & ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap & ) override;

	void
	set( core::Size const res1,
		std::string const & atom1,
		core::Size const res2,
		std::string const & atom2,
		bool const add_termini,
		bool const run_kic = false,
		core::Size const kic_res1 = 0,
		core::Size const kic_res2 = 0,
		bool const rebuild_fold_tree = false
	);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	set_selector_for_res1( core::select::residue_selector::ResidueSelectorCOP selector1 ){
		selector1_ = selector1;
	}

	void
	set_selector_for_res2( core::select::residue_selector::ResidueSelectorCOP selector2 ){
		selector2_ = selector2;
	}

private:
	core::select::residue_selector::ResidueSelectorCOP selector1_ = nullptr;
	core::Size res1_;
	std::string atom1_;

	core::select::residue_selector::ResidueSelectorCOP selector2_ = nullptr;
	core::Size res2_;
	std::string atom2_;

	bool add_termini_;
	bool run_kic_;
	core::Size kic_res1_;
	core::Size kic_res2_;


	/// @brief Should the foldtree be rebuilt?
	bool rebuild_fold_tree_;
};

} // moves
} // protocols

#endif
