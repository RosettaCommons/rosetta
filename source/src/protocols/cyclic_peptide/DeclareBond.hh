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

namespace protocols {
namespace cyclic_peptide {

class DeclareBond : public moves::Mover {
public:
	DeclareBond();
	virtual ~DeclareBond();

	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & );

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

private:
	Size res1_;
	std::string atom1_;
	Size res2_;
	std::string atom2_;

	bool add_termini_;
	bool run_kic_;
	Size kic_res1_;
	Size kic_res2_;


	/// @brief Should the foldtree be rebuilt?
	bool rebuild_fold_tree_;
};

} // moves
} // protocols

#endif
