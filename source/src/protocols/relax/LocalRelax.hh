// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/relax/LocalRelax.hh
/// @brief
/// @details
/// @author


#ifndef INCLUDED_protocols_relax_LocalRelax_hh
#define INCLUDED_protocols_relax_LocalRelax_hh

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>

#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

// option includes
#include <basic/options/keys/relax.OptionKeys.gen.hh>


namespace protocols {
namespace relax {


class LocalRelax : public protocols::moves::Mover {
private:
	core::Size NCYC_, NEXP_;
	core::Real K_, max_iter_;
	bool ramp_cart_;
	bool verbose_;

	core::scoring::ScoreFunctionOP pack_sfxn_, min_sfxn_;

public:
	LocalRelax();

	virtual std::string get_name() const { return "LocalRelax"; }

	/// @brief one cycle of local optimization
	void
	optimization_loop(
			Pose & pose,
			core::pack::task::PackerTaskOP ptask,
			core::kinematics::MoveMapOP mm,
			core::Real, core::Real);

	/// @brief get matrix of interacting residues
	void
	get_neighbor_graph(
		Pose & pose,
		utility::vector1< utility::vector1<bool> > &neighbor);

	/// @brief RS integration
	void parse_my_tag(
	  utility::tag::TagCOP tag,
	  basic::datacache::DataMap & data,
	  protocols::filters::Filters_map const &,
	  protocols::moves::Movers_map const &,
	  core::pose::Pose const &
	);

	/// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::moves::MoverOP clone() const {
		return protocols::moves::MoverOP( new LocalRelax(*this) );
	}

	void apply( core::pose::Pose & pose);
};


}
} // protocols

#endif
