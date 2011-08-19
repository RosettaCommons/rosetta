// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/qsar/qsarMover.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_qsarMover_hh
#define INCLUDED_protocols_qsar_qsarMover_hh

#include <protocols/qsar/qsarMover.fwd.hh>
#include <protocols/qsar/qsarMap.fwd.hh>
//#include <protocols/qsar/qsarTypeManager.fwd.hh>
#include <protocols/qsar/scoring_grid/GridManager.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace qsar {

class qsarMover : public protocols::moves::Mover
{
public:
	qsarMover(core::Real width, core::Real resolution);
	void set_chain(std::string chain_id);
	//qsarMover(core::Real width, core::Real resolution);
	void add_grid(std::string grid_name);
	void write_all_grids(std::string prefix);
	std::map<std::string,core::Real> get_cached_scores();
	virtual void apply(core::pose::Pose & pose);
//	virtual void parse_my_tag(utility::tag::TagPtr const tag,
//			moves::DataMap & data,
//			filters::Filters_map const& filters,
//			Movers_map const & movers,
//			core::pose::Pose const & pose)
	virtual std::string get_name() const;

private:
	scoring_grid::GridManagerOP grid_manager_;
	qsarMapOP qsar_map_;
	//core::conformation::ResidueOP current_ligand_;
	std::string chain_id_;
	utility::vector1<std::string> grids_to_make_;
	bool initialize_;
};

}
}

#endif /* QSARMOVER_HH_ */
