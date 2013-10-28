// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/AtomicContactCountFilter.hh
/// @brief
/// @author Alex Ford (fordas@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_AtomicContactCountFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_AtomicContactCountFilter_hh

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <string>

#include <protocols/protein_interface_design/filters/AtomicContactCountFilter.fwd.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

class AtomicContactCountFilter : public protocols::filters::Filter
{
	public:
		AtomicContactCountFilter();
		AtomicContactCountFilter(core::Real distance_cutoff);
    AtomicContactCountFilter( AtomicContactCountFilter const & copy );

    virtual protocols::filters::FilterOP clone() const;
    virtual protocols::filters::FilterOP fresh_instance() const;
    virtual ~AtomicContactCountFilter();

    // @brief Filter name
    virtual std::string name() const { return "AtomicContactCountFilter"; }

		void initialize_all_atoms(core::pack::task::TaskFactoryOP task_factory = NULL);

		void initialize_cross_jump(core::Size jump, std::string sym_dof_name = "", core::pack::task::TaskFactoryOP task_factory = NULL, bool normalize_by_sasa = false);

		void initialize_cross_chain(core::pack::task::TaskFactoryOP task_factory = NULL, bool normalize_by_sasa = false, bool detect_chains_for_interface_by_task = false);

    void parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

		/// @brief Returns true if the given pose passes the filter, false otherwise.
		virtual bool apply( core::pose::Pose const & /*pose*/ ) const { return true; }

		/// @brief used to report filter internals through a score or silent file
		// to determine that derived class has not overridden }
		virtual core::Real report_sm( core::pose::Pose const & pose ) const { return compute(pose); }

		core::Real compute( core::pose::Pose const &) const;

	private:
		//TODO Alex Ford Add support for specified chain mode, need to add chain-id resolution to tag parsing.
		//TODO Alex Ford Likely even better to abstract chain detection & jump code into an interface sasa calculator
		// and compute the sasa value there.
		enum Mode { ALL, CROSS_CHAIN_DETECTED, CROSS_CHAIN_ALL, CROSS_JUMP };

		core::pack::task::TaskFactoryOP task_factory_;
		core::Real distance_cutoff_;

		Mode filter_mode_;
		bool normalize_by_sasa_, ss_only_;

		core::Size jump_;
		std::string sym_dof_name_;
};

}
}
}

#endif
