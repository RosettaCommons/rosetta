// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/BindingStrainFilter.hh
/// @brief Reports the rotameric energy strain upon binding, by separating the monomers, repacking, minimizing, and measuring the energy difference from the bound but separated state
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_BindingStrainFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_BindingStrainFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/protein_interface_design/filters/BindingStrainFilter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>

// Unit headers

namespace protocols {
namespace protein_interface_design{
namespace filters {

class BindingStrainFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	BindingStrainFilter();
	/// @brief Constructor with a single target residue
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~BindingStrainFilter();
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	core::Size jump() const;
	void jump( core::Size const j );
	protocols::moves::MoverOP relax_mover() const;
	void relax_mover( protocols::moves::MoverOP const mover );
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task );
	core::Real threshold() const;
	void threshold( core::Real const t );
private:
	core::pack::task::TaskFactoryOP task_factory_; // what to repack?
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::moves::MoverOP relax_mover_; // potentially useful for minimizing after separation
	core::Size jump_; // dflt 1; which jump to separate
	core::Real threshold_; // dflt 3; how much strain to allow?
	void unbind( core::pose::Pose & ) const; //utility function for unbinding the pose
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_BindingStrainFilter_HH_
