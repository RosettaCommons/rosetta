// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/RelativePoseFilter.hh
/// @brief Computes a filter's value on a pose that is modified from one that is read from disk. Useful for computing values for the same sequence across many structures.
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_RelativePoseFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_RelativePoseFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/protein_interface_design/filters/RelativePoseFilter.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

//Auto Headers

// Unit headers

namespace protocols {
namespace protein_interface_design{
namespace filters {

class RelativePoseFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	RelativePoseFilter();
	///@brief Constructor with a single target residue
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~RelativePoseFilter();
	void parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	void filter( protocols::filters::FilterOP filter );
	protocols::filters::FilterOP filter() const;
	protocols::moves::MoverOP relax_mover() const;
	void relax_mover( protocols::moves::MoverOP const mover );
	core::pose::PoseOP pose() const;
	void pose( core::pose::PoseOP p );
	std::string dump_pose_fname() const;
	void dump_pose_fname( std::string const s );

	core::scoring::ScoreFunctionOP scorefxn() const;
  void scorefxn( core::scoring::ScoreFunctionOP const s );
	void packing_shell( core::Real const s );
	core::Real packing_shell() const;
	core::pose::PoseOP thread_seq( core::pose::Pose const & p) const; //key functionality of this filter: thread relevant sections, repack and relax
	bool thread() const;
	void thread( bool const t );
private:
	protocols::filters::FilterOP filter_; //which filter to use
	protocols::moves::MoverOP relax_mover_; // a mover to be called before evaluating the filter's value.
	std::map< core::Size, core::Size > alignment_; //alignment of active pose to the pose read from disk. Only the segments that are threaded need be specified. Assumes a mapping from disk sequence to input structure sequence
	std::string dump_pose_fname_; //filename of dumped pose. Empty means no dumping
	core::pose::PoseOP pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real packing_shell_; //dflt 8; after threading, what shell to repack around residue
	bool thread_;// dflt true; should we thread or only repack? Only repack is useful for computing baseline filter values at the start of the run, vs. the values for the actual mutated sequences later
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_RelativePoseFilter_HH_

