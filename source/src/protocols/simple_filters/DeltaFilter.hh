// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/DeltaFilter.hh
/// @brief Reports the average degree of connectivity of interface residues
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_DeltaFilter_hh
#define INCLUDED_protocols_simple_filters_DeltaFilter_hh

// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/simple_filters/DeltaFilter.fwd.hh>

#include <utility/vector1.hh>

// Unit headers

namespace protocols {
namespace simple_filters {

class DeltaFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	DeltaFilter();
	///@brief Constructor with a single target residue
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~DeltaFilter();
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	void ref_baseline( core::Real const rb );
	core::Real baseline() const;
	void baseline( core::Real const b );
	bool lower() const;
	void lower( bool const l );
	bool upper() const;
	void upper( bool const u );
	void filter( protocols::filters::FilterOP filter );
	protocols::filters::FilterOP filter() const;
	core::Real range() const;
	void range( core::Real const r );
	bool unbound() const;
	void unbound( bool const u );
	bool relax_unbound() const;
	void relax_unbound( bool const u );
	bool changing_baseline() const;
	void changing_baseline( bool const c );
	core::Size jump() const;
	void jump( core::Size const j );
	protocols::moves::MoverOP relax_mover() const;
	void relax_mover( protocols::moves::MoverOP const mover );
private:

	protocols::filters::FilterOP filter_; //which filter to use
	core::Real baseline_; // dflt 0.0; the baseline against which to compare
	core::Real range_; // dflt 0.0; how much above/below baseline to allow
  bool lower_, upper_; // dflt false, true, respectively; use a lower/upper cutoff
	bool unbound_; //dflt false; evaluate the filter in the unbound state? If so, activate jump, below
	bool relax_unbound_; //dflt false; call relax mover on unbound pose?
	bool changing_baseline_; //dflt false; Should the baseline be reset each time this delta filter is called? If so, then set changing baseline to true.
	core::Size jump_; //dflt 0, but defaults to 1 if unbound is true
	protocols::moves::MoverOP relax_mover_; // a mover to be called before evaluating the filter's value. Only called for computing the baseline at the start!
	core::pose::PoseOP reference_pose_; //the reference pose that the baseline will be calculated on. note: this will only get set if a pose saved in the middle of an RS protocol and not the starting structure is the reference
	mutable core::Real ref_baseline_; // The baseline from the reference pose

	void unbind( core::pose::Pose & ) const; //utility function for unbinding the pose
};

} // simple_filters
} // protocols

#endif //INCLUDED_protocols_simple_filters_DeltaFilter_HH_
