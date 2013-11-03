// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/AtomicDistanceFilter.hh
/// @brief Filter for looking at specific atom distances
/// @author Rocco Moretti (rmoretti@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_AtomicDistanceFilter_hh
#define INCLUDED_protocols_simple_filters_AtomicDistanceFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/ResId.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

/// @brief detects atomic contacts between two atoms of two residues
class AtomicDistanceFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	AtomicDistanceFilter();
	AtomicDistanceFilter( core::Size const res1, core::Size const res2, std::string atom_desig1="CB", std::string atom_desig2="CB", bool as_type1=false, bool as_type2=false, core::Real const distance=4.0);
	virtual bool apply( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual filters::FilterOP clone() const {
		return new AtomicDistanceFilter( *this );
	}
	virtual filters::FilterOP fresh_instance() const{
		return new AtomicDistanceFilter();
	}

	virtual ~AtomicDistanceFilter(){};
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
private:
	core::Size residue1_, residue2_;
	std::string atomdesg1_, atomdesg2_;
	bool astype1_, astype2_; // if desg should be interpreted as an atomtype
	core::Real distance_;
};

} // filters
} // protocols

#endif
