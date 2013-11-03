// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/HbondsToResidueFilter.hh
/// @brief definition of filter classes for iterations of docking/design.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_HbondsToResidueFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_HbondsToResidueFilter_hh


// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>


// C++ headers

// Unit headers
//#include <basic/datacache/DataMap.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

using protocols::filters::Filter;
using protocols::filters::FilterOP;
using protocols::filters::Filters_map;

/// @brief returns true if the number of hbonding partners to a particular residue exceeds a certain value
/// This filter is useful in conjunction with DesignMinimizeHbonds class
class HbondsToResidueFilter : public Filter
{
public:
	typedef core::Real Real;
	typedef core::Size Size;
public :
	HbondsToResidueFilter() : Filter( "HbondsToResidue" ) {}
	HbondsToResidueFilter( Size const resnum, Size const partners, Real const energy_cutoff=-0.5,
						   bool const backbone=false, bool const sidechain=true ) : Filter( "HbondsToResidue" ) {
		resnum_ = resnum; partners_ = partners; energy_cutoff_ = energy_cutoff; backbone_ = backbone;
		sidechain_ = sidechain;
		runtime_assert( backbone_ || sidechain_ );
		runtime_assert( partners_ );
		runtime_assert( energy_cutoff_ <= 0 );
	}
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new HbondsToResidueFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new HbondsToResidueFilter();
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Size compute( core::pose::Pose const & pose ) const;
	virtual ~HbondsToResidueFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	Size resnum_, partners_;
	Real energy_cutoff_;
	bool backbone_, sidechain_, bb_bb_;
};

}
} // protein_interface_design
} // devel


#endif /*INCLUDED_DOCK_DESIGN_FILTERS_H_*/
