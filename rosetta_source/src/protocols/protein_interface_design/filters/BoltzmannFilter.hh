// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/BoltzmannFilter.hh
/// @brief Reports the average degree of connectivity of interface residues
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_BoltzmannFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_BoltzmannFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/protein_interface_design/filters/BoltzmannFilter.fwd.hh>

//Auto Headers

// Unit headers

namespace protocols {
namespace protein_interface_design{
namespace filters {

class BoltzmannFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	BoltzmannFilter();
	///@brief Constructor with a single target residue
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~BoltzmannFilter();
	void parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	void add_positive_filter( protocols::filters::FilterOP f );
	void add_negative_filter( protocols::filters::FilterOP f );
	utility::vector1< protocols::filters::FilterOP > get_positive_filters() const;
	utility::vector1< protocols::filters::FilterOP > get_negative_filters() const;
	void temperature( core::Real const temp );
	core::Real temperature() const;
	core::Real fitness_threshold() const;
	void fitness_threshold( core::Real const f );
private:
	utility::vector1< protocols::filters::FilterOP > positive_filters_;
	utility::vector1< protocols::filters::FilterOP > negative_filters_;

	core::Real temperature_;
	core::Real fitness_threshold_;
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_BoltzmannFilter_HH_

