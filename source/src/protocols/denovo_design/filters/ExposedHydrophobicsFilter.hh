// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/filters/ExposedHydrophobicsFilter.hh
/// @brief Tom's Denovo Protocol. This is freely mutable and used for playing around with stuff
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_protocols_denovo_design_filters_ExposedHydrophobicsFilter_hh
#define INCLUDED_protocols_denovo_design_filters_ExposedHydrophobicsFilter_hh

// Unit headers
#include <protocols/denovo_design/filters/ExposedHydrophobicsFilter.fwd.hh>

// Project headers
#include <protocols/filters/Filter.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

//// C++ headers
#include <string>

#include <core/io/silent/silent.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace denovo_design {
namespace filters {

class ExposedHydrophobicsFilter : public protocols::filters::Filter {
public:

	/// @brief Initialize ExposedHydrophobicsFilter
	ExposedHydrophobicsFilter();

	/// @brief virtual constructor to allow derivation
	virtual ~ExposedHydrophobicsFilter();

	/// @brief Parses the ExposedHydrophobicsFilter tags
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	/// @brief Return the name of this mover.
	virtual std::string get_name() const;

	/// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::filters::FilterOP clone() const;

	/// @brief Apply the ExposedHydrophobicsFilter. Overloaded apply function from filter base class.
	virtual bool apply( core::pose::Pose const & pose ) const;
	protocols::filters::FilterOP fresh_instance() const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;

private:   // options
	/// @brief If total calculated filter score is <= theshold_, the filter passes, otherwise it fails.
	core::Real threshold_;
	/// @brief If a residue has sasa <= sasa_cutoff_, it is considered buried and ignored
	core::Real sasa_cutoff_;

private:   // other data
	static std::string const hydrophobic_residues_;

};


} // filters
} // denovo_design
} // protocols

#endif
