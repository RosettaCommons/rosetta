// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/ShapeComplementarityFilter.hh
/// @brief  header file for ShapeComplementarityFilter class
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)


#ifndef INCLUDED_protocols_simple_filters_ShapeComplementarityFilter_hh
#define INCLUDED_protocols_simple_filters_ShapeComplementarityFilter_hh

// Unit Headers
#include <protocols/simple_filters/ShapeComplementarityFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace simple_filters {

class ShapeComplementarityFilter : public protocols::filters::Filter {
public:
	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;

public:// constructor/destructor
	// @brief default constructor
	ShapeComplementarityFilter();

	// @brief constructor with arguments
	ShapeComplementarityFilter( Real const & filtered_sc, Real const & filtered_area,
		Size const & jump_id, Size const & quick, Size const & verbose);

	// @brief copy constructor
	ShapeComplementarityFilter( ShapeComplementarityFilter const & rval );

	virtual ~ShapeComplementarityFilter(){}

public:// virtual constructor
	// @brief make clone
	virtual filters::FilterOP clone() const { return new ShapeComplementarityFilter( *this ); }

	// @brief make fresh instance
	virtual filters::FilterOP fresh_instance() const { return new ShapeComplementarityFilter(); }

public:// accessor
	// @brief get name of this filter
	virtual std::string name() const { return "ShapeComplementarity"; }

public:// mutator
	void filtered_sc( Real const & filtered_sc );
	void filtered_area( Real const & filtered_area );
	void jump_id( Size const & jump_id );
	void quick( Size const & quick );
	void verbose( Size const & verbose );
	void residues1( utility::vector1< core::Size > const & residues );
	void residues2( utility::vector1< core::Size > const & residues );
	void sym_dof_name( std::string const & sym_dof_name );
	std::string sym_dof_name() const;
	void write_int_area( bool write_int_area );
	bool write_int_area( ) const;
	void multicomp( bool multicomp );
	bool multicomp( ) const;

public:// parser
	virtual void parse_my_tag( TagCOP const tag,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		Movers_map const &,
		Pose const & );
	void parse_def( utility::lua::LuaObject const & def,
					utility::lua::LuaObject const & score_fxns,
					utility::lua::LuaObject const & tasks );

public:// virtual main operation
	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	virtual bool apply( Pose const & pose ) const;

	/// @brief
	virtual Real report_sm( Pose const & pose ) const;

	/// @brief calc shape complementarity
	virtual Size compute( Pose const & pose ) const;

private:
	core::scoring::sc::ShapeComplementarityCalculator mutable scc_;

	Real filtered_sc_;
	Real filtered_area_;
	Size jump_id_;
	Size quick_;
	Size verbose_;
	utility::vector1<core::Size> residues1_;
	utility::vector1<core::Size> residues2_;
	bool write_int_area_;

	// symmetry-specific
	bool multicomp_;
	std::string sym_dof_name_;
};

} // filters
} // protocols

#endif
