// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/excn/EXCN_Base.hh>
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
	typedef core::scoring::sc::ShapeComplementarityCalculator ShapeComplementarityCalculator;
	typedef core::scoring::sc::RESULTS ShapeComplementarityCalculatorResults;

public:// constructor/destructor
	// @brief default constructor
	ShapeComplementarityFilter();

	// @brief constructor with arguments
	ShapeComplementarityFilter( Real const & filtered_sc, Real const & filtered_area,
		Size const & jump_id, Size const & quick, Size const & verbose);

	~ShapeComplementarityFilter() override= default;

public:// virtual constructor
	// @brief make clone
	filters::FilterOP clone() const override { return filters::FilterOP( new ShapeComplementarityFilter( *this ) ); }

	// @brief make fresh instance
	filters::FilterOP fresh_instance() const override { return filters::FilterOP( new ShapeComplementarityFilter() ); }

public:// accessor
	// @brief get name of this filter
	// XRW TEMP  std::string name() const override { return "ShapeComplementarity"; }

public:// mutator
	void filtered_sc( Real const & filtered_sc );
	void filtered_area( Real const & filtered_area );
	void jump_id( Size const & jump_id );
	void quick( Size const & quick );
	void verbose( Size const & verbose );
	void residues1( std::string const & res_string );
	void residues2( std::string const & res_string );
	void sym_dof_name( std::string const & sym_dof_name );
	std::string sym_dof_name() const;
	void write_int_area( bool write_int_area );
	bool write_int_area( ) const;
	void multicomp( bool multicomp );
	bool multicomp( ) const;

public:// parser
	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const & filters,
		Movers_map const & movers,
		Pose const & pose ) override;

public:// virtual main operation
	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	bool apply( Pose const & pose ) const override;

	/// @brief
	Real report_sm( Pose const & pose ) const override;

public:
	/// @brief calc shape complementarity, returns results of the ShapeComplementarityCalculator
	/// @param[in] pose Pose to be analyzed
	/// @exception EXCN_CalcInitFailed Thrown if calculator couldn't be initialized
	/// @exception EXCN_ResultsInvalid Thrown if computed results are invalid
	/// @returns ShapeComplementarityCalculator::RESULTS object
	ShapeComplementarityCalculatorResults
	compute( Pose const & pose ) const;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	/// @brief Uses residue selectors to set up the ShapeComplementarityCalculator
	/// @param[in]  pose Pose to be analyzed
	/// @param[out] scc Initialized, empty ShapeComplementarityCalculator, to which pose residues are added
	void
	setup_from_selectors( Pose const & pose, ShapeComplementarityCalculator & scc ) const;

	/// @brief Uses multi-component symmetric interfaces to set up the ShapeComplementarityCalculator
	/// @param[in]  pose              Pose to be analyzed
	/// @param[out] scc               Initialized, empty ShapeComplementarityCalculator, to which pose residues are added
	/// @param[out] nsubs_scalefactor Writes number of subunits, to be used as scaling factor for sc calculations
	void
	setup_multi_component_symm(
		Pose const & pose,
		ShapeComplementarityCalculator & scc,
		core::Real & nsubs_scalefactor ) const;

	/// @brief Uses single-component symmetric interfaces to set up the ShapeComplementarityCalculator
	/// @param[in]  pose              Pose to be analyzed
	/// @param[out] scc               Initialized, empty ShapeComplementarityCalculator, to which pose residues are added
	/// @param[out] nsubs_scalefactor Writes number of subunits, to be used as scaling factor for sc calculations
	void
	setup_single_component_symm(
		Pose const & pose,
		ShapeComplementarityCalculator & scc,
		core::Real & nsubs_scalefactor ) const;

	/// @brief prints results to given tracer in a human-readable format
	/// @param[out] tr std::ostream object to write to
	/// @param[in]  r  ShapeComplementarityCalculatorResults object containing results
	void
	print_sc_results(
		std::ostream & tr,
		ShapeComplementarityCalculatorResults const & r,
		core::Real const nsubs_scalefactor ) const;

	/// @brief writes area value to current jd2 job
	/// @param[in] pose     Pose being analyzed
	/// @param[in] area_val Area to be reported, before correcting for symmetry.  If area < 0,
	///                     the uncorrected value will be reported. If the pose isn't symmetric,
	///                     the uncorrected value will be reported.
	void
	write_area( Pose const & pose, core::Real const area_val ) const;

private:
	Real filtered_sc_;
	Real filtered_area_;
	Size jump_id_;
	Size quick_;
	Size verbose_;
	core::select::residue_selector::ResidueSelectorCOP selector1_;
	core::select::residue_selector::ResidueSelectorCOP selector2_;
	bool write_int_area_;

	// symmetry-specific
	bool multicomp_;
	std::string sym_dof_name_;
};

/// @brief Super-simple exception to be thrown when we can't initialize the SC calculator
class EXCN_InitFailed : public utility::excn::EXCN_Base {
public:
	void show( std::ostream & ) const override {}
};

/// @brief Super-simple exception to be thrown when the SC calculator fails to compute
class EXCN_CalcFailed : public utility::excn::EXCN_Base {
public:
	void show( std::ostream & ) const override {}
};

} // filters
} // protocols

#endif
