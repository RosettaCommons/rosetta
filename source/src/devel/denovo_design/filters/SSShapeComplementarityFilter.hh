// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/filters/SSShapeComplementarityFilter.hh
/// @brief Tom's Denovo Protocol. This is freely mutable and used for playing around with stuff
/// @detailed
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_devel_denovo_design_filters_SSShapeComplementarityFilter_hh
#define INCLUDED_devel_denovo_design_filters_SSShapeComplementarityFilter_hh

// Unit headers
#include <devel/denovo_design/filters/SSShapeComplementarityFilter.fwd.hh>

// Project headers
#include <protocols/filters/Filter.hh>
#include <protocols/fldsgn/topology/HelixPairing.fwd.hh>
#include <protocols/fldsgn/topology/HSSTriplet.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/sc/MolecularSurfaceCalculator.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.fwd.hh>

//// C++ headers
#include <string>

#include <core/io/silent/silent.fwd.hh>
#include <utility/vector1.hh>



namespace devel {
namespace denovo_design {
namespace filters {

class SSShapeComplementarityFilter : public protocols::filters::Filter {
public:

  /// @brief Initialize SSShapeComplementarityFilter
  SSShapeComplementarityFilter();

  /// @brief copy constructor
  SSShapeComplementarityFilter( SSShapeComplementarityFilter const & rval );

  /// @brief virtual constructor to allow derivation
	virtual ~SSShapeComplementarityFilter();

  /// @brief Parses the SSShapeComplementarityFilter tags
	void parse_my_tag(
	  utility::tag::TagCOP tag,
	  basic::datacache::DataMap & data,
	  protocols::filters::Filters_map const &,
	  protocols::moves::Movers_map const &,
	  core::pose::Pose const & );

  /// @brief Return the name of this mover.
  virtual std::string get_name() const;

  /// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::filters::FilterOP clone() const;

  /// @brief Apply the SSShapeComplementarityFilter. Overloaded apply function from filter base class.
	virtual protocols::filters::FilterOP fresh_instance() const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual bool apply( core::pose::Pose const & pose ) const;

	core::Real compute( core::pose::Pose const & pose ) const;

private:   // private functions
	/// @brief sets up the underlying filter to work based on a helix
	void
	setup_sc( core::pose::Pose const & pose,
						protocols::fldsgn::topology::SS_BaseCOP const ss ) const;

	/// @brief sets up the underlying shapecomplementarity filter to work based on secondary structure elements
	void
	setup_sc_hss( core::pose::Pose const & pose,
								protocols::fldsgn::topology::SS_Info2 const & ss_info,
								protocols::fldsgn::topology::HSSTripletCOP hss_triplet ) const;
	void
	setup_sc_hh( core::pose::Pose const & pose,
							 protocols::fldsgn::topology::SS_Info2 const & ss_info,
							 protocols::fldsgn::topology::HelixPairingCOP helix_pair ) const;

	/// @brief Runs the SC calculator to obtain an SC score and an interaction area. Returns a result in the format core::scoring::sc::RESULTS.  Assumes the SC calculator has been initialized and has the correct residues added.
	core::scoring::sc::RESULTS const &
	get_sc_and_area() const;

private:   // options
	/// @brief controls outputtting verbose information about SC
	bool verbose_;
	/// @brief should we calculate SC from each loop to the rest of the protein?
	bool calc_loops_;
	/// @brief should we calculate SC from each helix to the rest of the protein?
	bool calc_helices_;

private:   // other data
	/// @brief the blueprint file that contains secondary structure definitions
	protocols::jd2::parser::BluePrintCOP blueprint_;
	/// @brief the shape complementarity calculator
	core::scoring::sc::ShapeComplementarityCalculatorOP scc_;

};



} // filters
} // denovo_design
} // devel

#endif
