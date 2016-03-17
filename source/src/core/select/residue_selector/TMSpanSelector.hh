// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/TMSpanSelector.hh
/// @brief  Select residues within given transmembrane spans in a membrane protein
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_TMSpanSelector_HH
#define INCLUDED_core_select_residue_selector_TMSpanSelector_HH

// Unit headers
#include <core/select/residue_selector/TMSpanSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.hh> 
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace core {
namespace select {
namespace residue_selector {

/// @brief Select all transmembrane spans in a membrane protein
class TMSpanSelector : public core::select::residue_selector::ResidueSelector {

public:

	TMSpanSelector();
	
	TMSpanSelector(
		bool all_tm,
		utility::vector1< core::Size > tm_spans
	);
	
	virtual ~TMSpanSelector();

	virtual core::select::residue_selector::ResidueSubset apply( core::pose::Pose const & pose ) const;
	
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	virtual
	std::string
	get_name() const;
	
	virtual ResidueSelectorOP clone() const;

	static std::string class_name();
	static void provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd );
	
	/// @brief Check that the transmembrane spans selected from input are present
	/// in the pose's spanning topology object
	void check_valid_tm_selection( core::pose::Pose const & pose ) const;

public: // Accessors & Mutators

	/// @brief Use all transmembrane spans defined in the pose?
	void all( bool select_all ) { all_ = select_all; }
	bool all() const { return all_; }
	
	/// @brief Use the following transmembrane spans defined in the pose
	void add_span( core::Size span_no ) { tm_spans_.push_back( span_no ); }
	utility::vector1< core::Size > all_tm_spans() const { return tm_spans_; }

private:

	// Select all transmembrane spans
	bool all_;
	
	// Transmembrane spans to include (all or some)
	utility::vector1< core::Size > tm_spans_;
	
};


} //core
} //select
} //residue_selector

#endif //INCLUDED_core_select_residue_selector_TMSpanSelector_hh
