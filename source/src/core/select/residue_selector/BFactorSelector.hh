// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/BFactorSelector.hh
/// @brief  A residue selector dependent on b-factor values
/// @author AmeyaHarmalkar (harmalkar.ameya24@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_BFactorSelector_HH
#define INCLUDED_core_select_residue_selector_BFactorSelector_HH

// Unit headers
#include <core/select/residue_selector/BFactorSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief A residue selector dependent on b-factor values
class BFactorSelector : public ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	BFactorSelector();

	/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
	BFactorSelector(
		core::Size const lower_res,
		core::Size const upper_res,
		core::select::residue_selector::ResidueSelectorCOP const residue_selector,
		bool cross_chain_boundaries = false );

public:

	/// @brief Destructor.
	~BFactorSelector() override;

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.

	ResidueSelectorOP clone() const override;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").

	ResidueSubset apply( core::pose::Pose const & pose ) const override;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	/// @brief Get the mover class name.
	std::string
	get_name() const override;

	/// @brief Get the mover class name.
	static std::string
	class_name();

	//unit-specific
	void set_lower_residues( core::Size const nres );
	void set_upper_residues( core::Size const nres );
	void set_min_contiguous_res( core::Size const nres );
	void set_lower_bfactor_threshold( core::Real const lower_thresh );
	void set_upper_bfactor_threshold( core::Real const upper_thresh );
	void set_cross_chain_boundaries( bool cross );

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief This residue selector is unpublished.  It returns AmeyaHarmalkar as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

private:
	core::Real lower_threshold_ = 0.0 ;
	core::Real upper_threshold_ = 100.0;
	core::Size contiguous_res_ = 1;
	core::Size lower_res_;
	core::Size upper_res_;
	ResidueSelectorCOP residue_selector_;
	bool cross_chain_boundaries_;



#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


};


} //residue_selector
} //select
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_select_residue_selector_BFactorSelector )
#endif // SERIALIZATION

#endif //INCLUDEDcore_select_residue_selector_BFactorSelector_HH
