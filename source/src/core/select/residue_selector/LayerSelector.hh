// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/LayerSelector.hh
/// @brief  Header file for a ResidueSelector for selecting residues based
/// on layers defined by burial (core, boundary, surface, etc.).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_select_residue_selector_LayerSelector_HH
#define INCLUDED_core_select_residue_selector_LayerSelector_HH

// Unit headers
#include <core/select/residue_selector/LayerSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/select/util/SelectResiduesByLayer.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief The LayerSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions which match the given residue index. The index is read as comma-separated
/// list of either Rosetta indices (e.g. 10) or PDB numbers (e.g. 10A, residue 10 of chain A). Detection
/// and mapping from PDB to Rosetta residue numbers is done internally.
class LayerSelector : public ResidueSelector {
public:

	/// @brief Default constructor
	///
	LayerSelector();

	/// @brief Copy constructor
	///
	LayerSelector( LayerSelector const &src );

	/// @brief Destructor
	///
	~LayerSelector() override;

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	ResidueSelectorOP clone() const override;

	/// @brief Apply function: generate a ResidueSubset given a Pose.
	///
	ResidueSubset apply( core::pose::Pose const & pose ) const override;

	/// @brief Parse xml tag setting up this selector.
	///
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	/// @brief Get the class name.
	/// @details Calls class_name().
	std::string
	get_name() const override;

	/// @brief Get the class name.
	///
	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Set the layers that this selector will choose.
	///
	void set_layers( bool pick_core, bool pick_boundary, bool pick_surface );

	/// @brief Set the radius for the rolling ball algorithm used to determine burial.
	///
	void set_ball_radius( core::Real radius );

	/// @brief Set whether the sidechain neighbors algorithm is used to determine burial (as
	/// opposed to the rolling ball algorithm).
	void set_use_sc_neighbors( bool val );

	/// @brief Get whether the sidechain neighbors algorithm is used to determine burial (as
	/// opposed to the rolling ball algorithm).
	bool use_sc_neighbors() const;

	/// @brief Set whether to cache versus recompute the layer selection whenever it is accessed.
	void set_cache_selection( bool const val ) { cache_selection_ = val; }

	/// @brief Set the midpoint of the distance falloff if the sidechain neighbors method is used
	/// to define layers.
	void set_sc_neighbor_dist_midpoint( core::Real val );

	/// @brief Set the factor by which sidechain neighbor counts are divided if the sidechain
	/// neighbors method is used to define layers.
	void set_sc_neighbor_denominator( core::Real val );

	/// @brief Set the cutoffs for core and surface layers.
	/// @details Boundary is defined implicitly.  This can be a SASA cutoff or a neighbour count, depending
	/// on the algorithm.
	void set_cutoffs (core::Real core, core::Real surf);

	/// @brief Set the sidechain neighbor angle shift value.
	/// @details See the core::select::util::SelectResiduesByLayer class for details of the math.
	void set_angle_shift_factor( core::Real val );

	/// @brief Set the sidechain neighbor angle exponent.
	/// @details See the core::select::util::SelectResiduesByLayer class for details of the math.
	void set_angle_exponent( core::Real val );

	/// @brief Set the sidechain neighbor distance exponent.
	/// @details See the core::select::util::SelectResiduesByLayer class for details of the math.
	void set_dist_exponent( core::Real val );

private: // data members

	/// @brief Whether to cache the residue selection or recompute it when requested.
	bool cache_selection_;

	/// @brief Owning pointer to the calculator that determines what layer
	/// a residue is in.
	/// @details Object created and destroyed with LayerSelector objects.
	core::select::util::SelectResiduesByLayerOP srbl_;

#ifdef SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_select_residue_selector_LayerSelector )
#endif // SERIALIZATION

#endif
