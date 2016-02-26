// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/LayerSelector.cc
/// @brief  A ResidueSelector for selecting residues based on layers defined by burial (core, boundary, surface, etc.).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit headers
#include <core/select/residue_selector/LayerSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/select/residue_selector/util.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/select/util/SelectResiduesByLayer.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/exit.hh>

// C++ headers
#include <utility/assert.hh>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.LayerSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Default constructor
///
LayerSelector::LayerSelector() :
	cache_selection_( false ),
	srbl_( new core::select::util::SelectResiduesByLayer )
	//TODO -- initialize here
{}

/// @brief Copy constructor.
///
LayerSelector::LayerSelector( LayerSelector const &src ) :
	cache_selection_( src.cache_selection_ ),
	srbl_( src.srbl_->clone() ) //CLONE this -- don't copy it.
{}

/// @brief Destructor.
///
LayerSelector::~LayerSelector() {}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP LayerSelector::clone() const { return ResidueSelectorOP( new LayerSelector(*this) ); }

/// @brief Apply function: generate a ResidueSubset given a Pose.
///
ResidueSubset
LayerSelector::apply( core::pose::Pose const & pose ) const
{
	core::Size const nres( pose.n_residue() );

	// Return the cached value if requested by the user
	if( cache_selection_ && srbl_->is_selection_initialized() ) {
		ResidueSubset cached_subset( nres, false );
		const utility::vector1< Size >* layer_selections[3] = { &srbl_->selected_core_residues(), &srbl_->selected_boundary_residues(), &srbl_->selected_surface_residues() };
		for( int layer_idx= 0; layer_idx < 3; layer_idx++ ) {
			for( Size i = 1; i <= layer_selections[layer_idx]->size(); i++ ) {
				Size residue_index = (*layer_selections[layer_idx])[i];
				runtime_assert(residue_index <= nres);
				cached_subset[residue_index] = true;
			}
		}

		TR << "Returning cached selection result " << std::endl;
		return cached_subset;
	}

	utility::vector1<core::Size> selected_residues = srbl_->compute(pose, "", true);

	ResidueSubset subset( nres, false );
	for ( core::Size i=1, imax=selected_residues.size(); i<=imax; ++i ) {
		runtime_assert(selected_residues[i] > 0 && selected_residues[i] <= nres);
		subset[selected_residues[i]] = true;
	}

	return subset;
}

/// @brief Parse xml tag setting up this selector.
///
void
LayerSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &/*data*/)
{
	//Determine which layers we're picking:
	set_layers(
		tag->getOption<bool>("select_core", false),
		tag->getOption<bool>("select_boundary", false),
		tag->getOption<bool>("select_surface", false)
	);

	//Select algorithm:
	set_use_sc_neighbors( tag->getOption<bool>("use_sidechain_neighbors", true) );

	//Set algorithm details:
	if ( tag->hasOption("ball_radius") ) {
		set_ball_radius( tag->getOption<core::Real>("ball_radius") );
	}
	if ( tag->hasOption("sc_neighbor_dist_midpoint") ) {
		set_sc_neighbor_dist_midpoint( tag->getOption<core::Real>("sc_neighbor_dist_midpoint") );
	}
	if ( tag->hasOption("sc_neighbor_denominator") ) {
		set_sc_neighbor_denominator( tag->getOption<core::Real>("sc_neighbor_denominator") );
	}
	if ( tag->hasOption("sc_neighbor_angle_shift_factor") ) {
		set_angle_shift_factor( tag->getOption<core::Real>("sc_neighbor_angle_shift_factor") );
	}
	if ( tag->hasOption("sc_neighbor_angle_exponent") ) {
		set_angle_exponent( tag->getOption<core::Real>("sc_neighbor_angle_exponent") );
	}
	if ( tag->hasOption("sc_neighbor_dist_exponent") ) {
		set_dist_exponent( tag->getOption<core::Real>("sc_neighbor_dist_exponent") );
	}
	set_cutoffs(
		tag->getOption( "core_cutoff", (use_sc_neighbors() ? 5.2 : 20.0 ) ),
		tag->getOption( "surface_cutoff", (use_sc_neighbors() ? 2.0 : 40.0 ) )
	);
	cache_selection_ = tag->getOption<bool>("cache_selection", false);

	return;
}

/// @brief Get the class name.
/// @details Calls class_name().
std::string LayerSelector::get_name() const {
	return LayerSelector::class_name();
}

/// @brief Get the class name.
///
std::string
LayerSelector::class_name() {
	return "Layer";
}

void
LayerSelector::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute( "select_core",                    xs_boolean, "false" ));
	attributes.push_back( XMLSchemaAttribute( "select_boundary",                xs_boolean, "false" ));
	attributes.push_back( XMLSchemaAttribute( "select_surface",                 xs_boolean, "false" ));
	attributes.push_back( XMLSchemaAttribute( "cache_selection",                xs_boolean, "false" ));
	attributes.push_back( XMLSchemaAttribute( "use_sidechain_neighbors",        xs_boolean, "true" ));
	attributes.push_back( XMLSchemaAttribute( "ball_radius",                    xs_decimal ));
	attributes.push_back( XMLSchemaAttribute( "sc_neighbor_dist_midpoint",      xs_decimal ));
	attributes.push_back( XMLSchemaAttribute( "sc_neighbor_denominator",        xs_decimal ));
	attributes.push_back( XMLSchemaAttribute( "sc_neighbor_angle_shift_factor", xs_decimal ));
	attributes.push_back( XMLSchemaAttribute( "sc_neighbor_angle_exponent",     xs_decimal ));
	attributes.push_back( XMLSchemaAttribute( "sc_neighbor_dist_exponent",      xs_decimal ));
	attributes.push_back( XMLSchemaAttribute( "core_cutoff",                    xs_decimal ));
	attributes.push_back( XMLSchemaAttribute( "surface_cutoff",                 xs_decimal ));
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );
}


/// @brief Set the layers that this selector will choose.
///
void
LayerSelector::set_layers(
	bool const pick_core,
	bool const pick_boundary,
	bool const pick_surface
) {
	srbl_->set_design_layer(pick_core, pick_boundary, pick_surface);
	if ( TR.visible() ) {
		TR << "Setting core=" << (pick_core ? "true" : "false") << " boundary=" << (pick_boundary ? "true" : "false") << " surface=" << (pick_surface ? "true" : "false") << " in LayerSelector." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Set the radius for the rolling ball algorithm used to determine burial.
///
void
LayerSelector::set_ball_radius(
	core::Real const &radius
) {
	srbl_->pore_radius(radius);
	if ( TR.visible() ) {
		TR << "Setting radius for rolling ball algorithm to " << radius << " in LayerSelector.  (Note that this will have no effect if the sidechain neighbors method is used.)" << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Set whether the sidechain neighbors algorithm is used to determine burial (as
/// opposed to the rolling ball algorithm).
void LayerSelector::set_use_sc_neighbors( bool const val )
{
	srbl_->use_sidechain_neighbors( val );
	if ( TR.visible() ) {
		if ( val ) {
			TR << "Setting LayerSelector to use sidechain neighbors to determine burial." << std::endl;
		} else {
			TR << "Setting LayerSelector to use rolling ball-based occlusion to determine burial." << std::endl;
		}
		TR.flush();
	}
	return;
}

/// @brief Get whether the sidechain neighbors algorithm is used to determine burial (as
/// opposed to the rolling ball algorithm).
bool LayerSelector::use_sc_neighbors() const { return srbl_->use_sidechain_neighbors(); }

/// @brief Set the midpoint of the distance falloff if the sidechain neighbors method is used
/// to define layers.
void LayerSelector::set_sc_neighbor_dist_midpoint( core::Real const &val )
{
	srbl_->set_dist_midpoint(val);
	if ( TR.visible() ) {
		TR << "Set distance falloff midpoint for the LayerSelector to " << val << ".  Note that this has no effect if the sidechain neighbors method is not used." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Set the factor by which sidechain neighbor counts are divided if the sidechain
/// neighbors method is used to define layers.
void LayerSelector::set_sc_neighbor_denominator( core::Real const &val )
{
	srbl_->set_rsd_neighbor_denominator(val);
	if ( TR.visible() ) {
		TR << "Set denominator factor for the LayerSelector to " << val << ".  Note that this has no effect if the sidechain neighbors method is not used." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Set the cutoffs for core and surface layers.
/// @details Boundary is defined implicitly.  This can be a SASA cutoff or a neighbour count, depending
/// on the algorithm.
void LayerSelector::set_cutoffs (core::Real const &core, core::Real const &surf)
{
	srbl_->sasa_core(core);
	srbl_->sasa_surface(surf);
	if ( TR.visible() ) {
		TR << "Set cutoffs for core and surface to " << core << " and " << surf << ", respectively, in LayerSelector." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Set the sidechain neighbor angle shift value.
/// @details See the core::select::util::SelectResiduesByLayer class for details of the math.
void
LayerSelector::set_angle_shift_factor( core::Real const &val )
{
	srbl_->set_angle_shift_factor(val);
	if ( TR.visible() ) {
		TR << "Set angle shift value to " << val << " in LayerSelector.  Note that this has no effect if the sidechain neighbors method is not used." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Set the sidechain neighbor angle exponent.
/// @details See the core::select::util::SelectResiduesByLayer class for details of the math.
void
LayerSelector::set_angle_exponent( core::Real const &val )
{
	srbl_->set_angle_exponent(val);
	if ( TR.visible() ) {
		TR << "Set angle exponent to " << val << " in LayerSelector.  Note that this has no effect if the sidechain neighbors method is not used." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Set the sidechain neighbor distance exponent.
/// @details See the core::select::util::SelectResiduesByLayer class for details of the math.
void
LayerSelector::set_dist_exponent( core::Real const &val )
{
	srbl_->set_dist_exponent(val);
	if ( TR.visible() ) {
		TR << "Set distance exponent to " << val << " in LayerSelector.  Note that this has no effect if the sidechain neighbors method is not used." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Return an owning pointer to a new LayerSelector object.
///
ResidueSelectorOP
LayerSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new LayerSelector );
}

/// @brief Get the class name.
/// @details Also calls class_name.
std::string
LayerSelectorCreator::keyname() const {
	return LayerSelector::class_name();
}

void
LayerSelectorCreator::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) const {
	LayerSelector::provide_selector_xsd( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core
