// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJD_SASASummary.cc
/// @brief A class for storing information about solvent-exposed surface area and for transmitting it up the MPI hierarchy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_SASASummary.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.cyclic_peptide_predict.HierarchicalHybridJD_SASASummary" );


namespace protocols {
namespace cyclic_peptide_predict {

#define MIN_SASA 1.0e-10

/// @brief Options constructor.
HierarchicalHybridJD_SASASummary::HierarchicalHybridJD_SASASummary(
	int const originating_node_MPI_rank,
	core::Size const jobindex_on_originating_node,
	core::Real const &sasa,
	core::Real const &polar_sasa,
	core::Real const &hydrophobic_sasa
) :
	HierarchicalHybridJD_ResultsSummaryBase( originating_node_MPI_rank, jobindex_on_originating_node )
{
	set( sasa, polar_sasa, hydrophobic_sasa );
}

/// @brief Destructor.
HierarchicalHybridJD_SASASummary::~HierarchicalHybridJD_SASASummary(){}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
HierarchicalHybridJD_ResultsSummaryBaseOP
HierarchicalHybridJD_SASASummary::clone() const {
	return utility::pointer::make_shared< HierarchicalHybridJD_SASASummary >( *this );
}

/// @brief Set total sasa, polar sasa, and hydrophobic sasa.
/// @details Computes fraction_polar_sasa_ and fraction_hydrophobic_sasa_.
void
HierarchicalHybridJD_SASASummary::set(
	core::Real const &sasa,
	core::Real const &polar_sasa,
	core::Real const &hydrophobic_sasa
) {
	runtime_assert_string_msg( sasa > MIN_SASA, "Error in HierarchicalHybridJD_SASASummary::set(): The SASA cannot be zero or negative." );
	runtime_assert_string_msg( polar_sasa > 0, "Error in HierarchicalHybridJD_SASASummary::set(): The polar SASA cannot be negative." );
	runtime_assert_string_msg( hydrophobic_sasa > 0, "Error in HierarchicalHybridJD_SASASummary::set(): The hydrophobic SASA cannot be negative." );
	runtime_assert_string_msg( polar_sasa <= sasa, "Error in HierarchicalHybridJD_SASASummary::set(): The polar SASA cannot be greater than the total sasa." );
	runtime_assert_string_msg( hydrophobic_sasa <= sasa, "Error in HierarchicalHybridJD_SASASummary::set(): The hydrophobic SASA cannot be greater than the total sasa." );
	sasa_ = sasa;
	polar_sasa_ = polar_sasa;
	hydrophobic_sasa_ = hydrophobic_sasa;
	fraction_hydrophobic_sasa_ = hydrophobic_sasa / sasa;
	fraction_polar_sasa_ = polar_sasa / sasa;
}

} //cyclic_peptide_predict
} //protocols
