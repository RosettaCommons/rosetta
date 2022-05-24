// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/Thioether_Helper.cc
/// @brief A CrosslinkerMoverHelper for placing a terminal (or lariat-style) thioether linkage.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/Thioether_Helper.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

// Protocols headers
#include <protocols/cyclic_peptide/crosslinker/thioether_util.hh>
#include <protocols/simple_moves/DeclareBond.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "protocols.cyclic_peptide.crosslinker.Thioether_Helper" );

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

// Maximum Nterm CA - CYS CB distance is 8.75 A; slight overestimate of true maximum (which is closer to 7.8 A).
// We're subtracting off the cysteine neighbour radius here:
#define MAX_CA_CB_DIST_MINUS_NBR 5.3027

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
Thioether_Helper::Thioether_Helper() = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
Thioether_Helper::Thioether_Helper( Thioether_Helper const & /*src*/ ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
Thioether_Helper::~Thioether_Helper() = default;

//////////////////////
/// Public Methods ///
//////////////////////

/// @brief Provide an opportunity to provide a citation for this crosslinker type.
/// @details The base class implementation does nothing.  This override indicates that this helper is
/// unpublished.
void
Thioether_Helper::provide_citation_info(
	basic::citation_manager::CitationCollectionList & citations
) const {
	using namespace basic::citation_manager;
	UnpublishedModuleInfoOP unpub(
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		"Thioether_Helper", CitedModuleType::CrosslinkerMoverHelper,
		"Vikram K. Mulligan", "Systems Biology Group, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org", "Added support for thioether cyclization to the CrosslinkerMover."
		)
	);
	unpub->add_author( "Andy Watkins", "Prescient Design", "andy.watkins2@gmail.com", "Added thioether crosslink to Rosetta database."  );
	citations.add(unpub);
}

/// @brief Given a pose and a selection of residues, add the linker,
/// align it crudely to the selected residues, and set up covalent bonds.
/// @details Calls add_linker_bonds_asymmetric().  Version for asymmetric poses.
void
Thioether_Helper::add_linker_asymmetric(
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	core::Size firstres, secondres; //Firstres has the backbone linkage, secondres has the sidechain linkage.
	get_thioether_residues_from_selection( selection, pose, firstres, secondres );
	set_up_thioether_variants( pose, firstres, secondres ); //Defined in utility file, thioether_util.hh/cc.
	TR << "Added thioether variant types." << std::endl;
	add_linker_bonds_asymmetric( pose, utility::vector1< core::Size >{ firstres, secondres }, 0 );
}

/// @brief Given a pose and a linker, add bonds between the linker and the residues
/// that coordinate the linker.
/// @details Can be called by add_linker_asymmetric().  Must be defined by derived
/// classes (pure virtual).  Version for asymmetric poses.
void
Thioether_Helper::add_linker_bonds_asymmetric(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & res_indices,
	core::Size const /*linker_index*/
) const {
	std::string const errmsg( "Error in Thioether_Helper::add_linker_bonds_asymmetric(): " );
	protocols::simple_moves::DeclareBond declbond;
	runtime_assert_string_msg( res_indices.size() == 2 && res_indices[1] != res_indices[2], errmsg + "Exactly two residues must be selected to set up linker bonds." );
	runtime_assert_string_msg( res_indices[1] > 0 && res_indices[1] <= pose.total_residue(), errmsg + "The first residue index must be within the range of pose residues (1-" + std::to_string(pose.total_residue()) + "), but got " + std::to_string(res_indices[1]) + " instead." );
	runtime_assert_string_msg( res_indices[2] > 0 && res_indices[2] <= pose.total_residue(), errmsg + "The second residue index must be within the range of pose residues (1-" + std::to_string(pose.total_residue()) + "), but got " + std::to_string(res_indices[2]) + " instead." );
	set_up_thioether_bond_mover( declbond, pose, res_indices[1], res_indices[2] );
	declbond.apply(pose);
	correct_thioether_virtuals( pose, res_indices[1], res_indices[2] );
	TR << "Added thioether bond." << std::endl;
}

/// @brief Given a pose and a selection of residues, add the linker,
/// align it crudely to the selected residues, and set up covalent bonds.
/// @details Version for symmetric poses.  However, since the linker itself has
/// no symmetry, we can just call the asymmetric version.
void
Thioether_Helper::add_linker_symmetric(
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	add_linker_asymmetric( pose, selection );
}

/// @brief Given a pose and a linker, add bonds between the linker and the residues
/// that coordinate the linker.
/// @details Can be called by add_linker_symmetric().  Must be defined by derived
/// classes (pure virtual).  Version for symmetric poses.
void
Thioether_Helper::add_linker_bonds_symmetric(
	core::pose::Pose & pose,
	core::Size const /*res1*/,
	core::Size const linker_index1,
	core::Size const linker_index2
) const {
	add_linker_bonds_asymmetric( pose, utility::vector1< core::Size >{ linker_index1, linker_index2 }, 0 );
}


/// @brief Given a selection of residues that have already been connected to a crosslinker,
/// add constraints for the crosslinker.
/// @details Must be defined by derived classes.  Version for asymmetric poses.
void
Thioether_Helper::add_linker_constraints_asymmetric(
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	core::Size firstres, secondres; //Firstres has the backbone linkage, secondres has the sidechain linkage.
	get_thioether_residues_from_selection( selection, pose, firstres, secondres );
	set_up_thioether_constraints( pose, firstres, secondres );
}

/// @brief Given a selection of residues that have already been connected to a crosslinker,
/// add constraints for the crosslinker.
/// @details Version for symmetric poses, but this just calls the asymmetric version.
void
Thioether_Helper::add_linker_constraints_symmetric(
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & selection,
	bool const //linker_was_added
) const {
	add_linker_constraints_asymmetric( pose, selection );
}

/// @brief Given indices of residues that are already linked to a linker, get the index
/// of the linker.
/// @details Returns zero for this class, since no linker is added.
core::Size
Thioether_Helper::get_linker_index_asymmetric(
	core::pose::Pose const & /*pose*/,
	utility::vector1< core::Size > const & /*res_indices*/
) const {
	return 0;
}

/// @brief Given indices of residues that are already linked to pieces of a linker, get
/// of the indices of the symmetric pieces of the linker.
/// @details Returns an empty vector for linker indices for this class, since no linker is added.
void
Thioether_Helper::get_linker_indices_symmetric(
	core::pose::Pose const & /*pose*/,
	utility::vector1< core::Size > const & /*res_indices*/,
	utility::vector1< core::Size > & linker_indices
) const {
	linker_indices.clear();
}

/// @brief Given a pose with residues selected to be linked by a linker, determine whether
/// the residues are too far apart.
/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
/// @note Higher values of the filter multiplier make it more permissive.
bool
Thioether_Helper::filter_by_sidechain_distance_asymmetric(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & selection,
	core::Real const & filter_multiplier
) const {
	core::Size ntermres, sidechainres;
	get_thioether_residues_from_selection( selection, pose, ntermres, sidechainres );

	core::Real const ref_dist( MAX_CA_CB_DIST_MINUS_NBR + pose.residue_type(sidechainres).nbr_radius() );
	core::Real const ref_dist_sq( ref_dist*ref_dist );
	core::Real const comparison_val( filter_multiplier*filter_multiplier*ref_dist_sq );

	core::Real const distsq( pose.residue(ntermres).xyz("CA").distance_squared( pose.residue(sidechainres).xyz("CB") ) );
	if ( distsq > comparison_val ) return true; //Filter fails if distance greater than maximum CB-CB separation.

	return false;
}

/// @brief Given a pose with residues selected to be linked by a linker, determine whether
/// the residues are too far apart.
/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.  This
/// version is for symmetric poses.
/// @note Higher values of the filter multiplier make it more permissive.
bool
Thioether_Helper::filter_by_sidechain_distance_symmetric(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & selection,
	core::Real const & filter_multiplier
) const {
	return filter_by_sidechain_distance_asymmetric( pose, selection, filter_multiplier );
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
/// @note Higher values of the filter multiplier make it more permissive.
bool
Thioether_Helper::filter_by_constraints_energy_asymmetric(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & selection,
	core::Real const & filter_multiplier
) const {
	return filter_by_constraints_energy( pose, selection, false, false, filter_multiplier );
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints
/// score.  This version is for symmetric poses.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
/// @note Higher values of the filter multiplier make it more permissive.
bool
Thioether_Helper::filter_by_constraints_energy_symmetric(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & selection,
	bool const linker_was_added,
	core::Real const & filter_multiplier
) const {
	return filter_by_constraints_energy( pose, selection, true, linker_was_added, filter_multiplier );
}

/********************************************
PRIVATE FUNCTIONS
*********************************************/

/// @brief Given a selection of exactly two residues and a pose, return the two thioether
/// residues.  Ensures that the residue that can make an N-terminal connection is ntermres,
/// and that the residue that has a side-chain thiol and can make a side-chain connection
/// is cysres.
void
Thioether_Helper::get_thioether_residues_from_selection(
	core::select::residue_selector::ResidueSubset const & selection,
	core::pose::Pose const & pose,
	core::Size & ntermres,
	core::Size & cysres
) const {
	std::string const errmsg( "Error in Thioether_Helper::get_thioether_residues_from_selection(): " );
	ntermres = 0;
	cysres = 0;

	debug_assert( selection.size() == pose.total_residue() );
	debug_assert( selection.size() >= 2 );

	for ( core::Size i(1), imax(selection.size()); i<=imax; ++i ) {
		if ( selection[i] ) {
			if ( ntermres == 0 ) {
				ntermres = i;
			} else if ( cysres == 0 ) {
				cysres = i;
			} else {
				utility_exit_with_message( errmsg + "More than two residues were selected!  A thioether connection is specified by exactly two residues." );
			}
		}
	}

	runtime_assert_string_msg( ntermres != 0, errmsg + "No suitable residue that can accept a thioether bond at its N-terminus was selected!" );
	runtime_assert_string_msg( cysres != 0, "No suitable residue that can form a thioether connection from its side-chain was selected!" );

	if ( (!pose.residue_type( cysres ).is_sidechain_thiol()) || (pose.residue( ntermres ).connected_residue_at_lower() != 0) ) {
		runtime_assert_string_msg( pose.residue_type( ntermres ).is_sidechain_thiol(), errmsg + "Neither of the selected residues has a sidechain thiol, as is required for a thioether connection." );
		runtime_assert_string_msg( pose.residue( cysres ).connected_residue_at_lower() == 0, errmsg + "Neither of the selected residues has a free N-terminus, as is required for a thioether connection." );
		std::swap( cysres, ntermres );
	}
}

} //crosslinker
} //cyclic_peptide
} //protocols
