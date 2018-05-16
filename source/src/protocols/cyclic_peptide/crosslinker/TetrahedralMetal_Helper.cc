// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/TetrahedralMetal_Helper.cc
/// @brief A helper class for setting up tetrahedrally-coordinated metals
/// like Zn or Cu.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/TetrahedralMetal_Helper.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Protocols headers
#include <protocols/cyclic_peptide/CreateDistanceConstraint.hh>
#include <protocols/cyclic_peptide/CreateAngleConstraint.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static basic::Tracer TR( "protocols.cyclic_peptide.crosslinker.TetrahedralMetal_Helper" );

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

#define MAX_CB_CB_DIST_SQ 144.0

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
TetrahedralMetal_Helper::TetrahedralMetal_Helper(
	std::string const &metal_name_in
) :
	CrosslinkerMoverHelper(),
	metal_type_( TMH_Zn )
	//TODO initialize data here
{
	set_metal_type_from_name(metal_name_in);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
TetrahedralMetal_Helper::TetrahedralMetal_Helper( TetrahedralMetal_Helper const &src ) :
	CrosslinkerMoverHelper( src ),
	metal_type_( src.metal_type_ )
	//TODO copy data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
TetrahedralMetal_Helper::~TetrahedralMetal_Helper(){}

///////////////////////
/// Public Methods  ///
///////////////////////

/// @brief Given a pose and a selection of four residues, add the linker,
/// align it crudely to the selected residues, and set up covalent bonds.
/// @details Must be defined by derived classes.  Version for asymmetric poses.
void TetrahedralMetal_Helper::add_linker_asymmetric(
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	utility::vector1< core::Size > res_indices;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices );
	check_residue_indices_valid( res_indices, pose );

	for ( core::Size i(1), imax(res_indices.size()); i<=imax; ++i ) {
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_METAL_CONJUGATION, res_indices[i] );
		set_metal_bond_length( pose, res_indices[i] );
	}

	pose.update_residue_neighbors();
}

/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
/// @details Can be called by add_linker_asymmetric().  Must be defined by derived classes (pure virtual).  Version for asymmetric poses.
void
TetrahedralMetal_Helper::add_linker_bonds_asymmetric(
	core::pose::Pose &/*pose*/,
	utility::vector1< core::Size > const &/*res_indices*/,
	core::Size const /*linker_index*/
) const {
	/*Does nothing.  No bonds are needed, since the metals are virtual.*/
}

/// @brief Given a pose and a selection of four residues, add the linker,
/// align it crudely to the selected residues, and set up covalent bonds.
/// @details Must be defined by derived classes.  Version for symmetric poses.
void
TetrahedralMetal_Helper::add_linker_symmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	check_compatible_symmetry_type();
	runtime_assert_string_msg( core::pose::symmetry::is_symmetric( pose ), "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::add_linker_symmetric(): The pose is not symmetric." );
	runtime_assert_string_msg( CrosslinkerMoverHelper::selection_is_symmetric( selection, pose, symm_subunits_expected() ), "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::add_linker_symmetric(): The selection is not symmetric." );

	utility::vector1< core::Size > res_indices;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices );
	check_residue_indices_valid( res_indices, pose );

	core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const>( pose.conformation_ptr() ) );
	debug_assert( symmconf != nullptr ); //Should be true.
	core::conformation::symmetry::SymmetryInfoCOP symminfo( symmconf->Symmetry_Info() );
	debug_assert( symminfo != nullptr ); //Should also be true.

	for ( core::Size i(1), imax(res_indices.size()); i<=imax; ++i ) {
		if ( symminfo->bb_is_independent( res_indices[i] ) ) {
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_METAL_CONJUGATION, res_indices[i] );
			set_metal_bond_length( pose, res_indices[i] );
		}
	}

	pose.update_residue_neighbors();
}

/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
/// @details Can be called by add_linker_symmetric().  Must be defined by derived classes (pure virtual).  Version for symmetric poses.
void
TetrahedralMetal_Helper::add_linker_bonds_symmetric(
	core::pose::Pose &/*pose*/,
	core::Size const /*res1*/,
	core::Size const /*linker_index1*/,
	core::Size const /*linker_index2*/
) const {
	// GNDN
}

/// @brief Given a selection of four residues that have already been connected to a crosslinker,
/// add constraints for the crosslinker.
/// @details Must be defined by derived classes.  Version for asymmetric poses.
void
TetrahedralMetal_Helper::add_linker_constraints_asymmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	utility::vector1< core::Size > res_indices;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices );
	check_residue_indices_valid( res_indices, pose );

	for ( core::Size i(2), imax(res_indices.size()); i<=imax; ++i ) {
		for ( core::Size j(1); j<i; ++j ) {
			{ //Tether virtuals together:
				protocols::cyclic_peptide::CreateDistanceConstraint distcst;
				utility::vector1< core::Size > r1(1), r2(1);
				r1[1] = res_indices[i];
				r2[1] = res_indices[j];
				utility::vector1< std::string > a1(1), a2(1);
				a1[1] = a2[1] = std::string("VM1");
				utility::vector1< std::string > cstfxn(1);
				cstfxn[1] = std::string("HARMONIC 0.0 0.02");
				distcst.set(r1, a1, r2, a2, cstfxn);
				distcst.apply(pose);
			}
			{ //Angle constraints:
				add_angle_constraints(pose, res_indices, i, j);
				//add_angle_constraints(pose, res_indices, j, i);
			}
		}
	}
}

/// @brief Given a selection of four residues that have already been connected to a crosslinker,
/// add constraints for the crosslinker.
/// @details Must be defined by derived classes.  Version for symmetric poses.
void
TetrahedralMetal_Helper::add_linker_constraints_symmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	bool const /*linker_was_added*/
) const {
	check_compatible_symmetry_type();
	runtime_assert_string_msg( core::pose::symmetry::is_symmetric( pose ), "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::add_linker_symmetric(): The pose is not symmetric." );
	runtime_assert_string_msg( CrosslinkerMoverHelper::selection_is_symmetric( selection, pose, symm_subunits_expected() ), "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::add_linker_symmetric(): The selection is not symmetric." );

	utility::vector1< core::Size > res_indices;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices );
	check_residue_indices_valid( res_indices, pose );

	core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const>( pose.conformation_ptr() ) );
	debug_assert( symmconf != nullptr ); //Should be true.
	core::conformation::symmetry::SymmetryInfoCOP symminfo( symmconf->Symmetry_Info() );
	debug_assert( symminfo != nullptr ); //Should also be true.

	for ( core::Size i(2), imax(res_indices.size()); i<=imax; ++i ) {
		for ( core::Size j(1); j<i; ++j ) {
			if ( !symminfo->bb_is_independent(res_indices[j]) && !symminfo->bb_is_independent(res_indices[i]) ) continue;
			{ //Tether virtuals together:
				protocols::cyclic_peptide::CreateDistanceConstraint distcst;
				utility::vector1< core::Size > r1(1), r2(1);
				r1[1] = res_indices[i];
				r2[1] = res_indices[j];
				utility::vector1< std::string > a1(1), a2(1);
				a1[1] = a2[1] = std::string("VM1");
				utility::vector1< std::string > cstfxn(1);
				cstfxn[1] = std::string("HARMONIC 0.0 0.02");
				distcst.set(r1, a1, r2, a2, cstfxn);
				distcst.apply(pose);
			}
			{ //Angle constraints:
				add_angle_constraints(pose, res_indices, i, j);
				//add_angle_constraints(pose, res_indices, j, i);
			}
		}
	}
}

/// @brief Given indices of four residues that are already linked to a linker, get the index
/// of the linker.
/// @details Not applicable for this particular crosslinker.  A "GNDN" function -- goes nowhere, does nothing.
/// Only here because the base class has a pure virtual of this name.
core::Size
TetrahedralMetal_Helper::get_linker_index_asymmetric(
	core::pose::Pose const &/*pose*/,
	utility::vector1< core::Size > const &/*res_indices*/
) const {
	//GNDN.
	return 0;
}

/// @brief Given indices of four residues that are already linked to pieces of a linker, get
/// of the indices of the symmetric pieces of the linker.
/// @details Not applicable for this particular crosslinker.  A "GNDN" function -- goes nowhere, does nothing.
/// Only here because the base class has a pure virtual of this name.
void
TetrahedralMetal_Helper::get_linker_indices_symmetric(
	core::pose::Pose const &/*pose*/,
	utility::vector1< core::Size > const &/*res_indices*/,
	utility::vector1< core::Size > &linker_indices
) const {
	linker_indices.clear();
	//GNDN.  Deliberately goes nowhere, does nothing.
}

/// @brief Given a pose with residues selected to be linked by a linker, determine whether the residues are too far apart.
/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
/// @note Higher values of the filter multiplier make it more permissive.
bool
TetrahedralMetal_Helper::filter_by_sidechain_distance_asymmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	core::Real const &filter_multiplier
) const {
	utility::vector1< core::Size > res_indices;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices );
	check_residue_indices_valid( res_indices, pose );

	core::Real const comparison_val( filter_multiplier*filter_multiplier*MAX_CB_CB_DIST_SQ );

	for ( core::Size i(2); i<=4; ++i ) {
		for ( core::Size j(1); j<i; ++j ) {
			core::Real const distsq( pose.residue(res_indices[i]).xyz("CB").distance_squared( pose.residue(res_indices[j]).xyz("CB") ) );
			if ( distsq > comparison_val ) return true; //Filter fails if distance > 12
		}
	}

	return false;
}

/// @brief Given a pose with residues selected to be linked by a linker, determine whether the residues are too far apart.
/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.  This version is for symmetric poses.
/// @note Higher values of the filter multiplier make it more permissive.
bool
TetrahedralMetal_Helper::filter_by_sidechain_distance_symmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	core::Real const &filter_multiplier
) const {

	check_compatible_symmetry_type();
	runtime_assert_string_msg( core::pose::symmetry::is_symmetric( pose ), "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::filter_by_sidechain_distance_symmetric(): The pose is not symmetric." );
	runtime_assert_string_msg( CrosslinkerMoverHelper::selection_is_symmetric( selection, pose, symm_subunits_expected() ), "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::filter_by_sidechain_distance_symmetric(): The selection is not symmetric." );
	utility::vector1< core::Size > res_indices;
	CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices );
	check_residue_indices_valid( res_indices, pose );

	core::Real const comparison_val( filter_multiplier*filter_multiplier*MAX_CB_CB_DIST_SQ );

	//Just need to check one residue against the other three, thanks to symmetry.
	for ( core::Size i(2); i<=4; ++i ) {
		core::Real const distsq( pose.residue(res_indices[i]).xyz("CB").distance_squared( pose.residue(res_indices[1]).xyz("CB") ) );
		if ( distsq > comparison_val ) return true; //Filter fails if distance > 12
	}

	return false;
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
/// @note Higher values of the filter multiplier make it more permissive.
bool
TetrahedralMetal_Helper::filter_by_constraints_energy_asymmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	core::Real const &filter_multiplier
) const {
	return filter_by_constraints_energy( pose, selection, false, false, 20*filter_multiplier );
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
/// @note Higher values of the filter multiplier make it more permissive.
bool
TetrahedralMetal_Helper::filter_by_constraints_energy_symmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	bool const /*linker_was_added*/,
	core::Real const &filter_multiplier
) const {
	check_compatible_symmetry_type();
	runtime_assert_string_msg( core::pose::symmetry::is_symmetric( pose ), "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::filter_by_constraints_energy_symmetric(): The pose is not symmetric." );
	runtime_assert_string_msg( CrosslinkerMoverHelper::selection_is_symmetric( selection, pose, symm_subunits_expected() ), "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::filter_by_constraints_energy_symmetric(): The selection is not symmetric." );
	return filter_by_constraints_energy( pose, selection, true, false, 20*filter_multiplier );
}


/// @brief Given a metal name, set the type.
void
TetrahedralMetal_Helper::set_metal_type_from_name(
	std::string const & name_in
) {
	TetrahedralMetal_Helper_Metal const metal_type( metal_type_enum_from_string(name_in) );
	runtime_assert_string_msg( metal_type != TMH_unknown_metal, "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::set_metal_type_from_name(): The metal " + name_in + " is unsupported!" );
	metal_type_ = metal_type;
}

/// @brief Get the current metal type, as a string.
std::string
TetrahedralMetal_Helper::metal_type_string() const {
	return metal_type_string_from_enum( metal_type_ );
}

///////////////////////
/// Private Methods ///
///////////////////////

/// @brief Check that the correct number of residues have been selected, that they are within the pose, and that they are allowed residue types.
void
TetrahedralMetal_Helper::check_residue_indices_valid(
	utility::vector1< core::Size > const &indices,
	core::pose::Pose const &pose
) const {
	static const std::string errmsg( "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::check_residue_indices_valid(): " );
	runtime_assert_string_msg( indices.size() == 4, errmsg + "The number of side-chains needed for tetrahedral metal coordination is FOUR." );

	core::Size const nres( pose.total_residue() );

	for ( core::Size i(1), imax(indices.size()); i<=imax; ++i ) {
		runtime_assert_string_msg( indices[i] > 0 && indices[i] <= nres, errmsg + "All residue indices must be greater than zero and less than or equal to the length of the pose." );
		runtime_assert_string_msg( is_allowed_type( pose.residue_type( indices[i] ) ), errmsg + "A selected residue is not of an allowed type." );
		if ( i>1 ) {
			for ( core::Size j(1); j<i; ++j ) {
				runtime_assert_string_msg( indices[i] != indices[j], errmsg + "There are duplicates in the list of selected residues." );
			}
		}
	}
}

/// @brief Given a residue type, check whether it's an allowed residue type for tetrahedrally coordinating metals.
/// @details Returns "true" for pass (allowed type) and "false" for failure (prohibited type).  Currently, allowed types are L- and D-histidine,
/// L- or D-aspartate, L- or D-glutamate, L- or D-cysteine, L- or D-homocysteine, and the beta-3-amino acid equivalents.
bool
TetrahedralMetal_Helper::is_allowed_type(
	core::chemical::ResidueType const &type
) const {
	core::chemical::AA const type_aa( type.aa() );
	if ( type_aa == core::chemical::aa_his || type_aa == core::chemical::aa_dhi || type_aa == core::chemical::aa_b3h ||
			type_aa == core::chemical::aa_asp || type_aa == core::chemical::aa_das || type_aa == core::chemical::aa_b3d ||
			type_aa == core::chemical::aa_glu || type_aa == core::chemical::aa_dgu || type_aa == core::chemical::aa_b3e ||
			( (type_aa == core::chemical::aa_cys || type_aa == core::chemical::aa_dcs ||  type_aa == core::chemical::aa_b3c ||
			!type.base_name().compare("C26") || !type.base_name().compare("DC26") /*homocysteine*/ )
			&& !type.is_disulfide_bonded()
			)
			) {
		return true;
	}
	return false;
}

/// @brief Given a pose, a list of residues, and indices i and j in that list, add angle constraints between the two residues specified.
void
TetrahedralMetal_Helper::add_angle_constraints(
	core::pose::Pose &pose,
	utility::vector1< core::Size > const &res_indices,
	core::Size const i,
	core::Size const j
) const {

	core::chemical::ResidueType const &rtype1( pose.residue_type(res_indices[i]) );
	core::chemical::ResidueType const &rtype2( pose.residue_type(res_indices[j]) );

	bool const his_i(
		rtype1.aa() == core::chemical::aa_his ||
		rtype1.aa() == core::chemical::aa_dhi
	);
	bool const his_j(
		rtype2.aa() == core::chemical::aa_his ||
		rtype2.aa() == core::chemical::aa_dhi
	);

	std::stringstream cststring;

	static std::string const circularharmonic_string(" CIRCULARHARMONIC 1.9106332359 0.0349066\n"); //Static const data are threadsafe under the C++11 standard.

	if ( his_i || his_j ) cststring << "AmbiguousConstraint\n";
	std::string const parent_i( rtype1.atom_name( rtype1.icoor( rtype1.atom_index("VM1") ).stub_atom1().atomno() ) );
	std::string const parent_j( rtype2.atom_name( rtype2.icoor( rtype1.atom_index("VM1") ).stub_atom1().atomno() ) );
	cststring << "Angle " << parent_i << " " << res_indices[i] << " VM1 " << res_indices[i] << " " << parent_j << " " << res_indices[j] << circularharmonic_string;
	if ( his_i ) {
		std::string const alt_parent_i( (!rtype1.base_name().compare("HIS") || !rtype1.base_name().compare("DHIS") ) ? "NE2" : "ND1" ); //If this is "HIS" (as opposed to "HIS_D"), the parent of the virtual metal is ND1.  So the alternative is NE2.  Else the alternative is ND1.
		cststring << "Angle " << alt_parent_i << " " << res_indices[i] << " VM1 " << res_indices[i] << " " << parent_j << " " << res_indices[j] << circularharmonic_string;
	}
	if ( his_j ) {
		std::string const alt_parent_j( (!rtype2.base_name().compare("HIS") || !rtype2.base_name().compare("DHIS") ) ? "NE2" : "ND1" ); //If this is "HIS" (as opposed to "HIS_D"), the parent of the virtual metal is ND1.  So the alternative is NE2.  Else the alternative is ND1.
		cststring << "Angle " << parent_i << " " << res_indices[i] << " VM1 " << res_indices[i] << " " << alt_parent_j << " " << res_indices[j] << circularharmonic_string;
	}
	if ( his_i && his_j ) {
		std::string const alt_parent_i( (!rtype1.base_name().compare("HIS") || !rtype1.base_name().compare("DHIS") ) ? "NE2" : "ND1" ); //If this is "HIS" (as opposed to "HIS_D"), the parent of the virtual metal is ND1.  So the alternative is NE2.  Else the alternative is ND1.
		std::string const alt_parent_j( (!rtype2.base_name().compare("HIS") || !rtype2.base_name().compare("DHIS") ) ? "NE2" : "ND1" ); //If this is "HIS" (as opposed to "HIS_D"), the parent of the virtual metal is ND1.  So the alternative is NE2.  Else the alternative is ND1.
		cststring << "Angle " << alt_parent_i << " " << res_indices[i] << " VM1 " << res_indices[i] << " " << alt_parent_j << " " << res_indices[j] << circularharmonic_string;
	}

	if ( his_i || his_j ) cststring << "END\n";

	pose.add_constraints( core::scoring::constraints::ConstraintIO::get_instance()->read_constraints_new( cststring, core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet ), pose, false )->get_all_constraints() );
}

/// @brief Given a metal type and a metal-liganding atom type, return the ideal bond length.
/// @details This method must be updated if either enum is expanded.
core::Real const &
TetrahedralMetal_Helper::ideal_bond_length(
	TetrahedralMetal_Helper_Metal const metal_type,
	TetrahedralMetal_Helper_MetalLigand const ligand_type
) const {
	//This map must cover all combinatorial pairs.  Update it if you add to either enum:
	static const std::map< std::pair< TetrahedralMetal_Helper_Metal, TetrahedralMetal_Helper_MetalLigand >, core::Real> ideal_bond_lengths = {
		{ std::make_pair( TMH_Zn, TMHLigand_Nd_histidine ), 2.09 }, //From Tamames et al. (2007). Proteins 69(3): 466-75. DOI: 10.1002/prot.21536
		{ std::make_pair( TMH_Zn, TMHLigand_Ne_histidine ), 2.12 }, //From Tamames et al. (2007). Proteins 69(3): 466-75. DOI: 10.1002/prot.21536
		{ std::make_pair( TMH_Zn, TMHLigand_O_carboxyl ), 2.07 }, //From Tamames et al. (2007). Proteins 69(3): 466-75. DOI: 10.1002/prot.21536
		{ std::make_pair( TMH_Zn, TMHLigand_S_cysteine ), 2.32 } //From Tamames et al. (2007). Proteins 69(3): 466-75. DOI: 10.1002/prot.21536
		/*{ std::make_pair( TMH_Cu, TMHLigand_Nd_histidine ), xxx },
		{ std::make_pair( TMH_Cu, TMHLigand_Ne_histidine ), xxx },
		{ std::make_pair( TMH_Cu, TMHLigand_O_carboxyl ), xxx },
		{ std::make_pair( TMH_Cu, TMHLigand_S_cysteine ), xxx }*/
		}; //This initialization is threadsafe for const global data in C++11.

	std::map< std::pair< TetrahedralMetal_Helper_Metal, TetrahedralMetal_Helper_MetalLigand >, core::Real>::const_iterator data( ideal_bond_lengths.find( std::make_pair( metal_type, ligand_type ) ) );
	runtime_assert_string_msg( data != ideal_bond_lengths.end(), "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::ideal_bond_length(): The metal-metal ligand pair has no defined ideal bond length.  This should not happen.  Consult the protocol developer." );

	return data->second;
}

/// @brief Given a pose and a residue with the VIRTUAL_METAL_CONJUGATION variant type already added,
/// set the metal-metal liganding atom bond length appropriately for the metal in question.
/// @details Calls ideal_bond_length().
void
TetrahedralMetal_Helper::set_metal_bond_length(
	core::pose::Pose &pose,
	core::Size const res_index
) const {
	TetrahedralMetal_Helper_MetalLigand const ligand_type( liganding_atom_from_restype( pose.residue_type(res_index) ) );
	runtime_assert_string_msg( ligand_type != TMHLigand_unknown_ligand,
		"Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::set_metal_bond_length(): an unknown residue type was passed to this function." );
	core::Real const & bondlength( ideal_bond_length( metal_type(), ligand_type ) );
	core::Size const metalindex( pose.residue_type(res_index).atom_index("VM1") );
	core::Size const parentindex( pose.residue_type(res_index).icoor(metalindex).stub_atom1().atomno() );
	numeric::xyzVector< core::Real > const orig_metalposition( pose.residue(res_index).xyz(metalindex) );
	numeric::xyzVector< core::Real > const orig_parentposition( pose.residue(res_index).xyz(parentindex) );
	numeric::xyzVector< core::Real > const shifted_position( (orig_metalposition - orig_parentposition).normalize() * bondlength + orig_parentposition );
	pose.set_xyz( core::id::AtomID( metalindex, res_index ), shifted_position );
}

/// @brief Given a ResidueType with the VIRTUAL_METAL_CONJUGATION variant type already added, get the metal-liganding
/// atom enum.
TetrahedralMetal_Helper_MetalLigand
TetrahedralMetal_Helper::liganding_atom_from_restype(
	core::chemical::ResidueType const &restype
) const {
	if ( restype.is_sidechain_thiol() ) {
		return TMHLigand_S_cysteine;
	}
	if ( restype.aa() == core::chemical::aa_asp || restype.aa() == core::chemical::aa_glu ||
			restype.aa() == core::chemical::aa_das || restype.aa() == core::chemical::aa_dgu ||
			restype.aa() == core::chemical::aa_b3d || restype.aa() == core::chemical::aa_b3e
			) {
		return TMHLigand_O_carboxyl;
	}
	if ( restype.aa() == core::chemical::aa_his || restype.aa() == core::chemical::aa_dhi || restype.aa() == core::chemical::aa_b3h ) {
		std::string const & basename( restype.base_name() );
		if ( !basename.compare( "HIS" ) || !basename.compare( "DHIS" ) || !basename.compare( "B3H" ) ) {
			return TMHLigand_Nd_histidine;
		} else {
			return TMHLigand_Ne_histidine;
		}
	}
	return TMHLigand_unknown_ligand;
}

/// @brief Given a metal type enum, return the corresponding string.
std::string
TetrahedralMetal_Helper::metal_type_string_from_enum(
	TetrahedralMetal_Helper_Metal const metal_type
) const {
	switch( metal_type ) {
	case TMH_Zn :
		return "Zn";
		//case TMH_Cu:
		// return "Cu";
	case TMH_unknown_metal :
		utility_exit_with_message( "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::metal_type_string_from_enum(): The metal is unknown!" );
	}
	return ""; //To keep compiler happy.
}

/// @brief Given a metal type string, return the corresponding enum.
/// @details Returns TMH_unknown_metal if the string is unrecognized.
TetrahedralMetal_Helper_Metal
TetrahedralMetal_Helper::metal_type_enum_from_string(
	std::string const & metal_type_string
) const {
	for ( core::Size i(1); i < static_cast<core::Size>( TMH_end_of_list ); ++i ) {
		if ( !metal_type_string.compare( metal_type_string_from_enum( static_cast< TetrahedralMetal_Helper_Metal >(i) ) ) ) {
			return static_cast< TetrahedralMetal_Helper_Metal >(i);
		}
	}
	return TMH_unknown_metal;
}

/// @brief Check that the symmetry type is one of a few compatible types.
/// @details Allowed types are C2, S4, and D2.
void
TetrahedralMetal_Helper::check_compatible_symmetry_type() const {
	bool passed(false);
	if ( symm_type() == 'S' ) {
		if ( symm_count() == 4 ) {
			passed = true;
		}
	} else if ( symm_type() == 'C' || symm_type() == 'D' ) {
		if ( symm_count() == 2 ) {
			passed = true;
		}
	}
	runtime_assert_string_msg(passed, "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::check_compatible_symmetry_type(): Tetrahedral metal crosslinks are only compatible with C2, S4, or D2 symmetry.");
}



} //crosslinker
} //protocols
} //cyclic_peptide
