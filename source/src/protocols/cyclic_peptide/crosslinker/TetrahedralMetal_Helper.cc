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
	Metal_HelperBase( metal_name_in )
	//TODO initialize data here
{
	set_metal_type_from_name(metal_name_in);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
TetrahedralMetal_Helper::TetrahedralMetal_Helper( TetrahedralMetal_Helper const &src ) :
	Metal_HelperBase( src )
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
	Metal_HelperBase_Metal const metal_type,
	Metal_HelperBase_MetalLigand const ligand_type
) const {
	//This map must cover all combinatorial pairs.  Update it if you add to either enum:
	static const std::map< std::pair< Metal_HelperBase_Metal, Metal_HelperBase_MetalLigand >, core::Real> ideal_bond_lengths = {
		{ std::make_pair( MH_Zn, MHLigand_Nd_histidine ), 2.09 }, //From Tamames et al. (2007). Proteins 69(3): 466-75. DOI: 10.1002/prot.21536
		{ std::make_pair( MH_Zn, MHLigand_Ne_histidine ), 2.12 }, //From Tamames et al. (2007). Proteins 69(3): 466-75. DOI: 10.1002/prot.21536
		{ std::make_pair( MH_Zn, MHLigand_O_carboxyl ), 2.07 }, //From Tamames et al. (2007). Proteins 69(3): 466-75. DOI: 10.1002/prot.21536
		{ std::make_pair( MH_Zn, MHLigand_S_cysteine ), 2.32 } //From Tamames et al. (2007). Proteins 69(3): 466-75. DOI: 10.1002/prot.21536
		/*{ std::make_pair( MH_Cu, MHLigand_Nd_histidine ), xxx },
		{ std::make_pair( MH_Cu, MHLigand_Ne_histidine ), xxx },
		{ std::make_pair( MH_Cu, MHLigand_O_carboxyl ), xxx },
		{ std::make_pair( MH_Cu, MHLigand_S_cysteine ), xxx }*/
		}; //This initialization is threadsafe for const global data in C++11.

	std::map< std::pair< Metal_HelperBase_Metal, Metal_HelperBase_MetalLigand >, core::Real>::const_iterator data( ideal_bond_lengths.find( std::make_pair( metal_type, ligand_type ) ) );
	runtime_assert_string_msg( data != ideal_bond_lengths.end(), "Error in protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper::ideal_bond_length(): Tetrahedral coordination information is unavailable for the given ligand type when it binds to " + metal_type_string_from_enum(metal_type) + "." );

	return data->second;
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
