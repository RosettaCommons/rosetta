// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/TrigonalPlanarMetal_Helper.cc
/// @brief A helper class for setting up trigonal planarly-coordinated metals
/// like Zn or Cu.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/TrigonalPlanarMetal_Helper.hh>

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

static basic::Tracer TR( "protocols.cyclic_peptide.crosslinker.TrigonalPlanarMetal_Helper" );

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

#define MAX_CB_CB_DIST_SQ 144.0

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
TrigonalPlanarMetal_Helper::TrigonalPlanarMetal_Helper(
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
TrigonalPlanarMetal_Helper::TrigonalPlanarMetal_Helper( TrigonalPlanarMetal_Helper const &src ) :
	Metal_HelperBase( src )
	//TODO copy data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
TrigonalPlanarMetal_Helper::~TrigonalPlanarMetal_Helper(){}

///////////////////////
/// Public Methods  ///
///////////////////////


///////////////////////
/// Private Methods ///
///////////////////////

/// @brief Check that the correct number of residues have been selected, that they are within the pose, and that they are allowed residue types.
void
TrigonalPlanarMetal_Helper::check_residue_indices_valid(
	utility::vector1< core::Size > const &indices,
	core::pose::Pose const &pose
) const {
	static const std::string errmsg( "Error in protocols::cyclic_peptide::crosslinker::TrigonalPlanarMetal_Helper::check_residue_indices_valid(): " );
	runtime_assert_string_msg( indices.size() == 3, errmsg + "The number of side-chains needed for trigonal planar metal coordination is THREE." );

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

/// @brief Given a residue type, check whether it's an allowed residue type for trigonal planarly coordinating metals.
/// @details Returns "true" for pass (allowed type) and "false" for failure (prohibited type).  Currently, allowed types are L- and D-histidine,
/// L- or D-aspartate, L- or D-glutamate, and the beta-3-amino acid equivalents.
bool
TrigonalPlanarMetal_Helper::is_allowed_type(
	core::chemical::ResidueType const &type
) const {
	core::chemical::AA const type_aa( type.aa() );
	if ( type_aa == core::chemical::aa_his || type_aa == core::chemical::aa_dhi || type_aa == core::chemical::aa_b3h ||
			type_aa == core::chemical::aa_asp || type_aa == core::chemical::aa_das || type_aa == core::chemical::aa_b3d ||
			type_aa == core::chemical::aa_glu || type_aa == core::chemical::aa_dgu || type_aa == core::chemical::aa_b3e
			//|| ( (type_aa == core::chemical::aa_cys || type_aa == core::chemical::aa_dcs ||  type_aa == core::chemical::aa_b3c ||
			//!type.base_name().compare("C26") || !type.base_name().compare("DC26") /*homocysteine*/ )
			//&& !type.is_disulfide_bonded()
			//)
			) {
		return true;
	}
	return false;
}

/// @brief Given a pose, a list of residues, and indices i and j in that list, add angle constraints between the two residues specified.
void
TrigonalPlanarMetal_Helper::add_angle_constraints(
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

	static std::string const circularharmonic_string(" CIRCULARHARMONIC 2.09439510267 0.017453292522\n"); //Static const data are threadsafe under the C++11 standard.

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

/// @brief Add improper dihedral constraints to keep the metal in the plane of the liganding atoms.
void
TrigonalPlanarMetal_Helper::add_dihedral_constraints(
	core::pose::Pose &pose,
	utility::vector1< core::Size > const &res_indices
) const {
	static std::string const errmsg( "Error in TrigonalPlanarMetal_Helper::add_dihedral_constraints(): " );
	runtime_assert_string_msg( res_indices.size() == 3, errmsg + "The number of selected residues must be 3.");
	static std::string const circularharmonic_string(" CIRCULARHARMONIC 0.0 0.052359877566666664\n"); //Static const data are threadsafe under the C++11 standard.


	utility::fixedsizearray1< core::chemical::ResidueTypeCOP, 3 > const rtypes {  pose.residue_type_ptr(res_indices[1]),  pose.residue_type_ptr(res_indices[2]),  pose.residue_type_ptr(res_indices[3]) };
	utility::fixedsizearray1< bool, 3 > is_his;
	for ( core::Size i(1); i<=3; ++i ) {
		is_his[i] = ( rtypes[i]->aa() == core::chemical::aa_his || rtypes[i]->aa() == core::chemical::aa_dhi );
	}
	bool const any_his( is_his[1] || is_his[2] || is_his[3] );
	static utility::fixedsizearray1< utility::fixedsizearray1 < core::Size, 3 >, 3 > const res_permutations {
		utility::fixedsizearray1< core::Size, 3 > { 1, 2, 3 },
		utility::fixedsizearray1< core::Size, 3 > { 2, 3, 1 },
		utility::fixedsizearray1< core::Size, 3 > { 3, 1, 2 }
		};
	for ( core::Size permutation(1); permutation <= 3; ++permutation ) {
		std::stringstream cststring;
		if ( any_his ) {
			cststring << "AmbiguousConstraint\n";
		}
		for ( core::Size i(1), imax( is_his[res_permutations[permutation][1]] ? 2 : 1 ); i<=imax; ++i ) {
			core::chemical::ResidueType const rtype_i( *( rtypes[res_permutations[permutation][1]] ) );
			std::string const at1( is_his[res_permutations[permutation][1]] ? ( i == 1 ? "NE2" : "ND1" ) : rtype_i.atom_name( rtype_i.icoor( rtype_i.atom_index("VM1") ).stub_atom1().atomno() ) );
			for ( core::Size j(1), jmax( is_his[res_permutations[permutation][2]] ? 2 : 1 ); j<=jmax; ++j ) {
				core::chemical::ResidueType const rtype_j( *( rtypes[res_permutations[permutation][2]] ) );
				std::string const at2( is_his[res_permutations[permutation][2]] ? ( j == 1 ? "NE2" : "ND1" ) : rtype_j.atom_name( rtype_j.icoor( rtype_j.atom_index("VM1") ).stub_atom1().atomno() ) );
				for ( core::Size k(1), kmax( is_his[res_permutations[permutation][3]] ? 2 : 1 ); k<=kmax; ++k ) {
					core::chemical::ResidueType const rtype_k( *( rtypes[res_permutations[permutation][3]] ) );
					std::string const at3( is_his[res_permutations[permutation][3]] ? ( k == 1 ? "NE2" : "ND1" ) : rtype_k.atom_name( rtype_k.icoor( rtype_k.atom_index("VM1") ).stub_atom1().atomno() ) );
					cststring << "Dihedral " << at1 << " " << res_indices[res_permutations[permutation][1]] << " " << at2 << " " << res_indices[res_permutations[permutation][2]] << " "<< at3 << " " << res_indices[res_permutations[permutation][3]] << " VM1 " << res_indices[res_permutations[permutation][1]] << " " << circularharmonic_string;
				}
			}
		}
		if ( any_his ) {
			cststring << "END\n";
		}
		//std::cout << cststring.str(); //DELETE ME -- FOR DEBUGGING ONLY
		pose.add_constraints( core::scoring::constraints::ConstraintIO::get_instance()->read_constraints_new( cststring, core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet ), pose, false )->get_all_constraints() );
	}
}

/// @brief Given a metal type and a metal-liganding atom type, return the ideal bond length.
/// @details This method must be updated if either enum is expanded.
core::Real const &
TrigonalPlanarMetal_Helper::ideal_bond_length(
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
	runtime_assert_string_msg( data != ideal_bond_lengths.end(), "Error in protocols::cyclic_peptide::crosslinker::TrigonalPlanarMetal_Helper::ideal_bond_length(): Trigonal planar coordination information is unavailable for the given ligand type when it binds to " + metal_type_string_from_enum(metal_type) + "." );

	return data->second;
}

/// @brief Check that the symmetry type is one of a few compatible types.
/// @details Allowed types are C3.
void
TrigonalPlanarMetal_Helper::check_compatible_symmetry_type() const {
	bool passed(false);
	if ( symm_type() == 'C' ) {
		if ( symm_count() == 3 ) {
			passed = true;
		}
	}
	runtime_assert_string_msg(passed, "Error in protocols::cyclic_peptide::crosslinker::TrigonalPlanarMetal_Helper::check_compatible_symmetry_type(): Trigonal planar metal crosslinks are only compatible with C3 symmetry.");
}



} //crosslinker
} //protocols
} //cyclic_peptide
