// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/SquarePyramidalMetal_Helper.cc
/// @brief A helper class for setting up square pyramidally-coordinated metals
/// like Ni.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/SquarePyramidalMetal_Helper.hh>

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

static basic::Tracer TR( "protocols.cyclic_peptide.crosslinker.SquarePyramidalMetal_Helper" );

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

#define MAX_CB_CB_DIST_SQ 144.0

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SquarePyramidalMetal_Helper::SquarePyramidalMetal_Helper(
	std::string const &metal_name_in
) :
	SquarePlanarMetal_Helper( metal_name_in )
	//TODO initialize data here
{
	set_pyramidal(true);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
SquarePyramidalMetal_Helper::SquarePyramidalMetal_Helper( SquarePyramidalMetal_Helper const &src ) :
	SquarePlanarMetal_Helper( src )
	//TODO copy data here
{
	runtime_assert(is_pyramidal());
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
SquarePyramidalMetal_Helper::~SquarePyramidalMetal_Helper(){}

///////////////////////
/// Public Methods  ///
///////////////////////


///////////////////////
/// Private Methods ///
///////////////////////

/// @brief Check that the correct number of residues have been selected, that they are within the pose, and that they are allowed residue types.
void
SquarePyramidalMetal_Helper::check_residue_indices_valid(
	utility::vector1< core::Size > const &indices,
	core::pose::Pose const &pose
) const {
	static const std::string errmsg( "Error in protocols::cyclic_peptide::crosslinker::SquarePyramidalMetal_Helper::check_residue_indices_valid(): " );
	runtime_assert_string_msg( indices.size() == 5, errmsg + "The number of side-chains needed for square pyramidal metal coordination is FIVE." );

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

/// @brief Add improper dihedral constraints to keep the metal in the plane of the liganding atoms.
void
SquarePyramidalMetal_Helper::add_dihedral_constraints(
	core::pose::Pose &pose,
	utility::vector1< core::Size > const &res_indices
) const {
	runtime_assert_string_msg( res_indices.size() == 5, "Error in SquarePyramidalMetal_Helper::add_dihedral_constraints(): FIVE indices must be specified for square pyramidal metal coordination." );
	std::stringstream cststring;
	cststring << "AmbiguousConstraint\n";
	for ( core::Size ii(1); ii<=5; ++ii ) {
		core::Size counter(0);
		utility::vector1< core::Size > res_indices_abridged(4);
		for ( core::Size jj(1); jj<=5; ++jj ) {
			if ( jj == ii ) continue;
			++counter;
			res_indices_abridged[counter] = res_indices[jj];
		}
		add_dihedral_constraint_to_stream(cststring, pose, res_indices_abridged);
	}
	cststring << "END\n";
	//std::cout << cststring.str(); //DELETE ME -- FOR DEBUGGING ONLY
	pose.add_constraints( core::scoring::constraints::ConstraintIO::get_instance()->read_constraints_new( cststring, core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet ), pose, false )->get_all_constraints() );
}

/// @brief Check that the symmetry type is one of a few compatible types.
/// @details Square pyramidal metals are inherently asymmetric, so this always throws.
void
SquarePyramidalMetal_Helper::check_compatible_symmetry_type() const {
	runtime_assert_string_msg(false, "Error in protocols::cyclic_peptide::crosslinker::SquarePyramidalMetal_Helper::check_compatible_symmetry_type(): Square pyramidal metal crosslinks are incompatible with any symmetry.");
}



} //crosslinker
} //protocols
} //cyclic_peptide
