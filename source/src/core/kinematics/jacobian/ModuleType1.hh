// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/kinematics/jacobian/ModuleType1.hh
/// @brief Class that contains the linearized Cartesian motion/force space of a 1,2,3-amino acid module
/// @author teunhoevenaars (teunhoevenaars@gmail.com)


#ifndef INCLUDED_core_kinematics_jacobian_ModuleType1_hh
#define INCLUDED_core_kinematics_jacobian_ModuleType1_hh

#include <core/kinematics/jacobian/ModuleType1.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

// Protocol headers
#include <core/id/AtomID.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/MathMatrix.hh>

namespace core {
namespace kinematics {
namespace jacobian {

/// @brief The ModuleType1 class covers the Jacobian analysis at the lower-level for type-1 systems, which have internal
/// DoFs that can be organized as three sets of two torsion angles with intersecting axes.
/// @author teunhoevenaars (teunhoevenaars@gmail.com)

class ModuleType1 : public utility::VirtualBase {

public: // TYPEDEFS

	/// @brief typedef of a screw as a pair of xyzVectors
	typedef std::pair < numeric::xyzVector<core::Real>, numeric::xyzVector<core::Real> > Screw;

	/// @brief Struct to pass around all variations of Jacobian matrices
	/// @detail fw_dofs is the Jacobian matrix that maps differential torsion angles [rad/s] or [delta rad] on Cartesian reference frame
	/// fw_cons is the 6x6 Jacobian matrix that maps virtual (helper, constrained) torsion angles [rad/s] or [delta rad] on Cartesian reference frame
	/// fw_all is the 6x6 Jacobian matrix that maps all six torsion angles (real and virtual) [rad/s] or [delta rad] on Cartesian reference frame
	/// inv_dofs is the 6x6 Jacobian matrix that maps twist from Cartesian ref. frame onto torsion angles
	/// inv_cons is the 6x6 Jacobian matrix that maps twist from Cartesian ref. frame onto virtual torsion angles that represent the modules constraints
	/// inv_all is the 6x6 Jacobian matrix that maps twist from Cartesian ref. frame onto all torsion angles (real and virtual)
	struct jacobian_struct {
		numeric::MathMatrix<core::Real> fw_dofs{6,6, core::Real(0)};
		numeric::MathMatrix<core::Real> fw_cons{6,6, core::Real(0)};
		numeric::MathMatrix<core::Real> fw_all{6,6, core::Real(0)};
		numeric::MathMatrix<core::Real> inv_dofs{6,6, core::Real(0)};
		numeric::MathMatrix<core::Real> inv_cons{6,6, core::Real(0)};
		numeric::MathMatrix<core::Real> inv_all{6,6, core::Real(0)};
	};

	/// @brief struct with the different linear algebra objects that are needed to express the twists and wrenches to
	/// calculate the instantaneous Jacobian matrices of the module
	/// @details Hmatrices vector with the homogeneous matrices that describe the homogeneous transformations between the
	/// reference atom (ref_atom_) and the first atom in key_atoms_IDs_, followed by the homogeneous transformations between
	/// the subsequent atoms in key_atoms_IDs_.
	/// s_vectors is a vector1 with screw axis vectors
	/// r_vectorss is a vector1 with moment arm vectors for screw definitions
	struct screw_construct_struct {
		utility::vector1< numeric::HomogeneousTransform<core::Real> > Hmatrices{6};
		utility::vector1< numeric::xyzVector< core::Real > > s_vectors{6};
		utility::vector1< numeric::xyzVector< core::Real >  > r_vectors{6};
	};

public:
	/// @brief No default constructor because modules are not used on their own, but always as part of a SeriesJacobian
	ModuleType1() = delete;

	/// @brief Constructor for module with 1, 2 or 3 sets of two AtomIDs whose X-axes represent intersecting torsion axes.
	ModuleType1(core::Size const & dofs_module, utility::vector1< core::Size > const & res_numbers, core::id::AtomID const & ref_atom);

	/// @brief Copy constructor.
	ModuleType1(ModuleType1 const & src);

	/// @brief Destructor.
	~ModuleType1() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	ModuleType1OP clone() const;

	/// @brief Returns the residues that make up the Jacobian module
	utility::vector1< core::Size >
	get_residues(){
		return res_numbers_;
	}

	/// @brief Returns the number of DoFs of which this module represents the differential equations
	core::Size
	get_number_dofs_(){
		return number_dofs_;
	}

	/// @brief Returns the AtomID of the atom whose reference frame (expressed by its stub) is used to express all
	/// vectors and matrices of the module
	core::id::AtomID const &
	get_ref_atom_ID() {
		return ref_atom_ID_;
	}

	/// @brief calculate all Jacobian matrices for the instantaneous conformation
	ModuleType1::jacobian_struct
	update_jacobian_matrices(core::conformation::Conformation const & conformation);

	/// @brief apply vector of delta angles to the torsion angles that are associated to this instance
	/// @param[in] vector dq must have dimension 6 and all delta angles must be in radians
	void
	apply_delta_torsions_angles(core::conformation::Conformation & conformation, numeric::MathVector<core::Real> const & dq);

private: // FUNCTIONS
	/// @brief Updates the screw axes and moment arm vectors of the object
	ModuleType1::screw_construct_struct
	update_screw_vectors(core::conformation::Conformation const & conformation);

	/// @brief Returns the two wrenches that are dual to the input twists, which is used in Jacobian analysis
	utility::vector1< Screw >
	get_orthogonal_wrenches_from_twists(Screw const & twist_1,
		Screw const & twist_2,
		numeric::HomogeneousTransform<core::Real> const & H_CA0_CAa,
		numeric::HomogeneousTransform<core::Real> const & H_CAa_Ca,
		numeric::HomogeneousTransform<core::Real> const & H_Ca_CAb,
		numeric::HomogeneousTransform<core::Real> const & H_CAb_Cb);

private: // VARIABLES
	/// @brief number of dofs of the module
	core::Size number_dofs_;

	/// @brief residue numbers of the module
	utility::vector1<core::Size> res_numbers_;

	/// @brief the atom which acts as the reference frame for all vectors
	core::id::AtomID ref_atom_ID_;

private: // PARAMETERS
	/// @brief bound to check for singularity in Jacobian denominator
	static constexpr core::Real denom_singular_bound_ = 1e-2; //
	/// @brief bound to check for screw orthogonality
	static constexpr core::Real denom_orthogonal_bound_ = 1e-6; //
	/// @brief bound to check for singularity in atom position
	static constexpr core::Real pos_singular_bound_ = 1e-10; // [Ang], should be significantly smaller than denom_orthogonal_bound_
};

} //jacobian
} //kinematics
} //core

#endif //INCLUDED_core_kinematics_jacobian_ModuleType1_hh
