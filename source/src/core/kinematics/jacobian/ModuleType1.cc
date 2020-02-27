// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/kinematics/jacobian/ModuleType1.cc
/// @brief Class to perform Jacobian analyses on a series of atoms with mobility in three Universal joints (collocated torsions)
/// @author teunhoevenaars (teunhoevenaars@gmail.com)

/// @details
/// There are different methods to obtain a Jacobian, but here the method described in Huang et al. (2011) "Generalized
/// Jacobian analysis of lower mobility manipulators" is used. Idea towards future is to set up a general Module class,
/// where each ModuleTypeX class is a derived class.

// Project headers:
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/jacobian/ModuleType1.hh>
#include <core/types.hh>
#include <core/kinematics/AtomTree.hh>
#include <numeric/MathVector.hh>
#include <numeric/MathVector_operations.hh>
#include <numeric/MathMatrix_operations.hh>
#include <numeric/xyzVector.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>
#include <utility/vector1.hh>

static basic::Tracer TR( "core.kinematics.jacobian.ModuleType1" );

namespace core {
namespace kinematics {
namespace jacobian {

// allocate space for constexpr variables
constexpr core::Real ModuleType1::denom_singular_bound_;
constexpr core::Real ModuleType1::denom_orthogonal_bound_;
constexpr core::Real ModuleType1::pos_singular_bound_;

/// @brief Constructor for module with 1, 2 or 3 sets of two AtomIDs whose X-axes represent intersecting torsion axes.
/// @details Stores the basic information needed to perform the Jacobian analysis for this module for any conformation
/// @param[in] key_atoms Must have exactly six AtomIDs, whose associated stubs' X-axes identify the torsion axes
ModuleType1::ModuleType1(core::Size const & dofs_module, utility::vector1< core::Size > const & res_numbers, core::id::AtomID const & ref_atom){
	if ( res_numbers.size() != 3 ) {
		utility_exit_with_message("the res_numbers input vectors does not contain exactly three residues");
	}

	number_dofs_ = dofs_module;
	res_numbers_ = res_numbers;
	ref_atom_ID_ = ref_atom;
}

/// @brief Copy constructor.  Keep default unless deep copying is needed (and in that case,
/// consider using DeepCopyOPs.)
ModuleType1::ModuleType1(ModuleType1 const & ) = default;

/// @brief Destructor.
ModuleType1::~ModuleType1(){}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
ModuleType1OP
ModuleType1::clone() const {
	return utility::pointer::make_shared< ModuleType1 >(*this );
}

/// @brief Calculates all Jacobian matrices for the instantaneous conformation
/// @details Based on the definition of the torsion axes (provided by means of AtomIDs), this function calculates the set
/// of six twists and wrenches for this module. A twist is a dual Cartesian vector: one xyzVector expresses the
/// differential Cartesian rotational vector, and another expresses the differential Cartesian linear vector.
/// See Huang et al. (2011) Generalized Jacobian analysis of lower mobility manipulators for details on the definition
/// of twists and wrenches and their relations.
ModuleType1::jacobian_struct
ModuleType1::update_jacobian_matrices(core::conformation::Conformation const & conformation) {
	// update the vectors that are used in the Jacobian construction
	screw_construct_struct screw_construct = update_screw_vectors(conformation);

	// initialize empty construction twists and wrenches.
	utility::vector1< Screw > twists(6);
	utility::vector1< Screw > wrenches(6);

	// use the Screw vectors to define the six twists
	for ( core::Size i=1; i <= 6; ++i ) {
		twists[i].first = screw_construct.s_vectors[i];
		twists[i].second = screw_construct.s_vectors[i].cross_product(-screw_construct.r_vectors[i]);
	}

	// initialize empty Homogeneous helper matrices, which express vectors that are required to calculate the wrenches
	// that are dual to the above twists
	numeric::HomogeneousTransform<core::Real> H_CA0_CAa;
	numeric::HomogeneousTransform<core::Real> H_CAa_Ca;
	numeric::HomogeneousTransform<core::Real> H_Ca_CAb;
	numeric::HomogeneousTransform<core::Real> H_CAb_Cb;

	// calculate the various wrenches, which is done in sets of two. Each set of two has specific Homogeneous helper matrices
	for ( core::Size i=1; i <= 3; ++i ) {
		// define the two indices that make up the current set of two twists
		core::Size const index1 = 2*i-1;
		core::Size const index2 = 2*i;
		if ( i == 1 ) { // helper matrices to calculate wrenches that are dual to twists 1 and 2
			H_CA0_CAa = screw_construct.Hmatrices[1] * screw_construct.Hmatrices[2] * screw_construct.Hmatrices[3];
			H_CAa_Ca = screw_construct.Hmatrices[4];
			H_Ca_CAb = screw_construct.Hmatrices[5];
			H_CAb_Cb = screw_construct.Hmatrices[6];
		} else if ( i == 2 ) { // helper matrices to calculate wrenches that are dual to twists 3 and 4
			H_CA0_CAa = screw_construct.Hmatrices[1];
			H_CAa_Ca = screw_construct.Hmatrices[2];
			H_Ca_CAb = screw_construct.Hmatrices[3] * screw_construct.Hmatrices[4] * screw_construct.Hmatrices[5];
			H_CAb_Cb = screw_construct.Hmatrices[6];
		} else if ( i == 3 ) { // helper matrices to calculate wrenches that are dual to twists 5 and 6
			H_CA0_CAa = screw_construct.Hmatrices[1];
			H_CAa_Ca = screw_construct.Hmatrices[2];
			H_Ca_CAb = screw_construct.Hmatrices[3];
			H_CAb_Cb = screw_construct.Hmatrices[4];
		}
		// call the function that calculates the wrenches that are dual to each set of two twists.
		utility::vector1<Screw> const wrench_duo =
			get_orthogonal_wrenches_from_twists(twists[index1], twists[index2],
			H_CA0_CAa, H_CAa_Ca, H_Ca_CAb, H_CAb_Cb);
		// store the calculated wrenches at the same indices as the twists.
		wrenches[index1] = wrench_duo[1];
		wrenches[index2] = wrench_duo[2];
	}
	// Initialize core matrices for Jacobian analysis methods based on Screw Theory. A distinction is made between
	// Ta and Wa, which are the matrices that contain the twists and wrenches of actuation (i.e. DoFs, motions described
	// by phi/psi) and Tc and Tc, which are the matrices that contain the twists and wrenches of constraint (i.e. no
	// motion allowed along the associated axes).
	// N.B. These or 0-indexed!
	numeric::MathMatrix<core::Real> Ta(6,6,core::Real(0));
	numeric::MathMatrix<core::Real> Tc(6,6,core::Real(0));
	numeric::MathMatrix<core::Real> Wa(6,6,core::Real(0));
	numeric::MathMatrix<core::Real> Wc(6,6,core::Real(0));

	// copy all values from the screw object to the appropriate position in the MathMatrices
	for ( core::Size i=0; i < 6; ++i ) { // 6 DoFs in total
		for ( core::Size j=0; j < 3; ++j ) { // each xyzVector has 3 entries
			// Twists and wrenches associated with DoFs are copied into Ta and Wa
			if ( i < number_dofs_ ) {
				// in the matrix containing all twist, each twist populates one column
				Ta(j, i) = twists[i+1].first(j+1);
				Ta(j+3, i) = twists[i+1].second(j+1);
				// in the matrix containing all wrenches, each wrench populates one ROW (technically, wrenches are co-vectors)
				Wa(i, j) = wrenches[i+1].first(j + 1);
				Wa(i, j+3) = wrenches[i+1].second(j+1);
			} else {
				// twists and wrenches associated with constraints are copied into Tc and Wc
				// N.B. this matrix will only contain non-zeroes values if a module has <6 DoFs
				Tc(j,i) = twists[i+1].first(j+1);
				Tc(j+3, i) = twists[i+1].second(j+1);
				// in the matrix containing all wrenches, each wrench populates one ROW (technically, wrenches are co-vectors)
				Wc(i,j) = wrenches[i+1].first(j+1);
				Wc(i,j+3) = wrenches[i+1].second(j+1);
			}
		}
	}

	// Pre-calculate denominator for Jacobian analysis, which is also used to verify twists and wrenches are orthogonal
	// and not in singularity.
	// Wa*Ta should give diagonal matrix with values >0 on elements 0..(N_DoFs-1), Wa*Tc should only give zeros,
	// Wc*Tc should only give zeros, and Wc*Tc should give diagonal matrix with values >0 on elements N_DoFs..5.
	// Above matrices are checked in two groups: one group that is input to subsequent analysis (jac_denominators), and
	// one group that is only used to check.
	numeric::MathMatrix<core::Real> jac_denominators = (Wa * Ta) + (Wc * Tc);
	numeric::MathMatrix<core::Real> const jac_denominators_check = (Wa * Tc) + (Wc * Ta);

	for ( core::Size i=0; i < 6; ++i ) {
		for ( core::Size j = 0; j < 6; ++j ) {
			// the jac_denominators matrix should have finite values on its diagonal axes...
			if ( i == j && std::abs(jac_denominators(i, j)) < denom_singular_bound_ ) { // check for ill-defined matrix
				if ( jac_denominators(i,j) == 0 ) { // true singularity
					jac_denominators(i,j) = 1; // then take default value of 1
				} else { // ill-defined
					// get sign of the value by check for >0 and set denominator to either -0.01 or +0.01
					jac_denominators(i, j) = numeric::sign(jac_denominators(i, j)) * denom_singular_bound_; // ensure that denominators are not singular
					TR.Debug << "jacobian of module including residue " << res_numbers_[1]
						<< " is in a singular position and Jacobian has been capped to avoid large values"
						<< std::endl;
				}
				// ... but approximately zero values for the off-diagonal terms
			} else if ( i != j ) {
				runtime_assert_string_msg(std::abs(jac_denominators(i, j)) < denom_orthogonal_bound_, "Jacobian analysis failed"); // verify that all non-diagonal values are roughly zero
			}
			// the values in the jac_denominators_check matrix should all be approximately zero
			runtime_assert_string_msg(std::abs(jac_denominators_check(i, j)) < denom_orthogonal_bound_, "Jacobian analysis failed"); // verify orthogonality of twists and wrenches
		}
	}

	// Construct the various Jacobian matrices.
	jacobian_struct jacs;
	// Forward Jacobian matrices are direct copies of the twist matrices
	jacs.fw_dofs = Ta;
	jacs.fw_cons = Tc;
	jacs.fw_all = Ta + Tc;

	// Inverse Jacobian matrices need some processing, and are initialized empty
	numeric::MathMatrix<core::Real> inv_jac_dofs(6,6,core::Real(0) );
	numeric::MathMatrix<core::Real> inv_jac_constraints(6,6,core::Real(0) );
	// subsequently each row is computed according to Huang et al. 2011
	for ( core::Size i=0; i < 6; ++i ) { // 6 DoFs in total
		if ( i < number_dofs_ ) { // for allowed DoFs
			inv_jac_dofs.replace_row(i, numeric::operator/( Wa.get_row(i) , jac_denominators(i,i) ) );
		} else { // for constraints
			inv_jac_constraints.replace_row(i, numeric::operator/( Wc.get_row(i) , jac_denominators(i,i) ) );
		}
	}

	// inv_jac_dofs_ has the first number_dofs_ rows populated, while inv_jac_constraints_ has the remaining
	// (6-number_dofs_) populated. inv_jac_all_  is the combination of both Jacobians
	jacs.inv_dofs = inv_jac_dofs;
	jacs.inv_cons = inv_jac_constraints;
	jacs.inv_all = inv_jac_dofs + inv_jac_constraints;

	return jacs;
}

/// @brief Apply vector of delta angles to the torsion angles that are associated to this instance
/// @param[in] vector dq must have dimension 6 and all delta angles must be in radians
void
ModuleType1::apply_delta_torsions_angles(core::conformation::Conformation & conformation, numeric::MathVector<core::Real> const & dq){
	// Assert that the provided vector has size 6
	runtime_assert_string_msg(dq.size() == 6, "input vector 'dq' is not of size six, check input");

	// loop over dofs of the module
	for ( core::Size i=0; i < number_dofs_; ++i ) {
		// get torsion ID ('+ 1' because first two dofs are related to first residue)
		core::id::TorsionID torsion_id(res_numbers_[i / 2 + 1], core::id::BB, i % 2 + 1 );
		// apply torsion
		conformation.set_torsion( torsion_id, conformation.torsion(torsion_id) + dq(i) * numeric::constants::d::radians_to_degrees ); // in [deg]
	}
}

/// @brief Updates the screw axes and moment arm vectors of the object
/// @detail For the Jacobian analysis axis vectors and arm vectors are required
/// Handles to the six torsion axes are provided by means of AtomIDs. AtomID uniquely points to a Stub, whose X-axes
/// correspond to torsion axes. Thus, to identify the axes of the phi and psi torsions the stubs associated to respectively
/// the mainchain CA- and C-atoms needs to be provided.
ModuleType1::screw_construct_struct
ModuleType1::update_screw_vectors(core::conformation::Conformation const & conformation) {
	// Create temporary variables:
	// create empty vector of homogeneous matrices
	utility::vector1< numeric::HomogeneousTransform<core::Real> > Hmatrices;
	// the moment arm vectors are all initialized as zero vectors
	utility::vector1< numeric::xyzVector< core::Real > > r_vectors(6, numeric::xyzVector< core::Real >(0, 0, 0));
	// the screw axis vectors are initialized as unit vectors along the x-axis
	utility::vector1< numeric::xyzVector< core::Real > > s_vectors(6, numeric::xyzVector< core::Real >(0, 0, 0));
	// create empty RT instance
	core::kinematics::RT RT_temp;
	// create empty homogeneous transformation matrix instance
	numeric::HomogeneousTransform<core::Real> H_temp;

	// STEP 1: extract key atoms from these residues. There are always six key_atom_IDs in each module, but only the
	// first 'number_dofs_' are associated to DoFs, while the rest of the key_atom_IDs represent constraints
	utility::vector1< core::id::AtomID > key_atom_IDs;
	for ( core::Size i = 1; i <= 3; ++i ) {
		key_atom_IDs.push_back(core::id::AtomID(conformation.residue(res_numbers_[i]).mainchain_atom(2), res_numbers_[i]));
		key_atom_IDs.push_back(core::id::AtomID(conformation.residue(res_numbers_[i]).mainchain_atom(3), res_numbers_[i]));
	}


	// STEP 2: define the homogeneous matrices that ares used to define the axis vectors and arm vectors
	// start with homogeneous matrix between global reference atom and first atom of the key_atoms
	RT_temp = core::kinematics::RT(conformation.atom_tree().atom(ref_atom_ID_).get_stub(),
		conformation.atom_tree().atom(key_atom_IDs[1]).get_stub());
	Hmatrices.push_back( numeric::HomogeneousTransform<core::Real>(RT_temp.get_rotation(), RT_temp.get_translation()) );

	// update homogeneous matrices between all key atoms
	for ( core::Size i=1; i < key_atom_IDs.size(); ++i ) {
		// homogeneous matrix between two subsequent key atoms
		RT_temp = core::kinematics::RT(conformation.atom_tree().atom(key_atom_IDs[i]).get_stub(),
			conformation.atom_tree().atom(key_atom_IDs[i+1]).get_stub());
		Hmatrices.push_back ( numeric::HomogeneousTransform<core::Real>(RT_temp.get_rotation(), RT_temp.get_translation()) );
	}
	// check that size of Hmatrices is again 6, as expected
	runtime_assert_string_msg(Hmatrices.size() == 6, "attempt to update vectors while number of key atoms is (no longer) exactly six");

	// STEP 3: use the homogeneous matrices to construct the position vectors of the twist axes (r_vectors), as well as
	// axis vectors of the twist axes, expressed in global reference frame
	for ( core::Size i=1; i <= 6; ++i ) {
		// update H_temp to include transformation to next key atom (CA or C, i.e. atoms associated with free dihedral angles)
		H_temp = H_temp * Hmatrices[i];
		// define screw vector 's', which is local torison axis projected on defined global reference frame, i.e. R * [1;0;0]
		s_vectors[i] = H_temp.rotation_matrix().col_x();
		// moment arm from global reference frame to the CA-atom corresponds to the 'point' of the homogeneous matrix
		r_vectors[i] = H_temp.point();
	}

	// store vectors in struct
	screw_construct_struct screw_constructs;
	screw_constructs.s_vectors = s_vectors;
	screw_constructs.r_vectors = r_vectors;
	screw_constructs.Hmatrices = Hmatrices;

	return screw_constructs;
}

utility::vector1< ModuleType1::Screw >
ModuleType1::get_orthogonal_wrenches_from_twists(Screw const & twist_1,
	Screw const & twist_2,
	numeric::HomogeneousTransform<core::Real> const & H_CA0_CAa,
	numeric::HomogeneousTransform<core::Real> const & H_CAa_Ca,
	numeric::HomogeneousTransform<core::Real> const & H_Ca_CAb,
	numeric::HomogeneousTransform<core::Real> const & H_CAb_Cb){

	// STEP 1: prepare additional support variables
	// supporting homogeneous matrices to transform between various atom reference frames
	numeric::HomogeneousTransform< core::Real > const H_CAa_CAb = H_CAa_Ca * H_Ca_CAb;
	numeric::HomogeneousTransform< core::Real > const H_CA0_CAb = H_CA0_CAa * H_CAa_CAb;
	numeric::HomogeneousTransform< core::Real > const H_Ca_Cb = H_Ca_CAb * H_CAb_Cb;

	// STEP 2: define the screw vectors and moment arm vectors for the orthogonal screw axes
	// The first orthogonal screw axis is pointing from CA atom 'a' to CA atom 'b'
	numeric::xyzVector< core::Real > s_CAa_CAb = H_CAa_CAb.point(); // a relative vector that cannot be zero..
	s_CAa_CAb = H_CA0_CAa.rotation_matrix() * s_CAa_CAb; // ...that is expressed in IRF ...
	s_CAa_CAb /= s_CAa_CAb.norm(); // ...and subsequently normalized

	// The second orthogonal screw axis lies in the plane defined by the axes phi_a and psi_a, and intersects both the
	// phi_b and psi_b axes. First, the intersection point of axis phi_b with the plane is determined. The matrix
	// H_Ca_CAb transforms a vector from the ref frame of the C-atom of the b-input-residue to the reference frame connected
	// to the C-atom of the a-input-residue. If a vector [cx_phi_b; 0; 0] is defined in the former ref frame, then the projection
	// on the former ref frame has a component H_Ca_CAb.xz()*cx_phi_b + H_Ca_CAb.pz(). Setting this to zero, results in
	// below equation. However, first check that the Z-component is non-zero and artificially move it if so.
	core::Real cx_phi_b;
	if ( std::abs(H_Ca_CAb.xz()) > pos_singular_bound_ ) {
		cx_phi_b = -H_Ca_CAb.pz() / H_Ca_CAb.xz();
	} else {
		cx_phi_b = -H_Ca_CAb.pz() / (numeric::sign( H_Ca_CAb.xz() ) * pos_singular_bound_); // ensure that denominator is not singular
		TR.Warning << "Jacobian module containing residue " << res_numbers_[1] << " is near a singularity: atom temporarily moved " << pos_singular_bound_ << " Ang. to continue analysis." << std::endl;
	}
	// Above constant is used to construct a homogeneous transform, which goes from CA-atom of the b-residue to this construction point 'Pc1'
	numeric::HomogeneousTransform<core::Real> const H_CAb_Pc1(numeric::xyzMatrix<core::Real>::identity(),
		numeric::xyzVector<core::Real>(cx_phi_b,0,0));

	// Similarly, intersection point of axis psi_b with plane phi_a - psi_a expressed in CAb reference frame is labeled Pc2,
	// where Pc2 stand for construction point 2
	core::Real cx_psi_b;
	if ( std::abs(H_Ca_Cb.xz()) > pos_singular_bound_ ) {
		cx_psi_b = -H_Ca_Cb.pz() / H_Ca_Cb.xz();
	} else {
		cx_psi_b = -H_Ca_Cb.pz() / (numeric::sign( H_Ca_Cb.xz() ) * pos_singular_bound_); // ensure that denominator is not singular
		TR.Warning << "Jacobian module containing residue " << res_numbers_[1] << " is near a singularity: atom temporarily moved " << pos_singular_bound_ << " Ang. to continue analysis." << std::endl;
	}
	// Above constant is used to construct a homogeneous transform, which goes from C-atom of the b-residue to this construction point 'Pc2'
	numeric::HomogeneousTransform<core::Real> const H_Cb_Pc2(numeric::xyzMatrix<core::Real>::identity(),
		numeric::xyzVector<core::Real>(cx_psi_b,0,0));

	// Now, the screw axis pointing from Pc1 to Pc2 can be defined as
	numeric::xyzVector< core::Real > s_Pc1_Pc2 = ( H_CAb_Pc1.inverse() * H_CAb_Cb * H_Cb_Pc2 ).point();

	// Calculation of the arm vectors (from CA-atom of residue b to Pc1) acting on the vector s_Pc1_Pc2
	numeric::xyzVector< core::Real > const r_CAb_Pc1 = H_CAb_Pc1.point();

	// Express all vectors in the global reference frame
	numeric::HomogeneousTransform<core::Real> const H_CA0_Pc1 = H_CA0_CAb * H_CAb_Pc1;
	s_Pc1_Pc2 = H_CA0_Pc1.rotation_matrix() * s_Pc1_Pc2; // rotate the s_Pc1_Pc2 vector
	numeric::xyzVector< core::Real > const r_CA0_Pc1 = H_CA0_CAb * r_CAb_Pc1;
	// normalize screw axis vector s_Pc1_Pc2
	s_Pc1_Pc2 /= s_Pc1_Pc2.norm();

	// STEP 3: use the defined vectors to construct a set of two wrenches that are orthogonal to the other four wrenches
	// of the Jacobian module
	numeric::xyzVector<core::Real> const wrench_1_rot = s_CAa_CAb.cross_product(-H_CA0_CAa.point());
	numeric::xyzVector<core::Real> const wrench_1_lin = s_CAa_CAb;
	Screw const Wa_1(wrench_1_rot, wrench_1_lin);
	numeric::xyzVector<core::Real> const wrench_2_rot = s_Pc1_Pc2.cross_product(-r_CA0_Pc1);
	numeric::xyzVector<core::Real> const wrench_2_lin = s_Pc1_Pc2;
	Screw const Wa_2(wrench_2_rot, wrench_2_lin);

	// calculate the required ratio of these wrenches so that each wrench acts on only a single twist
	//utility::vector1 < numeric::xyzVector<core::Real> >
	core::Real const Wa_ratio1 = -(Wa_2.first.dot(twist_2.first) + Wa_2.second.dot(twist_2.second) ) /
		(Wa_1.first.dot(twist_2.first) + Wa_1.second.dot(twist_2.second) );
	core::Real const Wa_ratio2 = -(Wa_2.first.dot(twist_1.first) + Wa_2.second.dot(twist_1.second) ) /
		(Wa_1.first.dot(twist_1.first) + Wa_1.second.dot(twist_1.second) );

	// calculate the ratio to normalize based on the linear (force) vector
	core::Real const inv_norm_Wa1 = 1 / (Wa_ratio1 * Wa_1.second + Wa_2.second).norm();
	core::Real const inv_norm_Wa2 = 1 / (Wa_ratio2 * Wa_1.second + Wa_2.second).norm();

	// combine these calculations to construct the two wrenches of actuation, in disassembled format in order to allow them to be returned
	utility::vector1< Screw > wrench_set(2);
	wrench_set[1] = Screw (inv_norm_Wa1 * (Wa_ratio1 * Wa_1.first + Wa_2.first ), inv_norm_Wa1 * (Wa_ratio1 * Wa_1.second + Wa_2.second ) ); // linear part of first wrench
	wrench_set[2] = Screw (inv_norm_Wa2 * (Wa_ratio2 * Wa_1.first + Wa_2.first ), inv_norm_Wa2 * (Wa_ratio2 * Wa_1.second + Wa_2.second ) ); // linear part of 2nd wrench

	return(wrench_set);
}

} //jacobian
} //kinematics
} //core
