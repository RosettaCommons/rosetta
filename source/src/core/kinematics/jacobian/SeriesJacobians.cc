// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/kinematics/jacobian/SeriesJacobians.cc
/// @brief class that handles the Jacobian analysis of atoms connected in series
/// @author teunhoevenaars (teunhoevenaars@gmail.com)

/// @details This class handles the Jacobian analysis of a single series of atoms. A full Jacobian analysis consists of
/// one or more sets of such series, as described in @ref jacobian_structure_class. A series of atoms (represented
/// by a series of residues) in a protein can have any length, and therefore any number of internal Degrees of Freedom (DoFs)
/// However, the number of Cartesian DoFs of the final atom in the chain are by definition 6. Therefore, the maximum number
/// of independent internal DoFs to describe six differential Cartesian DoFs is 6. Jacobian analysis relies on orthogonal
/// vectors, which means that a Jacobian analysis can only be done for blocks of maximum six internal DoFs. To deal with
/// this, the SeriesJacobians class splits up a series of atoms (taken from input residues) into blocks that have six internal
/// DoFs.

// Project headers:
#include <core/kinematics/jacobian/SeriesJacobians.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

// Class headers
#include <core/kinematics/jacobian/ModuleType1.hh>
#include <core/conformation/Conformation.hh>

static basic::Tracer TR( "core.kinematics.jacobian.SeriesJacobians" );


namespace core {
namespace kinematics {
namespace jacobian {

/// @brief Constructor based on a vector containing a series of residues
/// @details This function constructs SeriesJacobians object based on a series of residues. First it is determined what
/// kind of residues the series consists of, because this determines which Jacobian modules need to be initialized.
/// @param[in] residue_set the provided vector with residue numbers needs to be in ascending order, but may contain gaps (e.g. a
/// vector with residues [16, 17, 19, 23, 24, 25, 26] is fine).
SeriesJacobians::SeriesJacobians(core::conformation::Conformation const & conformation, residue_series const & residue_set, core::id::AtomID const & ref_atom){
	using namespace core::kinematics::jacobian;

	residue_set_ = residue_set;
	ref_atom_ = ref_atom;

	// get residue series type (for now only single series)
	SeriesJacobianTypeEnum const residue_set_type = determine_residue_series_type( conformation, residue_set_ );

	// currently single option, but placeholder for extension to other types of residue series.
	if ( residue_set_type == SeriesJacobianTypeEnum::CANONICAL_AA ) {
		number_dofs_ = residue_set_.size() * 2;
		modules_ = init_modules_amino_acids();
	}
}

/// @brief Copy constructor.  Keep default unless deep copying is needed (and in that case,
/// consider using DeepCopyOPs.)
SeriesJacobians::SeriesJacobians(SeriesJacobians const & )=default;

/// @brief Destructor.
SeriesJacobians::~SeriesJacobians(){}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
SeriesJacobiansOP
SeriesJacobians::clone() const {
	return utility::pointer::make_shared< SeriesJacobians >(*this );
}

///@brief update all Jacobian matrices in the series chain
void
SeriesJacobians::update_Jacobian_matrices(core::conformation::Conformation const & conformation){
	// loop through modules and call Jacobian update functions
	for ( core::Size i=1; i <= modules_.size(); ++i ) {
		modules_[i]->update_jacobian_matrices(conformation);
	}
}

///@brief function that initializes the low-level Jacobian modules for a series of canonical amino acids
///@details For canonical amino acids the internal DoFs are the phi and psi torsion angles. That knowledge is used to
/// initialize the underlying Jacobian modules. First, the total number of DoFs in the residue series is calculated, and
/// subsequently the residue series is complemented with "assisting residues" such that the series can be split up in
/// modules of six torsion axes. This means that some of the torsion axes in the last module are actually constrained,
/// which is captured by the variable dofs_module.
utility::vector1< core::kinematics::jacobian::ModuleType1OP >
SeriesJacobians::init_modules_amino_acids(){
	// create empty vector
	utility::vector1< core::kinematics::jacobian::ModuleType1OP > modules;

	// determine modulus of free residues with 3, which indicates how many assisting residues need to be added to
	// get to a multiple of 3. Add residues that precede the first residue of the last module. E.g. if the last module
	// has one residue (41), then add 39 and 40. If the last module has residues 41 and 43, then add only residue 40.
	core::Size free_res_mod = residue_set_.size() % 3;
	while ( free_res_mod % 3 != 0 ) {
		residue_set_.push_back( residue_set_.end()[ - free_res_mod ] - 3 + free_res_mod ); // add residue to vector
		free_res_mod = residue_set_.size() % 3; // update modulus
	}

	// confirm that the number of residues in the vector (after adding assisting residues) is now a multiple of 3
	runtime_assert_string_msg(residue_set_.size() % 3 == 0, "initialization of Jacobian series failed to divide residue set in sets of 3 residues");

	// determine the number of modules
	core::Size const N_jac_modules = residue_set_.size()/3;

	// create empty vector of modules
	modules.reserve(N_jac_modules);

	// initialize individual members
	core::Size dofs_module_i = 6; // default number of dofs in each module

	for ( core::Size i=0; i < N_jac_modules; ++i ) {
		// Last modules not necessarily has 6 dofs
		if ( i == N_jac_modules-1 ) {
			// if it is >0 then the size of the last module is the modulus, otherwise it's the default 6
			if ( number_dofs_ % 6 != 0 ) {
				dofs_module_i = number_dofs_ % 6;
			}
		}
		// get residues that make up this module
		utility::vector1< core::Size > res_numbers_i(residue_set_.begin()+i*3, residue_set_.begin()+ (i+1)*3);

		// create new ModuleType1 pointer and add to vector
		modules.push_back(utility::pointer::make_shared< ModuleType1 >(dofs_module_i, res_numbers_i, ref_atom_) );
	}
	// verify that the expected number of Jacobian modules have been created
	runtime_assert_string_msg(N_jac_modules == modules.size(),"initialization of Jacobian series failed to create expected number of Jacobian modules");
	// verify that the residue number of the last atom of the last Jacobian module is the last residue in the provided input vector
	runtime_assert_string_msg(modules.back()->get_residues().back() == residue_set_.back(),"Last residue of last Jacobian module is not the same as last residue of input set");

	return modules;
}

SeriesJacobianTypeEnum
SeriesJacobians::determine_residue_series_type(core::conformation::Conformation const & conformation, residue_series const & res_numbers ){
	// N.B. currently, protocol can only work with canonical amino acids. Can be extended/generalized later.
	//      if function also needs to deal with series of residues of different types, then this could become complex
	for ( core::Size i=1; i <= res_numbers.size(); ++i ) {
		if ( !conformation.residue_type(res_numbers[i]).is_canonical_aa() ) {
			utility_exit_with_message("Non-canonical amino acid detected, but Jacobian structure analysis can currently only deal with canonical amino acids.");
		}
		//Later, put "return SeriesJacobianTypeEnum::NONCANONICAL_AA;" here.
	}

	// for now, if function completes, then the type is canonical_aa
	return SeriesJacobianTypeEnum::CANONICAL_AA;
}

} //jacobian
} //kinematics
} //core
