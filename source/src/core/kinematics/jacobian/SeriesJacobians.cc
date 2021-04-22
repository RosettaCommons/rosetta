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
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/TorsionID.hh>

static basic::Tracer TR( "core.kinematics.jacobian.SeriesJacobians" );


namespace core {
namespace kinematics {
namespace jacobian {

/// @brief Constructor based on a vector containing a series of residues
/// @details This function constructs SeriesJacobians object based on a series of residues. First it is determined what
/// kind of residues the series consists of, because this determines which Jacobian modules need to be initialized.
/// @param[in] residue_set the provided vector with residue numbers needs to be in ascending order, but may contain gaps (e.g. a
/// vector with residues [16, 17, 19, 23, 24, 25, 26] is fine).
SeriesJacobians::SeriesJacobians(core::conformation::Conformation const & conformation, residue_series const & residue_set, core::id::AtomID const & ref_atom_ID){
	using namespace core::kinematics::jacobian;

	residue_set_ = residue_set;
	ref_atom_ID_ = ref_atom_ID;

	// get residue series type (for now only single series)
	SeriesJacobianTypeEnum const residue_set_type = determine_residue_series_type( conformation, residue_set_ );

	// currently single option, but placeholder for extension to other types of residue series.
	if ( residue_set_type == SeriesJacobianTypeEnum::ALPHA_AA ) {
		number_dofs_ = residue_set_.size() * 2; // each residue has a phi and psi torsion angle
		// initialize the modules
		modules_ = init_modules_amino_acids();
	} else {
		utility_exit_with_message("Attempted to perform Jacobian analysis on residues that are not Alpha_AA residues. "
			"This is currently not supported.");
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

///@brief obtain vector of structs with all Jacobian matrices in the series chain
utility::vector1< ModuleType1::jacobian_struct >
SeriesJacobians::get_Jacobian_matrices(core::conformation::Conformation const & conformation) const{
	// loop through modules and call Jacobian update functions
	utility::vector1< ModuleType1::jacobian_struct > jacobian_structs(modules_.size());
	for ( core::Size i=1; i <= modules_.size(); ++i ) {
		jacobian_structs[i] = modules_[i]->get_jacobian_matrices(conformation);
	}

	return jacobian_structs;
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

	// determine variables that describes number of residue per module and number of free torsions per residue
	core::Size const N_res_per_module = 3;
	core::Size const N_dofs_per_res = 2;

	// create copy of residue_set for local use and manipulation
	utility::vector1< core::Size > residue_set_tmp = residue_set_;

	// determine modulus of residues that cannot be fitted in full module. 'left-over' residues have omega fixed,
	// which is necessary to create a module with sets of collocated axes that are fully contained within residues (i.e. phi and psi axes)
	// Assisting residues (to represent constraints of this final module) are taken as residues preceding the last module. E.g. if the last module
	// has one residue (41), then add 39 and 40 before entry 41. If the last module has residues 41 and 43, then add only residue 40.
	core::Size res_remain = residue_set_tmp.size() % N_res_per_module;
	while ( res_remain % 3 != 0 ) {
		residue_set_tmp.push_back( residue_set_tmp.end()[-res_remain] - 3 + res_remain ); // add residue to vector
		++res_remain; // increment remaining residues
	}

	// determine the number of modules
	core::Size const N_jac_modules = residue_set_tmp.size() / N_res_per_module; // using Size to round down inexact division

	// initialize individual members (modules) of the series as well as the move_map
	core::Size dofs_module_i = 6; // default number of dofs in each module
	utility::vector1< core::Size >::iterator res_module_i_start = residue_set_tmp.begin();
	core::kinematics::MoveMapOP mm(utility::pointer::make_shared<core::kinematics::MoveMap>());
	for ( core::Size i=0; i < N_jac_modules; ++i ) {
		// Last modules not necessarily has 6 dofs
		if ( i == N_jac_modules-1 ) {
			// if last module has <6 dofs, change the settings for that module
			if ( number_dofs_ % 6 != 0 ) {
				dofs_module_i = number_dofs_ % 6;
			}
		}
		// get residues that make up this module
		utility::vector1< core::Size > res_numbers_i(res_module_i_start,res_module_i_start + N_res_per_module);
		res_module_i_start += N_res_per_module; // update iterator to reflect residues accounted for

		// create vector of torsion IDs for this module and also create the related movemap
		utility::vector1< core::id::TorsionID > torsion_ids_i;
		// reset dof_counter
		core::Size dof_counter = 1;
		for ( core::Size j=1; j <= N_res_per_module; ++j ) { // loop over residues i module
			for ( core::Size k = 1; k <= N_dofs_per_res; ++k ) { // loop over phi and psi angles
				// add torsion id to this modules vector
				torsion_ids_i.push_back(core::id::TorsionID(res_numbers_i[j], core::id::TorsionType::BB, k));
				// add torsion to movemap if it is a free dof
				if ( dof_counter <= dofs_module_i ) {
					mm->set(torsion_ids_i.back(), true);
				}
				++dof_counter;
			}
		}

		// store the movemap
		move_map_ = mm;

		// create new ModuleType1 pointer and add to vector
		modules.push_back(utility::pointer::make_shared< ModuleType1 >(dofs_module_i, torsion_ids_i, ref_atom_ID_) );
	}
	// verify that the residue number of the last atom of the last Jacobian module is the last residue in the provided input vector
	runtime_assert_string_msg(modules.back()->get_torsion_ids().back().rsd() == residue_set_tmp.back(),"Residue of last torsionID of last Jacobian module is not the same as last torsionID of input set");

	return modules;
}

SeriesJacobianTypeEnum
SeriesJacobians::determine_residue_series_type(core::conformation::Conformation const & conformation, residue_series const & res_numbers ){
	// N.B. currently, protocol can only work with canonical amino acids. Can be extended/generalized later.
	//      if function also needs to deal with series of residues of different types, then this could become complex
	for ( core::Size i=1; i <= res_numbers.size(); ++i ) {
		if ( !conformation.residue_type(res_numbers[i]).is_alpha_aa() ) {
			utility_exit_with_message("Detected amino acid that is not an alpha amino acid, but Jacobian structure analysis can currently only deal with alpha amino acids.");
		}
		//Later, put other options here if applicable.
	}

	// for now, if function completes, then the type is canonical_aa
	return SeriesJacobianTypeEnum::ALPHA_AA;
}

} //jacobian
} //kinematics
} //core
