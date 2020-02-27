// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/kinematics/jacobian/JacobianStructure.cc
/// @brief class that defines structure of Jacobian modules
/// @author teunhoevenaars (teunhoevenaars@gmail.com)
/// @page jacobian_structure_class

/// @details A Jacobian is the first-order differential of a function that takes a vector as its inputs and outputs. In mechanism
/// analysis (and proteins can be considered complex mechanisms) Jacobian matrices are used to describe the tangent of
/// a conformation. A conformation relates Cartesian positions/orientations of atoms to internal coordinates (such as
/// torsion angles), which is also referred to as kinematics. The tangent (Jacobian) of a conformation is thus a first-order
/// approximation of how a change in torsion angles affects the positions/orientations of the various atoms.
///
/// The position/orientation of some atom with respect to some other atom dependss on all DoFs that lie between them.
/// In Rosetta the AtomTree uses Stubs to calculate the Cartesian position/orientation of the atom at the end of a series
/// of atoms. Here, JacobianStructure uses SeriesJacobians to calculate the tangent of the Cartesian position/orientation
/// of the atom at the end of a series of atoms with respect to the atom at the beginning of that series. Once expressed,
/// a Jacobian can be used to:
///
/// @li 1. iteratively solve the inverse kinematics problem, i.e. obtain the required torsion angles to achieve a specified
///     pose. This is for example of interest to loop closure.
/// @li 2. transform force/moment vectors that are expressed in a generalized space (e.g. as a function of bond angles, bond lengths, etc.)
///     into Cartesian force/moment vectors, and vice versa. This can be used to find the energy gradient, and thereby guide the
///     energy descent.
///
/// The first step in a Jacobian analysis is to obtain the kinematic relations. In most cases the atoms of a protein are
/// connected in series (the basis for the AtomTree), which leads to a single set of kinematic relations. However, in some
/// cases separate serial chains are connected in parallel (e.g. cyclical molecules, proline). In that case the Jacobian
/// analysis is a combination of sets of kinematic relations, one for each serial chain. This class handles that
/// upper-level structure that expresses how the Jacobian analysis of individual serial chains are related.

// Project headers:
#include <core/kinematics/jacobian/JacobianStructure.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

// Protocol headers:
#include <core/kinematics/jacobian/SeriesJacobians.hh>

static basic::Tracer TR( "core.kinematics.jacobian.JacobianStructure" );


namespace core {
namespace kinematics {
namespace jacobian {

/// @brief Constructor based on single chain of free residues.
/// @details In case of single serial chain atoms (provided as a residue series) the upper-level structure of a Jacobian
/// analysis is trivial, because it contains only one element.
JacobianStructure::JacobianStructure( core::conformation::Conformation const & conformation, utility::vector1<core::Size> const & free_residues, core::id::AtomID const & ref_atom ){

	// copy ref_atom
	ref_atom_ = ref_atom;
	// initialize Jacobian serial chain, which in turns initializes the underlying modules
	serial_chains_.push_back( utility::pointer::make_shared< core::kinematics::jacobian::SeriesJacobians >(conformation, free_residues, ref_atom_) );
}

/// @brief Copy constructor.  Keep default unless deep copying is needed (and in that case,
/// consider using DeepCopyOPs.)
JacobianStructure::JacobianStructure( JacobianStructure const & )=default;

/// @brief Destructor.
JacobianStructure::~JacobianStructure(){}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
JacobianStructureOP
JacobianStructure::clone() const {
	return utility::pointer::make_shared< JacobianStructure >( *this );
}

} //jacobian
} //kinematics
} //core
