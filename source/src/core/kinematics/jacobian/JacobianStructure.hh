// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/kinematics/jacobian/JacobianStructure.hh
/// @brief Class that defines structure of Jacobian modules
/// @author teunhoevenaars (teunhoevenaars@gmail.com)


#ifndef INCLUDED_core_kinematics_jacobian_JacobianStructure_hh
#define INCLUDED_core_kinematics_jacobian_JacobianStructure_hh

#include <core/kinematics/jacobian/JacobianStructure.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// protocol headers
#include <core/id/AtomID.hh>
#include <core/kinematics/jacobian/SeriesJacobians.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/types.hh>

namespace core {
namespace kinematics {
namespace jacobian {

/// @brief The JacobianStructure class is the upper-level wrapper of the Jacobian analysis of a protein's kinematics relations.
/// @author teunhoevenaars (teunhoevenaars@gmail.com)

class JacobianStructure : public utility::VirtualBase {

public:

	/// @brief No default constructor because JacobianStructure is currently not used on its own, but always as part of a mover
	JacobianStructure() = delete;

	/// @brief Constructor based on single chain of free residues.
	JacobianStructure( core::conformation::Conformation const & conformation, utility::vector1<core::Size> const & free_residues, core::id::AtomID const & ref_atom );

	/// @brief Copy constructor.
	JacobianStructure(JacobianStructure const & src);

	/// @brief Destructor.
	~JacobianStructure() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	JacobianStructureOP clone() const;

	///@brief get one serial chain of the structure
	core::kinematics::jacobian::SeriesJacobiansOP const
	get_single_chain( core::Size index){
		return serial_chains_[index];
	}

private: // VARIABLES

	/// @brief vector with Jacobian serial chains
	utility::vector1< core::kinematics::jacobian::SeriesJacobiansOP > serial_chains_;

	/// @brief atomID of the atom whose reference frame acts is reference frame for all vectors in the Jacobian structure
	core::id::AtomID ref_atom_;
};

} //jacobian
} //kinematics
} //core

#endif //INCLUDED_core_kinematics_jacobian_JacobianStructure_hh
