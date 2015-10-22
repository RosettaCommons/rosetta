// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/conformation/parametric/ParametersSet.cc
/// @brief  A class for holding sets of parameters for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <core/conformation/parametric/ParametersSet.hh>

// Package headers
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/Conformation.hh>

// Project headers

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <iostream>


namespace core {
namespace conformation {
namespace parametric {

static THREAD_LOCAL basic::Tracer TR( "core.conformation.parametric.ParametersSet" );

/// @brief Constructor.
///
ParametersSet::ParametersSet() :
	parameters_()
{
}

/// @brief Copy constructor.
///
ParametersSet::ParametersSet( ParametersSet const & src ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< ParametersSet >()
{
	parameters_.clear();
	//Make copies of the Parameters objects:
	if ( src.n_parameters()>0 ) {
		for ( core::Size i=1, imax=src.parameters_.size(); i<=imax; ++i ) {
			parameters_.push_back( src.parameters(i)->clone() );
		}
	}
}

ParametersSet::~ParametersSet() {}


/// @brief make a copy of this residue( allocate actual memory for it )
///
ParametersSetOP
ParametersSet::clone() const
{
	return ParametersSetOP( new ParametersSet( *this ) );
}

/// @brief Only for copying Conformation objects, this ensures that the new ParametersSet object's
/// Parameters objects have lists of ResidueOPs that point to residues in the new Conformation object,
/// rather than to residues that only exist in the Parameters objects.
void ParametersSet::update_residue_links( core::conformation::Conformation &new_conf ) {
	core::Size const nparams=parameters_.size();
	if ( nparams==0 ) return; //Do nothing if this object contains no Parameters objects.
	for ( core::Size i=1; i<=nparams; ++i ) {
		core::Size nresidue = parameters_[i]->n_residue();
		if ( nresidue==0 ) continue; //Do nothing to this Parameters object if it has no residues in its residue list.
		for ( core::Size j=1; j<=nresidue; ++j ) {
			core::Size seq_pos=parameters_[i]->residue(j)->seqpos();
			debug_assert( seq_pos > 0 && seq_pos <= new_conf.size() );
			parameters_[i]->set_residue(j, new_conf.residue_op(seq_pos));
		}
	}
	return;
}

/// @brief Get a summary of this ParametersSet object, for output to remark lines of a PDB file.
/// @details Default function can be overridden by derived classes.  The default message just says that the
/// parameterization does not support output to PDB remark lines.  Note that messages must have lines less than
/// about fifty characters, or they get truncated in REMARK output.
void ParametersSet::get_pdb_remark(std::stringstream &remark) const {
	//                                                                    X Anything from the X onward gets truncated
	remark << "This pose is parametric, but the parametrization does not" << std::endl;
	remark << "currently support output to PDB remark lines."             << std::endl;
	return;
}

} // namespace parametric
} // namespace conformation
} // namespace core

