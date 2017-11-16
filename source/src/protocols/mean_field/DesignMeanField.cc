// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    DesignMeanField.cc

/// @brief   Method definitions for DesignMeanField.
/// @author  arubenstein

// Unit headers
#include <protocols/mean_field/DesignMeanField.hh>

// Package headers
#include <protocols/mean_field/AAMatrix.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>
// Numeric headers


// C++ headers
#include <iostream>


// Construct tracer.
static basic::Tracer TR("protocols.mean_field.DesignMeanField");


namespace protocols {
namespace mean_field {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
DesignMeanField::DesignMeanField( Size option,
	pose::PoseOPs & poses,
	utility::vector1 < pack::task::PackerTaskOP > tasks,
	scoring::ScoreFunctionOP scfxn
) :
	MeanField::MeanField( option, poses, tasks, scfxn ),
	aa_matrix_( new AAMatrix() )
{}

DesignMeanField::DesignMeanField( Size option,
	pose::PoseOPs & poses,
	utility::vector1 < pack::task::PackerTaskOP > tasks,
	scoring::ScoreFunctionOP scfxn,
	Real lambda_mem,
	Real tolerance,
	Real temp,
	Real threshold ) :
	MeanField::MeanField( option, poses, tasks, scfxn, lambda_mem, tolerance, temp, threshold ),
	aa_matrix_( new AAMatrix() )
{}

// Destructor
DesignMeanField::~DesignMeanField() {}


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods

void
DesignMeanField::show( std::ostream & output ) const
{
	MeanField::show( output );
	output << "Specificity Profile (AAMatrix): " << std::endl;
	aa_matrix_->show( output );
}

void
DesignMeanField::process()
{
	MeanField::process();
	init_aa_matrix();
}

AAMatrixCOP
DesignMeanField::aa_matrix() const
{
	return aa_matrix_;
}

// Accessors/Mutators



// Private methods /////////////////////////////////////////////////////////////

// Friend methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that DesignMeanField can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, DesignMeanField const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

/// @details initalize a AAMatrix from the rot_matrix
void
DesignMeanField::init_aa_matrix()
{
	aa_matrix_->build_aa_matrix( *( rot_matrix() ), energies_matrix(), temperature() );
}


}  // namespace mean_field
}  // namespace protocols
