// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/forge/constraints/ConstraintFileRCG.cc
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009

// Unit header
#include <protocols/forge/constraints/ConstraintFileRCG.hh>

// Package headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/id/SequenceMapping.hh>

// Project headers
// AUTO-REMOVED #include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.forge.constraints.ConstraintFileRCG" );

namespace protocols{
namespace forge{
namespace constraints{

/// @brief
//ConstraintFileRCG::ConstraintFileRCG():
//	RemodelConstraintGenerator(),
//	filename_( "" )
//{
//	using namespace basic::options;
//	using namespace basic::options::OptionKeys;
//	if( option[ constraints::cst_file ].user() ){
//		filename_ = option[ constraints::cst_file ].value().at( 1 );
//	}
//}

/// @brief
ConstraintFileRCG::ConstraintFileRCG( String const & filename ):
	RemodelConstraintGenerator(),
	filename_( filename )
{}

/// @brief
ConstraintFileRCG::~ConstraintFileRCG() {}

/// @brief
void
ConstraintFileRCG::set_cstfile( String const & filename )
{
	filename_ = filename;
}

/// @brief
void
ConstraintFileRCG::generate_remodel_constraints( Pose const & pose )
{
	using namespace core::scoring::constraints;
	ConstraintSetOP constraints = ConstraintIO::get_instance()->read_constraints( filename_, new ConstraintSet, pose );
	this->add_constraints( constraints->get_all_constraints() );
} //generate constraints


} //namespace constraints
} //namespace forge
} //namespace protocols
