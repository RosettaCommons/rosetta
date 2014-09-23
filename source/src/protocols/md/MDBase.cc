// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/md/MDBase.cc
/// @brief   initialization for MD
/// @detailed
/// @author  Hahnbeom Park

#include <protocols/md/MDBase.hh>

#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Option
#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

// Constraints
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>

// parsing
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <fstream>

namespace protocols {
namespace md {

static thread_local basic::Tracer TR( "protocols.md" );

void
MDBase::set_constraint(	Real const sdev )
{
	constrained_ = true;
	cst_sdev_ = sdev;
	// starting coordinate constraint
	scorefxn_->set_weight( core::scoring::coordinate_constraint, 1.0 );
}

void
MDBase::cst_on_pose( pose::Pose &pose )
{
	// Remove all the constraints first
	pose.remove_constraints();

	// First, add cst_file info into pose
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file ].user() ) {
		TR << "Set constraints from input file..." << std::endl;
		scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose );
	}

	// Next, set coordinate constraint if specified
	if( constrained_ ){
		TR << "Set constraint on starting structure with sdev " << cst_sdev_ << std::endl;
		scoring::func::FuncOP fx( new scoring::func::HarmonicFunc( 0.0, cst_sdev_ ) );

		for( Size i_res = 1; i_res <= pose.total_residue(); ++i_res ){
			Size i_ca = pose.residue(i_res).atom_index(" CA ");
			id::AtomID atomID( i_ca, i_res );

			pose.add_constraint( new scoring::constraints::CoordinateConstraint(
				atomID, atomID, pose.residue(i_res).xyz(i_ca), fx
			) );
		}
	}
}

void
MDBase::parse_schfile( std::string const schfile ){

	std::vector< std::string > filelines;
	std::string line;

	std::ifstream infile( schfile.c_str() );
	TR.Debug << "================== Reading schedule file: ==================" << std::endl;
	if (!infile.good()) {
		utility_exit_with_message( "[ERROR] Error opening script file '" + schfile + "'" );
	}
	while( getline(infile,line) ) {
		filelines.push_back( line );
	}
	infile.close();

	// Clean first
	mdsch_.resize( 0 );

	// Put in by reading line by line
	for( Size i = 0; i < filelines.size(); i++ ){
		line = filelines[i];
		TR.Debug << line << std::endl;
		utility::vector1< std::string > tokens ( utility::split( line ) );
		// Format: sch nstep temp0
		if( tokens[1].compare("sch") == 0 ){
			MDscheduleData sch;
			sch.temp0 = atof( tokens[2].c_str() );
			sch.nstep = atoi( tokens[3].c_str() );
			mdsch_.push_back( sch );
		}
	}

}

}
}
