// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseBetaAntiParallelJumpSampleGenerator
/// @brief Subclass of StepWisePoseSampleGenerator
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/protein/beta/StepWiseBetaAntiParallelJumpSampleGenerator.hh>
#include <protocols/stepwise/sampling/protein/sample_generator/StepWiseProteinJumpSampleGenerator.hh>
#include <core/types.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/database/open.hh>

#include <utility/io/izstream.hh>

#include <iostream>
#include <string>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/io/mpistream.hh>



using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {
namespace beta {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
	StepWiseBetaAntiParallelJumpSampleGenerator::StepWiseBetaAntiParallelJumpSampleGenerator(
																		 core::pose::Pose const & pose,
																		 Size const moving_residue ):
		StepWiseProteinJumpSampleGenerator( 0, jumps_ ) //generic
	{
		Size which_jump = get_antiparallel_beta_jumps( pose, moving_residue );
		initialize( which_jump, jumps_ );
	}

	Size
	StepWiseBetaAntiParallelJumpSampleGenerator::get_antiparallel_beta_jumps( pose::Pose const & pose, int const sample_res ){

		using namespace core::kinematics;
		using ObjexxFCL::strip_whitespace;

		Size which_jump = 0;
		jumps_.clear();

		// Which jump connects into the new residue ("sample_res")? We're going to figure out nice beta-pairing jumps to it.
		bool downstream( false );
		FoldTree const & f( pose.fold_tree() );

		for ( Size i = 1; i <= f.num_jump(); i++ ) {
			if ( f.upstream_jump_residue( i )   ==  sample_res  ){
				downstream = true;
				which_jump = i;	 break;
			} else if (f.downstream_jump_residue( i ) == sample_res ){
				downstream = false;
				which_jump = i; break;
			}
		}
		if ( which_jump == 0 ) utility_exit_with_message( "Could not find the jump for sample_beta?!" );

		std::string atom_base, atom_sample;
		if ( downstream ){
			atom_base   = f.upstream_atom( which_jump );
			atom_sample = f.downstream_atom( which_jump );
		} else {
			atom_base   = f.downstream_atom( which_jump );
			atom_sample = f.upstream_atom( which_jump );
		}
		atom_base   = strip_whitespace( atom_base );
		atom_sample = strip_whitespace( atom_sample );

		std::string const jump_library_file( basic::database::full_name( "clustered_beta_pairs.dat" ) );
		utility::io::izstream data( jump_library_file.c_str() );
		if ( !data.good() ) utility_exit_with_message( "Unable to open file: " +jump_library_file + '\n' );

		std::string line, atom1_in, atom2_in, tag;
		Jump jump;
		while( getline(data, line) ) {
			std::istringstream is( line );
			is >> tag >> atom1_in >> atom2_in;
			is >> jump;

			std::cout << atom1_in << " " << atom_base << " " << atom2_in << " " << atom_sample << std::endl;

			if ( (atom1_in == atom_base) && (atom2_in == atom_sample ) ){
				if ( !downstream ) jump.reverse();
				jumps_.push_back( jump );
			}
		}
		data.close();

		std::cout << "FOUND " << jumps_.size() << " JUMPS!" << std::endl;

		return which_jump;

	}

} //beta
} //protein
} //sampling
} //stepwise
} //protocols

