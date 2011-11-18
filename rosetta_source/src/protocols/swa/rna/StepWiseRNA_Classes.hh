// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Classes.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_SWA_RNAClasses_HH
#define INCLUDED_protocols_swa_SWA_RNAClasses_HH


#include <core/pose/Pose.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <core/id/TorsionID.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/GreenPacker.fwd.hh>
#include <string>
#include <map>


typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace swa {
namespace rna {

enum PuckerState{ ALL, NORTH, SOUTH };
enum BaseState{ BOTH, ANTI, NONE };

class Jump_point{

		public:

		Jump_point():
		five_prime_seq_num( 0 ),
		five_prime_atom("O3*"),
		three_prime_seq_num( 0 ),
		three_prime_atom("P"),
		cut_point( 0 )
	{
	}

	~Jump_point(){};

	public:

	core::Size five_prime_seq_num;
	std::string five_prime_atom;
	core::Size three_prime_seq_num;
	std::string three_prime_atom;
	core::Size cut_point; //Choose the residue five_prime of the actual cutpoint position

};

// Why don't we just keep this in a "regular" kinematics::Stub?
struct base_stub{
  numeric::xyzVector<core::Real> centroid;
  Matrix base_coordinate_matrix;
};

struct pose_data_struct2{
	core::Real score;
	core::pose::PoseOP pose_OP;
	std::string tag;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Torsion_Info{
	core::id::TorsionID id;
	core::Real value;
};

/*
class Torsion_Info{
	public:
	 Torsion_Info():
			value(0.0)
		{
		}
		
		~Torsion_Info(){};

	core::id::TorsionID id;
	core::Real value;

};
*/



struct output_data_struct {

	core::Real rmsd_wrt_correct_minimized;
	core::Real rmsd_wrt_minimized;
	core::Real rmsd;
	core::Real loop_rmsd_wrt_correct_minimized;
	core::Real loop_rmsd_wrt_minimized;
	core::Real loop_rmsd;
	core::Real diff_torsions;
	core::Real current_score;
	core::Real O3_C5_distance;
};

class SillyCountStruct{

public:

  SillyCountStruct():
    output_pose_count( 0 ),
    good_rep_rotamer_count( 0 ),
    good_atr_rotamer_count( 0 ),
    good_angle_count( 0 ),
    good_distance_count( 0 ),
    chain_closable_count( 0 ),
    Near_Native_Rotamer_Count( 0 ),
    base_pairing_count( 0 ),
    base_stack_count( 0 ),
    both_count( 0 ),
    tot_rotamer_count( 0 ),
    fine_rmsd_count( 0 ),
    rmsd_count( 0 )
    {
  }


  ~SillyCountStruct(){};

public:

  core::Size output_pose_count;
  core::Size good_rep_rotamer_count;
  core::Size good_atr_rotamer_count;
  core::Size good_angle_count;
  core::Size good_distance_count;
  core::Size chain_closable_count;
  core::Size Near_Native_Rotamer_Count;
  core::Size base_pairing_count;
  core::Size base_stack_count;
  core::Size both_count;
  core::Size tot_rotamer_count;
  core::Size fine_rmsd_count;
  core::Size rmsd_count;
};

////////////////////////////////////////////////////////////
}
}
}

#endif
