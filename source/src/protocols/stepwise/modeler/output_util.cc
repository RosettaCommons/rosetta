// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/output_util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/output_util.hh>
#include <ObjexxFCL/format.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.output_util" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {

//////////////////////////////////////////////////////////////////////////////////////////////////////
void
output_boolean( std::string const & tag, bool boolean, std::ostream & outstream /* = std::cout */ ){

	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;
	outstream << tag;

	if ( boolean ) {
		outstream << A( 4, "T" );
	} else {
		outstream << A( 4, "F" );
	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
output_boolean( bool boolean, std::ostream & outstream /* = std::cout */ ){

	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;

	if ( boolean ) {
		outstream << A( 4, "T" );
	} else {
		outstream << A( 4, "F" );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void
output_rna_movemap_header( Size const & spacing, std::ostream & outstream  ){
	using namespace ObjexxFCL::format;
	outstream << A( spacing, "res_num" ) << A( spacing, " alpha  " ) << A( spacing, "  beta  " ) << A( spacing, " gamma  " );
	outstream << A( spacing, " delta  " ) << A( spacing, "epsilon " ) << A( spacing, "  zeta  " ) << A( spacing, " chi_1  " );
	outstream << A( spacing, "  nu_2  " ) << A( spacing, "  nu_1  " ) << A( spacing, "chi_O2' " ) << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void
output_protein_movemap_header( Size const & spacing, std::ostream & outstream  ){
	using namespace ObjexxFCL::format;
	outstream << A( spacing, "res_num" ) << A( spacing, " phi    " ) << A( spacing, "  psi   " ) << A( spacing, " omega   " );
	outstream << A( spacing, " chi_1 " ) << A( spacing, "  chi_2 " ) << A( spacing, " chi_3  " ) << A( spacing, " chi_4  " ) << std::endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
void
output_movemap( kinematics::MoveMap const & mm, core::pose::Pose const & pose, std::ostream & outstream  ){

	using namespace ObjexxFCL::format;
	using namespace core::id;
	using namespace core::kinematics;

	Size const total_residue = pose.total_residue();

	Size spacing = 10;

	outstream << "--------------------------------------------------------------------------------------" << std::endl;
	outstream << "Movemap ( in term of partial_pose seq_num ): " << std::endl;
	bool is_protein( false ), is_rna( false );
	for ( Size n = 1; n <= total_residue; n++ ) {
		if ( !is_protein && pose.residue_type( n ).is_protein() ) {
			output_protein_movemap_header( spacing, outstream ); is_rna = false; is_protein = true;
		} else  if ( !is_rna && pose.residue_type( n ).is_RNA() ) {
			output_rna_movemap_header( spacing, outstream ); is_rna = true; is_protein = false;
		}
		if ( is_rna ) {
			outstream << I( spacing, 3, n );
			for ( Size k = 1; k <= 6; k++ ) {
				outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::BB,  k ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			}
			for ( Size k = 1; k <= 4; k++ ) {
				outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::CHI, k ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			}
			outstream << std::endl;
		}
		if ( is_protein ) {
			outstream << I( spacing, 3, n );
			for ( Size k = 1; k <= 3; k++ ) {
				outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::BB,  k ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			}
			for ( Size k = 1; k <= pose.residue_type( n ).nchi(); k++ ) {
				outstream << A(  - 1 + ( spacing - 4 )/2, "" ); output_boolean( mm.get( TorsionID( n, id::CHI, k ) ), outstream ); outstream << A( 1 + ( spacing - 4 )/2, "" );
			}
			outstream << std::endl;
		}
	}
	outstream << "--------------------------------------------------------------------------------------" << std::endl;

	outstream << "print movemap jump_points [explicit method]: " << std::endl;
	for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ) {
		Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );

		outstream << "n = " << n << " | jump_pos1 = " << jump_pos1 << " | jump_pos2 = " << jump_pos2;
		outstream << " | mm.get_jump( n ) = "; output_boolean( mm.get_jump( n ), outstream );
		outstream << std::endl;

	}
	outstream << "--------------------------------------------------------------------------------------" << std::endl;

	//From core/kinematic/MoveMap.hh
	//typedef std::map< id::JumpID, bool > JumpID_Map
	outstream << "print movemap jump_points [iterator method]: " << std::endl;
	for ( std::map< id::JumpID, bool > ::const_iterator it = mm.jump_id_begin(), end = mm.jump_id_end(); it != end; ++it ) {
		outstream << "movemap jump == true for jump_pos1 = " << it->first << " | jump_pos2 = " << it->second << std::endl;
	}
	outstream << "--------------------------------------------------------------------------------------" << std::endl;
}


} //modeler
} //stepwise
} //protocols
