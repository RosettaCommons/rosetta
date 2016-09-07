// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
/// @author Chu Wang


// Unit Headers
#include <protocols/jumping/ResiduePairJumpSetup.hh>

// Package Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif

#include <core/fragment/OrderedFragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/SequenceMapping.hh>

#include <protocols/jumping/JumpSample.hh>
#include <protocols/jumping/ResiduePairJump.hh>

// Project Headers
#include <basic/Tracer.hh>

// Numeric headers

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/io/izstream.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.fwd.hh>


//numeric headers

// C++ headers


namespace protocols {
namespace jumping {

using namespace core;
using namespace fragment;

static THREAD_LOCAL basic::Tracer tr( "protocols.jumping" );

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResiduePairJumpSetup::read_file( std::string fname ) {


	utility::io::izstream data( fname.c_str() );
	tr.Info << "read ResiduePair jump-definitions from " << fname << std::endl;
	if ( !data ) {
		utility_exit_with_message( "Unable to open constraints file: " + fname +"\n");
	}

	core::chemical::ResidueTypeSetCOP residue_type_set(
		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
	);

	std::string line, tag;
	core::Size root = 1;
	ResiduePairJumpOP residue_pair_jump = nullptr;
	while ( getline(data, line ) ) {
		std::istringstream in( line );
		in >> tag;
		if ( tag == "BEGIN" ) {
			residue_pair_jump = ResiduePairJumpOP( new ResiduePairJump );
			continue;
		} else if ( tag == "END" ) {
			residue_pair_jump->init_mini_pose();
			add_residue_pair_jump( residue_pair_jump );
			residue_pair_jump = nullptr;
			continue;
		} else if ( tag == "ROOT" ) {
			in >> root;
		}

		if ( ! residue_pair_jump ) continue;

		if ( tag == "jump_def:" ) {
			Interval jump, cuts;
			in >> jump.start_ >> jump.end_ >> cuts.start_ >> cuts.end_;
			add_jump( jump, cuts );
		} else if ( tag == "aa:" ) {
			for ( int i = 1; i <= 2; ++i ) {
				std::string name;
				in >> name;
				if ( residue_type_set->name_map(name).is_protein() ) {
					name = name + ":CtermProteinFull:NtermProteinFull";
				}
				core::chemical::ResidueType const & res_type( residue_type_set->name_map(name) );
				residue_pair_jump->add_residue_single( res_type );
			}
		} else if ( tag == "cst_atoms:" ) {
			for ( int i = 1; i <= 2; ++i ) {
				std::string name;
				for ( int j = 1; j <= 3; ++j ) {
					in >> name;
					residue_pair_jump->set_cstAtoms(i,j,name);
				}
			}
		} else if ( tag == "jump_atoms:" ) {
			for ( int i = 1; i <= 2; ++i ) {
				std::string name;
				for ( int j = 1; j <= 3; ++j ) {
					in >> name;
					residue_pair_jump->set_jumpAtoms(i,j,name);
				}
			}
		} else if ( tag == "disAB:" ) {
			Real value;
			while ( in >> value ) {
				residue_pair_jump->set_cstInfo( disAB, value );
			}
		} else if ( tag == "angleA:" ) {
			Real value;
			while ( in >> value ) {
				residue_pair_jump->set_cstInfo( angleA, value );
			}
		} else if ( tag == "angleB:" ) {
			Real value;
			while ( in >> value ) {
				residue_pair_jump->set_cstInfo( angleB, value );
			}
		} else if ( tag == "dihedralA:" ) {
			Real value;
			while ( in >> value ) {
				residue_pair_jump->set_cstInfo( dihedralA, value );
			}
		} else if ( tag == "dihedralB:" ) {
			Real value;
			while ( in >> value ) {
				residue_pair_jump->set_cstInfo( dihedralB, value );
			}
		} else if ( tag == "dihedralAB:" ) {
			Real value;
			while ( in >> value ) {
				residue_pair_jump->set_cstInfo( dihedralAB, value );
			}
		}
	}
	set_root( root );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FragSetOP
ResiduePairJumpSetup::generate_jump_frags( JumpSample const& jumps, kinematics::MoveMap const& mm ) const {
	OrderedFragSetOP frags( new OrderedFragSet );
	FrameList jump_geometries;
	//runtime_assert( jumps.total_residue() == total_residue() );
	ObjexxFCL::FArray2D_int const & in_jumps ( jumps.jumps() );
	int ct = 1;
	for ( auto it=begin(), eit=end(); it!=eit; ++it, ct++ ) {
		Size jump_number = 0;
		for ( Size i = 1; i <= jumps.size(); ++i ) {
			if  ( ( in_jumps( 1, i ) == int( it->jump_.start_) ) && ( in_jumps( 2, i ) == int( it->jump_.end_ ) ) ) {
				jump_number = i;
				break;
			}
		}
		// did not find a match jump, so skip
		if ( jump_number == 0  ) continue;
		// this jump is not a flexible jump
		if ( ! mm.get_jump(jump_number) ) continue;
		// generate a frame based on this ResiduePairJump
		FrameOP jump_frame = ResiduePairJumps_[ct]->generate_frame();
		jump_frame->show( tr.Info );
		// map ResiduePairJump index 1 and 2 to the correct seqpos number defined in JumpSample
		core::id::SequenceMapping map;// 2, jumps.total_residue());
		map.insert_aligned_residue( 1, in_jumps(1,jump_number));
		map.insert_aligned_residue( 2, in_jumps(2,jump_number));
		// align the jump frame to the whole pose context so that they can be used in protocols level
		jump_frame->align( map );
		// save this frame to FrameList
		jump_geometries.push_back( jump_frame);
	}
	// add FrameList to FragSet
	frags->add( jump_geometries );

	return frags;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
JumpSample
ResiduePairJumpSetup::create_jump_sample() const
{
	ObjexxFCL::FArray2D_int jumps( 2, jumps_.size() );
	ObjexxFCL::FArray1D_int cuts( jumps_.size() );
	ObjexxFCL::FArray2D<std::string> jump_atoms(2, jumps_.size(),"");

	int ct = 1;
	Size total_residue = total_residue_;
	for ( auto it=begin(), eit=end(); it!=eit; ++it, ct++ ) {
		jumps( 1, ct ) = it->jump_.start_;
		jumps( 2, ct ) = it->jump_.end_;
		Size const crs ( it->cut_reg_.start_ );
		Size const cre ( it->cut_reg_.end_ );
		if ( crs > total_residue ) total_residue = crs;
		if ( cre > total_residue ) total_residue = cre;
		if ( it->jump_.end_ > total_residue ) total_residue = it->jump_.end_;
		if ( it->jump_.start_ > total_residue ) total_residue = it->jump_.start_;
		cuts( ct ) = crs+int( numeric::random::uniform()*( cre-crs ) + 0.5 );
		jump_atoms(1, ct) = ResiduePairJumps_[ct]->jumpAtoms(1)[1];
		jump_atoms(2, ct) = ResiduePairJumps_[ct]->jumpAtoms(2)[1];
	}
	return JumpSample( total_residue, jumps_.size(), jumps, jump_atoms, cuts, root() );
}

} //jumping
} //protocols
