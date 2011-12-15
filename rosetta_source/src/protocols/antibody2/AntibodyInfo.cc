// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)

// Rosetta Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/id/AtomID_Map.hh>

#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/loops/Loops.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>


// AUTO-REMOVED #include <utility/vector1.hh>

static basic::Tracer TR("antibody2.AntibodyInfo");

namespace protocols{
namespace antibody2{

/// default constructor
AntibodyInfo::AntibodyInfo() {
	set_default( false );

	for( core::Size i = 0; i <= 6; i++ )
		hfr_[i][0] = hfr_[i][1] = hfr_[i][2] = 0;
}


/// constructor with arguments
AntibodyInfo::AntibodyInfo( core::pose::Pose & pose ) {
	setup_loops( pose, false );
}

/// constructor with arguments
AntibodyInfo::AntibodyInfo( core::pose::Pose & pose, bool camelid ) {
	setup_loops( pose, camelid );
}

/// constructor with arguments
AntibodyInfo::AntibodyInfo( core::pose::Pose & pose, std::string cdr_name )
{
	set_default( false );

	if( !camelid_ ) {
		if( cdr_name == "l1" ) {
			L1_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', 24 ), pose.pdb_info()->pdb2pose( 'L', 34 ) );
			current_start = L1_->start();
			current_end = L1_->stop();
		}
		else if( cdr_name == "l2" ) {
			L2_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', 50 ), pose.pdb_info()->pdb2pose( 'L', 56 ) );
			current_start = L2_->start();
			current_end = L2_->stop();
		}
		else if( cdr_name == "l3" ) {
			L3_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', 89 ), pose.pdb_info()->pdb2pose( 'L', 97 ) );
			current_start = L3_->start();
			current_end = L3_->stop();
		}
	}
	if( cdr_name == "h1" ) {
		H1_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', 26 ), pose.pdb_info()->pdb2pose( 'H', 35 ) );
		current_start = H1_->start();
		current_end = H1_->stop();
	}
	else if( cdr_name == "h2" ) {
		H2_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', 50 ), pose.pdb_info()->pdb2pose( 'H', 65 ) );
		current_start = H2_->start();
		current_end = H2_->stop();
	}
	else if( cdr_name == "h3" ) {
		H3_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', 95 ), pose.pdb_info()->pdb2pose( 'H', 102 ) );
		current_start = H3_->start();
		current_end = H3_->stop();
	}
} // constructor with arguments

void
AntibodyInfo::set_default( bool camelid )
{
	camelid_ = camelid;
	current_start = 0;
	current_end = 0;
	kinked_ = false;
	extended_ = false;
}


void
AntibodyInfo::setup_loops( core::pose::Pose & pose, bool camelid ) {
	TR<<"sequence of the pose is "<< pose.sequence()<<std::endl;
	set_default( camelid );
	TR<<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<<std::endl;
	if( !camelid_ ) {
		L1_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', 24 ), pose.pdb_info()->pdb2pose( 'L', 34 ) );
		L1_->set_cut( L1_->start() + core::Size( ( ( L1_->stop() - L1_->start() ) + 1 ) / 2 ) );
		all_cdr_loops_.add_loop( *L1_ );
		loops_.insert( std::pair<std::string, loops::LoopOP>("l1", L1_) );
		L2_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', 50 ), pose.pdb_info()->pdb2pose( 'L', 56 ) );
		L2_->set_cut( L2_->start() + core::Size( ( ( L2_->stop() - L2_->start() ) + 1 ) / 2 ) );
		all_cdr_loops_.add_loop( *L2_ );
		loops_.insert( std::pair<std::string, loops::LoopOP>("l2", L2_) );
		L3_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', 89 ), pose.pdb_info()->pdb2pose( 'L', 97 ) );
		L3_->set_cut( L3_->start() + core::Size( ( ( L3_->stop() - L3_->start() ) + 1 ) / 2 ) );
		all_cdr_loops_.add_loop( *L3_ );
		loops_.insert( std::pair<std::string, loops::LoopOP>("l3", L3_) );
	}
	TR<<"PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"<<std::endl;
	H1_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', 26 ), pose.pdb_info()->pdb2pose( 'H', 35 ) );
	H1_->set_cut( H1_->start() + core::Size( ( ( H1_->stop() - H1_->start() ) + 1 ) / 2 ) );
	all_cdr_loops_.add_loop( *H1_ );
	loops_.insert( std::pair<std::string, loops::LoopOP>("h1", H1_) );
	H2_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', 50 ), pose.pdb_info()->pdb2pose( 'H', 65 ) );
	H2_->set_cut( H2_->start() + core::Size( ( ( H2_->stop() - H2_->start() ) + 1 ) / 2 ) );
	all_cdr_loops_.add_loop( *H2_ );
	loops_.insert( std::pair<std::string, loops::LoopOP>("h2", H2_) );
	H3_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', 95 ), pose.pdb_info()->pdb2pose( 'H', 102 )+1 );
	H3_->set_cut( H3_->start() + 1 );
	all_cdr_loops_.add_loop( *H3_ );
	loops_.insert( std::pair<std::string, loops::LoopOP>("h3", H3_) );
	TR<<"QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ"<<std::endl;
	hfr_[1][1] = pose.pdb_info()->pdb2pose( 'H', 5 );
	hfr_[1][2] = pose.pdb_info()->pdb2pose( 'H', 6 );
	hfr_[2][1] = pose.pdb_info()->pdb2pose( 'H', 10 );
	hfr_[2][2] = pose.pdb_info()->pdb2pose( 'H', 25 );
	hfr_[3][1] = pose.pdb_info()->pdb2pose( 'H', 36 );
	hfr_[3][2] = pose.pdb_info()->pdb2pose( 'H', 39 );
	hfr_[4][1] = pose.pdb_info()->pdb2pose( 'H', 46 );
	hfr_[4][2] = pose.pdb_info()->pdb2pose( 'H', 49 );
	hfr_[5][1] = pose.pdb_info()->pdb2pose( 'H', 66 );
	hfr_[5][2] = pose.pdb_info()->pdb2pose( 'H', 94 );
	hfr_[6][1] = pose.pdb_info()->pdb2pose( 'H', 103 );
	hfr_[6][2] = pose.pdb_info()->pdb2pose( 'H', 110 );

	all_cdr_loops_.sequential_order();

	all_cdr_fold_tree( pose );

	for( core::Size i = 1; i <= pose.total_residue(); ++i )
		Fv_sequence_.push_back( pose.residue(i).name1() );

	detect_CDR_H3_stem_type( pose );
	TR<<"RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"<<std::endl;
} // set_defaults

loops::LoopOP
AntibodyInfo::get_loop( std::string loop ) {
	LoopMap::iterator iter = loops_.begin();
	iter = loops_.find(loop);
	if ( iter != loops_.end() ) {return iter->second;}
}

void AntibodyInfo::align_to_native( core::pose::Pose & pose, antibody2::AntibodyInfo & native, core::pose::Pose & native_pose ) {

	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, core::id::BOGUS_ATOM_ID );

	for( core::Size j = 1; j <= 6; j++ ) {
		core::Size buffer_for_h3_end(0);
		if( j == 6 ) buffer_for_h3_end = 1;
		for( core::Size res_counter=hfr_[j][1]+buffer_for_h3_end,    nat_counter=native.hfr_[j][1]+buffer_for_h3_end;    res_counter<=hfr_[j][2];       res_counter++, nat_counter++ ) {
			for( core::Size atm_counter=1; atm_counter <= 4; atm_counter++ ) {
				core::id::AtomID const id1( atm_counter, res_counter );
				core::id::AtomID const id2( atm_counter, nat_counter );
				atom_map[ id1 ] = id2;
			}
		}
	}

	core::scoring::superimpose_pose( pose, native_pose, atom_map );

} // align_to_native()

void
AntibodyInfo::detect_CDR_H3_stem_type( core::pose::Pose & pose ) {
	if( camelid_ )
		detect_camelid_CDR_H3_stem_type();
	else
		detect_regular_CDR_H3_stem_type( pose );
	return;
} // detect_CDR_H3_stem_type

void
AntibodyInfo::detect_camelid_CDR_H3_stem_type() {
	TR << "AC Detecting Camelid CDR H3 Stem Type" << std::endl;

	// extract single letter aa codes for the chopped loop residues
	utility::vector1< char > cdr_h3_sequence;
	for( core::Size ii = H3_->start() - 2; ii <= (H3_->stop()); ++ii )
		cdr_h3_sequence.push_back( Fv_sequence_[ii] );

	// Rule for extended
	if( ( ( H3_->stop() - H3_->start() ) ) >= 12 ) {
		if( ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'Y' ) ||
					( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'W' ) ||
					( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'F' ) ) &&
				( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'H' ) &&
				( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'G' ) )
			extended_ = true;
	}

	if( !extended_ ) {
		kinked_ = true;
		if(           ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'R' ) ||
				      ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'Y' ) ||
				((     ( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'Y' ) || ( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'W' )    ) &&
				(     ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'Y' ) || ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'W' )    ) &&
				( 	  ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] != 'Y' ) || ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] != 'W' )    ))
		)
			kinked_ = false;
	}


	TR << "AC Finished Detecting Camelid CDR H3 Stem Type: "
		 << "Kink: " << kinked_ << " Extended: " << extended_ << std::endl;
} // detect_camelid_CDR_H3_stem_type()


void
AntibodyInfo::detect_regular_CDR_H3_stem_type( core::pose::Pose & pose ) {
	TR << "AC Detecting Regular CDR H3 Stem Type" << std::endl;

	bool is_H3( false );

	// extract single letter aa codes for the chopped loop residues
	utility::vector1< char > cdr_h3_sequence;
	for( core::Size ii = H3_->start() - 2; ii <= (H3_->stop() ); ++ii )
		cdr_h3_sequence.push_back( Fv_sequence_[ii] );

	// Rule 1a for standard kink
	if( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'D') {
		kinked_ = true;
		is_H3 = true;
	}

	// Rule 1b for standard extended form
	if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D')
			&& ( (cdr_h3_sequence[2] != 'K') &&
					 (cdr_h3_sequence[2] != 'R') ) && (is_H3 != true)) {
		extended_ = true;
		is_H3 = true;
	}

	if( !is_H3 ) {
		// Rule 1b extension for special kinked form
		bool is_basic( false ); // Special basic residue exception flag
		for(core::Size ii = 3; ii <= core::Size(cdr_h3_sequence.size() - 4);
				ii++) {
			if( cdr_h3_sequence[ii] == 'R' || cdr_h3_sequence[ii] == 'K') {
				is_basic = true;
				break;
			}
		}

		if( !is_basic ) {
			core::Size L49_pose_number = pose.pdb_info()->pdb2pose( 'L', 49 );
			char aa_code_L49 = pose.residue( L49_pose_number ).name1();
			if( aa_code_L49 == 'R' || aa_code_L49 == 'K')
				is_basic = true;
		}
		if( is_basic ) {
			kinked_ = true;
			is_H3 = true;
		}
	}

	// Rule 1c for kinked form with salt bridge
	if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D') &&
			( (cdr_h3_sequence[2] == 'K') ||
				(cdr_h3_sequence[2] == 'R') ) &&
			( (cdr_h3_sequence[1] != 'K') &&
				(cdr_h3_sequence[1] != 'R') ) && (is_H3 != true) ) {
		kinked_ = true;
		is_H3 = true;
		if( !is_H3 ) {
			bool is_basic( false ); // Special basic residue exception flag
			core::Size L46_pose_number = pose.pdb_info()->pdb2pose( 'L', 46 );
			char aa_code_L46 = pose.residue( L46_pose_number ).name1();
			if( aa_code_L46 == 'R' || aa_code_L46 == 'K')
				is_basic = true;
			if( is_basic ) {
				extended_ = true;
				is_H3 = true;
			}
		}
	}

	// Rule 1d for extened form with salt bridge
	if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D') &&
			( ( cdr_h3_sequence[ 2 ] == 'K') ||
				(cdr_h3_sequence[2] == 'R')) &&
			( (cdr_h3_sequence[1] == 'K') ||
				(cdr_h3_sequence[1] == 'R') ) && (is_H3 != true) ) {
		extended_ = true;
		is_H3 = true;
	}

	TR << "AC Finished Detecting Regular CDR H3 Stem Type: "
		 << "Kink: " << kinked_ << " Extended: " << extended_ << std::endl;
} // detect_regular_CDR_H3_stem_type()

void
AntibodyInfo::all_cdr_fold_tree( core::pose::Pose & pose ) {
	using namespace core::kinematics;

	all_cdr_loops_.sequential_order();

	FoldTree f;
	f.clear();

	core::Size jump_num = 0;
	for( loops::Loops::const_iterator it=all_cdr_loops_.begin(), it_end=all_cdr_loops_.end(), it_next; it < it_end; ++it ) {

		it_next = it;
		it_next++;

		if( it == all_cdr_loops_.begin() ) f.add_edge( 1, it->start()-1, Edge::PEPTIDE );

		jump_num++;
		f.add_edge( it->start()-1, it->stop()+1, jump_num );
		f.add_edge( it->start()-1, it->cut(),  Edge::PEPTIDE );
		f.add_edge( it->cut()+1, it->stop()+1, Edge::PEPTIDE );
		if( it == (it_end-1) )
			f.add_edge( it->stop()+1, pose.total_residue(), Edge::PEPTIDE);
		else
			f.add_edge( it->stop()+1, it_next->start()-1, Edge::PEPTIDE );
	}

	f.reorder(1);
	pose.fold_tree( f );

} // all_cdr_fold_tree()


} // namespace antibody2
} // namespace protocols

