// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody2/AntibodyInfo.cc
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)

#include <protocols/antibody2/AntibodyInfo.hh>

// Rosetta Headers
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <iostream>
#include <fstream>
#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <utility/exit.hh>
#include <protocols/antibody2/AntibodyUtil.hh>
#include <utility/tools/make_vector1.hh>




static basic::Tracer TR("antibody2.AntibodyInfo");

namespace protocols{
namespace antibody2{

/// default constructor
AntibodyInfo::AntibodyInfo() {
	set_default( false/*is_camelid*/ );

	for( core::Size i = 0; i <= 6; i++ ) hfr_[i][0] = hfr_[i][1] = hfr_[i][2] = 0;
}


/// constructor with arguments
AntibodyInfo::AntibodyInfo( core::pose::Pose & pose ) {
	setup_CDR_loops( pose, false/*is_camelid*/ );
}

/// constructor with arguments
AntibodyInfo::AntibodyInfo( core::pose::Pose & pose, bool is_camelid ) {
	setup_CDR_loops( pose, is_camelid );
}



/// constructor with one cdr loop arguments
AntibodyInfo::AntibodyInfo( core::pose::Pose & pose, std::string cdr_name )
{
    if(is_camelid_){
        if (cdr_name == "l1" || cdr_name == "l2" ||cdr_name == "l3") {
            utility_exit_with_message("This is Camelid antibody, No Light Chain !!!");
        }
    }

	set_default( false/*is_camelid*/ );

	if( !is_camelid_ ) {
		if( cdr_name == "l1" ) {
			L1_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', CDR_numbering_begin_["l1"] ), 
                                   pose.pdb_info()->pdb2pose( 'L', CDR_numbering_end_["l1"] ) );
			current_start = L1_->start();
			current_end = L1_->stop();
            
		}
		else if( cdr_name == "l2" ) {
			L2_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', CDR_numbering_begin_["l2"] ), 
                                   pose.pdb_info()->pdb2pose( 'L', CDR_numbering_end_["l2"] ) );
			current_start = L2_->start();
			current_end = L2_->stop();
		}
		else if( cdr_name == "l3" ) {
			L3_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', CDR_numbering_begin_["l3"] ), 
                                   pose.pdb_info()->pdb2pose( 'L', CDR_numbering_end_["l3"] ) );
			current_start = L3_->start();
			current_end = L3_->stop();
		}
	}
	if( cdr_name == "h1" ) {
		H1_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', CDR_numbering_begin_["h1"] ), 
                               pose.pdb_info()->pdb2pose( 'H', CDR_numbering_end_["h1"] ) );
		current_start = H1_->start();
		current_end = H1_->stop();
	}
	else if( cdr_name == "h2" ) {
		H2_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', CDR_numbering_begin_["h2"] ), 
                               pose.pdb_info()->pdb2pose( 'H', CDR_numbering_end_["h2"] ) );
		current_start = H2_->start();
		current_end = H2_->stop();
	}
	else if( cdr_name == "h3" ) {
		H3_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', CDR_numbering_begin_["h3"] ), 
                               pose.pdb_info()->pdb2pose( 'H', CDR_numbering_end_["h3"] ) );
		current_start = H3_->start();
		current_end = H3_->stop();
	}


} // constructor with arguments








void
AntibodyInfo::set_default( bool is_camelid )
{
	is_camelid_ = is_camelid;
	current_start = 0;
	current_end = 0;
	kinked_H3_ = false;
	extended_H3_ = false;
    get_CDRs_numbering();
    LH_dock_jumps_ .push_back(1);
    LH_dock_partners_ = "L_H";

}





// TODO:
// JQX:
// The code assumed that the input PDB has been been renumbered using the Chothia
// numbering scheme [see the "get_CDRs_numbering()" function below]: as a matter 
// of fact, since this code is desigend for the Rosetta Antibody Homology Modeling
// the input is always the structure made from different templates, and they are 
// always renumbered by the perl script from the Rosetta Antibody Server
// A smart way would be to use the "identify_CDR_from_a_sequence()" to automatically
// check this out. On my list!
void AntibodyInfo::setup_CDR_loops( core::pose::Pose & pose, bool is_camelid ) {

	set_default( is_camelid );

	if( !is_camelid_ ) {
		L1_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', CDR_numbering_begin_["l1"] ), 
                               pose.pdb_info()->pdb2pose( 'L', CDR_numbering_end_["l1"] ) );
		L1_->set_cut( L1_->start() + core::Size( ( ( L1_->stop() - L1_->start() ) + 1 ) / 2 ) );
		all_cdr_loops_.add_loop( *L1_ );
		loops_.insert( std::pair<std::string, loops::LoopOP>("l1", L1_) );
        L1_seq_ = get_seq_from_a_loop(pose, L1_);
        
		L2_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', CDR_numbering_begin_["l2"] ), 
                               pose.pdb_info()->pdb2pose( 'L', CDR_numbering_end_["l2"] ) );
		L2_->set_cut( L2_->start() + core::Size( ( ( L2_->stop() - L2_->start() ) + 1 ) / 2 ) );
		all_cdr_loops_.add_loop( *L2_ );
		loops_.insert( std::pair<std::string, loops::LoopOP>("l2", L2_) );
        L2_seq_ = get_seq_from_a_loop(pose, L2_);
        
		L3_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'L', CDR_numbering_begin_["l3"] ), 
                               pose.pdb_info()->pdb2pose( 'L', CDR_numbering_end_["l3"] ) );
		L3_->set_cut( L3_->start() + core::Size( ( ( L3_->stop() - L3_->start() ) + 1 ) / 2 ) );
		all_cdr_loops_.add_loop( *L3_ );
		loops_.insert( std::pair<std::string, loops::LoopOP>("l3", L3_) );
        L3_seq_ = get_seq_from_a_loop(pose, L3_);
	}

	H1_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', CDR_numbering_begin_["h1"] ), 
                           pose.pdb_info()->pdb2pose( 'H', CDR_numbering_end_["h1"] ) );
	H1_->set_cut( H1_->start() + core::Size( ( ( H1_->stop() - H1_->start() ) + 1 ) / 2 ) );
	all_cdr_loops_.add_loop( *H1_ );
	loops_.insert( std::pair<std::string, loops::LoopOP>("h1", H1_) );
    H1_seq_ = get_seq_from_a_loop(pose, H1_);
    
	H2_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', CDR_numbering_begin_["h2"] ), 
                           pose.pdb_info()->pdb2pose( 'H', CDR_numbering_end_["h2"] ) );
	H2_->set_cut( H2_->start() + core::Size( ( ( H2_->stop() - H2_->start() ) + 1 ) / 2 ) );
	all_cdr_loops_.add_loop( *H2_ );
	loops_.insert( std::pair<std::string, loops::LoopOP>("h2", H2_) );
    H2_seq_ = get_seq_from_a_loop(pose, H2_);

    
    // JQX:
    // One should always see 95-102 as the positions for your H3 in your FR02.pdb, but as a matter of fact,
    // the antibody script just copied h3.pdb (heavy atoms) into the FR02.pdb, sometimes one sees the stop
    // postition pdb number 98, not 102, if the h3.pdb is short. Therefore, one useing the pdb number 102 to 
    // define h3 fails! 
    // But in FR02.pdb, you always see 103, because 103 is on the framework. The idea is to find the pose number
    // of PDB number 103, then minus 1 will give you the last residue of h3.
    Size pose_num_end_plus_one = pose.pdb_info()->pdb2pose( 'H', CDR_numbering_end_["h3"]+1 ) ;  // PDB number 103 to get pose number
    Size pose_num_end = pose_num_end_plus_one - 1; // pose number -1
    H3_ = new loops::Loop( pose.pdb_info()->pdb2pose( 'H', CDR_numbering_begin_["h3"] ), 
                           pose_num_end );
    H3_->set_cut( H3_->start() + 1 );  // why this is different compared to other cuts of other loops? 
                                       // Aroop seems did this in his old R3 code, CHECK LATER !!!
	all_cdr_loops_.add_loop( *H3_ );
	loops_.insert( std::pair<std::string, loops::LoopOP>("h3", H3_) );
    H3_seq_ = get_seq_from_a_loop(pose, H3_);

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

	for( core::Size i = 1; i <= pose.total_residue(); ++i ){
		Fv_sequence_.push_back( pose.residue(i).name1() );
    }

	detect_and_set_CDR_H3_stem_type( pose );

} // set_defaults








loops::LoopOP AntibodyInfo::get_CDR_loop( std::string loop ) {
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







void AntibodyInfo::detect_and_set_CDR_H3_stem_type( core::pose::Pose & pose ) {
	if( is_camelid_ )
		detect_and_set_camelid_CDR_H3_stem_type();
	else
		detect_and_set_regular_CDR_H3_stem_type( pose );
	return;
} // detect_CDR_H3_stem_type




void AntibodyInfo::detect_and_set_camelid_CDR_H3_stem_type() {
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
			extended_H3_ = true;
	}

	if( !extended_H3_ ) {
		kinked_H3_ = true;
		if(           ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'R' ) ||
				      ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'Y' ) ||
				((     ( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'Y' ) || ( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'W' )    ) &&
				(     ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'Y' ) || ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'W' )    ) &&
				( 	  ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] != 'Y' ) || ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] != 'W' )    ))
		)
			kinked_H3_ = false;
	}


	TR << "AC Finished Detecting Camelid CDR H3 Stem Type: "
		 << "Kink: " << kinked_H3_ << " Extended: " << extended_H3_ << std::endl;
} 




void AntibodyInfo::detect_and_set_regular_CDR_H3_stem_type( core::pose::Pose & pose ) {
	TR << "AC Detecting Regular CDR H3 Stem Type" << std::endl;

	bool is_H3( false );

	// extract single letter aa codes for the chopped loop residues
	utility::vector1< char > cdr_h3_sequence;
	for( core::Size ii = H3_->start() - 2; ii <= (H3_->stop() ); ++ii )
		cdr_h3_sequence.push_back( Fv_sequence_[ii] );

	// Rule 1a for standard kink
	if( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'D') {
		kinked_H3_ = true;
		is_H3 = true;
	}

	// Rule 1b for standard extended form
	if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D')
			&& ( (cdr_h3_sequence[2] != 'K') &&
					 (cdr_h3_sequence[2] != 'R') ) && (is_H3 != true)) {
		extended_H3_ = true;
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
			kinked_H3_ = true;
			is_H3 = true;
		}
	}

	// Rule 1c for kinked form with salt bridge
	if( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D') &&
			( (cdr_h3_sequence[2] == 'K') ||
				(cdr_h3_sequence[2] == 'R') ) &&
			( (cdr_h3_sequence[1] != 'K') &&
				(cdr_h3_sequence[1] != 'R') ) && (is_H3 != true) ) {
		kinked_H3_ = true;
		is_H3 = true;
		if( !is_H3 ) {
			bool is_basic( false ); // Special basic residue exception flag
			core::Size L46_pose_number = pose.pdb_info()->pdb2pose( 'L', 46 );
			char aa_code_L46 = pose.residue( L46_pose_number ).name1();
			if( aa_code_L46 == 'R' || aa_code_L46 == 'K')
				is_basic = true;
			if( is_basic ) {
				extended_H3_ = true;
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
		extended_H3_ = true;
		is_H3 = true;
	}

	TR << "AC Finished Detecting Regular CDR H3 Stem Type: "
		 << "Kink: " << kinked_H3_ << " Extended: " << extended_H3_ << std::endl;
} // detect_regular_CDR_H3_stem_type()










void AntibodyInfo::all_cdr_fold_tree( core::pose::Pose & pose ) {
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


    
    
    
    
    // JQX:: assuming Chothia numbering
    //   setup_CDRs_numbering
void AntibodyInfo::get_CDRs_numbering(){
    CDR_numbering_begin_.insert( std::pair<std::string, core::Size>("l1", 24) );
    CDR_numbering_end_.insert( std::pair<std::string, core::Size>("l1", 34) );
        
    CDR_numbering_begin_.insert( std::pair<std::string, core::Size>("l2", 50) );
    CDR_numbering_end_.insert( std::pair<std::string, core::Size>("l2", 56) );
            
    CDR_numbering_begin_.insert( std::pair<std::string, core::Size>("l3", 89) );
    CDR_numbering_end_.insert( std::pair<std::string, core::Size>("l3", 97) );
            
    CDR_numbering_begin_.insert( std::pair<std::string, core::Size>("h1", 26) );
    CDR_numbering_end_.insert( std::pair<std::string, core::Size>("h1", 35) );
            
    CDR_numbering_begin_.insert( std::pair<std::string, core::Size>("h2", 50) );
    CDR_numbering_end_.insert( std::pair<std::string, core::Size>("h2", 65) );
            
    CDR_numbering_begin_.insert( std::pair<std::string, core::Size>("h3", 95) );
    CDR_numbering_end_.insert( std::pair<std::string, core::Size>("h3", 102) );
}

    
    
    
    
    
    



/// @details  Show the complete setup of the docking protocol
void AntibodyInfo::show( std::ostream & out ) {
    //      if ( !flags_and_objects_are_in_sync_ ){
    //              sync_objects_with_flags();
    //      }
    out << *this;
}

std::ostream & operator<<(std::ostream& out, const AntibodyInfo & ab_info )
{
        using namespace ObjexxFCL::fmt;
        // All output will be 80 characters - 80 is a nice number, don't you think?
        std::string line_marker = "///";
        out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
        out << line_marker << A( 47, "Rosetta Antibody Info" ) << space( 27 ) << line_marker << std::endl;
        out << line_marker << space( 74 ) << line_marker << std::endl;
        // Display the movable jumps that will be used in docking

        out << line_marker << " L1 info: "<<std::endl;
        out << line_marker << "           length:  "<<ab_info.L1_seq_.size() <<std::endl;
        out << line_marker << "         sequence:  "<<ab_info.L1_seq_ <<std::endl;
        out << line_marker << "        loop_info:  "<<*ab_info.L1_<<std::endl;

        out << line_marker << " L2 info: "<<std::endl;
        out << line_marker << "           length:  "<<ab_info.L2_seq_.size()<<std::endl;
        out << line_marker << "         sequence:  "<<ab_info.L2_seq_ <<std::endl;
        out << line_marker << "        loop_info:  "<<*ab_info.L2_<<std::endl;

        out << line_marker << " L3 info: "<<std::endl;
        out << line_marker << "           length:  "<<ab_info.L3_seq_.size()<<std::endl;
        out << line_marker << "         sequence:  "<<ab_info.L3_seq_<<std::endl;
        out << line_marker << "        loop_info:  "<<*ab_info.L3_<<std::endl;

        out << line_marker << " H1 info: "<<std::endl;
        out << line_marker << "           length:  "<<ab_info.H1_seq_.size()<<std::endl;
        out << line_marker << "         sequence:  "<<ab_info.H1_seq_<<std::endl;
        out << line_marker << "        loop_info:  "<<*ab_info.H1_<<std::endl;

        out << line_marker << " H2 info: "<<std::endl;
        out << line_marker << "           length:  "<<ab_info.H2_seq_.size()<<std::endl;
        out << line_marker << "         sequence:  "<<ab_info.H2_seq_<<std::endl;
        out << line_marker << "        loop_info:  "<<*ab_info.H2_<<std::endl;

        out << line_marker << " H3 info: "<<std::endl;
        out << line_marker << "           length:  "<<ab_info.H3_seq_.size()<<std::endl;
        out << line_marker << "         sequence:  "<<ab_info.H3_seq_<<std::endl;
        out << line_marker << "        loop_info:  "<<*ab_info.H3_<<std::endl;

        // Close the box I have drawn
        out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
        return out;
}


    
    
//TODO:
//JQX: make Daisuke's code compatible with my code
//
// Identify 3 CDRs from a sequence
// Automatically judge heavy or light chain (I hope!)
// Author: Daisuke Kuroda
// Last modified 03/16/2012 by Daisuke Kuroda
//
    
void AntibodyInfo::identify_CDR_from_a_sequence(std::string & querychain)
{
        using namespace std;

        int l1found = 0, l2found = 0, l3found = 1, h1found = 1, h2found = 0, h3found = 1; 
        // 0 if exst; otherwise 1.
        int lenl1 =0, lenl2 = 0, lenl3 = 0, lenh1 = 0, lenh2 = 0, lenh3 = 0;
        int posl1_s = 0, posl1_e = 0, posl2_s = 0, posl2_e = 0, posl3_s = 0, posl3_e = 0;
        int posh1_s = 0, posh1_e = 0, posh2_s = 0, posh2_e = 0, posh3_s = 0, posh3_e = 0;
        int i = 0, k = 0, l = 0, m = 0;
        
        string seql1, seql2, seql3, seqh1, seqh2, seqh3;
        string frl3, frh1, frh3;
        
        // For L3
        string p1_l3[] = {"G","A","P","C","D","E","Q","N","R","K","H","W","Y","F","M","T","V","I","S","L"};
        
        // For H1
        string p1_h1[] = {"I","V","F","Y","A","M","L","N"};
        string p2_h1[] = {"R","K","Q","V","N","C"};
        string p3_h1[] = {"Q","K","H","E","L","R"};
        
        // For H3
        string p1_h3[] = {"G","A","P","C","D","E","Q","N","R","K","H","W","Y","F","M","T","V","I","S","L"};
        string p2_h3[] = {"G","R","D"};
        
        
        cout << "*** Query sequence ***" << endl;
        cout << querychain << endl;
        cout << endl;
        
        /*****************************************************/
        /***************** Is it light chain? ****************/
        /*****************************************************/
        /* L1 search Start */
        if(querychain.find("WYL") != string::npos){
            posl1_e = querychain.find("WYL") - 1;
        }else if(querychain.find("WLQ") != string::npos){
            posl1_e = querychain.find("WLQ") - 1;
        }else if(querychain.find("WFQ") != string::npos){
            posl1_e = querychain.find("WFQ") - 1;
        }else if(querychain.find("WYQ") != string::npos){
            posl1_e = querychain.find("WYQ") - 1;
        }else if(querychain.find("WYH") != string::npos){
            posl1_e = querychain.find("WYH") - 1;
        }else if(querychain.find("WVQ") != string::npos){
            posl1_e = querychain.find("WVQ") - 1;
        }else if(querychain.find("WVR") != string::npos){
            posl1_e = querychain.find("WVR") - 1;
        }else if(querychain.find("WWQ") != string::npos){
            posl1_e = querychain.find("WWQ") - 1;
        }else if(querychain.find("WVK") != string::npos){
            posl1_e = querychain.find("WVK") - 1;
        }else if(querychain.find("WLL") != string::npos){
            posl1_e = querychain.find("WLL") - 1;
        }else if(querychain.find("WFL") != string::npos){
            posl1_e = querychain.find("WFL") - 1;
        }else if(querychain.find("WVF") != string::npos){
            posl1_e = querychain.find("WVF") - 1;
        }else if(querychain.find("WIQ") != string::npos){
            posl1_e = querychain.find("WIQ") - 1;
        }else if(querychain.find("WYR") != string::npos){
            posl1_e = querychain.find("WYR") - 1;
        }else if(querychain.find("WNQ") != string::npos){
            posl1_e = querychain.find("WNQ") - 1;
        }else if(querychain.find("WHL") != string::npos){
            posl1_e = querychain.find("WHL") - 1;
        }else{
            l1found = 1;
        }
        
        if(l1found != 1){
            posl1_s = querychain.find("C") + 1;
            lenl1   = posl1_e - posl1_s + 1;
            seql1   = querychain.substr(posl1_s,lenl1);
        }
        /* L1 search Finish */
        
        /* L2 search start */
        if(l1found != 1){
            posl2_s = posl1_e + 16;
            posl2_e = posl2_s + 6;
            lenl2   = posl2_e - posl2_s + 1;
            seql2   = querychain.substr(posl2_s,lenl2);
        }else{
            l2found = 1;
        }
        /* L2 search end */
        
        /* L3 search Start */
        for(k = 0;k < 20; k++){
            for(l = 0;l < 3; l++){
                frl3 = "FG" + p1_l3[k] + "G";
                
                if(querychain.find(frl3, 80) != string::npos){
                    posl3_e = querychain.find(frl3,80) - 1;
                    posl3_s = querychain.find("C",70) + 1;
                    lenl3   = posl3_e - posl3_s + 1;
                    seql3   = querychain.substr(posl3_s,lenl3);
                    l3found = 0;
                }
            }
        }
        /* L3 search Finish */
        
        /*****************************************************/
        /***************** Is it heavy chain? ****************/
        /*****************************************************/
        if(l1found == 1 || l2found == 1 || l3found == 1){
            /* H1 search Start */
            for(k = 0;k < 8; k++){
                for(l = 0; l < 6; l++){
                    for(m = 0; m < 6; m++){
                        frh1 = "W" + p1_h1[k] + p2_h1[l] + p3_h1[m];
                        
                        if(querychain.find(frh1, 0) != string::npos){
                            posh1_e = querychain.find(frh1, 0) - 1;
                            h1found = 0;
                        }
                    }
                }
            }
            
            if(h1found != 1){
                posh1_s = querychain.find("C") + 4;
                lenh1   = posh1_e - posh1_s + 1;
                seqh1   = querychain.substr(posh1_s, lenh1);
            }
            /* H1 search Finish */
            
            /* H3 search Start */
            for(k = 0;k < 20; k++){
                for(l = 0;l < 3; l++){
                    frh3 = "WG" + p1_h3[k] + p2_h3[l];
                    
                    if(querychain.find(frh3, 80) != string::npos){
                        posh3_e = querychain.find(frh3,80) - 1;
                        h3found = 0;
                    }
                }
            }
            
            if(querychain.find("C", 80) != string::npos){
                posh3_s = querychain.find("C", 80) + 3;
            }else{
                h3found = 1;
            }
            
            if(h3found != 1){
                lenh3 = posh3_e - posh3_s + 1;
                seqh3 = querychain.substr(posh3_s,lenh3);
            }
            /* H3 search Finish */
            
            /* H2 search start */
            if(h1found != 1 && h3found != 1){
                posh2_s = posh1_e + 15;
                posh2_e = posh3_s - 33;
                lenh2   = posh2_e - posh2_s + 1;
                seqh2   = querychain.substr(posh2_s,lenh2);
            }
            /* H2 search end */
        }
        
        if(l1found == 0 && l2found == 0 && l3found == 0){
            cout << "This is Light chain!" << endl;
            cout << "L1 info:" << "\t" << lenl1 << "\t" << posl1_s << "-" << posl1_e << "\t" << seql1 << endl;
            cout << "L2 info:" << "\t" << lenl2 << "\t" << posl2_s << "-" << posl2_e << "\t" << seql2 << endl;
            cout << "L3 info:" << "\t" << lenl3 << "\t" << posl3_s << "-" << posl3_e << "\t" << seql3 << endl;
        }else if(h1found == 0 && h2found == 0 && h3found == 0){
            cout << "This is Heavy chain!" << endl;
            cout << "H1 info:" << "\t" << lenh1 << "\t" << posh1_s << "-" << posh1_e << "\t" << seqh1 << endl;
            cout << "H2 info:" << "\t" << lenh2 << "\t" << posh2_s << "-" << posh2_e << "\t" << seqh2 << endl;
            cout << "H3 info:" << "\t" << lenh3 << "\t" << posh3_s << "-" << posh3_e << "\t" << seqh3 << endl;
        }else{
            cout << "Some CDRs seem to be missing!" << endl;
        }
        
 
}

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
void AntibodyInfo::load_CDR_query_info_to_check(){

        using namespace std;
        ifstream inf;
        std::string temp;

        inf.open("input/query.l1");
        if(!inf.is_open()) {utility_exit_with_message("Cannot open 'query.l1' file!!");}
        inf>>temp; inf>>L1_seq_;
        inf.close();

        inf.open("input/query.l2");
        if(!inf.is_open()) {utility_exit_with_message("Cannot open 'query.l2' file!!");}
        inf>>temp; inf>>L2_seq_;
        inf.close();

        inf.open("input/query.l3");
        if(!inf.is_open()) {utility_exit_with_message("Cannot open 'query.l3' file!!");}
        inf>>temp; inf>>L3_seq_;
        inf.close();

        inf.open("input/query.h1");
        if(!inf.is_open()) {utility_exit_with_message("Cannot open 'query.h1' file!!");}
        inf>>temp; inf>>H1_seq_;
        inf.close();

        inf.open("input/query.h2");
        if(!inf.is_open()) {utility_exit_with_message("Cannot open 'query.h2' file!!");}
        inf>>temp; inf>>H2_seq_;
        inf.close();

        inf.open("input/query.h3");
        if(!inf.is_open()) {utility_exit_with_message("Cannot open 'query.h3' file!!");}
        inf>>temp; inf>>H3_seq_;
        inf.close();
}




} // namespace antibody2
} // namespace protocols


