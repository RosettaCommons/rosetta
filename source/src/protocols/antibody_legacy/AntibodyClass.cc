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
/// @author Aroop Sircar

// Rosetta Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <protocols/antibody_legacy/AntibodyClass.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>


static THREAD_LOCAL basic::Tracer TR( "antibody" );

namespace protocols {
namespace antibody_legacy {

/// default constructor
Antibody::Antibody() {
	camelid_ = false;
	current_start = 0;
	current_end = 0;
	kinked_ = false;
	extended_ = false;

	for ( core::Size i = 1; i <= 3; i++ ) {
		cdrl_[i][1] = cdrl_[i][2] = cdrh_[i][1] = cdrh_[i][2] = 0;
	}
	for ( core::Size i = 1; i <= 7; i++ ) {
		lfr_[i][1] = lfr_[i][2] = 0;
	}
	for ( core::Size i = 1; i <= 6; i++ ) {
		hfr_[i][1] = hfr_[i][2] = 0;
	}
}


/// constructor with arguments
Antibody::Antibody( core::pose::Pose& pose_in ) {
	camelid_ = false;
	Fv = pose_in;
	set_defaults();
	current_start = 0;
	current_end = 0;

	update_sequence();
	kinked_ = false;
	extended_ = false;
	detect_CDR_H3_stem_type();
}

/// constructor with arguments
Antibody::Antibody( core::pose::Pose& pose_in, bool camelid ) {
	camelid_ = camelid;
	Fv = pose_in;
	set_defaults();
	current_start = 0;
	current_end = 0;

	update_sequence();
	kinked_ = false;
	extended_ = false;
	detect_CDR_H3_stem_type();
}

/// constructor with arguments
Antibody::Antibody(
	core::pose::Pose& pose_in,
	std::string cdr_name ) {
	camelid_ = false;
	Fv = pose_in;

	if ( !camelid_ ) {
		if ( cdr_name == "l1" ) {
			cdrl_[1][1] = Fv.pdb_info()->pdb2pose( 'L', 24 );
			cdrl_[1][2] = Fv.pdb_info()->pdb2pose( 'L', 34 );
			current_start = cdrl_[1][1];
			current_end = cdrl_[1][2];
		} else if ( cdr_name == "l2" ) {
			cdrl_[2][1] = Fv.pdb_info()->pdb2pose( 'L', 50 );
			cdrl_[2][2] = Fv.pdb_info()->pdb2pose( 'L', 56 );
			current_start = cdrl_[2][1];
			current_end = cdrl_[2][2];
		} else if ( cdr_name == "l3" ) {
			cdrl_[3][1] = Fv.pdb_info()->pdb2pose( 'L', 89 );
			cdrl_[3][2] = Fv.pdb_info()->pdb2pose( 'L', 97 );
			current_start = cdrl_[3][1];
			current_end = cdrl_[3][2];
		}
	}
	if ( cdr_name == "h1" ) {
		cdrh_[1][1] = Fv.pdb_info()->pdb2pose( 'H', 26 );
		cdrh_[1][2] = Fv.pdb_info()->pdb2pose( 'H', 35 );
		current_start = cdrh_[1][1];
		current_end = cdrh_[1][2];
	} else if ( cdr_name == "h2" ) {
		cdrh_[2][1] = Fv.pdb_info()->pdb2pose( 'H', 50 );
		cdrh_[2][2] = Fv.pdb_info()->pdb2pose( 'H', 65 );
		current_start = cdrh_[2][1];
		current_end = cdrh_[2][2];
	} else if ( cdr_name == "h3" ) {
		cdrh_[3][1] = Fv.pdb_info()->pdb2pose( 'H', 95 );
		cdrh_[3][2] = Fv.pdb_info()->pdb2pose( 'H', 102 );
		current_start = cdrh_[3][1];
		current_end = cdrh_[3][2];
	} else {
		current_start = 0;
		current_end = 0;
	}
} // constructor with arguments


void
Antibody::set_defaults() {
	if ( !camelid_ ) {
		lfr_[1][1] = Fv.pdb_info()->pdb2pose( 'L', 4 );
		lfr_[1][2] = Fv.pdb_info()->pdb2pose( 'L', 6 );
		lfr_[2][1] = Fv.pdb_info()->pdb2pose( 'L', 10 );
		lfr_[2][2] = Fv.pdb_info()->pdb2pose( 'L', 23 );
		cdrl_[1][1] = Fv.pdb_info()->pdb2pose( 'L', 24 );
		cdrl_[1][2] = Fv.pdb_info()->pdb2pose( 'L', 34 );
		lfr_[3][1] = Fv.pdb_info()->pdb2pose( 'L', 35 );
		lfr_[3][2] = Fv.pdb_info()->pdb2pose( 'L', 38 );
		lfr_[4][1] = Fv.pdb_info()->pdb2pose( 'L', 45 );
		lfr_[4][2] = Fv.pdb_info()->pdb2pose( 'L', 49 );
		cdrl_[2][1] = Fv.pdb_info()->pdb2pose( 'L', 50 );
		cdrl_[2][2] = Fv.pdb_info()->pdb2pose( 'L', 56 );
		lfr_[5][1] = Fv.pdb_info()->pdb2pose( 'L', 57 );
		lfr_[5][2] = Fv.pdb_info()->pdb2pose( 'L', 66 );
		lfr_[6][1] = Fv.pdb_info()->pdb2pose( 'L', 71 );
		lfr_[6][2] = Fv.pdb_info()->pdb2pose( 'L', 88 );
		cdrl_[3][1] = Fv.pdb_info()->pdb2pose( 'L', 89 );
		cdrl_[3][2] = Fv.pdb_info()->pdb2pose( 'L', 97 );
		lfr_[7][1] = Fv.pdb_info()->pdb2pose( 'L', 98 );
		lfr_[7][2] = Fv.pdb_info()->pdb2pose( 'L', 104 );
	}
	hfr_[1][1] = Fv.pdb_info()->pdb2pose( 'H', 5 );
	hfr_[1][2] = Fv.pdb_info()->pdb2pose( 'H', 6 );
	hfr_[2][1] = Fv.pdb_info()->pdb2pose( 'H', 10 );
	hfr_[2][2] = Fv.pdb_info()->pdb2pose( 'H', 25 );
	cdrh_[1][1] = Fv.pdb_info()->pdb2pose( 'H', 26 );
	cdrh_[1][2] = Fv.pdb_info()->pdb2pose( 'H', 35 );
	hfr_[3][1] = Fv.pdb_info()->pdb2pose( 'H', 36 );
	hfr_[3][2] = Fv.pdb_info()->pdb2pose( 'H', 39 );
	hfr_[4][1] = Fv.pdb_info()->pdb2pose( 'H', 46 );
	hfr_[4][2] = Fv.pdb_info()->pdb2pose( 'H', 49 );
	cdrh_[2][1] = Fv.pdb_info()->pdb2pose( 'H', 50 );
	cdrh_[2][2] = Fv.pdb_info()->pdb2pose( 'H', 65 );
	hfr_[5][1] = Fv.pdb_info()->pdb2pose( 'H', 66 );
	hfr_[5][2] = Fv.pdb_info()->pdb2pose( 'H', 94 );
	cdrh_[3][1] = Fv.pdb_info()->pdb2pose( 'H', 95 );
	cdrh_[3][2] = Fv.pdb_info()->pdb2pose( 'H', 102 );
	hfr_[6][1] = Fv.pdb_info()->pdb2pose( 'H', 103 );
	hfr_[6][2] = Fv.pdb_info()->pdb2pose( 'H', 110 );


	cdr_h3_cut_ = cdrh_[3][1] + 1;

	populate_all_cdrs();
	all_cdr_fold_tree();
} // set_defaults

void
Antibody::set_Fv( core::pose::Pose& pose_in ) {
	Fv = pose_in;
	camelid_ = false;
	set_defaults();
	current_start = 0;
	current_end = 0;

	update_sequence();
	kinked_ = false;
	extended_ = false;
	detect_CDR_H3_stem_type();
} // set_Fv

void
Antibody::set_Fv( core::pose::Pose& pose_in, bool camelid ) {
	Fv = pose_in;
	camelid_ = camelid;
	set_defaults();
	current_start = 0;
	current_end = 0;

	update_sequence();
	kinked_ = false;
	extended_ = false;
	detect_CDR_H3_stem_type();
} // set_Fv

void
Antibody::populate_all_cdrs() {

	core::Size begin(0), end(0), size(0), cut(0);

	std::string cdr_name[6] = { "l1", "l2", "l3", "h1", "h2", "h3" };

	core::Size cdr_total( 6 );

	for ( core::Size i = 0; i < cdr_total; i++ ) {
		if ( cdr_name[i] == "l1" ) {
			begin = cdrl_[1][1];
			end = cdrl_[1][2];
		} else if ( cdr_name[i] == "l2" ) {
			begin = cdrl_[2][1];
			end = cdrl_[2][2];
		} else if ( cdr_name[i] == "l3" ) {
			begin = cdrl_[3][1];
			end = cdrl_[3][2];
		} else if ( cdr_name[i] == "h1" ) {
			begin = cdrh_[1][1];
			end = cdrh_[1][2];
		} else if ( cdr_name[i] == "h2" ) {
			begin = cdrh_[2][1];
			end = cdrh_[2][2];
		} else if ( cdr_name[i] == "h3" ) {
			begin = cdrh_[3][1];
			end = cdrh_[3][2] + 1; // for the extra stem residue
		}

		//if( flank_relax && cdr_name[i] == "h3" ) {
		// begin = begin - h3_flank;
		// end = end + h3_flank;
		//}
		size = ( end - begin ) + 1;
		cut = begin + core::Size( size / 2 );
		if ( cdr_name[i] == "h3" ) {
			cut = cdr_h3_cut_;
		}
		if ( !camelid_ ) {
			all_cdr_loops.add_loop( begin, end, cut, 0, false);
		} else {
			if ( (cdr_name[i] == "h1") || (cdr_name[i] == "h2") ||
					(cdr_name[i] == "h3") ) {
				all_cdr_loops.add_loop( begin, end, cut, 0, false);
			}
		}
	}
	all_cdr_loops.sequential_order();
} // populate_all_cdrs

void
Antibody::update_sequence() {
	for ( core::Size i = 1; i <= Fv.size(); ++i ) {
		Fv_sequence_.push_back( Fv.residue(i).name1() );
	}
}

void
Antibody::detect_CDR_H3_stem_type() {
	if ( camelid_ ) {
		detect_camelid_CDR_H3_stem_type();
	} else {
		detect_regular_CDR_H3_stem_type();
	}
	return;
} // detect_CDR_H3_stem_type

void
Antibody::detect_camelid_CDR_H3_stem_type() {
	TR << "AC Detecting Camelid CDR H3 Stem Type" << std::endl;

	// extract single letter aa codes for the chopped loop residues
	utility::vector1< char > cdr_h3_sequence;
	for ( core::Size ii = cdrh_[3][1] - 2; ii <= (cdrh_[3][2] + 1); ++ii ) {
		cdr_h3_sequence.push_back( Fv_sequence_[ii] );
	}

	// Rule for extended
	if ( ( ( cdrh_[3][2] - cdrh_[3][1] ) + 1 ) >= 12 ) {
		if ( ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'Y' ) ||
				( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'W' ) ||
				( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'F' ) ) &&
				( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'H' ) &&
				( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'G' ) ) {
			extended_ = true;
		}
	}

	if ( !extended_ ) {
		kinked_ = true;
		if ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] == 'R' ) ||
				( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'Y' ) ||
				(( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'Y' ) ||
				( cdr_h3_sequence[ cdr_h3_sequence.size() - 1 ] != 'W' ) ) &&
				( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'Y' ) ||
				( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'W' ) ) &&
				( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] != 'Y' ) ||
				( cdr_h3_sequence[ cdr_h3_sequence.size() - 3 ] != 'W' ) )) ) {
			kinked_ = false;
		}
	}


	TR << "AC Finished Detecting Camelid CDR H3 Stem Type: "
		<< "Kink: " << kinked_ << " Extended: " << extended_ << std::endl;

	return;

} // detect_camelid_CDR_H3_stem_type()


void
Antibody::detect_regular_CDR_H3_stem_type() {
	TR << "AC Detecting Regular CDR H3 Stem Type" << std::endl;

	bool is_H3( false );

	// extract single letter aa codes for the chopped loop residues
	utility::vector1< char > cdr_h3_sequence;
	for ( core::Size ii = cdrh_[3][1] - 2; ii <= (cdrh_[3][2] + 1); ++ii ) {
		cdr_h3_sequence.push_back( Fv_sequence_[ii] );
	}

	// Rule 1a for standard kink
	if ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] != 'D' ) {
		kinked_ = true;
		is_H3 = true;
	}

	// Rule 1b for standard extended form
	if ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D')
			&& ( (cdr_h3_sequence[2] != 'K') &&
			(cdr_h3_sequence[2] != 'R') ) && (is_H3 != true) ) {
		extended_ = true;
		is_H3 = true;
	}

	if ( !is_H3 ) {
		// Rule 1b extension for special kinked form
		bool is_basic( false ); // Special basic residue exception flag
		for ( core::Size ii = 3; ii <= core::Size(cdr_h3_sequence.size() - 4);
				ii++ ) {
			if ( cdr_h3_sequence[ii] == 'R' || cdr_h3_sequence[ii] == 'K' ) {
				is_basic = true;
				break;
			}
		}

		if ( !is_basic ) {
			core::Size L49_pose_number = Fv.pdb_info()->pdb2pose( 'L', 49 );
			char aa_code_L49 = Fv.residue( L49_pose_number ).name1();
			if ( aa_code_L49 == 'R' || aa_code_L49 == 'K' ) {
				is_basic = true;
			}
		}
		if ( is_basic ) {
			kinked_ = true;
			is_H3 = true;
		}
	}

	// Rule 1c for kinked form with salt bridge
	if ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D') &&
			( (cdr_h3_sequence[2] == 'K') ||
			(cdr_h3_sequence[2] == 'R') ) &&
			( (cdr_h3_sequence[1] != 'K') &&
			(cdr_h3_sequence[1] != 'R') ) && (is_H3 != true) ) {
		kinked_ = true;
		is_H3 = true;
		if ( !is_H3 ) {
			bool is_basic( false ); // Special basic residue exception flag
			core::Size L46_pose_number = Fv.pdb_info()->pdb2pose( 'L', 46 );
			char aa_code_L46 = Fv.residue( L46_pose_number ).name1();
			if ( aa_code_L46 == 'R' || aa_code_L46 == 'K' ) {
				is_basic = true;
			}
			if ( is_basic ) {
				extended_ = true;
				is_H3 = true;
			}
		}
	}

	// Rule 1d for extened form with salt bridge
	if ( ( cdr_h3_sequence[ cdr_h3_sequence.size() - 2 ] == 'D') &&
			( ( cdr_h3_sequence[ 2 ] == 'K') ||
			(cdr_h3_sequence[2] == 'R')) &&
			( (cdr_h3_sequence[1] == 'K') ||
			(cdr_h3_sequence[1] == 'R') ) && (is_H3 != true) ) {
		extended_ = true;
		//is_H3 = true;
	}

	TR << "AC Finished Detecting Regular CDR H3 Stem Type: "
		<< "Kink: " << kinked_ << " Extended: " << extended_ << std::endl;

	return;

} // detect_regular_CDR_H3_stem_type()

void
Antibody::all_cdr_fold_tree() {
	using namespace core::kinematics;

	all_cdr_loops.sequential_order();

	FoldTree f;
	f.clear();

	core::Size jump_num = 0;
	for ( loops::Loops::const_iterator it=all_cdr_loops.begin(),
			it_end=all_cdr_loops.end(),
			it_next; it < it_end; ++it ) {

		it_next = it;
		it_next++;

		if ( it == all_cdr_loops.begin() ) {
			f.add_edge( 1, it->start()-1, Edge::PEPTIDE );
		}

		jump_num++;
		f.add_edge( it->start()-1, it->stop()+1, jump_num );
		f.add_edge( it->start()-1, it->cut(),  Edge::PEPTIDE );
		f.add_edge( it->cut()+1, it->stop()+1, Edge::PEPTIDE );
		if ( it == (it_end-1) ) {
			f.add_edge( it->stop()+1, Fv.size(), Edge::PEPTIDE);
		} else {
			f.add_edge( it->stop()+1, it_next->start()-1, Edge::PEPTIDE );
		}
	}

	f.reorder(1);
	Fv.fold_tree( f );

} // all_cdr_fold_tree()

void
Antibody::align_to_native( Antibody & native ) {

	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, Fv, core::id::AtomID::BOGUS_ATOM_ID() );

	for ( core::Size j = 1; j <= 6; j++ ) {
		core::Size buffer_for_h3_end(0);
		if ( j == 6 ) buffer_for_h3_end = 1;
		for ( core::Size res_counter=hfr_[j][1] + buffer_for_h3_end,
				nat_counter=native.hfr_[j][1] + buffer_for_h3_end;
				res_counter <= hfr_[j][2]; res_counter++, nat_counter++ ) {
			for ( core::Size atm_counter=1; atm_counter <= 4; atm_counter++ ) {
				core::id::AtomID const id1( atm_counter, res_counter );
				core::id::AtomID const id2( atm_counter, nat_counter );
				atom_map[ id1 ] = id2;
			}
		}
	}

	core::scoring::superimpose_pose( Fv, native.Fv, atom_map );

} // align_to_native()


} // namespace antibody
} // namespace protocols

