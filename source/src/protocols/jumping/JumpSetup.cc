// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @details
/// @author Oliver Lange


// Unit Headers
#include <protocols/jumping/JumpSetup.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/fragment/FrameList.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/OrderedFragSet.hh>
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif

//#include <core/scoring/constraints/ConstraintForest.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

//numeric headers
#include <numeric/random/random.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <fstream>


namespace protocols {
namespace jumping {

/// @details Auto-generated virtual destructor
BaseJumpSetup::~BaseJumpSetup() {}

using namespace core;
using namespace fragment;
using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer tr( "protocols.jumping" );


FragSetOP
BaseJumpSetup::generate_jump_frags( JumpSample const& jumps, kinematics::MoveMap const& mm ) const {
	OrderedFragSetOP frags( new OrderedFragSet );
	FrameList jump_geometries;
	jumps.generate_jump_frags( *jumping::StandardPairingLibrary::get_instance(),
		mm,
		true, /*bWithBBTorsions*/
		jump_geometries
	);
	frags->add( jump_geometries );
	return frags;
}

void
JumpSetup::read_file( std::string fname ) {
	utility::io::izstream data( fname.c_str() );
	tr.Info << "read jump-definitions from " << fname << std::endl;
	if ( !data ) {
		std::cerr << "ERROR:: Unable to open constraints file: "
			<< fname << std::endl;
		std::exit( 1 );
	}

	std::string line;

	while ( getline(data, line ) ) {
		std::istringstream in( line );
		Interval jump, cuts;
		in >> jump.start_ >> jump.end_ >> cuts.start_;
		if ( in ) {
			in >> cuts.end_;
		} else {
			cuts.end_=cuts.start_;
		}
		//  if ( in ) {
		//   tr.Warning << "Ignore remaining line: " << line << std::endl;
		//  }
		using namespace ObjexxFCL::format;
		tr.Debug << "read Jumps: " << RJ(4, jump.start_) << RJ( 4, jump.end_ ) << RJ( 4, cuts.start_ ) << RJ( 4, cuts.end_ )<< std::endl;
		add_jump( jump, cuts );
	}

}

//* ----------------------- JumpSelector ----------------------------*

JumpSelector::JumpSelector()
: total_weight_( 0 ), min_loop_length_( 5 ), loop_extension_(2), nr_jumps_min_(3), nr_jumps_max_(6) {}
JumpSelector::JumpSelector( std::string ss )
: secstruct_( ss), total_weight_( 0 ), min_loop_length_( 5 ), loop_extension_(2), nr_jumps_min_(3), nr_jumps_max_(6) {}

JumpSelector::~JumpSelector() {}


void
JumpSelector::read_file( std::string fname ) {
	utility::io::izstream data( fname.c_str() );
	tr.Info << "read jump-definitions from " << fname << std::endl;
	if ( !data ) {
		utility_exit_with_message("ERROR: could not open file " + fname );
	}

	std::string line;

	while ( getline(data, line ) ) {
		if ( line.find("nr_jumps") !=std::string::npos ) {
			std::istringstream in( line );
			std::string tag;
			in >> tag >> nr_jumps_min_ >> nr_jumps_max_;
			tr.Debug << "read nr_jumps: min " << nr_jumps_min_ << "... max " << nr_jumps_max_ << std::endl;
			continue;
		} else if ( line.find("loop_params") !=std::string::npos ) {
			std::istringstream in( line );
			std::string tag;
			in >> tag >> min_loop_length_ >> loop_extension_;
			tr.Debug << "read min_loop_length: " << min_loop_length_ << " loop_extension " << loop_extension_ << std::endl;
			continue;
		}
		std::istringstream in( line );
		Interval jump; Real weight;
		in >> jump.start_ >> jump.end_ >> weight;
		if ( jump.start_ > jump.end_ ) {
			Size dum = jump.start_;
			jump.start_ = jump.end_;
			jump.end_ = dum;
		};
		using namespace ObjexxFCL::format;
		tr.Debug << "read Jumps: " << RJ(7, jump.start_) << RJ( 7, jump.end_ ) << RJ( 7, weight ) <<std::endl;
		add_jump( jump, weight );
	}
}

void
dump_tags_( FArray1D_int const& tags, Size nres, std::ostream& out ) {
	out << "r_fold_cst: ";
	for ( Size i = 1; i<=nres; i++ ) {
		if ( (i-1)%10 == 0 ) { out << i; continue; }
		//large numbers take several characters... skip appropriate
		if ( (i>=10) && (i-2)%10 == 0 ) { continue; }
		if ( (i>=100) && (i-3)%10 == 0 ) { continue; }
		if ( (i>=1000) && (i-4)%10 == 0 ) { continue; }
		out << ".";
	}
	out << "\nr_fold_cst: ";
	for ( Size i = 1; i<=nres; i++ ) {
		if ( tags( i ) <= -10 ) out << 'C'; //cuttable
		if ( tags( i ) == -1 ) out << 'L'; //short loop
		if ( tags( i ) == 0 ) out << '.';
		if ( tags( i ) > 0 ) out << tags( i );
	}
	out << std::endl;
}

void
fill_tags_( FArray1D_int& tags, Size fill_pos, int tag, Size nres ) {
	Size pos = fill_pos;
	while ( tags( pos ) >= -1 ) {
		tags( pos ) = tag;
		if ( --pos  == 0 ) break;
		//dump_tags_( tags, nres, std::cout );
	}

	runtime_assert( pos == 0 || tags( pos ) == -10 ); //it should always be stopped by loop
	pos = fill_pos;
	while ( tags( pos ) >= -1 ) {
		tags( pos++ ) = tag;
		if ( pos > nres ) break;
		// dump_tags_( tags, nres, std::cout );
	}
	runtime_assert( pos > nres || tags( pos ) <= -10 );
}


// Utility headers

JumpSample
JumpSelector::create_jump_sample( ) const {
	typedef utility::vector1< Interval > CutList;
	typedef utility::vector1< Interval > JumpList;
	JumpList jump_list;
	//init tag-array
	// negative for loop regions
	// -1 for short loops(no jumps), -10 for long loops (cuttable)
	// 0 for free regions
	runtime_assert( min_loop_length_ > 0 );

	Size nr_jumps = nr_jumps_min_ + static_cast< int >( numeric::random::rg().uniform() * (nr_jumps_max_-nr_jumps_min_) );
	tr.Info << "generate " << nr_jumps << "jumps for following secstruct:\n";
	tr.Info << secstruct_ << std::endl;

	Size nres = secstruct_.size();
	FArray1D_int tags( nres , 0 );
	for ( Size i = 0; i < nres; i++ ) {
		if ( secstruct_[i] == 'L' ) {
			for ( int ii=std::max((int)i-(int)loop_extension_,0); ii<= std::min( (int) i + (int)loop_extension_, (int) nres-1); ii++ ) {
				tags(ii+1)= -10;
			}
		}
	}
	if ( tr.Trace.visible() ) dump_tags_( tags, nres, std::cout );

	// detect loops that are too short for cut-points
	for ( Size i = 1; i <= nres; i++ ) {
		if ( tags( i ) == -10 ) {
			bool bShort = false;

			// extend loops to each side
			for ( Size j = i; j <= i+min_loop_length_ && j <= nres; j++ ) {
				if ( tags( j ) > -10 ) {
					bShort = true;
					break;
				}
			}
			while ( tags( i ) == -10  ) { //go to end of current loop
				if ( bShort ) tags( i ) = -1;
				i++;
				if ( i > nres ) break;
			}
		}
	}

	// detect non-loop regions that are too small.
	Size min_free_size_ = 1;
	for ( Size i = 1; i<=nres; i++ ) {
		if ( tags( i ) == 0 ) {
			Size ct = 0;
			while ( tags( i+ct ) == 0 ) { ct++; if ( ct+i > nres ) break; }
			if ( ct <= min_free_size_ ) {
				for ( Size jj=0; jj<ct; jj++ ) {
					tags( i+jj ) = -1;
				}
			}
			i += (ct-1);
		}
	}

	if ( tr.Trace.visible() ) dump_tags_( tags, nres, std::cout );

	Size nr = 0;
	Size attempts = 1;
	while ( nr < nr_jumps  && attempts < 1000 ) {
		attempts += 1;
		Interval aJump = select_random();
		int start_tag = tags( aJump.start_ );
		int stop_tag = tags( aJump.end_ );
		//  bool bReTag = start_tag > 0 && stop_tag > 0;
		if ( start_tag == stop_tag && start_tag > 0 ) continue; //disallow jump
		if ( start_tag < 0 || stop_tag < 0 ) continue; //don't jump into loop regions
		++nr; // from now on we know that the jump will be inserted

		int new_tag =  start_tag;
		int replace_tag = stop_tag;
		if ( new_tag < stop_tag ) {
			new_tag = stop_tag;
			replace_tag = start_tag;
		}
		// now new_tag should be the larger of the start/stop-tags
		if ( new_tag == 0 ) new_tag = nr; //still zero: both start and stop are zero
		if ( replace_tag > 0 ) { //both start/stop > 0 needs to retag
			/// replace all < replace_tag > with < new_tag >
			for ( Size i = 1; i <= nres; i++ ) {
				if ( tags( i ) == replace_tag ) tags(i) = new_tag;
			}
		}
		tr.Info << "after " << attempts << " attempts, jump at " << aJump.start_ << " -- " << aJump.end_ << "  with tag " << new_tag << std::endl ;
		fill_tags_( tags, aJump.start_, new_tag, nres );
		fill_tags_( tags, aJump.end_, new_tag, nres );
		if ( tr.Trace.visible() ) dump_tags_( tags, nres, std::cout );
		jump_list.push_back( aJump );
	}
	if ( attempts >= 1000 ) {
		tr.Warning << "failed to find " << nr_jumps << " jumps in "
			<<attempts<<" attempts. Only " << nr << " possible jumps found\n";
	}
	typedef utility::vector1< CutList > MetaCutList; //take as many cuts from MetaCutList as jumps -- take one cutregion from each MetaCuts
	// determine cut-regions
	// not quite optimal, yet
	// a 1 CC 2 CC 2 CC 1  jump-setup will always be cut in first two loops, but it should sample both cut possibilities for jump 1


	MetaCutList all_cuts;
	std::map< int, bool > tag_list; // store all observed tags, only cut if one of these is found in downstream region (towards C-term)
	for ( Size pos=1; pos <= nres; pos ++ ) {
		if ( tags( pos ) > 0 ) {
			// add tag to tag_list
			tag_list[ tags( pos ) ] = true; // boolean doesn't really matter
			CutList cut_list;
			// move to end of jump region
			while ( tags( pos ) > 0 ) { ++pos; if ( pos > nres ) break; };
			if ( pos > nres ) break;

			runtime_assert ( tags( pos ) <= 0 );
			// go up until you find another jump region -- add all loops in between to possible cut-regions
			// if you find chain-break before next jump region no cut-point is needed.
			while ( tags( pos ) <= 0 ) {

				// skip possible white-space
				while ( tags( pos ) > -10 ) { ++pos; if ( pos > nres ) break; };
				if ( pos > nres ) break;

				// now we are at beginning of new cut region
				runtime_assert ( tags( pos ) == -10 );
				Interval cut;
				cut.start_ = pos;

				// skip to end of cutable-region
				while ( tags( pos ) == -10 ) { ++pos; if ( pos > nres ) break; };
				if ( pos > nres ) break;

				runtime_assert( tags( pos-1 ) == -10 && tags( pos ) > -10 );
				cut.end_ = pos-1;
				cut_list.push_back( cut );
				tr.Trace << "cut: possible " << cut.start_ << " " << cut.end_ << std::endl;
			}
			if ( pos > nres ) break;
			runtime_assert( tags( pos ) > 0 );
			// new jump region-- add recent collection of cut-regions to meta-list
			// only add if down-stream jumps are connected to up-stream
			bool bConnected = false;
			for ( std::map< int , bool>::const_iterator it = tag_list.begin(), eit = tag_list.end();
					it!=eit && !bConnected;
					++it ) {
				for ( Size ii = pos; ii<=nres && !bConnected; ii++ )  {
					if ( it->first == tags ( ii ) ) bConnected = true;
					if ( bConnected ) tr.Trace << "add cuts since tag " << it->first << " has been found at pos " << ii << std::endl;
				}
			}
			if ( bConnected ) all_cuts.push_back( cut_list );
		} // tags (pos) > 0
	} // for loop

	runtime_assert( all_cuts.size() >= jump_list.size() );
	MetaCutList::const_iterator cut_pool_it = all_cuts.begin();
	JumpSetup my_jumps( nres );
	for ( JumpList::const_iterator it = jump_list.begin(), eit = jump_list.end();
			it!=eit; ++it ) {
		int cut_frame_num = static_cast< int >( numeric::random::rg().uniform() * cut_pool_it->size() ) + 1;
		tr.Trace << "chosen jump-cut: " << it->start_ << " " << it->end_ << " cut:" << (*cut_pool_it)[cut_frame_num].start_
			<< " " << (*cut_pool_it)[cut_frame_num].end_ << std::endl;
		my_jumps.add_jump( *it, (*cut_pool_it)[cut_frame_num] );
		++cut_pool_it;
	}
	return JumpSample(my_jumps); //create random cut-points in given cut-regions
} // create_jump_sample

Interval
JumpSelector::select_random() const {
	// FArray1D_int freq(size(),0);
	// Size Nsample = 100000;
	// for ( int ii=1; ii < Nsample; ii++ ) {

	Real ran = numeric::random::rg().uniform() * total_weight_;
	Real cumsum = 0.0;
	const_iterator it = begin();
	const_iterator eit = end();

	int ct = 1;
	cumsum += it->weight_;
	while ( cumsum < ran  ) {
		runtime_assert ( it != eit );
		++it; ct++;
		cumsum += it->weight_;
	}
	// freq( ct ) += 1;
	runtime_assert( it != eit );

	// }
	// const_iterator it = begin();
	// for ( int i = 1; i<= size(); i++ ) {
	//  std::cout << 1.0*freq( i )/ ( 1.0*Nsample) << " " << it->weight_/total_weight_ << std::endl;
	//  ++it;
	// }
	return it->jump_;
}

//JumpsFromConstraintForest::JumpsFromConstraintForest(
// const core::Size total_residue,
// core::scoring::constraints::ConstraintForestOP forest,
// ObjexxFCL::FArray1D_float const& cut_probability
//) :
// total_residue_(total_residue),
// forest_(forest),
// cut_prob_(cut_probability)
//{}
//
//JumpSample
//JumpsFromConstraintForest::create_jump_sample() const {
// // STOPGAP SOLUTION: generate new jumps each time.  Need to take into account
// // that constraints will change too.
// // THIS HAS TO CHANGE.
// //forest_->generate_random_sample(); this is and has to be called where constraints are generated
// return JumpSample(total_residue_, forest_->get_pairings(), cut_prob_);
//}

//JumpsFromConstraintForest::~JumpsFromConstraintForest() { }
//
//JumpSample
//JumpsFromConstraintForest::clean_jumps( JumpSample const& js ) const
//{
// std::cerr << "ERROR: JumpSetup::clean_jumps() not implemented" << std::endl;
// return js;
//}


} //jumping
} //protocols
