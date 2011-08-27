// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/Loops.cc
/// @brief
/// @author Chu Wang
/// @author Mike Tyka
/// @author James Thompson

// Unit header
#include <protocols/loops/Loops.hh>

// C++ Headers
#include <iostream>
#include <string>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/NamedAtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

// Package headers
#include <protocols/loops/Loop.hh>

namespace protocols {
namespace loops {

using namespace core;
using namespace ObjexxFCL;

static basic::Tracer tr("loops");
static numeric::random::RandomGenerator RG(430);  // <- Magic number, do not change it (and dont try and use it anywhere else)

std::string get_loop_file_name() {
	using namespace basic::options;
	if ( option[ OptionKeys::loops::loop_file ].user() ) {
		utility::vector1< std::string>  loop_files = option[ OptionKeys::loops::loop_file]();
		if( loop_files.size() == 1 ) return loop_files[1];
		core::Size choice=core::Size( RG.random_range(1,(loop_files.size())  ));
		tr.Error << "Loop choice: " << loop_files[choice] << "  " << choice << std::endl;
		return loop_files[choice];
	}
	return std::string("");
}

Loops get_loops_from_file() {
	using namespace basic::options;

	Loops my_loops;

	if ( option[ OptionKeys::loops::loop_file ].user() ) {
		my_loops.read_loop_file( get_loop_file_name() );
		return my_loops;
	}

	return Loops(); // return empty loop definition if neither option is defined.
}

void Loops::center_of_mass(const core::pose::Pose& pose,
                           numeric::xyzVector<core::Real>* center) const {
  using core::Real;
  using core::Size;
  using core::id::NamedAtomID;

  assert(center);
  center->zero();

  Real count = 0;
  for (const_iterator i = begin(); i != end(); ++i) {
    for (Size j = i->start(); j <= i->stop(); ++j) {
      (*center) += pose.xyz(NamedAtomID("CA", j));
      ++count;
    }
  }

  (*center) /= count;
}

/// @brief switch DOF_Type for residues in loop. id::CHI, id::BB --- don't use
/// with id::JUMP
void
Loops::switch_movemap(
	 kinematics::MoveMap& movemap,
	 id::TorsionType id,
	 bool allow_moves
) const {
	for ( const_iterator it = begin(), eit = end(); it != eit; ++it ) {
		it->switch_movemap( movemap, id, allow_moves );
	}
}


//////////////////////////////////////////////////////////////////////
std::ostream & operator<< ( std::ostream & os, const Loops & loops ) {
	os << "LOOP  begin  end  cut  skip_rate  extended" << std::endl;
	for ( Loops::const_iterator it = loops.begin(), it_end = loops.end();
			 it != it_end; ++it ) {
		os << *it << std::endl;
	}
	return os;
}


void
Loops::write_loops_to_file(
  std::string const & filename,
	std::string token
) const
{

	utility::io::ozstream data;
	data.open( filename );
	if ( !data ) {
		utility_exit_with_message( "Couldn't write loops to file: "+filename );
	}

	write_loops_to_stream( data, token );

	data.close();
	data.clear();
}

void
Loops::write_loops_to_stream(
  std::ostream& data,
	std::string token
) const
{

	for( Loops::const_iterator it= this->begin(), it_end=this->end();
			 it != it_end; ++it ) {
		data << token << " " << it->start() << " " << it->stop() << " " << it->cut() << " "
				 << it->skip_rate() << " " << it->is_extended() << std::endl;
	}
}



void
Loops::add_loop( loops::Loop loop, core::Size minimal_gap ) {
	Size const start( loop.start() );
	Size const stop( loop.stop() );
	Size const cut( loop.cut() );
	tr.Trace << "adding loop " << loop << std::endl;
	if (  ( cut == 0 || ( cut>=start-1 && cut <= stop )) && start <= stop ) {
		for( iterator it = loops_.begin(), it_end = loops_.end();
				 it != it_end; ++it ) {
			// check for conflicts
			if( stop+minimal_gap >= it->start() && start <= it->stop() + minimal_gap ) {
				Loop new_loop( std::min( (int) start, (int) it->start() ), std::max( (int) it->stop(), (int) stop ), it->cut(), it->skip_rate() );
				loops_.erase( it );
				tr.Trace << "overlapping loop found: " << loop << " overlaps with " << *it << " create new loop " << new_loop << std::endl;
				add_loop( new_loop );
				return;
		// 		std::string msg;
// 				msg += "Loops::add_loop error -- overlapping loop regions\n";
// 				msg += "existing loop begin/end: " + string_of( it->start() ) + "/";
// 				msg += string_of( it->stop() ) + "\n";
// 				msg += "new loop begin/end: " + string_of(start) + "/" + string_of(stop);
// 				utility_exit_with_message( msg );
			}
		} // no overlaps
		loops_.push_back( loop );
	} else {
		std::string msg;
		msg += "Loops::add_loop error -- bad loop definition\n";
		msg += "begin/end/cut: " + string_of(start) + "/" + string_of(stop) + "/";
		msg += string_of(cut) + "\n";

		//		runtime_assert( false );
		utility_exit_with_message( msg );
	}
}

//////////////////////////////////////////////////////////////////////
void
Loops::add_loop(
	Size const start,
	Size const stop,
	Size const cut,
	core::Real skip_rate,
	bool const extended
)
{
	add_loop( Loop( start, stop, cut, skip_rate, extended ));
}
/////////////////////////////////////////////////////////////////////////////
void
Loops::add_loop( const Loops::const_iterator & it ) {
	add_loop( it->start(), it->stop(), it->cut(),
		it->skip_rate(), it->is_extended() );
}
/////////////////////////////////////////////////////////////////////////////
void
Loops::add_loop( const Loops::iterator & it ) {
	add_loop( it->start(), it->stop(), it->cut(),
		it->skip_rate(), it->is_extended() );
}
/////////////////////////////////////////////////////////////////////////////
void
Loops::add_overlap_loop( Loops loops ) {
	for( Loops::const_iterator it = loops.begin(), it_end = loops.end();
			 it != it_end; ++it ) {
		add_overlap_loop( *it );
	}
}
/////////////////////////////////////////////////////////////////////////////
void
Loops::add_overlap_loop( const Loop loop ) {

	//    if ( loop.cut() >=loop.start()-1 && loop.cut() <= loop.stop() && loop.start() < loop.stop() ) {
	Size temp_start = loop.start();
	Size temp_stop  = loop.stop();
	Size temp_cut   = loop.cut();
	std::vector<Size> loops_to_delete_start;
	std::vector<Size> loops_to_delete_stop;
	for( const_iterator it = loops_.begin(), it_end = loops_.end();
			 it != it_end; ++it ) {
		// check for conflicts
		if( temp_stop >= it->start() && temp_stop <= it->stop() ){
			temp_stop = it->stop();
			if( std::find( loops_to_delete_start.begin(), loops_to_delete_start.end(),
					it->start() ) == loops_to_delete_start.end() )
				{
					loops_to_delete_start.push_back( it->start() );
					loops_to_delete_stop.push_back( it->stop() );
				}
		}
		if( temp_start <= it->stop() && temp_start >= it->start()) {
			temp_start = it->start();
			if( std::find( loops_to_delete_start.begin(), loops_to_delete_start.end(),
					it->start() ) == loops_to_delete_start.end() )
				{
					loops_to_delete_start.push_back( it->start() );
					loops_to_delete_stop.push_back( it->stop() );
				}
		}
		if( temp_start <= it->start() && temp_stop >= it->stop() ){// include an existing loop
			if( std::find( loops_to_delete_start.begin(), loops_to_delete_start.end(),
					it->start() ) == loops_to_delete_start.end() )
				{
					loops_to_delete_start.push_back( it->start() );
					loops_to_delete_stop.push_back( it->stop() );
				}
		}
	}
	for( Size d = 0; d < loops_to_delete_start.size(); ++d )
		delete_loop(loops_to_delete_start[d], loops_to_delete_stop[d] );
	loops_.push_back( Loop( temp_start, temp_stop, temp_cut,  loop.skip_rate(), loop.is_extended() ) );

	//    } else {
	//      std::cerr << "Loops::add_loop error -- bad loop definition\n"
	//								<< "begin/cut/end: " << loop.start() << "/" << loop.stop() << "/"
	//								<< loop.cut() << std::endl;
	//			runtime_assert( false );
	//     utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	//    }
}

/////////////////////////////////////////////////////////////////////////////
void
Loops::delete_loop(
	Size const start,
	Size const stop
)
{
	runtime_assert( start < stop );

	for( iterator it=loops_.begin(), it_end=loops_.end();
			 it != it_end; ++it ) {
		if ( start == it->start() && stop == it->stop() ) {
			loops_.erase( it );
			break;
		}
	}
}
/////////////////////////////////////////////////////////////////////////////
Loops::const_iterator
Loops::one_random_loop() const {
	Size const size = loops_.size();
	runtime_assert( size > 0 );
	Size index =0;
	Size const end = static_cast< Size >( numeric::random::uniform()*size );
	const_iterator it = loops_.begin();
	while( index != end ) { ++index; ++it; }
	return it;
}
/////////////////////////////////////////////////////////////////////////////
Size
Loops::loop_size(
	Size const num
) const {
	runtime_assert( num > 0 && num <= loops_.size() );
	return loops_[num-1].size();
}
/////////////////////////////////////////////////////////////////////////////
Size
Loops::loop_size() const {
	Size size = 0;
	for( const_iterator it=loops_.begin(), it_end=loops_.end();
			 it != it_end; ++it ) {
		size += it->size();
	}
	return size;
}
/////////////////////////////////////////////////////////////////////////////
bool
Loops::is_loop_residue( Size const seqpos, int const offset ) const
{
	for( const_iterator it=loops_.begin(), it_end=loops_.end(); it != it_end; ++it ) {
		if ( seqpos >= (it->start()+offset) && seqpos <= (it->stop()-offset) ) return true;
	}
	return false;
}

bool
Loops::loop_of_residue( core::Size const seqpos, Loop& loop ) const {
	for ( const_iterator it=loops_.begin(), it_end=loops_.end(); it != it_end; ++it ) {
		if ( seqpos >= it->start() && seqpos <= it->stop() ) {
			loop = *it;
			return true;
		}
	}
	return false;
}

Size
Loops::loop_index_of_residue( core::Size const seqpos ) const {
	Size ct( 1 );
	for ( const_iterator it=loops_.begin(), it_end=loops_.end(); it != it_end; ++it,++ct ) {
		if ( seqpos >= it->start() && seqpos <= it->stop() ) {
			return ct;
		}
	}
	return 0;
}

void
Loops::remove_terminal_loops( pose::Pose const & pose ){
	LoopList new_loops_;
	iterator it_begin = loops_.begin();
	for ( const_iterator it = it_begin, it_end = loops_.end();
				it != it_end; ++it ) {
				if ( !it->is_terminal( pose ) ){
					new_loops_.push_back( *it );
				}
	}
	loops_ = new_loops_;
}

/////////////////////////////////////////////////////////////////////////////
void
Loops::sequential_order()
{
	if ( num_loop() <= 1 ) return;

	LoopList new_loops_;

	iterator it_begin = loops_.begin();
	new_loops_.push_back( *it_begin );

	for ( const_iterator it = ++it_begin, it_end = loops_.end();
				it != it_end; ++it ) {
		bool inserted = false;
		for( iterator it2 = new_loops_.begin(), it2_end = new_loops_.end();
				 it2 != it2_end; ++it2 ) {
			if ( it->start() < it2->start() ) {
				new_loops_.insert( it2, *it );
				inserted = true;
				break;
			}
		}
		if ( ! inserted ) new_loops_.push_back( *it );
	}
	runtime_assert( loops_.size() == new_loops_.size() );
	loops_ = new_loops_;
}
//////////////////////////////////////////////////////////////////////
void
Loops::clear(){
	loops_.clear();
}


void Loops::read_stream_to_END( std::istream &is, bool strict_looprelax_checks, std::string filename /*for error reports*/, std::string LOOP_token ) {
	std::string line;
	int linecount=0;
	int errcount=50; //if we reach 0 we bail!
	while( getline( is, line) ) {
		linecount++;
		std::vector< std::string > tokens ( utility::split( line ) );

		if( tokens.size() > 0 ) {
			if ( tokens[0].substr(0,3) == "END" ) break;
			if ( tokens[0] == LOOP_token ) {
				if ( tokens.size() < 3 ) {
					utility_exit_with_message( "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + " Minimum of 3 tokens necessary (begin, end, cutpoint)"  );
				}
				if ( tokens.size() > 6 ) {
					utility_exit_with_message( "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + " Maximum of 6 tokens allowed (LOOP begin end cutpoint skiprate extended)"  );
				}
				core::Size start_res = (core::Size) atoi(tokens[1].c_str());
				core::Size end_res   = (core::Size) atoi(tokens[2].c_str());
				core::Size cutpt = 0;        // default - let LoopRebuild choose cutpoint
				core::Real skip_rate = 0.0;  // default - never skip
				std::string extend_loop_str;
				bool extend_loop = false;

				if (tokens.size() > 3)
					cutpt = (core::Size) atoi(tokens[3].c_str());
				if (tokens.size() > 4)
					skip_rate = atof(tokens[4].c_str());
				if (tokens.size() > 5){
					if( tokens[5] == "X" ){
						tr.Error << "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + "[WARNING] DEPRECATED old style extended marker X is used" << std::endl;
						extend_loop = true;
						if ( errcount > 0 ) errcount--;
						else {
							utility_exit_with_message( "too many errors in loop-file " + filename );
						}
					}else{
						int extended_token = atoi(tokens[5].c_str());
						if( extended_token == 0 ) extend_loop = false;
						else                      extend_loop = true;
					}
				}
				if ( start_res > end_res || ( start_res==end_res && strict_looprelax_checks ) ) {
					utility_exit_with_message( "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + " Invalid loop definition (start residue " + ( strict_looprelax_checks ? ">=" : ">" )  + " end residue) - ERROR"  );
				} else {
					loops_.push_back( Loop(start_res, end_res, cutpt,  skip_rate, extend_loop) );
				}
			} else if ( tokens[0][0] != '#' ) {
				if (tokens.size() >= 2) {
					tr.Error << "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + "DEPRECATED r++ style loopfile" << std::endl;

					if ( errcount>0 ) errcount--;
					else {
						utility_exit_with_message( "too many errors in loop-file " + filename );
					}

					core::Size start_res = (core::Size) atoi(tokens[0].c_str());
					core::Size end_res   = (core::Size) atoi(tokens[1].c_str());
					core::Size cutpt = 0;        // default - let LoopRebuild choose cutpoint
					core::Real skip_rate = 0.0;  // default - never skip
					bool extend_loop = false;    // default - not extended
					if (tokens.size() > 2)
						cutpt = (core::Size) atoi(tokens[2].c_str());
					if (tokens.size() > 3)
						skip_rate = atof(tokens[3].c_str());
					if (tokens.size() > 4){
						if( tokens[4] == "X" ){
							tr.Error << "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + "[WARNING] DEPRECATED old style extended marker X is used" << std::endl;
							extend_loop = true;
						} else {
							int extended_token = atoi(tokens[4].c_str());
							if ( extended_token == 0 ) extend_loop = false;
							else                extend_loop = true;
						}
					}


					if ( start_res > end_res || ( start_res==end_res && strict_looprelax_checks ) ) {
						utility_exit_with_message( "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + " Invalid loop definition (start residue " + ( strict_looprelax_checks ? ">=" : ">" ) + "end residue) - ERROR"  );
					}

					loops_.push_back( Loop(start_res, end_res, cutpt,  skip_rate, extend_loop) );

				} else {
					tr.Warning << "[WARNING] Skipping line '" << line << "'" << std::endl;
				}
			}
 		}
	} //while
	// sort by start residue
	std::sort( loops_.begin(), loops_.end(), Loop_lt() );
}

void Loops::read_loop_file(
													 std::string filename,
													 bool strict_looprelax_checks,
													 std::string LOOP_token
) {
	clear();
	std::ifstream infile( filename.c_str() );

	if (!infile.good()) {
		utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + filename + "'" );
	}

	read_stream_to_END( infile, strict_looprelax_checks, filename, LOOP_token );

	tr.Warning << "LOOP formats were recently reconciled - with *some* backwards compatibility. Please check your definition files!" << std::endl;
	tr.Warning << "Please check that this is what you intended to read in: " << std::endl;
	tr.Warning << *this;
}

Loops Loops::invert( core::Size total_res ) const {
	Loops core; //the opposite of loops

	if ( loops_.size() == 0 ) {
		core.add_loop( 1, total_res, 0 );
		return core;
	}

	LoopList loops( loops_ ); //or can we assume order ?
	// sort by start residue
	std::sort( loops.begin(), loops.end(), Loop_lt() );
	for ( LoopList::const_iterator it = loops.begin(), eit = loops.end(); it!=eit; ++it )
		tr.Debug << (*it)  << std::endl;

	if ( loops.begin()->start() > 1 ) {
		core.add_loop( 1, loops.begin()->start()-1, 0 );
	}
	LoopList::const_iterator eit, it, next;
	for (  it = loops.begin(), next = loops.begin() + 1, eit = loops.end(); next!=eit; ++it,++next ) {
		if ( it->stop()+1 < next->start()-1 ) {
			tr.Debug << "add " << it->stop()+1 << " " << next->start()-1 << std::endl;
			core.add_loop( it->stop()+1, next->start()-1, 0 );
		} else {
			tr.Debug << "do not add " << it->stop()+1 << " " << next->start()-1 << std::endl;
		}
	}
	if ( it->stop() < total_res ) {
		core.add_loop( it->stop(), total_res, 0 );
	}
	return core;
}


///////////////////////////////////////////////////////////////////////
void
Loops::set_extended( bool input ) {
  for( Loops::iterator it=v_begin(), it_end=v_end(); it != it_end; ++it ) {
		it->set_extended( input );
	}
}



///////////////////////////////////////////////////////////////////////
void
Loops::auto_choose_cutpoints(	core::pose::Pose const & pose ) {
  for( Loops::iterator it=v_begin(), it_end=v_end(); it != it_end; ++it ) {
		it->auto_choose_cutpoint( pose );
	}
}


/////////////////////////////////////////////////////////////////////////
//void
//Loops::find_missing_density(	core::pose::Pose const & pose ) {
//	bool inloop = false;
//	for ( Size pos = 1; pos <= pose->total_residue(); pos++ ) {
//		numeric::xyzVector< core::Real> ca_pos = pose->residue( pos ).atom("CA").xyz();
//		bool good ( true );
//		for ( Size j=1; j<= pose->residue( pos ).natoms(); ++j ) {
//			if ( ( ca_pos - pose->residue( pos ).atom(j).xyz() ).length() > 20 ) {
//				good = false;
//			}
//		}
//		if ( good )
//	}
//	if ( tr.Trace.visible() ) {
//		tr.Trace << "selection of residues for rmsd of " << tag << std::endl;
//		for ( std::list< core::Size >::const_iterator it = selection_.begin(), eit = selection_.end();
//				it != eit; ++it ) {
//			tr.Trace << " " << *it;
//		}
//		tr.Trace << std::endl;
//	}
//}
//


///////////////////////////////////////////////////////////////////////
void
Loops::choose_cutpoints( core::pose::Pose const & pose ) {
  for( Loops::iterator it=v_begin(), it_end=v_end(); it != it_end; ++it ) {
		it->choose_cutpoint( pose );
	}
}


///////////////////////////////////////////////////////////////////////
// @Check loops are compatible with pose.

void
Loops::verify_against( core::pose::Pose const & pose ) const {
	using core::Size;
	Size nres = pose.total_residue();

  for( Loops::const_iterator it=begin(), it_end=end(); it != it_end; ++it ) {
		if ( it->start() <= 0 ){
			tr.Error << "ERROR invalid loop " << it->start() << " " << it->stop() << " " << it->cut() << ": Beginning less than 1" <<  std::endl;
			utility_exit_with_message("LoopRebuild::ERROR Loop definition out of boundary \n" );
		}
		if ( it->stop() > nres ){
			tr.Error << "ERROR invalid loop " << it->start() << " " << it->stop() << " " << it->cut() << ": End more than nres(" << nres << ")" << std::endl;
			utility_exit_with_message("LoopRebuild::ERROR Loop definition out of boundary \n" );
		}
		Size loopbegin_i = std::min(  it->start() , (Size)1 );
		Size loopend_i = std::max(  it->stop() , nres );
		Size cutpt_i = it->cut();
		if ( cutpt_i != 0 && ( cutpt_i > loopend_i || cutpt_i < loopbegin_i ) ) {
			tr.Error << "ERROR invalid loop " << loopbegin_i << " " << loopend_i << " " << cutpt_i << std::endl;
			utility_exit_with_message("LoopRebuild::ERROR Loop definition out of boundary \n" );
		}
	}
}

// @brief Extend a loop .. don't extend across cutpoints in the pose
void Loops::grow_all_loops( core::pose::Pose const & pose, core::Real magnitude ) {
	Loops &loops_ = *this;
	for ( core::Size i=1; i <= loops_.size(); i++ ) {
		grow_loop( pose, loops_[i], magnitude);
	}
}

// @brief Extend a loop
void Loops::grow_all_loops( core::Size nres, core::Real magnitude ) {
	Loops &loops_ = *this;
	for ( core::Size i=1; i <= loops_.size(); i++ ) {
		grow_loop( nres, loops_[i], magnitude);
	}
}

/// @brief Extend a loop .. don't extend across cutpoints in the pose
void Loops::grow_loop(
	core::pose::Pose const & pose,
	Loop & loop,
	core::Real magnitude
) {
	//fpd don't grow across chainbreaks
	//fpd do this by adjusting magnitude in both directions
	core::Real magL=magnitude, magR=magnitude;
	for (int i=0; i<=magnitude; ++i)
		if ( pose.fold_tree().is_cutpoint( loop.start()-i-1 ) ) {
			magL=i; break;
		}
	for (int i=0; i<=magnitude; ++i)
		if ( pose.fold_tree().is_cutpoint( loop.stop()+i ) ) {
			magR=i; break;
		}

	grow_loop( pose.total_residue(), loop, magL, magR );
}


/// @brief Extend a loop
void Loops::grow_loop(
	core::Size nres,
	Loop & loop,
	core::Real magnitude
) {
	grow_loop(nres,loop,magnitude,magnitude);
}

/// @brief Extend a loop
void Loops::grow_loop(
	core::Size nres,
	Loop & loop,
	core::Real magL,
	core::Real magR
) {
  Loop originalloop = loop;

	Loops &loops_ = *this;

	tr.Debug << "GrowLoop: " << loop << std::endl;

	core::Size extend_start = static_cast< core::Size >(
		numeric::random::uniform() * magL
	);
	core::Size extend_stop  = static_cast< core::Size >(
		numeric::random::uniform() * magR
	);

	if ( ( extend_start == 0 ) && ( extend_stop == 0 ) ) {
		if ( numeric::random::uniform() > 0.5 && magL > 0)
			extend_start = 1;
		else if (magR > 0)
			extend_stop  = 1;
		else // magR = magL = 0
			return;
	}

	Loop newloop( loop );
	core::Size new_start = static_cast< core::Size > (
		std::max( 1, (int)loop.start()  - (int)extend_start )
	);
	core::Size new_stop  = static_cast< core::Size > (
		std::min( (int)nres, (int)loop.stop()   + (int)extend_stop )
	);

	tr.Debug << "NewLoop: " << new_start << "  " << new_stop << std::endl;

	// make sure we dont eat into existing other loops
	//
	bool start_is_in_previous_loop = false;
	do {
		if ( new_stop <= new_start+1 ) break;
		if ( start_is_in_previous_loop ) new_start++;
		start_is_in_previous_loop = false;
		for ( Loops::iterator it=loops_.v_begin(), it_end=loops_.v_end();
			it != it_end; ++it
		) {
			if ( it->start() >= new_start ) continue; // ignore that start after this loop
			if ( (*it) == originalloop ) continue; // ignore the original loop
			if ( (it->stop()+1) >= new_start){
				tr.Warning << "Tried growing loop into previous loop" << *it << "  " << new_start << "  " << new_stop << std::endl;
				start_is_in_previous_loop = true;
				break;
			}
		}
	}
	while( start_is_in_previous_loop );


	bool stop_is_in_next_loop = false;
	do {
		if ( new_stop <= new_start+1 ) break;
		if ( stop_is_in_next_loop ) new_stop--;
		stop_is_in_next_loop = false;
		for ( Loops::iterator it=loops_.v_begin(), it_end=loops_.v_end();
			it != it_end; ++it ) {
			// ignore that stop before this loop
			if ( it->stop() <= new_stop ) continue;
			if ( (*it) == originalloop ) continue; // ignore the original loop
			if ( ( it->start()-1) <= new_stop ) {
				tr.Warning << "Tried growing loop into next loop:" << *it << "  " << new_start << "  " << new_stop << std::endl;
			 	stop_is_in_next_loop = true;
				break;
			}
		}
	} while( stop_is_in_next_loop );


	if ( new_stop < loop.stop() ) {
		tr.Info << "Loop stops earlier than before ???" << std::endl;
		new_stop = loop.stop();
	}
	if ( new_start > loop.start() ) {
		tr.Info << "Loop starts later than before ???" << std::endl;
		new_start = loop.start();
	}


	// make sure loop length is greater than 3
	//fpd don't extend across chainbreaks!!
	if ( new_stop - new_start < 2 && new_start != 1 && new_stop != nres ) {
		if (magL>0) new_start -= 1;
		if (magR>0) new_stop  += 1;
	}

	int final_extend_start = loop.start()-new_start;
	int final_extend_stop  = new_stop-loop.stop();

	loop.set_start( new_start );
	loop.set_stop( new_stop );
	tr.Info << "Extended:   (-" << final_extend_start << ",+" << final_extend_stop
		<< ")  " <<  loop << std::endl;
}

void Loops::get_residues( utility::vector1< Size>& selection ) const {
	selection.clear();
	for ( LoopList::const_iterator it = loops_.begin(); it != loops_.end(); ++it ) {
		it->get_residues( selection );
	}
}

bool Loops::operator== ( Loops const& other ) const {
	if ( other.size() != size() ) return false;
	LoopList::const_iterator other_it = other.loops_.begin();
	for ( LoopList::const_iterator it = loops_.begin(); it != loops_.end(); ++it, ++other_it ) {
		if ( *other_it != *it ) return false;
	}
	return true;
}

void
Loops::make_sequence_shift( int shift )
{
  for ( Loops::iterator it = v_begin(), it_end = v_end();
         it != it_end; ++it ) {
     it->set_start( core::Size( it->start() + shift ) );
     it->set_stop( it->stop() + shift );
     it->set_cut( it->cut() + shift );
  }
}

} // namespace loops
} // namespace protocols
