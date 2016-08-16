// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/Loops.cc
/// @brief
/// @author Chu Wang
/// @author Mike Tyka

// Unit header
#include <protocols/loops/Loops.hh>

// Package headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/LoopsFileIO.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/NamedAtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

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
#include <utility/vector1.hh>

// C++ Headers
#include <iostream>
#include <string>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace loops {

using namespace core;
using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer tr( "loops" );

Loops::Loops() : utility::pointer::ReferenceCount()
{
	init( LoopList() );
}

Loops::Loops( const Loops & src ) : utility::pointer::ReferenceCount()
{
	init( src.loops() );
}

Loops::Loops( SerializedLoopList const & src )
{
	init( setup_loops_from_data( src ) );
}

Loops::Loops( std::string const & loop_file_name )
{
	init( LoopList(), false, loop_file_name );
}

Loops::Loops( bool setup_loops_from_options_system )
{
	init( LoopList(), setup_loops_from_options_system );
}

Loops::Loops( utility::vector1< bool > const& selection,
	bool randomize_cutpoints ) {
	init( LoopList(), false );

	bool prev = false;
	Size start = 0;
	for ( Size i = 1; i <= selection.size(); ++i ) {
		if ( selection[i] && !prev ) {
			start = i;
		} else if ( !selection[i] && prev ) {
			assert( start != 0 );
			if ( randomize_cutpoints ) {
				this->add_loop( protocols::loops::Loop( start, i - 1, numeric::random::rg().random_range( (int) start+1, (int) i-1 ) ) );
			} else {
				this->add_loop( protocols::loops::Loop( start, i - 1, (i - start) / 2 ) );
			}
			start = 0;
		}

		prev = selection[i];
	}

	// if the last residue is true, the loop doesn't get added. For a regular loop, I imagine
	// that this is not useful, but this code gets used for other stuff...
	if ( start ) {
		core::Size end = selection.size();
		if ( randomize_cutpoints ) {
			this->add_loop( protocols::loops::Loop( start, end, numeric::random::rg().random_range( (int) start+1, (int) end ) ) );
		} else {
			this->add_loop( protocols::loops::Loop( start, end, (end - start + 1) / 2 ) );
		}
	}

}

// destructor
Loops::~Loops(){}

void Loops::init(
	LoopList const & loops_in,
	bool const read_loop_file_from_options,
	std::string const & passed_in_filename
)
{
	loops_ = loops_in;

	if ( read_loop_file_from_options ) {
		read_loops_options();
		return;
	}
	if ( passed_in_filename.compare( "" ) != 0 ) {
		set_loop_file_name_and_reset( passed_in_filename );
		return;
	}
	loop_filename_ = "";
}

void Loops::read_loops_options()
{
	using namespace basic::options;
	if ( option[ OptionKeys::loops::loop_file ].user() ) {
		utility::vector1< std::string >  loop_files = option[ OptionKeys::loops::loop_file]();
		if ( loop_files.size() == 1 ) {
			set_loop_file_name_and_reset( loop_files[ 1 ] );
			return;
		}
		core::Size choice = core::Size( numeric::random::rg().random_range(1,( loop_files.size() ) ));
		// Why is this tr.Error and not like... tr.Info?
		tr.Error << "Loop choice: " << loop_files[ choice ] << "  " << choice << std::endl;
		set_loop_file_name_and_reset( loop_files[ choice ] );
		return;
	}
	loop_filename_ = "";
}

bool Loops::empty() const { return num_loop() == 0; }
core::Size Loops::num_loop() const { return loops_.size(); }
Loops::const_iterator Loops::begin() const { return loops_.begin(); }
Loops::const_iterator Loops::end() const { return loops_.end(); }
Loops::iterator Loops::v_begin() { return loops_.begin(); }
Loops::iterator Loops::v_end() { return loops_.end(); }

void Loops::center_of_mass(const core::pose::Pose& pose,
	numeric::xyzVector<core::Real>* center) const {
	using core::Real;
	using core::Size;
	using core::id::NamedAtomID;

	assert(center);
	center->zero();

	Real count = 0;
	for ( const_iterator i = begin(); i != end(); ++i ) {
		for ( Size j = i->start(); j <= i->stop(); ++j ) {
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
	for ( Loops::const_iterator it = begin(), eit = end(); it != eit; ++it ) {
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

	for ( Loops::const_iterator it= this->begin(), it_end=this->end();
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
		for ( iterator it = loops_.begin(), it_end = loops_.end(); it != it_end; ++it ) {
			// check for conflicts
			if ( stop+minimal_gap >= it->start() && start <= it->stop() + minimal_gap ) {
				Loop new_loop(
					std::min( (int) start, (int) it->start() ),
					std::max( (int) it->stop(), (int) stop ),
					it->cut(),
					it->skip_rate() );
				tr.Trace << "overlapping loop found: " << loop << " overlaps with " << *it << " create new loop " <<
					new_loop << std::endl;
				loops_.erase( it );
				add_loop( new_loop );
				return;
			}
		} // no overlaps
		loops_.push_back( loop );
	} else {
		std::string msg;
		msg += "Loops::add_loop error -- bad loop definition\n";
		msg += "begin/end/cut: " + string_of(start) + "/" + string_of(stop) + "/";
		msg += string_of(cut) + "\n";

		//runtime_assert( false );
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
Loops::push_back( Loop loop ) {
	add_loop( loop );
}
/////////////////////////////////////////////////////////////////////////////
void
Loops::push_back(
	core::Size const start,
	core::Size const stop,
	core::Size const cut,
	core::Real skip_rate,
	bool const extended  ) {
	add_loop( start, stop, cut, skip_rate, extended );
}
/////////////////////////////////////////////////////////////////////////////
void
Loops::add_overlap_loop( Loops loops ) {
	for ( Loops::const_iterator it = loops.begin(), it_end = loops.end();
			it != it_end; ++it ) {
		add_overlap_loop( *it );
	}
}
/////////////////////////////////////////////////////////////////////////////
void
Loops::add_overlap_loop( Loop const & loop ) {

	//    if ( loop.cut() >=loop.start()-1 && loop.cut() <= loop.stop() && loop.start() < loop.stop() ) {
	Size temp_start = loop.start();
	Size temp_stop  = loop.stop();
	Size temp_cut   = loop.cut();
	std::vector<Size> loops_to_delete_start;
	std::vector<Size> loops_to_delete_stop;
	for ( const_iterator it = loops_.begin(), it_end = loops_.end();
			it != it_end; ++it ) {
		// check for conflicts
		if ( temp_stop >= it->start() && temp_stop <= it->stop() ) {
			temp_stop = it->stop();
			if ( std::find( loops_to_delete_start.begin(), loops_to_delete_start.end(),
					it->start() ) == loops_to_delete_start.end() ) {
				loops_to_delete_start.push_back( it->start() );
				loops_to_delete_stop.push_back( it->stop() );
			}
		}
		if ( temp_start <= it->stop() && temp_start >= it->start() ) {
			temp_start = it->start();
			if ( std::find( loops_to_delete_start.begin(), loops_to_delete_start.end(),
					it->start() ) == loops_to_delete_start.end() ) {
				loops_to_delete_start.push_back( it->start() );
				loops_to_delete_stop.push_back( it->stop() );
			}
		}
		if ( temp_start <= it->start() && temp_stop >= it->stop() ) { // include an existing loop
			if ( std::find( loops_to_delete_start.begin(), loops_to_delete_start.end(),
					it->start() ) == loops_to_delete_start.end() ) {
				loops_to_delete_start.push_back( it->start() );
				loops_to_delete_stop.push_back( it->stop() );
			}
		}
	}
	for ( Size d = 0; d < loops_to_delete_start.size(); ++d ) {
		delete_loop(loops_to_delete_start[d], loops_to_delete_stop[d] );
	}
	loops_.push_back( Loop( temp_start, temp_stop, temp_cut,  loop.skip_rate(), loop.is_extended() ) );

	//    } else {
	//      std::cerr << "Loops::add_loop error -- bad loop definition\n"
	//                << "begin/cut/end: " << loop.start() << "/" << loop.stop() << "/"
	//                << loop.cut() << std::endl;
	//      runtime_assert( false );
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

	for ( iterator it=loops_.begin(), it_end=loops_.end();
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
	while ( index != end ) { ++index; ++it; }
	return it;
}
/////////////////////////////////////////////////////////////////////////////
Size
Loops::loop_size(
	Size const num
) const {
	runtime_assert( num > 0 && num <= loops_.size() );
	return loops_[num].size();
}
/////////////////////////////////////////////////////////////////////////////
Size
Loops::loop_size() const {
	Size size = 0;
	for ( const_iterator it=loops_.begin(), it_end=loops_.end();
			it != it_end; ++it ) {
		size += it->size();
	}
	return size;
}
/////////////////////////////////////////////////////////////////////////////
core::Size Loops::size() const {
	return loops_.size();
}
/////////////////////////////////////////////////////////////////////////////
core::Size Loops::nr_residues() const {
	return loop_size();
}
/////////////////////////////////////////////////////////////////////////////
bool
Loops::is_loop_residue( Size const seqpos, int const offset ) const
{
	for ( const_iterator it=loops_.begin(), it_end=loops_.end(); it != it_end; ++it ) {
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
		if ( !it->is_terminal( pose ) ) {
			new_loops_.push_back( *it );
		}
	}
	loops_ = new_loops_;
}

/////////////////////////////////////////////////////////////////////////////
void Loops::sequential_order() {
	std::sort(loops_.begin(), loops_.end(), RationalLoopComparator());
}
//////////////////////////////////////////////////////////////////////
void
Loops::clear(){
	loops_.clear();
}

Loops::LoopList const & Loops::loops() const { return loops_; }

LoopsFileIOOP Loops::get_loop_file_reader()
{
	return LoopsFileIOOP( new LoopsFileIO );
}

Loops::LoopList Loops::setup_loops_from_data( SerializedLoopList const & loop_data )
{
	LoopList tmp_loops;
	for (  SerializedLoopList::const_iterator it=loop_data.begin(), it_end=loop_data.end(); it != it_end; ++it ) {
		tmp_loops.push_back( Loop( *it ) );
	}
	// sort by start residue
	std::sort( tmp_loops.begin(), tmp_loops.end(), Loop_lt() );
	return tmp_loops;
}

/// @details Soon to be deprecated.
void Loops::read_loop_file()
{
	clear();

	tr << "Loops object initializing itself from pose-numbered loop file named " << loop_file_name() << ". This functionality is soon to be deprecated." << std::endl;
	std::ifstream loopfile( loop_file_name().c_str() );
	PoseNumberedLoopFileReader reader;
	loops_ = setup_loops_from_data( reader.read_pose_numbered_loops_file( loopfile, loop_file_name() ));

	// TODO: Update these warnings/move them to the reader.
	tr.Warning << "LOOP formats were recently reconciled - with *some* backwards compatibility. Please check your definition files!" << std::endl;
	tr.Warning << "Please check that this is what you intended to read in: " << std::endl;
	tr.Warning << *this;
}

/// @detail Given the total number of residues, inverts this set of loops.
/// Note: the inline comments are geared toward comparative modeling, but
/// are not specific to it. Think of "unaligned" as this, "aligned" as its
/// inverse.
Loops Loops::invert(core::Size num_residues) const {
	Loops inv;

	// no unaligned regions
	if ( !num_loop() ) {
		inv.add_loop(1, num_residues);
		return inv;
	}

	Loops copy = *this;

	// aligned region preceding first unaligned region
	const Loop& first = copy[1];
	if ( first.start() != 1 ) {
		inv.add_loop(Loop(1, first.start() - 1));
	}

	// aligned region following last unaligned region
	const Loop& last = copy[copy.num_loop()];
	if ( last.stop() != num_residues ) {
		inv.add_loop(Loop(last.stop() + 1, num_residues));
	}

	// aligned regions in between unaligned regions
	for ( unsigned i = 2; i <= copy.num_loop(); ++i ) {
		const Loop& prev = copy[i - 1];
		const Loop& curr = copy[i];
		inv.add_loop(Loop(prev.stop() + 1, curr.start() - 1));
	}

	inv.sequential_order();
	return inv;
}

bool Loops::has( core::Size const seqpos, int const offset ) const {
	return is_loop_residue( seqpos, offset );
}

void Loops::set_extended( bool input ) {
	for ( Loops::iterator it=v_begin(), it_end=v_end(); it != it_end; ++it ) {
		it->set_extended( input );
	}
}

void Loops::auto_choose_cutpoints( core::pose::Pose const & pose ) {
	for ( Loops::iterator it=v_begin(), it_end=v_end(); it != it_end; ++it ) {
		it->auto_choose_cutpoint( pose );
	}
}

void Loops::choose_cutpoints( core::pose::Pose const & pose ) {
	for ( Loops::iterator it=v_begin(), it_end=v_end(); it != it_end; ++it ) {
		it->choose_cutpoint( pose );
	}
}

// @Check loops are compatible with pose.
void Loops::verify_against( core::pose::Pose const & pose ) const {
	using core::Size;
	Size nres = pose.total_residue();

	for ( Loops::const_iterator it=begin(), it_end=end(); it != it_end; ++it ) {
		if ( it->start() <= 0 ) {
			tr.Error << "ERROR invalid loop " << it->start() << " " << it->stop() << " " << it->cut() << ": Beginning less than 1" <<  std::endl;
			utility_exit_with_message("LoopRebuild::ERROR Loop definition out of boundary \n" );
		}
		if ( it->stop() > nres ) {
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
	for ( int i=0; i<=magnitude; ++i ) {
		if ( pose.fold_tree().is_cutpoint( loop.start()-i-1 ) ) {
			magL=i; break;
		}
	}
	for ( int i=0; i<=magnitude; ++i ) {
		if ( pose.fold_tree().is_cutpoint( loop.stop()+i ) ) {
			magR=i; break;
		}
	}

	grow_loop( pose.total_residue(), loop, magL, magR );
}

void Loops::grow_loop_away_from_sheets(
	core::pose::Pose const & pose,
	Loop & loop,
	core::Real magnitude
) {
	//fpd don't grow across chainbreaks
	//fpd do this by adjusting magnitude in both directions
	//adds a layer that if possible does not grow the loop into beta sheets
	//I allow growth into helices because sometimes
	//there are 2 or 3 residues helices in the middle of loops
	core::Real magL=magnitude, magR=magnitude;
	for ( int i=0; i<=magnitude; ++i ) {
		if ( ( pose.fold_tree().is_cutpoint( loop.start()-i-1 ) ) || (pose.secstruct(loop.start()-i-1) == 'E') ) {
			magL=i; break;
		}
	}
	for ( int i=0; i<=magnitude; ++i ) {
		if ( ( pose.fold_tree().is_cutpoint( loop.stop()+i ) ) ||  (pose.secstruct(loop.start()+i) == 'E') ) {
			magR=i; break;
		}
	}
	//if both sides hit extended structure than default to the cutpoint behavior
	if ( (magR == 0) && (magL == 0) ) {
		for ( int i=0; i<=magnitude; ++i ) {
			if ( pose.fold_tree().is_cutpoint( loop.start()-i-1 ) ) {
				magL=i; break;
			}
		}
		for ( int i=0; i<=magnitude; ++i ) {
			if ( pose.fold_tree().is_cutpoint( loop.stop()+i ) ) {
				magR=i; break;
			}
		}
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
		if ( numeric::random::uniform() > 0.5 && magL > 0 ) {
			extend_start = 1;
		} else if ( magR > 0 ) {
			extend_stop  = 1;
		} else { // magR = magL = 0
			return;
		}
	}

	Loop newloop( loop );
	core::Size new_start = static_cast< core::Size > (
		std::max( 1, (int)loop.start()  - (int)extend_start )
	);
	core::Size new_stop  = static_cast< core::Size > (
		std::min( (int)nres, (int)loop.stop()   + (int)extend_stop )
	);

	tr.Debug << "NewLoop before loop check: " << new_start << "  " << new_stop << std::endl;

	//grow loops to the start of previous or next loop
	for ( Loops::iterator it=loops_.v_begin(), it_end=loops_.v_end();
			it != it_end; ++it ) {
		if ( (*it) != originalloop ) {
			//case where the start has grown into the previous loop
			if ( (new_start >= it->start())&&(new_start <= it->stop()) ) {
				new_start = it->start();
				tr.Warning << "Tried growing loop into previous loop:" << *it << "  " << new_start << "  " << new_stop << std::endl;
			}
			//case where the stop has grown into the next loop
			if ( (new_stop <= it->stop())&&(new_stop >= it->start()) ) {
				new_stop = it->stop();
				tr.Warning << "Tried growing loop into next loop:" << *it << "  " << new_start << "  " << new_stop << std::endl;
			}
		}
	}

	tr.Debug << "NewLoop after loop check: " << new_start << "  " << new_stop << std::endl;

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
		if ( magL>0 ) new_start -= 1;
		if ( magR>0 ) new_stop  += 1;
	}

	int final_extend_start = loop.start()-new_start;
	int final_extend_stop  = new_stop-loop.stop();

	loop.set_start( new_start );
	loop.set_stop( new_stop );
	tr.Info << "Extended:   (-" << final_extend_start << ",+" << final_extend_stop
		<< ")  " <<  loop << std::endl;
}

void Loops::get_residues( utility::vector1< Size >& selection ) const {
	selection.clear();
	for ( const_iterator it = loops_.begin(); it != loops_.end(); ++it ) {
		it->get_residues( selection );
	}
}
std::string const & Loops::loop_file_name()
{
	return loop_filename_;
}
void Loops::set_loop_file_name_and_reset( std::string const & loop_filename )
{
	loop_filename_ = loop_filename;
	read_loop_file();
}

bool Loops::strict_looprelax_checks()
{
	return strict_looprelax_checks_on_file_reads_;
}

void Loops::set_strict_looprelax_checks( bool const check )
{
	strict_looprelax_checks_on_file_reads_ = check;
}

std::string const & Loops::file_reading_token()
{
	return file_reading_token_;
}

void Loops::set_file_reading_token( std::string const & token )
{
	file_reading_token_ = token;
}

const Loop & Loops::operator[] ( core::Size const i ) const
{
	return loops_[i];
}

Loop & Loops::operator[] ( core::Size const i )
{
	return loops_[i];
}

Loops & Loops::operator =( Loops const & src )
{
	loops_ = src.loops_;
	return *this;
}

bool Loops::operator== ( Loops const& other ) const
{
	if ( other.size() != size() ) return false;
	const_iterator other_it = other.loops_.begin();
	for ( const_iterator it = loops_.begin(); it != loops_.end(); ++it, ++other_it ) {
		if ( *other_it != *it ) return false;
	}
	return true;
}

bool Loops::operator!=( Loops const& other ) const
{
	return !( *this == other);
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

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::loops::Loops::save( Archive & arc ) const {
	arc( CEREAL_NVP( loops_ ) ); // LoopList
	//arc( CEREAL_NVP( loop_file_reader_ ) ); // LoopsFileIOOP
	arc( CEREAL_NVP( strict_looprelax_checks_on_file_reads_ ) ); // _Bool
	arc( CEREAL_NVP( loop_filename_ ) ); // std::string
	arc( CEREAL_NVP( file_reading_token_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::loops::Loops::load( Archive & arc ) {
	arc( loops_ ); // LoopList
	//arc( loop_file_reader_ ); // LoopsFileIOOP
	arc( strict_looprelax_checks_on_file_reads_ ); // _Bool
	arc( loop_filename_ ); // std::string
	arc( file_reading_token_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::loops::Loops );
CEREAL_REGISTER_TYPE( protocols::loops::Loops )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_loops_Loops )
#endif // SERIALIZATION
