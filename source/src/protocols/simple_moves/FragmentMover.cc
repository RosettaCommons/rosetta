// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Inserts a Fragment into a Pose, similar to old Rosetta++ main_frag_trial algorithm.
/// @author Oliver Lange

// Unit Headers
#include <protocols/simple_moves/FragmentMover.hh>

// Package Headers

// Project Headers
#include <core/fragment/Frame.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/kinematics/MoveMap.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>

// ObjexxFCL Headers

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>




// C++ headers

namespace protocols {
namespace simple_moves {

static numeric::random::RandomGenerator RG(489);  // <- Magic number, do not change it!

using namespace core;
using namespace fragment;
using namespace basic;

static basic::Tracer tr("protocols.simple_moves.FragmentMover");

FragmentMover::~FragmentMover() {}

/// @brief Empty constructor
FragmentMover::FragmentMover(std::string type)
{
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_bb(true); // standard MoveMap
	movemap_= movemap;
	protocols::moves::Mover::type(type);
	update_insert_map();
}

FragmentMover::FragmentMover(
	core::fragment::FragSetCOP fragset,
	std::string type
) :
	fragset_( fragset )
	//		movemap_( movemap )//,
	//bValidInsertMap_ ( false )
{
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_bb( true ); //standard movemap
	movemap_=movemap;
	protocols::moves::Mover::type( type );
	update_insert_map();
}

///@brief constructor
FragmentMover::FragmentMover(
	core::fragment::FragSetCOP fragset,
	core::kinematics::MoveMapCOP movemap,
	std::string type
) :
	fragset_( fragset ),
	movemap_( movemap )//,
	//bValidInsertMap_ ( false )
{
	protocols::moves::Mover::type( type );
	update_insert_map();
}

bool FragmentMover::apply( pose::Pose & pose, Size pos ) const {
	FrameList frames;
	if ( !fragset_->frames( pos, frames ) ) return false;
	return apply_frames( pose, frames );
}

std::string
FragmentMover::get_name() const {
	return "FragmentMover";
}

Size FragmentMover::apply_at_all_positions( core::pose::Pose& pose ) const {
	Size ct( 0 );
	for ( InsertMap::const_iterator it = insert_map().begin(), eit = insert_map().end(); it != eit; ++it ) {
		FrameList frames;
		if ( !fragset_->frames( *it, frames ) ) continue;
		for ( FrameList::const_iterator fit = frames.begin(); fit != frames.end(); ++fit ) {
			FrameList one_frame;
			one_frame.push_back( *fit );
			apply_frames( pose, one_frame );
		}
		//		ct += apply( pose, *it );
	}
	return ct;
}

///@brief accessor to the fragment set
core::fragment::FragSetCOP FragmentMover::fragments() const {
	return fragset_;
}

///@brief setter for the fragment set
void
FragmentMover::set_fragments( core::fragment::FragSetCOP new_frags_ ) {
	fragset_ = new_frags_;
  Size size_of_frags = fragset_->max_pos() - fragset_->min_pos() + 1;
  tr.Debug << " got new fragments with size " << size_of_frags << std::endl;
	on_new_fragments();
}

///@brief setter for the movemap
void
FragmentMover::set_movemap( core::kinematics::MoveMapCOP movemap ) {
	movemap_ = movemap;
	update_insert_map();
}

core::kinematics::MoveMapCOP
FragmentMover::movemap() const {
	return movemap_;
}

void
FragmentMover::update_insert_map() {
	//	if ( !bValidInsertMap_ ) {
	fragset_->generate_insert_map( *movemap_, insert_map_, insert_size_ );
	//	bValidInsertMap_ = true;
	//	}
}


// Empty constructor
ClassicFragmentMover::ClassicFragmentMover() : FragmentMover( "ClassicFragmentMover" )
{
	set_defaults();
}

// constructor
ClassicFragmentMover::ClassicFragmentMover(
	core::fragment::FragSetCOP fragset
)	: FragmentMover( fragset, "ClassicFragmentMover" )
{
	set_defaults();
}

// constructor
ClassicFragmentMover::ClassicFragmentMover(
	core::fragment::FragSetCOP fragset,
	core::kinematics::MoveMapCOP movemap
)	: FragmentMover( fragset, movemap, "ClassicFragmentMover" )
{
	set_defaults();
}

// constructor Temp work around for PyRosetta code, until we found a way how to handle owning pointers in this case
ClassicFragmentMover::ClassicFragmentMover(
	core::fragment::ConstantLengthFragSet const & fragset,
	core::kinematics::MoveMap const & movemap
)	: FragmentMover(fragset.clone(), new core::kinematics::MoveMap(movemap), "ClassicFragmentMover" )
{
	set_defaults();
}

ClassicFragmentMover::~ClassicFragmentMover()
{}

// alternative Constructor to be used by derived classes
ClassicFragmentMover::ClassicFragmentMover(
	core::fragment::FragSetCOP fragset,
	core::kinematics::MoveMapCOP movemap,
	std::string type
)	: FragmentMover( fragset, movemap, type )
{
	set_defaults();
}

// alternative Constructor to be used by derived classes
ClassicFragmentMover::ClassicFragmentMover(
	core::fragment::FragSetCOP fragset,
	std::string type
)	: FragmentMover( fragset, type )
{
	set_defaults();
}

std::string
ClassicFragmentMover::get_name() const {
	return "ClassicFragmentMover";
}

void
ClassicFragmentMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "-------------------Settings--------------------" << std::endl;
	output << "End bias:            " << end_bias_ << std::endl <<
					"Min overlap:         " << min_overlap_ << std::endl <<
					"Min fragment length: " << min_frag_length_ << std::endl <<
					"Check ss:            " << ( (check_ss_) ? "True" : "False" ) << std::endl <<
					"bApplyEndBias:       " << ( ( bApplyEndBias_ ) ? "True" : "False" ) << std::endl <<
					"Use predefined window start: " << ( (use_predefined_window_start_) ? "True": "False" ) << std::endl <<
					"Predefined window start:     " << predefined_window_start_ << std::endl;
	output << "-----------------------------------------------" << std::endl;
	output << "Movemap: " << std::endl;
	movemap()->show();
	output << "**Unless a movemap is specified above, all backbone torsion angles are set to TRUE**" << std::endl;
}

protocols::moves::MoverOP
ClassicFragmentMover::clone() const
{
	return new ClassicFragmentMover(*this);
}

protocols::moves::MoverOP
ClassicFragmentMover::fresh_instance() const
{
	return new ClassicFragmentMover();
}

//return a fragnum for given Frame, overload to make other choices
bool ClassicFragmentMover::choose_fragment(
	FrameList const& frames,
	pose::Pose const&,
	Size& frame_num,
	Size& frag_num
) const {
	// classically: choose randomly
	runtime_assert( frames.size() );
	for ( Size nfail = 1; nfail <= 100; nfail ++ ) {

		//choose frame
		frame_num = static_cast< int >( RG.uniform() * frames.size() ) + 1;
		Size N ( frames[ frame_num ]->nr_frags() );

		//choose frag_num in frame
		if ( N >= 1 ) { // nr_frags is indexed starting at 1
			frag_num = static_cast< int >( RG.uniform() * N ) + 1;
			return true;
		}
	}
	return false;
}

void
ClassicFragmentMover::set_defaults() {
	using namespace basic::options;
	check_ss_ = !option[ basic::options::OptionKeys::run::remove_ss_length_screen ]();

	bApplyEndBias_ = true;
	end_bias_ = 30.0; //classic is 60  // pose_simple_moves is 30
	min_overlap_ = 0;
	min_frag_length_ = 0;
	use_predefined_window_start_ = false;
}

// accept with probability 1 if the fragment window is centered on the center of the protein.
// accept with probability .3677 if the fragment window is centered end-bias residues away from
// the center of the protein
//		if ( total_insert+frag_length != pose.total_residue() ||  r <= std::exp( -( end_dist / end_bias ) ) ) {
// the question of total_insert+frag_length == pose.total_residue() doesn't make sense if different frag_lengths are involved
bool ClassicFragmentMover::end_bias_check( core::pose::Pose const& pose, Size begin ) const {
	Real r = RG.uniform();
	// classic bias
	// Real const end_bias ( 60.0 );
	// Real end_dist = std::abs( begin - ( pose.total_residue() / 2.0 ) );
	// return r <= std::exp( -( end_dist / end_bias ) );
	runtime_assert( begin > 0 && begin <= insert_size_.size() );
	Size size = insert_size_[ begin ];

	// the following assertion can happen if non-continuous (eg. Jump) fragments are inserted with bias-check, switch it off!

	if ( ( begin + size - 1 ) > pose.total_residue() ){
		tr.Error << "BEGIN: " << begin << " SIZE: " << size  << " TOTAL_RES: " << pose.total_residue() << std::endl;
		tr.Error << "Are the fragments compatible with the fasta or the input PDB used to extract the folding sequence ? " << std::endl;
		tr.Error << "It appears that the fragments go up to residue " << begin + size - 1 << " while the pose only has " << pose.total_residue() << " residues!" << std::endl;
		utility_exit_with_message("Assertion failure: runtime_assert( ( begin + size - 1 ) <= pose.total_residue() ); " );
	}


	Size min_fixed_residues;


	// THIS HAS TO TAKE MOVEMAP INTO ACCOUNT: if the core around jumps is fixed... basically no moves are accepted...
	Size const fixed_residues =
		pose.fold_tree().count_fixed_residues( begin, size, min_fixed_residues );

//if symmetric, we just consider a single subunit;
	Real factor ( 1.0 );
	if( core::pose::symmetry::is_symmetric( pose ) ) {
		factor = 1.0 / (Real) core::pose::symmetry::symmetry_info(pose)->subunits();
  }

	Real bias = std::exp( factor * ( (Real) min_fixed_residues - fixed_residues ) / end_bias() );
	//	tr.Trace << "biascheck: " << begin << " " << size << " " << factor << " " << min_fixed_residues << " " << fixed_residues << " " << end_bias() << " ("<<bias
	//					 << " => " << r << ")" << std::endl;
	return ( r <= bias );
}

void
ClassicFragmentMover::define_start_window( Size window_start ) {
	use_predefined_window_start_ = true;
	predefined_window_start_ = window_start;
	return;
}

bool ClassicFragmentMover::choose_window_start( pose::Pose const& pose, Size, Size &begin ) const {

	Size const total_insert ( insert_map_.size() );

	for ( Size i = 1; i<=total_insert; i++ ) tr.Trace << " " << insert_map_[ i ];

	if ( tr.Trace.visible() ) {
		tr.Trace << "size of insertmap: " << total_insert << " -- ";
		tr.Trace << "insert_size: ";
		if ( total_insert ) for ( Size i = 1; i<=insert_map_[ total_insert ]; i++ ) tr.Trace << " " << insert_size_[ i ];
		tr.Trace << std::endl;
	}

	if ( !total_insert ) {
		tr.Warning << "empty insert map ... no fragment insertion attempted" << std::endl;
		return false;
	}

	Size nfail ( 0 );
	while ( nfail < 100 ) {

		begin = insert_map_[ static_cast< int >( RG.uniform() * total_insert  ) + 1 ];
		//tr.Trace << "window start " << begin << std::endl;
		/// apl -- distance that the center of the fragment window is from the center of the protein
		/// apl -- SOON
		//Real end_dist = std::abs( (begin + frag_length / 2.0 ) - ( pose.total_residue() / 2.0 ) );
		if ( !bApplyEndBias_ || end_bias_check( pose, begin ) ) {
			return true;
		}
		//tr.Trace << "bias check failed for begin=" << begin << " " << std::endl;
		++nfail;
	}
	return false;
}

bool ClassicFragmentMover::apply_fragment(
	Frame const & frame,
	Size frag_num,
	kinematics::MoveMap const & movemap,
	pose::Pose & pose
) const {
	return frame.apply( movemap, frag_num, pose );
}

/// @brief DONT ALLOW HELICES OF LESS THAN 3 OR STRANDS OF LESS THAN 2
/// Fix this: inserting length two helices at the chain end is allowed
/// as is inserting length 1 strands... For the moment, preserving r++'s
/// incorrect behavior
bool ClassicFragmentMover::valid_ss( std::string const & new_ss ) const {
	bool valid = true;
	Size helix_len = 0;
	Size strand_len = 0;
	for ( Size i = 1; i < new_ss.size(); ++i ) {
		if ( new_ss[i-1] == 'H' ) ++helix_len;
		if ( new_ss[i-1] == 'E' ) ++strand_len;
		if ( new_ss[i-1] != new_ss[i] ) {
			if ( helix_len != 0 && helix_len < 3 ) {
				tr.Trace << "short helix_len" << std::endl;
				valid = false;
				break;
			}
			if ( strand_len != 0 && strand_len < 2 ) {
				tr.Trace << "short strand_len" << std::endl;
				valid = false;
				break;
			}
			helix_len = 0;
			strand_len = 0;
		} // if ( new_ss[i] != new_ss[i+1] )
	} // for ( Size i = 1; i < new_ss.length(); ++i )
	return valid;
}


/// @brief choose and insert a Fragment from the protocols::moves::Movers Fragment-Set into a Pose.
void ClassicFragmentMover::apply( core::pose::Pose & pose ) {
	PROF_START( basic::FRAGMENT_MOVER );
	//	update_insert_map( ); // checks if bValidInsertMap == false

	// If the insert map is empty dont attempt fragment insertions
	// to avoid corrupting the pose, memory and/or your grandmother.
	if( insert_size_.size() == 0 ){
		return;
	}

	// find a fragment
	bool success ( false );

	// since we have an OP of the movemap it might be changed from the outside:
	// this means the insertmap has to be regenerated... do this at most once
	bool insert_map_definitely_right( false );

	Size nfail = 0;
	while ( !success && ( nfail < 100 ) ) {
		Size frag_begin;
		Size window_length;
		FrameList frames;
		while ( nfail < 100 && frames.size() == 0 ) {
			// choose a fragment length
			if ( !choose_window_length( pose, window_length ) ) {
				nfail++;
				continue;
			}


			// choose an insertion point
			if( use_predefined_window_start_ )
				frag_begin = predefined_window_start_;
			else if ( !choose_window_start( pose, window_length, frag_begin ) ) {
				nfail++;
				continue;
			}

			// retrieve fragments for this position
			if ( !fragset_->region( *movemap(), frag_begin, frag_begin + window_length - 1,
					min_overlap_, min_frag_length_,	frames ) ) {
				// if we are here, we couldn't find fragments at this position --- maybe recompute insert_map
				if ( !insert_map_definitely_right ) {
					insert_map_definitely_right = true;
					tr.Debug << "didn't find fragment at predicted position ==> update_insert_map() " << std::endl;
					update_insert_map();
				} else {
					utility_exit_with_message(" couldn't find fragments --- inconsistency in the insert map " );
				}
				nfail++;
				continue;
			} // if fragset_->region
		} // while ( frames.size() == 0 )
			// if ssblock or random_frag only a subset of fragments was there to choose from, this decision could be made
			// outside of this code by supplying a Subset of the frag_set. Alternatively, we could do weight-based sampling
			// like the original design of mini fragments and just set some weights to zero.

		if ( frames.size() ) { // we got fragments
			success = apply_frames( pose, frames );
			//poss. reason for failure: ss-check, jump not present in fold-tree
		}
	} // while ( !success );

	if ( !success ) {
		tr.Error << "couldn't find fragment to insert!!" << std::endl;
		return;
	}
	PROF_STOP( basic::FRAGMENT_MOVER );

} // apply

bool ClassicFragmentMover::apply_frames( pose::Pose &pose, FrameList const& frames ) const {
	Size frame_num;
	Size frag_num;
	bool success( false );
	if ( !choose_fragment( frames, pose, frame_num /*output*/, frag_num /*output*/ ) ) return false;
	if ( tr.Trace.visible() ) tr.Trace
															<< "frag (" << frames[ frame_num ]->start() << ","
															<< frag_num << ","
															<< frames[ frame_num ]->nr_res_affected( *movemap_ )
															<< ")" << std::endl;
	if ( !check_ss() ) return apply_fragment( *frames[ frame_num ], frag_num, *movemap_, pose );

	// now do the ss-check!
	//	tr.Trace << "now do the ss-check!"<< std::endl;
	// get actual ss from pose
	std::string proposed_ss;
	proposed_ss.reserve( pose.total_residue() );
	proposed_ss = pose.secstruct();

	std::string old_ss = proposed_ss;

	// check if old ss is valid
	bool valid = !valid_ss( old_ss ); // if old_ss is not valid we can apply the fragment anyway

	// if old ss was fine ---> check fragments effect on ss
	if ( !valid ) { // if old_ss was valid we check if proposed_ss is still valid.
		frames[ frame_num ]->apply_ss( *movemap_, frag_num, proposed_ss );
		//		tr.Trace << !valid << " old_ss: " << old_ss << std::endl;
		valid = valid_ss( proposed_ss );
		//		tr.Trace << valid << "new_ss: " << proposed_ss << std::endl;
	}
	//	tr.Trace << "finished the ss-check! : " << valid << std::endl;
	if ( valid ) {
		success = apply_fragment( *frames[ frame_num ], frag_num, *movemap_, pose );
	} else {
		//		tr.Trace << "dissallow insertion due to short helix/strand " << std::endl;
	}
	return success;
}

std::ostream &operator<< ( std::ostream &os, ClassicFragmentMover const &cfmover )
{
	cfmover.show(os);
	return os;
}


LoggedFragmentMover::LoggedFragmentMover(
core::fragment::FragSetCOP fragset,
	core::kinematics::MoveMapCOP movemap
)	: ClassicFragmentMover( fragset, movemap, "ClassicFragmentMover" )
{}

LoggedFragmentMover::~LoggedFragmentMover()
{}

std::string
LoggedFragmentMover::get_name() const {
	return "LoggedFragmentMover";
}

bool
LoggedFragmentMover::apply_fragment(
		core::fragment::Frame const& frame,
		Size frag_num,
		core::kinematics::MoveMap const& movemap,
		core::pose::Pose &pose
) const {
	bool success = Parent::apply_fragment( frame, frag_num, movemap, pose );
	if ( success ) {
		logs_.push_back( Item( frame.start(), frag_num ) );
	}
	return success;
}

void
LoggedFragmentMover::show( std::ostream& out ) const {
	using namespace ObjexxFCL::format;
	for ( Storage::const_iterator it=logs_.begin(), eit=logs_.end(); it!=eit; ++it) {
		out << RJ(5, it->frame_pos) << ' ' << RJ(5, it->frag_num) << std::endl;
	}
}

void
LoggedFragmentMover::clear() {
	logs_.clear();
}

} // simple_moves
} // protocols
