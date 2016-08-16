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
/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <protocols/jumping/JumpSample.hh>

// Package Headers
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/util.hh>
#include <core/scoring/dssp/PairingsList.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/JumpingFrame.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/JumpSRFD.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/dssp/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <fstream>

#include <core/chemical/VariantType.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/util.hh>
#include <protocols/jumping/JumpSample.hh>
#include <utility/vector1.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace jumping {

using namespace core;
using namespace ObjexxFCL::format;

ChainbreakDistFunc::ChainbreakDistFunc( Real const x0_in ) : d2target_( x0_in*x0_in ) {}
ChainbreakDistFunc::~ChainbreakDistFunc() {}

core::scoring::func::FuncOP
ChainbreakDistFunc::clone() const { return core::scoring::func::FuncOP( new ChainbreakDistFunc( *this ) ); }

bool
ChainbreakDistFunc::operator == ( core::scoring::func::Func const & other ) const {
	if ( !same_type_as_me( other ) || !other.same_type_as_me(*this) ) return false;
	ChainbreakDistFunc const & other_downcast( static_cast< ChainbreakDistFunc const & > (other) );
	return d2target_ == other_downcast.d2target_;
}

bool
ChainbreakDistFunc::same_type_as_me( core::scoring::func::Func const & other ) const {
	return dynamic_cast< ChainbreakDistFunc const * > (&other);
}


Real ChainbreakDistFunc::func( Real const x ) const {
	return std::sqrt( std::abs( x*x - d2target_ ) );
}

Real ChainbreakDistFunc::dfunc( Real const x ) const {
	return std::sqrt( std::abs( x*x - d2target_ ) );
}



using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer tr( "protocols.jumping" );

JumpSample::JumpSample( JumpSetup const& def ) :
	total_residue_( def.total_residue() ),
	njump_( def.size() ),
	bValidTree_( false )
{
	resize( def.size() );
	int ct = 1;
	for ( JumpSetup::const_iterator it=def.begin(), eit=def.end(); it!=eit; ++it, ++ct ) {
		jumps_( 1, ct ) = it->jump_.start_;
		jumps_( 2, ct ) = it->jump_.end_;
		jump_atoms_( 1, ct ) = "";
		jump_atoms_( 2, ct ) = "";
		Size const crs ( it->cut_reg_.start_ );
		Size const cre ( it->cut_reg_.end_ );
		cuts_( ct ) = crs+int( numeric::random::uniform()*( cre-crs ) + 0.5 );
	}
	generate_tree_from_jumps_and_cuts();
	jumps2pairings();
	if ( !bValidTree_ && tr.Debug.visible() ) {
		tr.Debug << "invalid FoldTree in JumpSample: " << (*this) << std::endl;
	}
}

JumpSample::JumpSample ( Size total_residue, Size njump, FArray2D_int jumps, FArray1D_int cuts, Size root ) :
	total_residue_( total_residue ),
	njump_( njump ),
	jumps_ ( jumps ),
	cuts_ ( cuts ),
	bValidTree_( false )
{
	jump_atoms_.redimension(2, njump, "");
	generate_tree_from_jumps_and_cuts(root);
	jumps2pairings();
}

JumpSample::JumpSample ( Size total_residue, Size njump, FArray2D_int jumps, FArray2D<std::string> jump_atoms, FArray1D_int cuts, Size root) :
	total_residue_( total_residue ),
	njump_( njump ),
	jumps_ ( jumps ),
	jump_atoms_ ( jump_atoms ),
	cuts_ ( cuts ),
	bValidTree_( false )
{
	generate_tree_from_jumps_and_cuts(root);
	jumps2pairings();
}

JumpSample::JumpSample ( Size total_residue, core::scoring::dssp::PairingsList const& jump_pairings, core::fragment::SecondaryStructure const& ssdef, Size root ) :
	total_residue_( total_residue ),
	jump_pairings_( jump_pairings ),
	bValidTree_( false )
{
	resize( jump_pairings.size() );

	// first fill jumps_ array:
	njump_ = 0;
	for ( core::scoring::dssp::PairingsList::const_iterator it = jump_pairings.begin(),
			eit = jump_pairings.end(); it != eit; ++it ) {
		njump_++;
		jumps_(1,njump_) = it->Pos1();
		jumps_(2,njump_) = it->Pos2();
		jump_atoms_(1,njump_) = "";
		jump_atoms_(2,njump_) = "";
	}
	runtime_assert( njump_ == jump_pairings.size() ); //kind of trivial but something odd was going on
	if ( ssdef.total_residue() <= total_residue ) {
		core::fragment::SecondaryStructure new_ssdef = ssdef;
		new_ssdef.extend( total_residue );
		generate_random_tree_from_jumps( new_ssdef.loop_fraction(), root );
	} else {
		generate_random_tree_from_jumps( ssdef.loop_fraction(), root );
	}
}

JumpSample::JumpSample ( Size total_residue, core::scoring::dssp::PairingsList const& jump_pairings, FArray1D_float const& cut_probability, Size root) :
	total_residue_( total_residue ),
	jump_pairings_( jump_pairings ),
	bValidTree_( false )
{
	resize( jump_pairings.size() );

	// first fill jumps_ array:
	njump_ = 0;
	for ( core::scoring::dssp::PairingsList::const_iterator it = jump_pairings.begin(),
			eit = jump_pairings.end(); it != eit; ++it ) {
		njump_++;
		jumps_(1,njump_) = it->Pos1();
		jumps_(2,njump_) = it->Pos2();
		jump_atoms_(1,njump_) = "";
		jump_atoms_(2,njump_) = "";
	}
	runtime_assert( njump_ == jump_pairings.size() ); //kind of trivial but something odd was going on
	generate_random_tree_from_jumps( cut_probability, root );
}

/// we want N-CA-C as the stub atoms for the jump
/// to achieve that we have to set the upstream jump atom
/// according to folding direction N2C -->  JumpAtom C
///                                C2N --> JumpAtom N
/// since the AtomTree automatically uses the parent and grand-parent as second and third stub atom
void
JumpSample::correct_jump_atoms_for_fragments() const {
	//this method could also live in the FoldTree
	using namespace kinematics;
	if ( !bValidTree_ ) return;
	runtime_assert( fold_tree_ != 0 );
	for ( Size i = 1; i <= njump_; ++i ) {
		fold_tree_->set_jump_atoms( i, jump_atoms_(1,i), jump_atoms_(2,i) );
	}
	fold_tree_->put_jump_stubs_intra_residue();
}

void
JumpSample::generate_tree_from_jumps_and_cuts( Size root) {
	if ( total_residue_ == 0 ) total_residue_ = 2500; //don't make it too big.. it alloates an FArray1bool in fold-tree with this size
	fold_tree_ = core::kinematics::FoldTreeOP( new kinematics::FoldTree );
	bValidTree_ = fold_tree_->tree_from_jumps_and_cuts( total_residue_, njump_, jumps_, cuts_, root, true /* verbose */ );
	if ( bValidTree_ ) correct_jump_atoms_for_fragments();
}

void
JumpSample::jumps2pairings() {
	jump_pairings_.clear();
	for ( Size jump_nr = 1; jump_nr <= size(); jump_nr++ ) {
		jump_pairings_.push_back( core::scoring::dssp::Pairing( jumps_(1, jump_nr ), jumps_(2, jump_nr), 0, 0 ) ); //initialized with 0 = unknown orientation/pleating
	}
}

void
JumpSample::generate_random_tree_from_jumps( FArray1D_float const& prob, Size root ) {
	if ( total_residue_ == 0 ) total_residue_ = 2500;
	fold_tree_ = core::kinematics::FoldTreeOP( new kinematics::FoldTree );
	bValidTree_ = false;
	Size attempts( 10 );
	while ( !bValidTree_ && attempts-- > 0 )  {
		bValidTree_ =
			fold_tree_->random_tree_from_jump_points( total_residue_, njump_, jumps_, prob, root, true /*yes we allow 1 nres jumps*/ );
		if ( bValidTree_ ) {
			cuts_.dimension( njump_ );
			tr.Debug << "cut points ";// << std::endl;
			for ( Size i = 1; i <= njump_; i++ ) {
				cuts_( i ) = fold_tree_->cutpoint( i );
				//don't except trees that cut at jump point
				if ( fold_tree_->is_jump_point( cuts_( i ) ) ) {
					bValidTree_ = false;
					continue;
				}
				tr.Debug << cuts_( i ) << " ";
			}
			tr.Debug << std::endl;
		} // tree was connected
	} // while
	if ( bValidTree_ ) correct_jump_atoms_for_fragments();
}

JumpSample::JumpSample( kinematics::FoldTree const& f ) {
	apply_to( core::kinematics::FoldTreeOP( new kinematics::FoldTree( f ) ) );
	correct_jump_atoms_for_fragments();
}

JumpSample::JumpSample( kinematics::FoldTreeOP f ) {
	apply_to( f );
	correct_jump_atoms_for_fragments();
}


void JumpSample::apply_to ( kinematics::FoldTreeOP pf ) {
	fold_tree_ = pf;
	kinematics::FoldTree &f( *fold_tree_ );
	bValidTree_ = f.check_fold_tree();
	total_residue_ =  f.nres();
	if ( bValidTree_ ) {
		resize( f.num_jump() );

		// get jumps
		for ( Size ct = 1; ct <= njump_; ct++ ) {
			jumps_(1,ct) = f.upstream_jump_residue( ct );
			jumps_(2,ct) = f.downstream_jump_residue( ct );
		}

		// get cuts
		for ( Size i = 1; i <= njump_; i++ ) {
			cuts_( i ) = f.cutpoint( i );
		}

		jumps2pairings();

	} // bValidTree_
}

void
JumpSample::resize( Size njump ) {
	jumps_.redimension( 2, njump );
	jump_atoms_.redimension(2, njump);
	cuts_.redimension( njump );
	jump_pairings_.resize( njump );
	njump_=njump;
}

void
JumpSample::set_fold_tree_in_pose( pose::Pose &pose ) const {
	runtime_assert( fold_tree_ != 0 );
	runtime_assert( bValidTree_ != 0 );
	pose.fold_tree( *fold_tree_ );
}

void
JumpSample::safe_secstruct( pose::Pose &pose ) const {
	runtime_assert( fold_tree_ != 0 );
	runtime_assert( total_residue_ == pose.total_residue() );
	runtime_assert( *fold_tree_ == pose.fold_tree() );

	Size const num_jump ( size() );
	Size const nres( pose.total_residue() );

	for ( Size i = 1; i <= num_jump; ++i ) {
		for ( Size j = 1; j <= 2; ++j ) {
			Size const pos = jumps_(j,i);
			char const ss( pose.secstruct( pos ) );
			if ( ss != 'L' ) {
				for ( Size k = std::max( (Size) 1, pos-2 ), ke = std::min( nres, pos+2 );
						k <= ke; ++k ) {
					if ( pose.secstruct(k) != ss ) {
						pose.set_secstruct( k, ss );
					}
				}
			}
		}
	}
}

//@brief transfer native jump RT to poes, Sideeffect: changes fold-tree of native_pose.
void
JumpSample::transfer_jumps( core::pose::Pose &pose, core::pose::Pose &native_pose ) const {
	runtime_assert( pose.fold_tree() == *fold_tree_ );
	for ( Size i=1; i<=size(); i++ ) {
		kinematics::Jump aJump = native_pose.jump( i );
		tr.Info << "transfer jump " << i <<" " <<  aJump <<std::endl;
		pose.set_jump( i, aJump );
	}
}

//@brief transfer native jump RT to poes, Sideeffect: changes fold-tree of pose.
void
JumpSample::steal_jumps(
	core::pose::Pose &pose,
	core::fragment::FrameIterator const& begin,
	core::fragment::FrameIterator const& end
) const {
	using namespace core::fragment;
	set_fold_tree_in_pose( pose );
	for ( FrameIterator it = begin, eit = end;
			it!=eit; ++it ) {
		(*it)->steal( pose );
	}
}


//@brief generate list of frames ( one for each jump ) and steal RT from pose
void
JumpSample::generate_jump_frames(
	core::fragment::FrameList& all_frames,
	kinematics::MoveMap const& mm,
	bool bWithTorsion
) const
{
	all_frames.reserve( all_frames.size() + size() );// as many new frames as nr of jumps
	// find out how many different kind of fragments are we interested in:
	// max of four: A 1 , A 2, P 1, P 2

	for ( Size jump_nr = 1; jump_nr <= size(); jump_nr++ ) {
		int const startpos( jumps_( 1, jump_nr ) );
		int const endpos( jumps_( 2, jump_nr ) );

		tr.Debug << "jump " << jump_nr << " from: " << startpos << " --> " << endpos << std::endl;
		tr.Debug << "downstream res " <<  fold_tree().downstream_jump_residue( jump_nr )
			<< " upstream "    <<   fold_tree().upstream_jump_residue( jump_nr ) << std::endl;
		// I assume that jump_nr never change, runtime_assert that this is indeed the case
		int const down_jump_res( fold_tree().downstream_jump_residue( jump_nr ) );
		int const up_jump_res( fold_tree().upstream_jump_residue( jump_nr ) );
		runtime_assert( (down_jump_res == startpos) || (up_jump_res == startpos) );
		runtime_assert( (down_jump_res == endpos  ) || (up_jump_res == endpos  ) );
		runtime_assert( startpos != endpos );

		if ( mm.get_bb( up_jump_res ) || mm.get_bb( down_jump_res ) ) {
			using namespace core::fragment;
			FragDataOP frag_data( new FragData );

			if ( bWithTorsion && mm.get_bb( up_jump_res ) ) {
				BBTorsionSRFDOP start( new BBTorsionSRFD( 3, 'E', 'X' ) );
				frag_data->add_residue( start );
			}
			frag_data->add_residue( SingleResidueFragDataOP( new UpJumpSRFD ) );
			frag_data->add_residue( SingleResidueFragDataOP( new DownJumpSRFD ) );

			if ( bWithTorsion && mm.get_bb( down_jump_res ) ) {
				BBTorsionSRFDOP stop( new BBTorsionSRFD( 3, 'E', 'X' ) );
				frag_data->add_residue( stop );
			}

			JumpingFrameOP frame( new JumpingFrame( startpos, endpos, frag_data->size() ) );
			Size pos = 1;
			if ( bWithTorsion && mm.get_bb( up_jump_res ) ) frame->set_pos( pos++, startpos );
			frame->set_pos( pos++, startpos );
			frame->set_pos( pos++, endpos );
			if ( bWithTorsion && mm.get_bb( down_jump_res ) ) frame->set_pos( pos++, endpos );
			frame->add_fragment( frag_data );
			runtime_assert( frame->nr_frags() ); //adding has worked?
			all_frames.push_back( frame );
		}
	} // for Jumps iteration
} // method


void
JumpSample::steal_orientation_and_pleating( core::pose::Pose &native_pose ) {
	for ( Size jump_nr=1; jump_nr<=size(); jump_nr++ ) {
		// get native orientation to select the correct jump-geometries
		core::scoring::dssp::Pairing & p = jump_pairings_[ jump_nr ];
		tr.Info << "detect orientation and pleating for jump " << jump_nr <<" from " << p.Pos1() << " to " << p.Pos2() << std::endl;
		//compute_orientation_and_pleating
		core::Size orientation, pleating;
		core::scoring::dssp::get_pleating( native_pose, p.Pos1(), p.Pos2(), orientation, pleating );
		p.Orientation(orientation);
		p.Pleating(pleating);
		tr.Info << "orientation is " << p.Orientation() << " pleating is " << p.Pleating() << std::endl;
	}
}

bool JumpSample::has_orientation_and_pleating() const {
	return core::scoring::dssp::has_orientation_and_pleating( jump_pairings_ );
}

core::scoring::dssp::Pairing
JumpSample::get_pairing( Size res1, Size res2 ) const {
	runtime_assert ( has_orientation_and_pleating() );
	for ( core::scoring::dssp::PairingsList::const_iterator it = jump_pairings_.begin(), eit = jump_pairings_.end();
			it != eit; ++it ) {
		if ( it->Pos1() == res1 && it->Pos2() == res2 ) return *it;
		if ( it->Pos1() == res2 && it->Pos2() == res1 ) return it->generate_reversed();
	}
	return core::scoring::dssp::Pairing( 0, 0, 0, 0 );
}

//@brief generate fragset with RTs from library, take orientation and pleating from native_pose
void
JumpSample::generate_jump_frags(
	PairingLibrary const& lib,
	kinematics::MoveMap const& mm,
	bool bWithTorsion,
	core::fragment::FrameList& all_frames
) const {
	using namespace core::fragment;

	all_frames.reserve( all_frames.size() + size() );// as many new frames as nr of jumps

	// find out how many different kind of fragments are we interested in:
	// max of four: A 1 , A 2, P 1, P 2
	runtime_assert( has_orientation_and_pleating() );
	typedef utility::vector1< Size > JumpList;
	typedef std::map< std::pair< Size, Size >, JumpList > JumpOrientations;
	JumpOrientations jump_kind;
	Size jump_nr ( 1 );
	for ( core::scoring::dssp::PairingsList::const_iterator it = jump_pairings_.begin(), eit = jump_pairings_.end();
			it != eit; ++it ) {
		Size o_key ( it->Orientation() ); // < 0 ? 1 : 2 );
		Size p_key ( it->Pleating() ); // < 0 ? 1 : 2 );
		jump_kind[ std::make_pair( o_key, p_key ) ].push_back( jump_nr++ );
	}

	// now generate fragments for each of the maximum four JumpOrientations present
	for ( JumpOrientations::const_iterator it=jump_kind.begin(), eit=jump_kind.end();
			it!=eit;
			++it ) {
		Size o_key( it->first.first ); //orientation
		Size p_key( it->first.second ); //pleating ... believe me or not, it is in first.second
		FragDataOPs frag_data;
		lib.create_jump_fragments( o_key, p_key, bWithTorsion, frag_data );
		for ( JumpList::const_iterator jit=it->second.begin(), ejit=it->second.end();
				jit!=ejit; ++jit ) {
			int const jump_nr ( *jit );
			int const startpos( jumps_( 1, jump_nr ) );
			int const endpos( jumps_( 2, jump_nr ) );

			runtime_assert( o_key == jump_pairings_[ jump_nr ].Orientation() );
			runtime_assert( p_key == jump_pairings_[ jump_nr ].Pleating() );

			tr.Debug << "jump " << jump_nr << " from: " << startpos << " --> " << endpos << std::endl;
			tr.Debug << "downstream res " <<  fold_tree().downstream_jump_residue( jump_nr )
				<< " upstream "    <<   fold_tree().upstream_jump_residue( jump_nr ) << std::endl;

			// I assume that jump_nr never change, runtime_assert that this is indeed the case
			int const down_jump_res( fold_tree().downstream_jump_residue( jump_nr ) );
			int const up_jump_res( fold_tree().upstream_jump_residue( jump_nr ) );
			runtime_assert( (down_jump_res == startpos) || (up_jump_res == startpos) );
			runtime_assert( (down_jump_res == endpos  ) || (up_jump_res == endpos  ) );
			runtime_assert( startpos != endpos );
			if ( mm.get_bb( up_jump_res ) && mm.get_bb( down_jump_res ) ) {
				Size const length( bWithTorsion ? 4 : 2 );
				runtime_assert( length == frag_data.front()->size() );
				JumpingFrameOP frame = generate_empty_jump_frame( up_jump_res, down_jump_res, length );
				frame->add_fragment( frag_data );
				all_frames.push_back( frame );
			} else {
				utility_exit_with_message("need to implement this: make ss-library fragments that only contain those torsions for residues "
					"that can move according to movemap -- call this function with "
					"bWithTorsions = false ... and it works for now");
			}
		} // for JumpList iteration
	} // loop over orientations and pleatings
} // method

void
JumpSample::add_chainbreaks( pose::Pose &pose ) const {
	for ( Size i = 1; i<= njump_; i++ ) {
		if ( pose.residue_type( cuts_(i) ).has_variant_type( chemical::UPPER_TERMINUS_VARIANT ) ) continue;
		if ( pose.residue_type( cuts_(i)+1 ).has_variant_type( chemical::LOWER_TERMINUS_VARIANT ) ) continue;
		tr.Debug << "add chainbreak variant to residues " << cuts_(i) << " and " << cuts_(i)+1 << std::endl;
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cuts_(i) );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cuts_(i)+1 );
	}
}

void
JumpSample::add_chainbreaks( pose::Pose &pose, Size max_dist, core::kinematics::ShortestPathInFoldTree const& sp) const {
	//remove_chainbreaks( pose ); not necessary if max_dist is monotonoically increaseing
	for ( Size i = 1; i<= njump_; i++ ) {
		if ( sp.dist( cuts_(i), cuts_(i)+1 ) <= max_dist
				&& pose.residue( cuts_(i) ).is_polymer() && pose.residue( cuts_(i)+1 ).is_polymer() ) {
			if ( pose.residue_type( cuts_(i) ).has_variant_type( chemical::UPPER_TERMINUS_VARIANT ) ) continue;
			if ( pose.residue_type( cuts_(i)+1 ).has_variant_type( chemical::LOWER_TERMINUS_VARIANT ) ) continue;
			tr.Debug << "add chainbreak variant to residues " << cuts_(i) << " and " << cuts_(i)+1 << std::endl;
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cuts_(i) );
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cuts_(i)+1 );
		}
	}
}

using namespace scoring::constraints;
void
JumpSample::add_chainbreaks_as_distance_constraint( pose::Pose &pose ) const {
	for ( Size i = 1; i<= njump_; i++ ) {
		core::scoring::func::FuncOP f( new ChainbreakDistFunc( 1.7424 ) );
		pose.add_constraint(
			scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
			id::AtomID( pose.residue(cuts_(i)    ).atom_index("C"), cuts_(i)    ),
			id::AtomID( pose.residue(cuts_(i) + 1).atom_index("N"), cuts_(i) + 1),
			f
			) ) )
		);
	}
}

void
JumpSample::add_chainbreaks_as_distance_constraint(
	pose::Pose &pose,
	Size max_dist,
	core::kinematics::ShortestPathInFoldTree const& sp
) const {
	//remove_chainbreaks( pose ); not necessary if max_dist is monotonoically increaseing
	for ( Size i = 1; i<= njump_; i++ ) {
		if ( sp.dist( cuts_(i), cuts_(i)+1 ) <= max_dist ) {
			tr.Debug << "add chainbreak as distance constraint to residues " << cuts_(i) << " and " << cuts_(i)+1 << std::endl;
			core::scoring::func::FuncOP f( new ChainbreakDistFunc( 1.7424 ) );
			pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
				id::AtomID( pose.residue(cuts_(i)    ).atom_index("C"), cuts_(i)    ),
				id::AtomID( pose.residue(cuts_(i) + 1).atom_index("N"), cuts_(i) + 1),
				f
				) ) ) );
		}
	}
}

void
JumpSample::remove_chainbreaks( pose::Pose &pose ) const {
	for ( Size i = 1; i<= njump_; i++ ) {
		core::pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_LOWER, cuts_(i) );
		core::pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_UPPER, cuts_(i)+1 );
	}
}

std::ostream & operator <<(std::ostream & os, JumpSample const & t) {
	for ( Size i=1; i<=t.size(); i++ ) {
		os << RJ(4,t.jumps_(1,i)) << " " << RJ(4,t.jumps_(2,i)) << " | ";
	}
	os << "cuts:";
	for ( Size i=1; i<=t.size(); i++ ) {
		os << " " << RJ(4,t.cuts_(i));
	}
	return os;
}

void
JumpSample::dump_pymol( std::string fn ) const {
	utility::io::ozstream out( fn );
	if ( out ) {
		out << "from pymol import cmd\n";
		for ( Size i=1; i<=size(); i++ ) {
			out << "cmd.select( \"jumps"<<i<<"\", \"resi "<<jumps_(1,i)<<"+"<<jumps_(2,i)<<"\");" << std::endl;
			out << "cmd.show( \"sticks\", \"jumps"<<i<<"\");\n";
			out << "cmd.select( \"cut"<<i<<"\", \"resi " <<cuts_(i)<<"\");\n";
			out << "cmd.show( \"sphere\",\"cut"<<i<<"\");" << std::endl;
		}
	}
}

} //jumping
} //protocols

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::jumping::ChainbreakDistFunc::ChainbreakDistFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jumping::ChainbreakDistFunc::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::func::Func >( this ) );
	arc( CEREAL_NVP( d2target_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jumping::ChainbreakDistFunc::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::func::Func >( this ) );
	arc( d2target_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jumping::ChainbreakDistFunc );
CEREAL_REGISTER_TYPE( protocols::jumping::ChainbreakDistFunc )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jumping_JumpSample )
#endif // SERIALIZATION
