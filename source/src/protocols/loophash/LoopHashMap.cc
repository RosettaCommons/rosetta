// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/LoopHashMap.cc
/// @brief
/// @author Mike Tyka


// libRosetta headers
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/BackboneDB.hh>
#include <protocols/loophash/Exceptions.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallProvider.hh>

#include <boost/unordered_map.hpp>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

#include <numeric/HomogeneousTransform.hh>
//#include <protocols/match/Hit.fwd.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/pointer/owning_ptr.hh>

#include <numeric/angle.functions.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>

// C++ headers
#include <cstdio>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/lexical_cast.hpp>


#include <protocols/moves/Mover.fwd.hh>
#include <utility/vector1.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


using namespace protocols::moves;
using namespace core::scoring;
using namespace core;
using namespace core::pose;
using namespace conformation;
using namespace kinematics;
using namespace protocols::frag_picker;


namespace protocols {
namespace loophash {


static THREAD_LOCAL basic::Tracer TR( "LoopHashMap" );


/// @brief This takes a pose and two residue positions and determines the rigid body transform of the Leap described by those two residues.
///        Returns true is successful or false if something went haywire and one should just ignore this loop (this can happen at the ends)
bool get_rt_over_leap( const core::pose::Pose& orig_pose, core::Size ir, core::Size jr, numeric::geometry::hashing::Real6 &rt_6 ){
	using namespace core;
	using namespace core::pose;
	using namespace conformation;
	using namespace kinematics;
	using namespace numeric::geometry::hashing;

	core::pose::Pose pose = orig_pose;

	//fpd vrt/ligand trim
	core::Size nres = pose.total_residue();
	while ( !pose.residue_type(nres).is_polymer() ) nres--;

	// get current cutpoints; don't try to connect these
	utility::vector1< Size > cuts_in = pose.fold_tree().cutpoints();
	std::sort( cuts_in.begin(), cuts_in.end() );

	// bail if (ir,jr) crosses a cut
	for ( Size i=1; i<=cuts_in.size(); ++i ) {
		if ( cuts_in[i]<=jr && cuts_in[i]>=ir ) {
			TR.Error << "ERROR -- residue range crosses cut    IR: " << ir << "  JR: " << jr << "  CUT: " << cuts_in[i] << std::endl;
			return false;
		}
		//fpd insertions one position after the cut seem not to work ...
		//fpd perhaps if the foldtree for the local segment were reversed this might be ok
		if ( cuts_in[i]==ir-1 ) {
			TR.Error << "ERROR -- startres immediately follows cut    IR: " << ir << "  CUT: " << cuts_in[i] << std::endl;
			return false;
		}
	}

	// Create a fake foldtree with a jump from ir to jr, and the cutpoint just before jr. From that extract the
	// rigid body transfer. This is fairly hacky, but actually works v reliably and is easy to understand, without
	// having to mess deeply with the fold tree.
	// FoldTree f;
	// Size cutpoint= jr-1;
	// f.add_edge( 1, ir, Edge::PEPTIDE );
	// f.add_edge( ir, cutpoint, Edge::PEPTIDE );
	// f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
	// f.add_edge( jr, nres , Edge::PEPTIDE );
	// f.add_edge( ir, jr, 1 );  // this is the jump !!

	//fpd handle multiple chains/chainbreaks
	FoldTree f;
	core::Size last_cut=0, jump_num=2;
	Size cutpoint= jr-1;
	for ( Size i=1; i<=cuts_in.size(); ++i ) {
		if ( cuts_in[i] >= nres ) break;
		if ( cutpoint > last_cut && cutpoint < cuts_in[i] ) {
			f.add_edge( last_cut+1, ir, Edge::PEPTIDE );
			f.add_edge( ir, cutpoint, Edge::PEPTIDE );
			f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
			f.add_edge( jr, cuts_in[i] , Edge::PEPTIDE );
			f.add_edge( ir, jr, 1 );  // this is the jump !!
			if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
		} else {
			f.add_edge( last_cut+1, cuts_in[i], Edge::PEPTIDE );
			if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
		}
		last_cut = cuts_in[i];
	}
	if ( last_cut+1 <= nres ) {
		if ( cutpoint > last_cut && cutpoint < nres ) {
			f.add_edge( last_cut+1, ir, Edge::PEPTIDE );
			f.add_edge( ir, cutpoint, Edge::PEPTIDE );
			f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
			f.add_edge( jr, nres , Edge::PEPTIDE );
			f.add_edge( ir, jr, 1 );  // this is the jump !!
			if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
		} else {
			f.add_edge( last_cut+1, nres, Edge::PEPTIDE );
			if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
		}
	}
	for ( core::Size i=nres+1; i<=pose.total_residue(); ++i ) {
		f.add_edge( 1, i, jump_num++ );  // additional jumps
	}

	core::Size theroot = 1;
	if ( ir == 1 ) theroot = pose.total_residue();
	if ( orig_pose.residue_type( orig_pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) theroot = orig_pose.fold_tree().root();  //fpd
	if ( f.reorder(theroot) == false ) {
		TR.Error << "ERROR During reordering of fold tree - am ignoring this LOOP ! bailing: The root: " << theroot << " NRES " << pose.total_residue() << "   IR: " << ir << "  JR: " << jr << std::endl;
		return false; // continuing leads to a segfault - instead ignore this loop !
	}

	// Apply this new foldtree to the pose.
	pose.fold_tree( f );

	// Now extract the rigid body transform!
	Jump myjump;
	myjump = pose.jump( 1 );

	// Aha, now you have the RT (RigidbodyTransform)
	RT rt = myjump.rt();

	// Create a 6 value representation: (just change the data format)
	numeric::HomogeneousTransform< Real > ht( rt.get_rotation() , rt.get_translation() );
	numeric::xyzVector < Real > euler_angles =  ht.euler_angles_rad();

	rt_6[1] = rt.get_translation().x();
	rt_6[2] = rt.get_translation().y();
	rt_6[3] = rt.get_translation().z();
	rt_6[4] = euler_angles.x()*180.0/numeric::constants::d::pi;
	rt_6[5] = euler_angles.y()*180.0/numeric::constants::d::pi;
	rt_6[6] = euler_angles.z()*180.0/numeric::constants::d::pi;

	// indicate success
	return true;
}


/// @brief This takes a pose and two residue positions and determines the rigid body transform of the Leap described by those two residues
///       THe difference between this and the get_rt_over_leap function is that this version doesnt make a copy of the pose which makes it faster.
///     However this means that the pose passed cannot be a const pose, even though the function restores the fold tree afterwards..
bool get_rt_over_leap_fast( core::pose::Pose& pose, core::Size ir, core::Size jr, numeric::geometry::hashing::Real6 &rt_6 ){
	using namespace core;
	using namespace core::pose;
	using namespace conformation;
	using namespace kinematics;
	using namespace numeric::geometry::hashing;

	core::Size newroot=0;
	if ( pose.residue_type( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) newroot = pose.fold_tree().root();

	//fpd vrt/ligand trim
	core::Size nres = pose.total_residue();
	while ( !pose.residue_type(nres).is_polymer() ) nres--;

	// get current cutpoints; don't try to connect these
	utility::vector1< Size > cuts_in = pose.fold_tree().cutpoints();
	std::sort( cuts_in.begin(), cuts_in.end() );

	// bail if (ir,jr) crosses a cut
	for ( Size i=1; i<=cuts_in.size(); ++i ) {
		if ( cuts_in[i]<=jr && cuts_in[i]>=ir ) {
			TR.Error << "ERROR -- residue range crosses cut    IR: " << ir << "  JR: " << jr << "  CUT: " << cuts_in[i] << std::endl;
			return false;
		}
		//fpd insertions one position after the cut seem not to work ...
		//fpd perhaps if the foldtree for the local segment were reversed this might be ok
		if ( cuts_in[i]==ir-1 ) {
			TR.Error << "ERROR -- startres immediately follows cut    IR: " << ir << "  CUT: " << cuts_in[i] << std::endl;
			return false;
		}
	}

	// Create a fake foldtree with a jump from ir to jr, and the cutpoint just before jr. From that extract the
	// rigid body transfer. This is fairly hacky, but actually works v reliably and is easy to understand, without
	// having to mess deeply with the fold tree.
	// FoldTree f;
	// Size cutpoint= jr-1;
	// f.add_edge( 1, ir, Edge::PEPTIDE );
	// f.add_edge( ir, cutpoint, Edge::PEPTIDE );
	// f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
	// f.add_edge( jr, nres , Edge::PEPTIDE );
	// f.add_edge( ir, jr, 1 );  // this is the jump !!

	//fpd handle multiple chains/chainbreaks
	FoldTree f, f_orig=pose.fold_tree();
	core::Size last_cut=0, jump_num=2;
	Size cutpoint= jr-1;
	for ( Size i=1; i<=cuts_in.size(); ++i ) {
		if ( cuts_in[i] >= nres ) break;
		if ( cutpoint > last_cut && cutpoint < cuts_in[i] ) {
			f.add_edge( last_cut+1, ir, Edge::PEPTIDE );
			f.add_edge( ir, cutpoint, Edge::PEPTIDE );
			f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
			f.add_edge( jr, cuts_in[i] , Edge::PEPTIDE );
			f.add_edge( ir, jr, 1 );  // this is the jump !!
			if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
		} else {
			f.add_edge( last_cut+1, cuts_in[i], Edge::PEPTIDE );
			if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
		}
		last_cut = cuts_in[i];
	}
	if ( last_cut+1 <= nres ) {
		if ( cutpoint > last_cut && cutpoint < nres ) {
			f.add_edge( last_cut+1, ir, Edge::PEPTIDE );
			f.add_edge( ir, cutpoint, Edge::PEPTIDE );
			f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
			f.add_edge( jr, nres , Edge::PEPTIDE );
			f.add_edge( ir, jr, 1 );  // this is the jump !!
			if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
		} else {
			f.add_edge( last_cut+1, nres, Edge::PEPTIDE );
			if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
		}
	}
	for ( core::Size i=nres+1; i<=pose.total_residue(); ++i ) {
		f.add_edge( 1, i, jump_num++ );  // additional jumps
	}

	core::Size theroot = 1;
	if ( ir == 1 ) theroot = pose.total_residue();
	if ( newroot>0 ) theroot = newroot;  //fpd
	if ( f.reorder(theroot) == false ) {
		TR.Error << "ERROR During reordering of fold tree - am ignoring this LOOP ! bailing: The root: " << theroot << " NRES " << pose.total_residue() << "   IR: " << ir << "  JR: " << jr << std::endl;
		return false; // continuing leads to a segfault - instead ignore this loop !
	}

	// Apply this new foldtree to the pose.
	pose.fold_tree( f );

	// Now extract the rigid body transform!
	Jump myjump;
	myjump = pose.jump( 1 );
	RT rt = myjump.rt();

	// restore prev foldtree
	pose.fold_tree( f_orig );

	// TR.Info << "R: " << ir << " " << jr << rt << std::endl;

	// Create a 6 value representation:

	numeric::HomogeneousTransform< Real > ht( rt.get_rotation() , rt.get_translation() );
	numeric::xyzVector < Real > euler_angles =  ht.euler_angles_rad();


	rt_6[1] = rt.get_translation().x();
	rt_6[2] = rt.get_translation().y();
	rt_6[3] = rt.get_translation().z();
	rt_6[4] = euler_angles.x()*180.0/numeric::constants::d::pi;
	rt_6[5] = euler_angles.y()*180.0/numeric::constants::d::pi;
	rt_6[6] = euler_angles.z()*180.0/numeric::constants::d::pi;

	return true;
}

bool
get_rt_over_leap_without_foldtree_bs(
	core::pose::Pose const & pose,
	core::Size ir,
	core::Size jr,
	numeric::geometry::hashing::Real6 &rt_6
){
	using core::id::AtomID;

	if ( !pose.residue(ir).is_protein() ) return false;
	if ( !pose.residue(jr).is_protein() ) return false;

	core::id::StubID id1( AtomID(2,ir), AtomID(1,ir), AtomID(3,ir-1) );
	core::id::StubID id2;
	if ( pose.residue(jr).type().name1() == 'P' ) { //jumps to proline use Carbon as jump atom
		id2 = core::id::StubID( AtomID(2,jr), AtomID(1,jr), AtomID(3,jr) );
	} else {
		id2 = core::id::StubID( AtomID(2,jr), AtomID(1,jr), AtomID( pose.residue(jr).atom_index("H") ,jr) );
	}

	RT rt = pose.conformation().get_stub_transform(id1,id2);

	numeric::HomogeneousTransform< Real > ht( rt.get_rotation() , rt.get_translation() );
	numeric::xyzVector < Real > euler_angles =  ht.euler_angles_rad();

	rt_6[1] = rt.get_translation().x();
	rt_6[2] = rt.get_translation().y();
	rt_6[3] = rt.get_translation().z();
	rt_6[4] = euler_angles.x()*180.0/numeric::constants::d::pi;
	rt_6[5] = euler_angles.y()*180.0/numeric::constants::d::pi;
	rt_6[6] = euler_angles.z()*180.0/numeric::constants::d::pi;

	return true;
}


LoopHashMap::LoopHashMap( core::Size loop_size){
	loop_size_ = loop_size;
	setup(loop_size);
}

LoopHashMap::LoopHashMap(LoopHashMap const & ) = default;

LoopHashMap & LoopHashMap::operator=(LoopHashMap const & other)
{
	if ( this != & other ) {
		hash_ = other.hash_;
		backbone_index_map_ = other.backbone_index_map_;
		loopdb_ = other.loopdb_;
		loop_size_ = other.loop_size_;
	}
	return *this;

}

void LoopHashMap::mem_foot_print(){
	TR << "loopdb_: " << loopdb_.size() << " Size: " << loopdb_.size() * sizeof( LeapIndex ) << std::endl;
	TR << "BackboneIndexMap: " << backbone_index_map_.size() << " Size: " << backbone_index_map_.size() * (sizeof(boost::uint64_t) + sizeof(core::Size) ) << std::endl;
}

void LoopHashMap::sort() {
	std::sort( loopdb_.begin(), loopdb_.end(), by_index() );
}

void
LoopHashMap::setup( core::Size loop_size)
{
	loop_size_ = loop_size;

	if ( loop_size > 1 ) {
		// only show this if this class is being set up with a meaningful parameter. Loop_size == 1 means it was called without an argument to the constructor and
		// so its just an intermediate data structure.
		TR.Info << "Setting up hash_: Size:  " << loop_size << std::endl;

	}

	BoundingBox bounding_box( core::Vector(
		-HASH_POSITION_GRID_BASE*(int)loop_size,
		-HASH_POSITION_GRID_BASE*(int)loop_size,
		-HASH_POSITION_GRID_BASE*(int)loop_size
		),
		core::Vector(
		HASH_POSITION_GRID_BASE*(int)loop_size,
		HASH_POSITION_GRID_BASE*(int)loop_size,
		HASH_POSITION_GRID_BASE*(int)loop_size
		) );
	numeric::geometry::hashing::Size3 euler_offsets;
	euler_offsets[1] = 0;
	euler_offsets[2] = 0;
	euler_offsets[3] = 0;
	numeric::geometry::hashing::Real6 bin_widths;

	core::Real space_multiplier = 0.2;
	core::Real angle_multiplier =  15.0/6.0;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ lh::grid_space_multiplier].user() ) {
		space_multiplier = option[ lh::grid_space_multiplier]();
	}
	if ( option[ lh::grid_angle_multiplier].user() ) {
		angle_multiplier = option[ lh::grid_angle_multiplier]();
	}

	bin_widths[1] = space_multiplier*loop_size;
	bin_widths[2] = space_multiplier*loop_size;
	bin_widths[3] = space_multiplier*loop_size;
	bin_widths[4] = angle_multiplier*loop_size;
	bin_widths[5] = angle_multiplier*loop_size;
	bin_widths[6] = angle_multiplier*loop_size;

	hash_ = numeric::geometry::hashing::SixDCoordinateBinnerOP( new numeric::geometry::hashing::SixDCoordinateBinner( bounding_box, euler_offsets, bin_widths ) );
	// initialize the radial tree
	hash_->tree_init(option[lh::max_radius]());
}

void
LoopHashMap::add_legacyleap( const LegacyLeapIndex & legacyleap_index)
{
	numeric::geometry::hashing::Real6 rt_6;
	rt_6[1] = legacyleap_index.vecx;
	rt_6[2] = legacyleap_index.vecy;
	rt_6[3] = legacyleap_index.vecz;
	rt_6[4] = legacyleap_index.rotx;
	rt_6[5] = legacyleap_index.roty;
	rt_6[6] = legacyleap_index.rotz;
	LeapIndex leap_index;
	leap_index.index = 0;
	leap_index.offset = legacyleap_index.ba;
	add_leap( leap_index, rt_6 );
}

void LoopHashMap::add_leap( const LeapIndex &leap_index, boost::uint64_t key ) {
	core::Size cpindex = loopdb_.size();
	loopdb_.push_back( leap_index );
	// No check to see if it is in hash because this data is already processed
	backbone_index_map_.insert( std::make_pair( key, cpindex ));
}

void LoopHashMap::add_leap( const LeapIndex &leap_index, numeric::geometry::hashing::Real6 & transform ){
	core::Size cpindex = loopdb_.size();
	loopdb_.push_back( leap_index );
	numeric::geometry::hashing::Real6 rt_6;
	rt_6[1] = transform[1];
	rt_6[2] = transform[2];
	rt_6[3] = transform[3];
	rt_6[4] = transform[4];
	rt_6[5] = transform[5];
	rt_6[6] = transform[6];
	while (  rt_6[ 4 ] < 0.0 )    rt_6[ 4 ] += 360.0;
	while (  rt_6[ 4 ] >= 360.0 ) rt_6[ 4 ] -= 360.0;
	while (  rt_6[ 5 ] < 0.0 )    rt_6[ 5 ] += 360.0;
	while (  rt_6[ 5 ] >= 360.0 ) rt_6[ 5 ] -= 360.0;
	while (  rt_6[ 6 ] < 0.0 )  rt_6[ 6 ] += 180.0;
	while (  rt_6[ 6 ] > 180.0 ) rt_6[ 6 ] -= 180.0;
	if ( !( rt_6[ 4 ] >= 0.0 && rt_6[ 4 ] < 360.0 ) )  { std::cerr << "rt_6[4] out of bounds: " << rt_6[4] << std::endl;  utility_exit_with_message( "RANGE ERROR" ); }
	if ( !( rt_6[ 5 ] >= 0.0 && rt_6[ 5 ] < 360.0 ) )  { std::cerr << "rt_6[5] out of bounds: " << rt_6[5] << std::endl;  utility_exit_with_message( "RANGE ERROR" ); }
	if ( !( rt_6[ 6 ] >= 0.0 && rt_6[ 6 ] <= 180.0 ) ) { std::cerr << "rt_6[6] out of bounds: " << rt_6[6] << std::endl;  utility_exit_with_message( "RANGE ERROR" ); }
	hash_index( rt_6, cpindex  );
}

void LoopHashMap::hash_index( numeric::geometry::hashing::Real6 transform, core::Size data ){
	// Hash the transform
	if ( ! hash_->contains( transform ) ) {
		//TR.Error << "OutofBOIUNDS! " << std::endl;
		return;
	}

	boost::uint64_t bin_index = hash_->bin_index( transform );
	// update previously added leapindex with key
	loopdb_[data].key = bin_index;
	// Add Address to map with hashed index
	backbone_index_map_.insert( std::make_pair( bin_index, data ) );
}


void LoopHashMap::bbdb_range( std::pair< BackboneIndexMap::iterator, BackboneIndexMap::iterator > & range ){
	range.first = backbone_index_map_.begin();
	range.second = backbone_index_map_.end();
}

// append to a bucket of vectors in the appropriate bin
void
LoopHashMap::lookup(  numeric::geometry::hashing::Real6 transform, std::vector < core::Size > &result )
{
	// Hash the transform
	if ( ! hash_->contains( transform ) ) {
		//TR.Error << "OutofBOIUNDS! " << std::endl;
		return;
	}

	transform[4] = numeric::nonnegative_principal_angle_degrees(transform[4] );
	transform[5] = numeric::nonnegative_principal_angle_degrees(transform[5] );
	boost::uint64_t bin_index = hash_->bin_index( transform );

	// now get an iterator over that map entry

	TR.Info << "backbone bucket size:  " << backbone_index_map_.bucket_count() << std::endl;
	TR.Info << "backbone_index_map_.find(bin_index):  " << backbone_index_map_.count(bin_index) << std::endl;
	std::pair<  BackboneIndexMap::iterator,
		BackboneIndexMap::iterator> range = backbone_index_map_.equal_range( bin_index );

	for ( BackboneIndexMap::iterator it = range.first;
			it != range.second;
			++it ) {
		TR.Info << it->second << std::endl;
		TR.Info << "bucket_size:  " << result.size() << std::endl;
		result.push_back( it->second );
	}
}

void LoopHashMap::radial_lookup( core::Size radius,  numeric::geometry::hashing::Real6 center, std::vector < core::Size > &result )
{

	center[4] = numeric::nonnegative_principal_angle_degrees(center[4] );
	center[5] = numeric::nonnegative_principal_angle_degrees(center[5] );
	//TR.Info << "center:  " << center[4] << " " << center[5] << std::endl;
	std::vector< boost::uint64_t > bin_index_vec = hash_->radial_bin_index( radius, center );

	for ( auto & i : bin_index_vec) {
		//TR.Info << "bin_index_vec[i]:  " << bin_index_vec[i] << std::endl;
		// now get an iterator over that map entry
		std::pair< BackboneIndexMap::iterator,BackboneIndexMap::iterator> range = backbone_index_map_.equal_range( i );
		for ( BackboneIndexMap::iterator it = range.first; it != range.second; ++it ) {
			//TR.Info << "result:  " << result << std::endl;
			result.push_back( it->second );
		}
	}
}

Size LoopHashMap::radial_count( core::Size radius, numeric::geometry::hashing::Real6 center ) const
{
	center[4] = numeric::nonnegative_principal_angle_degrees(center[4] );
	center[5] = numeric::nonnegative_principal_angle_degrees(center[5] );
	std::vector< boost::uint64_t > bin_index_vec = hash_->radial_bin_index( radius, center );
	Size count = 0;
	for (auto & i : bin_index_vec) {
		count += backbone_index_map_.count( i );
	}
	return count;
}

void
LoopHashMap::lookup(
	core::Size index,
	std::vector < core::Size > &result
)
{
	if ( index > backbone_index_map_.bucket_count() ) {
		//TR.Error << "OutofBOIUNDS! " << std::endl;
		return;
	}

	BackboneIndexMap::local_iterator begin = backbone_index_map_.begin( index );
	BackboneIndexMap::local_iterator end = backbone_index_map_.end( index );

	for ( BackboneIndexMap::local_iterator it = begin; it != end; ++it ) {
		result.push_back( it->second );
	}
}

boost::uint64_t
LoopHashMap::return_key( core::Size bb_index )
{
	LeapIndex leap_index = get_peptide( bb_index );
	return leap_index.key;
}

void
LoopHashMap::radial_lookup_withkey( boost::uint64_t key, core::Size radius, std::vector < core::Size > &result )
{
	// get center of bin from key
	numeric::geometry::hashing::Real6 center = hash_->bin_center_point( hash_->bin_from_index( key ) );
	radial_lookup( radius, center, result );
}

void
LoopHashMap::lookup_withkey( boost::uint64_t key, std::vector < core::Size > &result )
{
	std::pair<  BackboneIndexMap::iterator,
		BackboneIndexMap::iterator> range = backbone_index_map_.equal_range( key );

	for ( BackboneIndexMap::iterator it = range.first; it != range.second; ++it ) {
		result.push_back( it->second );
	}
	return;
}

void LoopHashMap::read_legacydb(std::string filename )
{
	// use basic C input - C++ are too memory hungry to deal with these potentially v large files
	FILE *file = fopen( filename.c_str(), "r" );
	if ( file == nullptr ) throw EXCN_DB_IO_Failed( filename, "read" );

	loopdb_.clear();
	while ( !feof( file ) ) {
		const unsigned int bufsize = 128;
		LegacyLeapIndex bufferdata[bufsize];
		size_t readshorts = fread(&bufferdata[0],sizeof(LegacyLeapIndex),bufsize,file);
		for ( unsigned i = 0; i< readshorts; i ++ ) {
			add_legacyleap( bufferdata[i] );
		}
	}
	fclose( file );
	loopdb_.reserve( loopdb_.size() );

	TR.Info << "Rehashed dataset:" << backbone_index_map_.size() << "  " <<  loopdb_.size() << std::endl;
	TR.Info << "Rehashed dataset:  OUTOFBOUNDS: " <<  loopdb_.size() - backbone_index_map_.size() << std::endl;
}

void LoopHashMap::write_db( std::string filename ){
	std::ofstream file( filename.c_str() );
	if ( !file ) throw EXCN_DB_IO_Failed( filename, "write" );
	for (auto & i : loopdb_) {
		file << i.index << " " << i.offset << " " << i.key << std::endl;
	}
	file.close();
}

void
LoopHashMap::read_db(
	std::string filename, std::pair< core::Size, core::Size > loopdb_range,
	std::map< core::Size, bool > & homolog_index
){
	std::ifstream file( filename.c_str() );
	if ( !file ) throw EXCN_DB_IO_Failed( filename, "read" );
	std::string line;
	LeapIndex leap_index;
	while ( getline( file, line ) ) {
		std::vector<std::string> output;
		std::string t; std::istringstream sline(line);
		while ( sline >> t ) { output.push_back(t); }
		if ( output.size() != 3 ) throw EXCN_Wrong_DB_Format( filename );
		leap_index.index = boost::lexical_cast< core::Size > ( output[0] );
		if ( leap_index.index < loopdb_range.first ) continue;
		if ( leap_index.index >= loopdb_range.second && loopdb_range.second != 0 ) continue;
		// Since we're reading in partitions of the bbdb
		// We also have to offset where the leap_index.index points
		leap_index.index -= loopdb_range.first;
		// if the leap_index is in the homolog map, skip it
		if ( homolog_index.find( leap_index.index ) != homolog_index.end() ) continue;
		leap_index.offset = boost::lexical_cast< core::Size > ( output[1] );
		leap_index.key  = boost::lexical_cast< boost::uint64_t > ( output[2] );
		add_leap( leap_index, leap_index.key );
	}
	TR.Info << "Loophashmap range " << loopdb_[0].index << " " << loopdb_[loopdb_.size() - 1].index << std::endl;
	file.close();
}


} // namespace loops
} // namespace protocols
