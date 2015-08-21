// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragData.cc
/// @brief  a collection classes of the FragData and SingleResidueFragData class hirachy
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


// Unit Headers
#include <core/fragment/FragData.hh>

// Package Headers
#include <core/fragment/Frame.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {

using namespace kinematics;


//helper function
void
make_pose_from_sequence_(
	std::string sequence,
	chemical::ResidueTypeSet & residue_set,
	pose::Pose& pose
) {
	using namespace chemical;
	// clear all of the old data in the pose
	pose.clear();

	// setup the pose by appending the appropriate residues residues
	for ( Size seqpos = 1; seqpos <= sequence.length(); ++seqpos ) {
		char aa = sequence[seqpos-1]; // string indexing is zero-based!
		AA my_aa = aa_from_oneletter_code( aa );
		//ResidueTypeCOPs const & rsd_type_list( residue_set.aa_map( my_aa ) );
		//Size best_index = 1;
		//ResidueType const & rsd_type( *(rsd_type_list[ best_index ]) );
		ResidueType const & rsd_type( * residue_set.get_representative_type_aa( my_aa ) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );
	} // for seqpos
	// pose.conformation().insert_chain_ending( pose.total_residue() - 1 );   probably not necessary
} // make_pose_match_sequence_

FragData::FragData( SingleResidueFragDataOP SRFD, Size n) : valid_( false ), score_( 0.0 ) {
	for ( Size i = 1; i<=n; i++ ) {
		data_.push_back( SRFD->clone() );
	}
}


FragDataOP FragData::clone() const {
	FragDataOP fd( new FragData( size() ) );
	for ( Size pos = 1; pos<=size(); pos++ ) {
		fd->data_[pos] = data_[pos]->clone();
	};
	return fd;
}

void FragData::copy(FragData const& frag_data)
{
	data_ = SRFD_List(frag_data.data_.size());
	for ( Size pos = 1; pos<=frag_data.data_.size(); pos++ ) {
		data_[pos] = frag_data.data_[pos]->clone();
	};

	valid_ = frag_data.valid_;
	score_ = frag_data.score_;
}

Size FragData::apply( MoveMap const& mm, pose::Pose& pose, Size start, Size end ) const {
	if ( !is_valid() ) return 0;
	Size pos = start;
	Size ct( 0 );
	for ( SRFD_List::const_iterator it= data_.begin(), eit=data_.end(); it!=eit; ++it, ++pos ) {
		//success = success && data_[ pos-start+1 ]->apply( pose, pos );
		if ( pos > end ) break;
		if ( pos > pose.total_residue() ) break;
		if ( !(*it)->is_applicable( mm, pos ) ) continue; //don't apply for this residue
		if ( !(*it)->apply( mm, pose, pos ) ) continue; //don't apply for this residue
		ct++;
	}
	return ct;
}

Size FragData::apply( pose::Pose& pose, Size start, Size end ) const {
	Size pos = start;
	Size ct ( 0 );
	if ( !is_valid() ) return 0;
	for ( SRFD_List::const_iterator it= data_.begin(), eit=data_.end(); it!=eit; ++it, ++pos ) {
		if ( pos > end ) break;
		if ( !(*it)->apply( pose, pos ) ) continue;
		++ct;
	}
	return ct;
}

Size FragData::apply( MoveMap const& mm, pose::Pose& pose, Frame const& frame ) const {
	runtime_assert( size() == frame.length() );
	if ( !is_valid() ) return 0;
	if ( frame.is_continuous() ) {
		return apply( mm, pose, frame.start(), frame.end() );
	}
	Size ct ( 0 );
	Size ipos( 1 );
	for ( SRFD_List::const_iterator it= data_.begin(), eit=data_.end(); it!=eit; ++it, ++ipos ) {
		if ( frame.seqpos( ipos ) > pose.total_residue() ) continue;
		if ( !(*it)->is_applicable( mm, ipos, frame ) ) continue;
		if ( !(*it)->apply( mm, pose, ipos, frame ) ) continue;
		++ct;
	}
	return ct;
}

Size FragData::apply( pose::Pose& pose, Frame const& frame ) const {
	runtime_assert( size() == frame.length() );
	if ( !is_valid() ) return 0;
	Size ipos( 1 );
	Size ct ( 0 );
	for ( SRFD_List::const_iterator it= data_.begin(), eit=data_.end(); it!=eit; ++it, ++ipos ) {
		if ( frame.seqpos( ipos ) > pose.total_residue() ) continue;
		if ( !(*it)->apply( pose, ipos, frame ) ) continue;
		++ct;
	}
	return ct;
}

Size FragData::apply_ss( MoveMap const& mm, std::string& ss, Frame const& frame ) const {
	Size ipos ( 1 );
	Size ct( 0 );
	if ( !is_valid() ) return 0;
	for ( SRFD_List::const_iterator it= data_.begin(), eit=data_.end(); it!=eit; ++it, ++ipos ) {
		if ( !(*it)->is_applicable( mm, ipos, frame ) ) continue;
		if ( !(*it)->apply_ss( ss, ipos, frame ) ) continue;
		++ct;
	}
	return ct;
}


Size FragData::is_applicable( kinematics::MoveMap const& mm, Size start, Size end ) const {
	Size insert_size( 0 );
	// if ( !is_valid() ) return 0; //this might be necessary to have!!!
	for ( Size pos=start; pos<=end; pos++ ) {
		if ( !data_[ pos-start+1 ]->is_applicable( mm, pos ) ) continue;
		++insert_size;
	}
	return insert_size;
}

Size FragData::is_applicable( kinematics::MoveMap const& mm, Frame const& frame ) const {
	runtime_assert( size() == frame.length() );
	// if ( !is_valid() ) return 0;
	if ( frame.is_continuous() ) {
		return is_applicable( mm, frame.start(), frame.end() );
	}
	Size insert_size( 0 );
	for ( Size j=1; j<=size(); j++ ) {
		if ( !data_[ j ]->is_applicable( mm, j, frame ) ) continue;
		++insert_size;
	}
	return insert_size;
}

bool FragData::steal( pose::Pose const& pose, Size start, Size end ) {
	bool success( true );
	for ( Size pos=start; pos<=end; pos++ ) {
		runtime_assert( start > 0 );
		if ( pos > pose.total_residue() ) return false;
		success = success && data_[ pos-start+1 ]->steal( pose, pos );
	}
	if ( success ) set_valid();
	return success;
}


bool FragData::steal( pose::Pose const& pose, Frame const& frame ) {
	runtime_assert( size() == frame.length() );
	if ( frame.is_continuous() ) {
		return steal( pose, frame.start(), frame.end() );
	}
	bool success( true );
	for ( Size j=1; j<=size(); j++ ) {
		runtime_assert( frame.seqpos( j ) > 0 );
		if ( frame.seqpos( j ) > pose.total_residue() ) return false;
		success = success && data_[ j ]->steal( pose, j, frame );
	}
	if ( success ) set_valid();
	return success;
}

FragDataOP FragData::generate_sub_fragment( Size start, Size stop ) const {
	runtime_assert( stop >= start );
	runtime_assert( stop <= size() );
	FragDataOP new_frag( new FragData );
	for ( Size pos = start; pos<=stop; pos++ ) {
		new_frag->add_residue( data_[ pos ] ); //reuse of data
	}
	new_frag->set_valid();
	return new_frag;
}

bool FragData::is_compatible( FragData const& frag_data ) const {
	if ( frag_data.size() != size() ) return false;
	for ( SRFD_List::const_iterator it1= data_.begin(), it2=frag_data.data_.begin(), eit1=data_.end(); it1!=eit1; ++it1,++it2 ) {
		//  if ( typeid( *it1 ) != typeid( *it2 ) ) return false;
		if ( ! (*it1)->is_compatible( **it2 ) ) return false;
	}
	return true;
}

void FragData::show( std::ostream &os, Frame const& frame ) const {
	Size i = 1;
	if ( is_valid() ) {
		for ( SRFD_List::const_iterator it= data_.begin(), eit=data_.end(); it!=eit; ++it ) {
			os << ObjexxFCL::format::RJ( 10, frame.seqpos( i++ ) ) << " " << ObjexxFCL::format::RJ( 5, pdbpos() ) << " " << pdbid() <<" ";
			//   std::cerr << "FragData::show " << i-1 << std::endl;
			(*it)->show( os );
			os << std::endl;
			//   std::cerr << "FragData:show -- done" << std::endl;
		}
	} else {
		os << "EMPTY TEMPLATE" << std::endl;
	}
}

void FragData::show( std::ostream &os ) const {
	Size i = 1;
	if ( is_valid() ) {
		for ( SRFD_List::const_iterator it= data_.begin(), eit=data_.end(); it!=eit; ++it ) {
			os << ObjexxFCL::format::RJ( 10, i++ ) << " " << ObjexxFCL::format::RJ( 5, pdbpos() ) << " " << pdbid() << " " ;
			(*it)->show( os );
			os << std::endl;
		}
	} else {
		os << "EMPTY_TEMPLATE" << std::endl;
	}
}

void FragData::show_classic( std::ostream &os ) const {
	for ( SRFD_List::const_iterator it= data_.begin(), eit=data_.end(); it!=eit; ++it ) {
		os << " 1xxx X   111 ";
		(*it)->show( os );
		os << std::endl;
	}
}

FragDataOP AnnotatedFragData::clone() const {
	return FragDataOP( new AnnotatedFragData( pdbid_, startpos_, *Parent::clone() ) );
}

void AnnotatedFragData::copy(FragData const& frag_data)
{
	FragData::copy(frag_data);
	pdbid_ = frag_data.pdbid();
	chain_ = frag_data.chain();
	startpos_ = frag_data.pdbpos();
}

FragDataOP AnnotatedFragData::generate_sub_fragment( Size start, Size stop ) const {
	runtime_assert( stop >= start );
	runtime_assert( stop <= size() );
	FragDataOP new_frag( new AnnotatedFragData(pdbid_, startpos_ + start - 1, chain_) );

	for ( Size pos = start; pos<=stop; pos++ ) {
		new_frag->add_residue( data_[ pos] ); //reuse of data
	}
	new_frag->set_valid();
	return new_frag;
}


} // fragment
} // core
