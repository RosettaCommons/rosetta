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
/// @author

#include <core/types.hh>
#include <core/id/SequenceMapping.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

#include <algorithm>

#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace id {

using namespace ObjexxFCL;

/// @brief ctor
SequenceMapping::SequenceMapping( utility::vector1< Size > const & mapping ):
	size2_( mapping.size() ),
	mapping_( mapping )
{}

SequenceMapping::SequenceMapping(
	conformation::signals::LengthEvent const & event
)
: size2_( 0 )
{

	//int direction(0);
	bool longer( event.length_change > 0 );
	Size upstream_res(0), old_num_res(0);

	mapping_.clear();

	if ( event.tag == conformation::signals::LengthEvent::INVALIDATE ) {
		//utility_exit_with_message("invalidated LengthEvent passed");
		//i guess in this case it's better to return an empty sequence mapping and let
		//the observers deal with it
		return;
	} else if ( event.tag == conformation::signals::LengthEvent::RESIDUE_APPEND ) {
		runtime_assert( longer || event.length_change == 0);
		//direction = 1;
		upstream_res = event.position;
	} else if ( event.tag == conformation::signals::LengthEvent::RESIDUE_PREPEND ) {
		runtime_assert( longer || event.length_change == 0);
		//direction = 1;
		upstream_res = event.position - 1;
	} else if ( event.tag == conformation::signals::LengthEvent::RESIDUE_DELETE ) {
		runtime_assert( !longer || event.length_change == 0);
		//direction = -1;
		upstream_res = event.position - 1;
	} else {
		utility_exit_with_message(
			"unknown signal triggered by conformation length change. please update this file"
		);
	}

	runtime_assert( event.conformation_size != 0 );

	old_num_res = event.conformation_size - event.length_change;

	for ( Size i = 1; i <= upstream_res; ++i ) mapping_.push_back( i );

	//if ( direction == 1 ) mapping_.push_back( upstream_res + 1 + event.length_change );
	if ( longer ) {
		if ( upstream_res < old_num_res ) {
			mapping_.push_back( upstream_res + 1 + event.length_change );
		}
	} else {
		for ( int i = event.length_change; i < 0; ++i ) mapping_.push_back( 0 );
	}
	//for ( Size i = upstream_res + 2; i <= old_num_res; ++i ) mapping_.push_back( i + direction );
	Size downstream_res( mapping_.size() + 1 );
	runtime_assert( downstream_res <= old_num_res + 1);
	for ( Size i = downstream_res; i <= old_num_res; ++i ) mapping_.push_back( i + event.length_change );
}

SequenceMapping::~SequenceMapping() = default;


SequenceMapping::SequenceMapping( SequenceMapping const & src )
: ReferenceCount(src)
{
	*this = src;
}

SequenceMapping &
SequenceMapping::operator = ( SequenceMapping const & src ) {
	size2_      = src.size2();
	mapping_    = src.mapping();

	return *this;
}

void
SequenceMapping::resize( Size const s1, Size const s2 )
{
	size2_ = s2;
	mapping_.clear();
	mapping_.resize(s1,0);
}


void
SequenceMapping::reverse()
{
	// update size2!
	size2_ = *max_element( mapping_.begin(), mapping_.end() );

	utility::vector1< Size > new_mapping( size2_, 0 );
	for ( Size i = 1; i <= size1(); ++i ) {
		if ( mapping_[i] != 0 ) {
			new_mapping[ mapping_[i] ] = i;
		}
	}

	size2_ = mapping_.size();
	mapping_.swap( new_mapping );
}

void
SequenceMapping::downstream_combine( core::id::SequenceMapping const & smap_to_add )
{
	for ( core::Size i = 1; i <= mapping_.size(); ++i ) {
		if ( mapping_[i] != 0 ) {
			mapping_[i] = smap_to_add[ mapping_[i] ];
		}
	}
	size2_ = smap_to_add.size2();
}

void
SequenceMapping::upstream_combine( core::id::SequenceMapping const & smap_to_add )
{
	utility::vector1< Size > new_mapping( smap_to_add.mapping() );
	for ( core::Size i = 1; i <= new_mapping.size(); ++i ) {
		if ( new_mapping[i] != 0 && new_mapping[i] <= mapping_.size() ) {
			new_mapping[i] = mapping_[ new_mapping[i] ];
		} else {
			new_mapping[i] = 0;
		}
	}
	mapping_.swap( new_mapping );
	// size2_ stays the same
}

/// access
Size
SequenceMapping::size1() const
{
	return mapping_.size();
}

Size
SequenceMapping::size2() const
{
	return size2_;
}


bool
SequenceMapping::all_aligned() const
{
	bool aligned( true );
	for ( Size i=1; i<= size1(); ++i ) {
		if ( mapping_[i] == 0 ) {
			aligned = false;
			break;
		}
	}
	return aligned;
}


bool
SequenceMapping::is_identity() const {
	bool identity( true );
	for ( Size i=1; i<= size1(); ++i ) {
		if ( mapping_[i] != i ) {
			identity = false;
			break;
		}
	}
	return identity;
}

bool
SequenceMapping::is_identity_ignore_gaps() const {
	bool identity( true );
	for ( Size i = 1; i <= size1(); ++i ) {
		if ( mapping_[i] != 0 && mapping_[i] != i ) {
			identity = false;
			break;
		}
	}
	return identity;

}

void
SequenceMapping::size2( Size const s2 )
{
	size2_ = s2;
}

void
SequenceMapping::push_back( Size const al )
{
	mapping_.push_back( al );
}

void
SequenceMapping::delete_source_residue( Size const pos1 )
{
	mapping_.erase( mapping_.begin() + pos1-1 );
}

void
SequenceMapping::show() const
{
	show( basic::T("id.SequenceMapping") );
}

void
SequenceMapping::show( std::ostream & output ) const {
	for ( Size i=1; i<= size1(); ++i ) {
		output << ("id.SequenceMapping ") << i << " --> ";
		if ( mapping_[i] ) output << mapping_[i] << std::endl;
		else output << "----" << std::endl;
	}
}


void
SequenceMapping::set_offset( int setting ) {
	for ( Size i=1; i<= size1(); ++i ) {
		if ( mapping_[ i ] ) {
			int pos =  static_cast< int >(mapping_[ i ]) - setting;
			mapping_[ i ] = ( pos > 0 ) ? pos : 0;
		}
	}
}


void
SequenceMapping::insert_source_residue( Size const pos1 )
{
	mapping_.insert( mapping_.begin() + pos1-1, 0 );
}

void
SequenceMapping::insert_aligned_residue( Size const pos1, Size const pos2 )
{
	mapping_.insert( mapping_.begin() + pos1-1, pos2 );
}

// same as insert_aligned_residue, but a couple of extra checks on size1 and size2.
void
SequenceMapping::insert_aligned_residue_safe(
	Size const pos1,
	Size const pos2
) {
	if ( pos1 == 0 ) return;
	size2( std::max( pos2, size2() ) );
	if ( pos1 >= size1() ) {
		mapping_.resize( pos1, 0 );
	}

	mapping_[ pos1 ] = pos2;
}

void
SequenceMapping::insert_target_residue( Size const pos )
{
	++size2_;
	for ( Size i=1; i<= size1(); ++i ) {
		if ( mapping_[i] >= pos ) ++mapping_[i];
	}
}

void
SequenceMapping::delete_target_residue( Size const pos )
{
	--size2_;
	for ( Size i=1; i<= size1(); ++i ) {
		if ( mapping_[i] == pos ) mapping_[i] = 0;
		else if ( mapping_[i] > pos ) --mapping_[i];
	}
}

void
SequenceMapping::clear()
{
	mapping_.clear();
	size2_ = 0;
}

Size
SequenceMapping::operator[]( Size const pos1 ) const {
	if ( pos1 > mapping_.size() ) return 0;
	return mapping_[ pos1 ];
}

Size &
SequenceMapping::operator[]( Size const pos1 ) {
	if ( pos1 > mapping_.size() ) mapping_.resize( pos1, 0 );
	return mapping_[ pos1 ];
}

///@brief Get the corresponding (new) resnum from the old.
Size
SequenceMapping::get_corresponding_residue_in_current( Size original_resnum ) const {
	return mapping_[ original_resnum ];
}

//utility::vector1< Size >
//SequenceMapping::get_old_resnums() const {
// utility::vector1 < Size > mapped_residues;
// for (auto const kv : mapping_){
//  mapped_residues.push_back(kv.first);
// }
// return mapped_residues;
//}

bool SequenceMapping::operator == ( SequenceMapping const & rhs ) const
{
	if ( ! same_type_as_me( rhs ) || ! rhs.same_type_as_me( *this ) ) return false;
	return size2_ == rhs.size2_ && mapping_ == rhs.mapping_;
}

bool SequenceMapping::same_type_as_me( SequenceMapping const & other ) const
{
	SequenceMapping const * casted_other = dynamic_cast< SequenceMapping const *  > (&other);
	return casted_other;
}

std::string SequenceMapping::to_string() const {
	std::string retval( "" );
	for ( Size i = 1; i <= size1(); ++i ) {
		retval += string_of(i) + " --> ";
		if ( mapping_[i] ) {
			retval += string_of(mapping_[i]) + "\n";
		} else {
			retval += static_cast< std::string > ("----") + "\n";
		}
	}
	return retval;
}

/// @details combine all input sequence mappings into one.
/// sequentially, that is
core::id::SequenceMappingOP
combine_sequence_mappings(
	utility::vector1< core::id::SequenceMapping > const & smaps
){

	using namespace core::id;

	//gigo :)
	if ( smaps.size() == 0 ) return core::id::SequenceMappingOP( new SequenceMapping() ) ;

	SequenceMappingOP composite_smap( new SequenceMapping() );
	*composite_smap = smaps[1];

	for ( core::Size i = 2; i <= smaps.size(); ++i ) {
		composite_smap->downstream_combine( smaps[i] );
	}

	return composite_smap;

} //combine_sequence_mappings


/// @details combine smap_to_add into smap,
/// i.e. smap[j] becomes smap_to_add[ smap[j] ]
void
combine_sequence_mappings(
	core::id::SequenceMapping & smap,
	core::id::SequenceMapping const & smap_to_add )
{
	smap.downstream_combine( smap_to_add );
}

} // id
} // core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::SequenceMapping::save( Archive & arc ) const {
	arc( CEREAL_NVP( size2_ ) ); // Size
	arc( CEREAL_NVP( mapping_ ) ); // utility::vector1<Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::SequenceMapping::load( Archive & arc ) {
	arc( size2_ ); // Size
	arc( mapping_ ); // utility::vector1<Size>
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::SequenceMapping );
CEREAL_REGISTER_TYPE( core::id::SequenceMapping )

#endif // SERIALIZATION

#ifdef    SERIALIZATION
CEREAL_REGISTER_DYNAMIC_INIT( core_id_SequenceMapping )
#endif // SERIALIZATION
