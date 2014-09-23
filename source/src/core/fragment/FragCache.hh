// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/Frame.hh
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_FragCache_HH
#define INCLUDED_core_fragment_FragCache_HH

// Unit Headers
#include <core/fragment/FragCache.fwd.hh>

// Package Headers
// AUTO-REMOVED #include <core/fragment/FragData.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/BaseCacheUnit.hh>
#include <core/fragment/FragID_Iterator.hh>


// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/excn/Exceptions.hh>

// C++ STL Headers
#include <map>

#include <utility/vector1.hh>



namespace core {
namespace fragment {

template< class T>
//typedef core::Real T;
class MapCacheUnit : public core::fragment::BaseCacheUnit {
	typedef std::map< core::Size, T > TMap;
public:
	BaseCacheUnitOP clone() const {
		return BaseCacheUnitOP( new MapCacheUnit<T> );
	};

	void remap_value( BaseCacheUnit const& source, Size source_id, Size new_id ) {
		T value;
		dynamic_cast< MapCacheUnit<T> const& > (source).retrieve( source_id, value );
		store( new_id, value );
	};


	bool retrieve( core::Size frag_id, T& value ) const {
		typename TMap::const_iterator iter( map_.find( frag_id ) );
		if ( iter == map_.end() ) {
			return false;
		}	else {
			value = iter->second;
			return true;
		};
	};


	void store( Size frag_id, T const& value ) {
		map_[ frag_id ] = value;
	}

	void register_frag_id( Size ) { }; //do nothing -- cache has lazy evaluation
private:
	TMap map_;
};

template< class T>
//typedef core::Real T;
class VectorCacheUnit : public core::fragment::BaseCacheUnit {
	typedef utility::vector1< T > TVector;
public:
	BaseCacheUnitOP clone() const {
		return BaseCacheUnitOP( new VectorCacheUnit<T> );
	}

	void remap_value( BaseCacheUnit const& source, Size source_id, Size new_id ) {
		T value;
		dynamic_cast< VectorCacheUnit<T> const& > (source).retrieve( source_id, value );
		store( new_id, value );
	};


	bool retrieve( core::Size frag_id, T& value ) const {
		value = list_[ frag_id ];
		return true;
	}

	T const& retrieve( core::Size frag_id ) const {
		return list_[ frag_id ];
	}

	void store( Size frag_id, T const& value ) {
		if ( frag_id > list_.size() ) {
			list_.resize( frag_id );
		};
		list_[ frag_id ] = value;
	}

	void register_frag_id( Size frag_id ) {
		if ( frag_id > list_.size() ) {
			list_.resize( frag_id );
		};
	};
private:
	TVector list_;
};

template< class T, class XCacheUnit >
class CacheWrapper {
public:
	typedef utility::vector1< core::Size > IndexList;
	typedef utility::vector1< core::Real > ScoreList;
	typedef XCacheUnit TCacheUnit;
	typedef utility::pointer::weak_ptr< TCacheUnit > TCacheUnitAP;
	typedef utility::pointer::shared_ptr< TCacheUnit > TCacheUnitOP;
	typedef std::pair< FragID, T > ScoredFrag;
	typedef utility::vector1< ScoredFrag> ScoredList;

public:
	CacheWrapper( std::string tag ) :
		tag_ ( tag ),
		new_cache_ ( new TCacheUnit )
	{};

	~CacheWrapper() { };

	bool retrieve( Frame const& frame, core::Size frag_num, T& score) const {
		return cache( frame ).retrieve( frame.frag_id( frag_num ), score );
	};

	bool retrieve( FragID const& frag_id, T& score ) const {
		return cache( frag_id.frame() ).retrieve( frag_id.id(), score );
	}

	T retrieve( core::Size frag_id ) const {
		T val;
		if ( retrieve( frag_id, val ) ) {
				return val;
		} else {
			throw utility::excn::EXCN_RangeError( "no "+tag_+ "entry found for fragment: ");
		}
	}

	T retrieve( FragID const& frag_id ) const {
		T val;
		if ( retrieve( frag_id, val ) ) {
				return val;
		} else {
			throw utility::excn::EXCN_RangeError( "no "+tag_+ "entry found for fragment: ");
		}
	}

T retrieve( Frame const& frame, core::Size frag_num ) const {
		T val;
		if ( retrieve( frame, frag_num, val ) ) {
				return val;
		} else {
			throw utility::excn::EXCN_RangeError( "no "+tag_+ "entry found for fragment");
		}
	}


	void store( Frame const& frame, core::Size frag_num, T const& score ) {
		cache( frame ).store( frame.frag_id( frag_num ), score );
	};

	void store( FragID const& frag_id, T const& score ) {
		cache( frag_id.frame() ).store( frag_id.id() , score );
	};

	void scored_frag_ids( ScoredList &frag_ids, FragID_Iterator begin, FragID_Iterator end, T* empty = NULL ) const {
		for ( FragID_Iterator it = begin; it!=end; ++it ) {
			T score;
			if ( retrieve( *it, score ) )
				frag_ids.push_back( ScoredFrag( *it,  score ) );
			else if ( empty ) {
				frag_ids.push_back( ScoredFrag( *it, *empty ) );
			}
		};
	}

	void scored_frag_ids( ScoredList &frag_ids, FragID_Iterator begin, FragID_Iterator end, T empty  ) const {
		scored_frag_ids( frag_ids, begin, end, &empty );
	}

	TCacheUnit const& cache( Frame const& frame ) const {
		//    TCacheUnitAP ptr=utility::pointer::dynamic_pointer_cast< TCacheUnit >( frame.cache( tag_, new_cache_ ) );
		// assert(ptr);
		// better with refernce, throws exception automatic if it goes wrong
		return dynamic_cast< TCacheUnit const& >( frame.cache( tag_, new_cache_ ) );
	}

	TCacheUnit& cache( Frame const& frame ) {
		//    TCacheUnitAP ptr=utility::pointer::dynamic_pointer_cast< TCacheUnit >( frame.cache( tag_, new_cache_ ) );
		// assert(ptr);
		// better with refernce, throws exception automatic if it goes wrong
		return dynamic_cast< TCacheUnit& >( frame.cache( tag_, new_cache_ ) );
	}


	TCacheUnit& operator() ( Frame const& frame ) {
		return cache( frame );
		//    return utility::pointer::dynamic_pointer_cast< TCacheUnit >( frame.cache( tag_, new_cache_ ) );
	}

	std::string tag_;
	TCacheUnitOP new_cache_;
};

//FragCache uses MapCacheUnit, .i.e, some values might be not set --> retrieve returns false
//FragStore uses VectorCacheUnit, i.e., all values should be valid --> the user has to take care that every fragment in the frame

template < class T >
class FragCache : public CacheWrapper< T, MapCacheUnit< T> > {
	typedef CacheWrapper< T, MapCacheUnit< T> > Parent;
	//	typedef Parent::TCacheUnit TCacheUnit;
	typedef T ValueType;
public:
	FragCache( std::string tag ) : Parent( tag ) {};
};

template < class T >
class FragStore : public CacheWrapper< T, VectorCacheUnit< T> > {
	typedef CacheWrapper< T, VectorCacheUnit< T> > Parent;
	//	typedef Parent::TCacheUnit TCacheUnit;
	typedef T ValueType;
public:
	FragStore( std::string tag ) : Parent( tag ) {};
};


} //fragment
} //core

#endif
