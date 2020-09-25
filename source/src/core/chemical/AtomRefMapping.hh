// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/AtomRefMapping.hh
/// @brief A mapping of atom references (vds, indexes, names) to each other
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_AtomRefMapping_hh
#define INCLUDED_core_chemical_AtomRefMapping_hh

// Unit headers

// Project headers
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/types.hh>

// Utility headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace chemical {

/// @brief An internal implementation class, to get around the fact that functions cannot be partially specialized.
template<class OutRef>
class RefConvert {
	// C++ will ignore the "base" template for the specializations below.
	// The functions here are just to document the interface.
	// They deliberately don't have an implementation.

	/// @brief The value for invalid items of the class type
	static OutRef const &
	invalid();

	/// @brief For the given residue type, convert the given reference into the class reference type
	template<class ResTypeType, class InRef>
	static OutRef
	convert(ResTypeType const &, InRef const &);

	/// @brief Will pass through the value if InRef is the same type as OutRef,
	/// but will result in the second value (invalid) if it's not.
	template<class InRef>
	static OutRef const &
	identity_pass_through(InRef const &, OutRef const &);
};

///////////////////////////////////////
// The actual implementations, using template specialization

template<>
class RefConvert<VD> {
public:

	static
	VD const & invalid() {
		static VD const inval( ResidueGraph::null_vertex() );
		return inval;
	}

	static
	VD convert(MutableResidueType const &, VD const & ref ) { return ref; }

	template<class InRef>
	static
	VD convert(MutableResidueType const & restype, InRef const & ref ) { return restype.atom_vertex( ref ); }

	static
	VD const &
	identity_pass_through(VD const & ref, VD const & ) { return ref; }

	template<class InRef>
	static
	VD const &
	identity_pass_through(InRef const &, VD const & second) { return second; }
};

template<>
class RefConvert<core::Size> {
public:

	static
	core::Size const & invalid() {
		static core::Size const inval( 0 );
		return inval;
	}

	static
	core::Size convert(ResidueType const &, core::Size const & ref ) { return ref; }

	template<class InRef>
	static
	core::Size convert(ResidueType const & restype, InRef const & ref ) { return restype.atom_index( ref ); }

	static
	core::Size const &
	identity_pass_through(core::Size const & ref, core::Size const & ) { return ref; }

	template<class InRef>
	static
	core::Size const &
	identity_pass_through(InRef const &, core::Size const & second) { return second; }
};

template<>
class RefConvert<std::string> {
public:

	static
	std::string const & invalid() {
		static std::string const inval( "" );
		return inval;
	}

	static
	std::string convert(ResidueTypeBase const &, std::string const & ref ) { return ref; }

	template<class ResTypeType, class InRef>
	static
	std::string convert(ResTypeType const & restype, InRef const & ref ) { return restype.atom_name( ref ); }

	static
	std::string const &
	identity_pass_through(std::string const & ref, std::string const & ) { return ref; }

	template<class InRef>
	static
	std::string const &
	identity_pass_through(InRef const &, std::string const & second) { return second; }
};

/// @brief A class for mapping ResidueType atom references from one to another.
/// It's intended not only for intra ResidueType mapping, but also for mapping
/// corresponding references from one ResidueType to another.
/// @details This is a templated class to cut down on combinitorial complexity.

template< class FromRef, class ToRef >
class AtomRefMapping : public utility::VirtualBase {

private:
	typedef std::map< FromRef, ToRef > SubMap;

public:

	typedef typename SubMap::const_iterator const_iterator;

	AtomRefMapping():
		invalid_from_( RefConvert<FromRef>::invalid() ),
		invalid_to_( RefConvert<ToRef>::invalid() ),
		identity_(false)
	{}

	/// @brief Constructor with non-standard invalid values.
	AtomRefMapping(FromRef const & invalid_from, ToRef const & invalid_to):
		invalid_from_( invalid_from ),
		invalid_to_( invalid_to ),
		identity_(false)
	{}

	/// @brief Ctor with explicit identity setting.
	AtomRefMapping(bool identity):
		invalid_from_( RefConvert<FromRef>::invalid() ),
		invalid_to_( RefConvert<ToRef>::invalid() ),
		identity_(identity)
	{}

	/// @brief Make an explicit identity mapping for the ResType
	/// @details Not availible if either ref is a core::Size
	template< class F = FromRef, class T = ToRef, typename = typename std::enable_if< !std::is_same< F, core::Size >::value && !std::is_same< T, core::Size >::value, int >::type >
	AtomRefMapping(MutableResidueType const & restype):
		invalid_from_( RefConvert<FromRef>::invalid() ),
		invalid_to_( RefConvert<ToRef>::invalid() ),
		identity_(false)
	{
		VIter iter, iter_end;
		for ( boost::tie( iter, iter_end) = restype.atom_iterators(); iter != iter_end; ++iter ) {
			mapping_[ RefConvert<FromRef>::convert(restype,*iter) ] = RefConvert<ToRef>::convert(restype,*iter);
		}
	}

	/// @brief Make an explicit identity mapping for the ResType
	/// @details Not availible if either ref is a VD
	template< class F = FromRef, class T = ToRef, typename = typename std::enable_if< !std::is_same< F, VD >::value && !std::is_same< T, VD >::value, int >::type >
	AtomRefMapping(ResidueType const & restype):
		invalid_from_( RefConvert<FromRef>::invalid() ),
		invalid_to_( RefConvert<ToRef>::invalid() ),
		identity_(false)
	{
		for ( core::Size ii(1); ii <= restype.natoms(); ++ii ) {
			mapping_[ RefConvert<FromRef>::convert(restype,ii) ] = RefConvert<ToRef>::convert(restype,ii);
		}
	}

	/// @brief What's the invalid value for the keys?
	FromRef const & invalid_key() const { return invalid_from_; }

	/// @brief What's the invalid value for the mapped values?
	ToRef const & invalid_entry() const { return invalid_to_; }

	/// @brief Set the identity flag
	void identity( bool setting ) { identity_ = setting; }

	/// @brief Get the identity setting
	bool identity() const { return identity_; }

	/// @brief go from an A->B mapping to a B->A mapping
	AtomRefMapping< ToRef, FromRef > reverse() const {
		AtomRefMapping< ToRef, FromRef > reversed(invalid_to_,invalid_from_);
		typename AtomRefMapping< ToRef, FromRef >::SubMap & new_mapping( reversed.mapping_ );

		for ( typename SubMap::const_iterator iter(mapping_.begin()), iter_end(mapping_.end()); iter != iter_end; ++iter ) {
			if ( iter->second == invalid_to_ ) {
				continue; // We never want to map invalids to anything.
			}
			if ( new_mapping.count( iter->second ) ) {
				//TR.Warning << "Atom Reference " << iter->second << " found multiple times for inverted map - using first definition." << std::endl;
				continue;
			}
			new_mapping[ iter->second ] = iter->first;
		}
		// Identity stays the same.
		reversed.identity_ = identity_;
		return reversed;
	}

	/// @brief Apply a B->C mapping to the current A->B mapping to get an A->C mapping
	/// i.e. vmap[j] becomes vmap_to_add[ vmap[j] ]
	template< class OToRef >
	AtomRefMapping<FromRef,OToRef>
	downstream_combine( AtomRefMapping<ToRef,OToRef> const & vmap_to_add ) const {
		AtomRefMapping<FromRef,OToRef> combined(invalid_from_, vmap_to_add.invalid_to_);
		typename AtomRefMapping<FromRef,OToRef>::SubMap & new_mapping(combined.mapping_);
		// First map explicits
		for ( typename SubMap::const_iterator iter(mapping_.begin()), iter_end(mapping_.end()); iter != iter_end; ++iter ) {
			if ( iter->first == invalid_from_ ) continue;
			if ( vmap_to_add.identity_ || vmap_to_add.count( iter->second ) ) {
				new_mapping[ iter->first ] = vmap_to_add[ iter->second ];
			}
		}
		// Now add second map explicits where the first map is implicit.
		if ( identity_ ) {
			for ( typename AtomRefMapping<ToRef,OToRef>::SubMap::const_iterator iter(vmap_to_add.mapping_.begin()), iter_end(vmap_to_add.mapping_.end()); iter != iter_end; ++iter ) {
				// If we can't convert it into the appropriate type, ignore it.
				FromRef const & from( RefConvert<FromRef>::identity_pass_through(iter->first, invalid_from_) );
				if ( from == invalid_from_ ) continue;
				// Don't use first map implicit identity if it's already explicit
				if ( ! new_mapping.count(from) ) {
					new_mapping[from] = iter->second;
				}
			}
		}
		// The combined only has an implicit identity pass though if both do.
		combined.identity_ = identity_ && vmap_to_add.identity_;
		return combined;
	}

	/// @brief Apply a C->A mapping to the current A->B mapping to get a C->B mapping
	/// i.e. vmap[j] becomes vmap[ vmap_to_add[ j ] ]
	template< class OFromRef >
	AtomRefMapping<OFromRef, ToRef>
	upstream_combine( AtomRefMapping<OFromRef,FromRef> const & vmap_to_add ) const {
		return vmap_to_add.downstream_combine(*this);
	}

	// void show() const {
	//  show( TR );
	// }

	void show( std::ostream & output ) const {
		output << "core.chemical.AtomRefMapping  " << ( identity_ ? "identity items except for" : "" );
		output << "-----------------------> " << std::endl;
		for ( typename SubMap::const_iterator iter(mapping_.begin()), iter_end(mapping_.end()); iter != iter_end; ++iter ) {
			if ( iter->first == invalid_from_ ) continue;
			if ( identity_ && iter->first == RefConvert<FromRef>::identity_pass_through(iter->second,invalid_from_) ) continue;
			output << "\t" << iter->first << "\t" << iter->second << std::endl;
		}
		output << "----------------------- " << std::endl;
	}

	/// @brief Ignore all specializations on this map.
	void clear() {
		mapping_.clear();
	}

	/// @brief Make in act like the default
	/// (pay attention to identity_)
	void erase( FromRef const & in ) {
		mapping_.erase( in );
	}

	/// @brief Make the key value an invalid to value
	/// (Will work even with identity_ == true)
	void invalidate( FromRef const & in ) {
		mapping_[ in ]  = invalid_to_;
	}

	ToRef const operator[]( FromRef const & in ) const {
		if ( in == invalid_from_ ) {
			return invalid_to_; // Invalid always maps to invalid.
		}
		typename SubMap::const_iterator entry( mapping_.find( in ) );
		if ( entry == mapping_.end() ) {
			if ( identity_ ) {
				return RefConvert<ToRef>::identity_pass_through(in,invalid_to_);
			} else {
				return invalid_to_;
			}
		} else {
			return entry->second;
		}
	}

	ToRef & operator[]( FromRef const & in ) {
		if ( in == invalid_from_ ) {
			mapping_[in] = invalid_to_; //Always reset to invalid
		}
		if ( ! count(in) ) {
			if ( identity_ ) {
				mapping_[in] = RefConvert<ToRef>::identity_pass_through(in,invalid_to_);
			} else {
				mapping_[in] = invalid_to_;
			}
		}
		return mapping_[in];
	}

	/// @brief Gets the key which corresponds to the given mapped value, if any.
	FromRef const & reverse_lookup( ToRef const & out ) const {
		if ( out == invalid_to_ ) {
			return invalid_from_; // Invalid always maps to invalid.
		}
		for ( typename SubMap::const_iterator iter( mapping_.begin() ), iter_end( mapping_.end() ); iter != iter_end; ++iter ) {
			if ( iter->second == out ) {
				return iter->first;
			}
		}
		if ( identity_ ) {
			FromRef const & converted( RefConvert<FromRef>::identity_pass_through(out,invalid_from_) );
			if ( converted != invalid_from_ && count(converted) == 0 ) {
				// Can only be identity pass through if the from isn't explicitly mentioned for another conversion.
				return converted;
			}
		}
		return invalid_from_;
	}

	/// @brief Does this map have the given key?
	bool count( FromRef const & key ) const {
		if ( key == invalid_from_ ) { return true; } // Invalid entry is always present and always evaluates to invalid entry.
		return mapping_.count( key );
	}

	/// @brief Does this mapping contain any explicit (non-identity) mappings?
	bool empty() const {
		return mapping_.empty();
	}

	/// @brief Equality operator.
	bool operator==( AtomRefMapping<FromRef, ToRef> const & other ) const {
		if ( identity_ != other.identity_ ) return false;

		// Now we've checked the pass through, make sure all of the explicit maps match.
		// We can't just match the submaps, as some implicit pass-throughs might become explicit.
		for ( typename SubMap::const_iterator iter(mapping_.begin()), iter_end(mapping_.end()); iter != iter_end; ++iter ) {
			if ( iter->first == invalid_from_ ) { continue; }
			if ( ! other.count( iter->first ) || other[iter->first] != iter->second ) {
				return false;
			}
		}
		// We do this both ways, as one might might have explicit references the other might not.
		for ( typename SubMap::const_iterator iter(other.mapping_.begin()), iter_end(other.mapping_.end()); iter != iter_end; ++iter ) {
			if ( iter->first == other.invalid_from_ ) { continue; }
			if ( ! count( iter->first ) || (*this)[iter->first] != iter->second ) {
				return false;
			}
		}

		return true;
	}

	/// @brief Make a Size-to-Size vector, indexed by the from, with elements that specify the two entries.
	template< class F = FromRef, class T = ToRef, typename = typename std::enable_if< !std::is_same< F, VD >::value && !std::is_same< T, VD >::value, int >::type >
	utility::vector1< core::Size >
	vectorize(ResidueType const & from_res, ResidueType const & to_res ) const {
		utility::vector1< core::Size > out_vec( to_res.natoms(), 0 );
		for ( typename SubMap::const_iterator iter(mapping_.begin()), iter_end(mapping_.end()); iter != iter_end; ++iter ) {
			core::Size from_index( RefConvert<core::Size>::convert(from_res, iter->first ) );
			if ( from_index == 0 ) { continue; } // Explicitly code the invalid here, as we want zero-based invalids.
			if ( from_index > out_vec.size() ) {
				out_vec.resize(from_index,0);
			}
			out_vec[ from_index ] = RefConvert<core::Size>::convert(to_res, iter->second );
		}
		return out_vec;
	}

	/// @brief An iterator to std::pair(key,value) pairs for *explicit* entries
	const_iterator // AtomRefMapping::const_iterator
	begin() const {
		return mapping_.begin();
	}

	/// @brief An iterator to std::pair(key,value) pairs for *explicit* entries
	const_iterator // AtomRefMapping::const_iterator
	end() const {
		return mapping_.end();
	}

private:
	// Needed so AtomRefMappings with different template parameters can access each other's private data.
	// >>> The template parameters are *intentionally* different from FromRef and ToRef <<<
	template<class Foo, class Bar> friend class AtomRefMapping;

	/// @brief The actual mapping.
	SubMap mapping_;

	/// @brief What value to use for invalid from values
	FromRef invalid_from_;
	/// @brief What value to use for invalid to values
	ToRef invalid_to_;

	/// @brief Should non-explicitly stated atoms be treated like an identity relation? Only meaningful if the two types are the same.
	bool identity_;

}; // class AtomRefMapping

template<class FromRef, class ToRef>
inline std::ostream & operator<< (
	std::ostream & out,
	AtomRefMapping<FromRef, ToRef> const & map
) {
	map.show( out );
	return out;
}

template<class ARef, class BRef, class CRef>
AtomRefMapping<ARef,CRef>
combine( AtomRefMapping<ARef,BRef> const & first, AtomRefMapping<BRef,CRef> const & second ) {
	return first.downstream_combine(second);
}

typedef AtomRefMapping<VD,VD> VDVDMapping;
typedef AtomRefMapping<VD,core::Size> VDIndexMapping;
typedef AtomRefMapping<core::Size,VD> IndexVDMapping;
typedef AtomRefMapping<core::Size,core::Size> IndexIndexMapping;

typedef AtomRefMapping<std::string,VD> StringVDMapping;
typedef AtomRefMapping<VD,std::string> VDStringMapping;

// Additional aliases
typedef AtomRefMapping<std::string,VD> NameVDMapping;
typedef AtomRefMapping<VD,std::string> VDNameMapping;

typedef AtomRefMapping<std::string,core::Size> NameIndexMapping;
typedef AtomRefMapping<core::Size,std::string> IndexNameMapping;

/*
// PyRosetta on has issues with having the explicit instantiations here
template class AtomRefMapping<VD,VD>;
template class AtomRefMapping<VD,core::Size>;
template class AtomRefMapping<core::Size,VD>;
template class AtomRefMapping<core::Size,core::Size>;

template class AtomRefMapping<VD,std::string>;
template class AtomRefMapping<std::string,std::string>;
template class AtomRefMapping<std::string,VD>;

template class AtomRefMapping<core::Size,std::string>;
template class AtomRefMapping<std::string,core::Size>;
*/

} // chemical
} // core

#endif
