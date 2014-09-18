// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Entity.hh
/// @brief the unit employed/optimized by GeneticAlgorithm
/// @author ashworth

#ifndef INCLUDED_protocols_genetic_algorithm_Entity_hh
#define INCLUDED_protocols_genetic_algorithm_Entity_hh

// Unit headers
#include <protocols/genetic_algorithm/Entity.fwd.hh>

#include <core/types.hh>

// Utility headers
// AUTO-REMOVED #include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/factory/WidgetFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

// Boost headers
#include <boost/functional/hash.hpp> // hash_range

///#include <ObjexxFCL/format.hh> // apl needed?

// C++ headers
// AUTO-REMOVED #include <iostream>
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace protocols {
namespace genetic_algorithm {

/// Entity element

class EntityElement : public utility::pointer::ReferenceCount {
public:
	typedef core::Size Size;
	typedef utility::pointer::ReferenceCount parent;

public:
	EntityElement();
	EntityElement( Size index );
	EntityElement( std::string & word ); // This constructor nibbles away at the input string
	virtual ~EntityElement();

	virtual EntityElementOP clone() = 0;
	virtual EntityElementOP fresh_instance() = 0;

	Size index() const;
	void index( Size index );
	virtual Size hash() const = 0;

	virtual bool operator <  ( EntityElement const & rhs ) const;
	virtual bool operator == ( EntityElement const & rhs ) const;
	virtual EntityElement const & operator = ( EntityElement const & rhs );

	virtual std::string to_string() const;
	virtual std::string name() const = 0; // Each entity element must have a distinct name

private:
	Size index_;
};

/// Entity element creator

class EntityElementCreator : public utility::pointer::ReferenceCount {
public:
	virtual ~EntityElementCreator();
	virtual std::string widget_name() const = 0;
	virtual EntityElementOP new_entity( std::string const & word ) = 0;
};

/// Entity element factory

class EntityElementFactory : public utility::factory::WidgetFactory< EntityElementCreator > {
public:
	virtual ~EntityElementFactory() {};

	static EntityElementFactory * get_instance();

	EntityElementOP element_from_string( std::string const & );

	virtual std::string factory_name() const;

	//void register_creator( EntityElementCreatorOP creator );

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:
	EntityElementFactory();
	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static EntityElementFactory * create_singleton_instance();
private:
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< EntityElementFactory * > instance_;
#else
	static EntityElementFactory * instance_;
#endif

};

template < class T >
class EntityElementRegistrator : public utility::factory::WidgetRegistrator< EntityElementFactory, T > {
public:
	typedef typename utility::factory::WidgetRegistrator< EntityElementFactory, T > parent;
public:
	EntityElementRegistrator() : parent() {}
};


/// Entity: a vector of EntityElements used to describe the state for a system
/// under optimization.

class Entity : public utility::pointer::ReferenceCount {

public:
	typedef utility::pointer::owning_ptr< Entity > OP;
	typedef utility::pointer::owning_ptr< Entity const > COP;
	typedef utility::pointer::access_ptr< Entity const > CAP;
	typedef utility::vector1< COP > COPs;
	typedef utility::vector1< CAP > CAPs;

public:
	Entity();
	////@brief construct a duplicate Entity from another entity
	Entity( Entity const & entity );
	Entity const & operator = ( Entity const & );

	virtual ~Entity();

	////@brief construct Entity from std::string (e.g. from file)
	Entity( std::string const & line );

	virtual OP clone() const;
	virtual void set_traits_size( core::Size size /*, EntityElementOP element*/ );
	virtual void set_traits( EntityElements const & traits );
	virtual EntityElements const & traits() const;
	virtual void set_entity_element( core::Size index, EntityElementOP element );
	virtual void set_fitness( core::Real val );
	virtual core::Real fitness() const;
	virtual bool fitness_valid() const;
	virtual bool operator == ( Entity const & other ) const;
	virtual bool operator < ( Entity const & other ) const;
	virtual void show( std::ostream & os ) const;
	virtual std::string to_string() const;
	virtual std::string traits_string() const;

	virtual void write_checkpoint( std::ostream & os ) const;
	virtual bool read_checkpoint( std::istream & is );

private:
	EntityElements traits_;
	core::Real fitness_;
	bool fitness_valid_;
};

std::ostream & operator << ( std::ostream & os, Entity const & entity );


///@brief for sorting owning pointers by that to which they point
template <typename T>
bool lt_OP_deref(
	utility::pointer::owning_ptr<T> const & a,
	utility::pointer::owning_ptr<T> const & b
)
{
	if ( !a && !b ) return false;
	if ( !a ) return true;
	if ( !b ) return false;
	return *a < *b;
}

///@brief for assessing equality between owning pointers by that to which they point
template <typename T>
bool eq_OP_deref(
	utility::pointer::owning_ptr<T> const & a,
	utility::pointer::owning_ptr<T> const & b
)
{
	if ( !a && !b ) return true;
	if ( !a || !b ) return false;
	return *a == *b;
}

struct
Vec1Hash {
	std::size_t operator() ( EntityElements const & vec1 ) const {
		std::size_t seed = 0;
		for ( EntityElements::const_iterator
				iter = vec1.begin(), iter_end = vec1.end(); iter != iter_end; ++iter ) {
			boost::hash_combine( seed, (*iter)->hash() );
		}
		return seed;
	}
};

struct
EntityElementsEqual{
	bool operator() (
		EntityElements const & elems1,
		EntityElements const & elems2
	) const {
		if ( elems1.size() != elems2.size() ) return false;
		for ( core::Size ii = 1; ii <= elems1.size(); ++ii ) {
			if ( ! ((*elems1[ ii ]) == (*elems2[ ii ])) ) return false;
		}
		return true;
	}
};

} // namespace genetic_algorithm
} // namespace protocols

#endif
