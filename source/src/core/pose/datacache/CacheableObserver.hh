// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/datacache/CacheableObserver.hh
/// @brief  Base class for Pose/Conformation observers that are stored in
///         a Pose's ObserverCache.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pose_datacache_CacheableObserver_hh
#define INCLUDED_core_pose_datacache_CacheableObserver_hh

// unit headers
#include <core/pose/datacache/CacheableObserver.fwd.hh>
#include <core/pose/datacache/CacheableObserverType.hh>

// project headers
#include <core/pose/Pose.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace datacache {

/// @brief Base class for Pose/Conformation observers that are stored in
///  a Pose's DataCache.
/// @details Classes derived from CacheableObserver whose instances are
///  stored inside a Pose's ObserverCache will be cloned and re-attached to the
///  new Pose instance upon Pose copy construction/assignment.
///  CacheableObservers should not need to attach themselves independently
///  unless necessary -- attach/detach should be done through the ObserverCache's
///  interface.  For example, the ObserverCache can attach the observer upon
///  calling ObserverCache::set().
/// @warning When deriving from this class remember that Links or pointers to
///  subjects stored in classes are *not* to be copied when copy
///  constructor/assignment is called.
class CacheableObserver : public utility::pointer::ReferenceCount {


private: // typedefs


	typedef utility::pointer::ReferenceCount Super;


public: // construct/destruct


	/// @brief default constructor
	CacheableObserver();


	/// @brief copy constructor
	/// @warning Subject being observed (represented by Link/pointer) is not copied!
	CacheableObserver( CacheableObserver const & rval );


	/// @brief default destructor
	/// @warning Derived classes must remember to detach on destruction!
	virtual
	~CacheableObserver();


public: // assignment


	/// @brief copy assignment
	/// @warning Subject being observed (represented by Link/pointer) is not copied!
	CacheableObserver & operator =( CacheableObserver const & rval );


public: // virtual constructors


	/// @brief clone this object
	/// @warning Subject (represented by Link/pointer) is not copied!
	virtual
	CacheableObserverOP clone() = 0;


	/// @brief create a new instance of this object
	virtual
	CacheableObserverOP create() = 0;


public: // observer interface


	/// @brief attach to Pose/Conformation
	///  Derived classes do not overload this method -- see attach_impl()
	///  instead.
	void attach_to( Pose & pose );


	/// @brief detach from Pose/Conformation
	/// @remarks Derived classes do not overload this method -- see
	///  detach_impl() instead.
	void detach_from();


public: // virtual observer interface


	/// @brief is this observer attached to a Pose/Conformation?
	virtual
	bool is_attached() const = 0;
	
protected: // virtual observer interface


	/// @brief attach to Pose/Conformation
	virtual
	void attach_impl( Pose & pose ) = 0;


	/// @brief detach from Pose/Conformation
	virtual
	void detach_impl() = 0;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace pose
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_datacache_CacheableObserver )
#endif // SERIALIZATION


#endif /* INCLUDED_core_pose_datacache_CacheableObserver_HH */
