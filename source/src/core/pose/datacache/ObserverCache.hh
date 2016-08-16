// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/ObserverCache.hh
/// @brief  A DataCache storing objects derived from
///         basic::datacache::CacheableData.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pose_datacache_ObserverCache_hh
#define INCLUDED_core_pose_datacache_ObserverCache_hh

// unit headers
#include <core/pose/datacache/ObserverCache.fwd.hh>

// package headers

// project headers
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataCache.hh>

#include <core/pose/datacache/CacheableObserver.fwd.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace datacache {


/// @brief A DataCache storing Pose/Conformation observers derived from
///  core::pose::datacache::CacheableObserver.
class ObserverCache : public basic::datacache::DataCache< CacheableObserver > {


private: // typedefs


	typedef basic::datacache::DataCache< CacheableObserver > Super;


public: // typedefs


	typedef Super::Size Size;


protected: // typedefs


	typedef Super::DataOPs ObserverOPs;


public: // construct/destruct


	/// @brief constructor
	/// @param[in] n_types The number of slots for this ObserverCache.
	/// @param[in] pose The Pose that will be watched by the Observers in this cache.
	ObserverCache(
		Size n_slots,
		Pose & pose
	);


	/// @brief default destructor
	virtual
	~ObserverCache();


private: // disallowed constructors


	/// @brief copy constructor
	ObserverCache( ObserverCache const & rval );


public: // assignment


	/// @brief copy assignment
	ObserverCache & operator =( ObserverCache const & rval );


public: // methods


	/// @brief clear all the observers
	void clear();


	/// @brief clear the observer in a selected slot
	void clear( Size slot );


	/// @brief store a copy of the observer in the given slot and attach it to
	///  the Pose
	/// @param[in] The slot to use.
	/// @param[in] observer The Observer to clone() and store.
	/// @remarks this function exists to ensure the base class version is
	///  overridden
	void set(
		Size slot,
		CacheableObserverOP observer
	);


	/// @brief store a copy of the observer in the given slot
	/// @param[in] The slot to use.
	/// @param[in] observer The Observer to clone() and store.
	/// @param[in] auto_attach Attach the observer to the Pose?
	void set(
		Size slot,
		CacheableObserverOP observer,
		bool auto_attach
	);


public: // observer interface


	/// @brief is the observer in the slot attached to the Pose?
	/// @return true if attached, false if not attached or no observer
	///  exists in the slot
	bool is_attached( Size slot ) const;

	/// @brief attach all stored observers to the Pose
	void attach();


	/// @brief detach all observers from the Pose
	void detach();


	/// @brief attach an observer in a particular slot to the Pose
	/// @param[in] slot Attach the observer in this slot.
	void attach( Size slot );


	/// @brief detach an observer in a particular slot to the Pose
	/// @param[in] slot Detach the observer in this slot.
	void detach( Size slot );


private: // data

	/// @brief Which of the cacheable observers are attached to the pose, and which
	/// are just hanging out.
	utility::vector1< bool > attached_;

	/// @brief The Pose being watched by the observers in this cache.
	/// @remarks This must be a pointer and not an owning pointer in case
	///  the Pose is on the stack.
	Pose * pose_;


#ifdef    SERIALIZATION
protected:

	/// @brief default constructor -- used only during deserialization
	ObserverCache();
	friend class cereal::access;

public:
	/// @brief Attach the ObserverCache to a particular Pose; this is not the preferred
	/// method of attaching observers to Poses -- the preferred method is to set it the
	/// particular pose in the constructor.  It is, however, necessary for deserialization.
	/// It may be called only once and it requires that the (private) default constructor
	/// have been used to instantiate the class in the first place.
	void attach_pose( Pose & pose );

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace pose
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_datacache_ObserverCache )
#endif // SERIALIZATION


#endif /* INCLUDED_core_pose_datacache_ObserverCache_HH */
