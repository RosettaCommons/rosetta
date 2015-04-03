// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
		Size const n_slots,
		Pose & pose
	);


	/// @brief default destructor
	virtual
	~ObserverCache();


private: // disallowed constructors


	/// @brief default constructor
	ObserverCache();


	/// @brief copy constructor
	ObserverCache( ObserverCache const & rval );


public: // assignment


	/// @brief copy assignment
	ObserverCache & operator =( ObserverCache const & rval );


public: // methods


	/// @brief clear all the observers
	void clear();


	/// @brief clear the observer in a selected slot
	void clear( Size const slot );


	/// @brief store a copy of the observer in the given slot and attach it to
	///  the Pose
	/// @param[in] The slot to use.
	/// @param[in] observer The Observer to clone() and store.
	/// @remarks this function exists to ensure the base class version is
	///  overridden
	void set(
		Size const slot,
		CacheableObserverOP observer
	);


	/// @brief store a copy of the observer in the given slot
	/// @param[in] The slot to use.
	/// @param[in] observer The Observer to clone() and store.
	/// @param[in] auto_attach Attach the observer to the Pose?
	void set(
		Size const slot,
		CacheableObserverOP observer,
		bool const auto_attach
	);


public: // observer interface


	/// @brief is the observer in the slot attached to the Pose?
	/// @return true if attached, false if not attached or no observer
	///  exists in the slot
	bool is_attached( Size const slot ) const;


	/// @brief attach all stored observers to the Pose
	void attach();


	/// @brief detach all observers from the Pose
	void detach();


	/// @brief attach an observer in a particular slot to the Pose
	/// @param[in] slot Attach the observer in this slot.
	void attach( Size const slot );


	/// @brief detach an observer in a particular slot to the Pose
	/// @param[in] slot Detach the observer in this slot.
	void detach( Size const slot );


private: // data


	/// @brief The Pose being watched by the observers in this cache.
	/// @remarks This must be a pointer and not an owning pointer in case
	///  the Pose is on the stack.
	Pose * pose_;


};


} // namespace datacache
} // namespace pose
} // namespace core


#endif /* INCLUDED_core_pose_datacache_ObserverCache_HH */
