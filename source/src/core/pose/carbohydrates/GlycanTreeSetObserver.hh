// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/carbohydrates/GlycanTreeSetObserver.hh
/// @brief The CacheablePoseObserver version of GlycanTreeSet that will react to pose length changes..
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_core_pose_carbohydrates_GlycanTreeSetObserver_hh
#define INCLUDED_core_pose_carbohydrates_GlycanTreeSetObserver_hh

#include <core/pose/carbohydrates/GlycanTreeSetObserver.fwd.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.fwd.hh>
#include <core/pose/datacache/CacheableObserver.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/signals/Link.hh>


namespace core {
namespace pose {
namespace carbohydrates {


///@brief
///
///  Get the GlycanTreeSet from the pose.
///
///
///@details
///
/// This is so we can get the glycan_tree_set using a const pose!
///
conformation::carbohydrates::GlycanTreeSetCOP
get_glycan_tree_set(Pose const & pose);

///@brief The CacheablePoseObserver version of GlycanTreeSet that will react to pose length changes..
class GlycanTreeSetObserver : public core::pose::datacache::CacheableObserver {

public:

	GlycanTreeSetObserver();

	///@brief Construct the GlycanTreeSet, but do not attach it to any pose.
	GlycanTreeSetObserver( conformation::Conformation const & conf );

	///@Brief Construct the GlycanTreeSet and attach this object to the pose.
	GlycanTreeSetObserver( core::pose::Pose & pose );

	GlycanTreeSetObserver(GlycanTreeSetObserver const & src);

	virtual ~GlycanTreeSetObserver();

	GlycanTreeSetObserverOP
	clone() const;

public:

	///@brief Get the GlycanTreeSet that is maintained by this Observer.
	conformation::carbohydrates::GlycanTreeSetCOP
	get_glycan_tree_set() const;

public:

	//////////////////////////////////////////////////////////////////
	///                   ///
	///                Cacheable Observer Functions                ///
	///                  ///
	//////////////////////////////////////////////////////////////////


	virtual
	core::pose::datacache::CacheableObserverOP
	clone();

	virtual
	core::pose::datacache::CacheableObserverOP
	create();


public: //observer interface

	virtual
	bool is_attached() const;

protected: //observer interface

	virtual
	void attach_impl( core::pose::Pose & pose );

	virtual
	void detach_impl();

	void
	on_length_change( core::conformation::signals::LengthEvent const & event );

private:

	utility::signals::Link length_event_link_;
	conformation::carbohydrates::GlycanTreeSetOP glycan_tree_set_;

};


} //core
} //pose
} //carbohydrates



#endif //INCLUDED_core_pose_carbohydrates_GlycanTreeSetObserver_hh





