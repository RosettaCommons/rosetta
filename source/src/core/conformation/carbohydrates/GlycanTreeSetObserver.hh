// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/carbohydrates/GlycanTreeSetObserver.hh
/// @brief The CacheablePoseObserver version of GlycanTreeSet that will react to pose length changes..
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_core_conformation_carbohydrates_GlycanTreeSetObserver_hh
#define INCLUDED_core_conformation_carbohydrates_GlycanTreeSetObserver_hh

#include <core/conformation/carbohydrates/GlycanTreeSetObserver.fwd.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.fwd.hh>


#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/signals/Link.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace conformation {
namespace carbohydrates {


///@brief The CacheablePoseObserver version of GlycanTreeSet that will react to pose length changes..
class GlycanTreeSetObserver : public utility::pointer::ReferenceCount {

public:

	GlycanTreeSetObserver();

	///@brief Construct the GlycanTreeSet, but do not attach it to any pose.
	GlycanTreeSetObserver( Conformation const & conf );

	///@Brief Construct the GlycanTreeSet and attach this object to the pose.
	//GlycanTreeSetObserver( core::pose::Pose & pose );

	GlycanTreeSetObserver(GlycanTreeSetObserver const & src);

	virtual ~GlycanTreeSetObserver();

	GlycanTreeSetObserverOP
	clone() const;

public:

	///@brief Get the GlycanTreeSet that is maintained by this Observer.
	GlycanTreeSetCOP
	get_glycan_tree_set() const;


public: //observer interface

	bool
	is_attached() const;

public: //observer interface

	/// @brief Detach and attach to Conformation
	void
	attach_to( Conformation & conf );
	
	/// @brief Do the attachment to the length event signal
	void
	attach_impl( Conformation & conf );

	void
	detach_impl();

	void
	on_length_change( signals::LengthEvent const & event );

private:

	utility::signals::Link length_event_link_;
	GlycanTreeSetOP glycan_tree_set_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


}; // GlycanTreeSetObserver



} //core
} //conformation
} //carbohydrates

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_carbohydrates_GlycanTreeSetObserver )
#endif // SERIALIZATION

#endif //INCLUDED_core_conformation_carbohydrates_GlycanTreeSetObserver_hh





