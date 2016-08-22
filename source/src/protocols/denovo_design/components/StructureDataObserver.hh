// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/denovo_design/components/StructureDataObserver.hh
/// @brief  Observes a pose and updates StructureData accordingly
///
/// @author Tom Linsky (tlinsky at uw dot edu )

#ifndef INCLUDED_protocols_denovo_design_components_StructureDataObserver_hh
#define INCLUDED_protocols_denovo_design_components_StructureDataObserver_hh

// Unit headers
#include <protocols/denovo_design/components/StructureDataObserver.fwd.hh>
#include <core/pose/datacache/CacheableObserver.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.hh>

// Core headers
#include <core/conformation/signals/LengthEvent.hh>

// Utility headers
#include <utility/signals/Link.hh>
#include <utility/vector1.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace denovo_design {
namespace components {

/// @brief a cacheable observer that keeps track of what length events occured
class StructureDataObserver : public core::pose::datacache::CacheableObserver {

public: // construct/destruct

	// public for serialization only!  Don't construct StructureDataObservers without
	// a StructureData object!!
	StructureDataObserver();

	StructureDataObserver( StructureDataOP sd );

	/// @brief default destructor
	/// @remarks detaches during destruction
	virtual
	~StructureDataObserver();


public: // virtual constructors

	/// @brief clone this object
	virtual core::pose::datacache::CacheableObserverOP
	clone();

	/// @brief create a new instance of this object
	virtual core::pose::datacache::CacheableObserverOP
	create();


public: // interface

	/// @brief is this observer attached to a Pose/Conformation?
	virtual bool
	is_attached() const;

protected: // virtual observer interface

	/// @brief attach to Pose/Conformation
	virtual void
	attach_impl( core::pose::Pose & pose );


	/// @brief detach from Pose/Conformation
	virtual void
	detach_impl();

	void
	on_length_change( core::conformation::signals::LengthEvent const & event );

private:
	void
	on_residue_delete( core::conformation::signals::LengthEvent const & event );

	void
	on_residue_append( core::conformation::signals::LengthEvent const & event );

	void
	on_residue_prepend( core::conformation::signals::LengthEvent const & event );

private:
	StructureDataOP sd_;
	utility::signals::Link length_event_link_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace components
} // namespace denovo_design
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_denovo_design_components_StructureDataObserver )
#endif // SERIALIZATION

#endif // INCLUDED_protocols_denovo_design_components_StructuerDataObserver_hh
