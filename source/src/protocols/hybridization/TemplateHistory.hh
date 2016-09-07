// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   TemplateHistory.hh
/// @brief
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_hybridization_TemplateHistory_hh
#define INCLUDED_protocols_hybridization_TemplateHistory_hh

#include <protocols/hybridization/TemplateHistory.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/CacheableData.hh>
#include <core/types.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

///////////////////////////////////////////////////////////////////////////////
// ncs residue mapping
//   - stored in the pose and used by other movers (fragment insertion for example)
class TemplateHistory:  public basic::datacache::CacheableData {
public:
	TemplateHistory( core::pose::Pose &pose );

	basic::datacache::CacheableDataOP clone() const override {
		return basic::datacache::CacheableDataOP( new TemplateHistory(*this) );
	}

	void setall( int template_id );
	void set( core::Size res_start, core::Size res_stop, int template_id );
	int get( core::Size resid );
	core::Size size() { return history_.size(); }

private:
	utility::vector1< int > history_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	TemplateHistory();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // symmetry
//} // simple_moves
} // rosetta
#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_hybridization_TemplateHistory )
#endif // SERIALIZATION


#endif
