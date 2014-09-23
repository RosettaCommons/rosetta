// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   TemplateHistory.hh
/// @brief
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_hybridization_TemplateHistory_hh
#define INCLUDED_protocols_hybridization_TemplateHistory_hh

#include <protocols/hybridization/TemplateHistory.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/CacheableData.hh>
#include <core/types.hh>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

///////////////////////////////////////////////////////////////////////////////
// ncs residue mapping 
//   - stored in the pose and used by other movers (fragment insertion for example)
class TemplateHistory:  public basic::datacache::CacheableData {
public:
	TemplateHistory( core::pose::Pose &pose );

	basic::datacache::CacheableDataOP clone() const {
		return basic::datacache::CacheableDataOP( new TemplateHistory(*this) );
	}

	void setall( int template_id );
	void set( core::Size res_start, core::Size res_stop, int template_id );
	int get( core::Size resid );
	core::Size size() { return history_.size(); }

private:
	utility::vector1< int > history_;
};



} // symmetry
//} // simple_moves
} // rosetta
#endif
