// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/CstMinimizationData.hh
/// @brief  A cacheable data wrapper for ConstraintsOPs for use in minimization
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_constraints_CstMinimizationData_hh
#define INCLUDED_core_scoring_constraints_CstMinimizationData_hh

// Package headers
#include <core/scoring/constraints/Constraints.fwd.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace constraints {

class CstMinimizationData;
typedef utility::pointer::shared_ptr< CstMinimizationData > CstMinimizationDataOP;
typedef utility::pointer::shared_ptr< CstMinimizationData const > CstMinimizationDataCOP;

class CstMinimizationData : public basic::datacache::CacheableData
{
public:
	typedef basic::datacache::CacheableDataOP  CacheableDataOP;
	typedef basic::datacache::CacheableDataCOP CacheableDataCOP;

public:
	CstMinimizationData();
	CstMinimizationData( ConstraintsOP constraints );
	virtual ~CstMinimizationData();
	virtual CacheableDataOP clone() const;

	Constraints const & constraints() const { return *constraints_; }
	void set_constraints( ConstraintsOP constraints );

private:
	ConstraintsOP constraints_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_CstMinimizationData )
#endif // SERIALIZATION


#endif
