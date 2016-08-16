// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/CachedResidueSubset.hh
/// @brief
/// @author Tom Linsky ( tlinsky at uw dot edu )

#ifndef INCLUDED_core_select_residue_selector_CachedResidueSubset_hh
#define INCLUDED_core_select_residue_selector_CachedResidueSubset_hh

#include <core/select/residue_selector/CachedResidueSubset.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <map>
#include <vector>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

class CachedResidueSubset: public basic::datacache::CacheableData {
public:
	typedef std::vector< bool > BoolVector;
	typedef utility::pointer::shared_ptr< BoolVector > BoolVectorOP;
	typedef std::map< std::string, BoolVectorOP > ResidueSubsetMap;

public:
	// constructors
	CachedResidueSubset();
	virtual basic::datacache::CacheableDataOP clone() const;
	virtual basic::datacache::CacheableDataOP fresh_instance() const;

	void set_subset( ResidueSubsetCOP const subset, std::string const & name );
	ResidueSubsetCOP get_subset( std::string const & name ) const;

	bool has_subset( std::string const & name ) const;

private:
	ResidueSubsetMap subsets_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // residue_selector
} // select
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_select_residue_selector_CachedResidueSubset )
#endif // SERIALIZATION


#endif
