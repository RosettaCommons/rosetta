// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/MinimizationData.hh
/// @brief  A container class for use by certain EnergyMethods during derivative and
//          score function evaluation during minimization routines.
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_MinimizationData_hh
#define INCLUDED_core_scoring_MinimizationData_hh

// Unit headers
#include <core/scoring/MinimizationData.fwd.hh>

// Project headers
#include <basic/datacache/CacheableData.fwd.hh>
#ifdef WIN32
#include <basic/datacache/CacheableData.hh>
#endif

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {

enum min_single_data{
	etab_single_nblist = 1,
	etab_classic_intrares_single_nblist,
	mm_lj_intra_nblist,
	cst_res_data,
	lkb_res_data,
	vdw_res_data,
	mp_res_data,
	hbond_res_data,
	n_min_single_data = hbond_res_data // keep this guy last
};

enum min_pair_data {
	etab_pair_nblist = 1,
	etab_classic_intrares_pair_nblist,
	cst_respair_data,
	elec_pair_nblist,
	geom_solv_pair_nblist,
	lk_PolarNonPolar_pair_nblist,
	fa_dslf_respair_data,
	fa_custom_pair_dist_data,
	lkb_respair_data,
	vdw_respair_data,
	mp_respair_data,
	mg_pair_nblist,
	hbond_respair_data,
	n_min_pair_data = hbond_respair_data // keep this guy last
};


class ResSingleMinimizationData : public utility::pointer::ReferenceCount
{
public:
	typedef basic::datacache::CacheableData    CacheableData;
	typedef basic::datacache::CacheableDataOP  CacheableDataOP;
	typedef basic::datacache::CacheableDataCOP CacheableDataCOP;

public:
	ResSingleMinimizationData();
	virtual ~ResSingleMinimizationData();
	ResSingleMinimizationData( ResSingleMinimizationData const & ); // deep copy
	ResSingleMinimizationData & operator = ( ResSingleMinimizationData const & ); // deep copy

	void set_data( min_single_data index, CacheableDataOP data );
	CacheableDataOP get_data( min_single_data index );
	CacheableDataCOP get_data( min_single_data index ) const;
	CacheableData & get_data_ref( min_single_data index ) { return * data_cache_[ index ]; }
	CacheableData const & get_data_ref( min_single_data index ) const { return * data_cache_[ index ]; }

private:

private:
	utility::vector1< CacheableDataOP > data_cache_;

};

class ResPairMinimizationData : public utility::pointer::ReferenceCount
{
public:
	typedef basic::datacache::CacheableData    CacheableData;
	typedef basic::datacache::CacheableDataOP  CacheableDataOP;
	typedef basic::datacache::CacheableDataCOP CacheableDataCOP;

public:
	ResPairMinimizationData();
	virtual ~ResPairMinimizationData();
	ResPairMinimizationData( ResPairMinimizationData const & ); // deep copy
	ResPairMinimizationData & operator = ( ResPairMinimizationData const & ); // deep copy

	void set_data( min_pair_data index, CacheableDataOP );
	CacheableDataOP get_data( min_pair_data index );
	CacheableDataCOP get_data( min_pair_data index ) const;
	CacheableData & get_data_ref( min_pair_data index ) { return * data_cache_[ index ]; }
	CacheableData const & get_data_ref( min_pair_data index ) const { return * data_cache_[ index ]; }


private:
	utility::vector1< CacheableDataOP > data_cache_;
};

}
}

#endif
