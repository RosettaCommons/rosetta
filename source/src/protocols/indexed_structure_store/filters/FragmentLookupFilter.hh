// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/filters/FragmentLookupFilter.hh
/// @brief Header for FragmentLookupFilter class.
/// @details
/// @author Alex Ford (fordas@uw.edu) and Daniel Adriano Silva (dadriano@uw.edu)


#ifndef INCLUDED_protocols_indexed_structure_store_filters_FragmentLookupFilter_hh
#define INCLUDED_protocols_indexed_structure_store_filters_FragmentLookupFilter_hh

#include <protocols/indexed_structure_store/filters/FragmentLookupFilter.fwd.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
//#include <core/indexed_structure_store/FragmentLookup.fwd.hh>
#include <core/indexed_structure_store/FragmentLookup.hh>
#include <core/indexed_structure_store/FragmentStore.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>

#include <string>

namespace protocols {
namespace indexed_structure_store {
namespace filters {

// @brief Lookup mode used by filter.
//
// First - Returns the first fragment matching the lookup criteria.  Faster
// but not guarenteed to identify the closest potential library fragment.
//
// Closest - Returns the closest fragment with the target lookup.
// Significantly slower.
enum LookupMode {First, Closest};


// @brief Lookup pose sub-fragments within a target fragment library, filtering
// pose which contain fragments not present in the library.
//
// Targets a specified fragment lookup, which defines a fragment specification
// and store of candidate fragments. Filter extracts all fragments from the
// pose via the specification passes these fragments to the lookup. Filter
// fails if any pose fragment can not be found within the store.
//
// Lookup results, including a lookup success/failure flag and lookup
// confidence score for successful lookups cached within the filter after each
// apply call and exposed via the FragmentLookupFilter::lookup_result member.
//
// Fragment stores are identified by a lookup name and a target store path. The
// default store defined in the StructureStoreManager class via the options
// system.
class FragmentLookupFilter : public protocols::filters::Filter {

public:
	typedef core::indexed_structure_store::FragmentLookupOP FragmentLookupOP;
	typedef core::indexed_structure_store::FragmentLookupResult FragmentLookupResult;
	typedef core::indexed_structure_store::FragmentSpecification FragmentSpecification;
	typedef core::Real Real;
	typedef core::Size Size;

	typedef core::pose::Pose Pose;

	typedef utility::tag::TagPtr TagPtr;
	typedef protocols::filters::Filters_map Filters_map;
	typedef protocols::filters::FilterOP FilterOP;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;

	friend class FragmentLookupFilterCreator;

public:
	/// @brief Construct targeting the given named lookup in the default backend.
	//
	// lookup_name - Lookup name within store. See StructureStoreManager.
	// mode - Filter mode, see LookupMode.
	FragmentLookupFilter(std::string lookup_name,
		LookupMode mode=First);

	/// @brief Construct targeting the given lookup in the given backend.
	//
	// lookup_name - Lookup name within store. See StructureStoreManager.
	// store_path - Target store path. HDF5 file path (extension is optional) if
	//     HDF5 support is enabled. Root directory of binary store if not HDF5
	//     support.
	// mode - Filter mode, see LookupMode.
	FragmentLookupFilter(std::string lookup_name,
		std::string store_path,
		LookupMode mode=First,
		core::Size target_chain=0,
		core::Size threshold=0,
		bool b_target_chain=false);

	// @brief copy constructor
	FragmentLookupFilter( FragmentLookupFilter const & rval );

	virtual ~FragmentLookupFilter(){}


public:// virtual constructor

	// @brief make clone
	FilterOP clone() const override { return FilterOP( new FragmentLookupFilter( *this ) ); }

	// @brief make fresh instance
	FilterOP fresh_instance() const override { return FilterOP( new FragmentLookupFilter() ); }

public:

	// @brief Filter Name
	// XRW TEMP  virtual std::string name() const { return "FragmentLookupFilter"; }

	// @brief False if any pose fragment is outside the fragment lookup radius.
	//
	// Cached lookup results are available under
	// FragmentLookupFilter::lookup_result after apply call.
	bool apply( Pose const & pose ) const override;

	// @brief Parse arguments from rosettascripts XML
	void parse_my_tag( utility::tag::TagCOP tag,
		DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & ) override;

	//Main Filter computation routine
	core::Size compute( core::pose::Pose const & pose ) const;

	// @brief Write cached fragment lookup result as json to output stream.
	void report( std::ostream &, core::pose::Pose const & ) const override;

	// @Report back to JD2.
	core::Real report_sm( core::pose::Pose const & pose ) const override;

	// @brief Per-fragment lookup results of previous apply call, keyed by fragment start residue id.
	std::map<core::Size, FragmentLookupResult> const & lookup_result() { return cached_lookup_result_; }

	// @brief Fragment specification from underlying lookup, provides fragment length.
	FragmentSpecification const & fragment_specification();

	/// @brief sets chain
	inline void set_chain( core::Size const chainid ) { target_chain_ = chainid; }

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:
	// @brief Default constructor
	FragmentLookupFilter();

private:

	/// @brief
	FragmentLookupOP target_lookup_;
	mutable std::map<core::Size, FragmentLookupResult> cached_lookup_result_;

	LookupMode lookup_mode_;
	core::Size target_chain_;
	core::Size threshold_;
	bool b_target_chain_;
};

}
}
}

#endif
