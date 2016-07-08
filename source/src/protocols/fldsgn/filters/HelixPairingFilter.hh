// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/filters/HelixPairingFilter.hh
/// @brief header file for HelixPairingFilter class.
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_HelixPairingFilter_hh
#define INCLUDED_protocols_fldsgn_filters_HelixPairingFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/HelixPairingFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>
#include <protocols/fldsgn/topology/HelixPairing.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


// C++ headers
#include <set>

namespace protocols {
namespace fldsgn {
namespace filters {

class HelixPairingFilter : public protocols::filters::Filter {
public:


	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef std::string String;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::pose::Pose Pose;
	typedef protocols::fldsgn::topology::HelixPairing  HelixPairing;
	typedef protocols::fldsgn::topology::HelixPairings  HelixPairings;
	typedef protocols::fldsgn::topology::HelixPairingSet  HelixPairingSet;
	typedef protocols::fldsgn::topology::HelixPairingSetOP  HelixPairingSetOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	HelixPairingFilter();

	// @brief constructor with arguments
	HelixPairingFilter( String const & hf );

	// @brief constructor with arguments
	HelixPairingFilter( HelixPairings const & hpairs );

	// @brief copy constructor
	HelixPairingFilter( HelixPairingFilter const & rval );

	virtual ~HelixPairingFilter(){}


public:// virtual constructor


	// @brief make clone
	virtual FilterOP clone() const { return FilterOP( new HelixPairingFilter( *this ) ); }

	// @brief make fresh instance
	virtual FilterOP fresh_instance() const { return FilterOP( new HelixPairingFilter() ); }


public:// mutator


	void helix_pairings( String const & hpairs );

	void helix_pairings( HelixPairings const & hpairs );

	void secstruct( String const & ss );

	/// @brief sets distance cutoff
	void dist( Real const dist_val );

	/// @brief sets max helix bend
	void bend_angle( Real const bend_angle_val );

	/// @brief sets max cross angle
	void cross_angle( Real const cross_angle_val );

	/// @brief sets max alignment angle
	void align_angle( Real const align_angle_val );


public:// accessor


	// @brief get name of this filter
	virtual std::string name() const { return "HelixPairingFilter"; }


public:// parser


	virtual void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & );


public:// virtual main operation


	/// @brief
	Real report_sm( Pose const & pose ) const;

	/// @brief
	Real compute( Pose const & pose ) const;


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	virtual bool apply( Pose const & pose ) const;

private:
	/// @brief returns secondary structure to be used for finding helices
	/// @details If secstruct_ is set, returns that.
	///          If use_dssp_ is true, returns secstruct from DSSP
	///          Otherwise, returns the pose.secstruct()
	std::string
	secstruct( core::pose::Pose const & pose ) const;

	/// @brief Returns the Helix pairing set to filter based on.
	/// @details If hpairset_ is set, returns that.
	///          Otherwise, look up the pairing set using the StructureData cached in the pose
	topology::HelixPairingSet
	compute_helix_pairing_set( core::pose::Pose const & pose ) const;

	/// @brief Returns set of helix IDs that were not found in the pose.
	///        Set is empty if all helices are present.
	std::set< core::Size >
	find_missing_helices(
		topology::SS_Info2 const & ss_info,
		topology::HelixPairings const & hpairs ) const;

private:

	/// @brief if value is empty, dssp will run for ss definition ( default is emptry )
	String secstruct_;

	// TL: If this filter was used for two different poses sequentially, the DSSP secstruct from
	//     pose1 would be reused for pose2, because secstruct_ was being set in apply() if empty.
	//     Much better solution is to keep secstruct_ const.  Mutables should almost never be used.
	//mutable String secstruct_;

	/// @brief
	Real dist_cutoff_;

	/// @brief
	Real bend_angle_;

	/// @brief
	Real cross_angle_;

	/// @brief
	Real align_angle_;

	/// @brief helix pairings
	HelixPairingSetOP hpairset_;

	/// @brief HelixPairing id for output
	Size output_id_;

	/// @brief output type, dist or angle
	String output_type_;

	/// @brief if true, and secstruct_ is empty, DSSP will be used to determine secondary structure of the pose
	bool use_dssp_;

};

} // filters
} // fldsgn
} // protocols

#endif
