// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/filters/SheetTopologyFilter.hh
/// @brief header file for SheetTopologyFilter class.
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_SheetTopologyFilter_hh
#define INCLUDED_protocols_fldsgn_filters_SheetTopologyFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/SheetTopologyFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/fldsgn/topology/StrandPairing.fwd.hh>

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

class SheetTopologyFilter : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef std::string String;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::pose::Pose Pose;
	typedef protocols::fldsgn::topology::StrandPairingSet StrandPairingSet;
	typedef protocols::fldsgn::topology::StrandPairingSetOP StrandPairingSetOP;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::SS_Info2_OP SS_Info2_OP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;

	typedef std::pair< core::Size, core::Size > ResiduePairing;
	typedef std::set< ResiduePairing > ResiduePairingSet;
	typedef utility::vector1< ResiduePairingSet > ResiduePairingSets;

public:// constructor/destructor


	// @brief default constructor
	SheetTopologyFilter();

	// @brief constructor with arguments
	SheetTopologyFilter( StrandPairingSetOP const & sps );

	// @brief constructor with arguments
	SheetTopologyFilter( String const & sheet_topology );

	// @brief copy constructor
	SheetTopologyFilter( SheetTopologyFilter const & rval );

	virtual ~SheetTopologyFilter(){}


public:// virtual constructor


	// @brief make clone
	FilterOP clone() const override { return FilterOP( new SheetTopologyFilter( *this ) ); }

	// @brief make fresh instance
	FilterOP fresh_instance() const override { return FilterOP( new SheetTopologyFilter() ); }


public:// mutator


	// @brief set filtered sheet_topology by StrandPairingSetOP
	void filtered_sheet_topology( StrandPairingSetOP const & sps );

	// @brief set filtered sheet_topology by string
	void filtered_sheet_topology( String const & sheet_topology );

	/// @brief set user-specified pose secondary structure
	void set_secstruct( std::string const & ss );

	/// @brief if true, and secstruct is unset, dssp is used on the input.  Otherwise, the pose.secstruct() is used
	void set_use_dssp( bool const use_dssp );

public:// accessor

	// @brief get name of this filter
	// XRW TEMP  virtual std::string name() const { return "SheetTopologyFilter"; }


public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & ) override;


public:// virtual main operation

	/// @brief returns the fraction of pairings that pass the filter
	core::Real compute( Pose const & pose ) const;

	/// @brief returns the fraction of pairings that pass
	core::Real report_sm( Pose const & pose ) const override;

	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	bool apply( Pose const & pose ) const override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	/// @brief Returns the desired strand pairing topology string
	/// @details  Rules for selecting this topology string:
	///           1. If a user-specified filtered_sheet_topology_ is set, return that
	///           2. If StructureData is cached in the pose, determine pairings from that
	///           3. throw error
	std::string
	get_filtered_sheet_topology( core::pose::Pose const & pose ) const;

	/// @brief Returns the pose secondary structure to be used in computation
	/// @details  Rules for selecting the secondary structure:
	///           1. If a user-specified secstruct_input_ is set, return this
	///           2. If use_dssp is true, determine secondary structure by DSSP
	///           3. Return pose secondary stucture otherwise
	std::string
	get_secstruct( core::pose::Pose const & pose ) const;

	/// @brief Given the filtered strand pairings, compute the number of residue pairings possible
	/// @param[in] spairset The strand pairing set to be used to find residue pairings.  It is
	///                     non-const because the pairings are stored as OPs, so begin() and end()
	///                     are non-const
	/// @param[in] ss_info  SS_Info2 object describing the secondary structure of the pose
	/// @returns Vector of ResiduePairingSets, one for each strand pairing, in the same order as
	///          in spairset
	ResiduePairingSets
	compute_residue_pairings(
		topology::StrandPairingSet & spairset,
		topology::SS_Info2 const & ss_info ) const;

	/// @brief Computes number of pairings in the given StrandPairing
	/// @param[in] pairing  StrandPairing which contains residue pairing information
	/// @param[in] ss_info  SS_Info2 object describing the secondary structure of the pose
	/// @returns ResiduePairingSet containing pairs of residues
	ResiduePairingSet
	compute_paired_residues(
		topology::StrandPairing const & pairing,
		topology::SS_Info2 const & ss_info ) const;

	/// @brief Counts total number of residue pairings present in the ResiduePairingSets
	core::Size
	count_residue_pairings( ResiduePairingSets const & pair_sets ) const;

	/// @brief Counts number of residue pairs in the filtered_pair_set are present in the pose_pair_set
	core::Size
	count_good_pairings(
		ResiduePairingSet const & filtered_pair_set,
		ResiduePairingSet const & pose_pair_set ) const;

	/// @brief Replace register shift of pairings in pose_spairset with 99 if register shift in filtered_spairset
	///        is 99
	void
	replace_register_shifts(
		topology::StrandPairingSet & spairset,
		topology::StrandPairingSet & filtered_spairset ) const;

private:

	String filtered_sheet_topology_;

	String secstruct_input_;

	bool ignore_register_shift_;

	bool use_dssp_;

};

/// @brief Searches the StrandPairingSet for a pairing containing s1 and s2. Returns OP to it
topology::StrandPairingOP
find_pairing( topology::StrandPairingSet & spairset, core::Size const s1, core::Size const s2 );

/// @brief Searches the StrandPairingSet for a pairing containing s1 and s2. Returns its 1-based index
core::Size
find_pairing_idx( topology::StrandPairingSet & spairset, core::Size const s1, core::Size const s2 );

/// @brief helper function for replacing register shift of all pairs with 99
std::string
remove_register_shifts( std::string const & pair_str );

} // filters
} // fldsgn
} // protocols

#endif
