// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/filters/SecondaryStructureFilter.hh
/// @brief header file for SecondaryStructureFilter class.
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_SecondaryStructureFilter_hh
#define INCLUDED_protocols_fldsgn_filters_SecondaryStructureFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/SecondaryStructureFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace fldsgn {
namespace filters {

class SecondaryStructureFilter : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef std::string String;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	SecondaryStructureFilter();

	// @brief constructor with arguments
	SecondaryStructureFilter( String const & ss );

	virtual ~SecondaryStructureFilter(){}


public:// virtual constructor


	// @brief make clone
	virtual FilterOP clone() const { return FilterOP( new SecondaryStructureFilter( *this ) ); }

	// @brief make fresh instance
	virtual FilterOP fresh_instance() const { return FilterOP( new SecondaryStructureFilter() ); }


public:// mutator


	// @brief set filtered secondary structure
	void filtered_ss( String const & s );

	// @brief set filtered abego
	void filtered_abego( String const & s );

	/// @brief Should we use the secstruct in the pose (true), or recompute via dssp (false)? default=false
	void set_use_pose_secstruct( bool const use_ss );

public:// accessor


	// @brief get name of this filter
	virtual std::string name() const { return "SecondaryStructureFilter"; }

	/// @brief sets the blueprint file based on filename.  If a strand pairing is impossible (i.e. the structure has two strands, 5 and 6 residues, respectively, it sets the unpaired residues to 'h' so that they still match.
	void set_blueprint( std::string const & blueprint_file );

	/// @brief sets the residue selector used to choose which residues to scan.  Default is all protein residues
	void
	set_residue_selector( core::select::residue_selector::ResidueSelector const & selector );

public:// parser

	virtual void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & );


public:// virtual main operation


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	virtual bool apply( Pose const & pose ) const;
	virtual core::Real report_sm( Pose const & pose ) const;

private:
	/// @brief Returns the desired secondary structure
	/// @details  Rules for selecting this seconary structure:
	///           1. If a user-specified filtered_ss_ is set, return that
	///           2. If StructureData is cached in the pose, use the secondary structure of that
	///           3. Use pose secondary structure, throwing and error if use_pose_secstruct is true
	std::string
	get_filtered_secstruct( core::pose::Pose const & pose ) const;

	/// @brief Computes number of residues matching the desired secondary structure
	/// @param[in]  pose   The pose to scan
	/// @param[in]  subset The residue subset of the pose to include in counting. Must match
	///                    the pose size.
	/// @returns  Number of protein residues in the residue subset that match the desired
	///           secondary structure
	core::Size
	compute( Pose const & pose, core::select::residue_selector::ResidueSubset const & subset ) const;

private:
	core::select::residue_selector::ResidueSelectorCOP selector_;

	String filtered_ss_;

	utility::vector1< String > filtered_abego_;

	bool use_abego_;

	bool use_pose_secstruct_;

	core::Real threshold_;
};

} // filters
} // fldsgn
} // protocols

#endif
