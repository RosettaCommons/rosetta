// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/LongestContinuousPolarSegmentFilter.hh
/// @brief This filter computes the longest continuous stretch of polar residues within a pose or selection.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_LongestContinuousPolarSegmentFilter_hh
#define INCLUDED_protocols_simple_filters_LongestContinuousPolarSegmentFilter_hh

// Unit headers
#include <protocols/simple_filters/LongestContinuousPolarSegmentFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh

namespace protocols {
namespace simple_filters {

///@brief This filter computes the longest continuous stretch of polar residues within a pose or selection.
class LongestContinuousPolarSegmentFilter : public protocols::filters::Filter {

public:
	/// @brief Constructor.
	/// @details this constructor sets count_gly_as_polar to true
	LongestContinuousPolarSegmentFilter();

	/// @brief Constructor.
	/// @details this constructor sets count_gly_as_polar to false
	LongestContinuousPolarSegmentFilter( std::string const &filter_name );

	/// @brief destructor (important for properly forward-declaring smart-pointer members)
	///
	~LongestContinuousPolarSegmentFilter() override;

	/// @brief returns true if the structure passes the filter, false otherwise
	bool
	apply( core::pose::Pose const & pose ) const override;

	/// @brief required for reporting score values
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

	/// @brief allows printing data to a stream
	void
	report( std::ostream & os, core::pose::Pose const & pose ) const override;

public: //Setters and getters

	/// @brief Set whether I should exclude stretches of polars are the N- and C-termini of chains.
	inline void set_exclude_chain_termini( bool const setting ) { exclude_chain_termini_=setting; }

	/// @brief Get whether I should exclude stretches of polars are the N- and C-termini of chains.
	inline bool exclude_chain_termini() const { return exclude_chain_termini_; }

	/// @brief Set whether glycine is counted as a polar residue type.
	inline void set_count_gly_as_polar( bool const setting ) { count_gly_as_polar_ = setting; }

	/// @brief Get whether glycine is counted as a polar residue type.
	inline bool count_gly_as_polar() const { return count_gly_as_polar_; }

	/// @brief Set whether I should filter out high (true) or low (false) poses.
	inline void set_filter_out_high( bool const setting ) { filter_out_high_=setting; }

	/// @brief Get whether I should filter out high (true) or low (false) poses.
	inline bool filter_out_high() const { return filter_out_high_; }

	/// @brief Set the max (or min) tolerated number of polars.
	inline void set_cutoff( core::Size const cutoff ) { cutoff_ = cutoff; }

	/// @brief Get the max (or min) tolerated number of polars.
	inline core::Size cutoff() const { return cutoff_; }

	/// @brief Set the residue selector.
	/// @details Does not clone the input; uses it directly.  (Can be shared, or modified later).
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Get the residue selector.
	/// @details Can be shared, or modified later -- this is a true owning pointer to the stored selector, not to a clone.
	inline core::select::residue_selector::ResidueSelectorCOP residue_selector() const { return residue_selector_; };

public:
	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief parse XML tag (to use this Filter in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

	/// @brief Given a residue type, determine whether it's one of the types that this filter
	/// should count.
	/// @details Based on whether the residue type has the POLAR property.  Special-case exception
	/// is made for glycine, depending on whether count_gly_as_polar_ is true.
	virtual bool is_counted( core::chemical::ResidueType const &restype ) const;

	/// @brief returns type of counted residues (polar)
	virtual std::string counted_residue_description() const;

private: //Functions

	/// @brief Given a pose, actually compute the number of residues in the longest stretch, and return this value,
	/// along with the start and end indices (in Rosetta numbering).
	/// @details Returns 0, 0, 0 if no polar stretch could be found.
	/// @param[in] pose The pose to analyse.
	/// @param[in] selector An optional const-owning pointer to a ResidueSelector. If provided, only those stretches that have at least one residue selected by the ResidueSelector
	/// will be counted.
	/// @param[in] ignore_termini If true, stretches at the N- and C-termini of chains are ignored.  If false, they're counted.
	/// @param[out] longest_stretch_res_count The number of polar residues in the longest polar stretch.  Set to 0 if no polar stretch is found.
	/// @param[out] longest_stretch_start The index, in Rosetta numbering, of the first residue of the longest polar stretch.  Set to 0 if no polar stretch is found.
	/// @param[out] longest_stretch_end The index, in Rosetta numbering, of the last residue of the longest polar stretch.  Set to 0 if no polar stretch is found.
	void compute(
		core::pose::Pose const &pose,
		core::select::residue_selector::ResidueSelectorCOP selector,
		bool const ignore_termini,
		core::Size &longest_stretch_res_count,
		core::Size &longest_stretch_start,
		core::Size &longest_stretch_end
	) const;

private: //Variables

	/// @brief Should I exclude stretches of polars are the N- and C-termini of chains?
	/// @details Default true.
	bool exclude_chain_termini_;

	/// @brief If true, glycine will be considered "polar".  True by default.
	///
	bool count_gly_as_polar_;

	/// @brief Should I filter out designs with too many or too few polars
	/// in a stretch?
	/// @details Default is true (filter out designs with too many).
	bool filter_out_high_;

	/// @brief The cutoff value for filtering.
	/// @details Above (or below, if filtering out low) this value, designs are rejected.  If
	/// the number of polars is equal to or below this value (or equal to or above, if filtering
	/// out low), designs are accepted.
	core::Size cutoff_;

	/// @brief An optional residue selector to select a subset of a pose.
	/// @details Will be nullptr if not specified.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

};

} //protocols
} //simple_filters

#endif //INCLUDED_protocols_simple_filters_LongestContinuousPolarSegmentFilter_hh
