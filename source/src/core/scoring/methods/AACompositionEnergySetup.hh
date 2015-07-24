// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/methods/AACompositionEnergySetupSetup.hh
/// @brief Headers for a helper for the EnergyMethod that penalizes deviation from a desired amino acid composition.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).



#ifndef INCLUDED_core_scoring_methods_AACompositionEnergySetupSetup_hh
#define INCLUDED_core_scoring_methods_AACompositionEnergySetupSetup_hh

// Unit headers
#include <core/scoring/methods/AACompositionEnergySetup.fwd.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueProperty.hh>

// Project headers
#include <core/types.hh>
#include <map>
#include <string>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// @brief AACompositionPropertiesSet, a helper class that stores a set of properties that MUST
/// be present for a residue to be counted, and a set of properties that MUST NOT be present for
/// that residue to be counted.
class AACompositionPropertiesSet : public utility::pointer::ReferenceCount {
public:

	/// @brief Default constructor for AACompositionEnergySetupPropertiesSet.
	///
	AACompositionPropertiesSet();

	/// @brief Constructor for AACompositionEnergySetupPropertiesSet that takes lists of
	/// included and excluded properties.
	AACompositionPropertiesSet(
		utility::vector1< std::string > const &included_properties_strings,
		utility::vector1< std::string > const &excluded_properties_strings
	);
	
	/// @brief Copy constructor for AACompositionEnergySetupPropertiesSet.
	///
	AACompositionPropertiesSet( AACompositionPropertiesSet const &src );

	/// @brief Default destructor for AACompositionEnergySetupPropertiesSet.
	///
	virtual ~AACompositionPropertiesSet();
	
	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	virtual AACompositionPropertiesSetOP clone() const;

	/// @brief Check whether the properties of a residue type match those of this set.
	/// @details This checks (a) that all of the included properties are in the rsd_type, and
	/// (b) that none of the excluded properties is in the rsd_type.
	inline bool properties_match_this( core::chemical::ResidueType const &rsd_type ) const {
		for(core::Size i=1, imax=included_properties_.size(); i<=imax; ++i) {
			if(!rsd_type.has_property(included_properties_[i])) return false;
		}
		for(core::Size i=1, imax=excluded_properties_.size(); i<=imax; ++i) {
			if(rsd_type.has_property(excluded_properties_[i])) return false;
		}
		return true;
	}
	
	/// @brief Take a list of included property strings and parse it.
	/// @details Populates the included_properties_ vector based on the properties named in
	/// the list.
	void parse_included_properites( utility::vector1< std::string > const &proplist );

	/// @brief Take a list of excluded property strings and parse it.
	/// @details Populates the excluded_properties_ vector based on the properties named in
	/// the list.
	void parse_excluded_properites( utility::vector1< std::string > const &proplist );
	
	/// @brief Generate a one-line summary of the properties stored in this AACompositionPropertySet
	///
	std::string one_line_report() const;

private:
	/******************
	Private functions:
	******************/
	
	/// @brief Check whether a ResidueProperty is in a list.
	///
	inline bool is_in_list( core::chemical::ResidueProperty const property, utility::vector1 < core::chemical::ResidueProperty > const &list) const {
		for(core::Size i=1, imax=list.size(); i<=imax; ++i) {
			if(list[i]==property) return true;
		}
		return false;
	}

	/// @brief Add a property to the list of properties that must be present.
	///
	void add_included_property( core::chemical::ResidueProperty const property );

	/// @brief Add a property to the list of properties that must not be present.
	///
	void add_excluded_property( core::chemical::ResidueProperty const property );
	
	/// @brief Convert a property name to a ResidueProperty enum.
	///
	inline core::chemical::ResidueProperty parse_property( std::string const &name ) const {
		return core::chemical::ResidueProperties::get_property_from_string( name ); //Will throw an error if the property is not recognized.
	}
	
private:
	/******************
	Private variables:
	******************/
	
	/// @brief Properties that a residue MUST have.
	///
	utility::vector1 < core::chemical::ResidueProperty > included_properties_;

	/// @brief Properties that a residue MUST NOT have.
	///
	utility::vector1 < core::chemical::ResidueProperty > excluded_properties_;

}; // class AACompositionPropertiesSet


/// @brief AACompositionEnergySetup, a helper class for the AACompositionEnergy energy method
/// that stores all of its setup data.
class AACompositionEnergySetup : public utility::pointer::ReferenceCount {
public:

	/// @brief Default constructor for AACompositionEnergySetup.
	///
	AACompositionEnergySetup();
	
	/// @brief Copy constructor for AACompositionEnergySetup.
	///
	AACompositionEnergySetup( AACompositionEnergySetup const &src );

	/// @brief Default destructor for AACompositionEnergySetup.
	///
	virtual ~AACompositionEnergySetup();
	
	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	virtual AACompositionEnergySetupOP clone() const;

public:

	/**********************
	Public setup functions:
	**********************/

	/// @brief Reset all data in this data storage object.
	///
	void reset();
	
	/// @brief Initialize from a .comp file.
	///
	void initialize_from_file( std::string const &filename );

public:

	/*************************
	Public accessor functions:
	*************************/

	/// @brief Look up the penalty for a residue type of internal index res_type_index, given a deviation from
	/// the desired count given by delta_expected.
	inline core::Real type_penalty( signed long const delta_expected, core::Size const res_type_index) const {
		if(delta_expected < type_deviation_ranges_[res_type_index].first) return type_penalties_[res_type_index][1];
		if(delta_expected > type_deviation_ranges_[res_type_index].second) return type_penalties_[res_type_index][ type_penalties_[res_type_index].size() ];
		return type_penalties_[res_type_index][ delta_expected - type_deviation_ranges_[res_type_index].first + 1 ];
	}
	
	/// @brief Look up the penalty for a property set of internal index property_set_index, given a deviation from
	/// the desired count given by delta_expected.
	inline core::Real property_penalty( signed long const delta_expected, core::Size const property_set_index) const {
		if(delta_expected < property_deviation_ranges_[property_set_index].first) return property_penalties_[property_set_index][1];
		if(delta_expected > property_deviation_ranges_[property_set_index].second) return property_penalties_[property_set_index][ property_penalties_[property_set_index].size() ];
		return property_penalties_[property_set_index][ delta_expected - property_deviation_ranges_[property_set_index].first + 1 ];
	}
	
	/// @brief Get the number of types of residues that we'll be counting.
	///
	inline core::Size n_residue_types() const { return res_type_index_mappings_.size(); }
	
	/// @brief Get the number of sets of properties that we'll be counting.
	///
	inline core::Size n_property_sets() const { return property_sets_.size(); }
	
	/// @brief Return true if and only if a particular residue type name3 is in the map of residues to count.
	///
	inline bool has_type( std::string const &name ) const { return ( res_type_index_mappings_.count( name ) != 0 ); }
	
	/// @brief Get the internal index of a particular residue type that we're counting,
	/// given its name3 string.
	/// @details Does not check for whether the type exists in the map first!  Use has_type() to check before calling!
	inline core::Size res_type_index(std::string const &name) const { return res_type_index_mappings_.at( name ); }
	
	/// @brief Get the expected absolute number of residues for a given type, by type index.
	/// @details Warning!  For speed, there's no check that the index is in range!
	inline signed long expected_by_type_absolute( core::Size const index ) const { return expected_by_type_absolute_[index]; }
	
	/// @brief Get the expected absolute number of residues for a given set of properties, by property set index.
	/// @details Warning!  For speed, there's no check that the index is in range!
	inline signed long expected_by_properties_absolute( core::Size const index ) const { return expected_by_properties_absolute_[index]; }
	
	/// @brief Get the expected fractional number of residues for a given type, by type index.
	/// @details Warning!  For speed, there's no check that the index is in range!
	inline core::Real expected_by_type_fraction( core::Size const index ) const { return expected_by_type_fraction_[index]; }
	
	/// @brief Get the expected fractional number of residues for a given set of properties, by property set index.
	/// @details Warning!  For speed, there's no check that the index is in range!
	inline core::Real expected_by_properties_fraction( core::Size const index ) const { return expected_by_properties_fraction_[index]; }
	
	/// @brief Get the indices of all of the property sets corresponding to the current residue.
	/// @details Returns empty vector if none of the property sets that we're counting corresponds to the current residue.
	inline void property_set_indices_matching_residue( core::chemical::ResidueType const &rsd_type, utility::vector1 < core::Size > &indices_out ) const {
		indices_out.clear();
		for(core::Size i=1, imax=property_sets_.size(); i<=imax; ++i) {
			if(property_sets_[i]->properties_match_this( rsd_type ) ) indices_out.push_back(i);
		}
		return;
	}

	/// @brief Get a summary of the data stored in this object
	///
	std::string report() const;

private:

	/******************
	Private functions:
	******************/
	
	/// @brief Parse out penalty definition blocks from the data read from file.
	///
	void parse_penalty_definitions( utility::vector1 < std::string > const &lines );

	/// @brief Parse out penalty definition from a single block of lines from file.
	///
	void parse_a_penalty_definition( utility::vector1 < std::string > const &lines );
	
	/// @brief Do some final checks to ensure that data were loaded properly.
	///
	void check_data() const;

private:

	/******************
	Private variables:
	******************/

	/// @brief A map of the residue types that we will be scoring to their indices.
	/// @details Indicies are internal to this scoring function.
	std::map< std::string, core::Size > res_type_index_mappings_;
	
	/// @brief Penalties for each residue count and each residue type.
	/// @brief Outer vector is the residue type, by internal index.  Inner vector is
	/// the penalty value for a given deviation from the ideal count, by internal index.
	utility::vector1 < utility::vector1 < core::Real > > type_penalties_;
	
	/// @brief Deviation ranges for each residue type.
	/// @brief Outer vector is the residue type, by internal index.  Inner vector is pairs
	/// of deviation ranges (e.g. -10,5 indicating that we are storing penalties for having 10
	/// residues too few up to 5 residues too many for a given type).  Deviations outside of this
	/// range are assigned the value from the end of the range.
	utility::vector1 < std::pair < signed long, signed long > > type_deviation_ranges_;
	
	/// @brief The expected number of residues for each type (by type index), as a fraction (from 0 to 1) of total.
	///
	utility::vector1 < core::Real > expected_by_type_fraction_;

	/// @brief The expected number of residues for each type (by type index), as an absolute number.
	/// @details If this is negative, the fraction is used instead.
	utility::vector1 < signed long > expected_by_type_absolute_;
	
	/// @brief The residue property sets that we'll be counting.
	///
	utility::vector1 < AACompositionPropertiesSetOP > property_sets_;
	
	/// @brief The expected number of residues for each set of properties (by property set index), as a fraction (from 0 to 1) of total.
	///
	utility::vector1 < core::Real > expected_by_properties_fraction_;

	/// @brief The expected number of residues for each set of properties (by property set index), as an absolute number.
	/// @details If this is negative, the fraction is used instead.
	utility::vector1 < signed long > expected_by_properties_absolute_;
	
	/// @brief Penalties for each residue count (relative to expected) and each property set.
	/// @brief Outer vector is the property set, by internal index.  Inner vector is
	/// the penalty value for a given deviation from the ideal count, by internal index.
	utility::vector1 < utility::vector1 < core::Real > > property_penalties_;
	
	/// @brief Deviation ranges for each property set.
	/// @brief Outer vector is the property set, by internal index.  Inner vector is pairs
	/// of deviation ranges (e.g. -10,5 indicating that we are storing penalties for having 10
	/// residues too few up to 5 residues too many for a given property set).  Deviations outside
	/// of this range are assigned the value from the end of the range.
	utility::vector1 < std::pair < signed long, signed long > > property_deviation_ranges_;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
