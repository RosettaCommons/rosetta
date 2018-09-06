// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/residue_selectors/NativeSelector.hh
/// @brief  A ResidueSelector that applies a given residue selector to the native pose.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_residue_selectors_NativeSelector_HH
#define INCLUDED_protocols_residue_selectors_NativeSelector_HH

// Unit headers
#include <protocols/residue_selectors/NativeSelector.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Basic Headers
#include <basic/datacache/DataMap.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace residue_selectors {

class NativeSelector : public core::select::residue_selector::ResidueSelector {

public:

	/// @brief Constructor.
	NativeSelector();
	NativeSelector( NativeSelector const & src);

	/// @brief Constructor that calls set_residue_selector()
	NativeSelector( core::select::residue_selector::ResidueSelectorCOP inner_selector );

	/// @brief Destructor.
	///
	~NativeSelector() override;

	/// @brief Clone function.
	/// @details Copy this object and return owning pointer to the copy.
	core::select::residue_selector::ResidueSelectorOP clone() const override;

	/// @brief Calls apply on the inner residue selector but passes the native pose instead of pose
	core::select::residue_selector::ResidueSubset apply( core::pose::Pose const & pose ) const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override ;

	/// @brief Get the mover class name.
	std::string get_name() const override;

	/// @brief Get the mover class name.
	static std::string class_name();

	/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


public:
	/// @brief set the pose that will be passed to the inner residue selector
	void set_native_pose( core::pose::PoseCOP native ){
		native_ = std::move( native );
	}

	/// @brief set the inner residue selector. This selector's apply function will be called with the native pose
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector ){
		selector_ = std::move( selector );
	}

private: //data

	/// @brief This selector's apply function will be called with the native pose
	core::select::residue_selector::ResidueSelectorCOP selector_;

	/// @brief the pose that will be passed to the inner residue selector
	/// @details this element is marked as mutable so that we can load it at apply-time if a value is not given
	mutable core::pose::PoseCOP native_;


#ifdef    SERIALIZATION
protected:
	friend class cereal::access;

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //protocols
} //residue_selectors

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_residue_selectors_NativeSelector )
#endif // SERIALIZATION

#endif //INCLUDED
