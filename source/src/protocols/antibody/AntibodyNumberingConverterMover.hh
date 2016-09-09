// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/AntibodyNumberingConverterMover.hh
/// @brief Converts numbering schemes of an antibody.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_AntibodyNumberingConverterMover_hh
#define INCLUDED_protocols_antibody_AntibodyNumberingConverterMover_hh

// Unit headers
#include <protocols/antibody/AntibodyNumberingConverterMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/antibody/AntibodyEnum.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace antibody {

///@brief Converts numbering schemes of an antibody, independant of AntibodyInfo.
///        By default, works on a SINGLE antibody FAB with chains L and H, as in the rest of Rosetta.
///
///@details Usable Numbering Schemes:
///
///            Chothia_Scheme
///            Enhanced_Chothia_Scheme
///            AHO_Scheme
///            Kabat_Scheme
///            IMGT_Scheme
///
///
class AntibodyNumberingConverterMover : public protocols::moves::Mover {

public:

	AntibodyNumberingConverterMover();
	
	AntibodyNumberingConverterMover(AntibodyNumberingSchemeEnum const from, AntibodyNumberingSchemeEnum const to);
	
	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~AntibodyNumberingConverterMover();


	
public:

	///@brief Set the scheme to convert the pose into.
	void
	set_scheme_conversion(
	
		AntibodyNumberingSchemeEnum const from,
		AntibodyNumberingSchemeEnum const to
		
	);
	
	virtual void
	apply( core::pose::Pose & pose );


public:
	virtual void
	show( std::ostream & output = std::cout ) const;

	virtual std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//AntibodyNumberingConverterMover & operator=( AntibodyNumberingConverterMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

	static std::string
	class_name();
	
private:

	void
	set_defaults();
	
private:
	
	AntibodyNumberingSchemeEnum from_scheme_;
	AntibodyNumberingSchemeEnum to_scheme_;
	
};

std::ostream &
operator<<( std::ostream & os, AntibodyNumberingConverterMover const & mover );

} //protocols
} //antibody

#endif //protocols/antibody_AntibodyNumberingConverterMover_hh
