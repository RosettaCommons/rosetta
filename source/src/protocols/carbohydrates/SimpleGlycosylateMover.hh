// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/SimpleGlycosylateMover.hh
/// @brief A mover for glycosylation of common biological glycosylations.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_carbohydrates_SimpleGlycosylateMover_hh
#define INCLUDED_protocols_carbohydrates_SimpleGlycosylateMover_hh

// Unit headers
#include <protocols/carbohydrates/SimpleGlycosylateMover.fwd.hh>
#include <protocols/moves/Mover.hh>


#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace carbohydrates {



/// @brief A mover for glycosylation of biological glycosylations.
///  Currently glysolylation is done based on string, not from PDB.
///  Use GlycanRelax to model the resulting glycosylation!.
///
/// @details
///  Single Glycosylation:
///    If a single glycosylation is passed, it will glycosylate all positions set with this glycan.
///
///  Multiple Glycosylations:
///    If multiple glycosylations are passed, will randomly select from these on apply to sample glycan hetergenecity.
///    see glyco.set_weights to set weights for these glycans to sample non-random heterogenecity in the glycoform.
///
///  The site should be ASN for N-linked glycosylations OR SER/THR for O-linked glycosylations
///    If a glycan already exists, will strip off the current glycan by default.
///    set glyco.set_strip_existing_glycans( false ) to branch off existing glycans instead of deleting them.
///
///    see CreateGlycoSiteMover to create glyco sites (as N-linked glycosylations will need a specific motif to be biological)
///
///  Will randomly select from set positions and glycosylate all positions set.
///
///  Glycosylations:
///    1) If your name ends with .iupac or .gws (GlycoWorkBench), will try to load the file
///
///    2) Next, it will check the short names in the Rosetta database for your string.
///       If the string is in common_names.txt, will load the paired iupac sequence.
///       See database/chemical/carbohydrates/common_glycans/common_names.txt for accepted short names.
///     Names include man3, man5, and man9.
///
///    3)
///      If the name is not found, will attempt to build the glycan as an iupac sequence from the string.
///
///
///  TODO JAB: Add Glycosylate from pdb/cif files in addition to IUPAC and short names.
///   No current RosettaCarbohydrate functionality to load from PDB/CIF and attach to an existing pose..
///
class SimpleGlycosylateMover : public protocols::moves::Mover {

public:

	SimpleGlycosylateMover();

	// copy constructor
	SimpleGlycosylateMover( SimpleGlycosylateMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	~SimpleGlycosylateMover() override;

	void
	apply( core::pose::Pose & pose ) override;

public:

	/// @brief Set the glycosylation that will happen.
	///  See database/chemical/carbohydrates/common_glycans for common names.
	///  Names include man3, man5, and man9.
	///
	///  Can also use full iupac names.
	///
	///
	///  See also: glycosylate.set_glycosylation_weights
	///
	void
	set_glycosylation(std::string const & iupac_or_common_string);


	/// @brief Set possible glycosylations - the mover will randomly pick these on apply
	///  See database/chemical/carbohydrates/common_glycans for accepted names.
	///  Names include man3, man5, and man9.
	///
	///  Can also use full iupac names.
	///
	///
	///  See also: glycosylate.set_glycosylation_weights
	///
	void
	set_glycosylations(utility::vector1<std::string> const & iupac_or_common_strings);


	///@brief Set a single resnum position
	void
	set_position(core::Size position);

	///@brief Set multiple positions to glycosylate
	void
	set_positions(utility::vector1< bool > const & positions );

	///@brief Set mutliple positions to glycosylate
	void
	set_positions(utility::vector1< core::Size > const & positions);

	///@brief Set positions to glycosylate using a movemap.
	void
	set_positions_from_movemap( core::kinematics::MoveMapCOP mm);


public:


	///@brief Set weights for potential glycosylation if more than one is set.
	void
	set_glycosylation_weights( utility::vector1< core::Real > const & weights);

	///@brief This sets whether if we already have a glycan at a position, whether to extend it or delete the existing glycan.
	/// Advanced functionality - use with caution.
	///  Not yet work.
	void
	set_strip_existing_glycans( bool strip_existing );


public:


	void
	show( std::ostream & output=std::cout ) const override;

	// XRW TEMP  std::string
	// XRW TEMP  get_name() const override;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;



	//SimpleGlycosylateMover & operator=( SimpleGlycosylateMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	void
	remove_index(utility::vector1< core::Size > & current_vector, core::Size resnum) const;

	utility::vector1< std::string >
	setup_and_load_iupac_sequences() const;


private:

	utility::vector1< std::string > glycosylations_;
	utility::vector1< core::Real > glycosylation_weights_;

	utility::vector1< std::string > parsed_positions_;
	utility::vector1< core::Size > positions_;

	bool strip_existing_glycans_;

	std::string ref_pose_name_; //Used for RosettaScripts

	bool idealize_glycosylation_;

};

std::ostream &operator<< (std::ostream &os, SimpleGlycosylateMover const &mover);


} //protocols
} //carbohydrates


#endif //protocols/carbohydrates_SimpleGlycosylateMover_hh







