// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/relax/membrane/MPRangeRelaxMover.hh
/// @brief      Relaxes a membrane protein by relaxing in ranges
/// @details Relaxes a membrane protein by iteratively relaxing ranges of the protein;
///    No ramping required. Much faster than FastRelax and good for
///    large to very large proteins (tested up to 5250 residues)
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_relax_membrane_MPRangeRelaxMover_hh
#define INCLUDED_protocols_relax_membrane_MPRangeRelaxMover_hh

// Unit Headers
#include <protocols/relax/membrane/MPRangeRelaxMover.fwd.hh>
//#include <protocols/membrane/MPRangeRelaxMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace relax {
namespace membrane {

class MPRangeRelaxMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	MPRangeRelaxMover();

	/// @brief Copy Constructor
	MPRangeRelaxMover( MPRangeRelaxMover const & src );

	/// @brief Assignment Operator
	MPRangeRelaxMover & operator = ( MPRangeRelaxMover const & src );

	/// @brief Destructor
	virtual ~MPRangeRelaxMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Parse Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	/// @brief Provide xml schema for RosettaScripts compatibility
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (MPRangeRelaxMover)
	virtual std::string get_name() const override;

	/// @brief Do a RangeRelax of a membrane protein
	virtual void apply( Pose & pose ) override;

	/// @brief Optimize membrane
	void optimize_membrane( bool yesno );

	/////////////////////
	///  Get Methods  ///
	/////////////////////

	/// @brief Get native
	core::pose::PoseOP get_native() const;

	/// @brief Get scorefxn
	core::scoring::ScoreFunctionOP get_sfxn() const;

	/// @brief Get center residue number
	core::Size get_center_resnumber() const;

	/// @brief Get yes/no for set helical secondary structure in TM region
	bool get_set_tm_helical() const;

	/// @brief Get yes/no for optimize membrane option
	bool get_optmem() const;

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	void register_options();

	/// @brief Set default values
	void set_defaults();

	/// @brief Init from commandline
	void init_from_cmd();

private: // data

	/// @brief Native
	core::pose::PoseOP native_;

	/// @brief Scorefxn
	core::scoring::ScoreFunctionOP sfxn_;

	/// @brief Center residue number
	core::Size center_resnumber_;

	/// @brief Set helical secondary structure in TM region
	bool set_tm_helical_;

	/// @brief optimize membrane?
	bool optmem_;

	// Tell the compiler that we are not hiding the base
	// function with the parse_my_tag written above
	using protocols::moves::Mover::parse_my_tag;

};

} // membrane
} // relax
} // protocols

#endif // INCLUDED_protocols_relax_membrane_MPRangeRelaxMover_hh
