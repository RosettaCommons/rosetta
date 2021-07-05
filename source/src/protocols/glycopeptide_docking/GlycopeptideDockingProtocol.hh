// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/glycopeptide_docking/GlycopeptideDockingProtocol.hh
/// @brief A protocol to dock peptide or glycosylate peptide in the glycopeptide_docking site of
/// a glycosyltransferase. The protocol is designed for O-linked glycopeptide_docking in which
/// sugar are transferred one-at-a-time by specific glycosytransferases
/// to a peptide (or glycopeptide) or protein substrate from a nucleatide donor (eg. UDP-GalNAc)
/// The current version of the protocol has been tested on GalNAcT transferases that add the
/// first sugar in the O-linked pathway. The protocol is in active development for extension to
/// other enzymes in the O-glycopeptide_docking pathway.
/// Currently, the protocol requires the pose to contain the enzyme, activated sugar (or just UDP),
/// and a peptide substrate.
/// @author Yashes Srinivasan (yashess@gmail.com)

#ifndef INCLUDED_protocols_glycopeptide_docking_GlycopeptideDockingProtocol_HH
#define INCLUDED_protocols_glycopeptide_docking_GlycopeptideDockingProtocol_HH

// Unit headers
#include <protocols/glycopeptide_docking/GlycopeptideDockingProtocol.fwd.hh>
#include <protocols/glycopeptide_docking/GlycopeptideDockingFlags.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/constraint_movers/ConstraintSetMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>

// Utility headers

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

#include <string>

namespace protocols {
namespace glycopeptide_docking {

///@brief A protocol to predict glycosylation and elongation of glycans.
/// @details This protocol is intended to "sample" glycopeptide
/// poses/conformations in the context of the glycosyltransferase "ready"
/// to be glycosylated. Similar to enzyme-design protocols, constraints can be
/// specified to emulate the positioning of atoms/residues in the transition state
/// comingsoon(3-4 months): We are currently prototyping an extension to generalize the
/// protocol for elongation to core1, core2 and higher sugars with their respective
/// glcyosyltransferases. We hope to automate substrate generation from peptide-sequence
/// and glycosylation (for glycopeptide substrates), donor glycosylation (e.g. UDP->UDP-Sugar)
/// and initial setup of the peptide-substrate in the active site.
/// Additional options are defined in GlycopeptideDockingFlags.hh

class GlycopeptideDockingProtocol : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	GlycopeptideDockingProtocol();

	/// @brief Constructor with flags
	GlycopeptideDockingProtocol( GlycopeptideDockingFlags &flags);

	/// @brief Copy constructor (not needed unless you need deep copies)
	GlycopeptideDockingProtocol ( GlycopeptideDockingProtocol const & object_to_copy );

	/// @brief Assignment operator
	GlycopeptideDockingProtocol & operator=( GlycopeptideDockingProtocol const & object_to_copy );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~GlycopeptideDockingProtocol() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	/// @brief Show mover information
	void
	show( std::ostream & output = std::cout ) const override;

	/// @brief Set reference pose (native) from file name
	void
	set_ref_pose_from_filename( std::string const & filename );

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
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

private: // methods

	/// @brief Copy protocol data
	void
	copy_data( GlycopeptideDockingProtocol & object_to_copy_to, GlycopeptideDockingProtocol const & object_to_copy_from );

	/// @brief Prepare scorefunction. Currently the scorefunction
	/// is fixed by design. Benchmarking has not been done with other
	/// scorefunctions.
	void
	prepare_score_function();

private: // data

	protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags_;
	//////////////////////////
	//////// Scoring /////////
	//////////////////////////

	core::scoring::ScoreFunctionOP sf_;
	protocols::moves::MonteCarloOP mc_;


	//////////////////////////
	//////// Movers //////////
	//////////////////////////

	protocols::constraint_movers::ConstraintSetMoverOP constraint_setter_;


	//////////////////////////
	//////// Constants ///////
	//////////////////////////

	core::pose::PoseOP ref_pose_=nullptr;
	core::kinematics::FoldTreeOP ft_docking_=nullptr;
	core::kinematics::FoldTreeOP ft_substrate_=nullptr;

};



std::ostream &
operator<<( std::ostream & os, GlycopeptideDockingProtocol const & mover );

} //glycopeptide_docking
} //protocols

#endif //protocols_glycopeptide_docking_GlycopeptideDockingProtocol_HH
