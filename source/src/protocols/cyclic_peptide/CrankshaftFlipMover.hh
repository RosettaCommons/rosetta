// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/CrankshaftFlipMover.hh
/// @brief This mover performs a cyclic-aware crankshaft flip at a residue position.
/// @author P. Douglas Renfrew (doug.renfrew@gmail.com)
/// @author James Eastwood (jre318@nyu.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_CrankshaftFlipMover_HH
#define INCLUDED_protocols_cyclic_peptide_CrankshaftFlipMover_HH

// Unit headers
#include <protocols/cyclic_peptide/CrankshaftFlipMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

namespace protocols {
namespace cyclic_peptide {

///@brief This mover performs a cyclic-aware crankshaft flip at a residue position.

///@details In lattice-based polymer simulations, crankshaft flips are moves that involve two adjacent monomers in a loop moving together (http://aip.scitation.org/doi/abs/10.1063/1.431297). In peptides, the equivalent is two backbone atoms, C_i and N_(i+1) moving by correlated rotations in psi_i, phi_(i+1). This move has been observed in several experimental contexts, including in crystal structures of antibody loops(http://doi.org/10.1016/j.jmb.2010.10.030) and in simulations of cyclic peptides (https://pubs.acs.org/doi/10.1021/acs.jctc.6b00193) and in simulations of cyclic peptoids (https://pubs.acs.org/doi/full/10.1021/acs.jpcb.2c01669). It should be useful in loop modeling and cyclic-peptides/peptoids/hybrids. Depending on the residue type different rules are uswed to flip the backbone, outlined in the table below.

///
/// This mover is designed to be general and compatible with both homo- and hetero-polymers of the following type: peptides (not Gly or Pro), peptoids, glycine, proline, and N-methylated amino acids. The purtubations that are allowed for different combinations at the i and i+1 positions and compatible conformations are described in the table below.
///
/// | i            | i+1                    | Action                                                                               | Nickname  |
/// |--------------+------------------------+--------------------------------------------------------------------------------------+-----------|
/// | Peptoid, Gly | Peptoid, Pro, N-Methyl | phi_i negate, psi_i negate, omg_i flip by 180                                        | "Peptoid" |
/// | Any(?)       | Gly                    | psi_i 0<->(-30,180), phi_i+1 180<->(-90,+-60)                                        | "Pre-Gly" |
/// | Peptide      | Peptide (not P)        | if psi_i/phi_i+1 near -120/-60 or 120/60: psi_i negate, omg_i negate, phi_i+1 negate | "Peptide" |
///
///
class CrankshaftFlipMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	CrankshaftFlipMover();

	/// @brief Constructor
	CrankshaftFlipMover(
		core::select::residue_selector::ResidueSelectorCOP selector_in,
		bool allow_peptoid_flip_in,
		bool allow_peptide_flip_in,
		bool allow_pre_gly_flip_in,
		core::Real tolerance_in
	);

	/// @brief Copy Constructor
	CrankshaftFlipMover( CrankshaftFlipMover const &src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~CrankshaftFlipMover() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	/// @brief Check if peptoid flip can be done at position, if it can flip poistion and return true else return false
	bool
	try_peptoid_flip( core::pose::Pose & pose, core::Size res1_index, core::Size res2_index );

	/// @brief Check if peptide flip can be done at position, if it can flip poistion and return true else return false
	bool
	try_peptide_flip( core::pose::Pose & pose, core::Size res1_index, core::Size res2_index );

	/// @brief Check if pre-gly flip can be done at position, if it can flip poistion and return true else return false
	bool
	try_pre_gly_flip( core::pose::Pose & pose, core::Size res1_index, core::Size res2_index );

	/// @brief Setters and Getters
	void
	allow_peptoid_flips( bool af ){ allow_peptoid_flips_ = af; }

	void
	allow_peptide_flips( bool af ){ allow_peptide_flips_ = af; }

	void
	allow_pre_gly_flips( bool af ){ allow_pre_gly_flips_ = af; }

	bool
	allow_peptoid_flips(){ return allow_peptoid_flips_; }

	bool
	allow_peptide_flips(){ return allow_peptide_flips_; }

	bool
	allow_pre_gly_flips(){ return allow_pre_gly_flips_; }

	bool
	angle_degrees_in_range( core::Real angle, core::Real target );

	void
	selector( core::select::residue_selector::ResidueSelectorCOP selector );

	core::select::residue_selector::ResidueSelectorCOP
	selector(){ return selector_; }

	core::Real
	tolerance() const { return tolerance_; }

	void
	tolerance( core::Real tolerance ){ tolerance_ = tolerance; }


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	CrankshaftFlipMover & operator=( CrankshaftFlipMover const & src );

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

public: //Function overrides needed for the citation manager:

	/// @brief Sets up the citation info
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;

private: // methods

private: // data

	// what type of moves to allow
	bool allow_peptoid_flips_;
	bool allow_peptide_flips_;
	bool allow_pre_gly_flips_;

	// allowed tolerance in dihedral angle value for allowing flips
	core::Real tolerance_;

	//
	core::select::residue_selector::ResidueSelectorCOP selector_;

};

std::ostream &
operator<<( std::ostream & os, CrankshaftFlipMover const & mover );

} //cyclic_peptide
} //protocols

#endif //protocols_cyclic_peptide_CrankshaftFlipMover_HH
