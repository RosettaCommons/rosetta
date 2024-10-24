// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/simple_moves/BackboneMover.hh
/// @brief  Declarations and simple accessor/mutator definitions for SmallMover and ShearMover


#ifndef INCLUDED_protocols_simple_moves_BackboneMover_hh
#define INCLUDED_protocols_simple_moves_BackboneMover_hh

// Unit headers
#include <protocols/simple_moves/BackboneMover.fwd.hh>

// Package headers
#include <protocols/canonical_sampling/ThermodynamicMover.hh>

// Project headers
#include <core/types.hh>
#include <core/id/DOF_ID_Range.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/movemap/MoveMapFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>


// Utility Header
#include <utility/vector1.hh>

// Basic Header
#include <basic/datacache/DataMap.fwd.hh>

// C++ Headers
#include <map>
#include <string>
#include <tuple>

#include <core/id/TorsionID.hh> // MANUAL IWYU

namespace protocols {
namespace simple_moves {

///////////////////////////////////////////////////////////////////////////////
/// @brief BackboneMover class has elements of the MC temperature to do repetitions
/// of bb moves (small, shear, wobble, etc.).
class BackboneMover : public protocols::canonical_sampling::ThermodynamicMover {

public:
	typedef core::Real Real;


public:
	// empty constructor fills the mover with default values
	// default values from smallmoves.cc of Rosetta++ (small_move_param)
	BackboneMover();

	BackboneMover(
		core::kinematics::MoveMapOP movemap_in,
		core::Real temperature_in,
		core::Size nmoves_in
	);

	/// @brief Copy constructor.
	BackboneMover( BackboneMover const &src );

	// Destructor
	~BackboneMover() override;

	// Virtual functions that get overridden or called from the inheriting classes.
	void apply( core::pose::Pose & ) override;

	std::string get_name() const override;

	void show( std::ostream & output=std::cout ) const override;

	virtual void setup_list( core::pose::Pose & ) = 0;

	virtual void set_angles( core::Real ) = 0;

	bool make_move( core::pose::Pose & );

	virtual bool move_with_scorefxn( core::pose::Pose & pose ) = 0;
	virtual bool move_with_rama( core::pose::Pose & pose ) = 0;

	void test_move( core::pose::Pose & ) override;

	void clear();

	bool check_rama();


	// Properties set/get functions
	void temperature( core::Real const temperature_in ) { temperature_ = temperature_in; }

	core::Real temperature() const { return temperature_; }

	void nmoves( core::Size const nmoves_in ) { nmoves_ = nmoves_in; }

	core::Size nmoves() const { return nmoves_; }

	void movemap(core::kinematics::MoveMapOP new_movemap) { movemap_= new_movemap; }
	void movemap_factory(core::select::movemap::MoveMapFactoryCOP new_movemap_factory) { movemap_factory_= new_movemap_factory; }

	/// @brief Get the movemap for this pose.
	core::kinematics::MoveMapCOP movemap( core::pose::Pose const & pose );

	/// @brief Return the residue selector that can be used to dynamically select residues to perturb
	/// at runtime.
	core::select::residue_selector::ResidueSelectorCOP selector() const { return selector_; }

	void selector( core::select::residue_selector::ResidueSelectorCOP selector ) { selector_ = selector; }


	/// @brief If set, this scorefunction will be used instead of the original rama method (default = not set).
	core::scoring::ScoreFunctionOP scorefxn() const { return scorefxn_; }

	void scorefxn( core::scoring::ScoreFunctionOP sfxn ) { scorefxn_ = sfxn; }


	/// @brief Set the maximum angle of perturbation, independent of
	/// secondary structure. new_angle = old_angle +/- ( angle_max/2 )
	///
	/// Example:
	///     bbmover.angle_max(25)
	/// See also:
	///     ShearMover
	///     SmallMover
	void angle_max( core::Real const angle );

	/// @brief Sets the max angle of perturbation for residues with <type>
	/// secondary structure.  (<type> must be 'H', 'E', or 'L'.) new_angle = old_angle +/- ( angle_max/2 )
	///
	/// Example:
	///     bbmover.angle_max('H', 25)
	///
	/// See also:
	///     ShearMover
	///     SmallMover
	void angle_max( char const type, core::Real const angle );

	// Note: Pass in by value for one-direction assignment.
	/// @brief Sets the max angle of perturbation, for secondary structures
	/// 'H', 'E', and 'L'.
	/// new_angle = old_angle +/- ( angle_max/2 )
	///
	void angle_max( std::map< char, core::Real > angle_max_in );

	/// @brief Gets the max angle of perturbation for residues with <type>
	/// secondary structure.  (<type> must be 'H', 'E', or 'L'.)
	///
	/// Example:
	///     bbmover.angle_max('H')
	///
	/// See also:
	///     ShearMover
	///     SmallMover
	core::Real get_angle_max(char const type) const;

	/// @brief Set the ResidueSelector that this mover will use.
	/// @details Clones the input.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	static utility::tag::XMLSchemaComplexTypeGeneratorOP complex_type_generator_for_backbone_mover(
		utility::tag::XMLSchemaDefinition & xsd );

	static std::string backbone_mover_complex_type_namer( std::string tag_name );

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	/// @brief get whether detailed balance is preserved (i.e. no Ramachandran biasing)
	bool preserve_detailed_balance() const override { return preserve_detailed_balance_; }

	/// @brief set whether detailed balance is preserved (i.e. no Ramachandran biasing)
	void
	set_preserve_detailed_balance( bool preserve_detailed_balance ) override
	{
		preserve_detailed_balance_ = preserve_detailed_balance;
	}

	/// @brief get the DOF_IDs perturbed by the mover during moves, along with their ranges
	utility::vector1<core::id::DOF_ID_Range> dof_id_ranges(core::pose::Pose & pose) override = 0;

	///@brief Set to warn of an empty setup list. Default true.
	void
	set_warn_empty_list(bool warn_empty);

protected:
	/// @brief Set the ResidueSelector that this mover will use.
	core::select::residue_selector::ResidueSubset compute_selected_residues( core::pose::Pose const & pose ) const;

	/// @brief Return a list of TorsionIDs for the standard main-chain torsions of this residue.
	utility::vector1< core::id::TorsionID > get_mainchain_TorsionIDs(
		core::conformation::Conformation const & conf,
		core::uint const seqpos ) const;

	/// @brief Is this set of torsions allowed to move by the MoveMap?
	bool are_torsions_allowed(
		utility::vector1< core::id::TorsionID > torsions,
		core::kinematics::MoveMapCOP mm ) const;

private:
	core::select::movemap::MoveMapFactoryCOP movemap_factory_;
	core::kinematics::MoveMapOP movemap_;
	core::scoring::ScoreFunctionOP scorefxn_;

	// Controls bias w/ which uphill moves are accepted.
	core::Real temperature_;

	// Number of positions at which to make moves.
	core::Size nmoves_;
	bool warn_empty_list_ = true; //Warn if we try an empty list.  Typically, this is OK but for some ways to use these we can't tell ahead of time if we have residues or not (ShearMover, glycans especially).


protected:
	// Max allowed angle-change as a function of ss type.
	std::map< char, core::Real > angle_max_;

	// Variables for the apply.
	core::uint resnum_;
	core::Real big_angle_, small_angle_;
	utility::vector1< core::id::TorsionID > moving_torsions_;
	utility::vector1< std::tuple< core::uint, utility::vector1< core::id::TorsionID >, core::Real > > move_pos_list_;

protected:
	core::Real old_rama_score_, new_rama_score_;

private:
	bool preserve_detailed_balance_;

	// The ResidueSelector to use.
	core::select::residue_selector::ResidueSelectorCOP selector_;

};  // class BackboneMover

std::ostream & operator<< ( std::ostream & os, BackboneMover const & mover );


///////////////////////////////////////////////////////////////////////////////
/// @brief A mover that makes independent random perturbations of the phi and
/// psi torsion angles of residue i. It selects residue i at random among
/// movable residues (set by its MoveMap), and the final torsion angle
/// is subject to a metropolis criterion using the rama score to ensure that
/// only favorable backbone torsion angles are being selected. The number of
/// perturbations, and the magnitude of perturbations, and the temperature
/// in the rama check, can all be modified.
///
/// Common Methods:
///     SmallMover.apply
///     SmallMover.angle_max
class SmallMover : public BackboneMover {

public:
	// default constructor
	SmallMover();

	/// @brief Constructs a SmallMover
	/// smallmover = SmallMover( movemap , kT , n_moves )
	///
	/// MoveMap        movemap   /object storing BB torsion movability
	/// Real (float)   kT        /used in rama Metropolis Criterion
	/// core::Size (int)     n_moves   /the number of perturbations for one move
	SmallMover(
		core::kinematics::MoveMapOP movemap_in,
		core::Real temperature_in,
		core::Size nmoves_in );

	// destructor
	~SmallMover() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void setup_list( core::pose::Pose & pose ) override;

	void set_angles( core::Real angle_in ) override;

	/// @brief get the TorsionIDs perturbed by the mover during moves, along with their ranges
	utility::vector1<core::id::TorsionID_Range> torsion_id_ranges( core::pose::Pose & pose ) override;

	/// @brief get the DOF_IDs perturbed by the mover during moves, along with their ranges
	utility::vector1<core::id::DOF_ID_Range> dof_id_ranges( core::pose::Pose & pose ) override;

	std::string get_name() const override;

	static std::string mover_name();

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	bool move_with_scorefxn( core::pose::Pose & pose ) override;
	bool move_with_rama( core::pose::Pose & pose ) override;

};  // class SmallMover


///////////////////////////////////////////////////////////////////////////////
/// @brief A mover that perturbs the phi of residue i and the psi of residue
/// i-1 such that they create a 'shearing' effect, minimizing the downstream
/// consequences of this torsional perturbation. The final torsion angle
/// is subject to a metropolis criterion using the rama score to ensure that
/// only favorable backbone torsion angles are being selected. The number of
/// perturbations, and the magnitude of perturbations, and the temperature
/// in the rama check, can all be modified.
///
/// Common Methods:
///     ShearMover.apply
class ShearMover : public BackboneMover {

public:
	// default constructor
	ShearMover();

	/// @brief Constructs a ShearMover
	/// shearmover = ShearMover( movemap , kT , n_moves )
	///
	/// MoveMap        movemap   /object storing BB torsion movability
	/// Real (float)   kT        /used in rama Metropolis Criterion
	/// core::Size (int)     n_moves   /the number of perturbations for one move
	ShearMover(
		core::kinematics::MoveMapOP movemap_in,
		core::Real temperature_in,
		core::Size nmoves_in
	);

	//destructor
	~ShearMover() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void setup_list( core::pose::Pose & pose ) override;

	void set_angles( core::Real angle_in ) override;

	/// @brief get the TorsionIDs perturbed by the mover during moves, along with their ranges
	utility::vector1< core::id::TorsionID_Range > torsion_id_ranges( core::pose::Pose & pose ) override;

	/// @brief get the DOF_IDs perturbed by the mover during moves, along with their ranges
	utility::vector1< core::id::DOF_ID_Range > dof_id_ranges( core::pose::Pose & pose ) override;

	std::string get_name() const override;

	static std::string mover_name();

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	bool move_with_scorefxn( core::pose::Pose & pose ) override;
	bool move_with_rama( core::pose::Pose & pose ) override;

	std::pair< utility::vector1< core::id::TorsionID >, utility::vector1< core::id::TorsionID > >
	setup_list_for_saccharide_residue( core::pose::Pose & pose, core::uint seq_pos );
};  // class ShearMover

} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_BackboneMover_HH
