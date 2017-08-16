// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/potentials/SetSecStructEnergies.hh
/// @brief mover for setting centroid score of secondary structure through parser
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_fldsgn_potentials_SetSecStructEnergies_hh
#define INCLUDED_protocols_fldsgn_potentials_SetSecStructEnergies_hh

// unit headers
#include <protocols/fldsgn/potentials/SetSecStructEnergies.fwd.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasHelicesSheetPotential.fwd.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasHelixPairPotential.fwd.hh>

// type headers
#include <core/types.hh>

// project headers
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/parser/BluePrint.fwd.hh>
#include <protocols/moves/Mover.hh>

// C++ headers
#include <string>

#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {


class SetSecStructEnergies : public protocols::moves::Mover {


private: // typedefs


	typedef protocols::moves::Mover Super;


public: // typedefs


	typedef std::string String;

	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;

	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::parser::BluePrint BluePrint;
	typedef protocols::parser::BluePrintOP BluePrintOP;

	typedef core::conformation::symmetry::SymmetryInfoOP SymmetryInfoOP;
	typedef core::conformation::symmetry::SymmetryInfoCOP SymmetryInfoCOP;

	typedef protocols::fldsgn::potentials::sspot::NatbiasHelixPairPotential NatbiasHelixPairPotential;
	typedef protocols::fldsgn::potentials::sspot::NatbiasHelixPairPotentialOP NatbiasHelixPairPotentialOP;
	typedef protocols::fldsgn::potentials::sspot::NatbiasHelicesSheetPotential NatbiasHelicesSheetPotential;
	typedef protocols::fldsgn::potentials::sspot::NatbiasHelicesSheetPotentialOP NatbiasHelicesSheetPotentialOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public: // construct/destruct


	/// @brief default constructor
	SetSecStructEnergies();

	/// @brief value constructor
	SetSecStructEnergies( ScoreFunctionOP const sfx, String const & filename, bool const ss_from_blueprint=true );

	/// @brief value constructor
	SetSecStructEnergies( ScoreFunctionOP const sfx, BluePrintOP const blueprintOP, bool const ss_from_blueprint=true );

	/// @brief default destructor
	virtual ~SetSecStructEnergies();


private: // disallow assignment


	/// @brief copy assignment
	/// @remarks Mover base class prevents this from working properly...
	SetSecStructEnergies & operator =( SetSecStructEnergies const & rval );


public: // virtual constructors


	/// @brief clone this object
	MoverOP clone() const override;

	/// @brief create this type of object
	MoverOP fresh_instance() const override;


public: // mutators

	/// @brief access ptr to the modified score function
	ScoreFunctionOP
	scorefunction_ptr() const;

	/// @brief set the centroid level score function to be modified
	/// @details this also clones the object of ptr isn't null and stores the clone in sfx_orig_
	void set_scorefunction_ptr( ScoreFunctionOP sfx );

	/// @brief sets the secondary structure to be used for computation
	/// @param[in] secstruct Secondary structure to be forced on the pose
	/// @details If this is non-empty, it will be used as the pose secondary structure to
	///          determine secondary structure elements in the pose
	void
	set_secstruct( String const & secstruct );

public: // virtual main methods


	/// @brief apply defined moves to given Pose
	void apply( Pose & pose ) override;



public: //parser


	/// @brief parse xml file
	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const &,
		Movers_map const &,
		Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private: // helper functions


	String symmetric_secstruct( SymmetryInfoCOP const syminfo, String const & ss ) const;

	/// @brief chooses and return the secondary structure to be used in computation
	/// @param[in] pose Input pose
	/// @details  If secstruct_ is set, returns that.
	///           If use_dssp_ is set, returns DSSP secstruct string
	///           Otherwise, returns pose secondary structure
	std::string
	get_secstruct( Pose const & pose ) const;

	std::string
	get_helix_pairings( Pose const & pose ) const;

	/// @brief initializes pairings and secondary structure from blueprint file
	/// @param[in] bp                BluePrint object
	/// @param[in] ss_from_blueprint If true, secstruct_ will be set from the blueprint.  If false,
	///                              secstruct_ will be unchanged, and only the pairings will be set.
	void
	init_from_blueprint( BluePrint const & bp, bool const ss_from_blueprint );

private: // data

	/// @brief Secondary structure to be used
	std::string secstruct_;

	/// @brief If true, and secstruct is not set, DSSP will be used to determine secondary structure
	bool use_dssp_;

	/// @brief HHPAIR definition to be used
	std::string hh_pair_;

	/// @brief strand pairing definition to be used
	std::string ss_pair_;

	/// @brief helix-strand-strand triplet definition to be used
	std::string hss_triplet_;

	/// @brief Original copy of the centroid scorefunction to use
	ScoreFunctionCOP sfx_orig_;

	/// @brief the centroid scorefunction to use, default "remodel_cen"
	ScoreFunctionOP sfx_;

	/// @brief weight for helix and helix pairing potential
	Real hh_weight_;

	/// @brief weight for heilx and strand pairing potential
	Real hs_weight_;

	/// @brief weight for strand and strand pairing potential
	Real ss_weight_;

	/// @brief weight for sheet twist potential
	Real stwist_weight_;

	/// @brief weight for hs_pair
	Real hs_pair_weight_;

	/// @brief weight for ss_pair
	Real ss_pair_weight_;

	/// @brief weight for strand and strand pairing potential
	Real rsigma_weight_;

	/// @brief if this option is set, a copy of the input pose will have symmetry added to it, and this
	///        symmetry will be used to modify the secondary structure. The input pose will be
	///        unchanged. I can't think of any reason why anyone would want this to be true if the input pose
	///        is already symmetric.
	/// @details TL: this is added as a hack to work around behavior by SetupForSymmetryMover. SetupForSymmetry
	///          sets this command line flag to a bogus value in its parse_my_tag() function. Therefore,
	///          the presence of the option at runtime does not indicate the user's intentions. As a workaround,
	///          the option is queried set at parse time and this variable is used after than.
	///          WARNING: You should declare your SetupForSymmetry movers AFTER this one, or you will probably
	///          get crashes.
	///          Default is the value of option[ symmetry::symmetry_definition ].active() at parse-time
	bool add_symmetry_;

	NatbiasHelixPairPotentialOP hpairpot_;

	NatbiasHelicesSheetPotentialOP hspot_;

};

} // namespace potentials
} // namespace fldsgn
} // namespace protocols


#endif /* INCLUDED_protocols_forge_components_SetSecStructEnergies_HH */
