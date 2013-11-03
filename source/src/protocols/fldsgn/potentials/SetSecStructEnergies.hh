// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/jd2/parser/BluePrint.fwd.hh>
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

	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::jd2::parser::BluePrint BluePrint;
	typedef protocols::jd2::parser::BluePrintOP BluePrintOP;

  typedef core::conformation::symmetry::SymmetryInfoOP SymmetryInfoOP;

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

	/// @brief copy constructor
	SetSecStructEnergies( SetSecStructEnergies const & rval );

	/// @brief default destructor
	virtual	~SetSecStructEnergies();


private: // disallow assignment


	/// @brief copy assignment
	/// @remarks Mover base class prevents this from working properly...
	SetSecStructEnergies & operator =( SetSecStructEnergies const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	MoverOP clone() const;

	/// @brief create this type of object
	virtual
	MoverOP fresh_instance() const;


public: // mutators


	/// @brief define secondary structrue by blueprint
	inline void ss_from_blueprint( bool const flag ){ ss_from_blueprint_ = flag;	}

	/// @brief set the centroid level score function
	void scorefunction( ScoreFunction const & sfx );

	/// @brief set the centroid level score function
	void scorefunction( ScoreFunctionOP sfx );

	/// @brief set blueprint file by filename
	void set_blueprint( String const & filename );

	/// @brief set blueprint file
	void set_blueprint( BluePrintOP const blp );



public: // virtual main methods


	/// @brief apply defined moves to given Pose
	virtual
	void apply( Pose & pose );

	virtual
	std::string get_name() const;


public: //parser


	/// @brief parse xml file
	void parse_my_tag( TagCOP tag,
										 basic::datacache::DataMap & data,
										 Filters_map const &,
										 Movers_map const &,
										 Pose const & );


private: // helper functions


	String symmetric_secstruct( SymmetryInfoOP const syminfo, String const & ss );


private: // data


	/// @brief
	bool loaded_;

	/// @brief bluerprint file for setting build instruction
	BluePrintOP blueprint_;

	/// @brief
	bool ss_from_blueprint_;

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

	NatbiasHelixPairPotentialOP hpairpot_;

	NatbiasHelicesSheetPotentialOP hspot_;

};

} // namespace potentials
} // namespace fldsgn
} // namespace protocols


#endif /* INCLUDED_protocols_forge_components_SetSecStructEnergies_HH */
