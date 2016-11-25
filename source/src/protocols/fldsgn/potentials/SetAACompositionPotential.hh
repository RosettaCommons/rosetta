// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/components/SetAACompositionPotential.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_potentials_SetAACompositionPotential_HH
#define INCLUDED_protocols_fldsgn_potentials_SetAACompositionPotential_HH

// unit headers
#include <protocols/fldsgn/potentials/SetAACompositionPotential.fwd.hh>

// project headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>

// type headers
#include <core/types.hh>

// C++ headers
#include <map>
#include <string>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {

class SetAACompositionPotential : public protocols::moves::Mover {


private: // typedefs


	typedef protocols::moves::Mover Super;


public: // typedefs


	typedef std::string String;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef core::chemical::AA AA;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef protocols::moves::MoverOP MoverOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public: // construct/destruct


	/// @brief default constructor
	SetAACompositionPotential();

	/// @brief copy constructor
	SetAACompositionPotential( SetAACompositionPotential const & rval );

	/// @brief default destructor
	virtual ~SetAACompositionPotential();


private: // disallow assignment


	/// @brief copy assignment
	/// @remarks Mover base class prevents this from working properly...
	SetAACompositionPotential & operator =( SetAACompositionPotential const & rval );


public: // virtual constructors


	/// @brief clone this object
	MoverOP clone() const override;

	/// @brief create this type of object
	MoverOP fresh_instance() const override;


public: // virtual main methods


	/// @brief apply defined moves to given Pose
	void apply( Pose & ) override;

	// XRW TEMP  virtual
	// XRW TEMP  std::string get_name() const;


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


	bool set_parameters( String const & file );


private: // data


	/// @brief
	std::map< core::chemical::AA, std::pair< Real, Real > > comp_constraint_aas_;

	/// @brief weight for composition constraints
	Real weight_;

	/// @brief scorefunction to use
	ScoreFunctionOP sfx_;

	/// @brief
	bool loaded_;


};

} // namespace potentials
} // namespace fldsgn
} // namespace protocols


#endif /* INCLUDED_protocols_fldsgn_SetAACompositionPotential_HH */
