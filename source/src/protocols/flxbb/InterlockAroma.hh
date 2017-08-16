// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/flxbb/InterlockAroma.hh
/// @brief perform cycles of design and relax with filter
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_flxbb_InterlockAroma_hh
#define INCLUDED_protocols_flxbb_InterlockAroma_hh

// Unitt Header
#include <protocols/flxbb/InterlockAroma.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace flxbb {


///////////////////////////////////////////////////////////////////////////////////////////////////////
class InterlockAroma: public protocols::moves::Mover {
public:


	typedef protocols::moves::Mover Super;

	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef protocols::moves::MoverOP MoverOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public: // constructor/destructor


	/// @brief default constructor
	InterlockAroma();

	/// @brief copy constructor
	InterlockAroma( InterlockAroma const & rval );

	/// @brief destructor
	~InterlockAroma() override;


public: // virtual constructors


	/// @brief clone this object

	MoverOP clone() const override;


	/// @brief create this type of object

	MoverOP fresh_instance() const override;


public: // virtual main operation


	/// @brief mover apply
	void apply( Pose & pose ) override;



public:// parser


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



private:


	/// @brief score function
	ScoreFunctionOP scorefxn_;

	/// @brief input secondary structure information
	String input_ss_;

	/// @brief
	Real max_repulsion_energy_;

	/// @brief
	Real min_env_energy_;

	// @brief Exclude aromatic chi2 rotamers, of which angles are around 0
	bool limit_aroma_chi2_;

	// @brief
	bool output_pdbs_;

	// @brief
	bool verbose_;


};


} // namespace flxbb
} // namespace protocols

#endif
