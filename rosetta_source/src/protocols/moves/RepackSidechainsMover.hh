// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange


#ifndef INCLUDED_protocols_moves_RepackSidechainsMover_hh
#define INCLUDED_protocols_moves_RepackSidechainsMover_hh

// Unit headers
#include <protocols/moves/RepackSidechainsMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>



namespace protocols {
namespace moves {

/// @brief A Mover that packs the side-chains (very similar to pack_missing_sidechains()

class RepackSidechainsMover : public Mover {
public:
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;

public:
	///@brief default constructor
	RepackSidechainsMover();

	/// @brief constructor with typename
	//	RepackSidechainsMover( std::string const & );

	/// @brief Constructs a RepackSidechainsMover with PackerTask  <task>
	/// evaluated using  <scorefxn>
	///
	/// ScoreFunction  scorefxn   /function to minimize while changine rotamers
	/// PackerTask     task       /object specifying what to design/pack
	/// Size (int)     nloop      /number of loops in the Pose (???)
	RepackSidechainsMover(
		ScoreFunctionCOP scorefxn
	);

	// copy constructor
	RepackSidechainsMover( RepackSidechainsMover const & other );

	virtual void apply( Pose & pose );
	virtual std::string get_name() const;

	///@brief parse XML (specifically in the context of the parser/scripting scheme)
	virtual void parse_my_tag(
		TagPtr const,
		DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & );

	///@brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_score_function(
		TagPtr const,
		DataMap const &,
		Filters_map const &,
		Movers_map const &,
		Pose const & );

	///@brief required in the context of the parser/scripting scheme
	virtual MoverOP fresh_instance() const;

	///@brief required in the context of the parser/scripting scheme
	virtual MoverOP clone() const;

	void set_scorefxn( ScoreFunctionCOP sf );

	ScoreFunctionCOP scorefxn() const { return scorefxn_; };

protected:
private:
	// pointers to data that are passed in
	ScoreFunctionCOP scorefxn_;
};

// note: it is better to create new files, instead of adding additional classes here

} // moves
} // protocols

#endif
