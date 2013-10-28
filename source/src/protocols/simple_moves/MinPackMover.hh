// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/simple_moves/MinPackMover.hh
/// @brief  protocols::moves::Mover class to invoke core::pack::min_pack
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_MinPackMover_hh
#define INCLUDED_protocols_simple_moves_MinPackMover_hh

// Unit headers
#include <protocols/simple_moves/MinPackMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/pack/task/PackerTask.fwd.hh>

#ifdef __clang__
#include <core/pack/task/PackerTask.hh>
#endif

#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>



namespace protocols {
namespace simple_moves {

/// @brief a mover that packs and minimizes the side-chains.  It uses
/// a ScoreFunction for packing and either a PackerTask, or a TaskFactory that generates
/// a PackerTask for instructions on what rotamer sets are allowed at each residue
/// position during packing.
class MinPackMover : public protocols::moves::Mover {
public:
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::task::TaskFactoryCOP TaskFactoryCOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;

public:
	///@brief default constructor
	MinPackMover();

	///@brief constructor with typename
	MinPackMover( std::string const & );

	MinPackMover(
		ScoreFunctionCOP scorefxn
	);
	MinPackMover(
		ScoreFunctionCOP scorefxn,
		PackerTaskCOP task
	);

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~MinPackMover();

	// copy constructor
	MinPackMover( MinPackMover const & other );

	void
	init();

	// methods
	virtual void apply( Pose & pose );
	virtual std::string get_name() const;
	bool task_is_valid( Pose const & pose ) const; // should this be virtual?

	///@brief parse XML (specifically in the context of the parser/scripting scheme)
	virtual void parse_my_tag(
		TagCOP const,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	///@brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_score_function(
		TagCOP const,
		basic::datacache::DataMap const &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	///@brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_task_operations(
		TagCOP const,
		basic::datacache::DataMap const &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	///@brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP fresh_instance() const;

	///@brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP clone() const;

	// setters
	void score_function( ScoreFunctionCOP sf );
	void task_factory( TaskFactoryCOP tf );
	void task( PackerTaskCOP t );
	//void nloop( core::Size nloop_in );


	// accessors
	ScoreFunctionCOP score_function() const;
	PackerTaskCOP task() const;
	TaskFactoryCOP task_factory() const;

	void stochastic_pack( bool );
	bool stochastic_pack() const;

private:
	// pointers to data that are passed in
	ScoreFunctionCOP scorefxn_;
	PackerTaskCOP task_;
	TaskFactoryCOP task_factory_;

	bool stochastic_pack_; // calls stochastic_pack instead of min_pack if true
	bool nonideal_,cartesian_;
};


} // moves
} // protocols

#endif
