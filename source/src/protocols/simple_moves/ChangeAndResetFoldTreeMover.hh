// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/ChangeAndResetFoldTreeMover.hh
/// @brief Basic Mover used for setting up atomic protocols.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_simple_moves_ChangeAndResetFoldTreeMover_hh
#define INCLUDED_protocols_simple_moves_ChangeAndResetFoldTreeMover_hh

#include <protocols/simple_moves/ChangeAndResetFoldTreeMover.fwd.hh>

#include <protocols/moves/ChangeFoldTreeMover.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverApplyingMover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>

// Forward
namespace protocols {
namespace simple_moves {

/// @brief Basic mover used primarily for setting up atomic protocols, especially in generator functions.
/// 
/// As the Topology Broker system includes more claims and is more widely adopted,
/// This should not be used in favor of the Broker system.
///
/// @details
///  1) Holds the Pose's original FoldTree
///  2) Applies its ChangeFoldTreeMover
///  3) Applies its Main mover
///  4) Returns the Pose's original FoldTree
///
class ChangeAndResetFoldTreeMover : public protocols::moves::MoverApplyingMover {
public:

	/// @brief Basic Constructor.  
	ChangeAndResetFoldTreeMover();
	

	ChangeAndResetFoldTreeMover(protocols::moves::MoverOP main_mover);
	

	ChangeAndResetFoldTreeMover(
		protocols::moves::MoverOP main_mover,
		protocols::moves::ChangeFoldTreeMoverOP ft_mover);
	
	/// @brief Constructor with all needed classes as well as the optional scorefxn.
	/// Note that the scorefxn is ONLY used for printing score information in apply.
	///
	ChangeAndResetFoldTreeMover(
		protocols::moves::MoverOP main_mover,
		protocols::moves::ChangeFoldTreeMoverOP ft_mover,
		core::scoring::ScoreFunctionCOP scorefxn);
	
	ChangeAndResetFoldTreeMover(ChangeAndResetFoldTreeMover const & src);

	virtual ~ChangeAndResetFoldTreeMover();

public:
	
	/// @brief Set the main mover to apply
	virtual void
	set_mover(protocols::moves::MoverOP main_mover);
	
	/// @brief Get the main mover
	virtual protocols::moves::MoverOP
	mover() const;
	
	/// @brief Set the ChangeFoldTreeMover
	void
	set_ft_mover(protocols::moves::ChangeFoldTreeMoverOP ft_mover);
	
	/// @brief Get the ChangeFoldTreeMover
	protocols::moves::ChangeFoldTreeMoverOP
	ft_mover() const;
		

	/// @details
	///  1) Holds the Pose's original FoldTree
	///  2) Applies the ChangeFoldTreeMover
	///  3) Applies the Main mover
	///  4) Returns the Pose's original FoldTree
	virtual void
	apply(core::pose::Pose & pose);

public:
	
	/// @brief Set the ScoreFunction.
	/// Used ONLY for printing score information before and after the move!!
	void
	set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn);
	
	/// @brief Get the set ScoreFunction
	core::scoring::ScoreFunctionCOP
	scorefxn()  const;
	
public:

	std::string
	get_name() const;

	//void parse_my_tag(
	//	TagCOP tag,
	//	basic::datacache::DataMap & data,
	//	Filters_map const & filters,
	//	Movers_map const & movers,
	//	Pose const & pose
	//	);

	protocols::moves::MoverOP
	clone() const;

	//ChangeAndResetFoldTreeMover & operator=(ChangeAndResetFoldTreeMover const & src);

	virtual protocols::moves::MoverOP fresh_instance() const;

private:
	protocols::moves::MoverOP main_mover_;
	protocols::moves::ChangeFoldTreeMoverOP ft_mover_;
	core::scoring::ScoreFunctionCOP scorefxn_;

};

} //simple_moves
} //protocols


#endif //INCLUDED_protocols_simple_moves_ChangeAndResetFoldTreeMover_hh
