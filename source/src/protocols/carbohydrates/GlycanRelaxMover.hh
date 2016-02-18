// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/carbohydrates/GlycanRelaxMover.hh
/// @brief Main mover for Glycan Relax, which optimizes glycans in a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason W. Labonte (JWLabonte@jhu.edu)


#ifndef INCLUDED_protocols_carbohydrates_GlycanRelaxMover_hh
#define INCLUDED_protocols_carbohydrates_GlycanRelaxMover_hh

// Unit headers
#include <protocols/carbohydrates/GlycanRelaxMover.fwd.hh>
#include <protocols/carbohydrates/LinkageConformerMover.hh>
#include <protocols/moves/Mover.hh>


#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace carbohydrates {

///@brief Main mover for Glycan Relax, which optimizes glycans in a pose.
/// Each round optimizes either one residue for BB sampling, linkage, or multiple for minimization.
/// Currently uses a random sampler with a set of weights to each mover for sampling.
///
/// Weights are currently as follows:
///  .20 Phi Sugar BB Sampling
///  .20 Psi Sugar BB Sampling
///  .20 Linkage Conformer Sampling
///  .30 Small BB Sampling - equal weight to phi, psi, or omega
///    -> .17 +/- 15 degrees
///    -> .086 +/- 45 degrees
///    -> .044 +/- 90 degrees
///  .10 MinMover
///
class GlycanRelaxMover : public protocols::moves::Mover {

public:

	GlycanRelaxMover();
	
	//@brief constructor with arguments
	GlycanRelaxMover( core::kinematics::MoveMapCOP mm,
					  core::scoring::ScoreFunctionCOP scorefxn,
					  core::Size rounds = 50);
	
	// copy constructor
	GlycanRelaxMover( GlycanRelaxMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~GlycanRelaxMover();

	virtual void
	apply( core::pose::Pose & pose );

public:


	
	void
	set_movemap(core::kinematics::MoveMapCOP movemap);
	
	void
	set_scorefunction( core::scoring::ScoreFunctionCOP scorefxn);
	
	///@brief Each round applys a random mover to pose.
	/// This setting is multiplied by the number of glycan residues in the movemap for the total number of rounds
	void
	set_rounds( core::Size rounds);
	
	void
	set_kt( core::Real kt);
	
	void
	set_defaults();
	
public:
	virtual void
	show( std::ostream & output=std::cout ) const;

	std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );



	//GlycanRelaxMover & operator=( GlycanRelaxMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const;

private:
	
	void
	init_objects( core::pose::Pose const & pose );
	
	void
	set_cmd_line_defaults();
	
	void
	apply_to_res(
		core::pose::Pose & pose,
		core::Size resnum,
		core::kinematics::MoveMapOP mm,
		core::scoring::ScoreFunctionOP score,
		core::Size round);
	
private:

	core::kinematics::MoveMapOP full_movemap_;
	core::kinematics::MoveMapOP glycan_movemap_;
	
	moves::MonteCarloOP mc_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	
	LinkageConformerMoverOP linkage_mover_;
	moves::RandomMoverOP weighted_random_mover_;
	simple_moves::MinMoverOP min_mover_;
	core::Size rounds_;
	core::Real kt_;
	
	utility::vector1<std::string> accept_log_;
	
	bool test_;
	bool final_min_;
	bool pack_glycans_; //Problem packing glycans for now.
	
	core::Size total_glycan_residues_;
	bool pymol_movie_;
	
	std::string ref_pose_name_;
	bool use_branches_;
	
	utility::vector1< std::string > parsed_positions_;
	utility::vector1< core::Size > positions_;
};

std::ostream &operator<< (std::ostream &os, GlycanRelaxMover const &mover);


} //protocols
} //carbohydrates


#endif //protocols/carbohydrates_GlycanRelaxMover_hh







