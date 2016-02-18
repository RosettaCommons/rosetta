// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/simple_moves/BBDihedralSamplerMover.hh
/// @brief Mover interface to BBDihedralSampler.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason W. Labonte (JWLabonte@jhu.edu)


#ifndef INCLUDED_protocols_simple_moves_BBDihedralSamplerMover_hh
#define INCLUDED_protocols_simple_moves_BBDihedralSamplerMover_hh

// Unit headers
#include <protocols/simple_moves/BBDihedralSamplerMover.fwd.hh>
#include <protocols/simple_moves/bb_sampler/BBDihedralSampler.hh>
#include <protocols/moves/Mover.hh>


#include <core/pose/Pose.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace simple_moves {

///@brief Mover interface to BBDihedralSampler ( 1D ).
/// Randomly samples a particular torsion type set in the sampler using residues from a movemap.
class BBDihedralSamplerMover : public protocols::moves::Mover {

public:

	BBDihedralSamplerMover();
	BBDihedralSamplerMover( bb_sampler::BBDihedralSamplerOP sampler);
	BBDihedralSamplerMover( bb_sampler::BBDihedralSamplerOP sampler, core::kinematics::MoveMapCOP mm);
	
	
	
	// copy constructor
	BBDihedralSamplerMover( BBDihedralSamplerMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~BBDihedralSamplerMover();
	
	virtual void
	apply( core::pose::Pose & pose );
	
public:
	
	///@brief Sets residues to sample on FROM movemap.
	void
	set_movemap( core::kinematics::MoveMapCOP movemap);
	
	///@brief Set a single resnum instead of a movemap. 
	void
	set_single_resnum( core::Size resnum );
	
	void
	set_sampler( bb_sampler::BBDihedralSamplerOP sampler) {
		sampler_ = sampler;
	}
	
	
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



	//BBDihedralSamplerMover & operator=( BBDihedralSamplerMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const;


private:

	bb_sampler::BBDihedralSamplerOP sampler_;
	utility::vector1< core::Size > bb_residues_; //This is faster than holding a movemap, as we will need to randomly sample on these residues.
	
};

std::ostream &operator<< (std::ostream &os, BBDihedralSamplerMover const &mover);


} //simple_moves
} //protocols


#endif //protocols/carbohydrates_BBDihedralSamplerMover_hh







