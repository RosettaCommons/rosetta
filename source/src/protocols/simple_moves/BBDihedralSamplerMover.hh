// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/BBDihedralSamplerMover.hh
/// @brief Mover interface to BBDihedralSampler.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason W. Labonte (JWLabonte@jhu.edu)


#ifndef INCLUDED_protocols_simple_moves_BBDihedralSamplerMover_hh
#define INCLUDED_protocols_simple_moves_BBDihedralSamplerMover_hh

// Unit headers
#include <protocols/simple_moves/BBDihedralSamplerMover.fwd.hh>
#include <protocols/simple_moves/bb_sampler/BBDihedralSampler.hh>
#include <protocols/moves/Mover.hh>


#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace simple_moves {

///@brief Mover interface to BBDihedralSampler ( 1D ).
/// Samples using BBDihedralSamplers randomly using residues set in the movemap.
/// If multiple samplers are given, will randomly sample on them.
///
///@details Obeys Movemap BB Torsion ID settings if given in movemap.
///
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

	///@brief Sets residues (and optionally torsions) to sample on FROM movemap.
	void
	set_movemap( core::kinematics::MoveMapCOP movemap);

	///@brief Set a single resnum instead of a movemap.
	void
	set_single_resnum( core::Size resnum );

	///@brief Set a single sampler this mover will use.
	void
	set_sampler( bb_sampler::BBDihedralSamplerOP sampler );

	///@brief Add a sampler to this mover.
	void
	add_sampler( bb_sampler::BBDihedralSamplerOP sampler );


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

	//@brief Sets the union of residues available in movemap and sampler torsion ids.
	void
	setup_sampler_movemap_union( core::pose::Pose const & pose );

	///@brief Sets up bb_residues_ variable as all residues in the pose.
	void
	setup_all_bb_residues( core::pose::Pose const & pose );

private:

	//bb_sampler::BBDihedralSamplerOP sampler_;

	std::map< core::Size, utility::vector1< bb_sampler::BBDihedralSamplerOP > > samplers_;
	utility::vector1< core::Size > sampler_torsion_types_;

	core::kinematics::MoveMapCOP movemap_;
	utility::vector1< core::Size > bb_residues_; //This is faster than turning the movemap into a vector each apply, as we will need to randomly sample on these residues.

	std::map< core::Size, utility::vector1< core::Size > > sampler_movemap_union_; //Union of the torsion samplers we have and torsion IDs on in the movemap.

};

std::ostream &operator<< (std::ostream &os, BBDihedralSamplerMover const &mover);


} //simple_moves
} //protocols


#endif //protocols/carbohydrates_BBDihedralSamplerMover_hh







