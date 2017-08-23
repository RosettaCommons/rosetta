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

#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace simple_moves {

///@brief Mover interface to BBDihedralSampler ( 1D ).
/// Samples using BBDihedralSamplers randomly using residues set in the movemap.
/// If multiple samplers are given, will randomly sample on them.
///
class BBDihedralSamplerMover : public protocols::moves::Mover {

public:

	BBDihedralSamplerMover();
	BBDihedralSamplerMover( bb_sampler::BBDihedralSamplerOP sampler);
	BBDihedralSamplerMover( bb_sampler::BBDihedralSamplerOP sampler, core::select::residue_selector::ResidueSelectorCOP selector);



	// copy constructor
	BBDihedralSamplerMover( BBDihedralSamplerMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	~BBDihedralSamplerMover() override;

	void
	apply( core::pose::Pose & pose ) override;

public:

	///@brief Sets residues from a residue selector.
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector);

	///@brief Set a single resnum instead of a movemap.
	void
	set_single_resnum( core::Size resnum );

	///@brief Set a single sampler this mover will use.
	void
	set_sampler( bb_sampler::BBDihedralSamplerOP sampler );

	///@brief Add a sampler to this mover.
	void
	add_sampler( bb_sampler::BBDihedralSamplerOP sampler );


	///@brief Limit the torsions that are sampled using a mask.
	/// Mask is a vector of torsions we are allowed to sample at each residue.
	///  NOT every residue must be present - only the ones we need to mask.
	///  This helps for carbohydrates, since not every residue has all possible dihedrals, even if we have samplers to use.
	void
	set_dihedral_mask( std::map< core::Size, utility::vector1< core::Size >> mask );

public:

	void
	show( std::ostream & output=std::cout ) const override;

	std::string
	get_name() const override;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;



	//BBDihedralSamplerMover & operator=( BBDihedralSamplerMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

private:

	//@brief Sets the union of residues available in movemap and sampler torsion ids.
	void
	setup_samplers( core::pose::Pose const & pose );

	///@brief Sets up bb_residues_ variable as all residues in the pose.
	void
	setup_all_bb_residues( core::pose::Pose const & pose );

private:

	//bb_sampler::BBDihedralSamplerOP sampler_;

	std::map< core::Size, utility::vector1< bb_sampler::BBDihedralSamplerOP > > samplers_;
	utility::vector1< core::Size > sampler_torsion_types_;

	core::select::residue_selector::ResidueSelectorOP selector_ = nullptr;
	utility::vector1< core::Size > bb_residues_; //This is faster than turning the movemap into a vector each apply, as we will need to randomly sample on these residues.

	std::map< core::Size, utility::vector1< core::Size > > sampler_torsions_; //Resnum, vector of torsion IDs to sample.
	std::map< core::Size, utility::vector1< core::Size >> dihedral_mask_;

};

std::ostream &operator<< (std::ostream &os, BBDihedralSamplerMover const &mover);


} //simple_moves
} //protocols


#endif //protocols/carbohydrates_BBDihedralSamplerMover_hh







