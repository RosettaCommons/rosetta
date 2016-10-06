// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/ligand_docking/Transform.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_ligand_docking_Transform_hh
#define INCLUDED_protocols_ligand_docking_Transform_hh

#include <protocols/moves/Mover.hh>
#include <protocols/ligand_docking/Transform.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/UltraLightResidue.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {


struct Transform_info{ // including default values

public:
	std::string chain;
	core::Size chain_id;
	core::Size jump_id;
	core::Real move_distance;
	core::Real box_size;
	core::Real angle;
	core::Real rmsd;
	core::Size cycles;
	core::Real temperature;
	core::Size repeats;
	Transform_info():
		chain(""),
		chain_id(0),
		jump_id(0),
		move_distance(0.0),
		box_size(0.0),
		angle(0.0),
		rmsd(0.0),
		cycles(0),
		temperature(0.0),
		repeats(1){};
	Transform_info(Transform_info const & ) = default;
};

class Transform: public protocols::moves::Mover
{
public:
	Transform();
	Transform(
		std::string const & chain,
		core::Real const & box_size,
		core::Real const & move_distance,
		core::Real const & angle,
		core::Size const & cycles,
		core::Real const & temp
	);
	Transform(Transform const & other);
	~Transform() override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	std::string get_name() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	void apply(core::pose::Pose & pose) override;

private:

	/// @brief Performa a randomization of the ligand residue, translating by some random value within a uniform distribution with a max of distance,
	/// and rotating around a random axis with a uniformly random angle of between -angle/2 and +angle/2 (in degrees).
	void randomize_ligand(core::conformation::UltraLightResidue & residue, core::Real distance, core::Real angle);

	/// @brief translate and rotate a random value by the distances specified in the Transform_info object, using a gaussian distribution
	void transform_ligand(core::conformation::UltraLightResidue & residue);

	/// @brief randomly change the ligand conformation
	void change_conformer(core::conformation::UltraLightResidue & residue);

	/// @brief output the ligand residues to a pdb file
	void dump_conformer(core::conformation::UltraLightResidue & residue, utility::io::ozstream & output);

	/// @brief return true if the rmsd is within the specified cutoff
	bool check_rmsd(core::conformation::UltraLightResidue const & start, core::conformation::UltraLightResidue const& current) const;

private:
	//qsar::scoring_grid::GridManagerOP grid_manager_;
	Transform_info transform_info_;
	utility::vector1< core::conformation::ResidueOP >  ligand_conformers_;
	bool optimize_until_score_is_negative_ = false;
	bool output_sampled_space_ = false;
	bool check_rmsd_ = false;
	std::string sampled_space_file_;
	core::Real initial_perturb_ = 0.0;
	core::Real initial_angle_perturb_ = -360.0;

};

}
}

#endif /* TRANSFORM_HH_ */
