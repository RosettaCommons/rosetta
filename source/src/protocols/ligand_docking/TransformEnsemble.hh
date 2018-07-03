// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/ligand_docking/TransformEnsemble.hh
/// @author Thomas Willcock and Darwin Fu
/// Adapted from code by Sam Deluca

#ifndef INCLUDED_protocols_ligand_docking_TransformEnsemble_hh
#define INCLUDED_protocols_ligand_docking_TransformEnsemble_hh

#include <protocols/moves/Mover.hh>
#include <protocols/ligand_docking/TransformEnsemble.fwd.hh>
#include <protocols/qsar/scoring_grid/GridSet.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/UltraLightResidue.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/vector1.hh>
#include <vector>


namespace protocols {
namespace ligand_docking {


struct TransformEnsemble_info{ // including default values

public:
	//change first three into vectors (create a vector of chains)
	utility::vector1<std::string> chains;
	utility::vector1<core::Size> chain_ids;
	utility::vector1<core::Size> jump_ids;
	core::Real move_distance = 0;
	core::Real box_size = 0;
	core::Real angle = 0;
	core::Size cycles = 0;
	core::Real temperature; // no default
	core::Size repeats = 1;

	TransformEnsemble_info() = default;
	TransformEnsemble_info(TransformEnsemble_info const & ) = default;
};

class TransformEnsemble: public protocols::moves::Mover
{
public:
	TransformEnsemble();
	TransformEnsemble(
		protocols::qsar::scoring_grid::GridSetCOP grid_set_prototype,
		utility::vector1<std::string> const & chains,
		core::Real const & box_size,
		core::Real const & move_distance,
		core::Real const & angle,
		core::Size const & cycles,
		core::Real const & temp
	);
	TransformEnsemble(TransformEnsemble const & other);
	~TransformEnsemble() override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & pose
	) override;
	core::Vector weighted_center(utility::vector1<core::conformation::UltraLightResidue> & residues);
	void apply(core::pose::Pose & pose) override;

	void translate_ligand(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, utility::vector1<core::conformation::UltraLightResidue> & reference_residues, core::Real distance);
	void transform_ligand(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, utility::vector1<core::conformation::UltraLightResidue> & reference_residues);
	void change_conformers(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, const utility::vector1<core::conformation::UltraLightResidue> & reference_residues);
	void change_conformer(core::conformation::UltraLightResidue & ligand_residue, const core::conformation::UltraLightResidue & reference_residue, core::Size resid);
	void dump_conformer(core::conformation::UltraLightResidue const & residue, utility::io::ozstream & output);
	void print_xyz(core::Vector vector);
	bool monte_carlo(core::Real & current, core::Real & last);
	bool check_grid(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, core::Real distance = 0);

	/// @brief Calculate unconstrained score for given ligands
	core::Real score_ligands(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues);

	/// @brief Generate all extra grids
	void make_multi_pose_grids(core::Vector center);

	/// @brief copy the ligand into the desired receptor model and update conformation using the best model
	core::Real convert_to_full_pose(core::pose::Pose & pose, core::Size & best_pose_count);

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

private:
	protocols::qsar::scoring_grid::GridSetCOP grid_set_prototype_;
	utility::vector1<std::pair<std::string, qsar::scoring_grid::GridSetCOP>> grid_sets_;
	std::map<std::string, core::pose::Pose> grid_set_poses_;

	TransformEnsemble_info transform_info_;

	utility::vector1<core::conformation::UltraLightResidue> best_ligands_;
	utility::vector1<core::conformation::UltraLightResidue> ligand_residues_;
	utility::vector1<core::conformation::UltraLightResidue> reference_residues_;
	utility::vector1<core::conformation::UltraLightResidue> last_accepted_ligand_residues_;
	utility::vector1<core::conformation::UltraLightResidue> last_accepted_reference_residues_;

	std::map<core::Size, utility::vector1< core::conformation::ResidueOP > > ligand_conformers_;
	std::string ensemble_proteins_ = "";
	bool optimize_until_score_is_negative_ = false;
	bool output_sampled_space_ = false;
	bool use_conformers_ = true;
	bool use_main_model_ = false;


	std::string sampled_space_file_;
	core::Real initial_perturb_ = 0.0;


};

}
}

#endif /* TRANSFORM_HH_ */
