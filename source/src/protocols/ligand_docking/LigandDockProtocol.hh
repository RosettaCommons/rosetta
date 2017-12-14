// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/LigandDockProtocol.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_ligand_docking_LigandDockProtocol_hh
#define INCLUDED_protocols_ligand_docking_LigandDockProtocol_hh

#include <protocols/ligand_docking/LigandDockProtocol.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <protocols/ligand_docking/ResidueTorsionRestraints.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <map>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>
#include <utility/vector1.hh>


#ifdef WIN32
#include <protocols/ligand_docking/ResidueTorsionRestraints.hh>
#include <protocols/ligand_docking/MinimizeLigand.hh>
#endif

namespace protocols {
namespace ligand_docking {


/// @brief
///
/// @details
///
class LigandDockProtocol : public protocols::ligand_docking::LigandBaseProtocol
{
public:

	LigandDockProtocol();

	LigandDockProtocol(
		std::string const & protocol,
		bool const minimize_ligand,
		bool const minimize_backbone,
		bool const tether_ligand,
		bool const mutate_same_name3,
		core::Real const ligand_chi_stddev_deg,
		core::Real const protein_CA_stddev_Ang,
		core::Real const ligand_tether_stddev_Ang_=-1,
		core::Size const ligand_shear_moves=0
	);

	~LigandDockProtocol() override;


	LigandDockProtocolOP shared_from_this() { return utility::pointer::dynamic_pointer_cast<LigandDockProtocol>( Mover::shared_from_this() ); }

	void add_start_from(core::Real x, core::Real y, core::Real z);

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void
	append_ligand_docking_scores(
		core::pose::Pose const & before,
		core::pose::Pose const & after,
		core::scoring::ScoreFunctionCOP scorefxn,
		std::map< std::string, core::Real > & scores, //< appended to for output
		protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP constraint_io = nullptr
	) const;

private:

	LigandDockProtocol(LigandDockProtocol const & that);

	void
	classic_protocol(
		core::pose::Pose & pose,
		int jump_id,
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::moves::MonteCarloOP monteCarlo,
		core::Size num_cycles,
		core::Size repack_every_Nth
	) const;

	void
	shear_min_protocol(
		core::pose::Pose & pose,
		int jump_id,
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::moves::MonteCarloOP monteCarlo,
		core::Size num_cycles
	) const;

	void
	optimize_orientation3(
		core::pose::Pose & pose,
		int jump_id,
		core::Size lig_id
	);

	void
	random_conformer(
		core::pose::Pose & pose
	) const;

	protocols::moves::MoverOP
	make_dockmcm_mover(
		core::pose::Pose const & pose,
		int jump_id,
		protocols::moves::MoverOP repack_mover,
		protocols::moves::MoverOP rigbod_mover,
		core::kinematics::MoveMapOP movemap, //< would be COP but MinMover wants OP
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::moves::MonteCarloOP monteCarlo
	) const;

	void
	restrain_ligand_chis(
		core::pose::Pose & pose
	);

private:
	std::string protocol_; // i.e. 50 cycles instead of just 5
	bool minimize_ligand_;
	bool minimize_backbone_;
	bool tether_ligand_;
	bool ligand_protonation_;
	core::Real ligand_chi_stddev_deg_;
	core::Real protein_CA_stddev_Ang_;
	core::Real ligand_tether_stddev_Ang_;
	core::Size ligand_shear_moves_;
	bool minimize_all_rsds_;
	bool repack_all_rsds_;
	bool rottrials_all_rsds_;
	bool minimize_water_;
	utility::vector1< core::Vector > start_from_pts_; // may be empty
	utility::vector1< protocols::ligand_docking::ResidueTorsionRestraintsOP > ligand_torsion_restraints_;

}; // class LigandDockProtocol


} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_LigandDockProtocol_HH
