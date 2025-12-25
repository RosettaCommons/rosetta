// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon , adapted from the ResfileReader code
/// @author Steven Lewis (smlewi@gmail.com)
/// @author Andrew Leaver-Fay
/// @author Reworked a bit by Vikram K. Mulligan (vmulligan@flatironinstitue.org) during the 2019 PackerPalette mini-XRW to use PackerPalettes.

// Unit Headers
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/HighResDockerCreator.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/ligand_options/Interface.hh>
#include <protocols/ligand_docking/MinimizeLigand.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/TetherLigand.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>

#include <protocols/ligand_docking/UnconstrainedTorsionsMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/palette/PackerPalette.hh>
#include <core/pack/palette/CustomBaseTypePackerPalette.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResidueLevelTask.hh>

// Scripter Headers
#include <utility>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

//STL headers
#include <string>

#include <set>

//Auto Headers
#include <protocols/ligand_docking/LigandArea.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/CitationManager.hh>

#include <core/pack/task/ResidueLevelTask.hh> // AUTO IWYU For ResidueLevelTask

namespace protocols {
namespace ligand_docking {

static basic::Tracer TR( "protocols.ligand_docking.ligand_options.Protocol" );

/// @brief
HighResDocker::HighResDocker():
	Mover("HighResDocker"),
	num_cycles_(0),
	repack_every_Nth_(0),
	score_fxn_(/* NULL */),
	movemap_builder_(/* NULL */),
	resfile_("")
{
	resfile_.clear();
	// Now use cycles and repack_every_Nth to replicate these options...
	//meiler2006: 50, 8;
	//abbreviated: 5, 4;
	//abbrev2: 6, 3;
	initialize_from_options();
}

//With chains though chains do not appear to be used anywhere
HighResDocker::HighResDocker(
	core::Size num_cycles,
	core::Size repack_every_Nth,
	core::scoring::ScoreFunctionOP score_fxn,
	MoveMapBuilderOP movemap_builder,
	std::string const & resfile
): num_cycles_(num_cycles), repack_every_Nth_(repack_every_Nth), score_fxn_(std::move(score_fxn)), movemap_builder_(std::move(movemap_builder)), resfile_( resfile ){
	initialize_from_options();
}

HighResDocker::HighResDocker(HighResDocker const & that):
	//utility::VirtualBase(),
	protocols::moves::Mover( that ),
	num_cycles_(that.num_cycles_),
	repack_every_Nth_(that.repack_every_Nth_),
	score_fxn_(that.score_fxn_),
	movemap_builder_(that.movemap_builder_),
    rmsd_max_(that.rmsd_max_),
    check_rmsd_(that.check_rmsd_)
{
	initialize_from_options();
}

HighResDocker::~HighResDocker() = default;

protocols::moves::MoverOP HighResDocker::clone() const {
	return utility::pointer::make_shared< HighResDocker >( *this );
}

protocols::moves::MoverOP HighResDocker::fresh_instance() const {
	return utility::pointer::make_shared< HighResDocker >();
}

/// @brief function to be called by constructors to assign values to allow_minimization_, allow_repacking, and all_residues
void
HighResDocker::initialize_from_options(){
	bool allow_minimization_ = basic::options::option[ basic::options::OptionKeys::motifs::highresdocker_allow_minimization ];
	bool allow_repacking_ = basic::options::option[ basic::options::OptionKeys::motifs::highresdocker_allow_repacking ];
	bool all_residues_ = basic::options::option[ basic::options::OptionKeys::motifs::highresdocker_use_all_residues ];
	TR << "Initializing options: " << std::endl;
	TR << "allow_minimization_: " << allow_minimization_ << std::endl;
	TR << "allow_repacking_: " << allow_repacking_ << std::endl;
	TR << "all_residues_: " << all_residues_ << std::endl;
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
HighResDocker::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	if ( tag->getName() != "HighResDocker" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible");
	}

	// cycles and repack_every_Nth
	if ( ! tag->hasOption("cycles") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'HighResDocker' mover requires cycles tag");
	if ( ! tag->hasOption("repack_every_Nth") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'HighResDocker' mover requires repack_every_Nth tag");
	num_cycles_= tag->getOption<core::Size>("cycles");
	repack_every_Nth_= tag->getOption<core::Size>("repack_every_Nth");

	/// Score Function ///
	if ( ! tag->hasOption("scorefxn") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'HighResDocker' requires 'scorefxn' tag");
	std::string scorefxn_name= tag->getOption<std::string>("scorefxn");
	score_fxn_= datamap.get_ptr<core::scoring::ScoreFunction>( "scorefxns", scorefxn_name);

	/// MoveMapBuilder///
	if ( ! tag->hasOption("movemap_builder") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'HighResDocker' requires 'movemap_builder' tag");
	std::string movemap_builder_name= tag->getOption<std::string>("movemap_builder");
	movemap_builder_= datamap.get_ptr<MoveMapBuilder>( "movemap_builders", movemap_builder_name);

	/// Resfile ///
	if ( tag->hasOption("resfile") ) {
		resfile_= tag->getOption<std::string>("resfile");
	}

    /// RMSD ///
    if ( tag->hasOption("rmsd") ) {
        check_rmsd_ = true;
        rmsd_max_ = tag->getOption<core::Real>("rmsd");
    }
}


MinimizeLigandOPs
HighResDocker::setup_ligands_to_minimize(
	core::pose::Pose & pose,
	char chain // =0 This is for dealing with multiple ligand docking
){
	MinimizeLigandOPs minimize_ligands;

	LigandAreas const ligand_areas =
		movemap_builder_->get_sc_interface_builder()->get_ligand_areas();
	auto ligand_area_itr= ligand_areas.begin();
	LigandAreas::const_iterator const ligand_area_end= ligand_areas.end();

	//TODO Use BOOST_FOREACH


	//if chain is 0, then loop proceeds normally. If chain is given, then minimize ligand only runs for the matching chain.
	for ( ; ligand_area_itr != ligand_area_end; ++ligand_area_itr ) {
		if ( ligand_area_itr->first == chain || chain == 0 ) {
			char const & input_chain= ligand_area_itr->first;
			LigandAreaOP const ligand_area( ligand_area_itr->second );
			core::Real const & degrees( ligand_area->minimize_ligand_ );
			if ( degrees > 0 ) {
				MinimizeLigandOP minimize_ligand( new MinimizeLigand(input_chain, degrees) );
				minimize_ligand->apply(pose);
				minimize_ligands.push_back(minimize_ligand);
			}
		}
	}
	return minimize_ligands;
}

void
HighResDocker::remove_ligand_dihedral_restraints(core::pose::Pose & pose, MinimizeLigandOPs & minimized_ligands) const {
	for ( core::Size ii(1); ii <= minimized_ligands.size(); ++ii ) {
		minimized_ligands[ii]->remove_constraints( pose );
	}
}

TetherLigandOPs
HighResDocker::tether_ligands(
	core::pose::Pose & pose,
	char chain // =0 This is for dealing with multiple ligand docking
){
	TetherLigandOPs ligand_tethers;

	LigandAreas const ligand_areas =
		movemap_builder_->get_sc_interface_builder()->get_ligand_areas();
	auto ligand_area_itr= ligand_areas.begin();
	LigandAreas::const_iterator const ligand_area_end= ligand_areas.end();

	//if chain is 0, then loop proceeds normally. If chain is given, then minimize ligand only runs for the matching chain.

	for ( ; ligand_area_itr != ligand_area_end; ++ligand_area_itr ) {
		if ( ligand_area_itr->first == chain || chain == 0 ) {

			char const & input_chain= ligand_area_itr->first;
			LigandAreaOP const ligand_area( ligand_area_itr->second );
			core::Real const & tether_size( ligand_area->tether_ligand_ );
			if ( tether_size > 0 ) {
				TetherLigandOP tether_ligand( new TetherLigand(input_chain, tether_size) );

				tether_ligand->apply(pose);
				ligand_tethers.push_back(tether_ligand);
			}
		}
	}
	return ligand_tethers;
}

void
HighResDocker::remove_ligand_tethers(core::pose::Pose pose, TetherLigandOPs ligand_tethers){
	// TetherLigandOPs::const_iterator begin= ligand_tethers.begin(); // Unused variable causes warning.
	// TetherLigandOPs::const_iterator const end= ligand_tethers.end(); // Unused variable causes warning.

	for ( TetherLigandOP ligand_tether : ligand_tethers ) {
		ligand_tether->release(pose);
	}
}

void
HighResDocker::apply(core::pose::Pose & pose) {
	debug_assert(num_cycles_ > 0);

	MinimizeLigandOPs minimized_ligands( setup_ligands_to_minimize(pose) );

	TetherLigandOPs ligand_tethers= tether_ligands(pose);

	debug_assert(movemap_builder_ && score_fxn_ ); // make sure the pointers point
	core::kinematics::MoveMapOP movemap( movemap_builder_->build(pose) );

	protocols::moves::MonteCarloOP monteCarlo( new protocols::moves::MonteCarlo(pose, *score_fxn_, 2.0) );/* temperature, from RosettaLigand paper */
	score_fxn_->score( pose ); // without this neither of the movers below were working
	// I believe that this may have been related to adding constraints incorrectly at other places in my code.
	// Rigid body exploration
	utility::vector1<protocols::moves::MoverOP> rigid_body_movers = create_rigid_body_movers(pose);

    core::Size ligand_index = 0;
    for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
        core::conformation::Residue const & residue( pose.residue(i) );
        if ( residue.is_ligand() ) {
            ligand_index = i;
        }
    }
    core::conformation::Residue original_ligand(pose.residue(ligand_index));
    TR.Debug << "Save residue number " << ligand_index << " as ligand." << std::endl;

    core::Size rejected_rmsd_moves = 0;

	for ( core::Size cycle = 1; cycle <= num_cycles_; ++cycle ) {

        TR.Debug << "Start cycle " << cycle << ".";
        if ( repack_every_Nth_ != 0 && cycle % repack_every_Nth_ == 1 ) TR.Debug << " Repacking!";
        TR.Debug << std::endl;

		core::pack::task::PackerTaskOP packer_task;

		packer_task = make_packer_task(pose, all_residues_);

		protocols::moves::MoverOP pack_mover;
		//for no repacking, make sure to go to else
		if ( repack_every_Nth_ != 0 && cycle % repack_every_Nth_ == 1 ) {
			TR.Debug << "making PackRotamersMover" << std::endl;
			pack_mover = utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >(score_fxn_, packer_task);
		} else {
			TR.Debug << "making RotamerTrialsMover" << std::endl;
			pack_mover = utility::pointer::make_shared< protocols::minimization_packing::RotamerTrialsMover >(score_fxn_, *packer_task);
		}
		// Wrap it in something to disable the torsion constraints before packing!
		pack_mover = utility::pointer::make_shared< protocols::ligand_docking::UnconstrainedTorsionsMover >( pack_mover, minimized_ligands );

		protocols::minimization_packing::MinMoverOP min_mover( new protocols::minimization_packing::MinMover( movemap, score_fxn_, "lbfgs_armijo_nonmonotone_atol", 1.0, true /*use_nblist*/ ) );
		min_mover->min_options()->nblist_auto_update(true); // does this cost us lots of time in practice?

		core::Real const score1 = (*score_fxn_)( pose );
		apply_rigid_body_moves(pose, rigid_body_movers);

		//do not pack if not allowed
		if ( allow_repacking_ ) {
			pack_mover->apply(pose);
		}


		core::Real const score2 = (*score_fxn_)( pose );
		if ( score2 - score1 < 15.0 ) {
			min_mover->apply(pose);
		}

        if ( check_rmsd_ ) {
            core::Real rmsd = core::scoring::automorphic_rmsd( original_ligand, pose.residue(ligand_index), /*superimpose=*/false );
            TR.Debug << "Current rmsd is " << rmsd << ". Max rmsd: " << rmsd_max_ << std::endl;
            if (rmsd > rmsd_max_){
                pose = monteCarlo->last_accepted_pose();
                rejected_rmsd_moves++;
            }
        }
		monteCarlo->boltzmann( pose );

	}

    if (check_rmsd_) {
        TR << "Rejected " << rejected_rmsd_moves << " of " << num_cycles_ << " pose updates due to rmsd limitations." << std::endl;
    }

	// keep the best structure we found, not the current one
	monteCarlo->show_scores();
	monteCarlo->recover_low(pose);

	remove_ligand_tethers(pose, ligand_tethers);
	remove_ligand_dihedral_restraints(pose, minimized_ligands);
}

void
HighResDocker::apply(utility::vector1<core::pose::Pose> & poses, utility::vector1<core::Real> & current_scores, utility::vector1<char> qsar_chars, core::Size cycle) {

	debug_assert(num_cycles_ > 0);

	core::Size pose_counter = 1;

	for ( core::pose::Pose pose : poses ) {

		MinimizeLigandOPs minimized_ligands( setup_ligands_to_minimize(pose, qsar_chars[pose_counter]) );

		TetherLigandOPs ligand_tethers= tether_ligands(pose, qsar_chars[pose_counter]);

		debug_assert(movemap_builder_ && score_fxn_ ); // make sure the pointers point
		core::kinematics::MoveMapOP movemap( movemap_builder_->build(pose) );

		score_fxn_->score( pose ); // without this neither of the movers below were working

		// I believe that this may have been related to adding constraints incorrectly at other places in my code.
		// Rigid body exploration
		utility::vector1<protocols::moves::MoverOP> rigid_body_movers= create_rigid_body_movers(pose);

		core::pack::task::PackerTaskOP packer_task;

		packer_task = make_packer_task(pose, all_residues_);

		protocols::moves::MoverOP pack_mover;

		if ( repack_every_Nth_ != 0 && cycle % repack_every_Nth_ == 1 ) {
			TR.Debug << "making PackRotamersMover" << std::endl;
			pack_mover = utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >(score_fxn_, packer_task);
		} else {
			TR.Debug << "making RotamerTrialsMover" << std::endl;
			pack_mover = utility::pointer::make_shared< protocols::minimization_packing::RotamerTrialsMover >(score_fxn_, *packer_task);
		}

		// Wrap it in something to disable the torsion constraints before packing!
		pack_mover = utility::pointer::make_shared< protocols::ligand_docking::UnconstrainedTorsionsMover >( pack_mover, minimized_ligands );

		protocols::minimization_packing::MinMoverOP min_mover( new protocols::minimization_packing::MinMover( movemap, score_fxn_, "dfpmin_armijo_nonmonotone_atol", 1.0, true /*use_nblist*/ ) );
		min_mover->min_options()->nblist_auto_update(true); // does this cost us lots of time in practice?

		apply_rigid_body_moves(pose, rigid_body_movers);
		// not repack if not allowed
		if ( allow_repacking_ ) {
			pack_mover->apply(pose);
		}

		if ( allow_minimization_ ) {
			min_mover->apply(pose);  //okay to always apply MinMover?
		}


		current_scores[pose_counter] = (*score_fxn_)( pose );
		TR << "HighResDocker Pose Score For " << pose_counter << " is " << current_scores[pose_counter] << "\n";
		++pose_counter;

	}

}

void
HighResDocker::apply(core::pose::Pose & pose, core::Real & current_score, char qsar_char, core::Size cycle) {

	MinimizeLigandOPs minimized_ligands( setup_ligands_to_minimize(pose, qsar_char) );

	TetherLigandOPs ligand_tethers= tether_ligands(pose, qsar_char);

	debug_assert(movemap_builder_ && score_fxn_ ); // make sure the pointers point
	core::kinematics::MoveMapOP movemap( movemap_builder_->build(pose) );

	score_fxn_->score( pose ); // without this neither of the movers below were working

	// I believe that this may have been related to adding constraints incorrectly at other places in my code.
	// Rigid body exploration
	utility::vector1<protocols::moves::MoverOP> rigid_body_movers= create_rigid_body_movers(pose);

	core::pack::task::PackerTaskOP packer_task;

	packer_task = make_packer_task(pose, all_residues_);

	protocols::moves::MoverOP pack_mover;

	if ( repack_every_Nth_ != 0 && cycle % repack_every_Nth_ == 1 ) {
		TR.Debug << "making PackRotamersMover" << std::endl;
		pack_mover = utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >(score_fxn_, packer_task);
	} else {
		TR.Debug << "making RotamerTrialsMover" << std::endl;
		pack_mover = utility::pointer::make_shared< protocols::minimization_packing::RotamerTrialsMover >(score_fxn_, *packer_task);
	}

	// Wrap it in something to disable the torsion constraints before packing!
	pack_mover = utility::pointer::make_shared< protocols::ligand_docking::UnconstrainedTorsionsMover >( pack_mover, minimized_ligands );

	protocols::minimization_packing::MinMoverOP min_mover( new protocols::minimization_packing::MinMover( movemap, score_fxn_, "dfpmin_armijo_nonmonotone_atol", 1.0, true /*use_nblist*/ ) );
	min_mover->min_options()->nblist_auto_update(true); // does this cost us lots of time in practice?

	core::Real const score1 = (*score_fxn_)( pose );
	apply_rigid_body_moves(pose, rigid_body_movers);

	if ( allow_repacking_ ) {
		pack_mover->apply(pose);
	}


	//skip this if you do not want to minimize
	if ( allow_minimization_ ) {
		core::Real const score2 = (*score_fxn_)( pose );
		if ( score2 - score1 < 15.0 ) {
			min_mover->apply(pose);
		}
	}


	current_score = (*score_fxn_)( pose );


}


core::pack::task::PackerTaskOP
HighResDocker::make_packer_task_from_vector(
	core::pose::Pose const & pose,
	ligand_options::Interface const & allow_repack
) const{
	static bool pose_already_packed= false;

	// Initialize the PackerPalette:
	core::pack::palette::CustomBaseTypePackerPaletteOP palette( utility::pointer::make_shared< core::pack::palette::CustomBaseTypePackerPalette >() );
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		core::chemical::ResidueType const & this_res_type( pose.residue_type(i) );
		if ( this_res_type.is_ligand() ) {
			if ( !palette->has_type( this_res_type.base_name() ) ) {
				palette->add_type( this_res_type.base_name() );
			}

			//Add all other residue types with the same name3:
			core::chemical::ResidueTypeSetCOP rsd_type_set( pose.residue_type_set_for_pose( this_res_type.mode() ) );
			core::chemical::ResidueTypeCOPs allowed_types( core::chemical::ResidueTypeFinder( *rsd_type_set ).name3( this_res_type.name3() ).get_all_possible_residue_types() ); // a vector1
			for ( core::Size i(1), imax(allowed_types.size()); i<=imax; ++i ) {
				if ( !palette->has_type( allowed_types[i]->base_name() ) ) {
					palette->add_type( allowed_types[i]->base_name() );
				}
			}
		}
	}

	core::pack::task::PackerTaskOP pack_task( core::pack::task::TaskFactory::create_packer_task(pose, palette) );
	pack_task->initialize_from_command_line(); // -ex1 -ex2  etc.

	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot_( new core::pack::rotamer_set::UnboundRotamersOperation() );
	unboundrot_->initialize_from_command_line();
	pack_task->append_rotamerset_operation( unboundrot_ );

	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		/// If several params files have the same name, allow switching among them
		/// This was previously only enabled with mutate_same_name3.  Now default.
		if ( ! pose.residue(i).is_ligand() ) continue;
		TR.Debug<<  "enabling packing for ligand residue "<< i << std::endl;
		enable_ligand_rotamer_packing(pose, i, pack_task);
	}

	if ( resfile_.empty() ) {
		bool const use_resfile= basic::options::option[ basic::options::OptionKeys::packing::resfile ].user() ;
		if ( use_resfile ) {
			TR.Debug << "using OPTIONS resfile"<< std::endl;
			core::pack::task::parse_resfile(pose, *pack_task);
		} else {
			TR.Debug << "restricting to repack"<< std::endl;
			for ( core::Size i = 1; i <= pose.size(); ++i ) {
				if ( ! pose.residue(i).is_ligand() ) {
					pack_task->nonconst_residue_task( i ).restrict_to_repacking();
				}
			}
		}
	} else {
		TR.Debug << "using XML resfile"<< std::endl;
		core::pack::task::parse_resfile(pose, *pack_task, resfile_);
	}


	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( allow_repack[i].type == ligand_options::InterfaceInfo::non_interface  ) {
			pack_task->nonconst_residue_task( i ).prevent_repacking();
		}
	}

	// We always want the option (after the initial unbiased pack)
	// of sticking with our current nicely minimized conformation.
	if ( pose_already_packed ) {
		pack_task->or_include_current(true);
	} else {
		pose_already_packed=true;
	}

	return pack_task;
}

core::pack::task::PackerTaskOP
HighResDocker::make_packer_task(
	core::pose::Pose const & pose,
	bool all_residues
) const{

	//original value of all_residues was intended to be false
	if ( all_residues ) {
		ligand_options::Interface interface(pose.size(), ligand_options::InterfaceInfo(ligand_options::InterfaceInfo::is_interface)); // 0 is false, #
		return make_packer_task_from_vector(pose, interface);
	} else { // the packer task interface should match the movemap interface
		InterfaceBuilderOP sc_interface_builder= movemap_builder_->get_sc_interface_builder();
		ligand_options::Interface side_chain_interface= sc_interface_builder->build(pose);
		return make_packer_task_from_vector(pose, side_chain_interface);
	}
}

void
HighResDocker::enable_ligand_rotamer_packing(
	core::pose::Pose const & pose,
	core::Size const ligand_residue_id,
	core::pack::task::PackerTaskOP & pack_task
) const{
	core::conformation::Residue const & this_residue= pose.residue(ligand_residue_id);

	core::chemical::ResidueTypeSetCOP rsd_type_set = pose.residue_type_set_for_pose( this_residue.type().mode() );
	core::chemical::ResidueTypeCOPs allowed_types = core::chemical::ResidueTypeFinder( *rsd_type_set ).name3( this_residue.name3() ).get_all_possible_residue_types(); // a vector1


	debug_assert(allowed_types.size() > 0);
	/// TODO consider removing this so resfiles can specify ligand mutations to allow
	if ( allowed_types.size() == 1 ) {
		pack_task->nonconst_residue_task( ligand_residue_id ).restrict_to_repacking();
		return;
	}
	// else

	utility::vector1< std::string > types_to_allow;
	types_to_allow.reserve( allowed_types.size() + 1 ); //Max possible size
	types_to_allow.push_back( this_residue.type().base_name() );
	for ( core::Size j = 1; j <= allowed_types.size(); ++j ) {
		if ( types_to_allow.has_value( allowed_types[j]->base_name() ) ) continue; // already in the task's list
		///TODO figure out why this is nonconst.  Perhaps it could be const
		types_to_allow.push_back( allowed_types[j]->base_name() );
	}
	pack_task->nonconst_residue_task( ligand_residue_id ).restrict_restypes( types_to_allow );
}

utility::vector1<protocols::moves::MoverOP>
HighResDocker::create_rigid_body_movers(core::pose::Pose const & pose) const{
	utility::vector1<protocols::moves::MoverOP> rigid_body_movers;

	LigandAreas const ligand_areas =
		movemap_builder_->get_sc_interface_builder()->get_ligand_areas();

	for ( LigandAreas::value_type const & ligand_area_pair : ligand_areas ) {
		char const & chain= ligand_area_pair.first;
		utility::vector1<core::Size> jump_ids= core::pose::get_jump_ids_from_chain(chain, pose);
		for ( core::Size const jump_id : jump_ids ) {
			LigandAreaOP const ligand_area = ligand_area_pair.second;

			core::Real const & angstroms= ligand_area->high_res_angstroms_;
			core::Real const & degrees= ligand_area->high_res_degrees_;
			protocols::moves::MoverOP rigid_body_mover( new protocols::rigid::RigidBodyPerturbMover( jump_id, degrees, angstroms) );
			rigid_body_movers.push_back(rigid_body_mover);
		}

	}
	return rigid_body_movers;
}

void HighResDocker::apply_rigid_body_moves(
	core::pose::Pose & pose,
	utility::vector1<protocols::moves::MoverOP> & rigid_body_movers
){
	// utility::vector1<protocols::moves::MoverOP>::iterator rigid_body_mover= rigid_body_movers.begin(); // Unused variable causes warning.
	for ( protocols::moves::MoverOP rigid_body_mover : rigid_body_movers ) {
		rigid_body_mover->apply(pose);
	}
}

std::string HighResDocker::get_name() const {
	return mover_name();
}

std::string HighResDocker::mover_name() {
	return "HighResDocker";
}

void HighResDocker::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("scorefxn", xs_string, "Score function to be used.")
		+ XMLSchemaAttribute::required_attribute("movemap_builder", xs_string, "Name of a previously defined MoveMaoBuilder.")
		+ XMLSchemaAttribute("resfile", xs_string, "Name (path to) the resfile.")
		+ XMLSchemaAttribute("cycles", xsct_non_negative_integer, "Number of cycles to run.")
		+ XMLSchemaAttribute("repack_every_Nth", xsct_non_negative_integer, "Perform side chain repacking every Nth cycle.")
        + XMLSchemaAttribute("rmsd", xsct_real, "Maximum RMSD to be sampled away from the starting position.");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Randomly connects a fragment from the library to the growing ligand.", attlist );
}

/// @brief Provide the citation.
void
HighResDocker::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	basic::citation_manager::CitationCollectionOP cc(
		utility::pointer::make_shared< basic::citation_manager::CitationCollection >(
		"HighResDocker", basic::citation_manager::CitedModuleType::Mover
		)
	);

	cc->add_citation( basic::citation_manager::CitationManager::get_instance()->get_citation_by_doi( "10.1007/978-1-61779-465-0_10" ) );

	citations.add( cc );
}

std::string HighResDockerCreator::keyname() const {
	return HighResDocker::mover_name();
}

protocols::moves::MoverOP
HighResDockerCreator::create_mover() const {
	return utility::pointer::make_shared< HighResDocker >();
}

void HighResDockerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HighResDocker::provide_xml_schema( xsd );
}

void HighResDocker::set_allow_repacking(bool input)
{
	allow_repacking_ = input;
}

void HighResDocker::set_all_residues(bool input)
{
	all_residues_ = input;
}

void HighResDocker::set_rmsd_limit(core::Real rmsd_limit)
{
    check_rmsd_ = true;
    rmsd_max_ = rmsd_limit;
}

void HighResDocker::disable_rmsd_limit()
{
    check_rmsd_ = false;
}

} //namespace ligand_docking
} //namespace protocols