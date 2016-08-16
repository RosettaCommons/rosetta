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
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/Rotate.hh>
#include <protocols/ligand_docking/RotateCreator.hh>

#include <protocols/ligand_docking/grid_functions.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/rms_util.tmpl.hh>

#include <protocols/qsar/scoring_grid/GridManager.hh>

// Utility Headers
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <algorithm>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

//Auto Headers
using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer rotate_tracer( "protocols.ligand_docking.ligand_options.rotate", basic::t_debug );

std::string
RotateCreator::keyname() const
{
	return RotateCreator::mover_name();
}

protocols::moves::MoverOP
RotateCreator::create_mover() const {
	return protocols::moves::MoverOP( new Rotate );
}

std::string
RotateCreator::mover_name()
{
	return "Rotate";
}

Ligand_info::Ligand_info():residues(), atr(0), rep(0), jump(){}

Ligand_info::Ligand_info(core::conformation::ResidueCOPs const & residues, int atr, int rep):
	residues(residues), atr(atr), rep(rep), jump(){}

Ligand_info::Ligand_info(core::conformation::ResidueCOPs const & residues, std::pair<int,int> scores, core::kinematics::Jump jump):
	residues(residues), atr(scores.first), rep(scores.second), jump(jump){}


bool Ligand_info::operator<(Ligand_info const & ligand_info) const{
	return ( rep < ligand_info.rep || (rep == ligand_info.rep && atr < ligand_info.atr ) );
}
bool Ligand_info::operator<(std::pair<int,int> const & scores) const{
	return rep < scores.second || (rep == scores.second && atr < scores.first);
}
core::conformation::ResidueCOPs
const & Ligand_info::get_residues() const{
	return residues;
}

/// @brief
Rotate::Rotate(): Mover("Rotate")
{}

Rotate::Rotate(Rotate_info rotate_info): Mover("Rotate"), rotate_info_(rotate_info)
{}

Rotate::Rotate(Rotate const & that):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( that ),
	rotate_info_(that.rotate_info_)
{}

Rotate::~Rotate() {}

protocols::moves::MoverOP Rotate::clone() const {
	return protocols::moves::MoverOP( new Rotate( *this ) );
}

protocols::moves::MoverOP Rotate::fresh_instance() const {
	return protocols::moves::MoverOP( new Rotate );
}

std::string Rotate::get_name() const{
	return "Rotate";
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
Rotate::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data_map*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose
)
{
	if ( tag->getName() != "Rotate" ) {
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotate' mover requires 'chain' tag");
	if ( ! tag->hasOption("distribution") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotate' mover requires 'distribution' tag");
	if ( ! tag->hasOption("degrees") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotate' mover requires 'degrees' tag");
	if ( ! tag->hasOption("cycles") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotate' mover requires 'cycles' tag");

	rotate_info_.chain = tag->getOption<std::string>("chain");
	utility::vector1< core::Size > chain_ids( core::pose::get_chain_ids_from_chain(rotate_info_.chain, pose) );
	if ( chain_ids.size() == 0 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("'Rotate' mover: chain '"+rotate_info_.chain+"' does not exist.");
	} else if ( chain_ids.size() > 1 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("'Rotate' mover: chain letter '"+rotate_info_.chain+"' represents more than one chain. Consider using the 'Rotates' mover (with an 's') instead.");
	}
	rotate_info_.chain_id= chain_ids[1];
	rotate_info_.jump_id= core::pose::get_jump_id_from_chain_id(rotate_info_.chain_id, pose);
	std::string distribution_str= tag->getOption<std::string>("distribution");
	rotate_info_.distribution= get_distribution(distribution_str);
	rotate_info_.degrees = tag->getOption<core::Size>("degrees");
	rotate_info_.cycles = tag->getOption<core::Size>("cycles");

	if ( tag->hasOption("tag_along_chains") ) {
		std::string const tag_along_chains_str = tag->getOption<std::string>("tag_along_chains");
		utility::vector1<std::string> tag_along_chain_strs= utility::string_split(tag_along_chains_str, ',');
		BOOST_FOREACH ( std::string tag_along_chain_str, tag_along_chain_strs ) {
			utility::vector1<core::Size> chain_ids= get_chain_ids_from_chain(tag_along_chain_str, pose);
			BOOST_FOREACH ( core::Size chain_id, chain_ids ) {
				rotate_info_.tag_along_chains.push_back(chain_id);
				rotate_info_.tag_along_jumps.push_back( core::pose::get_jump_id_from_chain_id(chain_id, pose) );
				core::Size const chain_begin (pose.conformation().chain_begin(chain_id));
				assert( chain_begin == pose.conformation().chain_end(chain_id));
				rotate_info_.tag_along_residues.push_back( chain_begin );
			}
		}
	}
}

void Rotate::apply(core::pose::Pose & pose){
	core::Vector const center = protocols::geometry::downstream_centroid_by_jump(pose, rotate_info_.jump_id);

	qsar::scoring_grid::GridManager* grid_manager = qsar::scoring_grid::GridManager::get_instance();
	if ( grid_manager->size() == 0 ) {
		utility::vector1<core::Size> all_chain_ids = rotate_info_.tag_along_chains;
		all_chain_ids.push_back(rotate_info_.chain_id);
		utility::pointer::shared_ptr<core::grid::CartGrid<int> > const grid = make_atr_rep_grid_without_ligands(pose, center, all_chain_ids);
		rotate_ligand(grid, pose);// move ligand to a random point in binding pocket
	} else {
		// TODO refactor qsar map so it works properly
		/*
		if(grid_manager->is_qsar_map_attached())
		{
		//core::conformation::ResidueOP residue = new core::conformation::Residue(pose.residue(begin));
		//qsar::qsarMapOP qsar_map(new qsar::qsarMap("default",residue));
		//qsar_map->fill_with_value(1);
		//grid_manager->set_qsar_map(qsar_map);
		}
		*/
		//grid_manager->initialize_all_grids(center);
		//grid_manager->update_grids(pose,center);
	}
}

/// @brief for n random rotations, randomly pick one from among the best scoring set of diverse poses
void Rotate::rotate_ligand(
	utility::pointer::shared_ptr<core::grid::CartGrid<int> >  const & grid,
	core::pose::Pose & pose
) {
	if ( rotate_info_.degrees == 0 ) return;

	protocols::rigid::RigidBodyMoverOP mover;
	if ( rotate_info_.distribution == Uniform ) {
		mover = protocols::rigid::RigidBodyMoverOP( new protocols::rigid::RigidBodyRandomizeMover( pose, rotate_info_.jump_id, protocols::rigid::partner_downstream, rotate_info_.degrees, rotate_info_.degrees) );
	} else if ( rotate_info_.distribution == Gaussian ) {
		mover = protocols::rigid::RigidBodyMoverOP( new protocols::rigid::RigidBodyPerturbMover ( rotate_info_.jump_id, rotate_info_.degrees, 0 /*translate*/) );
	}

	core::Size chain_begin = pose.conformation().chain_begin(rotate_info_.chain_id);
	utility::vector1< Ligand_info> ligands= create_random_rotations(grid, mover, pose, chain_begin);

	core::Size const jump_choice=  (core::Size) numeric::random::rg().random_range(1, ligands.size());
	{
		pose.set_jump(rotate_info_.jump_id, ligands[jump_choice].jump);

		BOOST_FOREACH ( core::conformation::ResidueCOP residue, ligands[jump_choice].residues ) {
			pose.replace_residue(chain_begin, *residue, false /*orient backbone*/);// assume rotamers are oriented?
			++chain_begin;
		}
		for ( core::Size i=1; i <= rotate_info_.tag_along_residues.size(); ++i ) {
			// I cannot figure out what this assert is testing for; it seems to be comparing a ResidueOP to a Size.
			// In any case, it is causing a comparison warning, so I am commenting it out. ~Labonte
			//assert(rotate_info_.tag_along_residues.size() == ligands[jump_choice].tag_along_residues[i]);
			core::Size residue_id = rotate_info_.tag_along_residues[i];
			core::conformation::ResidueCOP residue = ligands[jump_choice].tag_along_residues[i];
			pose.replace_residue(residue_id, *residue, false /*orient backbone*/);// assume rotamers are oriented?
		}
	}
}

void Rotate::rotate_ligand(core::pose::Pose & pose)
{
	if ( rotate_info_.degrees == 0 ) return;

	protocols::rigid::RigidBodyMoverOP mover;
	if ( rotate_info_.distribution == Uniform ) {
		mover = protocols::rigid::RigidBodyMoverOP( new protocols::rigid::RigidBodyRandomizeMover( pose, rotate_info_.jump_id, protocols::rigid::partner_downstream, rotate_info_.degrees, rotate_info_.degrees) );
	} else if ( rotate_info_.distribution == Gaussian ) {
		mover = protocols::rigid::RigidBodyMoverOP( new protocols::rigid::RigidBodyPerturbMover ( rotate_info_.jump_id, rotate_info_.degrees, 0 /*translate*/) );
	}
	//core::Size chain_begin = pose.conformation().chain_begin(rotate_info_.chain_id);

}

utility::vector1< Ligand_info>
Rotate::create_random_rotations(
	utility::pointer::shared_ptr<core::grid::CartGrid<int> > const & grid,
	protocols::rigid::RigidBodyMoverOP const mover,
	core::pose::Pose & pose,
	core::Size begin
)const{
	core::Size const end = pose.conformation().chain_end(rotate_info_.chain_id);
	core::Size const heavy_atom_number= core::pose::num_heavy_atoms(begin, end, pose);
	core::pose::Pose local_pose= pose;
	local_pose.remove_constraints();
	core::Vector const center = protocols::geometry::downstream_centroid_by_jump(local_pose, rotate_info_.jump_id);

	utility::vector1< Ligand_info> ligands;  ///TODO make this a set.  The set should check for another pose with a similar RMSD.
	// "num_chi_angles" code comes from Ian Davis, who knows why. I added the max fxn so that waters can rotate (they have too few chi angles)
	core::Size const max_diversity= std::max(5, static_cast<int>(5*core::pose::num_chi_angles(begin, end, local_pose)+1) );

	Ligand_info best=create_random_rotation(grid, mover, center, begin, end, local_pose);// first case;
	add_ligand_conditionally(best, ligands, heavy_atom_number);
	for ( core::Size i=1; i<= rotate_info_.cycles && ligands.size() <= max_diversity ; ++i ) {
		Ligand_info current =create_random_rotation(grid, mover, center, begin, end, local_pose);
		if ( current < best ) {
			best= current;
		}
		add_ligand_conditionally(current, ligands, heavy_atom_number);
	}
	if ( ligands.empty() ) {
		ligands.push_back(best);
	}
	return ligands;
}

Ligand_info Rotate::create_random_rotation(
	utility::pointer::shared_ptr<core::grid::CartGrid<int> > const & grid,
	protocols::rigid::RigidBodyMoverOP const mover,
	core::Vector const & center,
	core::Size const begin,
	core::Size const end,
	core::pose::Pose & local_pose
) const{
	apply_rotate(mover, local_pose, center, rotate_info_.jump_id, rotate_info_.tag_along_jumps);
	rb_grid_rotamer_trials_atr_rep(*grid, local_pose, begin, end);
	core::kinematics::Jump jump= local_pose.jump(rotate_info_.jump_id);
	std::pair<int, int> const scores= get_rb_atr_and_rep_scores(*grid, local_pose, begin, end);
	Ligand_info ligand_info;
	ligand_info.jump= jump;
	ligand_info.atr= scores.first;
	ligand_info.rep= scores.second;

	ligand_info.residues = core::pose::get_chain_residues(local_pose, rotate_info_.chain_id);

	BOOST_FOREACH ( core::Size chain_id, rotate_info_.tag_along_chains ) {
		core::conformation::ResidueCOPs tag_along_residues = core::pose::get_chain_residues(local_pose, chain_id);
		assert(tag_along_residues.size() == 1);
		ligand_info.tag_along_residues.push_back(tag_along_residues[1]);
	}
	return ligand_info;
}

void add_ligand_conditionally(
	Ligand_info const & ligand_info,
	utility::vector1< Ligand_info> & ligands,
	core::Size const heavy_atom_number
){
	if (
			check_score(ligand_info, heavy_atom_number)
			&& check_RMSD(ligand_info, heavy_atom_number, ligands)
			) {
		ligands.push_back(ligand_info);
	}
}

void apply_rotate(
	protocols::rigid::RigidBodyMoverOP mover,
	core::pose::Pose & pose,
	core::Vector const & center,
	core::Size jump_id,
	utility::vector1<core::Size> const tag_along_jumps
){
	mover->rb_jump(jump_id);
	mover->apply(pose);
	pose.update_actcoords();///TODO Verify necessity
	mover->rot_center(center); // restore the old center so ligand doesn't walk away (really necessary?)

	mover->freeze();

	BOOST_FOREACH ( core::Size jump_id, tag_along_jumps ) {
		mover->rb_jump(jump_id);
		mover->apply(pose);
	}

	mover->unfreeze();

}

bool check_score(
	Ligand_info const & ligand,
	core::Size const heavy_atom_number
){
	int const rep_threshold=0;
	int const atr_threshold=-(int) (0.85 * heavy_atom_number);
	return ligand.atr <= atr_threshold && ligand.rep <= rep_threshold;
}

bool check_RMSD(
	Ligand_info const & ligand,
	core::Size const heavy_atom_number,
	utility::vector1< Ligand_info> const & ligands
){
	assert(heavy_atom_number > 0);

	// This next parameter is a wild heuristic guesses that seem OK for the Meiler x-dock set.
	core::Real const diverse_rms = 0.65 * std::sqrt((double) heavy_atom_number);

	core::conformation::ResidueCOPs const & these_residues= ligand.get_residues();

	BOOST_FOREACH ( Ligand_info ligand_info, ligands ) { // if ligands is empty we still return true so no need to check for this condition.
		core::conformation::ResidueCOPs const & compare_residues= ligand_info.get_residues();
		runtime_assert(these_residues.size() == compare_residues.size());

		core::Real const rms = (compare_residues.size() == 1) ///TODO write multi_residue automorphic fxn.
			? core::scoring::automorphic_rmsd(*these_residues[1], *compare_residues[1], false)
			:core::scoring::rmsd_no_super(these_residues, compare_residues, core::scoring::is_ligand_heavyatom_residues);

		if ( rms < diverse_rms ) return false;
	}
	return true;
}

} //namespace ligand_docking
} //namespace protocols
