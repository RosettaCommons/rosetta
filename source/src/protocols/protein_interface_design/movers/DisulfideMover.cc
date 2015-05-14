// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/DisulfideMover.hh
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date 4/30/2009

// Unit headers
#include <protocols/protein_interface_design/movers/DisulfideMover.hh>
#include <protocols/protein_interface_design/movers/DisulfideMoverCreator.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <protocols/scoring/Interface.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
//parsing
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/util/disulfide_util.hh>

//get_resnum crap

// Utility Headers
#include <utility/vector1.hh>

// Unit Headers

// C++ headers
//#include <map>
#include <algorithm>

#include <core/scoring/disulfides/CentroidDisulfidePotential.hh>
#include <utility/vector0.hh>
#include <basic/Tracer.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace protocols::moves;
using namespace core::chemical;
using namespace std;

using utility::vector1;
using utility::vector1_int;
using core::pose::Pose;
using core::pose::PoseOP;
using core::conformation::Residue;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.DisulfideMover" );

const core::scoring::disulfides::CentroidDisulfidePotential DisulfideMover::potential_;

std::string
DisulfideMoverCreator::keyname() const
{
	return DisulfideMoverCreator::mover_name();
}

protocols::moves::MoverOP
DisulfideMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DisulfideMover );
}

std::string
DisulfideMoverCreator::mover_name()
{
	return "DisulfideMover";
}

/// @brief default ctor
DisulfideMover::DisulfideMover() :
	parent(),
	rb_jump_(1)
{}

/// @brief copy ctor
DisulfideMover::DisulfideMover(DisulfideMover const& dm) :
	//utility::pointer::ReferenceCount(),
	parent( dm ),
	rb_jump_(dm.rb_jump_)
{}

/// @brief Constructor with a single target residue
DisulfideMover::DisulfideMover( core::Size targetResidue ) :
	parent(),
	rb_jump_(1)
{
	target_residues_.push_back(targetResidue);
}

/// @brief Constructor with multiple target residues
DisulfideMover::DisulfideMover( utility::vector1<core::Size> const& targetResidues ) :
	parent(),
	rb_jump_(1)
{
	target_residues_ = targetResidues;
}

DisulfideMover::~DisulfideMover() {}

/// @brief Modify the pose to define a disulfide bond between the two specified
///   residues.
/// @details Does not do the repacking & minimization required to place the
///   disulfide correctly.
void DisulfideMover::form_disulfide(Pose & pose, Size lower_res, Size upper_res)
{
	ResidueTypeSetCOP restype_set =
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	Residue const cyd( restype_set->name_map("CYS:disulfide"), true/*dummy*/ );

	//mutate the two residues
	pose.replace_residue(lower_res, cyd, true /*orient backbone*/);
	pose.replace_residue(upper_res, cyd, true /*orient backbone*/);

	//form the bond between the two residues
	pose.conformation().declare_chemical_bond(lower_res,"SG",upper_res,"SG");
}

void DisulfideMover::apply( Pose & pose ) {
	vector1< pair<Size,Size> > potential_disulfides;
	disulfide_list(pose, target_residues_, rb_jump_, potential_disulfides);
	TR.Info << "Found " << potential_disulfides.size()
		<< " potential disulfide bonds." << endl;

	allowed_aas_[ chemical::aa_cys ] = false;
	allowed_aas_[ chemical::aa_gly ] = false;
	allowed_aas_[ chemical::aa_pro ] = false;

	PoseOP best_pose = 0;

	for( vector1< pair<Size,Size> >::const_iterator
			disulf = potential_disulfides.begin(),
			end_disulf = potential_disulfides.end();
			disulf != end_disulf;
			++disulf )
	{
		TR.Debug << "Trying disulfide from " << disulf->first
			<< " to " << disulf->second << endl;

		PoseOP trial_pose( new Pose(pose) );

		//Introduce a disulfide bond
		form_disulfide(*trial_pose, disulf->first, disulf->second);

		trial_pose->update_residue_neighbors();
		protocols::scoring::Interface interface(rb_jump_);
		interface.calculate(*trial_pose);

		// Setup Packer
		task_ = core::pack::task::TaskFactory::create_packer_task( *trial_pose );
		task_->initialize_from_command_line().or_include_current( true );
		task_->restrict_to_repacking();

		// for each residue
		for( Size i(1); i <= trial_pose->total_residue(); ++i )
		{
			Residue const& res(trial_pose->residue(i));
			if( !res.is_protein() )
				continue;
			// Repack interface residues (except for disallowed aas)
			if( interface.is_interface(i) &&
				(	allowed_aas_[ res.aa() ] || //allowed or a target disulf
					i == disulf->first ||
					i == disulf->second
				) )
				continue;

			// Other residues are unchanged
			task_->nonconst_residue_task(i).prevent_repacking();
		}

		// Setup minimizer: allow interface to move
		kinematics::MoveMapOP mm( new kinematics::MoveMap );
		mm->set_bb( false );
		mm->set_jump( rb_jump_, false );
		for( core::Size i=1; i<=trial_pose->total_residue(); ++i ) {
			Residue const& res(trial_pose->residue(i));

			if ( !trial_pose->residue(i).is_protein() ) continue;
			if( interface.is_interface(i) &&
				(	allowed_aas_[ res.aa() ] || //allowed or a target disulf
					i == disulf->first ||
					i == disulf->second ))
			{
				mm->set_chi( i, true );
;
				//mm->set_bb( i, true );
			}
		}

		// Form disulfide bond
		core::util:: rebuild_disulfide(*trial_pose,disulf->first, disulf->second,
			task_, scorefxn_repack_, mm, scorefxn_minimize_);

		std::string name = trial_pose->residue(disulf->first).name();
		assert(name == "CYS:disulfide");

		// Is this pose better than the previous best pose?
		if( !best_pose ||
			best_pose->energies().total_energy() > trial_pose->energies().total_energy() )
			// && great disulfide bond?
		{
			best_pose = trial_pose;
		}
	}

	if( best_pose )
		pose = *best_pose;

}

std::string
DisulfideMover::get_name() const {
	return DisulfideMoverCreator::mover_name();
}

/// @brief Find all residues which could disulfide bond to a target
/// @return pairs of residues (target, host) from the target protein and the
///   docking protein.
void DisulfideMover::disulfide_list( Pose const& const_pose,
		vector1< Size > const& targets, Size rb_jump,
		vector1< pair<Size,Size> > & disulfides )
{
	disulfides.clear();

	//We make a nonconst copy just so that the interface can be updated first
	//I don't have a test case where this is necessary, but other code does it too
	Pose pose(const_pose);
	pose.update_residue_neighbors();

	protocols::scoring::Interface interface(rb_jump);
	interface.calculate(pose);

	vector1<vector1_int> pairs = interface.pair_list();

	// divide the targets by interface
	vector1<Size> lower_targets;
	vector1<Size> upper_targets;
	for(vector1<Size>::const_iterator target = targets.begin(),
			end_target = targets.end();
			target != end_target; ++target)
	{
		if( find(pairs[1].begin(), pairs[1].end(), static_cast<int>(*target)) != pairs[1].end() ) {
			//lower partner
			lower_targets.push_back(*target);
		}
		else if( find(pairs[2].begin(), pairs[2].end(), static_cast<int>(*target)) != pairs[2].end() ) {
			//upper partner
			upper_targets.push_back(*target);
		}
		else { //Target not in the interface
			TR.Debug << "Target "<< *target
				<< " does not lie on the protein interface." << std::endl;
			continue; //ignore this target
		}
	}

	// If no targets were specified on one partner, just use everything in the interface
	if( lower_targets.empty() )
		lower_targets = pairs[1];
	if( upper_targets.empty() )
		upper_targets = pairs[2];

	// loop through all targets. If one is found, return it.
	for(vector1<Size>::const_iterator lower_target = lower_targets.begin(),
			end_lower_target = lower_targets.end();
			lower_target != end_lower_target; ++lower_target)
	{
		for( vector1<Size>::const_iterator upper_target = upper_targets.begin(),
				end_host = upper_targets.end();
				upper_target != end_host; ++upper_target )
		{
			if( interface.is_pair(pose.residue(*lower_target), pose.residue(*upper_target)) &&
				potential_.is_disulfide(pose.residue(*lower_target), pose.residue(*upper_target)) )
			{
				disulfides.push_back(std::make_pair(*lower_target,*upper_target));
			}
		}
	}
}

/**
 * @brief Reinitialize this Mover with parameters from the specified tags.
 * @details Parameters recognized:
 *  - target_pdb_num or target_res_num. A single target residue to form disulfides to
 *  - target_pdb_nums or target_res_nums. A list of possible target residues
 */
void DisulfideMover::parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		Movers_map const &,
		Pose const & pose)
{

	// Set target to the residue specified by "target_pdb_num" or "target_res_num"
	if( tag->hasOption("targets") ){
		target_residues_ = core::pose::get_resnum_list(tag, "targets",pose);
	}

	using namespace core::scoring;
	scorefxn_repack_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_repack", data )->clone();
	scorefxn_minimize_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_minimize", data )->clone();

	TR<<"DisulfideMover targeting residues ";
	for(vector1<Size>::const_iterator target = target_residues_.begin();
			target != target_residues_.end();++target)
	{
		TR<<*target<<", ";
	}
	TR<<std::endl;
}


} // movers
} // protein_interface_design
} // protocols
