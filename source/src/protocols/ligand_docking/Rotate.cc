// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/Rotate.hh>
#include <protocols/ligand_docking/RotateCreator.hh>

#include <protocols/ligand_docking/grid_functions.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/palette/NoDesignPackerPalette.hh>

#include <core/scoring/rms_util.tmpl.hh>

#include <protocols/qsar/scoring_grid/GridSet.hh>
#include <protocols/qsar/scoring_grid/schema_util.hh>

// Utility Headers
#include <numeric/random/random.hh>
#include <utility>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <algorithm>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

//Auto Headers

namespace protocols {
namespace ligand_docking {

static basic::Tracer TR( "protocols.ligand_docking.ligand_options.rotate" );

//////////////////////////////////////////////////////////////////////////////////////

void
Rotate_info::set_chain_id( core::Size id ) {
	by_string_ = false;
	chain_number_ = id;
}

void
Rotate_info::set_chain_letter( std::string const & str) {
	by_string_ = true;
	chain_string_ = str;
}

core::Size
Rotate_info::chain_id( core::pose::Pose const & pose ) const {
	if ( ! by_string_ ) {
		return chain_number_;
	}
	utility::vector1< core::Size > chain_ids( core::pose::get_chain_ids_from_chain(chain_string_, pose) );
	if ( chain_ids.size() == 0 ) {
		utility_exit_with_message( "'Rotate' mover: chain '"+chain_string_+"' does not exist.");
	} else if ( chain_ids.size() > 1 ) {
		utility_exit_with_message( "'Rotate' mover: chain letter '"+chain_string_+"' represents more than one chain. Consider the 'Rotates' mover (with an 's') instead.");
	}
	return chain_ids[1];
}

char
Rotate_info::chain_letter( core::pose::Pose const & pose ) const {
	if ( by_string_ ) {
		if ( chain_string_.size() != 1 ) {
			utility_exit_with_message("'Translate' mover: chain designation '"+chain_string_+"' is not one character.");
		}
		return chain_string_[1];
	}

	return core::pose::get_chain_from_chain_id( chain_number_, pose );
}

core::Size
Rotate_info::jump_id( core::pose::Pose const & pose ) const {
	// This indirection existed in the previous version of the mover - we probably want to short-circuit it.
	return core::pose::get_jump_id_from_chain_id(chain_id(pose), pose);
}

utility::vector1<core::Size>
Rotate_info::tag_along_chain_ids( core::pose::Pose const & pose ) const {
	utility::vector1<core::Size> chain_ids;
	for ( auto const & tag_along_chain_str : tag_along_chains_ ) {
		utility::vector1<core::Size> chain_ids = get_chain_ids_from_chain(tag_along_chain_str, pose);
		for ( core::Size const chain_id : chain_ids ) {
			chain_ids.push_back(chain_id);
		}
	}
	return chain_ids;
}

utility::vector1<core::Size>
Rotate_info::tag_along_jumps( core::pose::Pose const & pose ) const {
	utility::vector1<core::Size> jump_ids;
	for ( core::Size chain_id: tag_along_chain_ids( pose ) ) {
		jump_ids.push_back( core::pose::get_jump_id_from_chain_id(chain_id, pose) );
	}
	return jump_ids;
}

utility::vector1<core::Size>
Rotate_info::tag_along_residues( core::pose::Pose const & pose ) const {
	utility::vector1<core::Size> residues;
	for ( core::Size chain_id: tag_along_chain_ids( pose ) ) {
		core::Size const chain_begin (pose.conformation().chain_begin(chain_id));
		debug_assert( chain_begin == pose.conformation().chain_end(chain_id) );
		residues.push_back( chain_begin );
	}
	return residues;
}

void
Rotate_info::set_tag_along_chains( utility::vector1<std::string> const & setting ) {
	tag_along_chains_ = setting;
}

//////////////////////////////////////////////////////////////////////////////////////

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

Rotate::Rotate(Rotate_info const & rotate_info): Mover("Rotate"), rotate_info_(rotate_info)
{}

protocols::moves::MoverOP Rotate::clone() const {
	return utility::pointer::make_shared< Rotate >( *this );
}

protocols::moves::MoverOP Rotate::fresh_instance() const {
	return utility::pointer::make_shared< Rotate >();
}


/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
Rotate::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map
)
{
	if ( ! tag->hasOption("chain") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Rotate' mover requires 'chain' tag");
	if ( ! tag->hasOption("distribution") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Rotate' mover requires 'distribution' tag");
	if ( ! tag->hasOption("degrees") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Rotate' mover requires 'degrees' tag");
	if ( ! tag->hasOption("cycles") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Rotate' mover requires 'cycles' tag");

	// Will return a nullptr if this XML didn't define any ScoringGrids
	grid_set_prototype_ = protocols::qsar::scoring_grid::parse_optional_grid_set_from_tag( tag, data_map );

	rotate_info_.set_chain_letter( tag->getOption<std::string>("chain") );

	std::string distribution_str= tag->getOption<std::string>("distribution");
	rotate_info_.distribution= get_distribution(distribution_str);
	rotate_info_.degrees = tag->getOption<core::Size>("degrees");
	rotate_info_.cycles = tag->getOption<core::Size>("cycles");

	runtime_assert_string_msg( rotate_info_.degrees >= 0.0, "Error in Rotate::parse_my_tag(): The number of degrees must be greater than or equal to zero.  Got " + std::to_string(rotate_info_.degrees) + ", which is invalid." );

	if ( tag->hasOption("tag_along_chains") ) {
		std::string const tag_along_chains_str = tag->getOption<std::string>("tag_along_chains");
		rotate_info_.set_tag_along_chains( utility::string_split(tag_along_chains_str, ',') );
	}
}

void Rotate::apply(core::pose::Pose & pose){
	core::Size chain_id = rotate_info_.chain_id( pose );
	core::Size jump_id = rotate_info_.jump_id( pose );

	core::Vector const center = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);

	if ( grid_set_prototype_ == nullptr ) {
		utility::vector1<core::Size> all_chain_ids = rotate_info_.tag_along_chain_ids( pose );
		all_chain_ids.push_back(chain_id);
		utility::pointer::shared_ptr<core::grid::CartGrid<int> > const grid = make_atr_rep_grid_without_ligands(pose, center, all_chain_ids);
		rotate_ligand(grid, pose);// move ligand to a random point in binding pocket
	} else {
		TR.Error << "############################################################" << std::endl << std::endl;
		TR.Error << "The Rotate mover will not work properly with defined GridSets - it will do nothing!" << std::endl;
		TR.Error << "You probably want to switch over to the Transform Mover, anyway." << std::endl;
		TR.Error << std::endl;
		TR.Error << "############################################################" << std::endl;
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

void
Rotate::set_chain( std::string const & chain ) {
	rotate_info_.set_chain_letter( chain );
}

void
Rotate::set_chain_id( core::Size chain_id ) {
	rotate_info_.set_chain_id( chain_id );
}

/// @brief for n random rotations, randomly pick one from among the best scoring set of diverse poses
void Rotate::rotate_ligand(
	utility::pointer::shared_ptr<core::grid::CartGrid<int> >  const & grid,
	core::pose::Pose & pose
) {
	if ( rotate_info_.degrees == 0 ) return;

	core::Size chain_id = rotate_info_.chain_id( pose );
	core::Size jump_id = rotate_info_.jump_id( pose );

	protocols::rigid::RigidBodyMoverOP mover;
	if ( rotate_info_.distribution == Uniform ) {
		mover = utility::pointer::make_shared< protocols::rigid::RigidBodyRandomizeMover >( pose, jump_id, protocols::rigid::partner_downstream, rotate_info_.degrees, rotate_info_.degrees);
	} else if ( rotate_info_.distribution == Gaussian ) {
		mover = utility::pointer::make_shared< protocols::rigid::RigidBodyPerturbMover > ( jump_id, rotate_info_.degrees, 0 /*translate*/);
	}

	core::Size chain_begin = pose.conformation().chain_begin(chain_id);
	utility::vector1< Ligand_info> ligands= create_random_rotations(grid, mover, pose, chain_begin);

	auto const jump_choice=  (core::Size) numeric::random::rg().random_range(1, ligands.size());
	{
		pose.set_jump(jump_id, ligands[jump_choice].jump);

		for ( core::conformation::ResidueCOP residue : ligands[jump_choice].residues ) {
			pose.replace_residue(chain_begin, *residue, false /*orient backbone*/);// assume rotamers are oriented?
			++chain_begin;
		}
		utility::vector1< core::Size > const & tag_along_residues = rotate_info_.tag_along_residues( pose );
		for ( core::Size i=1; i <= tag_along_residues.size(); ++i ) {
			debug_assert( i <= ligands[jump_choice].tag_along_residues.size() );
			core::Size residue_id = tag_along_residues[i];
			core::conformation::ResidueCOP residue = ligands[jump_choice].tag_along_residues[i];
			pose.replace_residue(residue_id, *residue, false /*orient backbone*/);// assume rotamers are oriented?
		}
	}
}

void Rotate::rotate_ligand(core::pose::Pose & pose)
{
	if ( rotate_info_.degrees == 0 ) return;

	core::Size jump_id = rotate_info_.jump_id( pose );

	protocols::rigid::RigidBodyMoverOP mover;
	if ( rotate_info_.distribution == Uniform ) {
		mover = utility::pointer::make_shared< protocols::rigid::RigidBodyRandomizeMover >( pose, jump_id, protocols::rigid::partner_downstream, rotate_info_.degrees, rotate_info_.degrees);
	} else if ( rotate_info_.distribution == Gaussian ) {
		mover = utility::pointer::make_shared< protocols::rigid::RigidBodyPerturbMover > ( jump_id, rotate_info_.degrees, 0 /*translate*/);
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
	core::Size chain_id = rotate_info_.chain_id( pose );
	core::Size jump_id = rotate_info_.jump_id( pose );

	core::Size const end = pose.conformation().chain_end(chain_id);
	core::Size const heavy_atom_number= core::pose::num_heavy_atoms(begin, end, pose);
	core::pose::Pose local_pose= pose;
	local_pose.remove_constraints();
	core::Vector const center = protocols::geometry::downstream_centroid_by_jump(local_pose, jump_id);

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
	core::Size chain_id = rotate_info_.chain_id( local_pose );
	core::Size jump_id = rotate_info_.jump_id( local_pose );
	utility::vector1< core::Size > tag_along_jumps = rotate_info_.tag_along_jumps( local_pose );

	apply_rotate(mover, local_pose, center, jump_id, tag_along_jumps);
	core::pack::task::PackerTaskCOP packertask( core::pack::task::TaskFactory::create_packer_task( local_pose, utility::pointer::make_shared< core::pack::palette::NoDesignPackerPalette >() ) );
	rb_grid_rotamer_trials_atr_rep(*grid, local_pose, begin, end, *packertask);
	core::kinematics::Jump jump= local_pose.jump(jump_id);
	std::pair<int, int> const scores= get_rb_atr_and_rep_scores(*grid, local_pose, begin, end);
	Ligand_info ligand_info;
	ligand_info.jump= jump;
	ligand_info.atr= scores.first;
	ligand_info.rep= scores.second;

	ligand_info.residues = core::pose::get_chain_residues(local_pose, chain_id);

	for ( core::Size const ta_chain_id : rotate_info_.tag_along_chain_ids( local_pose ) ) {
		core::conformation::ResidueCOPs tag_along_residues = core::pose::get_chain_residues(local_pose, ta_chain_id);
		debug_assert(tag_along_residues.size() == 1);
		ligand_info.tag_along_residues.push_back(tag_along_residues[1]);
	}
	return ligand_info;
}

std::string Rotate::get_name() const {
	return mover_name();
}

std::string Rotate::mover_name() {
	return "Rotate";
}

void Rotate::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	// Restriction for the distribution attribute
	XMLSchemaRestriction restriction_type;
	restriction_type.name( "distribution_string" );
	restriction_type.base_type( xs_string );
	restriction_type.add_restriction( xsr_pattern, "uniform|gaussian" );
	xsd.add_top_level_element( restriction_type );
	// Attributes
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("distribution", xs_string, "Sampling distribution; Either \"uniform\" or \"gaussian\"")
		+ XMLSchemaAttribute::required_attribute("degrees", xsct_real, "How degrees should be rotated around. Recommended=360")
		+ XMLSchemaAttribute::required_attribute("cycles", xsct_non_negative_integer, "Number of cycles. Recommended: 1000")
		+ XMLSchemaAttribute("chain", xs_string, "Chain ID. MUST be a completely connected single chain.")
		+ XMLSchemaAttribute("tag_along_chains", xs_string,
		"Comma separated list of chains to be moved together with the rotating chain. E.g. metals or water.");

	protocols::qsar::scoring_grid::attributes_for_parse_optional_grid_set_from_tag( attlist, "The ScoringGrid set to use for scoring the rotation. "
		"If no scoring grids (at all) are present in the XML, use a default classic grid." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Perform a course random rotation "
		" throughout all rotational degrees of freedom", attlist );
}

std::string RotateCreator::keyname() const {
	return Rotate::mover_name();
}

protocols::moves::MoverOP
RotateCreator::create_mover() const {
	return utility::pointer::make_shared< Rotate >();
}

void RotateCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	Rotate::provide_xml_schema( xsd );
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

	for ( core::Size const jump_id : tag_along_jumps ) {
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
	debug_assert(heavy_atom_number > 0);

	// This next parameter is a wild heuristic guesses that seem OK for the Meiler x-dock set.
	core::Real const diverse_rms = 0.65 * std::sqrt((core::Real) heavy_atom_number);

	core::conformation::ResidueCOPs const & these_residues= ligand.get_residues();

	for ( Ligand_info const & ligand_info : ligands ) { // if ligands is empty we still return true so no need to check for this condition.
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
