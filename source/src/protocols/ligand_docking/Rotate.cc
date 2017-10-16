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
using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer TR( "protocols.ligand_docking.ligand_options.rotate" );

// XRW TEMP std::string
// XRW TEMP RotateCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return Rotate::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RotateCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new Rotate );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP Rotate::mover_name()
// XRW TEMP {
// XRW TEMP  return "Rotate";
// XRW TEMP }

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

Rotate::Rotate(Rotate_info rotate_info): Mover("Rotate"), rotate_info_(std::move(rotate_info))
{}

Rotate::Rotate(Rotate const & that):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( that ),
	rotate_info_(that.rotate_info_)
{}

Rotate::~Rotate() = default;

protocols::moves::MoverOP Rotate::clone() const {
	return protocols::moves::MoverOP( new Rotate( *this ) );
}

protocols::moves::MoverOP Rotate::fresh_instance() const {
	return protocols::moves::MoverOP( new Rotate );
}

// XRW TEMP std::string Rotate::get_name() const{
// XRW TEMP  return "Rotate";
// XRW TEMP }

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
Rotate::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose
)
{
	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotate' mover requires 'chain' tag");
	if ( ! tag->hasOption("distribution") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotate' mover requires 'distribution' tag");
	if ( ! tag->hasOption("degrees") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotate' mover requires 'degrees' tag");
	if ( ! tag->hasOption("cycles") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotate' mover requires 'cycles' tag");

	// Will return a nullptr if this XML didn't define any ScoringGrids
	grid_set_prototype_ = protocols::qsar::scoring_grid::parse_optional_grid_set_from_tag( tag, data_map );

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
		for ( auto const & tag_along_chain_str : tag_along_chain_strs ) {
			utility::vector1<core::Size> chain_ids= get_chain_ids_from_chain(tag_along_chain_str, pose);
			for ( core::Size const chain_id : chain_ids ) {
				rotate_info_.tag_along_chains.push_back(chain_id);
				rotate_info_.tag_along_jumps.push_back( core::pose::get_jump_id_from_chain_id(chain_id, pose) );
				core::Size const chain_begin (pose.conformation().chain_begin(chain_id));
				debug_assert( chain_begin == pose.conformation().chain_end(chain_id));
				rotate_info_.tag_along_residues.push_back( chain_begin );
			}
		}
	}
}

void Rotate::apply(core::pose::Pose & pose){
	core::Vector const center = protocols::geometry::downstream_centroid_by_jump(pose, rotate_info_.jump_id);

	if ( grid_set_prototype_ == nullptr ) {
		utility::vector1<core::Size> all_chain_ids = rotate_info_.tag_along_chains;
		all_chain_ids.push_back(rotate_info_.chain_id);
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

		for ( core::conformation::ResidueCOP residue : ligands[jump_choice].residues ) {
			pose.replace_residue(chain_begin, *residue, false /*orient backbone*/);// assume rotamers are oriented?
			++chain_begin;
		}
		for ( core::Size i=1; i <= rotate_info_.tag_along_residues.size(); ++i ) {
			// I cannot figure out what this assert is testing for; it seems to be comparing a ResidueOP to a Size.
			// In any case, it is causing a comparison warning, so I am commenting it out. ~Labonte
			//debug_assert(rotate_info_.tag_along_residues.size() == ligands[jump_choice].tag_along_residues[i]);
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

	for ( core::Size const chain_id : rotate_info_.tag_along_chains ) {
		core::conformation::ResidueCOPs tag_along_residues = core::pose::get_chain_residues(local_pose, chain_id);
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
		+ XMLSchemaAttribute::required_attribute("degrees", xsct_non_negative_integer, "How degrees should be rotated around. Recommended=360")
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
	return protocols::moves::MoverOP( new Rotate );
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
	core::Real const diverse_rms = 0.65 * std::sqrt((double) heavy_atom_number);

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
