// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/ligand_docking/Transform.cc
/// @author Thomas Willcock and Darwin Fu
/// Adapted from code by Sam Deluca

// Testing

#include <protocols/ligand_docking/TransformEnsembleCreator.hh>
#include <protocols/ligand_docking/TransformEnsemble.hh>

#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>
#include <protocols/qsar/scoring_grid/schema_util.hh>
#include <protocols/ligand_docking/grid_functions.hh>
#include <protocols/jd2/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/ligand_docking/ligand_scores.hh>

#include <basic/Tracer.hh>

#include <sstream>
#include <utility/string_util.hh>

#include <numeric/numeric.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_xyz.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/map_util.hh>
#include <vector>

#include <ObjexxFCL/format.hh>
#include <utility/io/ozstream.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace ligand_docking {

static basic::Tracer transform_tracer("protocols.ligand_docking.TransformEnsemble", basic::t_debug);

TransformEnsemble::TransformEnsemble():
	Mover("TransformEnsemble")
	// & in declaration values
{}

TransformEnsemble::TransformEnsemble(TransformEnsemble const & ) = default;

TransformEnsemble::TransformEnsemble(
	protocols::qsar::scoring_grid::GridSetCOP grid_set_prototype,
	utility::vector1<std::string> const & chains,
	core::Real const & box_size,
	core::Real const & move_distance,
	core::Real const & angle,
	core::Size const & cycles,
	core::Real const & temp
) :
	Mover("TransformEnsemble"),
	grid_set_prototype_( grid_set_prototype )
	// & in declaration default values
{
	transform_info_.chains = chains;
	transform_info_.box_size = box_size;
	transform_info_.move_distance = move_distance;
	transform_info_.angle = angle;
	transform_info_.cycles = cycles;
	transform_info_.temperature = temp;
}

TransformEnsemble::~TransformEnsemble() = default;

protocols::moves::MoverOP TransformEnsemble::clone() const
{
	return protocols::moves::MoverOP( new TransformEnsemble (*this) );
}

protocols::moves::MoverOP TransformEnsemble::fresh_instance() const
{
	return protocols::moves::MoverOP( new TransformEnsemble);
}

void TransformEnsemble::parse_my_tag
(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose/*pose*/
)
{
	if ( tag->getName() != "TransformEnsemble" ) {
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("chains") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires chains tag");
	if ( ! tag->hasOption("move_distance") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires move_distance tag");
	if ( ! tag->hasOption("box_size") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires box_size tag");
	if ( ! tag->hasOption("angle") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires angle tag");
	if ( ! tag->hasOption("cycles") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires cycles tag");
	if ( !tag->hasOption("temperature") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires temperature tag");

	//Divides by root(3) so the center can only move a total equal to move_distance in each step
	transform_info_.move_distance = (tag->getOption<core::Real>("move_distance")) /sqrt(3);
	transform_info_.box_size = tag->getOption<core::Real>("box_size");
	transform_info_.angle = tag->getOption<core::Real>("angle");
	transform_info_.cycles = tag->getOption<core::Size>("cycles");
	transform_info_.temperature = tag->getOption<core::Real>("temperature");
	transform_info_.repeats = tag->getOption<core::Size>("repeats",1);
	optimize_until_score_is_negative_ = tag->getOption<bool>("optimize_until_score_is_negative",false);

	initial_perturb_ = (tag->getOption<core::Real>("initial_perturb",0.0));

	use_conformers_ = tag->getOption<bool>("use_conformers",true);
	optimize_until_ideal_ = tag->getOption<bool>("optimize_until_ideal",false);

	std::string const all_chains_str = tag->getOption<std::string>("chains");
	transform_info_.chains = utility::string_split(all_chains_str, ',');

	for ( core::Size i=1; i <= transform_info_.chains.size(); ++i ) {
		core::Size current_chain_id(core::pose::get_chain_id_from_chain(transform_info_.chains[i], pose));
		transform_info_.chain_ids.push_back(current_chain_id);
		transform_info_.jump_ids.push_back(core::pose::get_jump_id_from_chain_id(current_chain_id, pose));
	}

	if ( tag->hasOption("sampled_space_file") ) {
		output_sampled_space_ = true;
		sampled_space_file_ = tag->getOption<std::string>("sampled_space_file");
	}

	grid_set_prototype_ = protocols::qsar::scoring_grid::parse_grid_set_from_tag(tag, data_map);
}

void TransformEnsemble::apply(core::pose::Pose & pose)
{
	debug_assert( grid_set_prototype_ != nullptr );

	//Grid setup: Use centroid of all chains as center of grid
	core::Vector const center(protocols::geometry::centroid_by_chains(pose, transform_info_.chain_ids));

	// TODO: TransformEnsemble should be able to get the chains it's concerned about
	// itself, rather than relying on the GridSet.
	char chain = grid_set_prototype_->chain();
	// We pass exclude=false, so that only the passed chain is the one used in grid caching calculations
	qsar::scoring_grid::GridSetCOP grid_set( qsar::scoring_grid::GridManager::get_instance()->get_grids( *grid_set_prototype_, pose, center, chain, false ) );
	debug_assert(grid_set != nullptr); //something has gone hopelessly wrong if this triggers

	//Setting up ligands and ligand center
	utility::vector1<core::conformation::ResidueOP> single_conformers;
	core::Vector original_center(0,0,0);

	for ( core::Size i=1; i <= transform_info_.chains.size(); ++i ) {
		core::Size const begin(pose.conformation().chain_begin(transform_info_.chain_ids[i]));

		core::conformation::Residue original_residue = pose.residue(begin);
		original_center = original_center + original_residue.xyz(original_residue.nbr_atom());

		single_conformers.clear();
		rotamers_for_trials(pose,begin,single_conformers);
		ligand_conformers_.insert(std::pair<core::Size, utility::vector1< core::conformation::ResidueOP > >(i, single_conformers));

	}

	original_center = original_center/(transform_info_.chains.size());

	//Setup scores and move count

	core::Real last_score(10000.0);
	core::Real best_score(10000.0);
	core::Real current_score(10000.0);

	core::pose::Pose best_pose(pose);
	core::pose::Pose starting_pose(pose);
	core::Size accepted_moves = 0;
	core::Size rejected_moves = 0;


	for ( core::Size i=1; i <= transform_info_.chains.size(); ++i ) {
		transform_tracer << "Considering " << ligand_conformers_[i].size() << " conformers during sampling for chain " << i << std::endl;
	}

	utility::io::ozstream sampled_space;
	if ( output_sampled_space_ ) {
		sampled_space.open(sampled_space_file_);
	}


	for ( core::Size repeat = 1; repeat <= transform_info_.repeats; ++repeat ) {
		pose = starting_pose;
		core::Size cycle = 0;
		bool not_converged = true;
		best_ligands_.clear();
		ligand_residues_.clear();
		reference_residues_.clear();
		last_accepted_ligand_residues_.clear();
		last_accepted_reference_residues_.clear();
		for ( core::Size i=1; i <= transform_info_.chains.size(); ++i ) {
			core::Size begin(pose.conformation().chain_begin(transform_info_.chain_ids[i]));
			core::conformation::UltraLightResidue ligand_residue( pose.residue(begin).get_self_ptr() );
			ligand_residues_.push_back(ligand_residue);
			reference_residues_.push_back(ligand_residue);
		}

		last_accepted_ligand_residues_ = ligand_residues_;
		last_accepted_reference_residues_ = reference_residues_;

		//For benchmarking purposes it is sometimes desirable to translate/rotate the ligand
		//away from the starting point before beginning a trajectory.
		if ( initial_perturb_ > 0.0 ) {
			bool perturbed = false;
			while ( !perturbed )
					{
				translate_ligand(ligand_residues_,reference_residues_,initial_perturb_);

				//Also randomize starting conformer
				if ( ligand_conformers_.size() > 1 && use_conformers_ == true ) {

					for ( core::Size i=1; i <= ligand_residues_.size(); ++i ) {
						transform_tracer << "Doing a conformer change for ligand number: " << i << std::endl;
						change_conformer(ligand_residues_[i], reference_residues_[i], i);
					}
				}
				perturbed = true;

				//new_center calculated as the weighted atom count average of the center
				core::Vector new_center(0,0,0);
				new_center = weighted_center(ligand_residues_);
				core::Real distance = new_center.distance(original_center);

				//Not everything in grid, also checks center distance
				if ( !check_grid(grid_set, ligand_residues_, distance) ) {
					ligand_residues_ = last_accepted_ligand_residues_;
					reference_residues_ = last_accepted_reference_residues_;
					perturbed = false;
				}
			}

			last_accepted_ligand_residues_ = ligand_residues_;
			last_accepted_reference_residues_ = reference_residues_;
		}


		last_score = grid_set->average_score(ligand_residues_);
		//set post-perturb ligand as best model and score
		best_ligands_ = ligand_residues_;
		best_score = last_score;


		//Optimize until within 5 percent of the theoretical maximum score
		//Based on all atoms having an attractive score of 1
		core::Real ideal_limit = grid_set->ideal_score(ligand_residues_);
		ideal_limit = (core::Real)0.95 * ideal_limit;

		//Limit to 4 * number of cycles to prevent infinite loop
		core::Size max_cycles = (core::Size)4 * transform_info_.cycles;

		while ( not_converged )
				{
			if ( cycle >= max_cycles-1 ) {
				not_converged=false;
			} else if ( optimize_until_ideal_ ) {
				if ( cycle >= transform_info_.cycles-1 && last_score <= ideal_limit ) {
					not_converged= false;
				}
			} else if ( optimize_until_score_is_negative_ ) {

				if ( cycle >= transform_info_.cycles-1 && last_score <= 0.0 ) {
					not_converged= false;
				}

			} else {
				if ( cycle >= transform_info_.cycles-1 ) {
					not_converged= false;
				}

			}

			//Incrementer
			cycle++;
			transform_tracer << "Cycle Number " << cycle << std::endl;
			bool move_accepted = false;

			//during each move either move the ligand or try a new conformer (if there is more than one conformer)
			//Consider each conformer change scores individually but consider moves as the ensemble score.
			if ( ligand_conformers_.size() > 1 && use_conformers_ == true && numeric::random::uniform() >= 0.5 ) {

				for ( core::Size i=1; i <= ligand_residues_.size(); ++i ) {
					transform_tracer << "Doing a conformer change for ligand number: " << i << std::endl;

					change_conformer(ligand_residues_[i], reference_residues_[i], i);

					//Not everything in grid
					if ( !check_grid(grid_set, ligand_residues_) ) {
						move_accepted = false;
						ligand_residues_[i] = last_accepted_ligand_residues_[i];
						continue;
					}

					current_score = grid_set->average_score(ligand_residues_);

					//If monte_carlo rejected
					if ( !monte_carlo(current_score, last_score) ) {
						move_accepted = false;
						ligand_residues_[i]=last_accepted_ligand_residues_[i];
					} else {
						//If monte carlo accepted, check vs. best score
						move_accepted = true;
						last_accepted_ligand_residues_[i] = ligand_residues_[i];
					}

				}

				current_score = grid_set->average_score(ligand_residues_);
				if ( current_score < best_score ) {
					best_score = current_score;
					best_ligands_ = ligand_residues_;
				}

			} else {
				transform_ligand(ligand_residues_,reference_residues_);
				//new_center calculated as the weighted atom count average of the center
				core::Vector new_center(0,0,0);
				new_center = weighted_center(ligand_residues_);
				core::Real distance = new_center.distance(original_center);

				//Not everything in grid, also checks center distance
				if ( !check_grid(grid_set, ligand_residues_, distance) ) {
					ligand_residues_ = last_accepted_ligand_residues_;
					reference_residues_ = last_accepted_reference_residues_;
					move_accepted = false;
				} else {

					current_score = grid_set->average_score(ligand_residues_);

					//If monte_carlo rejected
					if ( !monte_carlo(current_score, last_score) ) {
						ligand_residues_ = last_accepted_ligand_residues_;
						reference_residues_ = last_accepted_reference_residues_;
						move_accepted = false;
					} else {
						//If monte carlo accepted, check vs. best score
						last_accepted_ligand_residues_ = ligand_residues_;
						last_accepted_reference_residues_ = reference_residues_;
						move_accepted = true;
						if ( current_score < best_score ) {
							best_score = current_score;
							best_ligands_ = ligand_residues_;
							transform_tracer << "accepting new best pose" << std::endl;
						}
					}
				}

			}

			if ( move_accepted ) {
				accepted_moves++;
			} else {
				rejected_moves++;
			}

			if ( output_sampled_space_ ) {
				for ( core::conformation::UltraLightResidue const & ligand_residue: ligand_residues_ ) {
					dump_conformer(ligand_residue,sampled_space);
				}
			}

		}


		core::Real accept_ratio =(core::Real)accepted_moves/((core::Real)accepted_moves
			+(core::Real)rejected_moves);
		transform_tracer <<"percent acceptance: "<< accepted_moves << " " << accept_ratio<<" " << rejected_moves <<std::endl;

		protocols::jd2::add_string_real_pair_to_current_job("Transform_accept_ratio", accept_ratio);

		std::string transform_ensemble = "TransformEnsemble";


		//best ligands normally
		for ( core::Size i=1; i <= best_ligands_.size(); ++i ) {
			best_ligands_[i].update_conformation(best_pose.conformation());

		}


		if ( output_sampled_space_ ) {
			sampled_space.close();
		}

		//Output reference poses for comparison
		core::pose::Pose reference_pose = starting_pose;
		for ( core::Size i=1; i <= best_ligands_.size(); ++i ) {
			reference_residues_[i].update_conformation(reference_pose.conformation());
		}
		std::string tag;
		//  reference_pose.dump_pdb("reference_pose.pdb",tag);


		//Reference poses for comparison above

		pose = best_pose;

		transform_tracer << "Accepted pose with grid score: " << best_score << std::endl;

		protocols::jd2::add_string_real_pair_to_current_job("Grid_score", best_score);


		std::map< std::string, core::Real > grid_scores;
		for ( core::Size i=1; i<=transform_info_.jump_ids.size(); ++i ) {
			core::Size jump = transform_info_.jump_ids[i];

			//Add ligand grid scores, hopefully to overall pose output
			utility::map_merge( grid_scores, get_ligand_grid_scores( *grid_set, jump, pose, "" ) );
		}

		for ( auto const & entry : grid_scores ) {
			protocols::jd2::add_string_real_pair_to_current_job( entry.first, entry.second );
		}

	}
}

bool TransformEnsemble::check_grid(qsar::scoring_grid::GridSetCOP grid, utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, core::Real distance) //distance=0 default
{

	if ( distance > transform_info_.box_size ) {
		transform_tracer << "Distance from original center: " << distance << std::endl;
		transform_tracer << "Pose rejected because the center has moved outside the box" << std::endl;
		return false;
	}

	//The score is meaningless if any atoms are outside of the grid
	if ( !grid->is_in_grid(ligand_residues) ) { //Reject the pose

		transform_tracer << "Pose rejected because atoms are outside of the grid" << std::endl;
		return false;
	}

	//Everything is awesome!
	return true;

}

bool TransformEnsemble::monte_carlo(core::Real & current, core::Real & last)
{

	core::Real const boltz_factor((last-current)/transform_info_.temperature);
	core::Real const probability = std::exp( boltz_factor);


	if ( probability < 1 && numeric::random::uniform() >= probability ) {  //reject the new pose
		transform_tracer << "Pose rejected because it didn't meet Metropolis criterion" << std::endl;
		return false;

	} else if ( probability < 1 ) {  // Accept the new pose
		last = current;
		transform_tracer << "Pose accepted because it did meet Metropolis criterion" << std::endl;
		return true;

	} else {  //Accept the new pose
		last = current;

		transform_tracer << "Pose accepted because it improved the score" << std::endl;
		return true;

	}

}

void TransformEnsemble::translate_ligand(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, utility::vector1<core::conformation::UltraLightResidue> & reference_residues, core::Real distance)
{
	//Random uniformly sampled translation sphere with radius equal to distance
	core::Vector translation = numeric::random::uniform_vector_sphere(distance);
	core::Real angle = 360;

	numeric::xyzMatrix<core::Real> rotation(
		numeric::z_rotation_matrix_degrees( angle*numeric::random::uniform() ) * (
		numeric::y_rotation_matrix_degrees( angle*numeric::random::uniform() ) *
		numeric::x_rotation_matrix_degrees( angle*numeric::random::uniform() ) ));

	core::Vector group_center = weighted_center(ligand_residues);
	core::Vector ref_group_center = weighted_center(reference_residues);

	for ( core::Size i=1; i <= ligand_residues.size(); ++i ) {
		ligand_residues[i].transform(rotation,translation,group_center);
	}

	for ( core::Size i=1; i <= reference_residues.size(); ++i ) {
		reference_residues[i].transform(rotation,translation,ref_group_center);
	}

}

void TransformEnsemble::transform_ligand(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, utility::vector1<core::conformation::UltraLightResidue> & reference_residues)
{
	if ( transform_info_.angle ==0 && transform_info_.move_distance == 0 ) {
		transform_tracer.Warning <<"angle and distance are both 0.  Transform will do nothing" <<std::endl;
		return;
	}

	core::Vector translation(
		transform_info_.move_distance*numeric::random::gaussian(),
		transform_info_.move_distance*numeric::random::gaussian(),
		transform_info_.move_distance*numeric::random::gaussian());

	numeric::xyzMatrix<core::Real> rotation(
		numeric::z_rotation_matrix_degrees( transform_info_.angle*numeric::random::gaussian() ) * (
		numeric::y_rotation_matrix_degrees( transform_info_.angle*numeric::random::gaussian() ) *
		numeric::x_rotation_matrix_degrees( transform_info_.angle*numeric::random::gaussian() ) ));

	core::Vector group_center = weighted_center(ligand_residues);
	core::Vector ref_group_center = weighted_center(reference_residues);

	for ( core::Size i=1; i <= ligand_residues.size(); ++i ) {

		ligand_residues[i].transform(rotation,translation,group_center);

	}

	for ( core::Size i=1; i <= reference_residues.size(); ++i ) {

		reference_residues[i].transform(rotation,translation,ref_group_center);
	}

}

void TransformEnsemble::change_conformer(core::conformation::UltraLightResidue & ligand_residue, const core::conformation::UltraLightResidue & reference_residue, core::Size resid)
{

	core::Size index_to_select;

	debug_assert(ligand_conformers_[resid].size());
	index_to_select = numeric::random::random_range(1,ligand_conformers_[resid].size());
	transform_tracer <<"Conformer is number: " << index_to_select << std::endl;
	core::conformation::UltraLightResidue new_residue(ligand_conformers_[resid][index_to_select]);
	new_residue.align_to_residue(reference_residue);
	ligand_residue = new_residue;

}

void TransformEnsemble::change_conformers(utility::vector1<core::conformation::UltraLightResidue> & ligand_residues, const utility::vector1<core::conformation::UltraLightResidue> & reference_residues)
{

	for ( core::Size i=1; i <= ligand_residues.size(); ++i ) {
		change_conformer(ligand_residues[i], reference_residues[i], i);

	}

}

void TransformEnsemble::dump_conformer(core::conformation::UltraLightResidue const & residue, utility::io::ozstream & output)
{
	using namespace ObjexxFCL::format;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms(); ++atom_index ) {
		core::PointPosition coords = residue[atom_index];
		std::string outline( "HETATM" + I( 5,  1 ) + "  V   AAA A"
			+ I( 4, 1 ) + "    "
			+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
			+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
		output <<outline <<std::endl;
	}
}

void TransformEnsemble::print_xyz(core::Vector vector)
{
	transform_tracer <<"X-coordinate is: " << vector[0] << std::endl;
	transform_tracer <<"Y-coordinate is: " << vector[1] << std::endl;
	transform_tracer <<"Z-coordinate is: " << vector[2] << std::endl;

}

core::Vector TransformEnsemble::weighted_center(utility::vector1<core::conformation::UltraLightResidue> & residues)
{
	core::Vector center(0,0,0);
	core::Size total_atoms = 0;

	for ( core::Size i=1; i <=residues.size(); ++i ) {
		center = center + residues[i].center()*residues[i].natoms();
		total_atoms = total_atoms + residues[i].natoms();
	}

	center = center / total_atoms;
	return center;
}

void TransformEnsemble::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	// AMW: In order to make non_negative_real possible with the least pain
	// it's easiest to make it come from decimal.
	// This means we're limited to 12 digits and no scientific notation.

	// To do: change this.
	XMLSchemaRestriction restriction;
	restriction.name( "non_negative_real" );
	restriction.base_type( xs_decimal );
	restriction.add_restriction( xsr_minInclusive, "0" );
	xsd.add_top_level_element( restriction );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("chains", xs_string, "Ligand chains, specified as the PDB chain IDs")
		+ XMLSchemaAttribute("sampled_space_file", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute::required_attribute("move_distance", xsct_real, "Maximum translation performed per step in the monte carlo search.")
		+ XMLSchemaAttribute::required_attribute("box_size", xsct_real, "Maximum translation that can occur from the ligand starting point.")
		+ XMLSchemaAttribute::required_attribute("angle", xsct_real, "Maximum rotation angle performed per step in the monte carlo search.")
		+ XMLSchemaAttribute::required_attribute("temperature", xsct_real, "Boltzmann temperature for the monte carlo simulation.")
		+ XMLSchemaAttribute::required_attribute("cycles", xsct_non_negative_integer,
		"Total number of steps to be performed in the monte carlo simulation.")
		+ XMLSchemaAttribute::attribute_w_default("repeats", xsct_non_negative_integer,
		"Total number of repeats of the monte carlo simulation to be performed.", "1")
		+ XMLSchemaAttribute::attribute_w_default("optimize_until_score_is_negative", xsct_rosetta_bool,
		"Continue sampling beyond \"cycles\" if score is positive", "false")
		+ XMLSchemaAttribute::attribute_w_default("optimize_until_ideal", xsct_rosetta_bool,
		"Continue sampling beyond \"cycles\" if score not close to minimum - all atoms has -1 score", "false")
		+ XMLSchemaAttribute::attribute_w_default("use_conformers", xsct_rosetta_bool,
		"Use ligand conformations while sampling", "true")
		+ XMLSchemaAttribute::attribute_w_default("initial_perturb", "non_negative_real" ,
		"Make an initial, unscored translation and rotation "
		"Translation will be selected uniformly in a sphere of the given radius (Angstrom)."
		"Automatically triggers 360 degrees randomized rotation", "0.0");

	protocols::qsar::scoring_grid::attributes_for_parse_grid_set_from_tag(attlist, "The Scoring Grid set to use with TransformEnsemble scoring");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Performs a monte carlo search of the ligand ensemble binding site using precomputed scoring grids. "
		"Replaces the Translate, Rotate, and SlideTogether movers.", attlist );
}


protocols::moves::MoverOP TransformEnsembleCreator::create_mover() const
{
	return protocols::moves::MoverOP( new TransformEnsemble );
}

std::string TransformEnsemble::mover_name()
{
	return "TransformEnsemble";
}

std::string TransformEnsemble::get_name() const {
	return mover_name();
}

std::string TransformEnsembleCreator::keyname() const {
	return TransformEnsemble::mover_name();
}

void TransformEnsembleCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TransformEnsemble::provide_xml_schema( xsd );
}


}
}
