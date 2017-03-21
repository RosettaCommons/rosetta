// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/ligand_docking/Transform.cc
/// @author Sam DeLuca

// Testing

#include <protocols/ligand_docking/TransformCreator.hh>
#include <protocols/ligand_docking/Transform.hh>

#include <protocols/ligand_docking/grid_functions.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>

#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/random/random_xyz.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/conversions.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/format.hh>
#include <utility/io/ozstream.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace ligand_docking {


static THREAD_LOCAL basic::Tracer transform_tracer( "protocols.ligand_docking.Transform" );

// XRW TEMP std::string TransformCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return Transform::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP TransformCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new Transform );
// XRW TEMP }

// XRW TEMP std::string Transform::mover_name()
// XRW TEMP {
// XRW TEMP  return "Transform";
// XRW TEMP }

Transform::Transform():
	Mover("Transform")
	// & in class defaults
{}

Transform::Transform(Transform const & ) = default;

Transform::Transform(
	std::string const & chain,
	core::Real const & box_size,
	core::Real const & move_distance,
	core::Real const & angle,
	core::Size const & cycles,
	core::Real const & temp
) : Mover("Transform")
	// & in class defaults

{
	transform_info_.chain = chain;
	transform_info_.box_size = box_size;
	transform_info_.move_distance = move_distance;
	transform_info_.angle = angle;
	transform_info_.cycles = cycles;
	transform_info_.temperature = temp;
}

Transform::~Transform() = default;

protocols::moves::MoverOP Transform::clone() const
{
	return protocols::moves::MoverOP( new Transform (*this) );
}

protocols::moves::MoverOP Transform::fresh_instance() const
{
	return protocols::moves::MoverOP( new Transform );
}

// XRW TEMP std::string Transform::get_name() const
// XRW TEMP {
// XRW TEMP  return "Transform";
// XRW TEMP }

void Transform::parse_my_tag
(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data_map*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "Transform" ) {
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires chain tag");
	if ( ! tag->hasOption("move_distance") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires move_distance tag");
	if ( ! tag->hasOption("box_size") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires box_size tag");
	if ( ! tag->hasOption("angle") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires angle tag");
	if ( ! tag->hasOption("cycles") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires cycles tag");
	if ( !tag->hasOption("temperature") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires temperature tag");

	transform_info_.chain = tag->getOption<std::string>("chain");

	//Divides by root(3) so the center can only move a total equal to move_distance in each step
	transform_info_.move_distance = (tag->getOption<core::Real>("move_distance")); // sqrt(3);

	transform_info_.box_size = tag->getOption<core::Real>("box_size");
	transform_info_.angle = tag->getOption<core::Real>("angle");
	transform_info_.cycles = tag->getOption<core::Size>("cycles");
	transform_info_.temperature = tag->getOption<core::Real>("temperature");
	transform_info_.repeats = tag->getOption<core::Size>("repeats",1);
	optimize_until_score_is_negative_ = tag->getOption<bool>("optimize_until_score_is_negative",false);

	use_conformers_ = tag->getOption<bool>("use_conformers",true);

	initial_perturb_ = tag->getOption<core::Real>("initial_perturb",0.0);
	if ( initial_perturb_ < 0 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("The initial_perturb option to the Transform mover must be positive.");
	}
	if ( tag->hasOption("initial_angle_perturb") ) {
		initial_angle_perturb_ = tag->getOption<core::Real>("initial_angle_perturb",0.0);
		if ( initial_angle_perturb_ < 0 ) {
			throw utility::excn::EXCN_RosettaScriptsOption("The initial_angle_perturb option to the Transform mover must be positive.");
		}
	} // else leave as default: 360 degree sampling

	if ( tag->hasOption("rmsd") ) {
		check_rmsd_ = true;
		transform_info_.rmsd = tag->getOption<core::Real>("rmsd");
	}

	if ( tag->hasOption("sampled_space_file") ) {
		output_sampled_space_ = true;
		sampled_space_file_ = tag->getOption<std::string>("sampled_space_file");
	}

}

void Transform::apply(core::pose::Pose & pose)
{
	qsar::scoring_grid::GridManager* grid_manager(qsar::scoring_grid::GridManager::get_instance());

	debug_assert(transform_info_.chain.size() == 1);
	transform_info_.chain_id = core::pose::get_chain_id_from_chain(transform_info_.chain, pose);
	transform_info_.jump_id = core::pose::get_jump_id_from_chain_id(transform_info_.chain_id, pose);
	core::Size const begin(pose.conformation().chain_begin(transform_info_.chain_id));
	core::Vector const center(protocols::geometry::downstream_centroid_by_jump(pose, transform_info_.jump_id));
	debug_assert(grid_manager != nullptr); //something has gone hopelessly wrong if this triggers

	core::conformation::Residue original_residue = pose.residue(begin);
	core::chemical::ResidueType residue_type = pose.residue_type(begin);

	grid_manager->initialize_all_grids(center);
	grid_manager->update_grids(pose,center);

	core::Real last_score(10000.0);
	core::Real best_score(10000.0);
	core::Size accepted_moves = 0;
	core::Size rejected_moves = 0;
	core::Size outside_grid_moves = 0;

	core::pose::Pose best_pose(pose);

	//Setup UltraLight residues for docking movements
	core::conformation::UltraLightResidue original_ligand(pose.residue(begin).get_self_ptr());
	core::conformation::UltraLightResidue ligand_residue = original_ligand;
	core::conformation::UltraLightResidue best_ligand = ligand_residue;

	core::Real temperature = transform_info_.temperature;
	core::Vector original_center(original_residue.xyz(original_residue.nbr_atom()));

	setup_conformers(pose, begin);

	utility::io::ozstream sampled_space;
	if ( output_sampled_space_ ) {
		sampled_space.open(sampled_space_file_);
	}


	for ( core::Size repeat = 1; repeat <= transform_info_.repeats; ++repeat ) {
		core::Size cycle = 1;
		bool not_converged = true;
		ligand_residue = original_ligand;
		core::conformation::UltraLightResidue last_accepted_ligand_residue = ligand_residue;

		//Initial Perturbation for benchmarking purposes or sampling a large binding volume
		// For benchmarking purposes it is sometimes desirable to translate the ligand
		// away from the starting point and randomize its orientation before beginning a trajectory.
		// The defaults of 0.0 & -360.0 won't trigger, but if either are set to positive they will.

		//Setting an initial perturb will also randomize the startign conformer

		if ( initial_perturb_ > 0.0 || initial_angle_perturb_ > 0 ) {
			bool perturbed = false;
			while ( !perturbed )
					{
				perturbed=true;
				randomize_ligand( ligand_residue, initial_perturb_, initial_angle_perturb_ );

				//Also randomize starting conformer to remove bias
				if ( ligand_conformers_.size() > 1 && use_conformers_ == true ) {
					change_conformer(ligand_residue);
				}
				core::Vector new_center(ligand_residue.center());
				core::Real distance = new_center.distance(original_center);


				//Not everything in grid, also checks center distance
				if ( !check_grid(grid_manager, ligand_residue, distance) ) {
					ligand_residue = last_accepted_ligand_residue;
					perturbed = false;
					continue;
				}
			}
			last_accepted_ligand_residue = ligand_residue;
		}

		last_score = grid_manager->total_score(ligand_residue);

		//Define starting ligand model (after perturb) as best model/score
		best_score = last_score;
		best_ligand = ligand_residue;

		while ( not_converged )
				{

			if ( optimize_until_score_is_negative_ ) {
				if ( cycle >= transform_info_.cycles && last_score <= 0.0 ) {

					not_converged= false;
				}
			} else {
				if ( cycle >= transform_info_.cycles ) {
					not_converged= false;
				}
			}

			cycle++;

			//during each move either move the ligand or try a new conformer (if there is more than one conformer)
			if ( ligand_conformers_.size() > 1 && use_conformers_ == true ) {
				if ( numeric::random::rg().uniform() >= 0.5 ) {
					transform_ligand(ligand_residue);
				} else {
					change_conformer(ligand_residue);
				}

			} else {
				transform_ligand(ligand_residue);
			}

			core::Vector new_center(ligand_residue.center());
			core::Real distance = new_center.distance(original_center);

			//Not everything in grid - Check distance too
			if ( !check_grid(grid_manager, ligand_residue, distance) ) {
				ligand_residue = last_accepted_ligand_residue;
				rejected_moves++;
				outside_grid_moves++;
				continue;
			}

			core::Real current_score = grid_manager->total_score(ligand_residue);
			core::Real const boltz_factor((last_score-current_score)/temperature);
			core::Real const probability = std::exp( boltz_factor ) ;


			if ( output_sampled_space_ ) {
				dump_conformer(ligand_residue,sampled_space);
			}

			if ( check_rmsd_ && !check_rmsd(original_ligand, ligand_residue) ) { //reject the new pose
				ligand_residue = last_accepted_ligand_residue;
				rejected_moves++;
			} else if ( probability < 1 && numeric::random::rg().uniform() >= probability ) {  //reject the new pose
				ligand_residue = last_accepted_ligand_residue;
				rejected_moves++;
				//  transform_tracer << "Move rejected- did not meet Monte Carlo probability " << std::endl;

			} else if ( probability < 1 ) {  // Accept the new pose
				last_score = current_score;
				last_accepted_ligand_residue = ligand_residue;
				accepted_moves++;

			} else {  //Accept the new pose
				last_score = current_score;
				last_accepted_ligand_residue = ligand_residue;
				accepted_moves++;
			}


			if ( last_score <= best_score ) {
				best_score = last_score;
				best_ligand = last_accepted_ligand_residue;
				//   transform_tracer << "accepting new best pose" << std::endl;
			} else {
				//   transform_tracer << "not accepting new best pose" << std::endl;

			}

		}


		core::Real accept_ratio =(core::Real)accepted_moves/((core::Real)accepted_moves+(core::Real)rejected_moves);
		transform_tracer <<"percent acceptance: "<< accepted_moves << " " << accept_ratio<<" " << rejected_moves <<std::endl;
		if ( outside_grid_moves > 0 ) {
			core::Real outside_grid_ratio = (core::Real)outside_grid_moves/((core::Real)accepted_moves+(core::Real)rejected_moves);
			transform_tracer << "Moves rejected for being outside of grid: " << outside_grid_moves << "  " << outside_grid_ratio << std::endl;
		}

		jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("Transform_accept_ratio", accept_ratio);
		best_ligand.update_conformation(best_pose.conformation());
	}

	if ( output_sampled_space_ ) {
		sampled_space.close();
	}
	pose = best_pose;

	transform_tracer << "Accepted pose with grid score: " << best_score << std::endl;
	jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("Grid_score", best_score);

}


void Transform::randomize_ligand(core::conformation::UltraLightResidue & residue, core::Real distance, core::Real angle)
{
	// Pick a random direction, then translate a random distance in that direction, up to the given maximum
	core::Vector trans_axis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );

	// sampling with sqrt( rnd * r^2 ) is to give an equal sampling distribution in volume (as opposed to equal sampling of distance
	core::Vector translation( std::sqrt( numeric::random::rg().uniform() * distance * distance ) * trans_axis );

	// Pick a (new) random axis, then rotate around that axis by up to the given maximum.
	// (Opposite direction rotation is handled by positive rotation around the opposite-direction vector)
	core::Vector axis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );

	core::Real chosen_angle( angle*numeric::random::rg().uniform() );

	numeric::xyzMatrix<core::Real> rotation(
		numeric::rotation_matrix( axis, numeric::conversions::to_radians( chosen_angle ) ) );

	residue.transform(rotation,translation,residue.center());

}

void Transform::transform_ligand(core::conformation::UltraLightResidue & residue)
{
	if ( transform_info_.angle ==0 && transform_info_.move_distance == 0 ) {
		transform_tracer.Warning << "angle and distance are both 0.  Transform will do nothing" <<std::endl;
		return;
	}

	core::Vector translation(
		transform_info_.move_distance*numeric::random::rg().gaussian(),
		transform_info_.move_distance*numeric::random::rg().gaussian(),
		transform_info_.move_distance*numeric::random::rg().gaussian());

	numeric::xyzMatrix<core::Real> rotation(
		numeric::z_rotation_matrix_degrees( transform_info_.angle*numeric::random::rg().gaussian() ) * (
		numeric::y_rotation_matrix_degrees( transform_info_.angle*numeric::random::rg().gaussian() ) *
		numeric::x_rotation_matrix_degrees( transform_info_.angle*numeric::random::rg().gaussian() ) ));

	residue.transform(rotation,translation,residue.center());

}

bool Transform::check_grid(qsar::scoring_grid::GridManager* grid, core::conformation::UltraLightResidue & ligand_residue, core::Real distance) //distance=0 default
{

	if ( distance > transform_info_.box_size ) {
		return false;
	}

	//The score is meaningless if any atoms are outside of the grid
	if ( !grid->is_in_grid(ligand_residue) ) { //Reject the pose
		return false;
	}

	//Everything is awesome!
	return true;

}

void Transform::setup_conformers(core::pose::Pose & pose, core::Size begin)
{
	ligand_conformers_.clear();
	rotamers_for_trials(pose,begin,ligand_conformers_);
	transform_tracer << "Considering " << ligand_conformers_.size() << " conformers during sampling" << std::endl;
}

void Transform::change_conformer(core::conformation::UltraLightResidue & residue)
{
	debug_assert(ligand_conformers_.size());
	core::Size index_to_select = numeric::random::rg().random_range(1,ligand_conformers_.size());
	core::conformation::UltraLightResidue new_residue(ligand_conformers_[index_to_select]);
	new_residue.align_to_residue(residue);
	residue = new_residue;

}

void Transform::dump_conformer(core::conformation::UltraLightResidue & residue, utility::io::ozstream & output)
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

bool Transform::check_rmsd(core::conformation::UltraLightResidue const & start, core::conformation::UltraLightResidue const& current) const
{
	debug_assert(start.natoms() == current.natoms());

	core::Real total_distance =0.0;
	for ( core::Size atomno = 1; atomno <= start.natoms(); ++atomno ) {
		total_distance += start[atomno].distance(current[atomno]);
	}

	core::Real rmsd = sqrt(total_distance/start.natoms());

	if ( rmsd <= transform_info_.rmsd ) {
		return true;
	} else {
		return false;
	}

}

std::string Transform::get_name() const {
	return mover_name();
}

std::string Transform::mover_name() {
	return "Transform";
}

void Transform::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
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
		+ XMLSchemaAttribute::required_attribute("chain", xs_string, "The ligand chain, specified as the PDB chain ID")
		+ XMLSchemaAttribute("sampled_space_file", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute::required_attribute("move_distance", xsct_real, "Maximum translation performed per step in the monte carlo search.")
		+ XMLSchemaAttribute::required_attribute("box_size", xsct_real, "Maximum translation that can occur from the ligand starting point.")
		+ XMLSchemaAttribute::required_attribute("angle", xsct_real, "Maximum rotation angle performed per step in the monte carlo search.")
		+ XMLSchemaAttribute::required_attribute("temperature", xsct_real, "Boltzmann temperature for the monte carlo simulation.")
		+ XMLSchemaAttribute("rmsd", xsct_real, "Maximum RMSD to be sampled away from the starting position.")
		+ XMLSchemaAttribute::required_attribute("cycles", xsct_non_negative_integer,
		"Total number of steps to be performed in the monte carlo simulation.")
		+ XMLSchemaAttribute::attribute_w_default("repeats", xsct_non_negative_integer,
		"Total number of repeats of the monte carlo simulation to be performed.", "1")
		+ XMLSchemaAttribute::attribute_w_default("optimize_until_score_is_negative", xsct_rosetta_bool,
		"Continue sampling beyond \"cycles\" if score is positive", "false")
		+ XMLSchemaAttribute::attribute_w_default("initial_perturb", "non_negative_real" ,
		"Make an initial, unscored translation and rotation "
		"Translation will be selected uniformly in a sphere of the given radius (Angstrom)."
		"Triggers 360 degress rotations are triggered.", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("initial_angle_perturb", "non_negative_real", "Control the size of the rotational perturbation by intitial_perturb. (Degrees)", "0.0");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Performs a monte carlo search of the ligand binding site using precomputed scoring grids. "
		"Replaces the Translate, Rotate, and SlideTogether movers.", attlist );
}

std::string TransformCreator::keyname() const {
	return Transform::mover_name();
}

protocols::moves::MoverOP
TransformCreator::create_mover() const {
	return protocols::moves::MoverOP( new Transform );
}

void TransformCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	Transform::provide_xml_schema( xsd );
}

}
}
