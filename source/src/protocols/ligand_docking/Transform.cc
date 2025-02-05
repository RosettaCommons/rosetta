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
#include <protocols/jd2/util.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>
#include <protocols/qsar/scoring_grid/schema_util.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/rms_util.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/chains_util.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/palette/NoDesignPackerPalette.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/conversions.hh>

#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/io/util.hh>
#include <ObjexxFCL/format.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/CitationManager.hh>

namespace protocols {
namespace ligand_docking {


static basic::Tracer TR( "protocols.ligand_docking.Transform" );

Transform::Transform():
	Mover("Transform")
	// & in class defaults
{}

Transform::Transform(Transform const & ) = default;

Transform::Transform(
	protocols::qsar::scoring_grid::GridSetCOP grid_set_prototype,
	std::string const & chain,
	core::Real const & box_size,
	core::Real const & move_distance,
	core::Real const & angle,
	core::Size const & cycles,
	core::Real const & temp
) : Mover("Transform"),
	grid_set_prototype_(std::move( grid_set_prototype ))
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
	return utility::pointer::make_shared< Transform > (*this);
}

protocols::moves::MoverOP Transform::fresh_instance() const
{
	return utility::pointer::make_shared< Transform >();
}

void Transform::parse_my_tag
(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	if ( tag->getName() != "Transform" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible");
	}

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
	use_constraints_ = tag->getOption<bool>("use_constraints",false);

	if ( use_constraints_ ) {
		//XML file takes precedence over command line option
		if ( tag->hasOption("cst_fa_file") ) {
			cst_fa_file_ = tag->getOption<std::string>("cst_fa_file");
		} else if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file].user() ) {
			cst_fa_file_ = basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file]().front();
		}

		if ( cst_fa_file_.empty() ) utility_exit_with_message("use_constraints=true requires cst_fa_file tag or cst_fa_file command-line option");
	}

	//XML file takes precedence over command line option
	if ( tag->hasOption("cst_fa_weight") ) {
		cst_fa_weight_ = tag->getOption<core::Real>("cst_fa_weight");
	} else if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_weight].user() ) {
		cst_fa_weight_ = basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_weight]();
	}

	use_main_model_ = tag->getOption<bool>("use_main_model",false);

	initial_perturb_ = tag->getOption<core::Real>("initial_perturb",0.0);
	if ( initial_perturb_ < 0 ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "The initial_perturb option to the Transform mover must be positive.");
	}
	if ( tag->hasOption("initial_angle_perturb") ) {
		initial_angle_perturb_ = tag->getOption<core::Real>("initial_angle_perturb",0.0);
		if ( initial_angle_perturb_ < 0 ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "The initial_angle_perturb option to the Transform mover must be positive.");
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

	grid_set_prototype_ = protocols::qsar::scoring_grid::parse_grid_set_from_tag(tag, data);

	if ( grid_set_prototype_ == nullptr ) {
		utility_exit_with_message( "The Transform mover needs to have the GridSet prototype set in order to work." );
	}

	if ( tag->hasOption("ensemble_proteins") ) {

		ensemble_proteins_ =  tag->getOption<std::string>("ensemble_proteins");
		utility::vector1<std::string> lines = utility::io::get_lines_from_file_data(ensemble_proteins_);

		for ( const auto & line : lines ) {
			core::pose::PoseOP extra_pose = core::import_pose::pose_from_file(line,false, core::import_pose::PDB_file);
			grid_set_poses_.insert(std::make_pair(line, *extra_pose));
		}
	}

}

void Transform::apply(core::pose::Pose & pose)
{
	grid_sets_.clear();

	debug_assert(transform_info_.chain.size() == 1);
	core::Size chain_id = core::pose::get_chain_id_from_chain(transform_info_.chain, pose);
	core::Size jump_id = core::pose::get_jump_id_from_chain_id(chain_id, pose);
	core::Vector const center(protocols::geometry::downstream_centroid_by_jump(pose, jump_id));

	core::Size const begin(pose.conformation().chain_begin(chain_id));

	core::conformation::Residue original_residue = pose.residue(begin);
	//core::chemical::ResidueType const & residue_type = pose.residue_type(begin);

	//Duplicate pose to avoid overriding existing constraints
	core::pose::PoseOP cst_pose;
	core::scoring::ScoreFunctionOP cst_function(new core::scoring::ScoreFunction);

	//Setup constraints
	if ( use_constraints_ ) {
		cst_pose = pose.clone();
		core::scoring::constraints::ConstraintSetOP cstset_ = core::scoring::constraints::ConstraintIO::get_instance()->read_constraints(cst_fa_file_, utility::pointer::make_shared< core::scoring::constraints::ConstraintSet >(), *cst_pose);
		cst_pose->constraint_set( cstset_ );

		cst_function->set_weight( core::scoring::atom_pair_constraint, cst_fa_weight_ );
		cst_function->set_weight( core::scoring::angle_constraint,  cst_fa_weight_ );
		cst_function->set_weight( core::scoring::dihedral_constraint,  cst_fa_weight_  );
		cst_function->set_weight( core::scoring::coordinate_constraint,  cst_fa_weight_ );
	}

	if ( grid_set_prototype_ == nullptr ) {
		utility_exit_with_message( "The Transform mover needs to have the GridSet prototype set in order to work." );
	}

	grid_sets_.push_back(std::make_pair("MAIN", qsar::scoring_grid::GridManager::get_instance()->get_grids( *grid_set_prototype_, pose, center, transform_info_.chain )));

	if ( ensemble_proteins_ != "" ) {
		make_multi_pose_grids(center);
	}

	core::Real last_score(10000.0);
	core::Real best_score(10000.0);
	core::Size accepted_moves = 0;
	core::Size rejected_moves = 0;
	core::Size outside_grid_moves = 0;
    core::Size rmsd_rejects = 0;

	// core::pose::Pose best_pose(pose);

	//Setup UltraLight residues for docking movements
	core::conformation::UltraLightResidue original_ligand(pose.residue(begin).get_self_ptr());
	core::conformation::UltraLightResidue ligand_residue = original_ligand;
	core::conformation::UltraLightResidue best_ligand = ligand_residue;

	core::Real temperature = transform_info_.temperature;
	core::Vector original_center(original_residue.xyz(original_residue.nbr_atom()));

	setup_conformers(pose, begin);

	// Check conformers, to make sure that they'll fit in the grid (at least in the initial position)
	if ( !check_conformers( original_ligand ) ) {
		// Already printed error message
		set_last_move_status( protocols::moves::FAIL_RETRY );
		return;
	}

	utility::io::ozstream sampled_space;
	if ( output_sampled_space_ ) {
		sampled_space.open(sampled_space_file_);
	}

	for ( core::Size repeat = 1; repeat <= transform_info_.repeats; ++repeat ) {
		TR.Trace << "Starting transform repeat " << repeat << " of " << transform_info_.repeats << std::endl;
		ligand_residue = original_ligand;
		core::conformation::UltraLightResidue last_accepted_ligand_residue = ligand_residue;

		//Initial Perturbation for benchmarking purposes or sampling a large binding volume
		// For benchmarking purposes it is sometimes desirable to translate the ligand
		// away from the starting point and randomize its orientation before beginning a trajectory.
		// The defaults of 0.0 & -360.0 won't trigger, but if either are set to positive they will.

		//Setting an initial perturb will also randomize the starting conformer

		if ( initial_perturb_ > 0.0 || initial_angle_perturb_ > 0 ) {
			bool perturbed = false;
			core::Size n_perturb_trials(0);
			while ( !perturbed && n_perturb_trials < 10*ligand_conformers_.size() ) {
				if ( initial_perturb_ > transform_info_.box_size ) {
					TR.Warning << "[ Warning ] In the Transform mover, the initial perturbation size is larger than the box size. This is highly inefficient." << std::endl;
				}

				++n_perturb_trials;
				perturbed=true;
				randomize_ligand( ligand_residue, initial_perturb_, initial_angle_perturb_ );

				//Also randomize starting conformer to remove bias
				if ( ligand_conformers_.size() > 1 && use_conformers_ == true ) {
					change_conformer(ligand_residue);
				}
				core::Vector new_center(ligand_residue.center());
				core::Real distance = new_center.distance(original_center);

				//Not everything in grid, also checks center distance
				if ( !check_grid(ligand_residue, distance) ) {
					TR.Debug << "In the initial perturbation, the ligand moved outside the grid - retrying." << std::endl;
					ligand_residue = last_accepted_ligand_residue;
					perturbed = false;
					continue;
				}
			}
			if ( n_perturb_trials >= 10*ligand_conformers_.size() ) {
				TR.Error << "Could not get a decent initial perturbation - check grid size, ligand size, and initial perturbation size." << std::endl;
				TR.Error << "    For this system, with a perturbation of " << initial_perturb_
					<< ", a grid size of at least " << utility::Real2string(recommended_grid_size(),1) << " is recommended." << std::endl;
				set_last_move_status( protocols::moves::FAIL_RETRY );
				return;
			} else if ( n_perturb_trials > 10 ) {
				TR.Warning << "It took more than 10 tries to get a decent initial perturbation. You likely want to check your grid size and initial perturbation size." << std::endl;
				TR.Warning << "    With the current perturbation setting of " << initial_perturb_ << ", a grid size of at least " << utility::Real2string(recommended_grid_size(),1) << " is recommended." << std::endl;
			}
			last_accepted_ligand_residue = ligand_residue;
		}

		last_score = score_ligand(ligand_residue);

		if ( use_constraints_ ) {
			last_score += score_constraints(*cst_pose, ligand_residue, cst_function);
		}

		if ( last_score < best_score ) {
			// best_ligand is accumulated across all repeats
			best_score = last_score;
			best_ligand = ligand_residue;
		}

		core::Size cycle = 1;
		bool not_converged = true;

		while ( not_converged ) {

			if ( optimize_until_score_is_negative_ ) {
				if ( cycle >= transform_info_.cycles && last_score <= 0.0 ) {
					not_converged= false;
				} else if ( cycle % 2*transform_info_.cycles == 0 ) { // Print every time we're twice the requested cycles.
					// Print a warning, so at least we can see if we're in an infinite loop.
					TR.Warning << "optimized for " << cycle << " cycles and the score (" << last_score << ") is still not negative." << std::endl;
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
			if ( !check_grid(ligand_residue, distance) ) {
				ligand_residue = last_accepted_ligand_residue;
				rejected_moves++;
				outside_grid_moves++;
				continue;
			}

			core::Real current_score = score_ligand(ligand_residue);

			if ( use_constraints_ ) {
				current_score += score_constraints(*cst_pose, ligand_residue, cst_function);
			}

			core::Real const boltz_factor((last_score-current_score)/temperature);
			core::Real const probability = std::exp( boltz_factor ) ;


			if ( output_sampled_space_ ) {
				dump_conformer(ligand_residue,sampled_space);
			}

			if ( check_rmsd_ && !check_rmsd(original_ligand, ligand_residue) ) { //reject the new pose
				ligand_residue = last_accepted_ligand_residue;
				rejected_moves++;
                rmsd_rejects++;
			} else if ( probability < 1 && numeric::random::rg().uniform() >= probability ) {  //reject the new pose
				ligand_residue = last_accepted_ligand_residue;
				rejected_moves++;
				//  TR << "Move rejected- did not meet Monte Carlo probability " << std::endl;

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
				//   TR << "accepting new best pose" << std::endl;
			} else {
				//   TR << "not accepting new best pose" << std::endl;
			}

		}


	}

	if ( output_sampled_space_ ) {
		sampled_space.close();
	}

	//Accept ratio over all repeats
	core::Real accept_ratio =(core::Real)accepted_moves/((core::Real)accepted_moves+(core::Real)rejected_moves);
	TR <<"percent acceptance: "<< accepted_moves << " " << accept_ratio <<" " << rejected_moves <<std::endl;
	if ( outside_grid_moves > 0 ) {
		core::Real outside_grid_ratio = (core::Real)outside_grid_moves/((core::Real)accepted_moves+(core::Real)rejected_moves);
		TR << "Moves rejected for being outside of grid: " << outside_grid_moves << "  " << outside_grid_ratio << std::endl;
		if ( outside_grid_ratio > 0.05 ) { // 5% is rather arbitrary here
			TR.Warning << "A large number of moves were rejected for being outside the grid. You likely want to reexamine your settings." << std::endl;
			TR.Warning << "    For the current settings, a grid size of at least " << utility::Real2string(recommended_grid_size( accept_ratio ),1);
			TR.Warning << "    and a box size of at least " << utility::Real2string(recommended_box_size( accept_ratio ),1) << " are recommended." << std::endl;
		}
	}
    if ( rmsd_rejects > 0 ) {
        core::Real rmsd_reject_ratio = (core::Real)rmsd_rejects / ( (core::Real)accepted_moves + (core::Real)rejected_moves );
        TR << "Moves rejected for going beyond max RMSD (" << transform_info_.rmsd << "): " << rmsd_rejects << "  " << rmsd_reject_ratio << std::endl;
    }

	protocols::jd2::add_string_real_pair_to_current_job("Transform_accept_ratio", accept_ratio);

	//Set pose to best scored model - report score of a single grid (the best one) if using multiple poses
	best_score = convert_to_full_pose(pose, best_ligand, chain_id);
	TR << "Accepted pose with grid score: " << best_score << std::endl;
	protocols::jd2::add_string_real_pair_to_current_job("Grid_score", best_score);
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
		TR.Warning << "angle and distance are both 0.  Transform will do nothing" <<std::endl;
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

bool Transform::check_grid(core::conformation::UltraLightResidue & ligand_residue, core::Real distance) //distance=0 default
{

	if ( distance > transform_info_.box_size ) {
		return false;
	}

	//The score is meaningless if any atoms are outside of the grid
	for ( const auto & grid_set : grid_sets_ ) {
		if ( !grid_set.second->is_in_grid(ligand_residue) ) { //Reject the pose
			return false;
		}
	}
	//Everything is awesome!
	return true;

}

void Transform::setup_conformers(core::pose::Pose & pose, core::Size begin)
{
	using namespace core::conformation;

	utility::vector1< core::conformation::ResidueOP > ligand_confs;
	core::pack::task::PackerTaskCOP the_task( core::pack::task::TaskFactory::create_packer_task(pose, utility::pointer::make_shared< core::pack::palette::NoDesignPackerPalette >() ) );
	rotamers_for_trials(pose,begin,ligand_confs, *the_task);

	ligand_conformers_.clear();
	for ( core::conformation::ResidueOP const & conf: ligand_confs ) {
		ligand_conformers_.push_back( utility::pointer::make_shared< UltraLightResidue >( conf ) );
	}
	if ( ligand_conformers_.empty() ) {
		// Add the starting conformer, so we at least have the one.
		ligand_conformers_.push_back( utility::pointer::make_shared< UltraLightResidue >( pose.residue(begin).get_self_ptr() ) );
	}
	TR << "Considering " << ligand_conformers_.size() << " conformers during sampling" << std::endl;
}

bool Transform::check_conformers(core::conformation::UltraLightResidue & starting_residue ) const
{
	core::Size n_outside( 0 );
	for ( core::conformation::UltraLightResidueOP conf: ligand_conformers_ ) {
		core::conformation::UltraLightResidue lig( *conf );
		lig.align_to_residue(starting_residue);
		for ( const auto & grid_set : grid_sets_ ) {
			if ( ! grid_set.second->is_in_grid( lig ) ) {
				++n_outside;
				break;
			}
		}
	}
	if ( n_outside == ligand_conformers_.size() ) {
		TR.Error << "All conformers start with atoms outside the grid -- Docking will be impossible. Increase the size of the grid." << std::endl;
		TR.Error << "(For the current system and settings, a grid size of at least " << utility::Real2string(recommended_grid_size(),1) << " is recommended.)" << std::endl;
		return false;
	} else if ( n_outside > ligand_conformers_.size()/2 ) {
		TR.Warning << n_outside << " of " << ligand_conformers_.size() << " conformers start with atoms outside the grid." << std::endl;
		TR.Warning << "     You likely want to increase the grid size to accomodate the size of the ligand. Recommended grid size for this run: at least " << utility::Real2string(recommended_grid_size(),1) << std::endl;
	}
	return true;
}

core::Real Transform::recommended_grid_size( core::Real success_rate ) const
{
	core::Real ligand_width( 0 );
	for ( core::conformation::UltraLightResidueOP conf: ligand_conformers_ ) {
		core::Real maxdist( conf->max_dist_to_center() );
		if ( maxdist > ligand_width ) {
			ligand_width = maxdist;
		}
	}

	core::Real perturb_size = initial_perturb_;
	core::Real estimated_travel = estimate_mc_travel( success_rate );
	core::Real total_travel = ligand_width + perturb_size + estimated_travel;
	return 2.0*total_travel; // total_travel is radius, recommended size is diameter
}

core::Real Transform::recommended_box_size( core::Real success_rate ) const
{
	// We don't need the ligand width here, as we're only concerned about the center.
	core::Real perturb_size = initial_perturb_;
	if ( success_rate < 0.33 ) { success_rate = 0.33; }
	core::Real estimated_travel = estimate_mc_travel( success_rate );
	core::Real total_travel = perturb_size + estimated_travel;
	return 2.0*total_travel; // total_travel is radius, recommended size is diameter
}

core::Real Transform::estimate_mc_travel( core::Real success_rate )  const
{
	if ( success_rate < 0.33 ) { // An accept ratio of 0.33 as a default/minimum is somewhat arbitrary.
		success_rate = 0.03;
	}
	// Estimate for how many translational moves we make
	core::Size effective_moves = success_rate * transform_info_.cycles;
	if ( ligand_conformers_.size() > 1 && use_conformers_ == true ) {
		effective_moves *= 0.5; // Reduce based on conformer switches.
	}
	// We can consider each axis separately, which is a gaussian random walk with sd = move_distance
	core::Real estimated_travel = std::sqrt( effective_moves ) * transform_info_.move_distance;
	return estimated_travel;
}

void Transform::change_conformer(core::conformation::UltraLightResidue & residue)
{
	debug_assert(ligand_conformers_.size());
	core::Size index_to_select = numeric::random::rg().random_range(1,ligand_conformers_.size());
	core::conformation::UltraLightResidue new_residue(*ligand_conformers_[index_to_select]);
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

	core::Real rmsd = core::scoring::automorphic_rmsd( *start.residue(), *current.residue(), /*superimpose=*/false );

    TR.Debug << "Current rmsd " << rmsd << ". Limit: " << transform_info_.rmsd << std::endl;

	if ( rmsd <= transform_info_.rmsd ) {
		return true;
	} else {
		return false;
	}

}

core::Real Transform::score_constraints(core::pose::Pose & pose, core::conformation::UltraLightResidue & residue, core::scoring::ScoreFunctionOP & sfxn)
{
	residue.update_conformation(pose.conformation());
	return (*sfxn)(pose);
}

core::Real Transform::score_ligand(core::conformation::UltraLightResidue & residue)
{
	core::Real score_total = 0;
	for ( const auto & grid_set : grid_sets_ ) {
		core::Real score = grid_set.second->total_score(residue);
		score_total = score_total + score;
		TR.Debug << "Scored ligand for " << grid_set.first << " and got " << score << std::endl;
	}
	return score_total / grid_sets_.size();
}


core::Real Transform::convert_to_full_pose(core::pose::Pose & pose, core::conformation::UltraLightResidue & residue, core::Size chain_id)
{

	core::Real best_score = 10000;
	std::string best_pose_hash = "";

	//Main model should be first grid set
	if ( use_main_model_ ) {
		best_score = grid_sets_[1].second->total_score(residue);
		TR.Debug << "Using main model with score " << best_score << std::endl;
		best_pose_hash = "MAIN";
	} else {
		for ( const auto & grid_set : grid_sets_ ) {
			core::Real current_score = grid_set.second->total_score(residue);

			TR.Debug << "Looked for best ligand at " << grid_set.first << " and got " << current_score << std::endl;

			//This has a bias of order you put poses in if scores equal. Could keep separate vector of keys to shuffle
			if ( current_score <= best_score ) {
				best_score = current_score;
				best_pose_hash = grid_set.first;
			}
		}
	}

	TR.Debug << "best pose is " << best_pose_hash << std::endl;

	if ( best_pose_hash == "MAIN" ) {
		residue.update_conformation(pose.conformation());
	} else {
		core::pose::PDBInfoCOP pdb_info_original(pose.pdb_info());
		core::conformation::Residue original_residue = pose.residue(pose.conformation().chain_begin(chain_id));
		core::Size original_residue_begin = pose.conformation().chain_begin(chain_id);
		core::Size original_residue_end = pose.conformation().chain_end(chain_id);

		//Set pose to be the best pose and then add ligand

		pose = grid_set_poses_[best_pose_hash];
		pose.append_residue_by_jump(original_residue, grid_set_poses_[best_pose_hash].size(),"","",true);
		pose.pdb_info()->obsolete(false);
		pose.pdb_info()->copy(*pdb_info_original, original_residue_begin, original_residue_end,pose.size());

		residue.update_conformation(pose.conformation());
	}

	return best_score;

}

void Transform::make_multi_pose_grids(core::Vector center)
{

	qsar::scoring_grid::GridManager* grid_manager = qsar::scoring_grid::GridManager::get_instance();
	for ( const auto & grid_set_pose : grid_set_poses_ ) {

		grid_sets_.push_back(std::make_pair(grid_set_pose.first, grid_manager->get_grids( *grid_set_prototype_, grid_set_pose.second, center, transform_info_.chain)));
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
		+ XMLSchemaAttribute::attribute_w_default("use_constraints", xsct_rosetta_bool,
		"Adjust scores based on constraint file input", "false")
		+ XMLSchemaAttribute::attribute_w_default("use_main_model", xsct_rosetta_bool,
		"Use the primary input protein model regardless of scores", "false")
		+ XMLSchemaAttribute::attribute_w_default("cst_fa_file", xs_string,
		"Full atom constraint file to read constraints from", "")
		+ XMLSchemaAttribute::attribute_w_default("ensemble_proteins", xs_string,
		"File to read protein ensemble options from", "")
		+ XMLSchemaAttribute::attribute_w_default("cst_fa_weight", "non_negative_real" ,
		"Weight for full atom constraints. Default of 1.0", "1.0")
		+ XMLSchemaAttribute::attribute_w_default("initial_perturb", "non_negative_real" ,
		"Make an initial, unscored translation and rotation "
		"Translation will be selected uniformly in a sphere of the given radius (Angstrom)."
		"Triggers 360 degress rotations are triggered.", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("initial_angle_perturb", "non_negative_real", "Control the size of the rotational perturbation by intitial_perturb. (Degrees)", "0.0");

	protocols::qsar::scoring_grid::attributes_for_parse_grid_set_from_tag(attlist, "The Scoring Grid set to use with Transform scoring");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Performs a monte carlo search of the ligand binding site using precomputed scoring grids. "
		"Replaces the Translate, Rotate, and SlideTogether movers.", attlist );
}

/// @brief Provide the citation.
void
Transform::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	basic::citation_manager::CitationCollectionOP cc(
		utility::pointer::make_shared< basic::citation_manager::CitationCollection >(
		"Transform", basic::citation_manager::CitedModuleType::Mover
		)
	);

	cc->add_citation( basic::citation_manager::CitationManager::get_instance()->get_citation_by_doi( "10.1371/journal.pone.0132508" ) );

	citations.add( cc );
}

std::string TransformCreator::keyname() const {
	return Transform::mover_name();
}

protocols::moves::MoverOP
TransformCreator::create_mover() const {
	return utility::pointer::make_shared< Transform >();
}

void TransformCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	Transform::provide_xml_schema( xsd );
}

}
}
