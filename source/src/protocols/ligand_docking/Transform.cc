// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/ligand_docking/Transform.cc
/// @author Sam DeLuca


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

#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/format.hh>
#include <utility/io/ozstream.hh>

namespace protocols {
namespace ligand_docking {

static numeric::random::RandomGenerator RG(23459);

static basic::Tracer transform_tracer("protocols.ligand_docking.Transform");

std::string TransformCreator::keyname() const
{
	return TransformCreator::mover_name();
}

protocols::moves::MoverOP TransformCreator::create_mover() const
{
	return new Transform;
}

std::string TransformCreator::mover_name()
{
	return "Transform";
}

Transform::Transform(): Mover("Transform"), transform_info_(),optimize_until_score_is_negative_(false), output_sampled_space_(false),check_rmsd_(false),initial_perturb_(0.0)
{

}


Transform::Transform(Transform const & other) :
	Mover(other),
	transform_info_(other.transform_info_),
	ligand_conformers_(other.ligand_conformers_),
	optimize_until_score_is_negative_(other.optimize_until_score_is_negative_),
	output_sampled_space_(other.output_sampled_space_),
	check_rmsd_(other.check_rmsd_),
	sampled_space_file_(other.sampled_space_file_),
	initial_perturb_(other.initial_perturb_)
{}

Transform::Transform(
	std::string const & chain,
	core::Real const & box_size,
	core::Real const & move_distance,
	core::Real const & angle,
	core::Size const & cycles,
	core::Real const & temp
) : Mover("Transform"), transform_info_(),optimize_until_score_is_negative_(false),output_sampled_space_(false),check_rmsd_(false)
{
	transform_info_.chain = chain;
	transform_info_.box_size = box_size;
	transform_info_.move_distance = move_distance;
	transform_info_.angle = angle;
	transform_info_.cycles = cycles;
	transform_info_.temperature = temp;
}

Transform::~Transform()
{
	//
}

protocols::moves::MoverOP Transform::clone() const
{
	return new Transform (*this);
}

protocols::moves::MoverOP Transform::fresh_instance() const
{
	return new Transform;
}

std::string Transform::get_name() const
{
	return "Transform";
}

void Transform::parse_my_tag
(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & /*data_map*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "Transform" )
	{
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires chain tag");
	if ( ! tag->hasOption("move_distance") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires move_distance tag");
	if (! tag->hasOption("box_size") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires box_size tag");
	if ( ! tag->hasOption("angle") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires angle tag");
	if ( ! tag->hasOption("cycles") ) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires cycles tag");
	if (!tag->hasOption("temperature")) throw utility::excn::EXCN_RosettaScriptsOption("'Transform' mover requires temperature tag");

	transform_info_.chain = tag->getOption<std::string>("chain");
	transform_info_.move_distance = tag->getOption<core::Real>("move_distance");
	transform_info_.box_size = tag->getOption<core::Real>("box_size");
	transform_info_.angle = tag->getOption<core::Real>("angle");
	transform_info_.cycles = tag->getOption<core::Size>("cycles");
	transform_info_.temperature = tag->getOption<core::Real>("temperature");
	transform_info_.repeats = tag->getOption<core::Size>("repeats",1);
	optimize_until_score_is_negative_ = tag->getOption<bool>("optimize_until_score_is_negative",false);
	initial_perturb_ = tag->getOption<core::Real>("initial_perturb",0.0);
    

	if(tag->hasOption("rmsd"))
	{
		check_rmsd_ = true;
		transform_info_.rmsd = tag->getOption<core::Real>("rmsd");
	}
	
	if(tag->hasOption("sampled_space_file"))
	{
		output_sampled_space_ = true;
		sampled_space_file_ = tag->getOption<std::string>("sampled_space_file");
	}

}

void Transform::apply(core::pose::Pose & pose)
{
	qsar::scoring_grid::GridManager* grid_manager(qsar::scoring_grid::GridManager::get_instance());

	assert(transform_info_.chain.size() == 1);
	transform_info_.chain_id = core::pose::get_chain_id_from_chain(transform_info_.chain, pose);
	transform_info_.jump_id = core::pose::get_jump_id_from_chain_id(transform_info_.chain_id, pose);
	core::Size const begin(pose.conformation().chain_begin(transform_info_.chain_id));
	core::Vector const center(protocols::geometry::downstream_centroid_by_jump(pose, transform_info_.jump_id));
	assert(grid_manager != 0); //something has gone hopelessly wrong if this triggers

	core::conformation::Residue original_residue = pose.residue(begin);
	core::chemical::ResidueType residue_type = pose.residue_type(begin);

	
	grid_manager->initialize_all_grids(center);
	grid_manager->update_grids(pose,center);

	core::Real last_score(10000.0);

    core::Real best_score(last_score);
    core::pose::Pose best_pose(pose);
    core::pose::Pose starting_pose(pose);
    core::conformation::UltraLightResidue best_ligand(&pose.residue(begin));
	core::conformation::UltraLightResidue original_ligand(best_ligand);

	core::Real temperature = transform_info_.temperature;
	core::Vector original_center(original_residue.xyz(original_residue.nbr_atom()));


	rotamers_for_trials(pose,begin,ligand_conformers_);
	transform_tracer << "Considering " << ligand_conformers_.size() << " conformers during sampling" << std::endl;
	core::Size accepted_moves = 0;
	core::Size rejected_moves = 0;


	utility::io::ozstream sampled_space;
	if(output_sampled_space_)
	{
		sampled_space.open(sampled_space_file_);
	}



	for(core::Size repeat = 1; repeat <= transform_info_.repeats; ++repeat)
	{
		pose = starting_pose;
		core::Size cycle = 1;
		bool not_converged = true;
		core::conformation::UltraLightResidue ligand_residue(&pose.residue(begin));
        
		//For benchmarking purposes it is sometimes desirable to translate the ligand
		//away from the starting point before beginning a trajectory.
		if(initial_perturb_ > 0.0)
		{
			translate_ligand(ligand_residue,initial_perturb_ );
		}
		last_score = grid_manager->total_score(ligand_residue);
		core::conformation::UltraLightResidue last_accepted_ligand_residue(ligand_residue);
		while(not_converged)
		{


			if(optimize_until_score_is_negative_)
			{
				if(cycle >= transform_info_.cycles && last_score <= 0.0)
				{
					not_converged= false;
				}
			}else
			{
				if(cycle >= transform_info_.cycles)
				{
					not_converged= false;
				}
			}

			cycle++;

			//during each move either move the ligand or try a new conformer (if there is more than one conformer)
			if(ligand_conformers_.size() > 1)
			{
				if(RG.uniform() >= 0.5)
				{
					transform_ligand(ligand_residue);
				}else
				{
					change_conformer(ligand_residue);
				}
			}else
			{
				transform_ligand(ligand_residue);
			}
			//The score is meaningless if any atoms are outside of the grid
			if(!grid_manager->is_in_grid(ligand_residue)) //Reject the pose
			{
				ligand_residue = last_accepted_ligand_residue;
				rejected_moves++;
				//transform_tracer.Trace << "probability: " << probability << " rejected (out of grid)"<<std::endl;
				continue;
			}

			if(output_sampled_space_)
			{
				dump_conformer(ligand_residue,sampled_space);
			}

			core::Real current_score = grid_manager->total_score(ligand_residue);
			core::Real const boltz_factor((last_score-current_score)/temperature);
			core::Real const probability = std::exp( boltz_factor ) ;
			core::Vector new_center(ligand_residue.center());

			if(new_center.distance(original_center) > transform_info_.box_size) //Reject the new pose
			{
				ligand_residue = last_accepted_ligand_residue;
				rejected_moves++;

			}else if(check_rmsd_ && !check_rmsd(original_ligand, ligand_residue)) //reject the new pose
			{
				ligand_residue = last_accepted_ligand_residue;
				rejected_moves++;
			}else if(probability < 1 && RG.uniform() >= probability)  //reject the new pose
			{
				ligand_residue = last_accepted_ligand_residue;
				rejected_moves++;

			}else if(probability < 1)  // Accept the new pose
			{
				last_score = current_score;
				last_accepted_ligand_residue = ligand_residue;
				accepted_moves++;

			}else  //Accept the new pose
			{
				last_score = current_score;
				last_accepted_ligand_residue = ligand_residue;
				accepted_moves++;

			}
			if(last_score <= best_score)
			{
				best_score = last_score;
				best_ligand = last_accepted_ligand_residue;
				//transform_tracer << "accepting new pose" << std::endl;
			}else
			{
				//transform_tracer << "not accepting new pose" << std::endl;
			}
			
		}
		core::Real accept_ratio =(core::Real)accepted_moves/((core::Real)accepted_moves
		+(core::Real)rejected_moves);
		transform_tracer <<"percent acceptance: "<< accepted_moves << " " << accept_ratio<<" " << rejected_moves <<std::endl;
		
		jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("Transform_accept_ratio", accept_ratio);
		best_ligand.update_conformation(best_pose.conformation());
	}

	if(output_sampled_space_)
	{
		sampled_space.close();
	}
	pose = best_pose;
	
	transform_tracer << "Accepted pose with grid score: " << best_score << std::endl;
	
}

    
void Transform::translate_ligand(core::conformation::UltraLightResidue & residue, core::Real distance)
{
	core::Vector translation(
		distance*RG.uniform(),
		distance*RG.uniform(),
		distance*RG.uniform());
    
	
	residue.slide(translation);
}
void Transform::transform_ligand(core::conformation::UltraLightResidue & residue)
{
	if(transform_info_.angle ==0 && transform_info_.move_distance == 0)
	{
		transform_tracer <<"WARNING: angle and distance are both 0.  Transform will do nothing" <<std::endl;
		return;
	}

	core::Vector translation(
		transform_info_.move_distance*RG.gaussian(),
		transform_info_.move_distance*RG.gaussian(),
		transform_info_.move_distance*RG.gaussian());

	numeric::xyzMatrix<core::Real> rotation(
		numeric::z_rotation_matrix_degrees( transform_info_.angle*RG.gaussian() ) * (
			numeric::y_rotation_matrix_degrees( transform_info_.angle*RG.gaussian() ) *
			numeric::x_rotation_matrix_degrees( transform_info_.angle*RG.gaussian() ) ));

	residue.transform(rotation,translation);
}

void Transform::change_conformer(core::conformation::UltraLightResidue & residue)
{
	assert(ligand_conformers_.size());
	core::Size index_to_select = RG.random_range(1,ligand_conformers_.size());
	core::conformation::UltraLightResidue new_residue(ligand_conformers_[index_to_select]);
	new_residue.align_to_residue(residue);
	residue = new_residue;

}

void Transform::dump_conformer(core::conformation::UltraLightResidue & residue, utility::io::ozstream & output)
{
	using namespace ObjexxFCL::format;
	for(core::Size atom_index = 1; atom_index <= residue.natoms();++atom_index)
	{
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
	assert(start.natoms() == current.natoms());
	
	core::Real total_distance =0.0;
	for(core::Size atomno = 1; atomno <= start.natoms();++atomno)
	{
		total_distance += start[atomno].distance(current[atomno]);
	}
	
	core::Real rmsd = sqrt(total_distance/start.natoms());
	
	if(rmsd <= transform_info_.rmsd)
	{
		return true;
	}else
	{
		return false;
	}
	
}
}
}
