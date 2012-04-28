// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/ddG.cc
/// @brief implementation of the ddG class for computing interface delta dGs
/// @author Sarel Fleishman (sarelf@u.washington.edu)


#include <core/types.hh>
//#include <protocols/moves/ResidueMover.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <numeric/xyzVector.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/simple_moves/ddGCreator.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/jd2/JobDistributor.hh>


#include <protocols/moves/DataMap.hh>

#include <ObjexxFCL/format.hh>

// C++ headers
#include <map>
#include <algorithm>

#include <basic/Tracer.hh>

#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>
#include <boost/functional/hash.hpp>
namespace protocols {
namespace simple_moves {

using core::conformation::symmetry::SymmetricConformation;
using core::conformation::symmetry::SymmetryInfoCOP;

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.protein_interface_design.movers.ddG" );

moves::MoverOP ddGCreator::create_mover() const
{
	return new ddG;
}

std::string ddGCreator::mover_name()
{
	return "ddG";
}

std::string ddGCreator::keyname() const
{
	return ddGCreator::mover_name();
}

using namespace core;
using namespace protocols::simple_moves;
using namespace core::scoring;

ddG::ddG() :
		simple_moves::DesignRepackMover(ddGCreator::mover_name()),
		bound_total_energy_(0.0),
		unbound_total_energy_(0.0),
		repeats_(0),
		rb_jump_(0),
		symmetry_(false),
		per_residue_ddg_(false),
		repack_(false),
		relax_mover_( NULL )
{
	bound_energies_.clear();
	unbound_energies_.clear();
	bound_per_residue_energies_.clear();
	unbound_per_residue_energies_.clear();
}

ddG::ddG( core::scoring::ScoreFunctionCOP scorefxn_in, core::Size const jump/*=1*/, bool const symmetry /*=false*/ ) :
		simple_moves::DesignRepackMover(ddGCreator::mover_name()),
		bound_total_energy_(0.0),
		unbound_total_energy_(0.0),
		repeats_(0),
		rb_jump_(0),
		symmetry_(false),
		per_residue_ddg_(false),
		repack_(false),
		relax_mover_( NULL )
{
	scorefxn_ = new core::scoring::ScoreFunction( *scorefxn_in );
	rb_jump_ = jump;
	symmetry_ = symmetry;
	per_residue_ddg_ = false;
	repack_ = true;
	repeats_ = 1;

	bound_energies_.clear();
	unbound_energies_.clear();
	bound_per_residue_energies_.clear();
	unbound_per_residue_energies_.clear();
}

void ddG::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap  & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const& pose)
{
	rb_jump_ = tag->getOption<core::Size>("jump", 1);
	symmetry_ = tag->getOption<bool>("symmetry",0);
	per_residue_ddg_ = tag->getOption<bool>("per_residue_ddg",0);
	repack_ = tag->getOption<bool>("repack",0);
	repeats_ = tag->getOption<Size>("repeats",1);

	if(tag->hasOption("chains") && symmetry_)
	{
		utility_exit_with_message("you cannot specify multiple chains and use symmetry mode in the ddG mover at the same time right now. Sorry");
	}

	if( ( tag->hasOption("chain_num") || tag->hasOption("chain_name") ) && tag->hasOption("jump"))
	{
		utility_exit_with_message("you can specify either chains or jump in the ddG mover, but not both");
	}

	if(tag->hasOption("chain_num"))
	{
		chain_ids_ = utility::string_split(tag->getOption<std::string>("chain_num"),',',core::Size());
	}

	if(tag->hasOption("chain_name"))
	{
		utility::vector1<std::string> chain_names = utility::string_split(tag->getOption<std::string>("chain_name"),',',std::string());
		for(utility::vector1<std::string>::iterator chain_name_it = chain_names.begin(); chain_name_it != chain_names.end(); ++chain_name_it)
		{
			chain_ids_.push_back(core::pose::get_chain_id_from_chain(*chain_name_it,pose));
		}
	}

	if(std::find(chain_ids_.begin(),chain_ids_.end(),1) != chain_ids_.end())
	{
		utility_exit_with_message("You can't move the first chain.  Moving chain 2 is the same as moving chain 1, so do that instead.");
	}

	std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn", "score12" ) );
	scorefxn_ = new ScoreFunction( *(data.get< ScoreFunction * >( "scorefxns", scorefxn_name )) );
}

ddG::~ddG() {}

void ddG::apply(Pose & pose)
{
	if(per_residue_ddg_)
	{
		core::scoring::methods::EnergyMethodOptionsOP energy_options(new core::scoring::methods::EnergyMethodOptions(scorefxn_->energy_method_options()));
		energy_options->hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn_->set_energy_method_options(*energy_options);

	}

	Real average_ddg = 0.0;
	std::map<Size, Real> average_per_residue_ddgs;

	for(Size repeat = 1; repeat <= repeats_; ++repeat)
	{
		calculate(pose);
		average_ddg += sum_ddG();
		report_ddG(TR);
		if(per_residue_ddg_)
		{
			for(core::Size i = 1; i <= pose.n_residue();++i)
			{
				core::Real bound_energy = bound_per_residue_energies_[i];
				core::Real unbound_energy = unbound_per_residue_energies_[i];
				core::Real residue_ddg = bound_energy - unbound_energy;
				if(average_per_residue_ddgs.find(i) == average_per_residue_ddgs.end())
				{
					average_per_residue_ddgs[i] = residue_ddg;
				}else
				{
					average_per_residue_ddgs[i] += residue_ddg;
				}
			}
		}
	}

	average_ddg /= repeats_;
	for(std::map<Size,Real>::iterator avg_it = average_per_residue_ddgs.begin(); avg_it != average_per_residue_ddgs.end();++avg_it)
	{
		avg_it->second /= repeats_;
	}

	jd2::JobOP job(jd2::JobDistributor::get_instance()->current_job());
	job->add_string_real_pair("ddg",average_ddg);
	if (per_residue_ddg_)
	{
		for (core::Size i = 1; i <= pose.n_residue(); ++i) {
			std::string residue_string(utility::to_string<core::Size>(i));
			job->add_string_real_pair("residue_ddg_"+residue_string,average_per_residue_ddgs[i]);
		}
	}

}

/// @details a private function for storing energy values
void
ddG::fill_energy_vector( pose::Pose const & pose, std::map< ScoreType, core::Real > & energy_map )
{
	energy_map.clear();
	using namespace core::scoring;
	for ( int i=1; i<= n_score_types; ++i ) {
		if ( (*scorefxn_)[ ScoreType(i) ] != 0.0 && ScoreType(i) != pro_close ) {
			energy_map.insert( std::make_pair( ScoreType( i ), (*scorefxn_)[ ScoreType(i)] * pose.energies().total_energies()[ ScoreType( i ) ] ) );
		}
	}
}

void ddG::fill_per_residue_energy_vector(pose::Pose const & pose, std::map<Size, Real> & energy_map)
{
	energy_map.clear();
	for(core::Size resid = 1; resid <= pose.total_residue(); ++resid)
	{
		core::Real energy = 0.0;
		for (int st = 1; st <= n_score_types; ++st)
		{
			if( (*scorefxn_)[ScoreType(st)] != 0.0 && ScoreType(st) != pro_close)
			{
				energy += (*scorefxn_)[ScoreType(st)] * pose.energies().residue_total_energies(resid)[ScoreType(st)];
			}
		}
		energy_map.insert(std::make_pair(resid,energy));
	}
}

/// @details output the ddG values
void
ddG::report_ddG( std::ostream & out ) const
{

	using ObjexxFCL::fmt::F;
	using ObjexxFCL::fmt::LJ;
	out << "-----------------------------------------\n";
	out << " Scores                       Wghtd.Score\n";
	out << "-----------------------------------------\n";
	std::map< ScoreType, Real >::const_iterator unbound_it=unbound_energies_.begin();
	for( std::map< ScoreType, Real >::const_iterator bound_it=bound_energies_.begin();
			 bound_it!=bound_energies_.end();
			 ++bound_it ) {
			 if( std::abs( unbound_it->second ) > 0.001 || std::abs( bound_it->second ) > 0.001 )
				out << ' ' << LJ( 24, bound_it->first ) << ' ' << F( 9,3, bound_it->second - unbound_it->second )<<'\n';
		++unbound_it;
	}
	out << "-----------------------------------------\n";
	out << "Sum ddg: "<< sum_ddG()<<std::endl;
}

/// @details returns the total ddG
Real
ddG::sum_ddG() const
{
	Real sum_energy(0.0);

	std::map< ScoreType, Real >::const_iterator unbound_it=unbound_energies_.begin();
	for( std::map< ScoreType, Real >::const_iterator bound_it=bound_energies_.begin();
			 bound_it!=bound_energies_.end();
			 ++bound_it ) {
		sum_energy += bound_it->second - unbound_it->second;
		++unbound_it;
	}
	return sum_energy;
}


/// @details compute the energy of the repacked complex in the bound and unbound states
void
ddG::calculate( pose::Pose const & pose_in )
{
	using namespace pack;
	using namespace protocols::moves;

	if ( core::pose::symmetry::is_symmetric( pose_in ) ) { //JBB 120423 changed from if (symmetry_)
		symm_ddG( pose_in );
		return;
	}

	if(!repack_)
	{
		no_repack_ddG(pose_in);
		return;
	}

	repack_partner1_ = true;
	repack_partner2_ = true;
	design_partner1_ = false;
	design_partner2_ = false;
	pose::Pose pose = pose_in;

	//setup_packer_and_movemap( pose );
	task_ = core::pack::task::TaskFactory::create_packer_task( pose );
	task_->initialize_from_command_line().or_include_current( true );
	core::pack::task::operation::RestrictToRepacking rpk;
	rpk.apply( pose, *task_ );
	core::pack::task::operation::NoRepackDisulfides nodisulf;
	nodisulf.apply( pose, *task_ );
	if(chain_ids_.size() > 0 )
	{
		//We want to translate each chain the same direction, though it doesnt matter much which one
		core::Size first_jump = core::pose::get_jump_id_from_chain_id(chain_ids_[1],pose);
		protocols::toolbox::task_operations::RestrictToInterface rti( first_jump, 8.0 /*interface_distance_cutoff_*/ );
		if(chain_ids_.size() > 1)
		{
			for(core::Size chain_index = 2; chain_index <= chain_ids_.size();++chain_index)
			{
				core::Size current_jump = core::pose::get_jump_id_from_chain_id(chain_ids_[chain_index],pose);
				rti.add_jump(current_jump);
			}
		}
		rti.apply(pose,*task_);
	}
	else
	{
		protocols::toolbox::task_operations::RestrictToInterface rti( rb_jump_, 8.0 /*interface_distance_cutoff_*/ );
		rti.apply( pose, *task_ );
	}

	pack::pack_rotamers( pose, *scorefxn_, task_ );
	(*scorefxn_)( pose );
	fill_energy_vector( pose, bound_energies_ );
	if(per_residue_ddg_)
	{
		fill_per_residue_energy_vector(pose, bound_per_residue_energies_);
	}

	if(chain_ids_.size() > 0)
	{
		//We want to translate each chain the same direction, though it doesnt matter much which one
		Vector translation_axis(1,0,0);
		for(utility::vector1<core::Size>::const_iterator chain_it = chain_ids_.begin(); chain_it != chain_ids_.end();++chain_it)
		{
			core::Size current_chain_id = *chain_it;
			core::Size current_jump_id = core::pose::get_jump_id_from_chain_id(current_chain_id,pose);
			rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( pose, current_jump_id) );
			translate->step_size( 1000.0 );
			translate->trans_axis(translation_axis);
			translate->apply( pose );
		}
	}else
	{
		rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( pose, rb_jump_ ) );
		translate->step_size( 1000.0 );
		translate->apply( pose );
	}

	pack::pack_rotamers( pose, *scorefxn_, task_ );
	if( relax_mover() )
		relax_mover()->apply( pose );
	(*scorefxn_)( pose );
	fill_energy_vector( pose, unbound_energies_ );
	if(per_residue_ddg_)
	{
		fill_per_residue_energy_vector(pose, unbound_per_residue_energies_);
	}
}

// @details compute the energy of a repacked symmetrical complex in bound and unbound states
void
ddG::symm_ddG( pose::Pose const & pose_in )
{
	using namespace pack;
	using namespace protocols::moves;


	pose::Pose pose = pose_in;

	assert( core::pose::symmetry::is_symmetric( pose ));
  SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

  std::map< Size, core::conformation::symmetry::SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	// convert to symetric scorefunction
	scorefxn_ = new scoring::symmetry::SymmetricScoreFunction( scorefxn_ );

	//setup_packer_and_movemap( pose );
	task_ = core::pack::task::TaskFactory::create_packer_task( pose );
	task_->initialize_from_command_line().or_include_current( true );
	core::pack::task::operation::RestrictToRepacking rpk;
	rpk.apply( pose, *task_ );
	core::pack::task::operation::NoRepackDisulfides nodisulf;
	nodisulf.apply( pose, *task_ );
	protocols::toolbox::task_operations::RestrictToInterface rti( rb_jump_, 8.0 /*interface_distance_cutoff_*/ );
	rti.apply( pose, *task_ );

	pack::symmetric_pack_rotamers( pose, *scorefxn_, task_ );
	(*scorefxn_)( pose );
	fill_energy_vector( pose, bound_energies_ );
	if(per_residue_ddg_)
	{
		fill_per_residue_energy_vector(pose, bound_per_residue_energies_);
	}
	rigid::RigidBodyDofSeqTransMoverOP translate( new rigid::RigidBodyDofSeqTransMover( dofs ) );
	translate->step_size( 1000.0 );
	translate->apply( pose );
	pack::symmetric_pack_rotamers( pose, *scorefxn_, task_ );
	if( relax_mover() )
		relax_mover()->apply( pose );
	(*scorefxn_)( pose );
	fill_energy_vector( pose, unbound_energies_ );
	if(per_residue_ddg_)
	{
		fill_per_residue_energy_vector(pose, unbound_per_residue_energies_);
	}
}

void
ddG::no_repack_ddG(Pose const & pose_in)
{
	pose::Pose pose = pose_in;

	(*scorefxn_)( pose );
	fill_energy_vector( pose, bound_energies_ );
	if(per_residue_ddg_)
	{
		fill_per_residue_energy_vector(pose, bound_per_residue_energies_);
	}

	if(chain_ids_.size()  > 0 )
	{
		//We want to translate each chain the same direction, though it doesnt matter much which one
		Vector translation_axis(1,0,0);

		for(utility::vector1<core::Size>::const_iterator chain_it = chain_ids_.begin(); chain_it != chain_ids_.end();++chain_it)
		{
			core::Size current_chain_id = *chain_it;
			core::Size current_jump_id = core::pose::get_jump_id_from_chain_id(current_chain_id,pose);
			rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( pose, current_jump_id) );
			translate->trans_axis(translation_axis);
			translate->step_size( 1000.0 );
			translate->apply( pose );
		}

	}else
	{
		rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( pose,rb_jump_ ) );
		translate->step_size( 1000.0 );
		translate->apply( pose );
	}

	if( relax_mover() )
		relax_mover()->apply( pose );
	(*scorefxn_)( pose );
	fill_energy_vector( pose,unbound_energies_ );
	if(per_residue_ddg_)
	{
		fill_per_residue_energy_vector(pose, unbound_per_residue_energies_);
	}
}

std::string
ddG::get_name() const {
	return "ddG";
}

protocols::moves::MoverOP
ddG::clone() const {
	return (protocols::moves::MoverOP) new ddG( *this );
}

} //movers
} //protocols
