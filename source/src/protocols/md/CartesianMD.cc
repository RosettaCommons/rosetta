// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/md/MDBase.hh
/// @brief   initialization for MD
/// @detailed
/// @author  Hahnbeom Park

#include <protocols/md/CartesianMDCreator.hh>
#include <protocols/md/CartesianMD.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/symmetry/util.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

//Optimization
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/cartesian_minimize.hh>
#include <core/optimization/types.hh>

//Mass setup
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Rosetta Scripts
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/util.hh>

//Constraints
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

//For reading native
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDB_Info.hh>

#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>

//Temporary
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/rms_util.hh>
#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

namespace protocols{
namespace md{

static basic::Tracer TR("protocols.md.Cartesian");

using namespace devel::md;
using namespace core::optimization;

static numeric::random::RandomGenerator RG( 51853218 ); //Magic number??

// creator
std::string
CartesianMDCreator::keyname() const
{
	return CartesianMDCreator::mover_name();
}

protocols::moves::MoverOP
CartesianMDCreator::create_mover() const {
	return new CartesianMD;
}

std::string
CartesianMDCreator::mover_name()
{
	return "CartesianMD";
}


// mover
CartesianMD::CartesianMD()
{
	init();
}

CartesianMD::CartesianMD( core::pose::Pose const & pose, 
													core::scoring::ScoreFunctionCOP sfxn,
													core::kinematics::MoveMapCOP movemap )
{
	set_movemap( pose, movemap );
	scorefxn_ = sfxn->clone();
	scorefxn_obj_ = scorefxn_->clone();
	init();
	if( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ){
		get_native_info( pose );
	}
}

void
CartesianMD::init()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// Default is to use Rattle
	use_rattle_ = true;
	dt_ = 0.002;

	// Access to these rather by set_reportstep & set_selectmode
	//md_report_stepsize_ = option[ md::report ]();
	//selectmode_ = option[ md::selectmode ]();
	md_report_stepsize_ = 50;
	selectmode_ = "final";

	nstep_ = 100;
	temp0_ = 300.0;

	if( option[ in::file::md_schfile ].user() ){
		scheduled_ = true;
		parse_schfile( option[ in::file::md_schfile ]() );
	} else {
		scheduled_ = false;
	}

	context_update_step_ = 10000000; // Default: Never update

	// Default
	ncyc_premin_ = 50;
  ncyc_postmin_ = 200;
	report_scorecomp_ = false;
	native_given_ = false;
	constrained_ = false;

	// Trajectory
	trj_.resize( 0 );
}

CartesianMD::~CartesianMD(){}

protocols::moves::MoverOP
CartesianMD::clone() const {
	return new CartesianMD(*this);
}

void CartesianMD::set_movemap(
								 core::pose::Pose const &,
								 core::kinematics::MoveMapCOP movemap )
{
	movemap_ = movemap->clone();
}

std::string CartesianMD::get_name() const
{
	return "CartesianMD";
}

void CartesianMD::get_native_info( pose::Pose const &pose )
{

	native_given_ = true;
	
	std::string nativepdb = basic::options::option[ basic::options::OptionKeys::in::file::native ]();
	core::chemical::ResidueTypeSetCAP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::import_pose::pose_from_pdb( native_, *rsd_set, nativepdb );
	
	// Set resmap
	std::map< Size, Size > resmap;
	for( Size ires = 1; ires <= pose.total_residue(); ++ires ){
		Size ii_pdb( pose.pdb_info()->number( ires ) );
		
		for( Size jres = 1; jres <= native_.total_residue(); ++jres ){
			Size jj_pdb( native_.pdb_info()->number( jres ) );
			if( ii_pdb == jj_pdb ){
				resmap[ires] = jres;
				break;
			}
		}
	}
	native_resmap_ = resmap;
}

void CartesianMD::do_initialize( core::pose::Pose &pose,
																 core::Real const &temp0 )
{

	// Steal utilities in CartesianMinimizer.cc
	CartesianMinimizerMap min_map;

	min_map.setup( pose, *movemap() );
	n_dof_ = min_map.ndofs();
	n_dof_temp_ = n_dof_ - 6;
	cummulative_time_ = 0.0;
	pose0_ = pose;

	// Check initial time
#ifndef WIN32
	gettimeofday(&inittime_, NULL );
#endif
	// Allocate
	xyz_.resize( n_dof() );
	vel_.resize( n_dof() );
	acc_.resize( n_dof() );
	mass_.resize( n_dof(), 0.0 );

	// Mass setup
	core::chemical::ElementSetCAP element_set
		( core::chemical::ChemicalManager::get_instance()->element_set("default") );
	core::chemical::ResidueTypeSetCAP rsdtype_set
		( core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );

	for( Size iatm = 1; iatm <= (Size)(min_map.natoms()); ++iatm ){
		id::AtomID AtomID = min_map.get_atom( iatm );
		Size resno = AtomID.rsd();
		Size atmno = AtomID.atomno();

		//std::string const element_name = pose.residue(resno).atom_type(atmno).element();
		//int const element_index = element_set->element_index(element_name);
		//mass_[iatm] = (*element_set)[element_index]->weight();
		mass_[iatm] = pose.residue_type(resno).atom(atmno).element_type()->weight();

	}

	// Velocity assignment
  initialize_velocity( temp0 );

}

void CartesianMD::apply( core::pose::Pose & pose ){
	using namespace core::optimization;

	//fpd we have to do this here since this the first time "seeing" the symm pose
	if (core::pose::symmetry::is_symmetric( pose ))
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap_ );

	// setup the map of the degrees of freedom
	const std::string minopt( "lbfgs_armijo_nonmonotone" );
	MinimizerOptions options_init( minopt, 0.0001, true, false, false );
	options_init.max_iter( ncyc_premin_ );
	MinimizerOptions options_final( minopt, 0.000001, true, false, false );
	options_final.max_iter( ncyc_postmin_ );

	const core::pose::Pose pose0( pose );

	// Set constraint
	cst_on_pose( pose );

	TR << "Reporting initial..." << std::endl;
	report_MD( pose );

	do_minimize( pose, options_init, true );

	TR << "Reporting after minimization..." << std::endl;
	report_MD( pose );

	if( report_scorecomp_ ) scorefxn_->show( TR, pose );

	// Set initial minobj_ after minimization
  // These are used to pick the best pose along trajectory (if selectmode == minobj )
	pose_minobj_ = pose;
	Emin_obj_ = 1.0e6;
	time_minobj_ = 0.0;

	// Main 
	if( !scheduled_ ){ // typical run
		do_MD( pose, nstep(), temp0(), true );

	} else { // scheduled run
		for( Size i_step = 1; i_step <= mdsch_.size(); ++i_step ){
			TR << "Changing schedule, Nstep/Temp:";
			TR << mdsch_[i_step].nstep << " " << mdsch_[i_step].temp0 << std::endl;
			if( i_step == 1 ){
				do_MD( pose, mdsch_[i_step].nstep, mdsch_[i_step].temp0, true );
			} else {
				do_MD( pose, mdsch_[i_step].nstep, mdsch_[i_step].temp0, false );
			}
		}
	}

	/// Selection for returning structure
	if( selectmode_.compare("final") == 0 ){
		// just return final pose
		TR << "Returning final structure for MD..." << std::endl;
	} else if ( selectmode_.compare("minobj") == 0 ){
		pose = pose_minobj_;
		TR << "Returning minimum objective function structure at ";
		TR << time_minobj_ << " in MD trajectory..." << std::endl;
	}

	do_minimize( pose, options_final, true );
	if( report_scorecomp_ ) scorefxn_->show( TR, pose );

	TR << "MD Done. " << std::endl;
}


void CartesianMD::do_minimize( core::pose::Pose &pose,
															 core::optimization::MinimizerOptions const &options,
															 bool const &show_energy )
{
	CartesianMinimizer minimizer;

	if( show_energy ){ 
		Real score_before = scorefxn_->score( pose ); 

		minimizer.run( pose, *movemap(), *scorefxn_, options );

		Real score_after = scorefxn_->score( pose ); 
		TR << "Energy before/after Min: " << score_before << " " << score_after << std::endl;

	} else {
		minimizer.run( pose, *movemap(), *scorefxn_, options );
	}
}

void CartesianMD::do_MD( core::pose::Pose & pose,
												 Size const &nstep,
												 Real const &temp0, 
												 bool const &initialize )
{
	if( initialize ){
		TR << "Running MD simulations for " << nstep*dt() << " ps at " << temp0 << " K... " << std::endl;
		do_initialize( pose, temp0 );
	}

	// Set dof variables

	CartesianMinimizerMap min_map;
	min_map.setup( pose, *movemap() );
	min_map.copy_dofs_from_pose( pose, xyz_ );

	// Setup RATTLE using min_map
	md::Rattle rattle( pose, min_map );
  if (use_rattle_) n_dof_temp_ = n_dof_ - 6 - rattle.ncst();

	// Set thermostat
  Thermostat thermostat( temp0, n_dof_temp_ );
 
  // Start MD integrator
	scorefxn_->setup_for_minimizing( pose, min_map );

  for ( Size istep = 1; istep <= nstep; istep++ ){
		cummulative_time_ += dt();

		// Report
		if ( istep%md_report_stepsize_ == 0 ){
			report_MD( pose );
		}

		if ( istep%context_update_step_ == 0 ){
			//TR << "update scorefunction..." << std::endl;
			scorefxn_->score( pose );
		}

		// Integrate
    VelocityVerlet_Integrator( pose, min_map, rattle );

		// Calculate/re-eval temperature
		temperature_ = thermostat.get_temperature( vel_, mass_ );

    if ( istep%thermostat.nstep_per_update() == 0 ){
      thermostat.rescale( vel_, dt(), mass_ );
			temperature_ = thermostat.get_temperature( vel_, mass_ );
    }
		kinetic_energy_ = 0.5*temperature_*n_dof()*GasConst;

  }

	min_map.copy_dofs_to_pose( pose, xyz_ );

}

void CartesianMD::VelocityVerlet_Integrator( pose::Pose &pose,
																						 CartesianMinimizerMap &min_map,
																						 md::Rattle &rattle ){

	Real dt2_2 = dt()*dt()*0.5;

	// Use previous acceleration here
	// and integrate first half of the velocity
	for( Size i_dof = 1; i_dof <= n_dof(); ++i_dof ){
		xyz_[i_dof] += vel_[i_dof]*dt() + acc_[i_dof]*dt2_2;
		vel_[i_dof] += 0.5*acc_[i_dof]*dt();
	}

	if (use_rattle_){
		rattle.run_rattle1( dt(), xyz_, vel_, mass_ );
	}

	// Reflect change in coordinates into pose
	min_map.copy_dofs_to_pose( pose, xyz_ );

	scorefxn_->score( pose );

  Multivec force;
	CartesianMultifunc f_ros( pose, min_map, *scorefxn_, false, false );
	f_ros.dfunc( xyz_, force );

	// Here, convert force into acceleration
	// and integrate remaining half of velocity
	for( Size i_dof = 1; i_dof <= n_dof(); ++i_dof ){
		Size i_atm = (i_dof+2)/3;
		// pass Virtual atoms
		if ( mass_[i_atm] < 1e-3 ) continue;

		acc_[i_dof] = -MDForceFactor*force[i_dof]/mass_[i_atm];
		vel_[i_dof] += 0.5*acc_[i_dof]*dt();
	}

	if (use_rattle_){
		rattle.run_rattle2( dt(), xyz_, vel_, mass_ );
	}

  //Stop rotation and translation
  //if((step%nrottrans)==0){
  //  stop_rot_trans(xyz, vel, size);
  //  thermostat_.get_temperature();
  //}
}

void CartesianMD::initialize_velocity( core::Real const &temperature )
{

  for( core::Size i_dof=1; i_dof<=n_dof(); i_dof++ ){
		core::Size i_atm = (i_dof+2)/3;
		// pass Virtual atoms
		if ( mass_[i_atm] < 1e-3 ) continue;

    vel_[i_dof] = sqrt(2.0*temperature*Boltzmann/mass_[i_atm])*RG.gaussian(); 
		acc_[i_dof] = 0.0;
  }
}

void CartesianMD::report_MD( core::pose::Pose &pose )
{
	core::Real const rmsd( core::scoring::CA_rmsd( pose0_, pose ));

	core::scoring::constraints::ConstraintSetCOP cstset( pose.constraint_set() );
	//TR << "Is there CST? " << std::endl;
	//cstset->show_numbers( TR );

	timeval currtime;
#ifndef WIN32
	gettimeofday(&currtime, NULL );
#endif
	Real elapsedTime = (currtime.tv_sec - inittime_.tv_sec) * 1000.0;
	elapsedTime += (currtime.tv_usec - inittime_.tv_usec) / 1000.0;
	elapsedTime /= 60000.0; // in minute

	core::Real Epot( scorefxn_->score( pose ) );
	//scorefxn_->show( pose );

	TR << "Time/E/Temp/RMSD/Elapsed(Min): " << std::setw(8) << cummulative_time();
	TR << " " << std::setw(12) << std::setprecision(6) << Epot;
	TR << " " << std::setw(6) << std::setprecision(4) << temperature_;
	TR << " " << std::setw(8) << std::setprecision(4) << rmsd;
	TR << " " << std::setw(8) << std::setprecision(4) << elapsedTime;
	if( native_given_ ){
		core::Real const rmsd_native( core::scoring::CA_rmsd( pose, native_, native_resmap_ ) );
		core::Real const gdtmm_native( core::scoring::CA_gdtmm( pose, native_, native_resmap_ ));
		TR << ", RMSD/GDTtoNative ";
		TR << std::setw(8) << std::setprecision(4) << rmsd_native;
		TR << std::setw(8) << std::setprecision(4) << gdtmm_native;
	}

	//core::Real Eobj( pose.energies().total_energies()[ core::scoring::elec_dens_fast ] );
	core::Real Eobj( scorefxn_obj_->score( pose ) );
	TR << " " << Eobj << " " << Emin_obj_;
	TR << std::endl;

	if( cummulative_time() > 0.1 && // Truncate initial 1ps to remove minimization memory 
		 selectmode_.compare("minobj") == 0 && Eobj < Emin_obj_ ){
		pose_minobj_ = pose;
		Emin_obj_ = scorefxn_obj_->score( pose );
		time_minobj_ = cummulative_time();
		TR << "Updating minimum objective score value / pose at time " << cummulative_time() << std::endl;
	}

	// Store trj
	// Can't we just make min_map as local variable?
	if( store_trj() ){
		CartesianMinimizerMap min_map;
		Multivec xyz;

		min_map.setup( pose, *movemap() );
		xyz.resize( min_map.ndofs() );
		min_map.copy_dofs_from_pose( pose, xyz );
		trj_.push_back( xyz );
	}

	/*
	std::stringstream ss;
	Size timesize( (int)(cummulative_time()*1000.0) );
	ss << "md" << timesize << ".pdb";
	pose.dump_pdb( ss.str() );
	*/

	/*
	//Check
	for ( Size i_atm = 1; i_atm <= n_dof()/3; ++i_atm ){
		Real vv = std::sqrt(vel_[3*i_atm-2]*vel_[3*i_atm-2] + vel_[3*i_atm-1]*vel_[3*i_atm-1] + vel_[3*i_atm]*vel_[3*i_atm]);
		Real av = std::sqrt(acc_[3*i_atm-2]*vel_[3*i_atm-2] + acc_[3*i_atm-1]*acc_[3*i_atm-1] + acc_[3*i_atm]*acc_[3*i_atm]);
		id::AtomID AtomID = min_map.get_atom( i_atm );
		Size resno = AtomID.rsd();
		Size atmno = AtomID.atomno();

		std::cout << "I/x/y/z/v/a: " << std::setw(6) << i_atm;
		std::cout << " " << std::setw(5) << pose.residue(resno).name();
		std::cout << " " << std::setw(4) << pose.residue(resno).atom_name(atmno);
		std::cout << " " << std::setw(10) << xyz_[3*i_atm-2];
		std::cout << " " << std::setw(10) << xyz_[3*i_atm-1];
		std::cout << " " << std::setw(10) << xyz_[3*i_atm];
		std::cout << " " << std::setw(10) << vv;
		std::cout << " " << std::setw(10) << av;
		std::cout << std::endl;
	}
	*/
}

// pose_ref should have exactly same molecule as stored in trj
utility::vector1< pose::Pose > 
CartesianMD::dump_poses( pose::Pose const &pose_ref ) const {

	utility::vector1< pose::Pose > poses;

	// Don't know why min_map doesn't support const pose,
	// so let's just copy for now
	pose::Pose pose_tmp( pose_ref );

	// note that min_map is not const function
	CartesianMinimizerMap min_map;
	min_map.setup( pose_tmp, *movemap() );

	pose::Pose pose;
	for( Size i_trj = 1; i_trj <= trj_.size(); ++i_trj ){
		min_map.copy_dofs_to_pose( pose, trj_[i_trj] );
		poses.push_back( pose );
	}

	return poses;
}

void  CartesianMD::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & pose )
{
	// Other options
	parse_opts( tag, data, pose );

	// Movemap
	parse_movemap( tag, data, pose );
}

void CartesianMD::parse_opts(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	//Filters_map const &,
	//protocols::moves::Movers_map const &,
	Pose const & pose )
{
	using namespace core;
	using namespace scoring;

	std::string const scorefxn_name( tag->getOption< std::string >( "scorefxn" ) );
	scorefxn_ = data.get< ScoreFunction * >( "scorefxns", scorefxn_name )->clone();

	std::string const scoreobj_name( tag->getOption< std::string >( "scorefxn_obj","" ) );

	if( scoreobj_name.compare("") == 0 ){
		scorefxn_obj_ = scorefxn_->clone();
	} else {
		scorefxn_obj_ = data.get< ScoreFunction * >( "scorefxns", scoreobj_name )->clone();
	}

	use_rattle_ = tag->getOption< bool >( "rattle", true );
	if( use_rattle_ ) dt_ = 0.002;

	nstep_ = tag->getOption< core::Size >( "nstep", 100 );
	temp0_ = tag->getOption<core::Real>("temp", 300.0);
	scheduled_ = false;

	ncyc_premin_ = tag->getOption< core::Size >( "premin", 50 );
	ncyc_postmin_ = tag->getOption< core::Size >( "postmin", 200 );

 	md_report_stepsize_ = tag->getOption< core::Size >( "report", 100 );
 	report_scorecomp_ = tag->getOption< bool >( "report_scorecomp", false );
 	selectmode_ = tag->getOption< std::string >( "selectmode", "final" );

	// Use parsed schedule file - this will overload nstep, temperature, etc. defined above
	std::string const schfile( tag->getOption< std::string >( "schfile","" ) );
	if( schfile.compare("") != 0 ){
		parse_schfile( schfile );
		scheduled_ = true;
	}

	//std::cout << "#######################scoreobj_/selectmode ?? " << scoreobj_name << " " << selectmode_ << std::endl;

	if( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ){
		get_native_info( pose );
	}
}

void CartesianMD::parse_movemap(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	//Filters_map const &,
	//protocols::moves::Movers_map const &,
	Pose const & pose )
{
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	bool const chi( tag->getOption< bool >( "chi", true ) ), bb( tag->getOption< bool >( "bb", true ) );
	movemap->set_chi( chi );
	movemap->set_bb( bb );
	set_movemap( pose, movemap );

	protocols::rosetta_scripts::parse_movemap( tag, pose, movemap_, data, false );
}

}
}
