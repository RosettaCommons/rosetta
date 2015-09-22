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
/// @details
/// @author  Hahnbeom Park

#include <protocols/md/CartesianMDCreator.hh>
#include <protocols/md/CartesianMD.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/symmetry/util.hh>

// Constraints
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

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
#include <core/pose/PDBInfo.hh>

#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <iostream>
#include <fstream>

//Temporary
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/rms_util.hh>
#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

using namespace ObjexxFCL::format;

namespace protocols {
namespace md {

static THREAD_LOCAL basic::Tracer TR( "protocols.md.Cartesian" );

using namespace devel::md;
using namespace core::optimization;


// creator
std::string
CartesianMDCreator::keyname() const
{
	return CartesianMDCreator::mover_name();
}

protocols::moves::MoverOP
CartesianMDCreator::create_mover() const {
	return protocols::moves::MoverOP( new CartesianMD );
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
	if ( movemap == 0 ) {
		movemap_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
		movemap_->set_jump( true ); movemap_->set_bb( true ); movemap_->set_chi( true );
		movemap_->set( core::id::THETA, true ); movemap_->set( core::id::D, true);
	} else {
		set_movemap( pose, movemap );
	}

	scorefxn_ = sfxn->clone();
	scorefxn_obj_ = scorefxn_->clone();
	init();
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		get_native_info( pose );
	}
}

CartesianMD::CartesianMD( core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &sfxn )
{
	core::kinematics::MoveMap movemap;
	movemap.set_jump( true ); movemap.set_bb( true ); movemap_->set_chi( true );
	movemap.set( core::id::THETA, true ); movemap.set( core::id::D, true );

	scorefxn_ = sfxn.clone();
	scorefxn_obj_ = sfxn.clone();
	init();
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		get_native_info( pose );
	}
}

CartesianMD::CartesianMD( core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &sfxn,
	core::kinematics::MoveMap const &movemap )
{
	set_movemap( pose, movemap.clone() );

	scorefxn_ = sfxn.clone();
	scorefxn_obj_ = sfxn.clone();
	init();
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
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
	md_report_stepsize_ = 500; // every 1 ps
	md_energy_report_stepsize_ = 50; // every 0.1 ps
	md_rsr_update_stepsize_ = 50; // every 0.1 ps
	selectmode_ = "final";

	nstep_ = 100;
	temp0_ = 300.0;

	if ( option[ in::file::md_schfile ].user() ) {
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
	store_trj_ = false;
	trj_.resize( 0 );
	report_as_silent_ = false;
	silentname_ = "";
	trj_score_only_ = true;

	// Adaptive restraint
	rsrfilename_ = "";
	write_dynamic_rsr_ = false;
	ref_xyz_.resize( 0 );
	trj_scratch_.resize( 0 );
	Kappa_ = 0.1;
	Gamma_ = 0.0;
}

CartesianMD::~CartesianMD(){}

protocols::moves::MoverOP
CartesianMD::clone() const {
	return protocols::moves::MoverOP( new CartesianMD(*this) );
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
	core::chemical::ResidueTypeSetCOP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::import_pose::pose_from_pdb( native_, *rsd_set, nativepdb );

	// Set resmap
	std::map< Size, Size > resmap;
	for ( Size ires = 1; ires <= pose.total_residue(); ++ires ) {
		if ( !pose.residue( ires ).is_protein() ) continue;
		Size ii_pdb( pose.pdb_info()->number( ires ) );

		for ( Size jres = 1; jres <= native_.total_residue(); ++jres ) {
			if ( !native_.residue( jres ).is_protein() ) continue;
			Size jj_pdb( native_.pdb_info()->number( jres ) );
			if ( ii_pdb == jj_pdb ) {
				resmap[ires] = jres;
				break;
			}
		}
	}
	native_resmap_ = resmap;
}

void CartesianMD::do_initialize( core::pose::Pose &pose )
{
	TR << "Initialize MD simulations" << std::endl;

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
	ref_xyz_.resize( n_dof() );

	// for adaptive rsr
	min_map.copy_dofs_from_pose( pose, ref_xyz_ );
	prv_eqxyz_ = ref_xyz_;

	// Mass setup
	core::chemical::ElementSetCAP element_set
		( core::chemical::ChemicalManager::get_instance()->element_set("default") );
	core::chemical::ResidueTypeSetCAP rsdtype_set
		( core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );

	for ( Size iatm = 1; iatm <= (Size)(min_map.natoms()); ++iatm ) {
		id::AtomID AtomID = min_map.get_atom( iatm );
		Size resno = AtomID.rsd();
		Size atmno = AtomID.atomno();

		//std::string const element_name = pose.residue(resno).atom_type(atmno).element();
		//int const element_index = element_set->element_index(element_name);
		//mass_[iatm] = (*element_set)[element_index]->weight();
		mass_[iatm] = pose.residue_type(resno).atom(atmno).element_type()->weight();

	}

}

void CartesianMD::use_rattle( bool const value )
{
	use_rattle_ = value;
	if ( use_rattle_ ) {
		dt_ = 0.002;
		TR << "Set Rattle on, changing dt as " << dt_ << std::endl;
	} else {
		dt_ = 0.001;
		TR << "Set Rattle off, changing dt as " << dt_ << std::endl;
	}
}

void
CartesianMD::update_restraint( pose::Pose & pose,
	CartesianMinimizerMap const &min_map )
{

	if ( Gamma_ == 0.0 || !constrained_ ) return;

	TR.Debug << "Update restraints with Kappa/Gamma = " << Kappa_ << "/" << Gamma_ << std::endl;
	Multivec curr_eqxyz = get_current_eqxyz();
	cst_on_pose_dynamic( pose, ref_xyz_, curr_eqxyz, prv_eqxyz_, min_map );

	// clear temporary trj
	trj_scratch_.resize( 0 );
}

Multivec
CartesianMD::get_current_eqxyz() const
{
	Multivec curr_eqxyz;
	curr_eqxyz.resize( n_dof() );
	core::Size const ntrj( trj_scratch_.size() );

	for ( core::Size i_dof = 1; i_dof <= n_dof(); ++i_dof ) {
		curr_eqxyz[i_dof] = 0.0;
		for ( core::Size i_trj = 1; i_trj <= ntrj; ++i_trj ) {
			curr_eqxyz[i_dof] += trj_scratch_[i_trj][i_dof];
		}

		if ( ntrj > 0 ) curr_eqxyz[i_dof] /= (core::Real)(ntrj);
	}
	return curr_eqxyz;
}

void
CartesianMD::cst_on_pose_simple( pose::Pose &pose ) const
{
	using namespace core::scoring::constraints;

	// Remove all the constraints first
	//pose.remove_constraints();  // fpd ... move this below!

	// First, add cst_file info into pose
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file ].user() ) {
		TR << "Set constraints from input file..." << std::endl;
		pose.remove_constraints();
		scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose );
	}

	// Next, set coordinate constraint if specified
	if ( constrained_ ) {
		TR << "Set constraints uniformly with stdev: " << cst_sdev_ << std::endl;
		pose.remove_constraints();
		Size nrsr( 0 );
		for ( Size i_res = 1; i_res <= pose.total_residue(); ++i_res ) {
			std::string resname = pose.residue(i_res).name();

			core::Size iatm;
			if ( pose.residue(i_res).has(" CA " ) ) {
				iatm = pose.residue( i_res ).atom_index(" CA ");
			} else if ( resname.compare("TP3") == 0 && pose.residue(i_res).has(" O  ") ) {
				iatm = pose.residue( i_res ).atom_index(" O  ");
			} else {
				continue;
			}
			id::AtomID atomID( iatm, i_res );
			core::Vector xyz = pose.xyz( atomID );
			scoring::func::FuncOP fx( new scoring::func::HarmonicFunc( 0.0, cst_sdev_ ) );
			pose.add_constraint(  ConstraintCOP( ConstraintOP(
				new CoordinateConstraint( atomID, atomID, xyz, fx )
				)));
			nrsr++;
		}
		TR << "Added " << nrsr << " coordinate constraints." << std::endl;
	}
}

void
CartesianMD::cst_on_pose_dynamic( pose::Pose &pose,
	Multivec const &ref_xyz,
	Multivec const &curr_eqxyz,
	Multivec &prv_eqxyz,
	CartesianMinimizerMap const &min_map ) const
{
	using namespace core::scoring::constraints;

	// not supporting cst_fa_file yet
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file ].user() ) return;

	// Remove all the constraints first
	pose.remove_constraints();

	//CartesianMinimizerMap min_map;
	//min_map.setup( pose, *movemap() );

	// Next, set coordinate constraint if specified
	Multivec prv_eqxyz0( prv_eqxyz );

	std::ofstream rsrfile( rsrfilename_.c_str(), std::ios_base::app );
	core::Size modality = (core::Size)(cummulative_time()*100+1)%100;
	bool write_dynamic_rsr = write_dynamic_rsr_ && ( modality <= 2);

	//TR << cummulative_time() << " " << ((core::Size)(cummulative_time()*100))%100 << " " << write_dynamic_rsr << std::endl;
	TR.Debug << "ref/curr/prv_eqxyz0? " << ref_xyz.size() << " " << curr_eqxyz.size() << " " << prv_eqxyz0.size() << std::endl;

	if ( constrained_ ) {
		if ( write_dynamic_rsr ) {
			rsrfile << "dynamic rsr: " << cummulative_time() << std::endl;
			rsrfile << "Res | X Y Z | prv_eq X Y Z | delta X Y Z" << std::endl;
		}
		core::Real rmsd( 0.0 );
		core::Real rmsd_rsr2ref( 0.0 );
		core::Real rmsd_rsr2crd( 0.0 );
		core::Size nrsr_dof( 0 );
		//for( Size i_res = 1; i_res <= pose.total_residue(); ++i_res ){
		for ( Size i_atm = 1; i_atm <= n_dof()/3; ++i_atm ) {
			id::AtomID atomID = min_map.get_atom( i_atm );
			Size resno = atomID.rsd();
			Size atmno = atomID.atomno();

			std::string resname = pose.residue(resno).name();
			std::string atmname = pose.residue(resno).atom_name(atmno);

			bool is_protein_ca = ( atmname.compare(" CA ") == 0 );
			bool is_water = ( resname.compare("TP3") == 0 && atmname.compare(" O  ") == 0 );
			if ( !(is_protein_ca || is_water ) ) continue;

			nrsr_dof++;
			TR.Debug << "on " << resno << " " << atmname << std::endl;

			core::Vector xyzmix( 0.0 );
			core::Real dist_res( 0.0 );
			core::Real dist_togo( 0.0 );
			core::Real dist_rsr( 0.0 );

			for ( core::Size i = 1; i <= 3; ++i ) {
				core::Size i_dof = (i_atm-1)*3 + i;
				xyzmix[i-1] = (1.0 - Kappa_)*prv_eqxyz[ i_dof ];
				xyzmix[i-1] += Kappa_*(Gamma_*curr_eqxyz[ i_dof ] +
					(1.0-Gamma_)*ref_xyz[ i_dof ] );

				// update prv_eq crd
				prv_eqxyz[ i_dof ] = xyzmix[i-1];
				dist_res  += (curr_eqxyz[ i_dof ] - ref_xyz[ i_dof ])  *(curr_eqxyz[ i_dof ] - ref_xyz[ i_dof ]);
				dist_togo += (curr_eqxyz[ i_dof ] - prv_eqxyz[ i_dof ])*(curr_eqxyz[ i_dof ] - prv_eqxyz[ i_dof ]);
				dist_rsr  += (prv_eqxyz[ i_dof ]  - ref_xyz[ i_dof ])  *(prv_eqxyz[ i_dof ] - ref_xyz[ i_dof ]);
			}
			rmsd += dist_res;
			rmsd_rsr2crd += dist_togo;
			rmsd_rsr2ref += dist_rsr;
			dist_res = std::sqrt(dist_res);
			//dist_togo = std::sqrt(dist_togo);

			core::Real sdev( cst_sdev_ );
			//if( dist_res ) sdev *=

			scoring::func::FuncOP fx( new scoring::func::HarmonicFunc( 0.0,  sdev ) );
			pose.add_constraint( ConstraintCOP( ConstraintOP(
				new CoordinateConstraint( atomID, atomID, xyzmix, fx )
				)));

			if ( write_dynamic_rsr ) {
				rsrfile << I(4,resno);
				rsrfile << " | " << F(8,3,xyzmix[0])  << " " << F(8,3,xyzmix[1]) << " " << F(8,3,xyzmix[2]);
				rsrfile << " | " << F(8,3,prv_eqxyz0[3*(i_atm-1)+1]);
				rsrfile << " "   << F(8,3,prv_eqxyz0[3*(i_atm-1)+2]);
				rsrfile << " "   << F(8,3,prv_eqxyz0[3*(i_atm-1)+3]);
				rsrfile << " | " << F(8,3,dist_res);
				rsrfile << " "   << F(8,3,prv_eqxyz[3*(i_atm-1)+1] - ref_xyz[3*(i_atm-1)+1]);
				rsrfile << " "   << F(8,3,prv_eqxyz[3*(i_atm-1)+2] - ref_xyz[3*(i_atm-1)+2]);
				rsrfile << " "   << F(8,3,prv_eqxyz[3*(i_atm-1)+3] - ref_xyz[3*(i_atm-1)+3]);
				rsrfile << std::endl;
			}

		} //iatm

		if ( write_dynamic_rsr ) {
			rmsd         /= (core::Real)(nrsr_dof);        rmsd = std::sqrt(rmsd);
			rmsd_rsr2ref /= (core::Real)(nrsr_dof); rmsd_rsr2ref = std::sqrt(rmsd_rsr2ref);
			rmsd_rsr2crd /= (core::Real)(nrsr_dof); rmsd_rsr2crd = std::sqrt(rmsd_rsr2crd);

			rsrfile << "At time " << cummulative_time() << ", RMSD of rsr to ref: " << F(8,3,rmsd_rsr2ref);
			rsrfile << " , crd to rsr: " << F(8,3,rmsd_rsr2crd);
			rsrfile << " , crd to ref:" << F(8,3,rmsd) << std::endl;
		}
	} //if uniform_coordinate_constrained

	rsrfile.close();

}

void CartesianMD::apply( core::pose::Pose & pose ){
	using namespace core::optimization;

	//fpd we have to do this here since this the first time "seeing" the symm pose
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap_ );
	}

	// setup the map of the degrees of freedom
	const std::string minopt( "lbfgs_armijo_nonmonotone" );
	MinimizerOptions options_init( minopt, 0.0001, true, false, false );
	options_init.max_iter( ncyc_premin_ );
	MinimizerOptions options_final( minopt, 0.000001, true, false, false );
	options_final.max_iter( ncyc_postmin_ );

	// Get initial crd info
	const core::pose::Pose pose0( pose );


	// Set constraint
	TR << "Set constraint on initial." << std::endl;
	cst_on_pose_simple( pose );

	TR << "Reporting initial..." << std::endl;
	CartesianMinimizerMap min_map;
	min_map.setup( pose, *movemap() );
	report_MD( pose, min_map, true );

	do_minimize( pose, options_init, true );

	TR << "Reporting after minimization..." << std::endl;
	report_MD( pose, min_map, true );

	if ( report_scorecomp_ ) scorefxn_->show( TR, pose );

	// Set initial minobj_ after minimization
	// These are used to pick the best pose along trajectory (if selectmode == minobj )
	pose_minobj_ = pose;
	Emin_obj_ = 1.0e6;
	time_minobj_ = 0.0;

	// Main
	if ( !scheduled_ ) { // typical run
		do_MD( pose, nstep(), temp0(), true );

	} else { // scheduled run
		for ( Size i_step = 1; i_step <= mdsch_.size(); ++i_step ) {
			std::string const runtype( mdsch_[i_step].type );
			if ( i_step == 1 ) do_initialize( pose );

			TR << "sch " << i_step << ", runtype: " << runtype << std::endl;
			if ( runtype.compare("sch") == 0 ) {
				TR << "Changing schedule, Nstep/Temp:";
				TR << mdsch_[i_step].nstep << " " << mdsch_[i_step].temp0 << std::endl;
				do_MD( pose, mdsch_[i_step].nstep, mdsch_[i_step].temp0, false );
			} else if ( runtype.compare("repack") == 0 ) {
				//
			} else if ( runtype.compare("min") == 0 ) {
				do_minimize( pose, options_init, true );

			} else if ( runtype.compare("set_weight") == 0 ) {
				scorefxn_->set_weight( mdsch_[i_step].scoretype, mdsch_[i_step].weight );
				TR << "Changing scoreweight: " << mdsch_[i_step].scorename << " " << mdsch_[i_step].weight << std::endl;
			}
		}
	}

	/// Selection for returning structure
	if ( selectmode_.compare("final") == 0 ) {
		// just return final pose
		TR << "Returning final structure for MD..." << std::endl;
	} else if ( selectmode_.compare("minobj") == 0 ) {
		pose = pose_minobj_;
		TR << "Returning minimum objective function structure at ";
		TR << time_minobj_ << " in MD trajectory..." << std::endl;
	}

	do_minimize( pose, options_final, true );
	if ( report_scorecomp_ ) scorefxn_->show( TR, pose );

	TR << "MD Done. " << std::endl;
}


void CartesianMD::do_minimize( core::pose::Pose &pose,
	core::optimization::MinimizerOptions const &options,
	bool const &show_energy )
{
	CartesianMinimizer minimizer;

	if ( show_energy ) {
		Real score_before = scorefxn_->score( pose );

		minimizer.run( pose, *movemap(), *scorefxn_, options );

		Real score_after = scorefxn_->score( pose );
		scorefxn_->show(TR, pose);
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
	if ( initialize ) {
		do_initialize( pose );
	}

	// Reporting about adaptive rsr
	if ( constrained_ ) {
		TR << "Run uniform restrained simulation ";
		if ( Gamma_ == 0.0 ) {
			TR << " with static restraints on starting pose." << std::endl;
		} else {
			TR << " with dynamic restraints with Kappa/Gamma = ";
			TR << Kappa_ << " " << Gamma_ << std::endl;
		}
	}

	// Set dof variables
	CartesianMinimizerMap min_map;
	min_map.setup( pose, *movemap() );
	min_map.copy_dofs_from_pose( pose, xyz_ );

	// Setup RATTLE using min_map
	md::Rattle rattle( pose, min_map );
	if ( use_rattle_ ) n_dof_temp_ = n_dof_ - 6 - rattle.ncst();

	// This should come after Rattle setup to get n_dof_temp_
	if ( initialize ) initialize_velocity( temp0 );

	// Set thermostat
	Thermostat thermostat( temp0, n_dof_temp_ );

	// Start MD integrator
	scorefxn_->setup_for_minimizing( pose, min_map );

	for ( Size istep = 1; istep <= nstep; istep++ ) {
		cummulative_time_ += dt();

		// Report
		if ( istep%md_report_stepsize_ == 0 ) {
			report_MD( pose, min_map, true ); // report including trajectory
		} else if ( istep%md_energy_report_stepsize_ == 0 ) {
			report_MD( pose, min_map, false ); // only report energy
		}

		bool update_score( false );
		if ( istep%context_update_step_ == 0 ) update_score = true;

		// For adaptive restraint
		if ( istep%md_rsr_update_stepsize_ == 0 ) {
			update_restraint( pose, min_map );
		} else if ( trj_scratch_.size() < 100 ) {
			// make sure scratch space doesn't use too much memory
			trj_scratch_.push_back( xyz_ );
		}

		// Integrate
		VelocityVerlet_Integrator( pose, min_map, rattle, update_score );

		// Calculate/re-eval temperature
		temperature_ = thermostat.get_temperature( vel_, mass_ );

		if ( istep%thermostat.nstep_per_update() == 0 ) {
			thermostat.rescale( vel_, dt(), mass_ );
			temperature_ = thermostat.get_temperature( vel_, mass_ );
		}
		kinetic_energy_ = 0.5*temperature_*n_dof()*GasConst;

	}

	min_map.copy_dofs_to_pose( pose, xyz_ );

}

void CartesianMD::VelocityVerlet_Integrator( pose::Pose &pose,
	CartesianMinimizerMap &min_map,
	md::Rattle &rattle,
	bool const update_score )
{
	Real dt2_2 = dt()*dt()*0.5;

	// Use previous acceleration here
	// and integrate first half of the velociy
	for ( Size i_dof = 1; i_dof <= n_dof(); ++i_dof ) {
		xyz_[i_dof] += vel_[i_dof]*dt() + acc_[i_dof]*dt2_2;
		vel_[i_dof] += 0.5*acc_[i_dof]*dt();
		if ( vel_[i_dof] > MaxVel ) vel_[i_dof] = MaxVel;
		if ( vel_[i_dof] < -MaxVel ) vel_[i_dof] = -MaxVel;
	}

	if ( use_rattle_ ) {
		rattle.run_rattle1( dt(), xyz_, vel_, mass_ );
	}

	// Reflect change in coordinates into pose
	min_map.copy_dofs_to_pose( pose, xyz_ );


	// Do we need for non-FACTS?
	if ( update_score ) scorefxn_->score( pose );

	Multivec force;
	CartesianMultifunc f_ros( pose, min_map, *scorefxn_, false, false );
	f_ros.dfunc( xyz_, force );

	// Here, convert force into acceleration
	// and integrate remaining half of velocity
	for ( Size i_dof = 1; i_dof <= n_dof(); ++i_dof ) {
		Size i_atm = (i_dof+2)/3;
		// pass Virtual atoms
		if ( mass_[i_atm] < 1e-3 ) continue;

		acc_[i_dof] = -MDForceFactor*force[i_dof]/mass_[i_atm];
		// safe boundary
		if ( acc_[i_dof] > MaxAccel ) acc_[i_dof] = MaxAccel;
		if ( acc_[i_dof] < -MaxAccel ) acc_[i_dof] = -MaxAccel;

		vel_[i_dof] += 0.5*acc_[i_dof]*dt();
		if ( vel_[i_dof] > MaxVel ) vel_[i_dof] = MaxVel;
		if ( vel_[i_dof] < -MaxVel ) vel_[i_dof] = -MaxVel;
	}

	if ( use_rattle_ ) {
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

	// Supposed to make Maxwell-Boltzmann distribution... is this really working?
	// To make sure we should use error-function, but too lazy to do that...
	for ( core::Size i_dof=1; i_dof<=n_dof(); i_dof++ ) {
		core::Size i_atm = (i_dof+2)/3;
		// pass Virtual atoms
		if ( mass_[i_atm] < 1e-3 ) continue;

		Real scalar = sqrt(2.0*temperature*Boltzmann/mass_[i_atm])*numeric::random::rg().gaussian();
		if ( numeric::random::rg().uniform() > 0.5 ) scalar *= -1.0;
		vel_[i_dof] = scalar;
		acc_[i_dof] = 0.0;
	}

	// Uniformly scale down to make sure init temperature assigned correctly
	Thermostat thermostat( temperature, n_dof_temp_ );
	Real init_temp = thermostat.get_temperature( vel_, mass_ );
	Real const scale( temperature/init_temp );
	for ( core::Size i_dof=1; i_dof<=n_dof(); i_dof++ ) vel_[i_dof] *= scale;

	TR << "Initial temperature assigned as " << init_temp;
	TR << ", scaling down by factor " << scale << std::endl;
}

void CartesianMD::report_MD( core::pose::Pose &pose,
	CartesianMinimizerMap const &min_map,
	bool const report_trj )
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
	core::Real rmsd_native( 0.0 ), gdttm_native( 0.0 ), gdtha_native( 0.0 );
	if ( native_given_ ) {
		rmsd_native = core::scoring::CA_rmsd( pose, native_, native_resmap_ );
		core::scoring::CA_gdttm( pose, native_, gdttm_native, gdtha_native, native_resmap_ );

		TR << ", RMSD/GDTtoNative ";
		TR << std::setw(8) << std::setprecision(4) << rmsd_native;
		TR << std::setw(8) << std::setprecision(4) << gdttm_native;
		TR << std::setw(8) << std::setprecision(4) << gdtha_native;
	}

	core::Real Eobj( scorefxn_obj_->score( pose ) );
	TR << " " << Eobj << " " << Emin_obj_ ;
	TR << std::endl;

	if ( cummulative_time() > 0.1 && // Truncate initial 1ps to remove minimization memory
			selectmode_.compare("minobj") == 0 && Eobj < Emin_obj_ ) {
		pose_minobj_ = pose;
		Emin_obj_ = scorefxn_obj_->score( pose );
		time_minobj_ = cummulative_time();
		TR << "Updating minimum objective score value / pose at time " << cummulative_time() << std::endl;
	}

	// Store trj
	if ( report_trj && store_trj() ) {
		//CartesianMinimizerMap min_map;
		//min_map.setup( pose, *movemap() );

		Multivec xyz;
		xyz.resize( min_map.ndofs() );
		min_map.copy_dofs_from_pose( pose, xyz );
		trj_.push_back( xyz );
		if ( report_as_silent_ ) {
			if ( native_given_ ) {
				report_silent( pose, rmsd, gdttm_native, gdtha_native );
			} else {
				report_silent( pose );
			}
		}
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

	for ( Size i_trj = 1; i_trj <= trj_.size(); ++i_trj ) {
		min_map.copy_dofs_to_pose( pose_tmp, trj_[i_trj] );
		poses.push_back( pose_tmp );
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
	scorefxn_ = data.get_ptr<ScoreFunction>( "scorefxns", scorefxn_name ); //->clone();

	std::string const scoreobj_name( tag->getOption< std::string >( "scorefxn_obj","" ) );

	if ( scoreobj_name.compare("") == 0 ) {
		scorefxn_obj_ = scorefxn_; //->clone();
	} else {
		scorefxn_obj_ = data.get_ptr<ScoreFunction>( "scorefxns", scoreobj_name ); //->clone();
	}

	use_rattle_ = tag->getOption< bool >( "rattle", true );
	if ( use_rattle_ ) dt_ = 0.002;

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
	if ( schfile.compare("") != 0 ) {
		parse_schfile( schfile );
		scheduled_ = true;
	}

	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
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
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	bool const chi( tag->getOption< bool >( "chi", true ) ), bb( tag->getOption< bool >( "bb", true ) );
	movemap->set_chi( chi );
	movemap->set_bb( bb );
	set_movemap( pose, movemap );

	protocols::rosetta_scripts::parse_movemap( tag, pose, movemap_, data, false );
}

}
}
