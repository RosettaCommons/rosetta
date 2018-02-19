// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/select/movemap/MoveMapFactory.hh>
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
#include <core/scoring/rms_util.hh>
#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace md {

static basic::Tracer TR( "protocols.md.Cartesian" );

using namespace devel::md;
using namespace core::optimization;
using namespace core;
using namespace ObjexxFCL::format;

// creator
// XRW TEMP std::string
// XRW TEMP CartesianMDCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return CartesianMD::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP CartesianMDCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new CartesianMD );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CartesianMD::mover_name()
// XRW TEMP {
// XRW TEMP  return "CartesianMD";
// XRW TEMP }


// mover
CartesianMD::CartesianMD():
	use_rattle_( true )
{
	init();
}

CartesianMD::CartesianMD( core::pose::Pose const & pose,
	core::scoring::ScoreFunctionCOP sfxn,
	core::kinematics::MoveMapCOP movemap ) :
	use_rattle_( true )
{
	if ( movemap == nullptr ) {
		core::kinematics::MoveMapOP mmloc( new core::kinematics::MoveMap );
		mmloc->set_jump( true ); mmloc->set_bb( true ); mmloc->set_chi( true );
		mmloc->set( core::id::THETA, true ); mmloc->set( core::id::D, true);
		set_movemap( pose, mmloc );
	} else {
		set_movemap( pose, movemap );
	}

	set_scorefxn( sfxn );
	set_scorefxn_obj( sfxn );
	init();
	get_native_info( pose );
}

CartesianMD::CartesianMD( core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &sfxn ) :
	use_rattle_( true )
{
	core::kinematics::MoveMapOP mmloc( new core::kinematics::MoveMap );
	mmloc->set_jump( true ); mmloc->set_bb( true ); mmloc->set_chi( true );
	mmloc->set( core::id::THETA, true ); mmloc->set( core::id::D, true );
	set_movemap( pose, mmloc );

	set_scorefxn( sfxn );
	set_scorefxn_obj( sfxn );
	init();
	get_native_info( pose );
}

CartesianMD::CartesianMD( core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &sfxn,
	core::kinematics::MoveMap const &movemap ) :
	use_rattle_( true )
{
	set_movemap( pose, movemap.clone() );

	set_scorefxn( sfxn );
	set_scorefxn_obj( sfxn );
	init(); // MDBase
	get_native_info( pose );
}

CartesianMD::~CartesianMD() = default;

protocols::moves::MoverOP
CartesianMD::clone() const {
	return protocols::moves::MoverOP( new CartesianMD(*this) );
}

// XRW TEMP std::string CartesianMD::get_name() const
// XRW TEMP {
// XRW TEMP  return "CartesianMD";
// XRW TEMP }

void CartesianMD::get_native_info( core::pose::Pose const &pose )
{

	native_given_ = false;
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		native_given_ = true;

		std::string nativepdb = basic::options::option[ basic::options::OptionKeys::in::file::native ]();
		core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );
		core::import_pose::pose_from_file( native_, *rsd_set, nativepdb, core::import_pose::PDB_file );

		// Set resmap
		std::map< Size, Size > resmap;
		for ( Size ires = 1; ires <= pose.size(); ++ires ) {
			if ( !pose.residue( ires ).is_protein() ) continue;
			Size ii_pdb( pose.pdb_info()->number( ires ) );

			for ( Size jres = 1; jres <= native_.size(); ++jres ) {
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
}

void CartesianMD::do_initialize( core::pose::Pose &pose )
{

	// Steal utilities in CartesianMinimizer.cc
	CartesianMinimizerMap min_map;

	min_map.setup( pose, *movemap() );

	set_n_dof( min_map.ndofs() );
	set_n_dof_temp( n_dof() - 6 );
	set_cummulative_time( 0.0 );
	set_pose0( pose );

	// Check initial time
#ifndef WIN32
	gettimeofday(&inittime_, nullptr );
#endif

	// Reallocate
	resize_natm_variables();

	// for adaptive rsr
	Multivec xyz_loc( n_dof() );
	min_map.copy_dofs_from_pose( pose, xyz_loc );
	set_ref_xyz( xyz_loc );
	set_prv_eqxyz( xyz_loc );

	// Mass setup
	core::chemical::ResidueTypeSetCOP rsdtype_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );
	core::chemical::ElementSetCOP element_set( rsdtype_set->element_set() );

	for ( Size iatm = 1; iatm <= (Size)(min_map.natoms()); ++iatm ) {
		core::id::AtomID AtomID = min_map.get_atom( iatm );
		Size resno = AtomID.rsd();
		Size atmno = AtomID.atomno();

		set_mass( iatm, pose.residue_type(resno).atom(atmno).element_type()->weight() );
	}

}

void CartesianMD::use_rattle( bool const value )
{
	use_rattle_ = value;
	if ( use_rattle_ ) {
		set_dt( 0.002 );
		TR << "Set Rattle on, changing dt as " << dt() << std::endl;
	} else {
		set_dt( 0.001 );
		TR << "Set Rattle off, changing dt as " << dt() << std::endl;
	}
}

void
CartesianMD::update_restraint( core::pose::Pose & pose,
	CartesianMinimizerMap const &min_map )
{

	if ( Gamma() == 0.0 || !uniform_coord_constrained() ) return;

	TR.Debug << "Update restraints with Kappa/Gamma = " << Kappa() << "/" << Gamma() << std::endl;
	Multivec curr_eqxyz = get_current_eqxyz();
	Multivec prv_eqxyz_loc( prv_eqxyz() );
	cst_on_pose_dynamic( pose, ref_xyz(), curr_eqxyz, prv_eqxyz_loc, min_map );
	set_prv_eqxyz( prv_eqxyz_loc );

	// clear temporary trj
	renew_trj_scratch();
}

Multivec
CartesianMD::get_current_eqxyz() const
{
	Multivec curr_eqxyz;
	curr_eqxyz.resize( n_dof() );
	utility::vector1< Multivec > const trj_tmp = trj_scratch();
	core::Size const ntrj( trj_tmp.size() );

	for ( core::Size i_dof = 1; i_dof <= n_dof(); ++i_dof ) {
		curr_eqxyz[i_dof] = 0.0;
		for ( core::Size i_trj = 1; i_trj <= ntrj; ++i_trj ) {
			curr_eqxyz[i_dof] += trj_tmp[i_trj][i_dof];
		}

		if ( ntrj > 0 ) curr_eqxyz[i_dof] /= (core::Real)(ntrj);
	}
	return curr_eqxyz;
}

void
CartesianMD::cst_on_pose_simple( core::pose::Pose &pose ) const
{
	using namespace core::scoring::constraints;
	using namespace core;

	// Remove all the constraints first
	//pose.remove_constraints(); // fpd ... move this below!

	// First, add cst_file info into pose
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file ].user() ) {
		TR << "Set constraints from input file..." << std::endl;
		pose.remove_constraints();
		scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose );
	}

	// Next, set coordinate constraint if specified
	if ( uniform_coord_constrained() ) {
		TR << "Set constraints uniformly with stdev: " << cst_sdev() << std::endl;
		pose.remove_constraints();

		for ( Size i_res = 1; i_res <= pose.size(); ++i_res ) {
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
			scoring::func::FuncOP fx( new scoring::func::HarmonicFunc( 0.0, cst_sdev() ) );
			pose.add_constraint(  ConstraintCOP( ConstraintOP(
				new CoordinateConstraint( atomID, atomID, xyz, fx )
				)));
		}
	}
}

void
CartesianMD::cst_on_pose_dynamic( core::pose::Pose &pose,
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

	std::ofstream rsrfile( rsrfilename().c_str(), std::ios_base::app );
	core::Size modality = (core::Size)(cummulative_time()*100+1)%100;
	bool write_rsr = write_dynamic_rsr() && ( modality <= 2);

	TR.Debug << cummulative_time() << " " << ((core::Size)(cummulative_time()*100))%100 << " " << write_rsr << std::endl;
	TR.Debug << "ref/curr/prv_eqxyz0? " << ref_xyz.size() << " " << curr_eqxyz.size() << " " << prv_eqxyz0.size() << std::endl;

	if ( uniform_coord_constrained() ) {
		if ( write_rsr ) {
			rsrfile << "dynamic rsr: " << cummulative_time() << std::endl;
			rsrfile << "Res | X Y Z | prv_eq X Y Z | delta X Y Z" << std::endl;
		}
		core::Real rmsd( 0.0 );
		core::Real rmsd_rsr2ref( 0.0 );
		core::Real rmsd_rsr2crd( 0.0 );
		core::Size nrsr_dof( 0 );
		//for ( Size i_res = 1; i_res <= pose.size(); ++i_res ) {
		for ( Size i_atm = 1; i_atm <= n_dof()/3; ++i_atm ) {
			core::id::AtomID atomID = min_map.get_atom( i_atm );
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
				xyzmix[i-1] = (1.0 - Kappa())*prv_eqxyz[ i_dof ];
				xyzmix[i-1] += Kappa()*(Gamma()*curr_eqxyz[ i_dof ] +
					(1.0-Gamma())*ref_xyz[ i_dof ] );

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

			core::Real sdev( cst_sdev() );

			core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0,  sdev ) );
			pose.add_constraint( ConstraintCOP( ConstraintOP(
				new CoordinateConstraint( atomID, atomID, xyzmix, fx )
				)));

			if ( write_rsr ) {
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

		if ( write_rsr ) {
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

void CartesianMD::apply( core::pose::Pose & pose ) {
	using namespace core::optimization;

	if ( movemap_factory_ ) {
		// reset the movemap if we have a valid MoveMapFactory
		set_movemap( pose, movemap_factory_->create_movemap_from_pose( pose ) );
	}

	//fpd we have to do this here since this the first time "seeing" the symm pose
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap() );
	}

	// setup the map of the degrees of freedom
	const std::string minopt( "lbfgs_armijo_nonmonotone" );
	MinimizerOptions options_init( minopt, 0.0001, true, false, false );
	options_init.max_iter( ncyc_premin() );
	MinimizerOptions options_final( minopt, 0.000001, true, false, false );
	options_final.max_iter( ncyc_postmin() );

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

	if ( report_scorecomp() ) scorefxn()->show( TR, pose );

	// Set initial minobj_ after minimization
	// These are used to pick the best pose along trajectory (if selectmode == minobj )
	set_pose_minobj( pose );
	set_Emin_obj( 1.0e6 );
	set_time_minobj( 0.0 );

	// Main
	if ( !scheduled() ) { // typical run
		do_MD( pose, nstep(), temp0(), true );

	} else { // scheduled run

		for ( Size i_step = 1; i_step <= mdsch().size(); ++i_step ) {
			std::string const runtype( mdsch(i_step).type );
			if ( runtype.compare("sch") == 0 ) {
				TR << "Changing schedule, Nstep/Temp:";
				TR << mdsch(i_step).nstep << " " << mdsch(i_step).temp0 << std::endl;
				if ( i_step == 1 ) {
					do_MD( pose, mdsch(i_step).nstep, mdsch(i_step).temp0, true );
				} else {
					do_MD( pose, mdsch(i_step).nstep, mdsch(i_step).temp0, false );
				}
			} else if ( runtype.compare("repack") == 0 ) {
				//
			}

		}
	}

	/// Selection for returning structure
	if ( selectmode().compare("final") == 0 ) {
		// just return final pose
		TR << "Returning final structure for MD..." << std::endl;
	} else if ( selectmode().compare("minobj") == 0 ) {
		pose = pose_minobj();
		TR << "Returning minimum objective function structure at ";
		TR << time_minobj() << " in MD trajectory..." << std::endl;
	}

	do_minimize( pose, options_final, true );
	if ( report_scorecomp() ) scorefxn()->show( TR, pose );

	TR << "MD Done. " << std::endl;
}


void CartesianMD::do_minimize( core::pose::Pose &pose,
	core::optimization::MinimizerOptions const &options,
	bool const &show_energy )
{
	CartesianMinimizer minimizer;

	if ( show_energy ) {
		core::Real score_before = scorefxn()->score( pose );

		minimizer.run( pose, *movemap(), *scorefxn(), options );

		core::Real score_after = scorefxn()->score( pose );
		TR << "Energy before/after Min: " << score_before << " " << score_after << std::endl;

	} else {
		minimizer.run( pose, *movemap(), *scorefxn(), options );
	}
}

void CartesianMD::do_MD( core::pose::Pose & pose,
	core::Size const &nstep,
	core::Real const &temp0,
	bool const &initialize )
{
	if ( initialize ) {
		TR << "Running MD simulations for " << nstep*dt() << " ps at " << temp0 << " K... " << std::endl;
		do_initialize( pose );
	}

	// Reporting about adaptive rsr
	if ( uniform_coord_constrained() ) {
		TR << "Run uniform restrained simulation ";
		if ( Gamma() == 0.0 ) {
			TR << " with static restraints on starting pose." << std::endl;
		} else {
			TR << " with dynamic restraints with Kappa/Gamma = ";
			TR << Kappa() << " " << Gamma() << std::endl;
		}
	}

	// Set dof variables
	CartesianMinimizerMap min_map;
	Multivec xyz_loc( n_dof() );
	min_map.setup( pose, *movemap() );
	min_map.copy_dofs_from_pose( pose, xyz_loc );
	set_xyz( xyz_loc );

	// Setup RATTLE using min_map
	md::Rattle rattle( pose, min_map );
	if ( use_rattle_ ) set_n_dof_temp( n_dof() - 6 - rattle.ncst() );

	// This should come after Rattle setup to get n_dof_temp_
	if ( initialize ) initialize_velocity( temp0 );

	// Set thermostat
	Thermostat thermostat( temp0, n_dof_temp() );

	// Start MD integrator
	scorefxn()->setup_for_minimizing( pose, min_map );

	for ( Size istep = 1; istep <= nstep; istep++ ) {
		set_cummulative_time( cummulative_time() + dt() );

		// Report
		if ( istep%md_report_stepsize() == 0 ) {
			report_MD( pose, min_map, true ); // report including trajectory
		} else if ( istep%md_energy_report_stepsize() == 0 ) {
			report_MD( pose, min_map, false ); // only report energy
		}

		bool update_score( false );
		if ( istep%context_update_step() == 0 ) update_score = true;

		// For adaptive restraint
		if ( istep%md_rsr_update_stepsize() == 0 ) {
			update_restraint( pose, min_map );
		} else if ( trj_scratch().size() < 100 ) {
			// make sure scratch space doesn't use too much memory
			add_trj_scratch( xyz() );
		}

		// Integrate
		VelocityVerlet_Integrator( pose, min_map, rattle, update_score );

		/*
		// let's get avrg velocities
		core::Real v2sum( 0.0 ), temperature( 0.0 );
		for (Size i_dof = 1; i_dof<=vel_.size(); ++i_dof) {
		int i_atm = (i_dof+2)/3;
		// pass Virtual atoms
		if ( mass_[i_atm] < 1e-3 ) continue;
		TR << "vel: " << i_atm << " " << vel_[i_dof] << std::endl;
		temperature += mass_[i_atm]*vel_[i_dof]*vel_[i_dof];
		v2sum += vel_[i_dof]*vel_[i_dof];
		}
		temperature /= n_dof_temp_*Boltzmann;
		*/

		// Calculate/re-eval temperature
		set_temperature( thermostat.get_temperature( vel(), mass() ) );

		if ( istep%thermostat.nstep_per_update() == 0 ) {
			Multivec vel_loc( vel() );
			thermostat.rescale( vel_loc, dt(), mass() );
			set_vel( vel_loc );
			set_temperature( thermostat.get_temperature( vel(), mass() ) );
		}
		//TR << "v2avrg/Temp/Temp2: " << std::sqrt(v2sum/vel_.size()) << " " << temperature << " " << temperature_ << std::endl;
		set_kinetic_energy( 0.5*temperature()*n_dof()*GasConst );

	}

	min_map.copy_dofs_to_pose( pose, xyz() );

}

void CartesianMD::VelocityVerlet_Integrator( core::pose::Pose &pose,
	CartesianMinimizerMap &min_map,
	md::Rattle &rattle,
	bool const update_score )
{
	core::Real dt2_2 = dt()*dt()*0.5;

	// Use previous acceleration here
	// and integrate first half of the velociy
	Multivec xyz_loc( xyz() ), vel_loc( vel() ), acc_loc( acc() );

	for ( Size i_dof = 1; i_dof <= n_dof(); ++i_dof ) {
		xyz_loc[i_dof] += vel_loc[i_dof]*dt() + acc_loc[i_dof]*dt2_2;
		vel_loc[i_dof] += 0.5*acc_loc[i_dof]*dt();
	}

	if ( use_rattle_ ) {
		rattle.run_rattle1( dt(), xyz_loc, vel_loc, mass() );
	}

	// Reflect change in coordinates into pose
	min_map.copy_dofs_to_pose( pose, xyz_loc );

	// Don't need this unless context needs to be updated
	if ( update_score ) scorefxn()->score( pose );

	Multivec force;
	CartesianMultifunc f_ros( pose, min_map, *scorefxn(), false, false );
	f_ros.dfunc( xyz_loc, force );

	// Here, convert force into acceleration
	// and integrate remaining half of velocity
	for ( Size i_dof = 1; i_dof <= n_dof(); ++i_dof ) {
		Size i_atm = (i_dof+2)/3;
		// pass Virtual atoms
		if ( mass(i_atm) < 1e-3 ) continue;

		//acc(i_dof) = -MDForceFactor*force[i_dof]/mass(i_atm);
		//vel(i_dof) += 0.5*acc(i_dof)*dt();
		acc_loc[i_dof] = -MDForceFactor*force[i_dof]/mass(i_atm);
		vel_loc[i_dof] += 0.5*acc_loc[i_dof]*dt();
	}

	if ( use_rattle_ ) {
		rattle.run_rattle2( dt(), xyz_loc, vel_loc, mass() );
	}
	set_xyz( xyz_loc );
	set_vel( vel_loc );
	set_acc( acc_loc );

	//Stop rotation and translation
	//if ((step%nrottrans)==0) {
	//  stop_rot_trans(xyz, vel, size);
	//  thermostat_.get_temperature();
	//}
}

void CartesianMD::initialize_velocity( core::Real const &temperature )
{

	TR.Debug << "Setting initial velocity with temp = " << temperature << std::endl;
	// Supposed to make Maxwell-Boltzmann distribution... is this really working?
	// To make sure we should use error-function, but too lazy to do that...
	Multivec vel_loc( n_dof() ), acc_loc( n_dof(), 0.0 );
	for ( core::Size i_dof=1; i_dof<=n_dof(); i_dof++ ) {
		core::Size i_atm = (i_dof+2)/3;
		// pass Virtual atoms
		if ( mass(i_atm) < 1e-3 ) continue;

		Real scalar = sqrt(2.0*temperature*Boltzmann/mass(i_atm))*numeric::random::rg().gaussian();
		if ( numeric::random::rg().uniform() > 0.5 ) scalar *= -1.0;

		vel_loc[i_dof] = scalar;
		//acc(i_dof) = 0.0;
		//printf("%4d %4d %8.3f %8.3f %8.3f\n",int(i_dof),int(i_atm),scalar,mass(i_atm),vel_loc[i_dof]);
	}

	// Uniformly scale down to make sure init temperature assigned correctly
	Thermostat thermostat( temperature, n_dof_temp() );
	Real init_temp = thermostat.get_temperature( vel_loc, mass() );
	Real const scale( temperature/init_temp );
	for ( core::Size i_dof=1; i_dof<=n_dof(); i_dof++ ) vel_loc[i_dof] *= scale;

	set_vel( vel_loc );

	TR << "Initial temperature assigned as " << init_temp;
	TR << ", scaling down by factor " << scale << std::endl;
}

void CartesianMD::report_MD( core::pose::Pose &pose,
	CartesianMinimizerMap const &min_map,
	bool const report_trj )
{
	core::Real const rmsd( core::scoring::CA_rmsd( pose0(), pose ));

	core::scoring::constraints::ConstraintSetCOP cstset( pose.constraint_set() );
	//TR << "Is there CST? " << std::endl;
	//cstset->show_numbers( TR );

	timeval currtime;
#ifndef WIN32
	gettimeofday(&currtime, nullptr );
#endif
	Real elapsedTime = (currtime.tv_sec - inittime_.tv_sec) * 1000.0;
	elapsedTime += (currtime.tv_usec - inittime_.tv_usec) / 1000.0;
	elapsedTime /= 60000.0; // in minute

	core::Real Epot( scorefxn()->score( pose ) );
	//scorefxn_->show( pose );

	TR << "Time/E/Temp/RMSD/Elapsed(Min): " << std::setw(8) << cummulative_time();
	TR << " " << std::setw(12) << std::setprecision(6) << Epot;
	TR << " " << std::setw(6) << std::setprecision(4) << temperature();
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

	core::Real Eobj( scorefxn_obj()->score( pose ) );
	TR << " " << Eobj << " " << Emin_obj() ;
	TR << std::endl;

	if ( cummulative_time() > 0.1 && // Truncate initial 1ps to remove minimization memory
			selectmode().compare("minobj") == 0 && Eobj < Emin_obj() ) {
		set_pose_minobj( pose );
		set_Emin_obj( scorefxn_obj()->score( pose ) );
		set_time_minobj(  cummulative_time() );
		TR << "Updating minimum objective score value / pose at time " << cummulative_time() << std::endl;
	}

	// Store trj
	if ( report_trj && store_trj() ) {
		//CartesianMinimizerMap min_map;
		//min_map.setup( pose, *movemap() );

		Multivec xyz_loc( min_map.ndofs() );
		min_map.copy_dofs_from_pose( pose, xyz_loc );
		add_trj( xyz_loc );
		if ( report_as_silent() ) {
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
	for ( Size i_atm = 1; i_atm <= n_dof()/3; ++i_atm ) {
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

	for ( Size i_trj = 1; i_trj <= trj().size(); ++i_trj ) {
		min_map.copy_dofs_to_pose( pose_tmp, trj(i_trj) );
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
	set_scorefxn( data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_name ) );

	std::string const scoreobj_name( tag->getOption< std::string >( "scorefxn_obj","" ) );

	if ( scoreobj_name.compare("") == 0 ) {
		set_scorefxn_obj( scorefxn() );
	} else {
		set_scorefxn_obj( data.get_ptr< ScoreFunction >( "scorefxns", scoreobj_name ) );
	}

	use_rattle_ = tag->getOption< bool >( "rattle", true );
	if ( use_rattle_ ) set_dt( 0.002 );

	set_nstep( tag->getOption< core::Size >( "nstep", 100 ) );
	set_temp0( tag->getOption<core::Real>("temp", 300.0) );
	set_scheduled( false );

	set_ncyc_premin( tag->getOption< core::Size >( "premin", 50 ) );
	set_ncyc_postmin( tag->getOption< core::Size >( "postmin", 200 ) );

	set_md_report_stepsize( tag->getOption< core::Size >( "report", 100 ) );
	set_report_scorecomp( tag->getOption< bool >( "report_scorecomp", false ) );
	set_selectmode( tag->getOption< std::string >( "selectmode", "final" ) );

	// Use parsed schedule file - this will overload nstep, temperature, etc. defined above
	std::string const schfile( tag->getOption< std::string >( "schfile","" ) );
	if ( schfile.compare("") != 0 ) {
		parse_schfile( schfile );
		set_scheduled( true );
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
	Pose const & )
{
	// set initial guess
	core::select::movemap::MoveMapFactoryOP mmf( new core::select::movemap::MoveMapFactory );
	bool const chi( tag->getOption< bool >( "chi", true ) ), bb( tag->getOption< bool >( "bb", true ) );
	mmf->all_chi( chi );
	mmf->all_bb( bb );

	movemap_factory_ = protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data, false, mmf );
}

std::string CartesianMD::get_name() const {
	return mover_name();
}

std::string CartesianMD::mover_name() {
	return "CartesianMD";
}

void CartesianMD::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "rattle", xsct_rosetta_bool,
		"Use Rattle algorithm to constraint hydrogen locations. "
		"This automatically sets integration step = 2fs. "
		"Otherwise uses integration step = 1fs" );
	attlist + XMLSchemaAttribute( "scorefxn", xs_string,
		"Specify a scorefunction to run MD simulation with" );
	attlist + XMLSchemaAttribute( "scorefxn_obj", xs_string,
		"Optional, identical to scorefxn unless specified. "
		"Specify a scorefunction to use as objective function "
		"for selecting a pose from trajectory. "
		"This will be used only when selectmode=\"minobj\"" );
	attlist + XMLSchemaAttribute( "nstep", xs_integer,
		"Number of steps to simulate. "
		"With Rattle on (default) each step is 2fs, "
		"and hence, nstep=10000 will be 20ps" );
	attlist + XMLSchemaAttribute( "temp", xsct_real,
		"Reference temperature for constant temperature simulation. "
		"Recommended values: "
		"150~200K for talaris2014_cart and ~250 for beta_nov15_cart" );
	attlist + XMLSchemaAttribute( "premin", xs_integer,
		"Steps of Cartesian minimization before MD simulation" );
	attlist + XMLSchemaAttribute( "postmin", xs_integer,
		"Steps of Cartesian minimization after MD simulation" );
	attlist + XMLSchemaAttribute( "report", xs_integer,
		"By how often the mover reports the simulation status to log" );
	attlist + XMLSchemaAttribute( "report_scorecomp", xsct_rosetta_bool,
		"Whether to report score components to log" );
	attlist + XMLSchemaAttribute( "selectmode", xs_string,
		"How to select single pose from the trajectory. "
		"\"final\" to take the final pose, \"minobj\" to take the "
		"lowest objective function (by scorefxn_obj) pose" );
	attlist + XMLSchemaAttribute( "schfile", xs_string,
		"Use user-defined schedule file. "
		"This overrides any other flags or options. "
		"Syntax: \"sch [temperature] [nsteps]\" to run simulation, "
		"or \"repack\" to repack side-chains" );
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"CartesianMD calls Molecular Dynamics simulation in Rosetta "
		"with user-defined energy function. Runs NVT simulation "
		"(constant volume and temperature) with Berendsen thermostat. "
		"Integrator uses Velocity Verlet algorithm", attlist );
}

std::string CartesianMDCreator::keyname() const {
	return CartesianMD::mover_name();
}

protocols::moves::MoverOP
CartesianMDCreator::create_mover() const {
	return protocols::moves::MoverOP( new CartesianMD );
}

void CartesianMDCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CartesianMD::provide_xml_schema( xsd );
}


}
}
