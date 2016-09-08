// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio

#include <protocols/cryst/refinable_lattice.hh>
#include <protocols/cryst/refinable_lattice_creator.hh>
#include <protocols/cryst/util.hh>
#include <protocols/cryst/wallpaper.hh>
#include <protocols/cryst/spacegroup.hh>

#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>

#include <core/scoring/electron_density/util.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymRotamerTrialsMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/relax/FastRelax.hh>



#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <core/scoring/constraints/util.hh>

#include <utility/excn/Exceptions.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


using basic::T;
using basic::Error;
using basic::Warning;

#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.295779513082323


namespace protocols {
namespace cryst {

static basic::Tracer TR("protocols.cryst.refinable_lattice");

using namespace protocols;
using namespace core;
using namespace kinematics;
using namespace scoring;
using namespace scoring::symmetry;
using namespace conformation;
using namespace conformation::symmetry;

//static int OUTCOUNTER=1;

////////////////////////////////////////////////////

// creators
std::string
UpdateCrystInfoCreator::keyname() const { return UpdateCrystInfoCreator::mover_name(); }

protocols::moves::MoverOP
UpdateCrystInfoCreator::create_mover() const { return protocols::moves::MoverOP(new UpdateCrystInfo); }

std::string
UpdateCrystInfoCreator::mover_name() { return "UpdateCrystInfo"; }

std::string
DockLatticeMoverCreator::keyname() const { return DockLatticeMoverCreator::mover_name(); }

protocols::moves::MoverOP
DockLatticeMoverCreator::create_mover() const { return protocols::moves::MoverOP(new DockLatticeMover); }

std::string
DockLatticeMoverCreator::mover_name() { return "DockLatticeMover"; }

std::string
MakeLatticeMoverCreator::keyname() const { return MakeLatticeMoverCreator::mover_name(); }

protocols::moves::MoverOP
MakeLatticeMoverCreator::create_mover() const { return protocols::moves::MoverOP(new MakeLatticeMover); }

std::string
MakeLatticeMoverCreator::mover_name() { return "MakeLatticeMover"; }

////////////////////////////////////////////////////

void
UpdateCrystInfo::apply( core::pose::Pose & pose ) {
	SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	Size Ajump_=0, Bjump_=0, Cjump_=0, SUBjump_=0;

	core::pose::Pose pose_asu;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);

	// find lattice jumps
	for ( Size i=1; i<=pose.fold_tree().num_jump(); ++i ) {
		std::string jumpname = SymmConf.Symmetry_Info()->get_jump_name( i );
		if ( jumpname == "A" ) Ajump_=i;
		else if ( jumpname == "B" ) Bjump_=i;
		else if ( jumpname == "C" ) Cjump_=i;
		else if ( jumpname == "SUB" ) SUBjump_=i;
	}
	runtime_assert( Ajump_!=0 &&  Bjump_!=0 ); // && Cjump_ != 0 );

	Vector Axform = pose.jump(Ajump_).get_translation();
	Vector Bxform = pose.jump(Bjump_).get_translation();
	Vector Cxform(-1000,0,0);

	if ( Cjump_ != 0 ) Cxform = pose.jump(Cjump_).get_translation();

	Bxform = Vector( Bxform[1], Bxform[2], Bxform[0] );
	Cxform = Vector( Cxform[2], Cxform[0], Cxform[1] );

	io::CrystInfo ci = pose.pdb_info()->crystinfo();

	bool need_angles=false;
	numeric::xyzVector<Size> grid;
	if ( Cjump_ != 0 ) {
		Spacegroup sg;
		sg.set_spacegroup(ci.spacegroup());
		if ( sg.setting() == TRICLINIC || sg.setting() == MONOCLINIC ) {
			need_angles = true;
		}
		grid = sg.get_nsubdivisions();
	} else {
		WallpaperGroup wg;
		wg.set_wallpaper_group(ci.spacegroup());
		if ( wg.setting() == wgMONOCLINIC ) {
			need_angles = true;
		}
		grid = wg.get_nsubdivisions();
	}

	Real A = Axform.length()*grid[0];
	Real B = Bxform.length()*grid[1];
	Real C = Cxform.length()*grid[2];
	ci.A(A); ci.B(B); ci.C(C);

	if ( need_angles ) {
		Axform = Vector(Axform[0], Axform[1], Axform[2]);
		Bxform = Vector(Bxform[1], Bxform[2], Bxform[0]);
		Cxform = Vector(-Cxform[2], Cxform[0], Cxform[1]);
		Real alpha = RAD2DEG*acos( Bxform.dot(Cxform) / (Bxform.length()*Cxform.length()) );
		Real beta = RAD2DEG*acos( Axform.dot(Cxform) / (Axform.length()*Cxform.length()) );
		Real gamma = RAD2DEG*acos( Axform.dot(Bxform) / (Axform.length()*Bxform.length()) );
		ci.alpha(alpha); ci.beta(beta); ci.gamma(gamma);
	}

	// origin jump root should be at (0,0,0)
	Vector k =  pose.residue( pose.fold_tree().jump_edge( SUBjump_ ).start() ).xyz("ORIG");
	pose_asu.apply_transform_Rx_plus_v( numeric::xyzMatrix<Real>::identity(),-k );

	// now delete the VRT
	pose = pose_asu;
	pose.conformation().delete_residue_slow( pose.size() );
	pose.pdb_info()->set_crystinfo(ci);
}

// parse_my_tag
void
UpdateCrystInfo::parse_my_tag(
	utility::tag::TagCOP const /*tag*/,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & ,
	moves::Movers_map const & ,
	core::pose::Pose const & /*pose*/ )
{

}

////////////////////////////////////////////////////

DockLatticeMover::DockLatticeMover(core::scoring::ScoreFunctionOP sf_in) {
	sf_ = sf_in->clone();
	sf_vdw_ = core::scoring::ScoreFunctionOP( new core::scoring::symmetry::SymmetricScoreFunction() );
	sf_vdw_->set_weight( core::scoring::vdw, 1.0 );

	SUBjump_ = 0;
	rot_mag_ = 0.5;
	trans_mag_ = 0.2;
	temp_ = 2.0;
	ncycles_ = 400;
	fullatom_ = false;
}

DockLatticeMover::DockLatticeMover() {
	sf_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth");
	sf_vdw_ = core::scoring::ScoreFunctionOP( new core::scoring::symmetry::SymmetricScoreFunction() );
	sf_vdw_->set_weight( core::scoring::vdw, 1.0 );

	SUBjump_ = 0;
	rot_mag_ = 0.5;
	trans_mag_ = 0.2;
	temp_ = 2.0;
	fullatom_ = false;
	design_ = false;
}

void
DockLatticeMover::init(core::pose::Pose & pose) {
	SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	symdofs_ = SymmConf.Symmetry_Info()->get_dofs();

	core::pose::Pose pose_asu;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);

	if ( pose.is_centroid() ) {
		monomer_bump_ = 2*((*sf_vdw_)( pose_asu )); // * SymmConf.Symmetry_Info()->subunits();
	}

	// find lattice jumps
	for ( Size i=1; i<=pose.fold_tree().num_jump(); ++i ) {
		std::string jumpname = SymmConf.Symmetry_Info()->get_jump_name( i );
		if ( jumpname == "SUB" ) SUBjump_=i;
	}
	runtime_assert( SUBjump_ != 0 );
}

void
DockLatticeMover::perturb_trial( core::pose::Pose & pose ) {
	std::map< Size, SymDof >::iterator it, it_begin = symdofs_.begin(), it_end = symdofs_.end();
	for ( it = it_begin; it != it_end; ++it ) {
		core::kinematics::Jump flexible_jump = pose.jump( it->first );
		SymDof const &dof = it->second;

		for ( Size i = 1; i<= 3; ++i ) {
			if ( dof.allow_dof(i) ) {
				flexible_jump.gaussian_move_single_rb( 1, trans_mag_, i );
			}
		}

		for ( Size i = 4; i<= 6; ++i ) {
			if ( dof.allow_dof(i) ) {
				flexible_jump.gaussian_move_single_rb( 1, rot_mag_, i );
			}
		}

		pose.set_jump( it->first, flexible_jump );
	}
}

void
DockLatticeMover::min_lattice( core::pose::Pose & pose ) {
	protocols::simple_moves::SwitchResidueTypeSetMover to_cen("centroid");
	protocols::simple_moves::SwitchResidueTypeSetMover to_fa("fa_standard");
	protocols::simple_moves::ReturnSidechainMoverOP restore_sc;

	if ( pose.is_fullatom() ) {
		restore_sc = protocols::simple_moves::ReturnSidechainMoverOP( new protocols::simple_moves::ReturnSidechainMover( pose ) );
		to_cen.apply(pose);
	}

	core::kinematics::MoveMapOP mm(new core::kinematics::MoveMap);
	mm->set_jump(true); mm->set_chi(false); mm->set_bb(false);
	core::pose::symmetry::make_symmetric_movemap( pose, *mm );

	core::scoring::ScoreFunctionOP sf = core::scoring::ScoreFunctionFactory::create_score_function("score1");

	protocols::simple_moves::symmetry::SymMinMoverOP min(
		new protocols::simple_moves::symmetry::SymMinMover(mm, sf, "lbfgs_armijo", 0.01, true) );
	min->apply(pose);

	if ( restore_sc ) {
		to_fa.apply(pose);
		restore_sc->apply(pose);
	}
}


void
DockLatticeMover::slide_lattice( core::pose::Pose & pose ) {
	Real min_step_size=0.25, max_step_size=1;  // keep max step small to stop subunits from jumping through one another
	Real max_sub_step_size=0.5;  // keep max step small to stop subunits from jumping through one another
	Real bump_threshold=1.0;
	Size tries=0, max_tries=100;
	Real score = (*sf_vdw_)( pose );

	std::map< Size, SymDof >::iterator it, it_begin = symdofs_.begin(), it_end = symdofs_.end();
	std::map< Size, Vector > step;
	numeric::xyzVector< Real > substep;

	std::map< Size, Vector > orig_trans;
	for ( it = it_begin; it != it_end; ++it ) {
		if ( it->first == SUBjump_ ) {
			step[it->first]=Vector(0.5,0.5,0.5);
		} else {
			step[it->first]=Vector(-1.0,0.0,0.0);
		}
		orig_trans[it->first] = pose.jump( it->first ).get_translation();
	}

	//static int xx=1;
	//pose.dump_pdb("slidepre"+utility::to_string(xx)+".pdb");

	// 1 slide out
	bool done = false;
	while ( tries < max_tries && !done ) {
		done = true;
		for ( it = it_begin; it != it_end; ++it ) {
			for ( int i = X_DOF; i <= Z_DOF; ++i ) {
				if ( !it->second.allow_dof(i) ) continue;

				core::Size ii = 0;
				if ( i==X_DOF ) ii=0;
				if ( i==Y_DOF ) ii=1;
				if ( i==Z_DOF ) ii=2;

				core::kinematics::Jump jump_i = pose.jump( it->first );
				numeric::xyzVector<core::Real> T_curr = jump_i.get_translation();
				numeric::xyzVector<core::Real> move_dir(0,0,0);
				move_dir[ii] = step[it->first][ii];
				if ( it->first == SUBjump_ ) {
					numeric::xyzVector<core::Real> T_unit = T_curr; T_unit.normalize();
					move_dir[ii] *= T_unit[ii];
				}
				core::Real score_old = (*sf_vdw_)( pose );
				jump_i.set_translation( T_curr + move_dir );
				pose.set_jump( it->first, jump_i );
				score = (*sf_vdw_)( pose );

				// if this move didnt change the vdw energy, reset
				if ( std::fabs(score-score_old) < 1e-2 ) {
					jump_i.set_translation( T_curr );
					pose.set_jump( it->first, jump_i );
				}

				// if we are still above bump thresh, we are not done
				if ( (score-monomer_bump_) > bump_threshold ) { // BUMP
					done = false;
				}

				//std::cerr << "slide-out " << move_dir << " " << score-monomer_bump_ << std::endl;
				//pose.dump_pdb("slideout"+utility::to_string(xx++)+".pdb");

				tries++;
			}
		}
	}

	//pose.dump_pdb("slideout"+utility::to_string(xx)+".pdb");

	if ( !done ) {
		std::cerr << "slide-out fail after " << max_tries << " attempts" << std::endl;
		//pose.dump_pdb("slidefail.pdb");
		for ( it = it_begin; it != it_end; ++it ) {
			core::kinematics::Jump jump_i = pose.jump( it->first );
			jump_i.set_translation( orig_trans[it->first] );
			pose.set_jump( it->first, jump_i );
		}
		return;
	}

	numeric::xyzVector<core::Real> T_unit = pose.jump(SUBjump_).get_translation();
	T_unit.normalize();
	step[SUBjump_] = T_unit;

	// 2 slide in
	tries=0;
	done = false;
	while ( tries < max_tries && !done ) {
		done = true;
		it_begin = symdofs_.begin();
		it_end = symdofs_.end();
		for ( it = it_begin; it != it_end; ++it ) {
			//Size ndims=(it->first == SUBjump_)?3:1;

			for ( int i = X_DOF; i <= Z_DOF; ++i ) {
				if ( !it->second.allow_dof(i) ) continue;

				core::Size ii = 0;
				if ( i==X_DOF ) ii=0;
				if ( i==Y_DOF ) ii=1;
				if ( i==Z_DOF ) ii=2;

				core::kinematics::Jump jump_i = pose.jump( it->first );
				numeric::xyzVector<core::Real> T_curr = jump_i.get_translation();
				numeric::xyzVector<core::Real> move_dir(0,0,0);
				move_dir[ii] = step[it->first][ii];
				//if (it->first == SUBjump_) {
				// numeric::xyzVector<core::Real> T_unit = T_curr; T_unit.normalize();
				// move_dir[ii] *= T_unit[ii];
				//}
				//core::Real score_old = (*sf_vdw_)( pose );
				jump_i.set_translation( T_curr - move_dir );
				pose.set_jump( it->first, jump_i );
				score = (*sf_vdw_)( pose );

				if ( (score-monomer_bump_) > bump_threshold ) { // BUMP
					jump_i.set_translation( T_curr );
					pose.set_jump( it->first, jump_i );
					if ( std::fabs(step[it->first][ii]) > std::fabs(min_step_size+1e-4) ) {
						done = false;
						step[it->first][ii] /= 2.0; // take smaller steps
					}
				} else { // NO BUMP
					done = false;
					if ( it->first == SUBjump_ ) {
						if ( std::fabs(step[it->first][ii]) < std::fabs(max_sub_step_size-1e-4) ) {
							step[it->first][ii] *= 2.0;
						}
					} else {
						if ( std::fabs(step[it->first][ii]) < std::fabs(max_step_size-1e-4) ) {
							step[it->first][ii] *= 2.0;
						}
					}
				}

				//std::cerr << "slide-in " << move_dir << " " << score-monomer_bump_ << std::endl;
				//pose.dump_pdb("slidein"+utility::to_string(xx++)+".pdb");

				tries++;

			}
		}
	}


	if ( !done ) {
		std::cerr << "slide-in fail after " << max_tries << " attempts" << std::endl;
		for ( it = it_begin; it != it_end; ++it ) {
			if ( it->first == SUBjump_ ) continue;
			core::kinematics::Jump jump_i = pose.jump( it->first );
			jump_i.set_translation( orig_trans[it->first] );
			pose.set_jump( it->first, jump_i );
		}
	}

}


void
DockLatticeMover::modify_lattice( core::pose::Pose & pose, core::Real mag ) {
	std::map< Size, SymDof >::iterator it, it_begin = symdofs_.begin(), it_end = symdofs_.end();
	for ( it = it_begin; it != it_end; ++it ) {
		if ( it->first == SUBjump_ ) continue;
		core::kinematics::Jump jump_i = pose.jump( it->first );
		jump_i.set_translation( jump_i.get_translation() - mag*numeric::xyzVector<core::Real>(1,0,0) );  // neg. direction expands lattice
	}
}

void
DockLatticeMover::apply( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::simple_moves::SwitchResidueTypeSetMover to_cen("centroid");
	protocols::simple_moves::SwitchResidueTypeSetMover to_fa("fa_standard");

	if ( !fullatom_ ) {
		if ( pose.is_fullatom() ) to_cen.apply(pose);
		init(pose); // find which are lattice jumps && which are dof jumps
		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo(pose, *sf_, temp_ ) );

		for ( int i=1; i<(int)ncycles_; ++i ) {
			if ( i%24 == 0 ) {
				mc->show_scores();
				mc->show_counters();
			}
			perturb_trial(pose);
			if ( i%24 == 0 ) {
				slide_lattice(pose);
			}
			//pose.dump_pdb("dock_"+utility::to_string( OUTCOUNTER++ )+".pdb");
			mc->boltzmann( pose );
		}
		mc->recover_low( pose );
		//to_fa.apply(pose);
	} else {
		if ( !pose.is_fullatom() ) to_fa.apply(pose);

		// set up packer
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		using namespace protocols::toolbox::task_operations;
		TaskFactoryOP tf (new TaskFactory);
		tf->push_back( TaskOperationCOP(new InitializeFromCommandline) );
		tf->push_back( TaskOperationCOP(new IncludeCurrent) );
		if ( !design_ ) {
			tf->push_back( TaskOperationCOP(new RestrictToRepacking) );
			tf->push_back( TaskOperationCOP(new NoRepackDisulfides) );
		}
		tf->push_back( TaskOperationCOP(new RestrictToInterface( SUBjump_ ) ) );
		protocols::simple_moves::PackRotamersMoverOP pack_interface_repack( new protocols::simple_moves::symmetry::SymPackRotamersMover( sf_ ) );
		pack_interface_repack->task_factory(tf);
		protocols::simple_moves::RotamerTrialsMoverOP pack_interface_rtrials( new protocols::simple_moves::symmetry::SymRotamerTrialsMover( sf_, tf ) );

		// set up minimizer
		core::kinematics::MoveMapOP mm(new core::kinematics::MoveMap);
		mm->set_jump(true); mm->set_chi(true); mm->set_bb(false);
		core::pose::symmetry::make_symmetric_movemap( pose, *mm );
		protocols::simple_moves::symmetry::SymMinMoverOP min(
			new protocols::simple_moves::symmetry::SymMinMover(mm, sf_, "lbfgs_armijo", 0.01, true) );

		init(pose); // find which are lattice jumps && which are dof jumps
		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo(pose, *sf_, 2.0 ) );
		for ( int i=1; i<(int)ncycles_; ++i ) {
			if ( i%24 == 0 ) {
				mc->show_scores();
				mc->show_counters();
			}

			perturb_trial(pose);

			if ( i%24 == 1 ) {
				pack_interface_repack->apply( pose );
			} else {
				pack_interface_rtrials->apply( pose );
			}

			min->apply(pose);
			mc->boltzmann( pose );
		}
		mc->recover_low( pose );
	}
}


// parse_my_tag
void
DockLatticeMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const & ,
	moves::Movers_map const & ,
	core::pose::Pose const & /*pose*/ )
{
	fullatom_ = tag->getOption<bool>( "fullatom", true );
	if ( fullatom_ ) {
		sf_ = core::scoring::ScoreFunctionFactory::create_score_function("talaris2013");
	} else {
		sf_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth");
	}

	if ( tag->hasOption( "scorefxn" ) ) {
		sf_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	}

	if ( tag->hasOption( "ncycles" ) ) {
		ncycles_ = tag->getOption<core::Size>( "ncycles" );
	}
	if ( tag->hasOption( "trans_step" ) ) {
		trans_mag_ = tag->getOption<core::Real>( "trans_step" );
	}
	if ( tag->hasOption( "rot_step" ) ) {
		rot_mag_ = tag->getOption<core::Real>( "rot_step" );
	}


}

////////////////////////////////////////////////////

void
MakeLatticeMover::apply( core::pose::Pose & pose ) {
	// initialize sg_ from pose CRYST1 line
	io::CrystInfo ci = pose.pdb_info()->crystinfo();
	runtime_assert(ci.A()*ci.B()*ci.C() != 0);  // TODO: allow these to be randomized

	sg_.set_spacegroup(ci.spacegroup());
	sg_.set_parameters(ci.A(),ci.B(),ci.C(), ci.alpha(), ci.beta(), ci.gamma());

	// get principal axes
	Vector Ax, Bx, Cx;
	Ax = sg_.f2c()*Vector(1,0,0); Ax.normalize();
	Bx = sg_.f2c()*Vector(0,1,0); Bx.normalize();
	Cx = sg_.f2c()*Vector(0,0,1); Cx.normalize();

	Size rootres = place_near_origin( pose );

	core::pose::Pose posebase;
	utility::vector1<Size> Ajumps, Bjumps, Cjumps, monomer_jumps, monomer_anchors;
	Size base_monomer;

	Vector max_extent(0,0,0);
	for ( Size i=1; i<= pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		for ( Size j=1; j<=pose.residue(i).natoms(); ++j ) {
			Vector f_ij = sg_.c2f()*pose.residue(i).xyz(j);
			for ( int k=0; k<3; ++k ) max_extent[k] = std::max( std::fabs(f_ij[k]) , max_extent[k] );
		}
	}
	max_extent += sg_.c2f()*Vector(6,6,6);

	// in each direction, find closest symmcenter we are not considering
	numeric::xyzVector< Size > grid = sg_.get_nsubdivisions();
	numeric::xyzVector<core::Real> closest_nongen_symmcenter(1.0/grid[0], 1.0/grid[1], 1.0/grid[2]);

	numeric::xyzVector<int> EXTEND(
		std::max( 1, (int)std::ceil( 2*max_extent[0]-closest_nongen_symmcenter[0] )),
		std::max( 1, (int)std::ceil( 2*max_extent[1]-closest_nongen_symmcenter[1] )),
		std::max( 1, (int)std::ceil( 2*max_extent[2]-closest_nongen_symmcenter[2] )) );

	TR.Debug << "using extent " << EXTEND[0] << "," << EXTEND[1] << "," << EXTEND[2] << std::endl;
	TR.Debug << "with max_extent = " << max_extent[0] << "," << max_extent[1] << "," << max_extent[2] << std::endl;
	TR.Debug << "with closest_nongen_symmcenter = " << closest_nongen_symmcenter[0] << "," << closest_nongen_symmcenter[1] << "," << closest_nongen_symmcenter[2] << std::endl;

	build_lattice_of_virtuals( posebase, EXTEND, Ajumps, Bjumps, Cjumps, monomer_anchors, base_monomer);

	Size nvrt = posebase.size();
	Size nres_monomer = pose.size();

	// only connecting ones!
	// updates monomer_anchors
	detect_connecting_subunits( pose, posebase, monomer_anchors, base_monomer );
	add_monomers_to_lattice( pose, posebase, monomer_anchors, monomer_jumps, rootres );

	Size nsubunits = monomer_anchors.size();

	core::pose::PDBInfoOP pdbinfo_old=pose.pdb_info(), pdbinfo_new;
	pose = posebase;

	conformation::symmetry::SymmetryInfo syminfo;
	setup_xtal_symminfo( pose, nsubunits, nvrt, base_monomer, nres_monomer, Ajumps, Bjumps, Cjumps, monomer_jumps, syminfo );

	bool symmdetectdisulf = basic::options::option[ basic::options::OptionKeys::symmetry::detect_bonds ]();
	basic::options::option[ basic::options::OptionKeys::symmetry::detect_bonds ].value(false);
	pose::symmetry::make_symmetric_pose( pose, syminfo );
	basic::options::option[ basic::options::OptionKeys::symmetry::detect_bonds ].value(symmdetectdisulf);

	pdbinfo_new = pose::PDBInfoOP( new pose::PDBInfo( pose, true ) );
	core::pose::symmetry::make_symmetric_pdb_info( pose, pdbinfo_old, pdbinfo_new );
	pdbinfo_new->set_crystinfo(ci);
	pose.pdb_info( pdbinfo_new );

	//fpd foldtree above assumes no monomer jumps
	//fpd here, we reinitialize the foldtree to let the SymmetryInfo machinery take care of
	//    monomer jumps
	//fpd this also has the advantage of renumbering jumps in a consistent way so that
	//    all symmetric foldtree manipulations work well with symm poses made from this function
	core::kinematics::FoldTree ft = core::conformation::symmetry::get_asymm_unit_fold_tree(pose.conformation());
	core::conformation::symmetry::symmetrize_fold_tree(pose.conformation(), ft);
	pose.fold_tree( ft );

	// force symmetrization
	core::conformation::symmetry::SymmetryInfoCOP symminfo_new =
		dynamic_cast<core::conformation::symmetry::SymmetricConformation const & >( pose.conformation()).Symmetry_Info();
	for ( core::Size i=pose.fold_tree().num_jump(); i>=1; --i ) {
		if ( symminfo_new->jump_is_independent(i) && symminfo_new->jump_clones(i).size() > 0 ) {
			core::kinematics::Jump j_i = pose.jump( i );
			pose.set_jump( i, j_i );
		}
	}

	// if we are mirror symmetric, update restypes
	if ( core::conformation::symmetry::is_mirror_symmetric( pose.conformation() ) ) {
		core::conformation::symmetry::MirrorSymmetricConformation & mirror_conf(
			dynamic_cast< core::conformation::symmetry::MirrorSymmetricConformation& >( pose.conformation() ) );
		mirror_conf.update_residue_identities();
	}

	// update disulf info
	pose.conformation().detect_disulfides();

}

Size
MakeLatticeMover::place_near_origin (
	Pose & pose
) {
	Size rootpos=0;
	Size nres = pose.size();

	Vector com(0,0,0);
	for ( Size i=1; i<= nres; ++i ) {
		com += pose.residue(i).xyz("CA");
		if ( pose.residue(i).is_upper_terminus() ) break;
	}
	com /= nres;

	Real mindis2(1e6);
	for ( Size i=1; i<= nres; ++i ) {
		Real const dis2( com.distance_squared(  pose.residue(i).xyz("CA") ) );
		if ( dis2 < mindis2 ) {
			mindis2 = dis2;
			rootpos = i;
		}
		if ( pose.residue(i).is_upper_terminus() ) break;
	}

	Size nsymm = sg_.nsymmops();
	Size bestxform=0;
	Vector bestoffset(0,0,0);
	mindis2=1e6;
	com = sg_.c2f()*com;
	for ( Size i=1; i<=nsymm; ++i ) {
		Vector foffset = sg_.symmop(i).get_rotation()*com + sg_.symmop(i).get_translation(), rfoffset;
		rfoffset[0] = min_mod( foffset[0], 1.0 );
		rfoffset[1] = min_mod( foffset[1], 1.0 );
		rfoffset[2] = min_mod( foffset[2], 1.0 );
		Real dist = (sg_.f2c()*rfoffset).length_squared();
		if ( dist<mindis2 ) {
			mindis2=dist;
			bestxform=i;
			bestoffset = foffset - rfoffset;
		}
	}

	numeric::xyzMatrix<Real> R = sg_.f2c()*sg_.symmop(bestxform).get_rotation()*sg_.c2f();
	numeric::xyzVector<Real> T = sg_.f2c()*(sg_.symmop(bestxform).get_translation() - bestoffset);
	pose.apply_transform_Rx_plus_v( R,T );

	return rootpos;
}

void
MakeLatticeMover::detect_connecting_subunits(
	Pose const & monomer_pose,
	Pose const & pose,
	utility::vector1<Size> & monomer_anchors,
	Size &basesubunit
) {
	utility::vector1<Size> new_monomer_anchors;
	Size new_basesubunit=0;

	utility::vector1< numeric::xyzMatrix<core::Real> > new_allRs;
	utility::vector1< numeric::xyzVector<core::Real> > new_allTs;


	// get pose radius
	Size const num_monomers( monomer_anchors.size() ), nres_monomer( monomer_pose.size () );

	Vector com(0,0,0);
	Real radius = 0;
	utility::vector1<Vector> monomer_cas(nres_monomer);

	Vector T0 = pose.residue(monomer_anchors[basesubunit]).xyz("ORIG");
	runtime_assert( T0.length() < 1e-6);

	for ( Size i=1; i<= nres_monomer; ++i ) {
		Vector ca_i(0,0,0);
		if ( pose.residue(i).is_protein() ) {
			ca_i = monomer_pose.residue(i).xyz("CA");
		} else {
			// com
			for ( Size j=1; j<= monomer_pose.residue(i).natoms(); ++j ) {
				ca_i += monomer_pose.residue(i).xyz(j);
			}
			ca_i = ca_i/((core::Real)monomer_pose.residue(i).natoms());
		}
		monomer_cas[i] = ca_i;
		radius = std::max( (ca_i).length_squared() , radius );
	}
	radius = sqrt(radius);

	// make master first
	new_monomer_anchors.push_back(monomer_anchors[basesubunit]);
	new_allRs.push_back( allRs_[basesubunit] );
	new_allTs.push_back( allTs_[basesubunit] );
	new_basesubunit = 1;

	for ( Size i=1; i<=num_monomers; ++i ) {
		if ( i==basesubunit ) continue;

		// pass 1 check vrt-vrt dist to throw out very distant things
		Vector T = pose.residue(monomer_anchors[i]).xyz("ORIG");
		Real disVRT = T.length();
		if ( disVRT>contact_dist_+2*radius ) continue;

		// pass 2 check ca-ca dists
		Vector X = pose.residue(monomer_anchors[i]).xyz("X") - T;
		Vector Y = pose.residue(monomer_anchors[i]).xyz("Y") - T;
		Vector Z = X.cross( Y );

		numeric::xyzMatrix<core::Real> R = numeric::xyzMatrix<core::Real>::cols(X,Y,Z);

		bool contact=false;
		for ( Size j=1; j<= nres_monomer && !contact; ++j ) {
			Vector Ri = R*monomer_cas[j] + T;
			for ( Size k=1; k<= nres_monomer && !contact; ++k ) {
				contact = ((Ri-monomer_cas[k]).length_squared() < contact_dist_*contact_dist_);
			}
		}
		if ( contact ) {
			new_monomer_anchors.push_back(monomer_anchors[i]);
			new_allRs.push_back( allRs_[i] );
			new_allTs.push_back( allTs_[i] );
		}
	}
	basesubunit = new_basesubunit;
	monomer_anchors = new_monomer_anchors;
	allRs_ = new_allRs;
	allTs_ = new_allTs;
}


void
MakeLatticeMover::add_monomers_to_lattice(
	Pose const & monomer_pose,
	Pose & pose,
	utility::vector1<Size> const & monomer_anchors,
	utility::vector1<Size> & monomer_jumps,
	Size rootpos
) {
	Size const num_monomers( monomer_anchors.size() ), nres_monomer( monomer_pose.size () );

	monomer_jumps.clear();
	Size n_framework_jumps = pose.fold_tree().num_jump();

	Size nres_protein(0);
	for ( Size i=1; i<= num_monomers; ++i ) {
		//std::cerr << "inserting " << i << " of " << num_monomers << std::endl;
		Size const anchor( monomer_anchors[i] + nres_protein ); // since we've already done some insertions
		Size const old_nres_protein( nres_protein );
		pose.insert_residue_by_jump( monomer_pose.residue(rootpos), old_nres_protein+1, anchor ); ++nres_protein;
		for ( Size j=rootpos-1; j>=1; --j ) {
			pose.prepend_polymer_residue_before_seqpos( monomer_pose.residue(j), old_nres_protein+1, false ); ++nres_protein;
		}
		for ( Size j=rootpos+1; j<= monomer_pose.size(); ++j ) {
			if ( monomer_pose.residue(j).is_lower_terminus() ) {
				pose.insert_residue_by_jump( monomer_pose.residue(j), nres_protein+1, nres_protein ); ++nres_protein;
			} else {
				pose.append_polymer_residue_after_seqpos( monomer_pose.residue(j), nres_protein, false ); ++nres_protein;
			}
		}
		monomer_jumps.push_back( n_framework_jumps+i );
	}
	for ( Size i=1; i<= num_monomers; ++i ) {
		pose.conformation().insert_chain_ending( nres_monomer*i );
	}

	core::kinematics::FoldTree f( pose.fold_tree() );
	f.reorder( nres_protein+1 );
	pose.fold_tree(f);
}


void
MakeLatticeMover::build_lattice_of_virtuals(
	core::pose::Pose & posebase,
	numeric::xyzVector<int> EXTEND,
	utility::vector1<Size> &Ajumps,
	utility::vector1<Size> &Bjumps,
	utility::vector1<Size> &Cjumps,
	utility::vector1<Size> &subunit_anchors,
	Size &basesubunit
) {
	numeric::xyzVector<int> nvrts_dim(2*EXTEND[0]+1, 2*EXTEND[1]+1, 2*EXTEND[2]+1);

	posebase.clear();
	allRs_.clear();
	allTs_.clear();

	// get gridspace in each dimension
	numeric::xyzVector< Size > grid = sg_.get_nsubdivisions();
	numeric::xyzVector< Size > trans_dofs = sg_.get_trans_dofs();

	Vector O, Ax, Bx, Cx, Ay, By, Cy;
	if ( sg_.setting() == HEXAGONAL ) {
		Ax = sg_.f2c()*Vector(1,0,0); Ax.normalize();
		Bx = sg_.f2c()*Vector(0,1,0); Bx.normalize();
		Cx = sg_.f2c()*Vector(0,0,1); Cx.normalize();
	} else {
		Ax = Vector(1,0,0); Ax.normalize();
		Bx = Vector(0,1,0); Bx.normalize();
		Cx = Vector(0,0,1); Cx.normalize();
	}

	// we don't care where y is as long as it is perpendicular
	Ay = Bx - Ax.dot(Bx)*Ax; Ay.normalize();
	By = Cx - Bx.dot(Cx)*Bx; By.normalize();
	Cy = Ax - Cx.dot(Ax)*Cx; Cy.normalize();

	ObjexxFCL::FArray3D<int> vrtX, vrtY, vrtZ;
	vrtX.dimension(nvrts_dim[0]*grid[0]+1,1,1); vrtX=0;
	vrtY.dimension(nvrts_dim[0]*grid[0]+1,nvrts_dim[1]*grid[1]+1,1); vrtY=0;
	vrtZ.dimension(nvrts_dim[0]*grid[0]+1,nvrts_dim[1]*grid[1]+1,nvrts_dim[2]*grid[2]+1); vrtZ=0;

	// 1 expand A (j==k==1)
	for ( int i=1; i<=(int)(nvrts_dim[0]*grid[0]+1); ++i ) {
		Vector fX( (Real)(i-1-(Real)(EXTEND[0]*grid[0]))/(Real)grid[0], -EXTEND[1], -EXTEND[2] );
		O = sg_.f2c()*fX;

		// now add 3 virtuals to the pose, with X pointing toward A,B,C, respectively
		ResidueOP vrt_x = make_vrt(O,Ax,Ay);
		ResidueOP vrt_y = make_vrt(O,Bx,By);
		ResidueOP vrt_z = make_vrt(O,Cx,Cy);

		if ( i==1 ) {
			posebase.append_residue_by_bond( *vrt_x );    vrtX(1,1,1) = 1;
			posebase.append_residue_by_jump( *vrt_y, 1);  vrtY(1,1,1) = 2;
			posebase.append_residue_by_jump( *vrt_z, 2);  vrtZ(1,1,1) = 3;
		} else {
			posebase.append_residue_by_jump( *vrt_x, vrtX(i-1,1,1));  vrtX(i,1,1) = posebase.size();
			Ajumps.push_back(posebase.fold_tree().num_jump());
			posebase.append_residue_by_jump( *vrt_y, vrtX(i  ,1,1));  vrtY(i,1,1) = posebase.size();
			posebase.append_residue_by_jump( *vrt_z, vrtY(i  ,1,1));  vrtZ(i,1,1) = posebase.size();
		}
	}

	// expand B (k==1)
	for ( int i=1; i<=(int)(nvrts_dim[0]*grid[0]+1); ++i ) {
		for ( int j=2; j<=(int)(nvrts_dim[1]*grid[1]+1); ++j ) {
			Vector fX( (i-1-(Real)(EXTEND[0]*grid[0]))/(Real)grid[0], (j-1-(Real)(EXTEND[1]*grid[1]))/(Real)grid[1], -EXTEND[2] );
			O = sg_.f2c()*fX;

			// now add 3 virtuals to the pose, with X pointing toward A,B,C, respectively
			ResidueOP vrt_y = make_vrt(O,Bx,By);
			ResidueOP vrt_z = make_vrt(O,Cx,Cy);

			posebase.append_residue_by_jump( *vrt_y, vrtY(i,j-1,1)); vrtY(i,j,1) = posebase.size();
			Bjumps.push_back(posebase.fold_tree().num_jump());
			posebase.append_residue_by_jump( *vrt_z, vrtY(i,j,1));   vrtZ(i,j,1) = posebase.size();
		}
	}

	// expand C
	for ( int i=1; i<=(int)(nvrts_dim[0]*grid[0]+1); ++i ) {
		for ( int j=1; j<=(int)(nvrts_dim[1]*grid[1]+1); ++j ) {
			for ( int k=2; k<=(int)(nvrts_dim[2]*grid[2]+1); ++k ) {
				Vector fX( (i-1-(Real)(EXTEND[0]*grid[0]))/(Real)grid[0], (j-1-(Real)(EXTEND[1]*grid[1]))/(Real)grid[1], (k-1-(Real)(EXTEND[2]*grid[2]))/(Real)grid[2] );
				O = sg_.f2c()*fX;

				// now add 3 virtuals to the pose, with X pointing toward A,B,C, respectively
				ResidueOP vrt_z = make_vrt(O,Cx,Cy);

				posebase.append_residue_by_jump( *vrt_z, vrtZ(i,j,k-1));  vrtZ(i,j,k) = posebase.size();
				Cjumps.push_back(posebase.fold_tree().num_jump());
			}
		}
	}

	// add "hanging" virtuals
	for ( int s=1; s<=(int)sg_.nsymmops(); ++s ) {
		numeric::xyzMatrix<Real> R_i = sg_.symmop(s).get_rotation(), R_i_cart;
		numeric::xyzVector<Real> T_i = sg_.symmop(s).get_translation();

		// T_i -> indices
		for ( int i=-(int)EXTEND[0]; i<=(int)EXTEND[0]; ++i ) {
			for ( int j=-(int)EXTEND[1]; j<=(int)EXTEND[1]; ++j ) {
				for ( int k=-(int)EXTEND[2]; k<=(int)EXTEND[2]; ++k ) {
					// find lattice anchor
					int x_i = (int)std::floor( (i+EXTEND[0]+T_i[0])*grid[0] + 1.5 );
					int y_i = (int)std::floor( (j+EXTEND[1]+T_i[1])*grid[1] + 1.5 );
					int z_i = (int)std::floor( (k+EXTEND[2]+T_i[2])*grid[2] + 1.5 );

					O = posebase.residue(vrtZ(x_i,y_i,z_i)).xyz("ORIG");

					if ( sg_.setting() == HEXAGONAL ) {
						R_i_cart = sg_.f2c()*R_i*sg_.c2f();
						Ax = R_i_cart*Vector(1,0,0); Ax.normalize();
						Ay = R_i_cart*Vector(0,1,0); Ay.normalize();
					} else {
						R_i_cart = R_i;
						Ax = R_i*Vector(1,0,0); Ax.normalize();
						Ay = R_i*Vector(0,1,0); Ay.normalize();
					}
					ResidueOP vrt_z = make_vrt(O,Ax,Ay, (R_i_cart.det()<0));

					posebase.append_residue_by_jump( *vrt_z, vrtZ(x_i,y_i,z_i));

					allRs_.push_back(R_i);
					allTs_.push_back(numeric::xyzVector<Real>(i+T_i[0],j+T_i[1],k+T_i[2]));

					subunit_anchors.push_back(posebase.size());
					if ( s==1 && i==0 && j==0 && k==0 ) {
						basesubunit = subunit_anchors.size();
					}
				}
			}
		}
	}
}


void
MakeLatticeMover::setup_xtal_symminfo(
	Pose & pose,
	Size const num_monomers,
	Size const num_virtuals,
	Size const base_monomer,
	Size const nres_monomer,
	utility::vector1<Size> const &Ajumps,
	utility::vector1<Size> const &Bjumps,
	utility::vector1<Size> const &Cjumps,
	utility::vector1<Size> const &monomer_jumps,
	conformation::symmetry::SymmetryInfo & symminfo
) {

	// bb clones
	for ( Size i=1; i<= num_monomers; ++i ) {
		if ( i != base_monomer ) {
			Size const offset( (i-1)*nres_monomer ), base_offset( (base_monomer-1)*nres_monomer);
			for ( Size j=1; j<= nres_monomer; ++j ) {
				symminfo.add_bb_clone ( base_offset+j, offset+j );
				symminfo.add_chi_clone( base_offset+j, offset+j );
			}
		}
	}

	// subunit base jump clones
	Size const base_monomer_jump( monomer_jumps[ base_monomer ] );
	for ( Size i=1; i<=monomer_jumps.size(); ++i ) {
		if ( monomer_jumps[i]!= base_monomer_jump ) {
			symminfo.add_jump_clone( base_monomer_jump, monomer_jumps[i], 0.0 );
		}
	}

	// unit cell clones
	using core::conformation::symmetry::SymDof;
	std::map< Size, SymDof > symdofs;

	Size Amaster=Ajumps[1], Bmaster=Bjumps[1], Cmaster=Cjumps[1];
	numeric::xyzVector<core::Size> linked_dofs=sg_.get_trans_dofs();
	if ( linked_dofs[1]==1 ) { Bmaster=Ajumps[1]; /*std::cerr << "clone B->A" << std::endl;*/ }
	if ( linked_dofs[2]==1 ) { Cmaster=Ajumps[1]; /*std::cerr << "clone C->A" << std::endl;*/ }

	for ( Size i=2; i<=Ajumps.size(); ++i ) {
		symminfo.add_jump_clone( Amaster, Ajumps[i], 0.0 );
	}
	for ( Size i=1; i<=Bjumps.size(); ++i ) {
		if ( Bmaster!=Bjumps[i] ) {
			symminfo.add_jump_clone( Bmaster, Bjumps[i], 0.0 );
		}
	}
	for ( Size i=1; i<=Cjumps.size(); ++i ) {
		if ( Cmaster!=Cjumps[i] ) {
			symminfo.add_jump_clone( Cmaster, Cjumps[i], 0.0 );
		}
	}

	SymDof symdof_a;
	SymDof symdof_b;
	SymDof symdof_c;

	if ( refinable_lattice_ ) {
		core::Size nrot_dofs=sg_.get_nrot_dofs();
		if ( nrot_dofs == 3 ) {
			symdof_c.read( "x y z" );
			symdof_b.read( "x z" );
		} else if ( nrot_dofs == 1 ) {
			symdof_c.read( "x y" );
			symdof_b.read( "x" );
		} else {
			symdof_c.read( "x" );
			symdof_b.read( "x" );
		}
		symdof_a.read( "x" );

		symdofs[ Ajumps[1] ] = symdof_a;
		if ( Bmaster==Bjumps[1] ) {
			symdofs[ Bjumps[1] ] = symdof_b;
		}
		if ( Cmaster==Cjumps[1] ) {
			symdofs[ Cjumps[1] ] = symdof_c;
		}
	}

	// jump names
	TR << "Initializing " << pose.num_jump() << " jumps." << std::endl;
	for ( Size v=1; v<=pose.num_jump(); ++v ) symminfo.set_jump_name(v, "v_"+utility::to_string(v));

	symminfo.set_jump_name(Ajumps[1], "A");
	symminfo.set_jump_name(Bjumps[1], "B");
	symminfo.set_jump_name(Cjumps[1], "C");

	SymDof symdof_m;
	symdof_m.read( sg_.get_moveable_dofs()+" angle_x angle_y angle_z");
	symdofs[ base_monomer_jump ] = symdof_m;
	symminfo.set_jump_name(base_monomer_jump, "SUB");

	symminfo.set_dofs( symdofs );

	symminfo.num_virtuals( num_virtuals );
	symminfo.set_use_symmetry( true );
	symminfo.set_flat_score_multiply( pose.size(), 0 );
	symminfo.set_nres_subunit( nres_monomer );

	Size const nres_protein( num_monomers * nres_monomer );
	for ( Size i=1; i<= nres_protein; ++i ) {
		if ( symminfo.bb_is_independent( i ) ) symminfo.set_score_multiply( i, 2 );
		else symminfo.set_score_multiply( i, 1 );
	}
	symminfo.update_score_multiply_factor();
}

// parse_my_tag
void
MakeLatticeMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & ,
	moves::Movers_map const & ,
	core::pose::Pose const & /*pose*/ )
{
	if ( tag->hasOption( "contact_dist" ) ) {
		contact_dist_ = tag->getOption<core::Real>( "contact_dist" );
	}
}

}
}
