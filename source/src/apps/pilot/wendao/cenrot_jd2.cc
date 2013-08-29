// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file cenrot_jd2.cc
/// @brief test centroid rot model
/// @author Yuan Liu


// Core Headers
#include <core/chemical/AA.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// 
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/Minimizer.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/Mover.fwd.hh>

//sampling
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/moves/NullMover.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <protocols/jd2/JobDistributor.hh>

#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <map>

/////////////////////////////////////////////////////////////////
using namespace core;
using namespace core::import_pose;
using namespace core::io::pdb;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core::pose;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::scoring;
using namespace core::kinematics;
using namespace protocols::moves;

using namespace core::pack::task;
using namespace core::pack::dunbrack;
using namespace core::pack::dunbrack::cenrot;

//////////////////////////////////////////////////////////////////
static numeric::random::RandomGenerator RG(62331911);
basic::Tracer TR("pilot.wendao.cenrotjd2");

utility::vector1<core::Size> nrecovery(20);
utility::vector1<core::Size> n_total(20);

void* my_main( void* );

///////////////////////////////////////////////////////////////////
OPT_KEY(Boolean, fit_best_rotamer)
OPT_KEY(Boolean, switch_to_centroid)
OPT_KEY(Boolean, input_cenrot_pdb)
OPT_KEY(Boolean, output_cenrot_intcoord)
OPT_KEY(Boolean, output_bestrot_err)
OPT_KEY(Boolean, output_cenrot_pdb)
OPT_KEY(Boolean, repack_cenrot)
OPT_KEY(Boolean, relax_cenrot)
OPT_KEY(Boolean, opt_after_relax)
OPT_KEY(Integer, repack_buried_cutoff)
OPT_KEY(Real, repack_bfactor_cutoff)
OPT_KEY(Integer, repack_ncycle)
OPT_KEY(String, output_cenrot_score)
OPT_KEY(String, output_cenrot_dir)
OPT_KEY(String, output_cenrot_prefix)

OPT_KEY(String, cenrot_score)

OPT_KEY(Real, relax_temp)
OPT_KEY(Integer, relax_step_per_cycle)
OPT_KEY(Integer, relax_cycle_number)

OPT_1GRP_KEY(Boolean, min, cenrot) //add min mover alone
OPT_1GRP_KEY(Boolean, min, debug)
OPT_1GRP_KEY(Boolean, min, cartesian)

int main( int argc, char * argv [] ) {
	NEW_OPT(fit_best_rotamer, "fit the exact centroid to the closest in lib", false);
	NEW_OPT(switch_to_centroid, "switch the fa pdb to the old centroid", false);
	NEW_OPT(output_cenrot_intcoord, "output the internal coordinate, for building lib", false);
	NEW_OPT(output_bestrot_err, "output the distance between the native rot and best fitted one", false);
	NEW_OPT(input_cenrot_pdb, "input centroid pdbs for scoring", false);
	NEW_OPT(output_cenrot_pdb, "output centroid pdbs for building database", false);
	NEW_OPT(output_cenrot_score, "score the centroid pdbs", "cenrot_score.out");
	NEW_OPT(output_cenrot_dir, "dir for output centroid pdbs", ".");
	NEW_OPT(output_cenrot_prefix, "prefix for pdbs", "idealized_");
	NEW_OPT(repack_cenrot, "repack the centroid rotamer model", false);
	NEW_OPT(repack_bfactor_cutoff, "count repack side-chain with Bfactor lower than default 100", 100.0);
	NEW_OPT(repack_buried_cutoff, "count repack side-chain with buried cutoff", 0);
	NEW_OPT(repack_ncycle, "how many times to repack", 1);
	NEW_OPT(cenrot_score, "cenrot score weight file", "test.wts");

	NEW_OPT(relax_cenrot, "relax the centroid rotamer model", false);
	NEW_OPT(opt_after_relax, "opt after relax", false);
	NEW_OPT(relax_temp, "temp", 1.0);
	NEW_OPT(relax_step_per_cycle, "step", 100);
	NEW_OPT(relax_cycle_number, "cycle", 1);

	NEW_OPT(min::cenrot, "add min mover alone", false);
	NEW_OPT(min::debug, "debug derivs?", false);
	NEW_OPT(min::cartesian, "cartesian minimization?", false);

	try {
		devel::init(argc, argv);
		protocols::viewer::viewer_main( my_main );
	}
	catch ( utility::excn::EXCN_Base const & e) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}

///////////////////////////////////////////////////////////////////
class OutputCenrotIntCoord : public protocols::moves::Mover {
public:
	OutputCenrotIntCoord() {}
	virtual std::string get_name() const {return "OutputCenrotIntCoordMover";}

	void apply( Pose & p ) {
		for ( core::Size i=1; i<= p.total_residue(); ++i )
		{
			Residue const & rsd( p.residue(i) );
			/// for each residue, find out the dof-id of CEN
			id::DOF_ID id_dis(id::AtomID(p.residue(i).atom_index("CEN"), i), id::D);
			id::DOF_ID id_ang(id::AtomID(p.residue(i).atom_index("CEN"), i), id::THETA);
			id::DOF_ID id_dih(id::AtomID(p.residue(i).atom_index("CEN"), i), id::PHI);

			//output internal coordinates of centroids
			//as well as phi/psi angle 
			if ( !rsd.is_terminus() ) {
				TR << "CEN-INT: " << rsd.name3() << " " << i << " "
				<< p.dof(id_dis) << " "
				<< p.dof(id_ang) << " "
				<< p.dof(id_dih) << " "
				<< p.psi(i) << " "
				<< p.phi(i) << std::endl;
			}
		}
	}
};

class RescoreCenrot : public protocols::moves::Mover {
private:
	core::scoring::ScoreFunctionOP scfxn_;
	pose::PoseOP native_pose_;

public:
	RescoreCenrot(){}
	virtual std::string get_name() const {return "RescoreCenrotMover";}

	void set_scfxn(core::scoring::ScoreFunctionOP const &in_score) {
		scfxn_ = in_score;
	}
	void set_native(pose::PoseOP const &in_pose) {
		native_pose_ = in_pose;
	}

	void apply( Pose & p ) {
		//setup the ss correctly
		core::scoring::dssp::Dssp dssp( p );
		dssp.insert_ss_into_pose( p );

		(*scfxn_)(p);

		// get rmsd
		//if (native_pose_) {
		//	Real rmsd = core::scoring::CA_rmsd(p, *native_pose_);
		//	setPoseExtraScores(p, "rmsd", rmsd);
		//}
	}
};

class RepackCenrotMover : public protocols::moves::Mover {
private:
	core::scoring::ScoreFunctionOP scfxn_;

public:
	virtual std::string get_name() const {return "RepackCenrotMover";}
	void set_scfxn(core::scoring::ScoreFunctionOP const &in_score) {
		scfxn_ = in_score;
	}

	void apply( Pose & p ) {
		using namespace core::pack::task;
		TaskFactoryOP main_task_factory = new TaskFactory;
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory->push_back( rtrop );

		protocols::simple_moves::PackRotamersMover packrotamersmover;
		packrotamersmover.task_factory(main_task_factory);
		packrotamersmover.score_function(scfxn_);

		packrotamersmover.apply(p);

		p.dump_pdb("repacked.pdb");
	}
};

class MinCenrotMover : public protocols::moves::Mover {
public:
	MinCenrotMover(){}
	virtual std::string get_name() const { return "MinCenrotMover"; }

	void apply( core::pose::Pose & pose) {
		using namespace core;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::id;

		bool debug_derivs = option[ OptionKeys::min::debug ]();

		//movemap
		core::kinematics::MoveMap mm;
		mm.set_bb  ( false ); mm.set_chi ( true ); mm.set_jump( false );

		// min option
		core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 1e-3, true, debug_derivs, debug_derivs );
		minoptions.nblist_auto_update( true );
		minoptions.max_iter( 100 );

		// scorefunction1
		core::scoring::ScoreFunctionOP scorefxn;
		if ( option[ score::weights ].user() ) {
			scorefxn = core::scoring::getScoreFunction();
		}
		else {
			scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("score_cenrot");
		}

		if (option[OptionKeys::min::cartesian]()) {
			core::optimization::CartesianMinimizer minimizer;

			minimizer.run( pose, mm, *scorefxn, minoptions );
		}
		else {
			core::optimization::AtomTreeMinimizer minimizer;

			//setup movemap for sidechain
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				if ( res_i.aa()!=chemical::aa_gly 
				     && res_i.aa()!=chemical::aa_ala 
				     && res_i.type().has_atom_name( "CEN")) {
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::D ), true );
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::THETA ), true );
				}
			}

			minimizer.run( pose, mm, *scorefxn, minoptions );
		}

		scorefxn->show( TR.Debug, pose  );
		pose.dump_pdb("minimized.pdb");
	}
};

class CenRotRelaxMover : public protocols::moves::Mover {
public:
	CenRotRelaxMover(){}
	virtual std::string get_name() const { return "CenRotRelaxMover"; }

	void apply( core::pose::Pose & pose) {
		using namespace core;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::id;

		bool debug_derivs = option[ OptionKeys::min::debug ]();

		////////////////////////
		//general setup
		////////////////////////

		//movemap
		core::kinematics::MoveMap mm;
		mm.set_bb  ( true ); mm.set_chi ( true ); mm.set_jump( true );

		// min option
		core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 1e-2, true, debug_derivs, debug_derivs );
		minoptions.nblist_auto_update( true );
		minoptions.max_iter( 100 );

		// scorefunction1
		core::scoring::ScoreFunctionOP scorefxn1;
		if ( option[ score::weights ].user() ) {
			scorefxn1 = core::scoring::getScoreFunction();
		}
		else {
			scorefxn1 = core::scoring::ScoreFunctionFactory::create_score_function("score_cenrot");
		}

		// ramp vdw term
		core::Real max_vdw = scorefxn1->get_weight( core::scoring::vdw );
		int ncyc = option[ relax::default_repeats ]();

		//repack
		RepackCenrotMover repack;
		repack.set_scfxn(scorefxn1);
		repack.apply(pose);

		//min type
		if (option[OptionKeys::min::cartesian]()) {
			core::optimization::CartesianMinimizer minimizer;

			// only for cart_min, set cart_bonded
			// scorefunction0 -- fix bad chainbreaks
			core::scoring::ScoreFunctionOP scorefxn0 = core::scoring::getScoreFunction();
			scorefxn0->reset();
			scorefxn0->set_weight( core::scoring::vdw, 0.1 );
			scorefxn0->set_weight( core::scoring::cart_bonded, 0.1 );
			core::scoring::methods::EnergyMethodOptions options0(scorefxn0->energy_method_options());
			options0.set_cartesian_bonded_linear(true);
			options0.set_cartesian_bonded_parameters(10,1,0,0,0);
			scorefxn0->set_energy_method_options(options0);

			if ( scorefxn1->get_weight(core::scoring::cart_bonded)==0 ) {
				scorefxn1->set_weight( core::scoring::cart_bonded, 0.1 );
			}
			core::scoring::methods::EnergyMethodOptions options1(scorefxn1->energy_method_options());
			options1.set_cartesian_bonded_linear(true);
			scorefxn1->set_energy_method_options(options1);

			minimizer.run( pose, mm, *scorefxn0, minoptions );
			TR << "cart_min " << ncyc << std::endl;

			for (int i=1; i<=ncyc; ++i) {
				scorefxn1->set_weight( core::scoring::vdw, max_vdw/(pow(2.0,Real(ncyc-i))) );
				(*scorefxn1)(pose);
				repack.apply(pose);
				minimizer.run( pose, mm, *scorefxn1, minoptions );
				scorefxn1->show( TR.Debug, pose  );
			}
		}
		else {
			core::optimization::AtomTreeMinimizer minimizer;

			//setup movemap for sidechain
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				if ( res_i.aa()!=chemical::aa_gly 
				     && res_i.aa()!=chemical::aa_ala 
				     && res_i.type().has_atom_name( "CEN")) {
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::D ), true );
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::THETA ), true );
				}
			}

			if ( scorefxn1->get_weight(core::scoring::cart_bonded)>0 ) {
				scorefxn1->set_weight( core::scoring::cart_bonded, 0.0 );
			}

			TR << "atomtree_min " << ncyc << std::endl;

			for (int i=1; i<=ncyc; ++i) {
				scorefxn1->set_weight( core::scoring::vdw, max_vdw/(pow(2.0,Real(ncyc-i))) );
				(*scorefxn1)(pose);
				repack.apply(pose);
				minimizer.run( pose, mm, *scorefxn1, minoptions );
				scorefxn1->show( TR.Debug, pose  );
			}
		}
	}
};


///////////////////////////////////////////////////////////////////
void*
my_main( void* ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::moves;
	using namespace protocols::simple_moves;

	SequenceMoverOP do_cenrot( new SequenceMover() );

	//switch
	if (!option[input_cenrot_pdb]() && option[in::file::residue_type_set]()!="CENTROID_ROT") {
		if (option[switch_to_centroid]()) {
			do_cenrot->add_mover(new SwitchResidueTypeSetMover("centroid"));
		}
		else {
			TR << "Switch to CenRot model" << std::endl;
			do_cenrot->add_mover(new SwitchResidueTypeSetMover("centroid_rot"));
		}
	}

	//output intcoord
	if (option[output_cenrot_intcoord]()) {
		do_cenrot->add_mover(new OutputCenrotIntCoord());
	}
	else {
		//setup the scorefxn
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::getScoreFunction();

		//repack
		if (option[repack_cenrot]()) {
			utility::pointer::owning_ptr< RepackCenrotMover > repack = new RepackCenrotMover();
			repack->set_scfxn(score_fxn);
			do_cenrot->add_mover(repack);
		}

		if (option[min::cenrot]()) {
			utility::pointer::owning_ptr< MinCenrotMover > mincenrot = new MinCenrotMover();
			do_cenrot->add_mover(mincenrot);
		}

		//relax
		if (option[relax_cenrot]()) {
			utility::pointer::owning_ptr< CenRotRelaxMover > relax = new CenRotRelaxMover();
			do_cenrot->add_mover(relax);
		}

		//rescore anyway
		utility::pointer::owning_ptr< RescoreCenrot > rescore = new RescoreCenrot();
		rescore->set_scfxn(score_fxn);
		if (option[in::file::native].user()) {
			ResidueTypeSetCAP rsd_set;
			rsd_set=ChemicalManager::get_instance()->residue_type_set( "centroid" );
			core::pose::PoseOP native_pose;
			native_pose = new Pose();
			core::import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
			rescore->set_native(native_pose);
		}
		do_cenrot->add_mover(rescore);
	}

	protocols::jd2::JobDistributor::get_instance()->go( do_cenrot );

	return 0;
}


