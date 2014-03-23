// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file fa_to_cenrot_score.cc
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
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
//#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/scoring/dssp/Dssp.hh>

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
basic::Tracer TR("pilot.wendao.cenrot");
std::map<std::string, Real> masslst;

utility::vector1<core::Size> nrecovery(20);
//utility::vector1<core::Real> err_buried(20);
utility::vector1<core::Size> n_total(20);

void *my_main( void* );
//void switch_to_residue_type_set_cenrot( core::pose::Pose & pose );
void fit_centroid_to_the_best_rot( core::pose::Pose & pose, utility::vector1<Size> &rotlist );
void process_the_pose(core::pose::PoseOP &native_pose, 
	core::pose::Pose & p, std::string const &tag);
void relax_cenrot_pose(core::pose::PoseOP &native_pose, 
	core::pose::Pose & p, std::string const &tag);
void rescore_pose(
	core::scoring::ScoreFunctionOP score_fxn,
	core::pose::PoseOP &native_pose,
	core::pose::Pose & p,
	std::string const &tag);

///////////////////////////////////////////////////////////////////
// class SwitchResidueTypeSetCenrot : public protocols::moves::Mover {
// public:
// 	SwitchResidueTypeSetCenrot(){}
// 	void apply (core::pose::Pose & pose) {
// 		//clear old energy cache
// 		pose.energies().clear();
// 		//get restype
// 		ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( "centroid_rot" );

// 		//loop for each residue
// 		for ( core::Size i=1; i<= pose.total_residue(); ++i ) {
// 			//get the residue
// 			Residue const & rsd( pose.residue(i) );
// 			//check current restype
// 			std::string const & current_type_set_name ( rsd.type().residue_type_set().name() );
// 			if ( current_type_set_name == "centroid_rot" ) {
// 				TR.Warning << "core::util::switch_to_residue_type_set: residue " << i 
// 				<< " already in centroid_rot residue_type_set" << std::endl;
// 				continue;
// 			}

// 			//get temperature
// 			core::Real maxB=0.0;
// 			for (Size k=rsd.first_sidechain_atom(); k<=rsd.nheavyatoms(); ++k) {
// 				if (rsd.is_virtual(k)) continue;
// 				Real B = pose.pdb_info()->temperature( i, k );
// 				maxB = std::max( B, maxB );
// 			}
// 			if (maxB==0.0) {maxB = pose.pdb_info()->temperature( i, rsd.atom_index("CA") );}

// 			//gen new residue
// 			core::conformation::ResidueOP new_rsd( 0 );
// 			if( ( rsd.aa() == aa_unk ) ){
// 				//skip
// 			}
// 			else if ( rsd.name().substr(0,3)=="CYD" ) {
// 				core::chemical::ResidueTypeCOPs const & rsd_types( rsd_set->name3_map( "CYS" ) );
// 				for (core::Size j=1; j<=rsd_types.size(); ++j ) {
// 					core::chemical::ResidueType const & new_rsd_type( *rsd_types[j] );
// 					if ( new_rsd_type.name3()=="CYS" ) {
// 						new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
// 						break;
// 					}
// 				}
// 			}
// 			else if ( rsd.name().substr(0,5)=="HIS_D" ) {
// 				core::chemical::ResidueTypeCOPs const & rsd_types( rsd_set->name3_map( rsd.name3() ) );
// 				for (core::Size j=1; j<=rsd_types.size(); ++j ) {
// 					core::chemical::ResidueType const & new_rsd_type( *rsd_types[j] );
// 					if ( new_rsd_type.name3()=="HIS" ) {
// 						new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
// 						break;
// 					}
// 				}
// 			}
// 			else  if (rsd.is_terminus()) {
// 				//get the terminal type (maybe no need, but to consist with reading a cenrot pdb)
// 				core::chemical::ResidueType const & new_rsd_type( rsd_set->name_map(rsd.name()) );
// 				new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
// 			}
// 			else {
// 				//just find the standard aa restype
// 				core::chemical::ResidueType const & new_rsd_type( rsd_set->name_map(rsd.name3()) );
// 				new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
// 			}

// 			if ( ! new_rsd ) {
// 				TR.Warning << "Did not find perfect match for residue: "  << rsd.name()
// 				<< " at position " << i << ". Trying to find acceptable match. " << std::endl;
				
// 				core::chemical::ResidueTypeCOPs const & rsd_types( rsd_set->name3_map( rsd.name3() ) );
// 				for ( core::Size j=1; j<= rsd_types.size(); ++j ) {
// 					core::chemical::ResidueType const & new_rsd_type( *rsd_types[j] );
// 					if ( rsd.type().name3()  == new_rsd_type.name3()  ) {
// 						new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
// 						break;
// 					}
// 				}
// 				if (  new_rsd ) {
// 					TR.Warning << "Found an acceptable match: " << rsd.type().name() << " --> " << new_rsd->name() << std::endl;
// 				}
// 				else {
// 					//bug here?
// 					utility_exit_with_message( "switch_to_cenrot_residue_type_set fails\n" );
// 				}
// 			}

// 			//find the centroid postion
// 			PointPosition cenrotxyz(0,0,0);
// 			Real mass = 0.0;
// 			if (rsd.name3()=="GLY") {
// 				cenrotxyz = rsd.atom("CA").xyz();
// 			}
// 			else if (rsd.name3()=="ALA") {
// 				cenrotxyz = rsd.atom("CB").xyz();
// 			}
// 			else {
// 				//std::cout<<rsd.name()<<std::endl;
// 				for (	Size na=rsd.type().first_sidechain_atom(); 
// 						na<=rsd.type().nheavyatoms();
// 						na++) {
// 					if (rsd.atom_name(na)==" CB ") continue;
// 					std::string elem = rsd.atom_name(na).substr(1,1);
// 					cenrotxyz += (rsd.atoms()[na].xyz()*masslst[elem]);
// 					mass += masslst[elem];
// 					TR.Debug << "|" << rsd.atom_name(na) << "|"
// 					<< " (" << masslst[elem] << ") "
// 					<< rsd.atoms()[na].xyz().x() << ","
// 					<< rsd.atoms()[na].xyz().y() << ","
// 					<< rsd.atoms()[na].xyz().z() << std::endl;
// 				}
// 				cenrotxyz = cenrotxyz/mass;
// 			}

// 			//replace
// 			if ( ! new_rsd ) {
// 				std::cerr << pose.sequence() << std::endl;
// 				std::cerr  << "can not find a residue type that matches the residue " << rsd.name()
// 				<< " at position " << i << std::endl;
// 				utility_exit_with_message( "switch_to_cenrot_residue_type_set fails\n" );
// 				//continue;
// 			}

// 			//replace it
// 			pose.replace_residue( i, *new_rsd, false );
// 			//set centroid_rot xyz
// 			pose.set_xyz(id::AtomID(pose.residue(i).atom_index("CEN"), i), cenrotxyz);
// 			//set temperature
// 			pose.pdb_info()->temperature(i, pose.residue(i).nbr_atom(), maxB);
// 		}
// 	}

// 	virtual std::string get_name() const {
// 		return "SwitchResidueTypeSetCenrot";
// 	}
// };

///////////////////////////////////////////////////////////////////
OPT_KEY(Boolean, fit_best_rotamer)
OPT_KEY(Boolean, switch_to_centroid)
OPT_KEY(Boolean, input_cenrot_pdb)
OPT_KEY(Boolean, output_cenrot_intcoord)
OPT_KEY(Boolean, output_bestrot_err)
OPT_KEY(Boolean, output_cenrot_pdb)
OPT_KEY(Boolean, repack_cenrot)
OPT_KEY(Boolean, repack_high_vdw)
OPT_KEY(Boolean, min_after_repack)
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

int main( int argc, char * argv [] ) {
	NEW_OPT(fit_best_rotamer, "fit the exact centroid to the closest in lib", false);
	NEW_OPT(switch_to_centroid, "switch the fa pdb to the old centroid", false);
	NEW_OPT(output_cenrot_intcoord, "output the internal coordinate, for building lib", false);
	NEW_OPT(output_bestrot_err, "output the distance between the native rot and best fitted one", false);
	NEW_OPT(input_cenrot_pdb, "input centroid pdbs for scoring", false);
	NEW_OPT(output_cenrot_pdb, "output centroid pdbs for building database", false);
	NEW_OPT(output_cenrot_score, "score the centroid pdbs", "cenrot_score.out");
	NEW_OPT(cenrot_score, "cenrot score weight file", "test.wts");
	NEW_OPT(output_cenrot_dir, "dir for output centroid pdbs", ".");
	NEW_OPT(output_cenrot_prefix, "prefix for pdbs", "idealized_");
	NEW_OPT(repack_cenrot, "repack the centroid rotamer model", false);
	NEW_OPT(repack_high_vdw, "repack the centroid rotamer model with high vdw (4.0)", false);
	NEW_OPT(min_after_repack, "min after relack", false);
	NEW_OPT(repack_bfactor_cutoff, "count repack side-chain with Bfactor lower than default 100", 100.0);
	NEW_OPT(repack_buried_cutoff, "count repack side-chain with buried cutoff", 0);
	NEW_OPT(repack_ncycle, "how many times to repack", 1);

	NEW_OPT(relax_cenrot, "relax the centroid rotamer model", false);
	NEW_OPT(opt_after_relax, "opt after relax", false);
	NEW_OPT(relax_temp, "temp", 1.0);
	NEW_OPT(relax_step_per_cycle, "step", 100);
	NEW_OPT(relax_cycle_number, "cycle", 1);

	devel::init(argc, argv);

	//element mass
	masslst["C"]=12.0107;
	masslst["O"]=15.9994;
	masslst["S"]=32.066;
	masslst["N"]=14.00674;
	masslst["H"]=1.00794;

	protocols::viewer::viewer_main( my_main );
	return 0;
}

////////////////////////////////////////////////////////////////////
void * my_main( void* ) {
	ResidueTypeSetCAP rsd_set;
	if (option[input_cenrot_pdb])
		rsd_set=ChemicalManager::get_instance()->residue_type_set( "centroid_rot" );
	else {
		rsd_set=ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	}
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
	//SwitchResidueTypeSetCenrot to_cenrot;
	protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot("centroid_rot");
	
	//load native pdb
	PoseOP native_pose;
	if (option[in::file::native].user()) {
		native_pose = new Pose();
		core::import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
	}

	if (option[in::file::silent].user()) {
		//TR.Error << "No pdb file found, use -in::file::l to specify a list of pdb!" << std::endl;
		// try load a silent file
		ResidueTypeSetCAP fa_rsd_set;
		fa_rsd_set=ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

		core::io::silent::SilentFileData silent_file_data_in;
		std::string const infile  = option[in::file::silent][1];

		TR << "Reading Silent file" << std::endl;
		silent_file_data_in.read_file( infile );

		for ( core::io::silent::SilentFileData::iterator 
			iter = silent_file_data_in.begin(), 
			end = silent_file_data_in.end(); 
			iter != end; ++iter ) {
			Pose p;
			iter->fill_pose( p, *fa_rsd_set ); // reading fa silent structure as decoys
			std::string const tag = iter->decoy_tag();
			clearPoseExtraScores( p ); // clean all the old score

			if (option[switch_to_centroid]) {
				//for centroid model, rescore only
				to_centroid.apply(p);
				// core::scoring::ScoreFunctionOP score_fxn = new core::scoring::ScoreFunction();
				// score_fxn->set_weight(core::scoring::pair, 1.0);
				// score_fxn->set_weight(core::scoring::env, 1.0);
				//core::scoring::ScoreFunctionOP score_fxn = core::scoring::getScoreFunction();
				core::scoring::ScoreFunctionOP score_fxn = core::scoring::ScoreFunctionFactory::create_score_function(option[cenrot_score]);
				rescore_pose(score_fxn, native_pose, p, tag);
			}
			else {
				//switch_to_residue_type_set_cenrot(p);
				to_cenrot.apply(p);
				process_the_pose(native_pose, p, tag);
			}
		}
	}

	if (option[ in::file::l ].user()) {
		Size npdbs = option[ in::file::l ]().size();
		// -in:file:l
		for (Size npdb=1; npdb<=npdbs; npdb++)
		{
			PoseOP pose = new Pose();
			Pose &p(*pose);
			pose_from_pdb( p, *rsd_set, option[ in::file::l ]()[npdb] );

			//std::cerr << option[ in::file::l ]()[npdb] << " start..." << std::endl;
			if (option[switch_to_centroid]) {
				to_centroid.apply(p);
				// core::scoring::ScoreFunctionOP score_fxn = new core::scoring::ScoreFunction();
				// score_fxn->set_weight(core::scoring::pair, 1.0);
				// score_fxn->set_weight(core::scoring::env, 1.0);
				//core::scoring::ScoreFunctionOP score_fxn = core::scoring::getScoreFunction();
				core::scoring::ScoreFunctionOP score_fxn = core::scoring::ScoreFunctionFactory::create_score_function(option[cenrot_score]);

				rescore_pose(score_fxn, native_pose, p, option[ in::file::l ]()[npdb]);
				continue;
			}
			
			if (!option[input_cenrot_pdb]) to_cenrot.apply(p); //switch_to_residue_type_set_cenrot(p);

			if (option[relax_cenrot]) {
				relax_cenrot_pose(native_pose, p, option[ in::file::l ]()[npdb]);
			}
			else {
				process_the_pose(native_pose, p, option[ in::file::l ]()[npdb]);
			}
		}
	}

	//output for repack
	//if "fit" output recovery rate
	//else output distance errors
	if(option[repack_cenrot]) {
		// if (option[fit_best_rotamer]) {
		// 	TR << "Recovery Rate: " << std::endl;
		// 	for ( int i=core::chemical::aa_ala;
		// 		i <= core::chemical::num_canonical_aas; i++ ) {
		// 		core::chemical::AA aa=static_cast<core::chemical::AA>(i);
		// 		TR << aa << ": " << std::fixed 
		// 		<< std::setw(6) << std::setprecision(3)
		// 		<< core::Real(nrecovery[aa])/n_total[aa]*100
		// 		<< "% of " << n_total[aa] << std::endl;
		// 	}
		// }
		// else {
		// 	TR << "Repacking error: " << std::endl;
		// 	for ( int i=core::chemical::aa_ala;
		// 		i <= core::chemical::num_canonical_aas; i++ ) {
		// 		core::chemical::AA aa=static_cast<core::chemical::AA>(i);
		// 		TR << aa << ": " << err_buried[aa]/n_total[aa] << std::fixed 
		// 		<< std::setw(4) << std::setprecision(3) << std::endl;
		// 	}
		// }

		Size sum_re=0, sum_all=0;
		TR << "Recovery Rate: " << std::endl;
		for ( int i=core::chemical::aa_ala;
			i <= core::chemical::num_canonical_aas; i++ ) {
			core::chemical::AA aa=static_cast<core::chemical::AA>(i);
			if (n_total[aa]==0) continue;
			TR << aa << ": " << std::fixed 
			<< std::setw(6) << std::setprecision(3)
			<< core::Real(nrecovery[aa])/n_total[aa]*100
			<< " % of " << n_total[aa] << std::endl;
			sum_re += nrecovery[aa];
			sum_all += n_total[aa];
		}
		TR << "_ALL_: " << std::fixed 
			<< std::setw(6) << std::setprecision(3)
			<< core::Real(sum_re)/sum_all*100
			<< " % of " << sum_all << std::endl;
	}

	return 0;
}

void relax_cenrot_pose(core::pose::PoseOP &native_pose, core::pose::Pose & p, std::string const &tag)
{
	using namespace core::pack::task;
	using namespace protocols;

	pose::Pose native_p(p);

	//score function
	core::scoring::ScoreFunctionOP score_bb_fxn;
	core::scoring::ScoreFunctionOP score_sc_fxn;
	//score_bb_fxn = core::scoring::getScoreFunction();
	score_bb_fxn = core::scoring::ScoreFunctionFactory::create_score_function(option[cenrot_score]);
	//score_bb_fxn->set_weight(core::scoring::cen_rot_dun, 0.6);
	//score_sc_fxn = core::scoring::getScoreFunction();
	score_sc_fxn = core::scoring::ScoreFunctionFactory::create_score_function(option[cenrot_score]);
	//score_sc_fxn->set_weight(core::scoring::cen_rot_dun, 0.6);

	//repack mover
	TaskFactoryOP main_task_factory = new TaskFactory;
	operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
	main_task_factory->push_back( rtrop );
	protocols::simple_moves::PackRotamersMoverOP packrotamersmover(new protocols::simple_moves::PackRotamersMover());
	packrotamersmover->task_factory(main_task_factory);
	packrotamersmover->score_function(score_sc_fxn);

	//gaussian mover
	simple_moves::BBG8T3AMoverOP bbgmover(new protocols::simple_moves::BBG8T3AMover());

	//monte carlo
	MonteCarloOP mc = new MonteCarlo(p, *score_bb_fxn, option[relax_temp]);

	//combo
	moves::SequenceMoverOP combo( new moves::SequenceMover() );
	combo->add_mover(bbgmover);
	combo->add_mover(packrotamersmover);
	moves::TrialMoverOP trial ( new moves::TrialMover(combo, mc) );
	moves::RepeatMoverOP run( new moves::RepeatMover(trial, option[relax_step_per_cycle]) );

	protocols::simple_moves::MinMoverOP minmover = new protocols::simple_moves::MinMover();
	minmover->score_function(*score_bb_fxn);
	minmover->min_type("lbfgs_armijo_nonmonotone");
	minmover->tolerance(1e-6);
	core::kinematics::MoveMapOP final_mm = new core::kinematics::MoveMap();
	final_mm->set_bb( true );
	minmover->movemap(final_mm);

	moves::SequenceMoverOP combo2( new moves::SequenceMover() );
	combo2->add_mover(packrotamersmover);
	combo2->add_mover(minmover);

	//this is not real relax
	for (int i=1; i<=option[relax_cycle_number]; i++) {
		run->apply(p);

		score_bb_fxn->show(TR,p);
		TR << "score: " << (*score_bb_fxn)(p) << " rmsd: " << core::scoring::CA_rmsd(p, native_p) << std::endl;
		TR.flush();

		std::ostringstream outputfn;
		mc->show_counters();
		mc->reset_counters();

		if (option[output_cenrot_pdb]) {
			outputfn << "traj_" << i << ".pdb";
			p.dump_pdb(outputfn.str());
		}

		if (option[opt_after_relax]) {
			pose::Pose bestpose(mc->lowest_score_pose());
			mc->reset(bestpose);
			p = bestpose;

			moves::RepeatMover(combo2, 3).apply(bestpose);
			(*score_bb_fxn)(bestpose);
			score_bb_fxn->show(TR,bestpose);
			TR.flush();

			std::ostringstream outputmin;
			if (option[output_cenrot_pdb]) {
				outputmin << "minimized_" << i << ".pdb";
				bestpose.dump_pdb(outputmin.str());
			}
		}
	}
}

//should be cenrot
void process_the_pose(core::pose::PoseOP &native_pose, core::pose::Pose & p, std::string const &tag) {
	using namespace core::id;

	core::scoring::ScoreFunctionOP score_fxn;
	if ( option[cenrot_score].user() ) {
		//score_fxn = core::scoring::getScoreFunction();
	  score_fxn = core::scoring::ScoreFunctionFactory::create_score_function(option[cenrot_score]);
	}
	else {
		score_fxn = new core::scoring::ScoreFunction();
		score_fxn->set_weight(core::scoring::cen_rot_pair, 1.0);
		score_fxn->set_weight(core::scoring::cen_rot_env, 1.0);
	}

	(*score_fxn)(p); //score, update nbr, fixing the tenA graph bug
	if (option[output_cenrot_score].user()) {
		rescore_pose( score_fxn, native_pose, p, tag);
	}

	if (option[output_cenrot_intcoord]) {
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

	//save close native rot
	utility::vector1<Size> old_rotlist;
	if (option[fit_best_rotamer]) fit_centroid_to_the_best_rot(p, old_rotlist);

	if (option[repack_cenrot]) {
		//save the old pose
		Pose::Pose ref(p);

		//repack
		using namespace core::pack::task;
		TaskFactoryOP main_task_factory = new TaskFactory;
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory->push_back( rtrop );
		
		// C-beta atoms should not be altered during packing because branching atoms are optimized
		// main_task_factory->push_back( new operation::PreserveCBeta );

		core::scoring::ScoreFunctionOP sf_pack = score_fxn->clone();

		if (option[repack_high_vdw]) {
			sf_pack->set_weight(core::scoring::vdw, 4.0);
		}

		protocols::simple_moves::PackRotamersMover packrotamersmover;
		packrotamersmover.task_factory(main_task_factory);
		packrotamersmover.score_function(sf_pack);

		for (int np=0; np<option[repack_ncycle]; np++) {
			p=ref;
			packrotamersmover.apply(p);

			if (option[min_after_repack]) {
				core::kinematics::MoveMap mm;
				mm.set_bb  ( false ); 
				mm.set_chi ( false );

				Size const n_res( p.n_residue() );
				for (Size i=1; i<=n_res; i++) {
					core::conformation::Residue const &res_i = p.residue(i);
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), i ), core::id::RB1 ), true );
				}

				core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 1e-2, true, false, false );
				minoptions.max_iter( 20 );
				core::optimization::CartesianMinimizer minimizer;

				minimizer.run( p, mm, *score_fxn, minoptions );
			}

			//rescore
			if (option[output_cenrot_score].user()) {
				rescore_pose( score_fxn, native_pose, p, "packed_"+tag);
			}

			for ( core::Size i=1; i<= ref.n_residue(); ++i ) {
				if (ref.pdb_info()->temperature(i, p.residue(i).nbr_atom()) > option[repack_bfactor_cutoff]) continue;
				if (ref.energies().tenA_neighbor_graph().get_node(i)->num_edges() < Size(option[repack_buried_cutoff])) continue;

				Real dis= ( p.residue(i).atom("CEN").xyz()-ref.residue(i).atom("CEN").xyz()).length();
				Real r = (ref.residue(i).atom(ref.residue(i).nbr_atom()).xyz()-ref.residue(i).atom("CEN").xyz()).length();

				if (r < 0.5) continue; //gly and ala should be discard
				//save this
				n_total[p.residue(i).aa()]++;

				if (dis/r > 0.6427876) continue; //sin(20)=0.342, sin(40)=0.6427876
				if (p.residue(i).aa()==chemical::aa_pro && (dis/r > 0.156)) continue; //for pro sin(18/2)
				nrecovery[p.residue(i).aa()]++;
			}

			// the logic here is: only if "fit" first, then cal recovery rate
			// ** fit and repack: output recovery rate
			// ** repack only: output err
			// ** fit only: out put err distribution(?)
			//stat
			// if (option[fit_best_rotamer]) {
			// 	utility::vector1<Size> new_rotlist;
			// 	fit_centroid_to_the_best_rot(p, new_rotlist);
			// 	for ( core::Size i=1; i<= old_rotlist.size(); ++i ) {
			// 		//only count the buried residue !!!
			// 		//std::cout << "res: " << ref.residue(i).name() << std::endl;
			// 		if (ref.pdb_info()->temperature(i, p.residue(i).nbr_atom()) < option[repack_bfactor_cutoff]) {
			// 			if ( p.pdb_info()->temperature(i, p.residue(i).nbr_atom())>0 &&
			// 				p.energies().tenA_neighbor_graph().get_node(i)
			// 				->num_edges()<option[repack_buried_cutoff]) continue;
			// 			//save the number of buried AA
			// 			n_total[p.residue(i).aa()]++;
			// 			//count right number
			// 			if (old_rotlist[i]==new_rotlist[i])
			// 				nrecovery[p.residue(i).aa()]++;
			// 		}
			// 	}
			// }
			// else {
			// 	//save err
			// 	for ( core::Size i=1; i<=p.n_residue(); ++i ) {
			// 		//cutoff Bfactor
			// 		if (ref.pdb_info()->temperature(i, p.residue(i).nbr_atom()) < option[repack_bfactor_cutoff]) {
			// 			if ( p.pdb_info()->temperature(i, p.residue(i).nbr_atom())>0 &&
			// 				p.energies().tenA_neighbor_graph().get_node(i)
			// 				->num_edges()<option[repack_buried_cutoff]) continue;
			// 			//save the number of buried AA
			// 			n_total[p.residue(i).aa()]++;
			// 			err_buried[p.residue(i).aa()]+=(p.residue(i).atom("CEN").xyz()
			// 			-ref.residue(i).atom("CEN").xyz()).length();
			// 		}
			// 	}
			// }
		}
	}

	if (option[output_cenrot_pdb]) {
		//output
		std::ostringstream outfn;
		std::string outdir(option[output_cenrot_dir]());
		std::string outpre(option[output_cenrot_prefix]());
		outfn << outdir << "/" << outpre << tag;
		//std::cout << outfn.str() <<std::endl;
		p.dump_pdb(outfn.str());
		//p.dump_scored_pdb(outfn.str(), *score_fxn, "demo");
	}
}

void rescore_pose(
	core::scoring::ScoreFunctionOP score_fxn,
	core::pose::PoseOP &native_pose,
	core::pose::Pose & p,
	std::string const &tag) {
	//silentfile
	static core::io::silent::SilentStructOP ss(
		core::io::silent::SilentStructFactory::get_instance() \
		->get_silent_struct("score"));
	static core::io::silent::SilentFileData sfd;

	//setup the ss correctly
	core::scoring::dssp::Dssp dssp( p );
	dssp.insert_ss_into_pose( p );
	(*score_fxn)(p);

	ss->fill_struct(p, tag);
	// get rmsd
	if (native_pose) {
		Real rmsd = core::scoring::CA_rmsd(p, *native_pose);
		ss->add_energy("rmsd", rmsd);
	}

	// header
	// if ( npdb==1 ) {
	// 	std::string fn(option[output_cenrot_score]());
	// 	std::ofstream scoreos(fn.c_str());
	// 	ss->print_header( scoreos );
	// }
	//write
	sfd.write_silent_struct(*ss, option[output_cenrot_score](), true );
}

void fit_centroid_to_the_best_rot( core::pose::Pose & p, utility::vector1<Size> &rotlist)
{
	//don't use this any more
	//std::cout << "since we use continuous rot, this func is no use, just return" << std::endl;
	//return;

	rotlist.clear();
	//get rotamerlib
	RotamerLibrary const & rlcap = RotamerLibrary::get_instance();
	
	for ( core::Size i=1; i<= p.total_residue(); ++i ) {
		//fit to the closest rotamer
		SingleResidueRotamerLibraryCAP residue_rotamer_library(
			rlcap.get_rsd_library(p.residue(i).type()));

		SingleResidueCenrotLibraryCAP residue_cenrot_library(
			dynamic_cast< SingleResidueCenrotLibrary const * >(residue_rotamer_library.get()));

		if (residue_rotamer_library==0){
			rotlist.push_back(1);
			continue;
		}

		Size closest_rot;
		Real closest_dis;
		CentroidRotamerSampleData const sampledata (
			residue_cenrot_library->get_closest_rotamer(
				p.residue(i), closest_rot, closest_dis));

		//get the best rot int
		id::DOF_ID id_dis(id::AtomID(p.residue(i).atom_index("CEN"), i), id::D);
		id::DOF_ID id_ang(id::AtomID(p.residue(i).atom_index("CEN"), i), id::THETA);
		id::DOF_ID id_dih(id::AtomID(p.residue(i).atom_index("CEN"), i), id::PHI);
		p.set_dof(id_dis, sampledata.distance());
		p.set_dof(id_ang, sampledata.angle());
		p.set_dof(id_dih, sampledata.dihedral());

		if (option[output_bestrot_err]) {
			TR << "FIT-ROT: res: " << p.residue(i).name3()
			<< " dis= "<< sqrt(closest_dis) << std::endl;
		}

		rotlist.push_back(closest_rot);
	}
}
