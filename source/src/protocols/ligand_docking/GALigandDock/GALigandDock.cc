// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/ga_dock/GALigandDock.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#include <protocols/ligand_docking/GALigandDock/GALigandDock.hh>
#include <protocols/ligand_docking/GALigandDock/GALigandDockCreator.hh>
#include <protocols/ligand_docking/GALigandDock/GAOptimizer.hh>
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/selection.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/rtmin.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/TopOutFunc.hh>
#include <core/scoring/constraints/util.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/import_pose/import_pose.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/relax/FastRelax.hh>

#include <basic/datacache/DataMap.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/mover_schemas.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <basic/options/option.hh> // HACK
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <ctime>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

using namespace core;

static basic::Tracer TR( "protocols.ligand_docking.GALigandDock.GALigandDock" );

GALigandDock::GALigandDock() {
	scfxn_ = core::scoring::get_score_function();

	ngen_ = 20;
	npool_ = 50;
	pmut_ = 0.2;
	rmsdthreshold_ = 1.0;
	smoothing_ = 0.0;

	grid_ = 0.25;
	padding_ = 5.0;
	hashsize_ = 8.0;
	subhash_ = 3;

	final_exact_minimize_ = "sc";
	debug_ = exact_ = false;
	sidechains_ = "none";
	sc_edge_buffer_ = 2.0; // in A
	max_rot_cumulative_prob_ = 0.99;
	rot_energy_cutoff_ = 1000;
	altcrossover_= false;
	favor_native_ = 0;

	maxiter_ = 100;
	packer_cycles_ = 20;
	init_oversample_ = 10;
	initial_pool_ = "";
	fast_relax_script_file_ = "";
	rtmin_nonideal_ = false;
	report_all_samples_="";

	runmode_ = "dock";
}

void
GALigandDock::apply( pose::Pose & pose ) {
	std::string prefix( basic::options::option[ basic::options::OptionKeys::out::prefix ]() );
	core::Size lig_resno = pose.total_residue();
	if ( ligid_.length() > 0 ) {
		lig_resno = core::pose::parse_resnum( ligid_, pose, false );
	}

	// if no cycles are specified, _AND_ final refinement is on,
	// do not do grid calculations
	bool no_grid_score = ( exact_ );

	// [[1]] setup grid scoring
	GridScorerOP gridscore(new GridScorer( scfxn_ ));
	gridscore->set_voxel_spacing( grid_ );
	gridscore->set_bbox_padding( padding_ );
	gridscore->set_hash_gridding( hashsize_ );
	gridscore->set_hash_subgridding( subhash_ );
	gridscore->set_exact( exact_ );
	gridscore->set_debug( debug_ );

	// prepare the grid but don't calculate scores yet
	TR << "Preparing grid using input pose " << std::endl;
	gridscore->prepare_grid( pose, lig_resno );

	// now figure out movable sidechains (using sidechains_ flag)
	utility::vector1< core::Size > movable_scs;
	if ( sidechains_.length() == 0 || sidechains_ == "none" || sidechains_ == "NONE" ) {
		; // do nothing
	} else if ( sidechains_.substr(0,4) == "auto" || sidechains_.substr(0,4) == "AUTO" ) {
		TR << "Detect flexible sidechains with mode: " << sidechains_ << ", "
			<< " edge_buffer = " << sc_edge_buffer_ << std::endl;

		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			if ( i == lig_resno ) continue;
			if ( pose.residue(i).aa() == core::chemical::aa_ala || pose.residue(i).aa() == core::chemical::aa_gly ) continue;
			if ( pose.residue(i).type().has_variant_type( core::chemical::DISULFIDE ) ) continue;
			if ( !pose.residue(i).is_protein() ) continue; //skip anything not amino-acid

			if ( gridscore->is_residue_in_grid(pose.residue(i), 75.0, sc_edge_buffer_ ) ) { // 75 degree "angle buffer" check
				movable_scs.push_back( i );
			}
		}
	} else if ( sidechains_.substr(0,5) == "aniso" || sidechains_.substr(0,5) == "ANISO" ) {
		TR << "Detect flexible sidechains with mode: " << sidechains_ << ", "
			<< " edge_buffer = " << sc_edge_buffer_ << std::endl;

		// model ligand as ellipsoid
		// first compute mean and covariance of ligand position
		numeric::xyzVector< core::Real > mean(0.0,0.0,0.0);
		for ( core::Size iatm=1; iatm<=pose.residue(lig_resno).natoms(); ++iatm ) {
			mean += pose.residue(lig_resno).xyz(iatm);
		}
		mean /= pose.residue(lig_resno).natoms();

		numeric::xyzMatrix< core::Real > covariance(0.0);
		for ( core::Size iatm=1; iatm<=pose.residue(lig_resno).natoms(); ++iatm ) {
			numeric::xyzVector< core::Real > diff = pose.residue(lig_resno).xyz(iatm) - mean;
			covariance += numeric::outer_product( diff , diff );
		}
		covariance /= pose.residue(lig_resno).natoms();

		// compute sorted eigenvectors
		numeric::xyzVector< core::Real > eigval, eigvalS;
		numeric::xyzMatrix< core::Real > eigvec, eigvecS;
		eigval = numeric::eigenvector_jacobi( covariance, 1e-4, eigvec );
		utility::vector1<core::Size> idx = {0,1,2};
		sort(idx.begin(), idx.end(), [&](core::Size i1, core::Size i2) {return eigval[i1] > eigval[i2];});
		TR << "idx: " << idx[1] << "," << idx[2] << "," << idx[3] << std::endl;
		for ( core::Size i=1; i<=3; ++i ) {
			eigvalS[i-1] = eigval[idx[i]];
			eigvecS.col(i-1, eigvec.col(idx[i]));
		}

		TR << "eigenvals = " << "[ " << eigvalS[0] << "," << eigvalS[1] << "," << eigvalS[2] << " ]" << std::endl;
		TR << "eigenvecs = " << "[ [ "
			<< eigvecS.xx() << "," << eigvecS.yx() << "," << eigvecS.zx() << "] ; ["
			<< eigvecS.xy() << "," << eigvecS.yy() << "," << eigvecS.zy() << "] ; ["
			<< eigvecS.xz() << "," << eigvecS.yz() << "," << eigvecS.zz() << "] ; ]"
			<< std::endl;

		// for each residue, look and see if sidechain sphere intersects ellipsoid
		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			if ( i == lig_resno ) continue;
			if ( pose.residue(i).aa() == core::chemical::aa_ala || pose.residue(i).aa() == core::chemical::aa_gly ) continue;
			if ( pose.residue(i).type().has_variant_type( core::chemical::DISULFIDE ) ) continue;
			if ( !pose.residue(i).is_protein() ) continue; //skip anything not amino-acid

			if ( gridscore->is_residue_in_grid(pose.residue(i), sc_edge_buffer_, eigvalS, eigvecS) ) {
				movable_scs.push_back( i );
			}
		}
	} else {
		// parse as residue numbers
		std::set<core::Size> resnums = core::pose::get_resnum_list( sidechains_, pose );
		for ( auto res_i=resnums.begin(); res_i!=resnums.end(); ++res_i ) {
			if ( *res_i == lig_resno ) continue;
			if ( pose.residue(*res_i).aa() == core::chemical::aa_ala
					|| pose.residue(*res_i).aa() == core::chemical::aa_gly ) continue; // don't warn
			if ( pose.residue(*res_i).type().has_variant_type( core::chemical::DISULFIDE ) ) continue; // don't warn

			if ( gridscore->is_residue_in_grid(pose.residue(*res_i), 0.0, 0.0 ) ) { // for user-defined, ignore angle/dist checks
				movable_scs.push_back( *res_i );
			} else {
				TR.Warning << "Residue " << *res_i << " is not within grid!" << std::endl;
				TR.Warning << "Increase grid padding to include it!" << std::endl;
			}
		}
	}

	// prepare the input pose
	// [a] optimize hydrogens in apo protein
	// do this only if flex sc case
	if ( movable_scs.size() > 0 ) {
		for ( core::Size i=1; i<=pose.residue(lig_resno).natoms(); ++i ) {
			pose.set_xyz( id::AtomID( i,lig_resno ), pose.residue(lig_resno).xyz(i)+250 );
		}
		core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
		task->initialize_from_command_line();
		task->or_optimize_h_mode( true );
		task->or_include_current( true );
		task->or_flip_HNQ( true );
		task->or_multi_cool_annealer( true );
		core::pack::pack_rotamers( pose, *scfxn_, task );

		// [b] idealize and minimize moving SCs
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );

		// b.1 idealize
		core::Real score0 = (*scfxn_)(pose);
		for ( core::Size i=1; i<=movable_scs.size(); ++i ) {
			utility::vector1< core::Real > chis_i = pose.residue(movable_scs[i]).chi();
			core::conformation::Residue newres( pose.residue_type(movable_scs[i]) , false);
			pose.replace_residue(movable_scs[i], newres, true);
			for ( core::Size j=1; j<=chis_i.size(); ++j ) {
				pose.set_chi( j, movable_scs[i], chis_i[j] );
			}
			mm->set_chi( movable_scs[i], true );
		}

		// b.2 minimize
		core::optimization::AtomTreeMinimizer minimizer;
		core::optimization::MinimizerOptions options( "lbfgs_armijo", 0.0001, true , false );
		options.max_iter(50);
		core::Real score1 = (*scfxn_)(pose);
		minimizer.run( pose, *mm, *scfxn_, options );
		core::Real score2 = (*scfxn_)(pose);

		TR << "Sidechain idealize score: " << score0 << "->" << score1 << "->" << score2 << std::endl;
	}

	for ( core::Size i=1; i<=pose.residue(lig_resno).natoms(); ++i ) {
		pose.set_xyz( id::AtomID( i,lig_resno ), pose.residue(lig_resno).xyz(i)-250 );
	}

	// compute the grid
	if ( !no_grid_score ) {
		TR << "Build grid for ligand and for " << movable_scs.size() << " sidechains" << std::endl;
		TR << lig_resno;
		for ( core::Size i=1; i<=movable_scs.size(); ++i ) {
			TR << "+" << movable_scs[i];
		}
		TR << std::endl;
		gridscore->calculate_grid( pose, lig_resno, movable_scs );
	} else {
		TR << "Skipping grid calculation! (exact=1)" << std::endl;
	}

	// [[2]] perturb
	core::pose::PoseOP pose_working( new core::pose::Pose( pose ) );
	//Debug gz:
	//(*pose_working).dump_pdb(prefix+"before_perturb.pdb");
	LigandConformer gene_initial( pose_working, lig_resno, movable_scs );

	// override default params; currently very close to default params...
	// use these vals for init pert & mutation
	if ( runmode_ == "refine" ) {
		gene_initial.set_rotwidth( 30.0 );
		gene_initial.set_transwidth( 1.0 );
		gene_initial.set_chiwidth( 120.0 );
		gene_initial.set_torsmutrate( 0.1 ); //mutation by 10% chance
	} else if ( runmode_ == "optligand" ) {
		if ( !exact_ ) {
			TR.Error << "runmode=optligand should be defined with exact=1 option together";
			utility_exit();
		}
		// no trans/rot
		gene_initial.set_rotwidth( 0.0 );
		gene_initial.set_transwidth( 0.0 );
		gene_initial.set_chiwidth( 120.0 );
		gene_initial.set_torsmutrate( 0.2 );
	}

	// load inputs & generate randomized starting points
	LigandConformers genes = generate_perturbed_structures( gene_initial, gridscore, protocol_[1].pool );
	//Debug gz:
	//genes[1].dump_pose(prefix+"after_perturb.pdb");

	// [[3]] main optimization cycle
	GAOptimizerOP optimizer = get_optimizer( gene_initial, gridscore );
	optimizer->run( genes );

	//for( core::Size igene = 1; igene <= genes.size(); ++igene )
	//genes[igene].dump_pose( "beforerelax"+std::to_string(igene)+".pdb" );

	// [[4]] optionally minimize a final generation with exact scores
	if ( final_exact_minimize_.substr(0,4) == "bbsc" ||  final_exact_minimize_ == "sc" || final_exact_minimize_ == "rtmin" ) {
		utility::vector1< core::pose::PoseOP > poses;
		utility::vector1<std::string> pose_names;
		std::string pose_name;
		if ( report_all_samples_ == "final_pose" ) {
			for ( core::Size i=1; i<=genes.size(); ++i ) {
				pose_name = prefix+"final."+std::to_string(i)+".pre_scmin.pdb";
				genes[i].dump_pose( pose_name );
				pose_names.push_back( pose_name );
			}
		}

		if ( final_exact_minimize_.substr(0,4) == "bbsc" ) {
			core::Size N = 0;
			if ( final_exact_minimize_.length() > 4 ) {
				N = std::atoi( final_exact_minimize_.substr(4).c_str() );
			}
			final_exact_cartmin( N, genes, poses );
		} else if ( final_exact_minimize_ == "sc" ) {
			final_exact_scmin( genes, poses );
		} else {
			final_exact_rtmin( genes, poses );
		}

		genes.clear(); // out of date now

		utility::vector1< std::pair<core::Real, core::pose::PoseOP> > poses_scored;
		utility::vector1< std::pair<std::pair<core::Real, core::pose::PoseOP>, std::string> > pose_scored_premin_index;
		for ( core::Size i=1; i<= poses.size(); ++i ) {
			core::Real score = (*scfxn_)(*poses[i]);
			TR << "apply: score after relax: " << score <<std::endl;
			poses[i]->energies().clear();
			TR << "apply: score after energies 1st clear: " << (*scfxn_)(*poses[i]) <<std::endl;
			poses_scored.push_back( std::make_pair( score, poses[i] ) );
			if ( report_all_samples_ == "final_pose" ) {
				pose_scored_premin_index.push_back( std::make_pair( std::make_pair( score, poses[i] ), pose_names[i] ));
			}
		}

		std::sort(poses_scored.begin(), poses_scored.end(),
			[&](std::pair<core::Real, core::pose::PoseOP> const &ii, std::pair<core::Real, core::pose::PoseOP> const &jj)
			{ return ii.first < jj.first; } );
		if ( report_all_samples_ == "final_pose" ) {
			std::sort(pose_scored_premin_index.begin(), pose_scored_premin_index.end(),
				[&](std::pair<std::pair<core::Real, core::pose::PoseOP>, std::string> const &ii,
				std::pair<std::pair<core::Real, core::pose::PoseOP>, std::string> const &jj)
				{ return ii.first.first < jj.first.first; } );
		}
		for ( core::Size i=1; i<=poses_scored.size(); ++i ) {
			poses_scored[i].second->energies().clear();
			if ( pose_native_ ) {
				core::Real rms = core::scoring::automorphic_rmsd(
					pose_native_->residue( lig_resno ), poses_scored[i].second->residue( lig_resno ), false );
				core::pose::setPoseExtraScore( *poses_scored[i].second, "lig_rms", rms);
			}
			if ( report_all_samples_ == "final_pose" ) {
				core::pose::setPoseExtraScore(
					*pose_scored_premin_index[i].first.second, "final_premin_pose_description", pose_scored_premin_index[i].second);
			}
			remaining_expanded_outputs_.push( poses_scored[i].second );
			TR << "apply: score after energies 2nd clear: " << (*scfxn_)(*poses_scored[i].second) <<std::endl;
		}
	} else {
		for ( core::Size i=1; i<=genes.size(); ++i ) {
			remaining_compact_outputs_.push( genes[i] );
		}
	}

	// also append input structure for reference
	if ( runmode_ == "optligand" ) {
		remaining_compact_outputs_.push(gene_initial);
	}

	pose_working = get_additional_output();
	pose = *pose_working;
}


// final optimziation cycle with sidechain flexibility
void
GALigandDock::final_exact_cartmin(
	core::Size nneigh,
	LigandConformers & genes,
	utility::vector1< core::pose::PoseOP > &poses
) {
	poses.resize( genes.size() );

	for ( core::Size i = 1; i <= genes.size(); ++i ) {
		poses[i] = core::pose::PoseOP( new core::pose::Pose );
		genes[i].to_pose( poses[i] );

		// setup movemap
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
		core::Size lig_jump = poses[i]->fold_tree().get_jump_that_builds_residue( genes[i].ligand_id() );
		mm->set_chi( genes[i].ligand_id(), true );
		mm->set_jump( lig_jump, true );
		for ( core::Size j=1; j<=genes[i].moving_scs().size(); ++j ) {
			int resid = (int)genes[i].moving_scs()[j];
			for ( int k=-((int)nneigh); k<=((int)nneigh); ++k ) {
				if ( resid+k < 1 || resid+k >= (int)poses[i]->total_residue() ) continue;
				if ( !poses[i]->residue(resid+k).is_protein() ) continue;
				mm->set_chi( resid+k, true );
				mm->set_bb( resid+k, true );
			}
		}

		// setup constraints if cst weight is on
		if ( scfxn_->get_weight( core::scoring::atom_pair_constraint ) != 0 ) {
			core::Size MINSEQSEP=2;
			core::Real MAXDIST=5.0;
			core::Size addedCsts = 0;

			for ( core::Size j=1; j<poses[i]->size(); ++j ) {
				for ( core::Size k=j+MINSEQSEP; k<poses[i]->size(); ++k ) {
					if ( !mm->get_bb(j) && !mm->get_bb(k) ) continue;

					for ( core::Size jatm=1; jatm<poses[i]->residue(j).nheavyatoms(); ++jatm ) {
						for ( core::Size katm=1; katm<poses[i]->residue(k).nheavyatoms(); ++katm ) {
							core::Real dist = poses[i]->residue(j).xyz(jatm).distance(
								poses[i]->residue(k).xyz(katm)
							);

							if ( dist > MAXDIST ) continue;
							using namespace core::scoring::func;
							FuncOP fx( new ScalarWeightedFunc( 1.0, FuncOP( new TopOutFunc( 0.25, dist, 2.0 ) ) ) );
							poses[i]->add_constraint(
								core::scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP(
								new core::scoring::constraints::AtomPairConstraint(
								core::id::AtomID(jatm,j), core::id::AtomID(katm,k), fx ) ) )
							);
							addedCsts++;
						}
					}
				}
			}
			TR << "Adding " << addedCsts << " constraints to the pose." << std::endl;
		}

		protocols::relax::FastRelax relax;
		std::vector< std::string > lines;
		if ( fast_relax_script_file_ != "" ) {
			TR << "==== Use FastRelax script: " << fast_relax_script_file_ << std::endl;
			relax=protocols::relax::FastRelax( scfxn_, fast_relax_script_file_ );
		} else {
			TR << "==== Use FastRelax hardcoded. "<< std::endl;
			lines.push_back( "switch:cartesian" );
			lines.push_back( "repeat 1" );
			lines.push_back( "ramp_repack_min 0.02 0.01 1.0 100" );
			lines.push_back( "ramp_repack_min 1.0  0.001 0.0 100" );
			lines.push_back( "accept_to_best" );
			lines.push_back( "endrepeat" );
			relax.set_script_from_lines( lines );
			relax.set_scorefxn( scfxn_ );
		}
		relax.set_movemap( mm );
		relax.apply( *poses[i] );
		TR << "final_cartmin: score after relax: " << (*scfxn_)(*poses[i]) <<std::endl;
	}
}


// final optimziation cycle with sc flexibility only
void
GALigandDock::final_exact_scmin( LigandConformers & genes, utility::vector1< core::pose::PoseOP > &poses ) {
	poses.resize( genes.size() );

	for ( core::Size i = 1; i <= genes.size(); ++i ) {
		poses[i] = core::pose::PoseOP( new core::pose::Pose );
		genes[i].to_pose( poses[i] );

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
		core::Size lig_jump = poses[i]->fold_tree().get_jump_that_builds_residue( genes[i].ligand_id() );
		mm->set_jump( lig_jump, true );
		for ( core::Size j=0; j<=genes[i].moving_scs().size(); ++j ) {
			core::Size resid = (j==0 ? genes[i].ligand_id() : genes[i].moving_scs()[j]);
			mm->set_chi( resid, true );
		}

		protocols::relax::FastRelax relax;
		std::vector< std::string > lines;
		if ( fast_relax_script_file_ != "" ) {
			TR << "==== Use FastRelax script: " << fast_relax_script_file_ << std::endl;
			relax=protocols::relax::FastRelax( scfxn_, fast_relax_script_file_ );
		} else {
			lines.push_back( "repeat 1" );
			lines.push_back( "ramp_repack_min 0.02 0.01 1.0 50" );
			lines.push_back( "ramp_repack_min 1.0  0.001 0.0 50" );
			lines.push_back( "accept_to_best" );
			lines.push_back( "endrepeat" );
			relax.set_script_from_lines( lines );
			relax.set_scorefxn( scfxn_ );
		}
		relax.set_movemap( mm );
		relax.apply( *poses[i] );
		TR << "final_scmin: score after relax: " << (*scfxn_)(*poses[i]) <<std::endl;

	}
}


//
void
GALigandDock::final_exact_rtmin(
	LigandConformers & genes,
	utility::vector1< core::pose::PoseOP > &poses
) {
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;

	poses.resize( genes.size() );
	for ( core::Size i = 1; i <= genes.size(); ++i ) {
		poses[i] = core::pose::PoseOP( new core::pose::Pose );
		genes[i].to_pose( poses[i] );

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
		core::Size lig_jump = poses[i]->fold_tree().get_jump_that_builds_residue( genes[i].ligand_id() );
		mm->set_jump( lig_jump, true );
		for ( core::Size j=0; j<=genes[i].moving_scs().size(); ++j ) {
			core::Size resid = (j==0 ? genes[i].ligand_id() : genes[i].moving_scs()[j]);
			mm->set_chi( resid, true );
		}
		core::pack::task::PackerTaskOP rtmin_task = core::pack::task::TaskFactory::create_packer_task( *poses[i] );

		TaskFactoryOP local_tf( new TaskFactory() );
		local_tf->push_back(TaskOperationCOP( new InitializeFromCommandline() ));
		local_tf->push_back(TaskOperationCOP( new RestrictToRepacking() ));
		PreventRepackingOP turn_off_packing( new PreventRepacking() );
		for ( core::Size pos = 0; pos <= (*poses[i]).size(); ++pos ) {
			core::Size resid = (pos==0 ? genes[i].ligand_id() : pos);
			if ( pos ==0 ) {
				turn_off_packing->include_residue(resid);
			} else if ( ! mm->get_chi(resid) ) {
				turn_off_packing->include_residue(resid);
			}
		}

		local_tf->push_back(turn_off_packing);
		local_tf->push_back(TaskOperationCOP( new IncludeCurrent() ));
		core::pack::task::PackerTaskOP packer( local_tf->create_task_and_apply_taskoperations( *poses[i] ) );
		//packer->task_factory(local_tf);
		(*scfxn_)(*poses[i]); // Ensure scorefunction data is appropriately initialized
		TR << "Score before rtmin: " << (*scfxn_)(*poses[i]) << std::endl;
		core::pack::RTMin RTMin;
		RTMin.set_nonideal(rtmin_nonideal_); // Default: False, change through xml
		RTMin.rtmin( *poses[i], *scfxn_, packer );
		TR << "Score after rtmin: " << (*scfxn_)(*poses[i]) << std::endl;

		core::Real tolerance(0.01);
		protocols::minimization_packing::MinMoverOP min_mover;
		if ( core::pose::symmetry::is_symmetric( *poses[i] ) )  {
			min_mover = protocols::minimization_packing::MinMoverOP( new minimization_packing::symmetry::SymMinMover( mm, scfxn_, "linmin", tolerance, true ) );
		} else {
			min_mover = protocols::minimization_packing::MinMoverOP( new protocols::minimization_packing::MinMover( mm, scfxn_, "linmin", tolerance, true ) );
		}

		min_mover->cartesian( false );
		min_mover->max_iter( 200 );

		min_mover->apply( *poses[i] );
		TR << "Score after min_mover: " << (*scfxn_)(*poses[i]) << std::endl;

	}
}


core::pose::PoseOP
GALigandDock::get_additional_output() {
	if ( remaining_expanded_outputs_.size() > 0 ) {
		core::pose::PoseOP retval = remaining_expanded_outputs_.front();
		remaining_expanded_outputs_.pop();
		(*scfxn_)(*retval);
		return retval;
	}

	if ( remaining_compact_outputs_.size() == 0 ) {
		return nullptr;
	}

	core::pose::PoseOP retval( new core::pose::Pose() );
	LigandConformer gene_i = remaining_compact_outputs_.front();
	remaining_compact_outputs_.pop();
	gene_i.to_pose( retval );
	(*scfxn_)(*retval);

	// add rms to file
	core::pose::setPoseExtraScore( *retval, "lig_rms", gene_i.rms());

	return retval;
}


// initial perturbation
LigandConformers
GALigandDock::generate_perturbed_structures(
	LigandConformer const &gene_initial,
	GridScorerOP gridscorer,
	core::Size npool
) const {
	LigandConformers genes_sel, genes_gen;

	core::Real rmscut = protocol_[1].rmsthreshold;

	// load inputs
	// fd maybe a bit wasteful to read in every time?  Should probably move this to parse_my_tag!
	if ( initial_pool_ != "" ) {
		core::Real score_initial = -9999; // ensure structures are kept though 1st filter round
		utility::vector1<std::string> input_pdbs = utility::string_split( initial_pool_, ',' );

		for ( core::Size ipdb = 1; ipdb <= input_pdbs.size(); ++ipdb ) {
			std::string tag = input_pdbs[ipdb];

			if ( (tag.length() >= 3) && (tag.substr( tag.length()-3 )=="pdb") ) {
				LigandConformer gene=gene_initial;
				gene.score( score_initial );
				core::pose::PoseOP pose = core::import_pose::pose_from_file( tag, false, core::import_pose::PDB_file );
				gene.update_conf( pose );
				genes_sel.push_back( gene );
			} else if ( tag != "INPUT" && tag != "Input" && tag != "input" ) {
				// silent file
				core::io::silent::SilentFileOptions opts; // initialized from the command line
				core::io::silent::SilentFileData sfd( opts );
				sfd.read_file( tag );
				for ( auto iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
					LigandConformer gene=gene_initial;
					gene.score( score_initial );
					core::pose::PoseOP pose(new core::pose::Pose);
					iter->fill_pose( *pose );
					gene.update_conf( pose );

					bool is_similar( false );
					for ( core::Size jj = 1; jj <= genes_sel.size() && !is_similar; ++jj ) {
						core::Real d = distance_fast( gene, genes_sel[jj] );
						is_similar = ( d < rmscut );
					}
					if ( !is_similar ) {
						genes_sel.push_back( gene );
					}

					//genes_sel.push_back( gene );
					//if( genes_inp.size() >= npool ) break;
				}
			} else {
				// "input" file
				LigandConformer gene=gene_initial;
				gene.score( score_initial );
				genes_sel.push_back( gene );
			}

			//if( genes_sel.size() >= npool ) break;
		}
	}

	int nleft = (int)npool - (int)genes_sel.size();
	int npert = debug_? nleft : (int)(nleft*init_oversample_);

	TR << "Input conformations = " << genes_sel.size() << std::endl;
	if ( npert <= 0 ) {
		return genes_sel;
	}

	TR << "Randomly generating " << npert << " structures and selecting best " << npool << std::endl;

	// if there is a smoothing schedule, use iter 1 smoothness to generate structures
	core::Real smoothing = protocol_[1].smoothing;
	gridscorer->set_smoothing( smoothing );

	// using softer repulsion for input generation
	Real const w_rep_org( gridscorer->get_w_rep() );
	gridscorer->set_w_rep( 0.02 );

	// perturb
	for ( int i=1; i<=npert; ++i ) {
		LigandConformer gene( gene_initial );
		if ( runmode_ == "refine" ) {
			gene = mutate( gene_initial ); //mutation parameters are set as small in advance
		} else {
			gene.randomize( gridscorer->get_padding() - 1.0 );  // radius of search
		}

		Real score_soft = gridscorer->score( gene, true ); // score with soft repulsive
		gene.score( score_soft );
		genes_gen.push_back( gene );
	}

	// recover original weight
	gridscorer->set_w_rep( w_rep_org );

	// fd: select lowest by energy while ensuring diversity
	// do not apply diversity criteria to inputs
	std::sort(genes_gen.begin(), genes_gen.end(),
		[&](LigandConformer const &lig_i, LigandConformer const &lig_j){ return lig_i.score() < lig_j.score(); } );
	utility::vector1< bool > selected( genes_gen.size(), false );

	for ( core::Size ii = 1; ii <= genes_gen.size(); ++ii ) {
		LigandConformer &gene_gen = genes_gen[ii];
		bool is_similar( false );
		for ( core::Size jj = 1; jj <= genes_sel.size() && !is_similar; ++jj ) {
			core::Real d = distance_fast( gene_gen, genes_sel[jj] );
			is_similar = ( d < rmscut );
		}
		if ( !is_similar ) {
			genes_sel.push_back( gene_gen );
			selected[ii] = true;
		}

		if ( genes_sel.size() == npool ) break;
	}
	for ( core::Size ii = 1; ii <= genes_gen.size() && genes_sel.size()<npool; ++ii ) {
		if ( !selected[ii] ) {
			genes_sel.push_back( genes_gen[ii] );
			selected[ii] = true;
		}
	}

	TR << "Finished generating initial pool." << std::endl;

	return genes_sel;
}

void
GALigandDock::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
) {

	scfxn_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );

	if ( tag->hasOption("ngen") ) { ngen_ = tag->getOption<int>("ngen"); }

	if ( tag->hasOption("runmode") ) { runmode_ = tag->getOption<std::string>("runmode"); }

	// allowed movement
	if ( tag->hasOption("ligand") ) { ligid_ = tag->getOption<std::string>("ligand"); }
	if ( tag->hasOption("sidechains") ) { sidechains_ = tag->getOption<std::string>("sidechains"); }
	if ( tag->hasOption("sc_edge_buffer") ) { sc_edge_buffer_ = tag->getOption<core::Real>("sc_edge_buffer"); }
	if ( tag->hasOption("rotprob") ) {  max_rot_cumulative_prob_ = tag->getOption<core::Real>("rotprob"); }
	if ( tag->hasOption("rotEcut") ) {  rot_energy_cutoff_ = tag->getOption<core::Real>("rotEcut"); }
	if ( tag->hasOption("favor_native") ) { favor_native_ = tag->getOption<core::Real>("favor_native"); }

	if ( tag->hasOption("altcrossover") ) { altcrossover_ = tag->getOption<core::Real>("altcrossover"); }

	// input params
	if ( tag->hasOption("init_oversample") ) { init_oversample_ = tag->getOption<core::Real>("init_oversample"); }
	if ( tag->hasOption("initial_pool") ) { initial_pool_ = tag->getOption<std::string>("initial_pool"); }
	if ( tag->hasOption("nativepdb") ) {
		std::string nativepdb = tag->getOption<std::string>("nativepdb");
		pose_native_ = core::import_pose::pose_from_file( nativepdb, false, core::import_pose::PDB_file );
	}

	// post-processing
	if ( tag->hasOption("final_exact_minimize") ) {
		final_exact_minimize_ = tag->getOption<std::string>("final_exact_minimize");
		if ( final_exact_minimize_ != "sc" && final_exact_minimize_.substr(0,4) != "bbsc"
				&& final_exact_minimize_ != "rtmin" && final_exact_minimize_ != "none" && final_exact_minimize_ != "aniso" ) {
			TR.Error << "The tag 'final_exact_minimize' must be one of: sc, bbscX, rtmin, aniso, none" << std::endl;
			utility_exit();
		}
	}
	//protocols
	if ( tag->hasOption("fastrelax_script") ) {
		fast_relax_script_file_ = tag->getOption<std::string>("fastrelax_script");
	}
	if ( tag->hasOption("report_all_samples") ) {
		report_all_samples_ = tag->getOption<std::string>("report_all_samples");
	}

	if ( tag->hasOption("rtmin_nonideal") ) {
		rtmin_nonideal_ = tag->getOption<bool>("rtmin_nonideal");
	}

	// grid params
	if ( tag->hasOption("exact") ) { exact_ = tag->getOption<bool>("exact"); }
	if ( tag->hasOption("debug") ) { debug_ = tag->getOption<bool>("debug"); }
	if ( tag->hasOption("grid_step") ) { grid_ = tag->getOption<core::Real>("grid_step"); }
	if ( tag->hasOption("padding") ) { padding_ = tag->getOption<core::Real>("padding"); }
	if ( tag->hasOption("hashsize") ) { hashsize_ = tag->getOption<core::Real>("hashsize"); }
	if ( tag->hasOption("subhash") ) { subhash_ = tag->getOption<core::Real>("subhash"); }

	// per-cycle defaults
	if ( tag->hasOption("npool") ) { npool_ = tag->getOption<core::Size>("npool"); }
	if ( tag->hasOption("pmut") ) { pmut_ = tag->getOption<core::Real>("pmut"); }
	if ( tag->hasOption("smoothing") ) { smoothing_ = tag->getOption<core::Real>("smoothing"); }
	if ( tag->hasOption("rmsdthreshold") ) { rmsdthreshold_ = tag->getOption<core::Real>("rmsdthreshold"); }
	if ( tag->hasOption("maxiter") ) { maxiter_ = tag->getOption<core::Size>("maxiter"); }
	if ( tag->hasOption("pack_cycles") ) { packer_cycles_ = tag->getOption<core::Size>("pack_cycles"); }
	if ( tag->hasOption("ramp_schedule") ) {
		std::string ramp_schedule_string = tag->getOption<std::string>("ramp_schedule");
		utility::vector1<std::string> ramp_schedule_stringV( utility::string_split( ramp_schedule_string , ',' ) );
		for ( std::string & scale : ramp_schedule_stringV ) {
			ramp_schedule_.push_back( atof( scale.c_str() ) );
		}
	}

	// detailed per-cycle controls
	utility::vector1< utility::tag::TagCOP > const stage_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for ( tag_it = stage_tags.begin(); tag_it != stage_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "Stage" ) {
			GADockStageParams stage_i( 1, npool_, rmsdthreshold_, pmut_, maxiter_, packer_cycles_, smoothing_, ramp_schedule_ );
			if ( (*tag_it)->hasOption("repeats") ) { stage_i.repeats = (*tag_it)->getOption<core::Size>("repeats"); }
			if ( (*tag_it)->hasOption("npool") ) { stage_i.pool = (*tag_it)->getOption<core::Size>("npool"); }
			if ( (*tag_it)->hasOption("pmut") ) { stage_i.pmut = (*tag_it)->getOption<core::Real>("pmut"); }
			if ( (*tag_it)->hasOption("smoothing") ) { stage_i.smoothing = (*tag_it)->getOption<core::Real>("smoothing"); }
			if ( (*tag_it)->hasOption("rmsdthreshold") ) { stage_i.rmsthreshold = (*tag_it)->getOption<core::Real>("rmsdthreshold"); }
			if ( (*tag_it)->hasOption("maxiter") ) { stage_i.maxiter = (*tag_it)->getOption<core::Size>("maxiter"); }
			if ( (*tag_it)->hasOption("pack_cycles") ) { stage_i.packcycles = (*tag_it)->getOption<core::Size>("pack_cycles"); }
			if ( (*tag_it)->hasOption("ramp_schedule") ) {
				std::string ramp_schedule_string = (*tag_it)->getOption<std::string>("ramp_schedule");
				utility::vector1<std::string> ramp_schedule_stringV( utility::string_split( ramp_schedule_string , ',' ) );
				stage_i.ramp_schedule.clear();
				for ( std::string & scale : ramp_schedule_stringV ) {
					stage_i.ramp_schedule.push_back( atof( scale.c_str() ) );
				}
			}

			protocol_.push_back( stage_i );
		}
	}

	// sanity checks
	if ( tag->hasOption("ngen") && protocol_.size() > 0 ) {
		TR.Error << "ngen specified but detailed schedule also given.  Aborting!" << std::endl;
		utility_exit();
	}

	// default protocol
	if ( protocol_.size() == 0 ) {
		TR << "Using default protocol." << std::endl;
		GADockStageParams stage_i( ngen_, npool_, rmsdthreshold_, pmut_, maxiter_, packer_cycles_, smoothing_, ramp_schedule_ );
		protocol_.push_back( stage_i );
	} else {
		TR << "Using custom " << protocol_.size() << "-stage protocol." << std::endl;
	}
}


GAOptimizerOP
GALigandDock::get_optimizer(
	LigandConformer const &gene_initial,
	GridScorerOP gridscorer
) const {
	GAOptimizerOP optimizer(new GAOptimizer(gridscorer));
	if ( pose_native_ ) {
		LigandConformer gene_native( pose_native_, gene_initial.ligand_id(), gene_initial.moving_scs() );
		optimizer->set_native( gene_native );
	}

	optimizer->set_protocol( protocol_ );
	optimizer->set_max_rot_cumulative_prob( max_rot_cumulative_prob_ );
	optimizer->set_rot_energy_cutoff( rot_energy_cutoff_ );  // at some point make this a parameter?
	optimizer->set_favor_native( favor_native_ );
	optimizer->set_altcrossover( altcrossover_ );
	return optimizer;
}

std::string GALigandDock::get_name() const {
	return mover_name();
}

std::string GALigandDock::mover_name() {
	return "GALigandDock";
}

// xml stuff
std::string gadock_subelement_ct_name( std::string const & name ) {
	return "GALigandDock_subelement_" + name + "Type";
}

void GALigandDock::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "scorefxn", xs_string, "weights file");

	attlist + XMLSchemaAttribute( "runmode", xs_string, "run mode [dock/refine/optligand]");

	attlist + XMLSchemaAttribute( "altcrossover", xsct_rosetta_bool, "Use alternate xover.");
	attlist + XMLSchemaAttribute( "init_oversample", xsct_real, "scale factor to ntrial of initial pool generation");
	attlist + XMLSchemaAttribute( "rotprob", xsct_real, "max cumulative rotamer probability");
	attlist + XMLSchemaAttribute( "rotEcut", xsct_real, "rotamer 1b energy");
	attlist + XMLSchemaAttribute( "ligand", xs_string, "ligand residue id (if not specified default to last residue)");
	attlist + XMLSchemaAttribute( "nativepdb", xs_string, "name of native pdb");
	attlist + XMLSchemaAttribute( "favor_native", xsct_real, "give a bonus score to the input rotamer");
	attlist + XMLSchemaAttribute( "grid_step", xsct_real, "Grid step (A) for grid-based scoring");
	attlist + XMLSchemaAttribute( "padding", xsct_real, "Padding (A) step for grid-based scoring");
	attlist + XMLSchemaAttribute( "hashsize", xsct_real, "Width of hash bins (A)");
	attlist + XMLSchemaAttribute( "subhash", xsct_non_negative_integer, "When scanning gridspace, subhash to this level");
	attlist + XMLSchemaAttribute( "final_exact_minimize", xs_string, "Minimize the Genes by exact score after GA.");
	attlist + XMLSchemaAttribute( "rtmin_nonideal", xsct_rosetta_bool, "Set nonideal to true in rtmin. Default: false");
	attlist + XMLSchemaAttribute( "fastrelax_script", xs_string, "FastRelax script file for exact minimize.");
	attlist + XMLSchemaAttribute( "exact", xsct_rosetta_bool, "Use exact scoring.");
	attlist + XMLSchemaAttribute( "debug", xsct_rosetta_bool, "Debug grid scoring: report both exact and grid scores.");
	attlist + XMLSchemaAttribute( "initial_pool", xs_string, "Include these structures in the initial pool.");
	attlist + XMLSchemaAttribute( "report_all_samples", xs_string, "Dump intermediate structures/scores sampled.");
	attlist + XMLSchemaAttribute( "sidechains", xs_string, "Sidechains to move: none, auto, or residue IDs.");
	attlist + XMLSchemaAttribute( "sc_edge_buffer", xsct_real, "Scaling factor of maxdistance when deciding to include sc as movable");

	// per-cycle parameters (defaults)
	attlist + XMLSchemaAttribute( "ngen", xs_integer, "number of generations");
	attlist + XMLSchemaAttribute( "npool", xsct_non_negative_integer, "(default) pool size");
	attlist + XMLSchemaAttribute( "pmut", xsct_real, "(default) probability of mutation");
	attlist + XMLSchemaAttribute( "smoothing", xsct_real, "(default) grid smoothing");
	attlist + XMLSchemaAttribute( "rmsdthreshold", xsct_real, "(default) RMSD threshold between pool structures");
	attlist + XMLSchemaAttribute( "ramp_schedule", xs_string, "(default) During minimization, ramp fa_rep according to this schedule");
	attlist + XMLSchemaAttribute( "maxiter", xsct_non_negative_integer, "(default) maxiter for minimizer");
	attlist + XMLSchemaAttribute( "pack_cycles", xsct_non_negative_integer, "(default) pack for (N x #res) cycles");

	// attributes for "Stage" subelement
	AttributeList stage_subelement_attributes;
	stage_subelement_attributes
		+ XMLSchemaAttribute( "repeats", xs_integer, "number of generations in this stage")
		+ XMLSchemaAttribute( "npool", xsct_non_negative_integer, "pool size in this stage" )
		+ XMLSchemaAttribute( "smoothing", xsct_real, "Grid smoothing in this stage" )
		+ XMLSchemaAttribute( "pmut", xsct_real, "Sampling frequency weight for this template" )
		+ XMLSchemaAttribute( "rmsdthreshold", xsct_real, "symmdef file associated with this template (only if using symmetry)" )
		+ XMLSchemaAttribute( "ramp_schedule", xs_string, "comma-seprated list of chains to randomize - not documented" )
		+ XMLSchemaAttribute( "maxiter", xsct_non_negative_integer, "maxiter for minimizer" )
		+ XMLSchemaAttribute( "pack_cycles", xsct_non_negative_integer, "pack for (N x #res) cycles");

	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & gadock_subelement_ct_name );
	subelements.add_simple_subelement( "Stage", stage_subelement_attributes, "Per-stage parameters");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(), "This mover runs ligand docking using a GA (with gridded scoring) to optimize ligand-protein interaction energies.",
		attlist, subelements );
}

protocols::moves::MoverOP
GALigandDock::clone() const {
	return( protocols::moves::MoverOP( new GALigandDock( *this ) ) );
}

protocols::moves::MoverOP
GALigandDock::fresh_instance() const {
	return protocols::moves::MoverOP( new GALigandDock );
}

std::string GALigandDockCreator::keyname() const {
	return GALigandDock::mover_name();
}

protocols::moves::MoverOP
GALigandDockCreator::create_mover() const {
	return protocols::moves::MoverOP( new GALigandDock );
}

void GALigandDockCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	GALigandDock::provide_xml_schema( xsd );
}

} // ga_dock
} // ligand_docking
} // protocols


