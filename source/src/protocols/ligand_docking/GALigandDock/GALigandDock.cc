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
#include <protocols/ligand_docking/GALigandDock/util.hh>
#include <protocols/ligand_docking/GALigandDock/EntropyEstimator.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/util.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>
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
#include <protocols/simple_moves/WaterBoxMover.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/TopOutFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/import_pose/import_pose.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>

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
#include <numeric/random/random.hh>
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
	scfxn_relax_ = core::scoring::get_score_function();

	ngen_ = 20;
	npool_ = 50;
	nrelax_ = 20;
	nreport_ = 20;
	pmut_ = 0.2;
	rmsdthreshold_ = 1.0;
	smoothing_ = 0.0;
	sample_ring_conformers_ = true; // still only works if params contains proper ring definition

	// grid
	grid_ = 0.25;
	padding_ = 4.5;
	hashsize_ = 8.0;
	subhash_ = 3;
	fa_rep_grid_ = 0.2;
	grid_bound_penalty_ = 100.0;

	debug_ = exact_ = false;
	sidechains_ = "none";
	sc_edge_buffer_ = 2.0; // in A
	altcrossover_= false;
	optimize_input_H_ = true;

	// final relaxation
	final_exact_minimize_ = "sc";
	min_neighbor_ = false;
	cartmin_lig_ = true;
	full_repack_before_finalmin_ = false;
	final_solvate_ = false;
	fast_relax_script_file_ = "";
	fast_relax_lines_ = std::vector<std::string>();
	redefine_flexscs_at_relax_ = false;

	// estimate and report dG
	estimate_dG_ = false;

	// packing behavior
	max_rot_cumulative_prob_ = 0.9;
	rot_energy_cutoff_ = 100;
	maxiter_ = 50;
	packer_cycles_ = 100;
	favor_native_ = 2;

	// input generation
	premin_ligand_ = false;
	random_oversample_ = 10;
	reference_oversample_ = 2;
	reference_frac_ = 0.5;
	reference_frac_auto_ = true;
	use_pharmacophore_ = true; // turn on by default!

	multiple_ligands_ = utility::vector1< std::string >();
	initial_pool_ = "";
	reference_pool_ = "";

	move_water_ = true;

	runmode_ = "dockflex"; // high-resolution pharmacophore docking
}

void
GALigandDock::apply( pose::Pose & pose )
{
	std::string prefix( basic::options::option[ basic::options::OptionKeys::out::prefix ]() );
	core::Size lig_resno = pose.total_residue();
	if ( ligid_.length() > 0 ) {
		lig_resno = core::pose::parse_resnum( ligid_, pose, false );
	}

	// if no cycles are specified, _AND_ final refinement is on,
	// do not do grid calculations
	bool no_grid_score = ( exact_ );

	// [[1]] setup grid scoring
	GridScorerOP gridscore( new GridScorer( scfxn_ ));
	gridscore->set_voxel_spacing( grid_ );
	gridscore->set_bbox_padding( padding_ );
	gridscore->set_hash_gridding( hashsize_ );
	gridscore->set_hash_subgridding( subhash_ );
	gridscore->set_exact( exact_ );
	gridscore->set_debug( debug_ );
	gridscore->set_out_of_bound_e( grid_bound_penalty_ ); // should perhaps make gridscore parse tags...

	// prepare the grid but don't calculate scores yet
	TR << "Preparing grid using input pose " << std::endl;
	gridscore->prepare_grid( pose, lig_resno );

	// now figure out movable sidechains (using sidechains_ flag)
	utility::vector1< core::Size > movable_scs = get_movable_scs( pose, gridscore, lig_resno );

	// pass all the ligand residues for grid construction
	utility::vector1< core::conformation::Residue > rsds_to_build_grids;
	rsds_to_build_grids.push_back( pose.residue(lig_resno) );
	for ( core::Size ires = 1; ires <= movable_scs.size(); ++ires ) {
		rsds_to_build_grids.push_back( pose.residue(movable_scs[ires]) );
	}
	if ( multiple_ligands_.size() > 0 ) {
		for ( core::Size ilig = 1; ilig <= multiple_ligands_.size(); ++ilig ) {
			core::conformation::ResidueOP ligand = core::conformation::get_residue_from_name( multiple_ligands_[ilig] );
			rsds_to_build_grids.push_back( *ligand );
		}
	}
	gridscore->get_grid_atomtypes( rsds_to_build_grids );

	// prepare the input pose
	idealize_and_repack_pose( pose, movable_scs, lig_resno );

	// compute the grid
	if ( !no_grid_score ) {
		TR << "Build grid for ligand and for " << movable_scs.size() << " sidechains" << std::endl;
		TR << "Residue nums: " << lig_resno;
		for ( core::Size i=1; i<=movable_scs.size(); ++i ) {
			TR << "+" << movable_scs[i];
		}
		TR << std::endl;
		gridscore->calculate_grid( pose, lig_resno, movable_scs );

	} else {
		TR << "Skipping grid calculation! (exact=1)" << std::endl;
		TR << "Mobile sidechains: ";
		TR << lig_resno;
		for ( core::Size i=1; i<=movable_scs.size(); ++i ) {
			TR << "+" << movable_scs[i];
		}
		TR << std::endl;
	}

	// Prepare ligand aligner
	LigandAligner aligner;
	if ( reference_pool_ != "" ) {
		scfxn_->score( pose ); // make sure scored to get neighbor graph
		TR << "Setting up Ligand aligner." << std::endl;
		aligner = setup_ligand_aligner( pose, lig_resno, movable_scs );
	}

	if ( multiple_ligands_.size() > 0 ) {
		/*
		core::pose::PoseOP pose_apo( new core::pose::Pose( pose ) );
		pose_apo->energies().clear();
		pose_apo->data().clear();
		core::conformation::Residue const &ligand_ref = pose.residue(lig_resno);
		pose_apo->delete_residue_slow( lig_resno );
		lig_resno = pose.size(); // reset it to the last one
		*/

		for ( core::Size ilig = 1; ilig <= multiple_ligands_.size(); ++ilig ) {
			TR << "===============================================================================" << std::endl;
			TR << " Starting " <<  multiple_ligands_[ilig]
				<< " (" << ilig << "/" << multiple_ligands_.size() << ")" << std::endl;
			TR << "===============================================================================" << std::endl;

			core::pose::PoseOP pose_working =
				make_starting_pose_for_virtual_screening( pose, lig_resno, multiple_ligands_[ilig] );
			//pose_working->dump_pdb("start."+multiple_ligands_[ilig]+".pdb");

			if ( premin_ligand_ ) {
				TR << "Preminimize ligand." << std::endl;
				premin_ligand( pose, lig_resno );
			}

			LigandConformer gene_initial( pose_working, lig_resno, movable_scs );
			gene_initial.set_sample_ring_conformers( sample_ring_conformers_ );

			OutputStructureStore temporary_outputs;
			// take lowest score pose from each ligand
			pose = run_docking( gene_initial, gridscore, aligner, temporary_outputs );

			// store to remaining outputs
			core::Real score, rms, ligscore, recscore;
			std::string ligandname;
			score = (*scfxn_relax_)(pose);
			core::pose::getPoseExtraScore( pose, "ligscore", ligscore );
			core::pose::getPoseExtraScore( pose, "recscore", recscore );
			core::pose::getPoseExtraScore( pose, "lig_rms", rms );
			core::pose::getPoseExtraScore( pose, "ligandname", ligandname );
			//ignore ranking_prerelax
			remaining_outputs_.push( pose, score, rms, ligscore, recscore, 0, ligandname );
		}

	} else {
		core::pose::PoseOP pose_working( new core::pose::Pose( pose ) );
		pose_working->energies().clear();
		pose_working->data().clear();
		if ( premin_ligand_ ) {
			TR << "Preminimize ligand." << std::endl;
			premin_ligand( pose, lig_resno );
		}

		LigandConformer gene_initial( pose_working, lig_resno, movable_scs );
		gene_initial.set_sample_ring_conformers( sample_ring_conformers_ );

		pose = run_docking( gene_initial, gridscore, aligner, remaining_outputs_ );
	}
}

core::pose::Pose
GALigandDock::run_docking( LigandConformer const &gene_initial,
	GridScorerOP gridscore,
	LigandAligner &aligner,
	OutputStructureStore &outputs )
{
	core::Size const lig_resno( gene_initial.ligand_id() );

	// load inputs & generate randomized starting points
	LigandConformers genes = generate_perturbed_structures( gene_initial, gridscore, protocol_[1].pool,
		aligner );

	// [[3]] main optimization cycle

	GAOptimizerOP optimizer = get_optimizer( gene_initial, gridscore );
	optimizer->run( genes );

	// trim genes for final min & report
	if ( nrelax_ <= genes.size() ) {
		genes.resize( nrelax_ );
	} else {
		TR << "Warning: Nrelax(" << nrelax_ << ") >= Npool(" << genes.size() << "); ignore." << std::endl;
	}

	// [[4]] optionally minimize a final generation with exact scores
	// FD this is starting to get ugly...
	bool finalbbscmin = (final_exact_minimize_.substr(0,4) == "bbsc");
	bool finalscmin = (final_exact_minimize_.substr(0,2) == "sc");
	//bool dualrelax = (finalbbscmin && final_exact_minimize_.length() > 8 && final_exact_minimize_.substr(5,9) == "dual");

	utility::vector1< core::Real > dEs;
	for ( core::Size i=1; i<=genes.size(); ++i ) {
		core::pose::PoseOP pose_tmp( new core::pose::Pose );
		genes[i].to_pose( pose_tmp );
		//pose_tmp->dump_pdb("premin."+std::to_string(i)+".pdb");

		utility::vector1< core::Size > const &movable_scs = genes[i].moving_scs();
		// idealize again... why do we need this again here?
		for ( core::Size i=1; i<=movable_scs.size(); ++i ) {
			utility::vector1< core::Real > chis_i = pose_tmp->residue(movable_scs[i]).chi();
			core::conformation::Residue newres( pose_tmp->residue_type(movable_scs[i]) , false);
			pose_tmp->replace_residue(movable_scs[i], newres, true);
			for ( core::Size j=1; j<=chis_i.size(); ++j ) {
				pose_tmp->set_chi( j, movable_scs[i], chis_i[j] );
			}
		}

		if ( final_solvate_ ) {
			final_solvate( genes[i], *pose_tmp );
		}

		if ( cartmin_lig_ ) {
			final_cartligmin(genes[i], *pose_tmp );
		}
		//pose_tmp->dump_pdb("cartmin1."+std::to_string(i)+".pdb");

		if ( finalbbscmin ) {
			core::Size N = 0;
			if ( final_exact_minimize_.length() > 4 ) {
				N = std::atoi( final_exact_minimize_.substr(4).c_str() );
			}
			final_exact_cartmin( N, genes[i], *pose_tmp );

			// post relax
			if ( cartmin_lig_ ) {
				final_cartligmin(genes[i], *pose_tmp);
			}
		} else if ( finalscmin ) {
			final_exact_scmin( genes[i], *pose_tmp );
			// post relax
			if ( cartmin_lig_ ) {
				final_cartligmin(genes[i], *pose_tmp );
			}
		}

		pose_tmp->energies().clear();
		pose_tmp->data().clear();
		core::Real score = (*scfxn_relax_)(*pose_tmp);

		//TR << "FINAL score for " << i << "-th"; scfxn_relax_->show(TR,*pose_tmp);

		core::Real rms = 0.0;
		if ( pose_native_ ) {
			//fd  if the input ligand is the last residue, use the last residue of the _native_ as the ligand
			//fd  otherwise, match residue IDs
			core::Size native_lig = lig_resno;
			if ( lig_resno == pose_tmp->total_residue() ) native_lig = pose_native_->total_residue();
			rms = core::scoring::automorphic_rmsd(
				pose_native_->residue( native_lig ), pose_tmp->residue( lig_resno ), false );
		}

		// report ligand-only energy
		core::Real ligscore = calculate_free_ligand_score( pose_tmp->residue( lig_resno ) );
		core::Real recscore = calculate_free_receptor_score( *pose_tmp, lig_resno, movable_scs, true );
		core::Real dE( score - recscore - ligscore );
		std::string ligandname( pose_tmp->residue(lig_resno).name() );
		outputs.push( *pose_tmp, score, rms, ligscore, recscore, i, ligandname );

		dEs.push_back( dE );
	}

	// lowest energy; use output class function instead of overrided one
	//core::pose::PoseOP pose = get_additional_output();
	core::pose::PoseOP pose = outputs.pop();

	if ( estimate_dG_ ) {
		EntropyEstimator entropy_estimator( scfxn_relax_, *pose, lig_resno );
		//entropy_estimator.simple( runmode_=="VSX" ); // TODO
		if ( runmode_ == "VSX" ) entropy_estimator.set_niter( 1000 );
		core::Real TdS = entropy_estimator.apply( *pose ); //comes out in energy unit; sign is opposite

		std::sort( dEs.begin(), dEs.end() ); // default comparator
		core::Real mindE = dEs[1];

		core::Real dG = mindE + TdS;

		TR << "Estimated Binding Free Energy (arbitrary energy unit, just for relative ranking)" << std::endl;
		TR << "dH: " << std::setw(6) << mindE << std::endl;
		TR << "-T*dS: " << std::setw(6) << TdS << std::endl;
		TR << "dG (dH-T*dS): " << dG << std::endl;
	}

	return *pose; // return lowest energy one
}

utility::vector1< core::Size >
GALigandDock::get_movable_scs( core::pose::Pose const &pose,
	GridScorerCOP gridscore,
	core::Size const lig_resno ) const
{
	utility::vector1< core::Size > movable_scs;

	if ( sidechains_.length() == 0 || sidechains_ == "none" || sidechains_ == "NONE" ) {
		; // do nothing
	} else if ( sidechains_.substr(0,4) == "auto" || sidechains_.substr(0,4) == "AUTO" ) {
		TR << "Detect flexible sidechains with mode: " << sidechains_ << ", "
			<< " edge_buffer=" << sc_edge_buffer_ << std::endl;

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
			eigvecS.col(i, eigvec.col(idx[i]+1));  // .col is 1-indexed (for some reason)
		}

		// special case for planar ligands
		eigvalS[2] /= eigvalS[0];
		eigvalS[1] /= eigvalS[0];
		eigvalS[0] = 1.0;

		eigvalS[2] = std::max( eigvalS[2], 0.05 );
		eigvalS[1] = std::max( eigvalS[1], 0.05 );

		TR << "eigenvals = " << "[ " << eigvalS[0] << "," << eigvalS[1] << "," << eigvalS[2] << " ]" << std::endl;
		//TR << "eigenvecs = " << "[ [ "
		// << eigvecS.xx() << "," << eigvecS.yx() << "," << eigvecS.zx() << "] ; ["
		// << eigvecS.xy() << "," << eigvecS.yy() << "," << eigvecS.zy() << "] ; ["
		// << eigvecS.xz() << "," << eigvecS.yz() << "," << eigvecS.zz() << "] ; ]"
		// << std::endl;

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

	return movable_scs;
}

void
GALigandDock::idealize_and_repack_pose( core::pose::Pose &pose,
	utility::vector1< core::Size > const &movable_scs,
	core::Size const lig_resno ) const
{
	// do this only if flex sc case
	if ( movable_scs.size() == 0 ) return;

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
	core::Real score0 = (*scfxn_)(pose);

	for ( core::Size i=1; i<=pose.residue(lig_resno).natoms(); ++i ) {
		pose.set_xyz( id::AtomID( i,lig_resno ), pose.residue(lig_resno).xyz(i)+250 );
	}

	// a. idealize
	for ( core::Size i=1; i<=movable_scs.size(); ++i ) {
		utility::vector1< core::Real > chis_i = pose.residue(movable_scs[i]).chi();
		core::conformation::Residue newres( pose.residue_type(movable_scs[i]) , false);
		pose.replace_residue(movable_scs[i], newres, true);
		for ( core::Size j=1; j<=chis_i.size(); ++j ) {
			pose.set_chi( j, movable_scs[i], chis_i[j] );
			mm->set_chi( movable_scs[i], true );
		}
	}

	// opt hydrogen as apo + minimize; only for "docking" case
	// b.2 optimize hydrogen
	if ( optimize_input_H_ ) {
		TR << "Re-optimizing hydrogens in whole structure." << std::endl;
		core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
		task->initialize_from_command_line();
		task->or_optimize_h_mode( true );
		task->or_include_current( true );
		task->or_flip_HNQ( true );
		task->or_multi_cool_annealer( true );
		core::pack::pack_rotamers( pose, *scfxn_, task );
	}

	// b.3 minimize
	core::optimization::AtomTreeMinimizer minimizer;
	core::optimization::MinimizerOptions options( "lbfgs_armijo", 0.0001, true , false );
	options.max_iter(50);
	core::Real score1 = (*scfxn_)(pose);
	minimizer.run( pose, *mm, *scfxn_, options );
	core::Real score2 = (*scfxn_)(pose);
	TR << "Sidechain idealize score: " << score0 << "->" << score1 << "->" << score2 << std::endl;

	for ( core::Size i=1; i<=pose.residue(lig_resno).natoms(); ++i ) {
		pose.set_xyz( id::AtomID( i,lig_resno ), pose.residue(lig_resno).xyz(i)-250 );
	}
}

// generate a pose for virtual screening with the ligand corresponding to ligand_name
core::pose::PoseOP
GALigandDock::make_starting_pose_for_virtual_screening( core::pose::Pose const &pose,
	core::Size const &lig_resno,
	std::string const ligand_name
) const
{
	// will it die here if ligand name not exist in residue type set??
	core::conformation::ResidueOP ligand = core::conformation::get_residue_from_name( ligand_name );
	/*
	core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()
	->residue_type_set( "fa_standard" ) );
	core::chemical::ResidueTypeCOP rsd_type = residue_set->name_map( name ).get_self_ptr();
	ligand = conformation::ResidueFactory::create_residue( *rsd_type );
	*/

	core::Size jumpid = pose.fold_tree().get_jump_that_builds_residue( lig_resno );
	core::kinematics::Jump ligjump = pose.jump( jumpid );
	numeric::xyzVector< core::Real > T = ligjump.get_translation();

	core::pose::PoseOP pose_working( new core::pose::Pose( pose ) );
	pose_working->replace_residue( lig_resno, *ligand, false );

	/*
	// refine placement by delta com?
	numeric::xyzVector< core::Real > com_ref( 0.0 ), com( 0.0 );
	for( core::Size iatm = 1; iatm <= ligand_ref.nheavyatoms(); ++iatm )
	com_ref += ligand_ref.xyz(iatm);
	com_ref /= ligand_ref.nheavyatoms();

	for( core::Size iatm = 1; iatm <= ligand->nheavyatoms(); ++iatm )
	com += ligand->xyz(iatm);
	com /= ligand->nheavyatoms();
	*/

	//numeric::xyzVector< core::Real > T2 = ligjump.get_translation();
	ligjump.set_translation( T ); // no change in orientation; is this line necessary?
	pose_working->set_jump( jumpid, ligjump );

	return pose_working;
}

Real
GALigandDock::calculate_free_receptor_score( core::pose::Pose pose, // call by value
	core::Size const lig_resno,
	utility::vector1< core::Size > const& moving_scs,
	bool simple
) const
{
	if ( pose.size() == 1 ) return 0.0; // ligand-only

	// necessary for memory issue?
	pose.energies().clear();
	pose.data().clear();

	// delete ligand
	pose.delete_residue_slow( lig_resno );

	if ( simple ) {
		return (*scfxn_relax_)(pose);
	}

	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
	task->initialize_from_command_line();
	task->or_optimize_h_mode( false ); // to reduce noise
	task->or_flip_HNQ( false ); // to reduce noise

	utility::vector1< core::Size > moving_his_scs;
	for ( core::Size ires = 1; ires <= moving_scs.size(); ++ires ) {
		if ( pose.residue(moving_scs[ires]).aa() == core::chemical::aa_his ) moving_his_scs.push_back( moving_scs[ires] );
	}
	task->or_fix_his_tautomer( moving_his_scs, true ); // to reduce noise
	//task->or_multi_cool_annealer( true );

	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->modify_task( pose, task );

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
	for ( core::Size j=1; j<=moving_scs.size(); ++j ) {
		core::Size resid = moving_scs[j];
		mm->set_chi( resid, true );
	}

	// let's use hard-coded scopt schedule for now...
	protocols::relax::FastRelax relax;
	relax.set_task_factory( tf );

	// sometimes just stuck here wasting huge memory... why is it?

	std::vector< std::string > lines;
	if ( fast_relax_script_file_ != "" ) {
		TR << "==== Use FastRelax script: " << fast_relax_script_file_ << std::endl;
		relax=protocols::relax::FastRelax( scfxn_relax_, fast_relax_script_file_ );
	} else if ( fast_relax_lines_.size() > 0 ) {
		lines = fast_relax_lines_;
	} else {
		TR << "==== Use FastRelax hardcoded. "<< std::endl;
		lines.push_back( "switch:torsion" );
		lines.push_back( "repeat 3" );
		lines.push_back( "ramp_repack_min 0.02 0.01 1.0 50" );
		lines.push_back( "ramp_repack_min 1.0  0.00001 0.0 50" );
		lines.push_back( "accept_to_best" );
		lines.push_back( "endrepeat" );
		relax.set_script_from_lines( lines );
		relax.set_scorefxn( scfxn_relax_ );
	}
	relax.set_movemap( mm );
	relax.apply( pose );

	return (*scfxn_relax_)(pose);
}

Real
GALigandDock::calculate_free_ligand_score( core::conformation::Residue const ligand ) const
{
	// make a ligand-only pose; root ligand on virtual if no residues to anchor jump
	core::pose::PoseOP pose( new core::pose::Pose );
	pose->append_residue_by_jump( ligand, 0 );
	core::pose::addVirtualResAsRoot(*pose);

	// optimize slightly...
	{
		core::scoring::ScoreFunctionOP scfxn_ligmin( scfxn_relax_ );
		scfxn_ligmin->set_weight( core::scoring::coordinate_constraint, 1.0 );

		core::pose::Pose pose_premin( *pose );
		for ( Size iatm = 1; iatm <= ligand.natoms(); ++iatm ) {
			core::id::AtomID atomid( iatm, 1 );
			core::Vector const &xyz = pose->xyz( atomid );
			core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
			pose->add_constraint( core::scoring::constraints::ConstraintCOP
				( core::scoring::constraints::ConstraintOP
				( new core::scoring::constraints::CoordinateConstraint( atomid, atomid, xyz, fx ) )));
		}

		//Real score0 = scfxn_ligmin->score( *pose );

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( true ); mm->set_chi( true ); mm->set_jump( false );

		protocols::minimization_packing::MinMoverOP min_mover =
			protocols::minimization_packing::MinMoverOP
			( new protocols::minimization_packing::MinMover( mm, scfxn_ligmin, "linmin", 0.01, true ) );
		min_mover->max_iter( 30 );
		min_mover->apply( *pose );

		// make sure structure hasn't changed much by minimization
		//core::Real rmsd_by_min = core::scoring::automorphic_rmsd( pose->residue(1), pose_premin.residue(1), false );
		//Real score = scfxn_ligmin->score( *pose );
		//std::cout << "RMSD: " << rmsd_by_min << ", score drop" << score0 << " -> " << score
		//     << std::endl;
	}

	Real ligandscore = scfxn_relax_->score( *pose );

	return ligandscore;
}

void
GALigandDock::premin_ligand( core::pose::Pose &pose, core::Size const lig_resno ) const
{
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
	mm->set_chi( lig_resno, true );

	core::Real w_cart = (*scfxn_relax_)[ core::scoring::cart_bonded ];
	core::Real w_proclose = (*scfxn_relax_)[ core::scoring::pro_close ];
	core::Real w_atr = (*scfxn_relax_)[ core::scoring::fa_atr ];
	core::Real w_rep = (*scfxn_relax_)[ core::scoring::fa_rep ];
	if ( w_cart < 1.0e-6 ) {
		scfxn_relax_->set_weight( core::scoring::cart_bonded, 0.5 );
		scfxn_relax_->set_weight( core::scoring::pro_close, 0.0 );
	}
	// turn non-bonded terms off
	scfxn_relax_->set_weight( core::scoring::fa_atr, 0.0 );
	scfxn_relax_->set_weight( core::scoring::fa_rep, 0.0 );

	core::optimization::CartesianMinimizer minimizer;
	core::optimization::MinimizerOptions options( "lbfgs_armijo", 0.0001, true , false );
	options.max_iter(50);
	minimizer.run( pose, *mm, *scfxn_relax_, options );

	// reset
	scfxn_relax_->set_weight( core::scoring::cart_bonded, w_cart );
	scfxn_relax_->set_weight( core::scoring::pro_close, w_proclose );
	scfxn_relax_->set_weight( core::scoring::fa_atr, w_atr );
	scfxn_relax_->set_weight( core::scoring::fa_rep, w_rep );
}

void
GALigandDock::apply_coord_cst_to_sctip( core::pose::PoseOP pose,
	utility::vector1< core::Size > const& moving_scs
) const
{

	TR << "Applying coordinate cst to sidechain tips:" << std::endl;
	std::string reportline( "CSTATMS " );
	utility::vector1< core::id::AtomID > cstatoms;
	for ( core::Size isc = 1; isc <= moving_scs.size(); ++isc ) {
		std::string aname;
		core::Size ires( moving_scs[isc] );
		core::chemical::AA aa( pose->aa( ires ) );
		switch ( aa ){
		case core::chemical::aa_cys : { aname = "SG"; break; }
		case core::chemical::aa_asp : { aname = "CG"; break; }
		case core::chemical::aa_glu : { aname = "CD"; break; }
		case core::chemical::aa_phe : { aname = "CZ"; break; }
		case core::chemical::aa_his : { aname = "NE2"; break; }
		case core::chemical::aa_ile : { aname = "CG2"; break; }
		case core::chemical::aa_lys : { aname = "NZ"; break; }
		case core::chemical::aa_leu : { aname = "CD1"; break; }
		case core::chemical::aa_met : { aname = "CE"; break; }
		case core::chemical::aa_asn : { aname = "OD1"; break; }
		case core::chemical::aa_gln : { aname = "OE1"; break; }
		case core::chemical::aa_arg : { aname = "CZ"; break; }
		case core::chemical::aa_ser : { aname = "OG"; break; }
		case core::chemical::aa_thr : { aname = "OG1"; break; }
		case core::chemical::aa_val : { aname = "CG1"; break; }
		case core::chemical::aa_trp : { aname = "CZ2"; break; }
		case core::chemical::aa_tyr : { aname = "OH"; break; }
		default : continue;
		}

		reportline += " "+std::to_string(ires)+"."+aname;
		core::Size iatm( pose->residue(ires).atom_index(aname) );
		core::id::AtomID atomid( iatm, ires );
		cstatoms.push_back( atomid );
	}
	TR << reportline << std::endl;

	for ( Size i = 1; i <= cstatoms.size(); ++i ) {
		core::id::AtomID const& atomid = cstatoms[i];
		core::Vector const &xyz = pose->xyz( atomid );
		core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
		pose->add_constraint( core::scoring::constraints::ConstraintCOP(
			core::scoring::constraints::ConstraintOP(
			new core::scoring::constraints::CoordinateConstraint( atomid, atomid, xyz, fx )
			)));
	}
}

// final optimziation cycle with sidechain flexibility
void
GALigandDock::final_exact_cartmin(
	core::Size nneigh,
	LigandConformer & gene,
	core::pose::Pose &pose
	//bool dualrelax
) {

	/////
	// (1) setup movemap
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
	if ( pose.size() == 1 ) {
		mm->set_chi( true ); // ligand-only; virtual root
	} else {
		core::Size lig_jump = pose.fold_tree().get_jump_that_builds_residue( gene.ligand_id() );
		mm->set_chi( gene.ligand_id(), true );
		mm->set_jump( lig_jump, true );
		for ( core::Size j=1; j<=gene.moving_scs().size(); ++j ) {
			int resid = (int)gene.moving_scs()[j];
			mm->set_chi( resid, true ); // chi on just moving scs to reduce noise
			for ( int k=-((int)nneigh); k<=((int)nneigh); ++k ) {
				if ( resid+k < 1 || resid+k >= (int)pose.total_residue() ) continue;
				if ( !pose.residue(resid+k).is_protein() ) continue;
				//mm->set_chi( resid+k, true );
				mm->set_bb( resid+k, true );
			}
		}
	}

	core::kinematics::MoveMapOP mm2 = mm->clone(); // for torsion space

	// waters, hetmol
	for ( core::Size j=1; j<=pose.total_residue(); ++j ) {
		if ( (pose.residue_type(j).is_water() && move_water_)
				|| !pose.residue(j).is_protein() ) {
			core::Size wjump = pose.fold_tree().get_jump_that_builds_residue( j );
			mm->set_jump( wjump, true );
			mm2->set_jump( wjump, true ); // mm2 only jump
			mm->set_chi( j, true );
		}
	}

	std::vector< std::string > lines;

	// Torsion-relax in dual relax: run before any cst
	// hard coded internal coord part
	/*
	if ( dualrelax ) {
	TR << "Dual-relax mode: calling torsion-relax before cartrelax" << std::endl;
	protocols::relax::FastRelax relax_dual;
	relax_dual.set_scorefxn( scfxn_ ); // use grid score with softer farep
	mm2->set_bb( false ); // turn off bb
	lines.push_back( "switch:torsion" );
	lines.push_back( "repeat 2" );
	lines.push_back( "ramp_repack_min 0.02 0.01 1.0 50" );
	lines.push_back( "ramp_repack_min 1.0  0.00001 0.0 100" );
	lines.push_back( "accept_to_best" );
	lines.push_back( "endrepeat" );
	relax_dual.set_script_from_lines( lines );

	relax_dual.set_movemap( mm2 );
	relax_dual.apply( pose );
	lines.resize( 0 );
	}
	*/

	/////
	// (2) setup constraints if cst weight is on
	core::Size MINSEQSEP=2;
	core::Real MAXDIST=5.0, TOPOUT_WIDTH=2.0;
	core::Real TOPOUT_WT=1/(TOPOUT_WIDTH*TOPOUT_WIDTH); // max penalty per cst = 1
	core::Size addedCsts = 0;

	if ( scfxn_relax_->get_weight( core::scoring::atom_pair_constraint ) != 0 ) {
		TR << "Adding constraints to pose before final relax..." << std::endl;

		for ( core::Size j=1; j<pose.size(); ++j ) {
			for ( core::Size k=j+MINSEQSEP; k<pose.size(); ++k ) {
				if ( !mm->get_bb(j) && !mm->get_bb(k) ) continue;

				for ( core::Size jatm=1; jatm<pose.residue_type(j).nheavyatoms(); ++jatm ) {
					for ( core::Size katm=1; katm<pose.residue_type(k).nheavyatoms(); ++katm ) {

						core::Real dist = pose.residue(j).xyz(jatm).distance(
							pose.residue(k).xyz(katm)
						);

						if ( dist > MAXDIST ) continue;
						using namespace core::scoring::func;
						FuncOP fx( new ScalarWeightedFunc( 1.0, FuncOP( new TopOutFunc( TOPOUT_WT, dist, TOPOUT_WIDTH ) ) ) );
						pose.add_constraint(
							core::scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP(
							new core::scoring::constraints::AtomPairConstraint(
							core::id::AtomID(jatm,j), core::id::AtomID(katm,k), fx ) ) )
						);
						addedCsts++;
					}
				}
			}
		}
		TR << "Added " << addedCsts << " atmpair constraints to the pose." << std::endl;
	}
	if ( scfxn_relax_->get_weight( core::scoring::coordinate_constraint ) != 0 ) {
		TR << "Adding constraints to pose before final relax..." << std::endl;
		for ( core::Size j=1; j<pose.size(); ++j ) {
			if ( !mm->get_bb(j) ) continue;

			core::id::AtomID atomid(pose.residue(j).atom_index("CA"),j);
			core::Vector const &xyz = pose.xyz( atomid );
			using namespace core::scoring::func;
			FuncOP fx( new ScalarWeightedFunc( 1.0, FuncOP( new HarmonicFunc( 0.0, 1.0 ) ) ) );
			pose.add_constraint(
				core::scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP(
				new core::scoring::constraints::CoordinateConstraint
				( atomid, atomid, xyz, fx ) ) ) );
			addedCsts++;
		}
		TR << "Added " << addedCsts << " coord constraints to the pose." << std::endl;
	}
	if ( scfxn_relax_->get_weight( core::scoring::cart_bonded ) == 0 ) {
		scfxn_relax_->set_weight( core::scoring::cart_bonded, 0.5 );
		scfxn_relax_->set_weight( core::scoring::pro_close, 0.0 );
		TR << "scfxn_relax is not properly set for cartmin! setting cart_bonded=0.5 pro_close=0.0." << std::endl;
	}

	protocols::relax::FastRelax relax;
	if ( fast_relax_script_file_ != "" ) {
		TR << "==== Use FastRelax script: " << fast_relax_script_file_ << std::endl;
		relax=protocols::relax::FastRelax( scfxn_relax_, fast_relax_script_file_ );
	} else if ( fast_relax_lines_.size() > 0 ) {
		lines = fast_relax_lines_;
		relax.set_script_from_lines( lines );
		relax.set_scorefxn( scfxn_relax_ );
	} else {
		TR << "==== Use FastRelax hardcoded. "<< std::endl;
		lines.push_back( "switch:cartesian" );
		lines.push_back( "repeat 3" );
		lines.push_back( "ramp_repack_min 0.02 0.01 1.0 50" );
		lines.push_back( "ramp_repack_min 1.0  0.00001 0.0 50" );
		lines.push_back( "accept_to_best" );
		lines.push_back( "endrepeat" );
		relax.set_script_from_lines( lines );
		relax.set_scorefxn( scfxn_relax_ );
	}

	relax.set_movemap( mm );
	relax.apply( pose );
	TR << "final_cartmin: score after relax: ";
	TR << (*scfxn_relax_)(pose) <<std::endl;
	//scfxn_relax_->show( TR, pose );
	//pose.dump_pdb("dualcart.pdb");
}

// final optimziation cycle with sc flexibility only
void
GALigandDock::final_exact_scmin(
	LigandConformer const & gene,
	core::pose::Pose &pose
) {
	utility::vector1< core::Size > repack_scs = gene.moving_scs();

	// movemap for repacking
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( false );
	if ( pose.size() == 1 ) {
		mm->set_chi( true ); // ligand-only; virtual root
	} else {
		core::Size lig_jump = pose.fold_tree().get_jump_that_builds_residue( gene.ligand_id() );
		mm->set_jump( lig_jump, true );
		for ( core::Size j=0; j<=repack_scs.size(); ++j ) {
			core::Size resid = (j==0 ? gene.ligand_id() : repack_scs[j]);
			mm->set_chi( resid, true );
		}
	}

	// waters, hetmol
	for ( core::Size j=1; j<=pose.total_residue(); ++j ) {
		if ( (pose.residue_type(j).is_water() && move_water_)
				|| !pose.residue(j).is_protein() ) {
			//core::Size wjump = pose.fold_tree().get_jump_that_builds_residue( j );
			//mm->set_jump( wjump, true );
			// let's be more conservative...
			mm->set_chi( j, true );
		}
	}

	// opt-H: re-optimize full pose!
	//if( optimize_input_H_ ){
	TR << "Re-optimizing hydrogens in whole structure." << std::endl;
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
	task->initialize_from_command_line();
	task->or_optimize_h_mode( true );
	task->or_include_current( true );
	task->or_flip_HNQ( true );
	task->or_multi_cool_annealer( true );
	core::pack::pack_rotamers( pose, *scfxn_relax_, task );
	//}

	// main relax
	protocols::relax::FastRelax relax;
	std::vector< std::string > lines;
	if ( fast_relax_script_file_ != "" ) {
		TR << "==== Use FastRelax script: " << fast_relax_script_file_ << std::endl;
		relax = protocols::relax::FastRelax( scfxn_relax_, fast_relax_script_file_ );
	} else if ( fast_relax_lines_.size() > 0 ) {
		lines = fast_relax_lines_;
		relax.set_script_from_lines( lines );
		relax.set_scorefxn( scfxn_relax_ );
	} else {
		TR << "==== Use FastRelax hardcoded. " << std::endl;
		lines.push_back( "switch:torsion" );
		lines.push_back( "repeat 3" );
		lines.push_back( "ramp_repack_min 0.02 0.01 1.0 50" );
		lines.push_back( "ramp_repack_min 1.0  0.00001 0.0 50" );
		lines.push_back( "accept_to_best" );
		lines.push_back( "endrepeat" );
		relax.set_script_from_lines( lines );
		relax.set_scorefxn( scfxn_relax_ );
	}
	relax.set_movemap( mm );

	TR << "final_scmin: score after relax: " << (*scfxn_relax_)(pose) <<std::endl;
	//TR << "final_scmin:"; scfxn_relax_->show( TR, pose );
}


// final optimziation cycle with ligand flexibility only
void
GALigandDock::final_cartligmin(
	LigandConformer const & gene,
	core::pose::Pose &pose
) {
	// Nov08!! (turn this off for pre-Nov08)
	utility::vector1< core::Size > contact_scs;
	if ( redefine_flexscs_at_relax_ ) {
		contact_scs = get_atomic_contacting_sidechains( pose, gene.ligand_id(), 4.5 );
		TR << "Redefined flexible sidechains: ";
		for ( core::Size ires = 1; ires < contact_scs.size(); ++ires ) TR << contact_scs[ires] << "+";
		if ( contact_scs.size() > 0 ) TR << contact_scs[contact_scs.size()];
		TR << std::endl;
	} else {
		contact_scs = gene.moving_scs();
	}

	// setup for lig-only cartmin
	core::Real w_cart = (*scfxn_relax_)[ core::scoring::cart_bonded ];

	core::kinematics::MoveMapOP mmlig( new core::kinematics::MoveMap );
	mmlig->set_bb( false ); mmlig->set_chi( false ); mmlig->set_jump( false );
	if ( pose.size() == 1 ) {
		mmlig->set_chi( true ); // ligand-only; virtual root
	} else {
		core::Size lig_jump = pose.fold_tree().get_jump_that_builds_residue( gene.ligand_id() );
		mmlig->set_jump( lig_jump, true );
		mmlig->set_chi( gene.ligand_id(), true ); // ligand-chi only
		if ( min_neighbor_ ) {
			for ( core::Size j=1; j<=contact_scs.size(); ++j ) {
				mmlig->set_chi( contact_scs[j], true );
			}
		}
	}
	core::optimization::CartesianMinimizer minimizer;
	core::optimization::MinimizerOptions options( "lbfgs_armijo", 0.0001, true , false );
	options.max_iter(50);

	core::scoring::ScoreFunctionOP scfxn_cartmin = scfxn_relax_->clone();
	if ( w_cart < 1.0e-6 ) {
		scfxn_cartmin->set_weight( core::scoring::cart_bonded, 0.5 );
		scfxn_cartmin->set_weight( core::scoring::cart_bonded_ring, -0.5 );
		scfxn_cartmin->set_weight( core::scoring::pro_close, 0.0 );
	}
	//core::Real scorepremin = scfxn_cartmin->score( pose );
	minimizer.run( pose, *mmlig, *scfxn_cartmin, options );

	//core::Real scoremin = scfxn_cartmin->score( pose );
	//TR << "lig-cartmin, " << scorepremin << " -> " << scoremin << std::endl;
	//scfxn_cartmin->show( TR, pose );
}

// use WaterBoxMover to solvate ligand
void
GALigandDock::final_solvate(
	LigandConformer & gene,
	core::pose::Pose & pose
) {
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;

	core::pack::task::PackerTaskOP task_new = core::pack::task::TaskFactory::create_packer_task( pose );
	task_new->or_include_current(true);

	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		auto i_it = std::find( gene.moving_scs().begin(), gene.moving_scs().end(), i);
		if ( i == gene.ligand_id() || i_it != gene.moving_scs().end() ) {
			task_new->nonconst_residue_task(i).restrict_to_repacking();
		} else {
			task_new->nonconst_residue_task(i).prevent_repacking();
		}
	}

	core::scoring::ScoreFunctionOP sfwater = scfxn_->clone();
	sfwater->set_weight( core::scoring::pointwater, 1.0 );

	protocols::simple_moves::WaterBoxMover wb(sfwater);
	wb.set_taskop( task_new );
	wb.apply( pose );
}


// for multi-outputting, get the next pose
core::pose::PoseOP
GALigandDock::get_additional_output() {
	core::pose::PoseOP retval = remaining_outputs_.pop();
	if ( retval == nullptr ) return retval;

	(*scfxn_relax_)(*retval);
	return retval;
}

// load the initial inputs specified by initial_pool_
//   - reads all files ending in .pdb as PDBs
//   - reads special tag input as the input pose
//   - reads everything else as a silent file
void
GALigandDock::load_initial_pool(
	LigandConformer const &gene_initial,
	LigandConformers &genes_sel
) const {
	utility::vector1<std::string> input_pdbs = utility::string_split( initial_pool_, ',' );

	for ( core::Size ipdb = 1; ipdb <= input_pdbs.size(); ++ipdb ) {
		std::string tag = input_pdbs[ipdb];

		if ( (tag.length() >= 3) && (tag.substr( tag.length()-3 )=="pdb") ) {
			LigandConformer gene=gene_initial;
			gene.score( 0.0 );
			core::pose::PoseOP pose = core::import_pose::pose_from_file( tag, false, core::import_pose::PDB_file );
			gene.update_conf( pose );
			genes_sel.push_back( gene );

		} else if ( tag != "INPUT" && tag != "Input" && tag != "input" ) {
			// silent file
			core::io::silent::SilentFileOptions opts; // initialized from the command line
			core::io::silent::SilentFileData sfd( opts );
			sfd.read_file( tag );

			//fd: temporary hack for terminal type mismatches
			//core::pose::PoseOP pose_ref(new core::pose::Pose);
			//gene_initial.to_pose( pose_ref );

			for ( auto iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
				LigandConformer gene=gene_initial;
				gene.score( 0.0 );
				core::pose::PoseOP pose(new core::pose::Pose);
				iter->fill_pose( *pose );

				//fd: temporary hack for terminal type mismatches
				//pose->replace_residue( 1, pose_ref->residue(1), false );
				gene.update_conf( pose );

				genes_sel.push_back( gene );
			}
		} else {
			// tag == input
			LigandConformer gene=gene_initial;
			gene.score( 0.0 );
			genes_sel.push_back( gene );
		}
	}
}

// load the initial inputs specified by reference_pool_
//   - reads all files ending in .pdb as PDBs
//   - reads everything else as a silent file
//  we assume the reference pose ligand is the last residue (or should we use _all_ ligands?)
void
GALigandDock::load_reference_pool(
	LigandConformer const &/*gene_initial*/,
	utility::vector1< ConstraintInfo > & ref_ligs
) const {
	utility::vector1<std::string> ref_pdbs = utility::string_split( reference_pool_, ',' );

	for ( core::Size ipdb = 1; ipdb <= ref_pdbs.size(); ++ipdb ) {
		std::string tag = ref_pdbs[ipdb];

		if ( tag == "input" || tag == "Input" || tag == "INPUT" ) {
			// MOVED TO ALIGNER setup in order to share VS-info; do nothing in this case
			/*
			core::pose::PoseOP pose = gene_initial.receptor();
			scfxn_->score( *pose ); // make sure scored to get neighbor graph
			ref_ligs.push_back( ConstraintInfo(*pose, gridscorer, use_pharmacophore_,
			(ipdb==1) ) // report only at first case
			);
			*/

		} else if ( (tag.length() >= 3) && (tag.substr( tag.length()-3 )=="pdb") ) {
			core::pose::PoseOP pose = core::import_pose::pose_from_file( tag, false, core::import_pose::PDB_file );
			core::Size lig = pose->total_residue();
			if ( !pose->residue(lig).is_ligand() ) {
				utility_exit_with_message("error!  No ligand found in pdb file "+tag);
			}
			ref_ligs.push_back( ConstraintInfo(pose->residue(lig), use_pharmacophore_, (ipdb==1) ) );

		} else {
			// silent file
			core::io::silent::SilentFileOptions opts; // initialized from the command line
			core::io::silent::SilentFileData sfd( opts );
			sfd.read_file( tag );

			for ( auto iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
				core::pose::PoseOP pose(new core::pose::Pose);
				iter->fill_pose( *pose );
				core::Size lig = pose->total_residue();
				while ( pose->residue(lig).aa() == core::chemical::aa_vrt ) lig--;
				if ( !pose->residue(lig).is_ligand() ) {
					utility_exit_with_message("error!  No ligand found in silent file "+tag);
				}
				ref_ligs.push_back( ConstraintInfo(pose->residue(lig), use_pharmacophore_, (ipdb==1) ) );
			}
		}
	}
}

// LigandAligner
LigandAligner
GALigandDock::setup_ligand_aligner( core::pose::Pose const & pose,
	core::Size const lig_resno,
	utility::vector1< core::Size > movable_scs_in_ref // call by value
) const
{
	GridScorerOP gridscore_ref(new GridScorer( scfxn_ ));
	gridscore_ref->set_voxel_spacing( 0.5 ); // hard-coded

	// In dockPH mode, 0.02 prefers collapsing inside receptor
	if ( use_pharmacophore_ ) {
		gridscore_ref->set_w_rep( 0.1 );
		movable_scs_in_ref.resize( 0 ); // drop these from grid construction
	} else {
		gridscore_ref->set_w_rep( 0.02 ); // softer for reference-docking mode
		// keep input movable scs
	}

	gridscore_ref->set_smoothing( 0.75 ); // hard-coded
	gridscore_ref->set_bbox_padding( padding_ );
	gridscore_ref->set_hash_gridding( hashsize_ );
	gridscore_ref->set_hash_subgridding( subhash_ );
	gridscore_ref->set_exact( false );
	gridscore_ref->set_debug( false );
	gridscore_ref->set_out_of_bound_e( grid_bound_penalty_ );
	gridscore_ref->prepare_grid( pose, lig_resno );

	// pass all the ligand residues for grid construction
	utility::vector1< core::conformation::Residue > rsds_to_build_grids;
	rsds_to_build_grids.push_back( pose.residue(lig_resno) );
	if ( multiple_ligands_.size() > 0 ) {
		for ( core::Size ilig = 1; ilig <= multiple_ligands_.size(); ++ilig ) {
			core::conformation::ResidueOP ligand = core::conformation::get_residue_from_name( multiple_ligands_[ilig] );
			rsds_to_build_grids.push_back( *ligand );
		}
	}
	gridscore_ref->get_grid_atomtypes( rsds_to_build_grids );

	gridscore_ref->calculate_grid( pose, lig_resno, movable_scs_in_ref );

	LigandAligner aligner( use_pharmacophore_, movable_scs_in_ref, (runmode_ == "VSX"));
	aligner.set_sf( gridscore_ref );
	aligner.refine_input( (runmode_ == "refine") );

	if ( use_pharmacophore_ ) {
		// make sure reference_pool is not set to something else
		if ( !(reference_pool_ == "input" || reference_pool_ == "Input" ||
				reference_pool_ == "INPUT" ) ) {
			utility_exit_with_message("error!  pharmacophore docking requires reference_pool to be as 'INPUT'!" );
		}
		core::pose::PoseOP receptor( new core::pose::Pose( pose ) );
		receptor->delete_residue_slow( lig_resno );
		scfxn_->score( *receptor ); // for energygraph!!
		aligner.set_pharmacophore_reference( *receptor );
	}

	return aligner;
}

// initial perturbation
LigandConformers
GALigandDock::generate_perturbed_structures(
	LigandConformer const &gene_initial,
	GridScorerOP gridscorer,
	core::Size npool,
	LigandAligner aligner //call by value
) const {
	LigandConformers genes_sel, genes_ref, genes_rand;

	core::Real rmscut = protocol_[1].rmsthreshold;

	auto start = std::chrono::steady_clock::now();

	/////
	// 1: (optionally) load inputs and add to pool
	if ( initial_pool_ != "" ) {
		load_initial_pool(gene_initial, genes_sel);
	}
	core::Size nstruct_input = genes_sel.size(), nstruct_ref=0;
	int nleft = (int)npool - (int)nstruct_input;
	if ( nleft <= 0 ) {
		if ( reference_pool_ != "" ) {
			TR << "WARN! Reference pool provided but will not be used.  Increase pool size!" << std::endl;
		}
		return genes_sel;
	}

	// set up grid scorer for steps 2 & 3
	// if there is a smoothing schedule, use iter 1 smoothness to generate structures
	Real const w_rep_org( gridscorer->get_w_rep() );
	core::Real smoothing = protocol_[1].smoothing;
	gridscorer->set_smoothing( smoothing );
	gridscorer->set_w_rep( 0.02 ); // using softer repulsion for input generation

	/////
	// 2: (optionally) load reference structures and randomly generate conformers
	if ( reference_pool_ != "" ) {

		// make a temporary gridscorer for align docking; keep movable scs fixed in grid
		TR << "Construct a separate grid for LigandAligner, with 0.5 grid step and no movable scs." << std::endl;
		core::pose::PoseOP pose( new core::pose::Pose );
		gene_initial.to_pose( pose );

		utility::vector1< ConstraintInfo > ref_poses;
		load_reference_pool(gene_initial, ref_poses);

		// assign num structures generating from reference
		core::Size nstruct_ref = (core::Size)(nleft*reference_frac_+0.5);
		core::Size nrefgen = debug_? nstruct_ref : int(reference_oversample_)*nstruct_ref;

		// re-assign numbers based on Nmatches if using pharmacophore
		if ( reference_frac_auto_ && use_pharmacophore_ ) {
			nrefgen = aligner.estimate_nstruct_sample( pose->residue( gene_initial.ligand_id() ), nrefgen );
			nstruct_ref = (core::Size)(nrefgen/reference_oversample_);

			TR << "Automatically set nstruct-from-reference to " << nstruct_ref
				<< " (from " << nrefgen << " trials) of total " << nleft << " left to sample." << std::endl;
			TR << "Est. time for matching: " << nrefgen << "~" << nrefgen*2 << " seconds..." << std::endl;
		}

		for ( core::Size i=1; i<=nrefgen; ++i ) {
			LigandConformer gene( gene_initial );

			if ( ref_poses.size() > 0 ) {
				ConstraintInfo const & selected_ref =
					ref_poses[ numeric::random::rg().random_range( 1, ref_poses.size() ) ];
				aligner.set_target( selected_ref ); // random reference from pool
			}

			gene.randomize( 0.0 );  // for torsion&ring; fix trans
			aligner.apply( gene );

			// rescore with orignal gridscorer to match scale
			Real score_soft = gridscorer->score( gene, true );
			gene.score( score_soft );
			genes_ref.push_back( gene );
		}

		// take best scoring subset
		std::sort(genes_ref.begin(), genes_ref.end(),
			[&](LigandConformer const &lig_i, LigandConformer const &lig_j){ return lig_i.score() < lig_j.score(); } );
		for ( core::Size i=1; i<=nstruct_ref; ++i ) {
			genes_sel.push_back( genes_ref[i] );
		}
		TR << "Ref scorecut = " << genes_ref[nstruct_ref].score() << std::endl;
		nleft -= nstruct_ref;
	} // if


	if ( nleft <= 0 ) {
		return genes_sel;
	}

	/////
	// 3: random structures
	core::Size nrand = debug_? nleft : (int)(nleft*random_oversample_);
	for ( core::Size i=1; i<=nrand; ++i ) {
		LigandConformer gene( gene_initial );
		if ( runmode_ == "refine" ) {
			gene = mutate( gene_initial ); //mutation parameters are set as small in advance
		} else {
			gene.randomize( gridscorer->get_padding() - 1.0 );  // radius of search
		}

		Real score_soft = gridscorer->score( gene, true ); // score with soft repulsive
		gene.score( score_soft );
		genes_rand.push_back( gene );
	}

	// select lowest random structures by energy while ensuring diversity
	std::sort(genes_rand.begin(), genes_rand.end(),
		[&](LigandConformer const &lig_i, LigandConformer const &lig_j){ return lig_i.score() < lig_j.score(); } );
	utility::vector1< bool > selected( genes_rand.size(), false );
	for ( core::Size ii = 1; ii <= genes_rand.size(); ++ii ) {
		LigandConformer &gene_gen = genes_rand[ii];
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
	// no random structures left that are unique.  Fill with lowest-scoring
	for ( core::Size ii = 1; ii <= genes_rand.size() && genes_sel.size()<npool; ++ii ) {
		if ( !selected[ii] ) {
			genes_sel.push_back( genes_rand[ii] );
			selected[ii] = true;
		}
	}
	//std::chrono::duration<double> pack_time, min_time;
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> gen_time = end-start;

	TR << "Finished generating initial pool in " << gen_time.count() << " seconds." << std::endl;
	TR << "   # from inputs: " << nstruct_input << std::endl;
	TR << "   # aligned to ref: " << nstruct_ref << std::endl;
	TR << "   # random placement: " << nleft << std::endl;

	// recover original weight
	gridscorer->set_w_rep( w_rep_org );

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
	if ( tag->hasOption("scorefxn_relax") ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn_relax" ) );
		scfxn_relax_ = (datamap.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	} else {
		scfxn_relax_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );
	}

	// First read runmode; then allow options to change by user-specified options
	if ( tag->hasOption("runmode") ) {
		runmode_ = tag->getOption<std::string>("runmode");
		setup_params_for_runmode( runmode_ );
	}

	// Below are detailed controls
	if ( tag->hasOption("ngen") ) { ngen_ = tag->getOption<int>("ngen"); }

	// allowed movement
	if ( tag->hasOption("ligand") ) { ligid_ = tag->getOption<std::string>("ligand"); }
	if ( tag->hasOption("sidechains") ) { sidechains_ = tag->getOption<std::string>("sidechains"); }
	if ( tag->hasOption("sc_edge_buffer") ) { sc_edge_buffer_ = tag->getOption<core::Real>("sc_edge_buffer"); }
	if ( tag->hasOption("fa_rep_grid") ) { fa_rep_grid_ = tag->getOption<core::Real>("fa_rep_grid"); }
	if ( tag->hasOption("grid_bound_penalty") ) { grid_bound_penalty_ = tag->getOption<core::Real>("grid_bound_penalty"); }
	if ( tag->hasOption("rotprob") ) {  max_rot_cumulative_prob_ = tag->getOption<core::Real>("rotprob"); }
	if ( tag->hasOption("rotEcut") ) {  rot_energy_cutoff_ = tag->getOption<core::Real>("rotEcut"); }
	if ( tag->hasOption("favor_native") ) { favor_native_ = tag->getOption<core::Real>("favor_native"); }
	if ( tag->hasOption("optimize_input_H") ) { optimize_input_H_ = tag->getOption<bool>("optimize_input_H"); }

	if ( tag->hasOption("sample_ring_conformers") ) { sample_ring_conformers_ = tag->getOption<core::Real>("sample_ring_conformers"); }
	if ( tag->hasOption("altcrossover") ) { altcrossover_ = tag->getOption<core::Real>("altcrossover"); }

	// input params
	if ( tag->hasOption("use_pharmacophore") ) { use_pharmacophore_ = tag->getOption<bool>("use_pharmacophore"); }
	if ( tag->hasOption("initial_pool") ) { initial_pool_ = tag->getOption<std::string>("initial_pool"); }
	if ( tag->hasOption("reference_oversample") ) { reference_oversample_ = tag->getOption<core::Real>("reference_oversample"); }
	if ( tag->hasOption("reference_pool") ) { reference_pool_ = tag->getOption<std::string>("reference_pool"); }
	if ( tag->hasOption("reference_frac") ) { reference_frac_ = tag->getOption<core::Real>("reference_frac"); }
	if ( tag->hasOption("reference_frac_auto") ) { reference_frac_auto_ = tag->getOption<bool>("reference_frac_auto"); }
	if ( tag->hasOption("random_oversample") ) { random_oversample_ = tag->getOption<core::Real>("random_oversample"); }
	if ( tag->hasOption("premin_ligand") ) { premin_ligand_ = tag->getOption<bool>("premin_ligand"); }

	if ( tag->hasOption("multiple_ligands") ) {
		std::string ligands_string = tag->getOption<std::string>("multiple_ligands");
		multiple_ligands_ = utility::string_split( ligands_string, ',' );
		if ( multiple_ligands_.size() > 100 ) {
			TR.Error << "multiple_ligands arguments cannot be more than 100! " << std::endl;
			utility_exit();
		}
	}

	// reporting
	if ( tag->hasOption("nativepdb") ) {
		std::string nativepdb = tag->getOption<std::string>("nativepdb");
		pose_native_ = core::import_pose::pose_from_file( nativepdb, false, core::import_pose::PDB_file );
	}


	// post-processing
	if ( tag->hasOption("final_exact_minimize") ) {
		final_exact_minimize_ = tag->getOption<std::string>("final_exact_minimize");
		if ( final_exact_minimize_.substr(0,2) != "sc" && final_exact_minimize_.substr(0,4) != "bbsc"
				&& final_exact_minimize_ != "none" ) {
			TR.Error << "The tag 'final_exact_minimize' must be one of: sc, bbscX, rtmin, none" << std::endl;
			utility_exit();
		}
	}

	if ( tag->hasOption("cartmin_lig") ) { cartmin_lig_ = tag->getOption<bool>("cartmin_lig"); }
	if ( tag->hasOption("min_neighbor") ) { min_neighbor_ = tag->getOption<bool>("min_neighbor"); }
	if ( min_neighbor_ && !cartmin_lig_ ) {
		TR.Error << "min_neighbor may only be specified if cartmin_lig is enabled!" << std::endl;
		utility_exit();
	}

	if ( tag->hasOption("estimate_dG") ) { estimate_dG_ = tag->getOption<bool>("estimate_dG"); }

	if ( tag->hasOption("final_solvate") ) {
		final_solvate_ = tag->getOption<bool>("final_solvate");
		if ( final_exact_minimize_ == "none" ) {
			TR.Error << "The option 'final_solvate' requires a final minimize set!" << std::endl;
			utility_exit();
		}
	}

	if ( tag->hasOption("full_repack_before_finalmin") ) { full_repack_before_finalmin_ = tag->getOption<bool>("full_repack_before_finalmin"); }

	//protocols
	if ( tag->hasOption("fastrelax_script") ) {
		fast_relax_script_file_ = tag->getOption<std::string>("fastrelax_script");
	}
	if ( tag->hasOption("move_water") ) { move_water_ = tag->getOption<bool>("move_water"); }
	if ( tag->hasOption("redefine_flexscs_at_relax") ) { redefine_flexscs_at_relax_ = tag->getOption<bool>("redefine_flexscs_at_relax"); }

	// grid params
	if ( tag->hasOption("exact") ) { exact_ = tag->getOption<bool>("exact"); }
	if ( tag->hasOption("debug") ) { debug_ = tag->getOption<bool>("debug"); }
	if ( tag->hasOption("grid_step") ) { grid_ = tag->getOption<core::Real>("grid_step"); }
	if ( tag->hasOption("padding") ) { padding_ = tag->getOption<core::Real>("padding"); }
	if ( tag->hasOption("hashsize") ) { hashsize_ = tag->getOption<core::Real>("hashsize"); }
	if ( tag->hasOption("subhash") ) { subhash_ = tag->getOption<core::Real>("subhash"); }

	// per-cycle defaults
	if ( tag->hasOption("npool") ) { npool_ = tag->getOption<core::Size>("npool"); }

	if ( tag->hasOption("nrelax") ) nrelax_ = tag->getOption<core::Size>("nrelax");
	if ( tag->hasOption("nreport") ) nreport_ = tag->getOption<core::Size>("nreport");

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
	bool stage_specified( false );
	for ( tag_it = stage_tags.begin(); tag_it != stage_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "Stage" ) stage_specified = true;
	}

	if ( stage_specified ) {
		// remove info in case specified from runmode
		if ( protocol_.size() > 0 ) protocol_.resize( 0 );

		utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
		for ( tag_it = stage_tags.begin(); tag_it != stage_tags.end(); ++tag_it ) {
			if ( (*tag_it)->getName() == "Stage" ) {
				GADockStageParams stage_i( 1, npool_, rmsdthreshold_, pmut_, maxiter_, packer_cycles_, smoothing_, ramp_schedule_ );
				if ( (*tag_it)->hasOption("repeats") ) { stage_i.repeats = (*tag_it)->getOption<core::Size>("repeats"); }
				if ( (*tag_it)->hasOption("npool") ) { stage_i.pool = (*tag_it)->getOption<core::Size>("npool"); }
				if ( (*tag_it)->hasOption("pmut") ) { stage_i.pmut = (*tag_it)->getOption<core::Real>("pmut"); }
				if ( (*tag_it)->hasOption("rb_maxrank") ) { stage_i.rb_maxrank = (*tag_it)->getOption<core::Size>("rb_maxrank"); }
				if ( (*tag_it)->hasOption("smoothing") ) { stage_i.smoothing = (*tag_it)->getOption<core::Real>("smoothing"); }
				if ( (*tag_it)->hasOption("elec_scale") ) { stage_i.elec_scale = (*tag_it)->getOption<core::Real>("elec_scale"); }
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
	}

	// sanity checks
	if ( tag->hasOption("ngen") && protocol_.size() > 0 ) {
		TR.Error << "ngen specified but detailed schedule also given.  Aborting!" << std::endl;
		utility_exit();
	}

	// default protocol
	if ( protocol_.size() == 0 ) {
		TR << "Using default protocol." << std::endl;
		GADockStageParams stage_i( ngen_, npool_, rmsdthreshold_, pmut_, maxiter_, packer_cycles_, smoothing_, ramp_schedule_ ); // elec scale SHOULD be controled by per-stage schedule
		protocol_.push_back( stage_i );
	} else {
		TR << "Using custom " << protocol_.size() << "-stage protocol." << std::endl;
	}
}

// note that this only works through rosetta_scripts
void
GALigandDock::setup_params_for_runmode( std::string runmode )
{
	utility::vector1< core::Real > ramp_schedule;
	ramp_schedule.push_back( 0.1 ); ramp_schedule.push_back( 1.0 );
	protocol_.resize( 0 );
	fast_relax_lines_.resize( 0 );
	bool die( false );

	if ( runmode == "dockrigid" ) {
		GADockStageParams stage1( 10, 100, 1.0, // repeats, npool, rmsdthreshold
			0.2, 100, 100, 0.375, ramp_schedule ); // pmut, maxiter, packcycles, smoothing, ramp_schedule
		protocol_.push_back( stage1 );

		sidechains_ = "none";
		final_exact_minimize_ = "sc"; //default
		reference_pool_ = "input";
		nrelax_ = nreport_ = 20;
		// use default fast_relax schedule

	} else if ( runmode == "dockflex" ) {
		// hi-res cross docking scenario
		GADockStageParams stage1( 10, 100, 2.0, // repeats, npool, rmsdthreshold
			0.2, 100, 100, 0.375, ramp_schedule ); // pmut, maxiter, packcycles, smoothing, ramp_schedule
		protocol_.push_back( stage1 );

		sidechains_ = "aniso";
		final_exact_minimize_ = "bbsc1";
		sc_edge_buffer_ = 0.0; // in A
		reference_pool_ = "input";
		reference_frac_auto_ = true;
		nrelax_ = nreport_ = 20;

		fast_relax_lines_.push_back("switch:cartesian");
		fast_relax_lines_.push_back("repeat 1");
		fast_relax_lines_.push_back("ramp_repack_min 1.0   0.001  1.0 200");
		fast_relax_lines_.push_back("accept_to_best");
		fast_relax_lines_.push_back("endrepeat");

	} else if ( runmode.substr(0,2) == "VS" ) {
		// virtual screening
		final_exact_minimize_ = "sc"; //default
		reference_pool_ = "input";
		estimate_dG_ = true;
		reference_frac_auto_ = true;
		premin_ligand_ = true;

		if ( runmode.substr(2) == "H" ) { // VS-Hires
			sidechains_ = "aniso";
			nrelax_ = nreport_ = 20;
			//1-round minimization
			min_neighbor_ = true;
			/* // == default
			fast_relax_lines_.push_back("switch:torsion");
			fast_relax_lines_.push_back("repeat 3");
			fast_relax_lines_.push_back("ramp_repack_min 0.02  0.01   1.0 50");
			fast_relax_lines_.push_back("ramp_repack_min 1.0   0.000001  1.0 50");
			fast_relax_lines_.push_back("accept_to_best");
			fast_relax_lines_.push_back("endrepeat");
			*/

			// repeats, npool, rmsdthreshold, pmut, maxiter, packcycles, smoothing, ramp_schedule
			GADockStageParams stage1( 5, 100, 2.0, 0.2, 100, 100, 0.75, ramp_schedule );
			stage1.elec_scale = 3.0; // upweight at early stages
			stage1.maxiter = 25;
			protocol_.push_back( stage1 );
			GADockStageParams stage2( 5, 100, 2.0, 0.2, 100, 100, 0.375, ramp_schedule );
			stage2.maxiter = 25;
			protocol_.push_back( stage2 );

		} else if ( runmode.substr(2) == "X" ) { // VS-eXpress
			sidechains_ = "none";
			nrelax_ = 5;
			nreport_ = 1;

			//1-round minimization
			min_neighbor_ = false;
			fast_relax_lines_.push_back("switch:torsion");
			fast_relax_lines_.push_back("repeat 2");
			fast_relax_lines_.push_back("ramp_repack_min 0.02  0.01   1.0 50");
			fast_relax_lines_.push_back("ramp_repack_min 1.0   0.001  1.0 50");
			fast_relax_lines_.push_back("accept_to_best");
			fast_relax_lines_.push_back("endrepeat");

			//ramp_schedule.push_back( 0.1 ); ramp_schedule.push_back( 1.0 );
			GADockStageParams stage1( 5, 50, 2.0, 0.2, 50, 25, 0.375, ramp_schedule );
			protocol_.push_back( stage1 );

		} else {
			die = true;
		}
	} else if ( runmode == "refine" ) {
		//
	} else {
		die = true;
	}

	if ( die ) utility_exit_with_message("error!  No runmode defined for "+runmode);
}


GAOptimizerOP
GALigandDock::get_optimizer(
	LigandConformer const &gene_initial,
	GridScorerOP gridscorer
) const {
	GAOptimizerOP optimizer(new GAOptimizer(gridscorer));
	if ( pose_native_ ) {
		//fd  if the input ligand is the last residue, use the last residue of the _native_ as the ligand
		//fd  otherwise, match residue IDs
		core::Size native_lig = gene_initial.ligand_id();
		if ( gene_initial.is_ligand_terminal() ) native_lig = pose_native_->total_residue();
		LigandConformer gene_native( pose_native_, native_lig, gene_initial.moving_scs() );
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
	attlist + XMLSchemaAttribute( "scorefxn_relax", xs_string, "weights file");

	attlist + XMLSchemaAttribute( "runmode", xs_string, "run mode [dock/dockPH/refine/optligand]");

	attlist + XMLSchemaAttribute( "altcrossover", xsct_rosetta_bool, "Use alternate xover.");
	attlist + XMLSchemaAttribute( "sample_ring_conformers", xsct_rosetta_bool, "Allow ring conformer sampling if defined in params.");
	attlist + XMLSchemaAttribute( "rotprob", xsct_real, "max cumulative rotamer probability");
	attlist + XMLSchemaAttribute( "rotEcut", xsct_real, "rotamer 1b energy");
	attlist + XMLSchemaAttribute( "ligand", xs_string, "ligand residue id (if not specified default to last residue)");
	attlist + XMLSchemaAttribute( "nativepdb", xs_string, "name of native pdb");
	attlist + XMLSchemaAttribute( "favor_native", xsct_real, "give a bonus score to the input rotamer");
	attlist + XMLSchemaAttribute( "optimize_input_H", xsct_rosetta_bool, "do not optimize H at the begining (which is used for grid construction)");
	attlist + XMLSchemaAttribute( "grid_step", xsct_real, "Grid step (A) for grid-based scoring");
	attlist + XMLSchemaAttribute( "padding", xsct_real, "Padding (A) step for grid-based scoring");
	attlist + XMLSchemaAttribute( "hashsize", xsct_real, "Width of hash bins (A)");
	attlist + XMLSchemaAttribute( "subhash", xsct_non_negative_integer, "When scanning gridspace, subhash to this level");
	attlist + XMLSchemaAttribute( "nrelax", xsct_non_negative_integer, "Num. structs to run final minimize");
	attlist + XMLSchemaAttribute( "nreport", xsct_non_negative_integer, "Num. structs to report");
	attlist + XMLSchemaAttribute( "final_exact_minimize", xs_string, "Minimize the Genes by exact score after GA.");
	attlist + XMLSchemaAttribute( "cartmin_lig", xsct_rosetta_bool, "Cartmin ligand-only before and after final relax");
	attlist + XMLSchemaAttribute( "premin_ligand", xsct_rosetta_bool, "Cartmin ligand-only at the beginning");
	attlist + XMLSchemaAttribute( "min_neighbor", xsct_rosetta_bool, "If cartmin is enabled, also cartmin SCs before and after final relax.");
	attlist + XMLSchemaAttribute( "full_repack_before_finalmin", xsct_rosetta_bool, "Full repack before final relax.");
	attlist + XMLSchemaAttribute( "final_solvate", xsct_rosetta_bool, "Solvate pose (via WaterBoxMover) in final optimize. Default: false");
	attlist + XMLSchemaAttribute( "fastrelax_script", xs_string, "FastRelax script file for exact minimize.");
	attlist + XMLSchemaAttribute( "move_water", xsct_rosetta_bool, "Move water at final relaxation.");
	attlist + XMLSchemaAttribute( "redefine_flexscs_at_relax", xsct_rosetta_bool, "Redefine movable residues at final relaxation.");
	attlist + XMLSchemaAttribute( "exact", xsct_rosetta_bool, "Use exact scoring.");
	attlist + XMLSchemaAttribute( "debug", xsct_rosetta_bool, "Debug grid scoring: report both exact and grid scores.");
	attlist + XMLSchemaAttribute( "use_pharmacophore", xsct_rosetta_bool, "Use pharmacophore info at initial pool generation.");
	attlist + XMLSchemaAttribute( "initial_pool", xs_string, "Include these structures in the initial pool.");
	attlist + XMLSchemaAttribute( "multiple_ligands", xs_string, "Scan ligands with these residue types.");
	attlist + XMLSchemaAttribute( "random_oversample", xsct_real, "scale factor to ntrial of initial random pool generation");
	attlist + XMLSchemaAttribute( "reference_oversample", xsct_real, "scale factor to ntrial of initial reference pool generation");
	attlist + XMLSchemaAttribute( "reference_pool", xs_string, "Use this structures as _references_ to generate the initial pool.");
	attlist + XMLSchemaAttribute( "reference_frac", xsct_real, "If reference pool is provided, the fraction of structures from the reference pool.");
	attlist + XMLSchemaAttribute( "reference_frac_auto", xsct_rosetta_bool, "Select Nstruct to sample by reference automatically");
	attlist + XMLSchemaAttribute( "sidechains", xs_string, "Sidechains to move: none, auto, or residue IDs.");
	attlist + XMLSchemaAttribute( "sc_edge_buffer", xsct_real, "Scaling factor of maxdistance when deciding to include sc as movable");
	attlist + XMLSchemaAttribute( "fa_rep_grid", xsct_real, "Repulsion weight at grid scoring stage");
	attlist + XMLSchemaAttribute( "grid_bound_penalty", xsct_real, "Penalty factor when ligand atm gets out of boundary");
	attlist + XMLSchemaAttribute( "estimate_dG", xsct_rosetta_bool, "Estimate dG of binding on lowest-energy docked pose. Default: false");

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
		+ XMLSchemaAttribute( "elec_scale", xsct_real, "Scale of elec and hbond terms at this stage")
		+ XMLSchemaAttribute( "pmut", xsct_real, "Sampling frequency weight for this template" )
		+ XMLSchemaAttribute( "rb_maxrank", xs_integer, "superimpose sampled pose to motifs of topN parent structures" )
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

