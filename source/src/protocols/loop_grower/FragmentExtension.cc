// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Brandon Frenz
/// @author Frank DiMaio

//keep first
#include <protocols/loop_grower/FragmentExtensionCreator.hh>
#include <protocols/loop_grower/FragmentExtension.hh>
#include <protocols/loop_grower/FragmentExtension.fwd.hh>
#include <protocols/loop_grower/LoopGrower.hh>
#include <fstream>
#include <iostream>

#include <protocols/hybridization/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/relax/util.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>

#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/comparative_modeling/coord_util.hh>
#include <basic/datacache/DataMap.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/Energies.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/selection.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/id/AtomID.hh>

#include <core/pose/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/USOGFunc.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/conformation/util.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/BBTorsionSRFD.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>


#include <utility/excn/Exceptions.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/graph/Graph.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/WeightedSampler.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>



namespace protocols {
namespace loop_grower {

using namespace core;

static THREAD_LOCAL basic::Tracer TR("protocols.loop_grower.FragmentExtension");
//static numeric::random::RandomGenerator RG(8403178);

/////////////
// creator
std::string
FragmentExtensionCreator::keyname() const {
	return FragmentExtension::mover_name();
}

std::string fragext_subelement_ct_name( std::string const & name ) {
	return "FragExt_subelement_" + name + "Type";
}

protocols::moves::MoverOP
FragmentExtensionCreator::create_mover() const {
	return protocols::moves::MoverOP( new FragmentExtension );
}

std::string
FragmentExtension::mover_name() {
	return "FragmentExtension";
}
////////////

FragmentExtension::FragmentExtension( ) {
	beamwidth_ = 1;
	fragtrials_ = 200;
	chainbreak_ = 1.5;
	continuous_weight_ = 1.0;
	rmscutoff_ = 0.0;
	fragcluster_ = 0.0;
	master_beam_cutoff_ = 0.0;
	sheetbonus_ = 0.5;
	sheet_tolerance_ = 0.7;
	sc_scale_ = 1;
	fa_bonus_ = 0.0;
	window_dens_weight_ = 30;
	master_beam_width_ = beamwidth_;
	rmswindow_ = 20;
	steps_ = 1e4;
	parallelcount_ = 0;
	montecarlorounds_ = 1;
	sheetcriteria_ = 2;
	beamscorecutoff_ = 0.85;
	debug_ = false;
	fragmelt_ = 2;
	minmelt_ = 15;
	pack_min_cycles_ = 2;
	looporder_ = 0;
	maxloopsize_ = 8;
	sf_ = core::scoring::get_score_function();
	cen_sf_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth");
	cenrot_sf_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_cenrot_relax");
	direction_ = 1;
	parametercheck_= false;
	minimize_= true;
	nativegrow_= false;
	dumpbeam_ = false;
	dumperrors_ = false;
	dumpfinalbeam_ = false;
	dumprms_ = false;
	cenrot_ = false;
	greedy_ = true;
	readbeams_ = false;
	writebeams_ = true;
	writelps_ = true;
	fafilter_ = true;
	famin_ = false;
	samplesheets_ = true;
	trackfragments_ = false;
	filterprevious_ = false;
	rephasemap_ = false;
	checksymm_ = false;
	continuous_sheets_ = true;
	auto_stop_ = false;
	clustercheck_ = false;
	rescorebeams_ = false;
	storedbeams_ = "";
	filterbeams_ = "";
	coordfile_ = "";
	skeleton_file_ = "";
	read_from_file_ = false;
	assign_incomplete_ = false;

}

protocols::moves::MoverOP FragmentExtension::clone() const { return protocols::moves::MoverOP( new FragmentExtension( *this ) ); }
protocols::moves::MoverOP FragmentExtension::fresh_instance() const { return protocols::moves::MoverOP( new FragmentExtension ); }

void
FragmentExtension::apply( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;

	// native
	if ( option[ in::file::native ].user() ) {
		native_ = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_file( *native_, option[ in::file::native ](), true );

	}
	if ( !native_ && nativegrow_ == true ) {
		TR << " attemping to use custom model torsions when no model is provided" << std::endl;
		runtime_assert( native_ && nativegrow_ != true);
	}
	//autogenerate fragments
	if ( fragments_.size() == 0 || dumpfragments_ ) {
		if ( frag_sizes_.size() == 0 ) frag_sizes_.push_back(3);
		for ( core::uint i = 1; i <= frag_sizes_.size(); ++i ) {
			fragments_.push_back( protocols::hybridization::create_fragment_set_no_ssbias(fullength_seq_->sequence(), frag_sizes_[i], nfrags_) );
		}
	}

	if ( fragcluster_ != 0 ) {
		cluster_fragments(fragments_, fragcluster_);
	}

	core::pose::Pose pose_for_seq;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_for_seq, false);
	core::sequence::SequenceOP t_pdb_seq( new core::sequence::Sequence( pose_for_seq.sequence(), "pose_seq" ));
	core::sequence::SWAligner sw_align;
	core::sequence::ScoringSchemeOP ss(  new core::sequence::SimpleScoringScheme(120, 0, -100, 0));
	core::sequence::SequenceAlignment fasta2template_;

	fasta2template_ = sw_align.align(fullength_seq_, t_pdb_seq, ss);
	core::id::SequenceMapping sequencemap = fasta2template_.sequence_mapping(1,2);

	TR << "Aligning:" << std::endl;
	TR << fasta2template_ << std::endl;
	for ( core::Size ii = 1; ii <= fullength_seq_->sequence().size(); ++ii ) {
		if ( sequencemap[ii] ) {
			TR << fullength_seq_->at(ii);
		} else {
			TR << "-";
		}
	}
	TR << std::endl;

	utility::vector1<protocols::loops::Loop> loops_from_aln = get_unaligned(sequencemap), loops_to_build;


	// split chainbreaks
	for ( int i=1; i<=(int)loops_from_aln.size() ; ++i ) {
		bool loop_was_split=false;
		for ( int j=1; j<=(int)cbreaks_.size(); ++j ) {
			if ( loops_from_aln[i].start() <= cbreaks_[j] && loops_from_aln[i].stop() > cbreaks_[j] ) {
				runtime_assert( !loop_was_split );
				loop_was_split = true;
				loops_to_build.push_back(protocols::loops::Loop(loops_from_aln[i].start(), cbreaks_[j]));
				loops_to_build.push_back(protocols::loops::Loop(cbreaks_[j]+1, loops_from_aln[i].stop()));
			}
		}
		if ( !loop_was_split ) {
			loops_to_build.push_back(loops_from_aln[i]);
		}
	}

	//fix alignment errors if the first and last residue are the same and it assigns it incorrectly
	for ( core::Size i=1; i<=loops_to_build.size(); i++ ) {
		bool loopagain = false;
		core::Size start_fasta = loops_to_build[i].start();
		core::Size stop_fasta = loops_to_build[i].stop();
		bool done = false;
		while ( !done ) {
			if ( start_fasta != 1 && fullength_seq_->at(start_fasta-1) == fullength_seq_->at(stop_fasta) ) {
				core::Vector startxyz = pose.residue(sequencemap[start_fasta-1]).atom("N").xyz();
				core::Vector start_negone_xyz = pose.residue(sequencemap[start_fasta-2]).atom("C").xyz();
				Real distance = (startxyz-start_negone_xyz).length();
				TR << "distance = " << distance << std::endl;

				if ( distance > 2 ) {
					sequencemap[start_fasta-1] = 0;
					sequencemap[stop_fasta] = sequencemap[stop_fasta+1]-1;
					loops_to_build[i].set_start(start_fasta-1);
					loops_to_build[i].set_stop(stop_fasta-1);
					loopagain = true;
					TR << "loops to build start and stop in loop " << loops_to_build[i].start() << " " << loops_to_build[i].stop() << std::endl;
				}

			}

			start_fasta = loops_to_build[i].start();
			stop_fasta = loops_to_build[i].stop();
			if ( loopagain ) {
				done=false;
			} else {
				done=true;
			}
			loopagain = false;
		}
	}

	if ( dumpfragments_ ) {
		for ( Size i=1; i<=frag_sizes_.size(); i++ ) {
			std::string totalfrags = utility::to_string(nfrags_);
			std::string fraglength = utility::to_string(frag_sizes_[i]);
			std::string fragname = totalfrags+"."+fraglength+"mers";
			core::fragment::FragmentIO().write_data( fragname, *(fragments_[i]) );
			exit(0);
		}
	}


	if ( loops_to_build.size() <= 1 ) greedy_ = true;
	if ( !greedy_ ) looporder_ = -2; //this insures the ordering of the loops is from N term to C term which is currently required for the proper growing in loop comparator
	if ( looporder_ == 0 || looporder_ == -1 ) numeric::random::random_permutation( loops_to_build.begin(), loops_to_build.end(), numeric::random::rg() );

	int start_fasta, stop_fasta, start_pose, stop_pose;
	bool singleloop = false;
	int i=0;
	TR << "loops_to_build.size() = \t" << loops_to_build.size() << std::endl;
	LoopComparator comparator;
	comparator.set_read_from_file(read_from_file_);
	comparator.set_assign_incomplete(assign_incomplete_);
	comparator.set_loops( loops_to_build );
	core::pose::Pose originalpose = pose;
	core::scoring::ScoreFunctionOP sf_clashing;
	core::scoring::ScoreFunctionOP sf_dens;
	if ( read_from_file_ ) {
		greedy_ = false;
		pack_min_cycles_ = 0;
	}
	if ( !greedy_ ) {
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			sf_dens = core::scoring::ScoreFunctionOP( new core::scoring::symmetry::SymmetricScoreFunction );
			sf_clashing = core::scoring::ScoreFunctionOP( new core::scoring::symmetry::SymmetricScoreFunction );
			core::scoring::symmetry::SymmetricScoreFunctionOP symmdens( new core::scoring::symmetry::SymmetricScoreFunction );
			sf_dens = utility::pointer::dynamic_pointer_cast< core::scoring::symmetry::SymmetricScoreFunction >(symmdens);
			core::scoring::symmetry::SymmetricScoreFunctionOP symmclash( new core::scoring::symmetry::SymmetricScoreFunction );
			sf_clashing = utility::pointer::dynamic_pointer_cast< core::scoring::symmetry::SymmetricScoreFunction >(symmclash);
		} else {
			sf_dens = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
			sf_clashing = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
		}
		sf_dens->set_weight( core::scoring::elec_dens_fast , cen_sf_->get_weight(core::scoring::elec_dens_fast) );
		//sf_dens->set_weight( core::scoring::vdw , cen_sf_->get_weight(core::scoring::vdw) ); //vdw needs to be included if we want to get 1 body and 2 body in a single pass.
		sf_clashing->set_weight( core::scoring::vdw , cen_sf_->get_weight(core::scoring::vdw) );
		comparator.set_density_sf(sf_dens);
		comparator.set_vdw_sf(sf_clashing);
		comparator.set_sequence(fullength_seq_);
	}

	//hacky way to turn off growing if reading in loops from text file
	if ( read_from_file_ ) singleloop = true;
	while ( i<=(int)loops_to_build.size()-1 && !singleloop ) {
		i++;
		if ( looporder_ > 0 ) {
			i = looporder_;
			singleloop = true;
			int loopsizecheck = loops_to_build.size();
			if ( looporder_ > loopsizecheck ) {
				TR << "You have assigned an invalid loop. Loops are numbered from N to C terminus please use a different looporder value" << std::endl;
				runtime_assert( looporder_ <= loopsizecheck );
			}
		}
		TR << "reseting start & stop" << std::endl;
		start_fasta = loops_to_build[i].start();
		stop_fasta = loops_to_build[i].stop();

		if ( looporder_ == -1 ) {
			if ( (unsigned)stop_fasta - (unsigned)start_fasta > maxloopsize_ ) continue;
		}

		start_pose = stop_pose = 0;

		TR << "finished checking " << std::endl;
		bool nterm=(start_fasta == 1);
		for ( int j=1; j<=(int)cbreaks_.size() && !nterm; ++j ) {
			nterm = (start_fasta == (int)cbreaks_[j]+1);
		}
		bool cterm=(stop_fasta == (int)fullength_seq_->sequence().size());
		for ( int j=1; j<=(int)cbreaks_.size() && !nterm; ++j ) {
			cterm = (stop_fasta == (int)cbreaks_[j]);
		}

		if ( !nterm ) start_pose = sequencemap[start_fasta-1];
		if ( !cterm ) stop_pose = sequencemap[stop_fasta+1];
		if ( (unsigned)start_pose > pose.total_residue() && stop_pose == 0 ) {
			TR << "not growing from virtual" << std::endl;
			continue;
		}
		Size loopnumber = i;
		if ( looporder_ > 0 ) loopnumber  = looporder_;

		TR << " Growing loop between residues " << start_fasta-1 << " and " << stop_fasta+1 << std::endl;
		TR << " loops to build start and stop " << loops_to_build[i].start() << " " << loops_to_build[i].stop() << " start and stop pose " << start_pose << " " << stop_pose << std::endl;
		LoopGrower loopgrow( fullength_seq_, loops_to_build[i], start_pose, stop_pose,  fragments_ );
		loopgrow.set_scorefunction( sf_ );
		loopgrow.set_cen_scorefunction( cen_sf_ );
		loopgrow.set_cenrot_scorefunction( cenrot_sf_);
		loopgrow.set_beamwidth( beamwidth_ );
		loopgrow.set_fragtrials( fragtrials_ );
		loopgrow.set_chainbreak( chainbreak_ );
		loopgrow.set_continuous_weight( continuous_weight_ );
		loopgrow.set_debug( debug_ );
		loopgrow.set_fragmelt( fragmelt_ );
		loopgrow.set_minmelt( minmelt_ );
		loopgrow.set_dumpbeam ( dumpbeam_ );
		loopgrow.set_dumperrors( dumperrors_ );
		loopgrow.set_dumpfinalbeam( dumpfinalbeam_ );
		loopgrow.set_dumprms( dumprms_ );
		loopgrow.set_pack_min_cycles( pack_min_cycles_ );
		loopgrow.set_direction ( direction_ );
		loopgrow.set_minimize ( minimize_ );
		loopgrow.set_native( native_ );
		loopgrow.set_nativegrow ( nativegrow_ );
		loopgrow.set_parametercheck( parametercheck_ );
		loopgrow.set_cenrot( cenrot_ );
		loopgrow.set_writebeams( writebeams_ );
		loopgrow.set_readbeams( readbeams_ );
		loopgrow.set_clustercheck( clustercheck_ );
		loopgrow.set_rescorebeams( rescorebeams_ );
		loopgrow.set_writelps( writelps_ );
		loopgrow.set_fafilter( fafilter_ );
		loopgrow.set_famin( famin_ );
		loopgrow.set_samplesheets( samplesheets_ );
		loopgrow.set_trackfragments( trackfragments_ );
		loopgrow.set_filterprevious( filterprevious_ );
		loopgrow.set_rephasemap( rephasemap_ );
		loopgrow.set_checksymm( checksymm_ );
		loopgrow.set_continuous_sheets( continuous_sheets_ );
		loopgrow.set_auto_stop( auto_stop_ );
		loopgrow.set_rmscutoff ( rmscutoff_ );
		loopgrow.set_master_beam_cutoff ( master_beam_cutoff_ );
		loopgrow.set_sheetbonus ( sheetbonus_ );
		loopgrow.set_sheet_tolerance ( sheet_tolerance_ );
		loopgrow.set_sc_scale ( sc_scale_ );
		loopgrow.set_windowdensweight( window_dens_weight_ );
		loopgrow.set_master_beam_width( master_beam_width_ );
		loopgrow.set_rmswindow( rmswindow_ );
		loopgrow.set_steps( steps_ );
		loopgrow.set_pcount( parallelcount_ );
		loopgrow.set_sheetcriteria( sheetcriteria_ );
		loopgrow.set_beamscorecutoff ( beamscorecutoff_ );
		loopgrow.set_greedy( greedy_ );
		loopgrow.set_storedbeams( storedbeams_ );
		loopgrow.set_filterbeams( filterbeams_ );
		loopgrow.set_coordfile( coordfile_ );
		loopgrow.set_skeletonfile( skeleton_file_ );
		loopgrow.set_loopnumber( loopnumber );

		loopgrow.apply(pose);

		if ( !greedy_ && master_beam_width_ != 1337 ) {
			LoopPartialSolutionStore solutionset = loopgrow.getSolutionSet();
			comparator.push(solutionset);
			pose = originalpose;
		}

		// update sequencemap
		core::pose::symmetry::extract_asymmetric_unit(pose, pose_for_seq, false);
		t_pdb_seq = core::sequence::SequenceOP( new core::sequence::Sequence( pose_for_seq.sequence(), "pose_seq" ));
		fasta2template_ = sw_align.align(fullength_seq_, t_pdb_seq, ss);
		sequencemap = fasta2template_.sequence_mapping(1,2); //creates the sequence map
	}
	if ( !greedy_ && master_beam_width_ != 1337 ) {
		// The comparator uses the loop scores from the grower and it's modified score function. Clashing is determined via vdw scores in centroid hence the conversion here.
		protocols::simple_moves::SwitchResidueTypeSetMover tocen("centroid");
		tocen.apply(pose);
		Real pre_grow_clash = (*sf_clashing)(pose);
		TR << "calling fill pose" << std::endl;
		comparator.fill_pose(pose);
		comparator.set_scores(pose);

		protocols::moves::MoverOP restore_sc;
		restore_sc = protocols::simple_moves::ReturnSidechainMoverOP( new protocols::simple_moves::ReturnSidechainMover( pose ));
		TR << "applying comparator " << std::endl;

		//the creation of comparator pose is only necessary if you want to use the option to remove a clashing loop. This is currrently disabled but I'm leaving the logic here for now
		pose.conformation().reset_chain_endings();
		core::pose::Pose comparatorpose = pose;
		core::pose::Pose bestpose = pose;
		Real least_clashing = pre_grow_clash + 99999;
		for ( Size i=1; i<=montecarlorounds_; i++ ) {
			TR << " comparator pose size is " << comparatorpose.total_residue() << std::endl;
			comparator.apply(comparatorpose);
			(*cen_sf_)(comparatorpose);
			Real clashscore = (*sf_clashing)(comparatorpose);
			if ( clashscore < least_clashing ) {
				least_clashing = clashscore;
				bestpose = comparatorpose;
			}
			core::pose::PDBInfoOP newpdbinfo(new core::pose::PDBInfo( comparatorpose, true ));
			comparatorpose.pdb_info( newpdbinfo );
			comparatorpose.dump_pdb("aftercomparator_"+utility::to_string(i)+"_"+utility::to_string(clashscore)+".pdb");
			comparatorpose = pose;
		}
		Real finalclash = least_clashing-pre_grow_clash;
		TR << " finalclash " << finalclash << " pre_grow_clash " << pre_grow_clash << std::endl;
		bool needs_another_round = false;
		if ( finalclash > 100.0 ) needs_another_round = true;
		std::ofstream outfile;
		outfile.open("recommendation.txt");
		if ( needs_another_round ) {
			outfile << " I suggest another round of growing the clash score is " << finalclash << std::endl;
		} else {
			outfile << " Another round probably isn't necessary. Clash score is " << finalclash << std::endl;
		}
		protocols::simple_moves::SwitchResidueTypeSetMover to_full(chemical::FA_STANDARD);
		to_full.apply(pose);

	}
}


utility::vector1<protocols::loops::Loop>
FragmentExtension::get_unaligned( core::id::SequenceMapping const & sequencemap ) const
{
	utility::vector1<protocols::loops::Loop> loops;

	int startpos=0;
	bool inloop=false;
	for ( unsigned int currentpos = 1; currentpos <= sequencemap.size1(); ++currentpos ) {

		bool ismapped=false;
		if ( sequencemap[currentpos]==0 ) {
			ismapped = false;
		} else {
			ismapped = true;
		}
		if ( !inloop && !ismapped ) {
			startpos = currentpos;
		}

		if ( inloop && ismapped ) {
			loops.push_back(protocols::loops::Loop(startpos,currentpos-1));
		}
		inloop = !ismapped;
	}
	if ( inloop ) {
		loops.push_back(protocols::loops::Loop(startpos, sequencemap.size1()));
	}
	return loops;
}



// parse_my_tag
void
FragmentExtension::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const & ,
	moves::Movers_map const & ,
	core::pose::Pose const & /*pose*/ )
{
	std::string fastaname = tag->getOption<std::string>("fasta");
	utility::vector1<std::string> seq = core::sequence::read_fasta_file_str(fastaname);
	TR << "Read "  << seq.size() << " sequences" << std::endl;

	// read chainbreaks
	std::string fulllength_clean;
	for ( int i=0; i<(int)seq[1].length(); ++i ) {
		if ( seq[1][i] == '/' ) {
			cbreaks_.push_back(fulllength_clean.length());
		} else {
			fulllength_clean += seq[1][i];
		}
	}
	fullength_seq_ = core::sequence::SequenceOP(new core::sequence::Sequence( fulllength_clean, "target" ));

	// fragments
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "Fragments" ) {
			using namespace core::fragment;
			fragments_.push_back( FragmentIO().read_data( (*tag_it)->getOption<std::string>( "fragfile" )  ) );
		}
	}
	//runtime_assert( fragments_.size()<=1 );  // only 1 fragset for now
	dumpfragments_ = false;

	// scorefunctions
	if ( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		sf_ = ( data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ) )->clone();
	}
	if ( tag->hasOption( "censcorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "censcorefxn" ) );
		cen_sf_ = ( data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ) )->clone();
	}
	if ( tag->hasOption( "cenrotscorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "cenrotscorefxn" ) );
		cenrot_sf_ = ( data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ) )->clone();
	}
	if ( tag->hasOption( "beamwidth" ) ) {
		beamwidth_ = tag->getOption<core::Size  >( "beamwidth" );
	}
	if ( tag->hasOption( "fragtrials" ) ) {
		fragtrials_ = tag->getOption<core::Size>( "fragtrials" );
	}
	if ( tag->hasOption( "debug" ) ) {
		debug_ = tag->getOption<bool>( "debug" );
	}
	if ( tag->hasOption( "chainbreak" ) ) {
		chainbreak_ = tag->getOption<core::Real>( "chainbreak" );
	}
	if ( tag->hasOption( "continuous_weight" ) ) {
		continuous_weight_ = tag->getOption<core::Real>( "continuous_weight" );
	}
	if ( tag->hasOption( "rmscutoff" ) ) {
		rmscutoff_ = tag->getOption<core::Real>( "rmscutoff" );
	}
	if ( tag->hasOption( "fragcluster" ) ) {
		fragcluster_ = tag->getOption<core::Real>( "fragcluster" );
	}
	if ( tag->hasOption( "masterbeamcutoff" ) ) {
		master_beam_cutoff_ = tag->getOption<core::Real>( "masterbeamcutoff" );
	}
	if ( tag->hasOption( "sheetbonus" ) ) {
		sheetbonus_ = tag->getOption<core::Real>( "sheetbonus" );
	}
	if ( tag->hasOption( "sheet_tolerance" ) ) {
		sheet_tolerance_ = tag->getOption<core::Real>( "sheet_tolerance" );
	}
	if ( tag->hasOption( "sc_scale" ) ) {
		sc_scale_ = tag->getOption<core::Real>( "sc_scale" );
	}
	if ( tag->hasOption( "fabonus" ) ) {
		fa_bonus_ = tag->getOption<core::Real>( "fabonus" );
	}
	if ( tag->hasOption( "windowdensweight" ) ) {
		window_dens_weight_ = tag->getOption<core::Real>( "windowdensweight" );
	}
	if ( tag->hasOption( "masterbeamwidth" ) ) {
		master_beam_width_ = tag->getOption<core::Size>( "masterbeamwidth" );
	}
	if ( tag->hasOption( "rmswindow" ) ) {
		rmswindow_ = tag->getOption<core::Size>( "rmswindow" );
	}
	if ( tag->hasOption( "steps" ) ) {
		steps_ = tag->getOption<core::Size>( "steps" );
	}
	if ( tag->hasOption( "pcount" ) ) {
		parallelcount_ = tag->getOption<core::Size>( "pcount" );
	}
	if ( tag->hasOption( "sheetcriteria" ) ) {
		sheetcriteria_ = tag->getOption<core::Size>( "sheetcriteria" );
	}
	if ( tag->hasOption( "comparatorrounds" ) ) {
		montecarlorounds_ = tag->getOption<core::Size>( "comparatorrounds" );
	}
	if ( tag->hasOption( "beamscorecutoff" ) ) {
		beamscorecutoff_ = tag->getOption<core::Real>( "beamscorecutoff" );
	}
	if ( tag->hasOption( "pack_min_cycles" ) ) {
		pack_min_cycles_ = tag->getOption<core::Size>( "pack_min_cycles" );
	}
	if ( tag->hasOption( "dumpbeam" ) ) {
		dumpbeam_ = tag->getOption<bool>( "dumpbeam" );
	}
	if ( tag->hasOption( "dumperrors" ) ) {
		dumperrors_ = tag->getOption<bool>( "dumperrors" );
	}
	if ( tag->hasOption( "dumpfinalbeam" ) ) {
		dumpfinalbeam_ = tag->getOption<bool>( "dumpfinalbeam" );
	}
	if ( tag->hasOption( "dumprms" ) ) {
		dumprms_ = tag->getOption<bool>( "dumprms" );
	}
	if ( tag->hasOption( "dumpfragments" ) ) {
		dumpfragments_ = tag->getOption<bool>( "dumpfragments" );
	}
	if ( tag->hasOption( "direction" ) ) {
		direction_ = tag->getOption<core::Size>( "direction" );
	}
	if ( tag->hasOption( "minimize" ) ) {
		minimize_ = tag->getOption<bool>( "minimize" );
	}
	if ( tag->hasOption( "nativegrow" ) ) {
		nativegrow_ = tag->getOption<bool>( "nativegrow" );
	}
	if ( tag->hasOption( "parametercheck" ) ) {
		parametercheck_ = tag->getOption<bool>( "parametercheck" );
	}
	if ( tag->hasOption( "cenrot" ) ) {
		cenrot_ = tag->getOption<bool>( "cenrot" );
	}
	if ( tag->hasOption( "writebeams" ) ) {
		writebeams_ = tag->getOption<bool>( "writebeams" );
	}
	if ( tag->hasOption( "readbeams" ) ) {
		readbeams_ = tag->getOption<bool>( "readbeams" );
	}
	if ( tag->hasOption( "clustercheck" ) ) {
		clustercheck_ = tag->getOption<bool>( "clustercheck" );
	}
	if ( tag->hasOption( "rescorebeams" ) ) {
		rescorebeams_ = tag->getOption<bool>( "rescorebeams" );
	}
	if ( tag->hasOption( "writelps" ) ) {
		writelps_ = tag->getOption<bool>( "writelps" );
	}
	if ( tag->hasOption( "fafilter" ) ) {
		fafilter_ = tag->getOption<bool>( "fafilter" );
	}
	if ( tag->hasOption( "famin" ) ) {
		famin_ = tag->getOption<bool>( "famin" );
	}
	if ( tag->hasOption( "samplesheets" ) ) {
		samplesheets_ = tag->getOption<bool>( "samplesheets" );
	}
	if ( tag->hasOption( "filterprevious" ) ) {
		filterprevious_ = tag->getOption<bool>( "filterprevious" );
	}
	if ( tag->hasOption( "rephasemap" ) ) {
		rephasemap_ = tag->getOption<bool>( "rephasemap" );
	}
	if ( tag->hasOption( "checksymm" ) ) {
		checksymm_ = tag->getOption<bool>( "checksymm" );
	}
	if ( tag->hasOption( "continuous_sheets" ) ) {
		continuous_sheets_ = tag->getOption<bool>( "continuous_sheets" );
	}
	if ( tag->hasOption( "auto_stop" ) ) {
		auto_stop_ = tag->getOption<bool>( "auto_stop" );
	}
	greedy_ = true;
	if ( tag->hasOption( "greedy" ) ) {
		greedy_ = tag->getOption<bool>( "greedy" );
	}
	if ( tag->hasOption( "read_from_file" ) ) {
		read_from_file_ = tag->getOption<bool>( "read_from_file" );
	}
	if ( tag->hasOption( "assign_incomplete" ) ) {
		assign_incomplete_ = tag->getOption<bool>( "assign_incomplete" );
	}
	if ( tag->hasOption( "looporder" ) ) {
		looporder_ = tag->getOption<int>( "looporder" );
	}
	if ( tag->hasOption( "maxloopsize" ) ) {
		maxloopsize_ = tag->getOption<int>( "maxloopsize" );
	}

	nfrags_ = tag->getOption< core::Size >( "nfrags", 25 );
	// autofragments
	if ( tag->hasOption( "fraglens" ) ) {
		utility::vector1<std::string> fraglens = utility::string_split( tag->getOption< std::string >( "fraglens" ), ',' );
		for ( core::Size i=1; i<=fraglens.size(); ++i ) {
			frag_sizes_.push_back( atoi(fraglens[i].c_str()) );
		}
	}
	//runtime_assert( frag_sizes_.size()<=1 );  // only 1 fragset for now

	// melt
	if ( tag->hasOption( "fragmelt" ) ) {
		fragmelt_ = tag->getOption<core::Size  >( "fragmelt" );
	}
	if ( tag->hasOption( "minmelt" ) ) {
		minmelt_ = tag->getOption<core::Size  >( "minmelt" );
	}
	if ( tag->hasOption( "storedbeams" ) ) {
		storedbeams_ = tag->getOption<std::string  >( "storedbeams" );
	}
	if ( tag->hasOption( "filterbeams" ) ) {
		filterbeams_ = tag->getOption<std::string  >( "filterbeams" );
	}
	if ( tag->hasOption( "coordfile" ) ) {
		coordfile_ = tag->getOption<std::string  >( "coordfile" );
	}
	if ( tag->hasOption( "skeletonfile" ) ) {
		skeleton_file_ = tag->getOption<std::string  >( "skeletonfile" );
	}
}

void FragmentExtensionCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FragmentExtension::provide_xml_schema( xsd );
}

void FragmentExtension::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//basic direct attributes
	attlist
		//
		+ XMLSchemaAttribute::attribute_w_default( "fasta", xs_string, "The name of the fasta file", "")
		+ XMLSchemaAttribute::attribute_w_default( "beamwidth", xsct_non_negative_integer, "Controls the maximize size of the beam search", "32" )
		+ XMLSchemaAttribute::attribute_w_default( "fragtrials", xsct_non_negative_integer, "The maximum number of fragments to apply at each extension.", "200" )
		+ XMLSchemaAttribute::attribute_w_default( "debug", xsct_rosetta_bool, "Triggers all the debuging options", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "chainbreak", xsct_real, "The weight on the chainbreak term when closing the loop", "1.5" )
		+ XMLSchemaAttribute::attribute_w_default( "continuous_weight", xsct_real, "The weight on the non-continuous density penalty", "0.3" )
		+ XMLSchemaAttribute::attribute_w_default( "rms_cutoff", xsct_real, "The rms cutoff for two structures to be considered different", "1.5" )
		+ XMLSchemaAttribute::attribute_w_default( "fragcluster", xsct_real, "The torsional RMS for two fragments to be considered different. Default doesn't filter", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "masterbeamcutoff", xsct_real, "The RMS different for two structures to be considered in a different path", "3.0" )
		+ XMLSchemaAttribute::attribute_w_default( "sheetbonus", xsct_real, "The percentage of the score to add for the sampled sheets that fit the density", "0.5" )
		+ XMLSchemaAttribute::attribute_w_default( "sheet_tolerance", xsct_real, "How good do the sheets need to be in order to be accepted", "0.7" )
		+ XMLSchemaAttribute::attribute_w_default( "sc_scale", xsct_real, "Scales down the side chain density score", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "windowdensweight", xsct_real, "Density weight on the windowdens", "30" )
		+ XMLSchemaAttribute::attribute_w_default( "masterbeamwidth", xsct_non_negative_integer, "Cap on the number of beams within a cluster", "200" )
		+ XMLSchemaAttribute::attribute_w_default( "rmswindow", xsct_real, "The number of residues to include in the ", "20" )
		+ XMLSchemaAttribute::attribute_w_default( "rmscutoff", xsct_real, "How close two conformations must be to be considered identical", "1.5" )
		+ XMLSchemaAttribute::attribute_w_default( "steps", xsct_non_negative_integer, "The number of steps to take before exiting.", "10000" )
		+ XMLSchemaAttribute::attribute_w_default( "pcount", xsct_non_negative_integer, "The tag for the parallelcount. Purely used for managing jobs through the python script", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "sheetcriteria", xsct_non_negative_integer, "The strategy to use for deciding whether or not to grow sheets", "2" )
		+ XMLSchemaAttribute::attribute_w_default( "comparatorrounds", xsct_non_negative_integer, "The number of models to produce with the loop comparator", "100" )
		+ XMLSchemaAttribute::attribute_w_default( "beamscorecutoff", xsct_real, "The number of models to produce with the loop comparator", "0.85" )
		+ XMLSchemaAttribute::attribute_w_default( "pack_min_cycles", xsct_non_negative_integer, "The number of cycles of packing and minimization to use in full atom mode", "2" )
		+ XMLSchemaAttribute::attribute_w_default( "dumpbeam", xsct_rosetta_bool, "Dump the structures after filtering", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "dumperrors", xsct_rosetta_bool, "Dump the poses when the beam deviates from the native", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "dumpfinalbeam", xsct_rosetta_bool, "Dump the poses after the last step", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "dumprms", xsct_rosetta_bool, "Dumps the rms of the beam for debugging purposes", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "dumpfragments", xsct_rosetta_bool, "Writes the fragments as a fragfile after picking", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "direction", xsct_non_negative_integer, "Which side of the missing segment to start growing from.", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize", xsct_rosetta_bool, "Whether or not to minimize the structures", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "nativegrow", xsct_rosetta_bool, "Use the native structure as a fragment when growing", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "parametercheck", xsct_rosetta_bool, "Uses the native structures to figure out what the settings would need to be to capture the native", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "cenrot", xsct_rosetta_bool, "Use the cenrot representation", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "writebeams", xsct_rosetta_bool, "Writes 'beamfiles' used in parallelization through python", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "readbeams", xsct_rosetta_bool, "Reads the beams from the file", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "clustercheck", xsct_rosetta_bool, "Dumps all the structures in the neighborhood of the native for debugging", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rescorebeams", xsct_rosetta_bool, "Rescores all the beams in the beamfile. Useful when testing score functions", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "writelps", xsct_rosetta_bool, "Writes an lps file when the loop is finished that can be input for the monte carlo assembler", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "fafilter", xsct_rosetta_bool, "Rescores the top 2N structures with the full atom representation before filtering down to N.", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "famin", xsct_rosetta_bool, "Toggles minimization in full atom", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "samplesheets", xsct_rosetta_bool, "Toggles the sheet sampler", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "filterprevious", xsct_rosetta_bool, "Filter against the beams in the filterbeams file. Used in taboo search", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rephasemap", xsct_rosetta_bool, "Triggers rephasing of the map when working with xray data, currently distabled", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "checksymm", xsct_rosetta_bool, "Checks all the symmetric versions of the connecting segments to see if an alternative is a better match", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "continuous_sheets", xsct_rosetta_bool, "Toggles a one time vs maintained bonus for the sheets", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "auto_stop", xsct_rosetta_bool, "Toggles autostoping when running out of density. Currently unfinished", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "greedy", xsct_rosetta_bool, "Only runs on the loop specificed by the loop order variable", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "read_from_file", xsct_rosetta_bool, "Read the partial models directly from the file and run the monte carlo assembly", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "assign_incomplete", xsct_rosetta_bool, "Have the option to assign no segment to a given peice in the structure when using Monte Carlo Assembly", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "looporder", xsct_non_negative_integer, "Which loop to build. Labelted 1-N from N-C termini", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "maxloopsize", xsct_non_negative_integer, "If loop order set to -1 all loops shorter than this will be skipped", "8" )
		+ XMLSchemaAttribute( "nfrags", xsct_nnegative_int_cslist, "The number of fragments to pick")
		+ XMLSchemaAttribute( "fraglens", xsct_nnegative_int_cslist, "Fragment Lengths To pick")
		+ XMLSchemaAttribute::attribute_w_default( "fragmelt", xsct_non_negative_integer, "How overlap of the old fragment to replace at every step. fraglen-fragmelt = # of new residues added", "2" )
		+ XMLSchemaAttribute::attribute_w_default( "minmelt", xsct_non_negative_integer, "How overlap of the old fragment to replace at every step", "15" )
		+ XMLSchemaAttribute::attribute_w_default( "storedbeams", xs_string, "The name of the file containing the input beams", "")
		+ XMLSchemaAttribute::attribute_w_default( "filterbeams", xs_string, "The name of the file containing the beams just used in filtering", "")
		+ XMLSchemaAttribute::attribute_w_default( "coordfile", xs_string, "The name of the file with coordinates and score weights, currently not fully developed", "")
		+ XMLSchemaAttribute::attribute_w_default( "skeletonfile", xs_string, "The name of the skeleton density map", "")
		//score functions
		+ XMLSchemaAttribute( "scorefxn", xs_string, "Full Atom Scorefunction")
		+ XMLSchemaAttribute( "censcorefxn", xs_string, "Centroid Scorefunction")
		+ XMLSchemaAttribute( "cenrotscorefxn", xs_string, "Scorefunction for use with cenrot. This is typically not used");

	// attributes for Fragments subelement
	AttributeList fragment_subelement_attributes;
	fragment_subelement_attributes
		+ XMLSchemaAttribute( "fragfile", xs_string, "The name of the fragment file");

	//rosetta_scripts::attributes_for_parse_task_operations(DetailedControls_subelement_attributes);

	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & fragext_subelement_ct_name );
	subelements.add_simple_subelement( "Fragments", fragment_subelement_attributes, "Name of the fragment file");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "This is XML that calls RosettaES (aka the loop grower)", attlist, subelements );
}

void LoopComparator::set_scores(core::pose::Pose & pose){


	//This code calculates and stores the 1 body (electron density) and 2 body (van der walls) energies for all the possible solutions.
	//each gap in the protein is refered to as a "loop" and each individual solution for that gap is the refered to as the "beam".
	//The LoopPartialSolutionStore contains all the beams for an individual loop.

	TR << "setting scores" << std::endl;
	//Store 1 body loopscores
	for ( Size i=1; i<= solutionsets_.size(); i++ ) {
		LoopPartialSolutionStore loop = solutionsets_[i];
		utility::vector1<Real> one_body_vector(loop.size());
		Size upper = loop.get_upper_fasta();
		Size lower = loop.get_lower_fasta();
		for ( Size ii=1; ii<=loop.size(); ii++ ) {
			LoopPartialSolution beam = loop[ii];
			Real beamscore = beam.get_score();
			//Note one body score is the per residue score
			Real score = beamscore/(upper-lower);
			TR << "one body score is " << score << " beam size is " << upper-lower << std::endl;

			//This code calculates the one body scores directly. It's not being used as we are using the scores stored by the loop grower instead.
			/*beam.set_coordinates(pose);
			Real dens_score = 0.0;
			(*sf_dens_)(pose);
			for(Size j=lower; j<=upper; j++){
			core::scoring::ScoreType elec_dens = core::scoring::score_type_from_name( "elec_dens_fast" );
			core::Real residue_dens = (sf_dens_->get_weight(core::scoring::elec_dens_fast) * pose.energies().residue_total_energies(j)[elec_dens]);
			dens_score += residue_dens;
			}
			//TR << "one body " << dens_score << std::endl;*/

			one_body_vector[ii] = score;
		}
		one_body_energies_.push_back(one_body_vector);
	}
	TR << "finished 1 body scores" << std::endl;
	//Store energies
	//// using this method calculates 1 body and 2 body energies at the same time. This requires the electron density be scored over additional residues but is probably faster than scoring
	//an additional time for each loop to change this behavior set first_round = to false and uncomment the above code.
	for ( Size i=1; i<=solutionsets_.size(); i++ ) {
		LoopPartialSolutionStore loop = solutionsets_[i];
		//iterate over all the n+1 loops
		Size z = i+1;
		Size upper_one = loop.get_upper_fasta();
		Size lower_one = loop.get_lower_fasta();
		TR << "starting solutionset_" << i << "'s two body scores " << std::endl;
		while ( z<=solutionsets_.size() ) {
			utility::vector1<utility::vector1< Real > > two_body_beam_vector(loop.size());
			LoopPartialSolutionStore compared_loop = solutionsets_[z];
			Size upper_two = compared_loop.get_upper_fasta();
			Size lower_two = compared_loop.get_lower_fasta();
			for ( Size ii=1; ii<=loop.size(); ii++ ) {
				utility::vector1<Real> two_body_score_vector(compared_loop.size());
				LoopPartialSolution beam_one = loop[ii];
				beam_one.set_coordinates(pose);
				TR << " solutionset_"<< i << " compared to solutionset " << z << " using beam " << ii << std::endl;
				//iterate over all the beams in the second loop
				for ( Size iii=1; iii<=compared_loop.size(); iii++ ) {
					LoopPartialSolution beam_two = compared_loop[iii];
					//beam_two.apply(pose, lower_one, upper_one);
					beam_two.set_coordinates(pose);
					(*sf_vdw_)(pose);

					//get second loop score
					core::scoring::EnergyGraph const & energy_graph( pose.energies().energy_graph() );
					Real score = 0;
					for ( Size j=lower_one; j<=upper_one; ++j ) {
						for ( utility::graph::Graph::EdgeListConstIter
								iru = energy_graph.get_node(j)->const_upper_edge_list_begin(),
								irue = energy_graph.get_node(j)->const_upper_edge_list_end();
								iru != irue; ++iru ) {
							core::scoring::EnergyEdge const *edge( static_cast< core::scoring::EnergyEdge const *> (*iru) );
							Size const e1( edge->get_first_node_ind() );
							Size const e2( edge->get_second_node_ind() );
							if ( (e1 >=lower_one && e1<=upper_one) || (e2 >= lower_one && e2 <= upper_one) ) {
								if ( (e1 >=lower_one && e1 <=upper_one) || (e2 >=lower_one && e2 <= upper_one) ) {
									if ( (e1 >= lower_two && e1<= upper_two) || (e2 >=lower_two && e2<= upper_two) ) {
										score += (*edge)[core::scoring::vdw];
									}
								}
							}
						}
					}
					two_body_score_vector[iii] = score;
					//TR << "two_body_score " << two_body_score_vector[iii] << std::endl;
				}
				two_body_beam_vector[ii] = two_body_score_vector;
			}
			std::pair<Size,Size> mapkey(i,z);
			two_body_energies_map_[mapkey] = two_body_beam_vector;
			z++;
		}
	}
}

void LoopComparator::apply( core::pose::Pose & pose ){

	utility::vector1<Size> working_solution;
	utility::vector1<Size> best_solution;
	utility::vector1<Size> current_solution;
	Real working_score = 0;
	Real best_score = 99999;

	for ( Size i=1; i<=solutionsets_.size(); i++ ) {
		Size beam_number = solutionsets_[i].size();
		Size beam = numeric::random::random_range(1, beam_number);
		working_solution.push_back(beam);
	}
	for ( Size i=1; i<=working_solution.size(); i++ ) {
		Size beam = working_solution[i];
		working_score += one_body_energies_[i][beam];
		Size z = i+1;
		while ( z<=working_solution.size() ) {
			std::pair<Size,Size> mapkey(i,z);
			Real two_body_score = two_body_energies_map_[mapkey][working_solution[i]][working_solution[z]];
			working_score += two_body_score;
			z++;
		}
	}
	best_score = working_score;
	best_solution = working_solution;
	current_solution = working_solution;


	Size montecount = 0;
	Real kt = 200;
	for ( Size i=1; i<=6; i++ ) {
		kt = kt/2;
		for ( Size j=1; j<=250; j++ ) {
			montecount++;
			Size loop = numeric::random::random_range(1, solutionsets_.size());
			//change the range of the following function to start from 0 if you want the possibility of no loop being assigned
			Size start_range = 1;
			if ( assign_incomplete_ ) {
				start_range = 0;
			}
			Size beam = numeric::random::random_range(start_range, solutionsets_[loop].size());
			current_solution[loop] = beam;
			Real current_score = 0;
			Real two_body_score = 0;
			for ( Size ii=1; ii<=current_solution.size(); ii++ ) {
				Size storedbeam = current_solution[ii];
				if ( storedbeam == 0 ) continue;
				Size looplen = solutionsets_[ii].get_upper_fasta() - solutionsets_[ii].get_lower_fasta();
				current_score += one_body_energies_[ii][storedbeam]/looplen;
				Size z = ii;
				while ( z<=current_solution.size()-1 ) {
					z++;
					std::pair<Size,Size> mapkey(ii,z);
					if ( current_solution[z] == 0 ) continue;
					two_body_score = two_body_energies_map_[mapkey][current_solution[ii]][current_solution[z]];
					current_score += two_body_score;
				}
			}
			if ( current_score < best_score ) {
				best_score =  current_score;
				best_solution = current_solution;
			}
			bool accept = false;
			Real const arg = (current_score-working_score)/kt;
			Real probability = exp(-(arg));
			if ( probability > 1.0 ) probability = 1.0;
			if ( probability >= numeric::random::rg().uniform() ) {
				accept = true;
			}
			if ( accept ) {
				working_score = current_score;
				working_solution = current_solution;
			} else {
				current_score = working_score;
				current_solution = working_solution;
			}
		}
	}
	//TR << "best solution size " << best_solution.size() << std::endl;
	//TR << best_solution << std::endl;
	Size deletecount = 0;

	one_body_score_ = 0;
	for ( Size i=1; i<=best_solution.size(); i++ ) {
		LoopPartialSolutionStore loop = solutionsets_[i];
		if ( best_solution[i] != 0 ) {
			LoopPartialSolution beam = loop[best_solution[i]];
			beam.set_coordinates(pose);
			one_body_score_ += one_body_energies_[i][best_solution[i]];
		}
	}
	for ( Size i=1; i<=best_solution.size(); i++ ) {
		LoopPartialSolutionStore loop = solutionsets_[i];
		if ( best_solution[i] == 0 ) {
			Size start = loops_[i].start()-deletecount;
			Size end = loops_[i].stop()-deletecount;
			TR << "deleting loop " << start << " to " << end << std::endl;
			pose.delete_residue_range_slow(start, end);
			deletecount += end-start+1;
		}
	}
	TR << " best score " << best_score << " one body score part is " << one_body_score_ << std::endl;

}
void LoopComparator::read_from_file(){
	std::ifstream lpsfile("lpsfile.txt");
	Size loopcount;
	lpsfile >> loopcount;
	for ( Size i=1; i<=loopcount; i++ ) {
		Size totalbeams;
		lpsfile >> totalbeams;
		LoopPartialSolutionStore loop;
		Size lower_fasta, upper_fasta, lower_pose, upper_pose, cut_point;
		lpsfile >> lower_fasta >> upper_fasta >> lower_pose >> upper_pose >> cut_point;
		if ( lpsfile.eof() ) break;
		loop.set_fastas(lower_fasta, upper_fasta);
		loop.set_poses(lower_pose, upper_pose);
		loop.set_cutpoint(cut_point);
		for ( Size ii=1; ii<=totalbeams; ii++ ) {
			LoopPartialSolution beam;
			utility::vector1< core::id::AtomID > ids;
			utility::vector1< numeric::xyzVector<core::Real> > positions;
			Size m;
			Real score;
			lpsfile >> score;
			beam.set_score(score);
			lpsfile >> m;
			if ( lpsfile.eof() ) break;
			for ( Size j=1; j<=m; j++ ) {
				Size atomno, rsd;
				lpsfile >> atomno >> rsd;
				core::id::AtomID atomid(atomno,rsd);
				ids.push_back(atomid);
				if ( lpsfile.eof() ) break;
			}
			numeric::xyzVector<core::Real> xyzvector;
			for ( Size j=1; j<=m; j++ ) {
				lpsfile >> xyzvector[0] >> xyzvector[1] >> xyzvector[2];
				positions.push_back(xyzvector);
				if ( lpsfile.eof() ) break;
			}
			beam.set_positions(positions);
			beam.set_ids(ids);
			loop.store(beam);
			if ( lpsfile.eof() ) break;
		}
		//  std::cout << "testing shit " << n << " lower_fasta  " << lower_fasta << " upper_fasta " << upper_fasta << " lower_pose " << lower_pose << " upper_pose " << upper_pose << std::endl;
		solutionsets_.push_back(loop);
	}
}

void LoopComparator::fill_pose(core::pose::Pose & pose){

	//create residues from sequence
	core::chemical::ResidueTypeSetCOP restypeset_cen = pose.residue_type_set_for_pose();
	core::chemical::ResidueTypeCOPs restypes_pose_cen = core::pose::residue_types_from_sequence(fullength_seq_->sequence(), *restypeset_cen);
	core::chemical::ResidueTypeSetCOP restypeset = pose.residue_type_set_for_pose();
	core::chemical::ResidueTypeCOPs restypes_pose = core::pose::residue_types_from_sequence(fullength_seq_->sequence(), *restypeset);

	bool is_centroid;
	if ( pose.is_centroid() ) {
		is_centroid = true;
	} else {
		is_centroid = false;
	}
	if ( read_from_file_ ) read_from_file();
	Size cutpoint;
	//Place Loops into Pose
	bool is_cterm = false;
	bool is_nterm = false;
	for ( Size i=1; i<=solutionsets_.size(); i++ ) {
		is_nterm = false;
		is_cterm = false;
		LoopPartialSolutionStore loop = solutionsets_[i];
		Size range_hi = loop.get_upper_fasta();
		//This lower pose logic is probably broken it should be 1 lower than it is
		Size lower_pose = loop.get_lower_pose();
		Size upper_pose = loop.get_upper_pose();
		cutpoint = loop.get_cutpoint();
		if ( range_hi == fullength_seq_->sequence().size() ) is_cterm = true;
		if ( lower_pose == 1 ) is_nterm = true;
		if ( !is_nterm && !is_cterm ) {
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, lower_pose-1  );
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, lower_pose );
			core::conformation::remove_upper_terminus_type_from_conformation_residue( pose.conformation(), lower_pose-1 );
			core::conformation::remove_lower_terminus_type_from_conformation_residue( pose.conformation(), lower_pose );
		}
		for ( Size j=upper_pose; j>cutpoint; j-- ) {
			if ( is_centroid ) {
				core::conformation::ResidueOP newres_cen( core::conformation::ResidueFactory::create_residue(*restypes_pose_cen[j] ));
				pose.conformation().safely_prepend_polymer_residue_before_seqpos(*newres_cen, lower_pose, true);
			} else {
				core::conformation::ResidueOP newres( core::conformation::ResidueFactory::create_residue(*restypes_pose[j] ));
				pose.conformation().safely_prepend_polymer_residue_before_seqpos(*newres, lower_pose, true);
			}
		}
		for ( Size j=lower_pose; j<=cutpoint; j++ ) {
			if ( is_centroid ) {
				core::conformation::ResidueOP newres_cen( core::conformation::ResidueFactory::create_residue(*restypes_pose_cen[j] ));
				pose.conformation().safely_append_polymer_residue_after_seqpos(*newres_cen, j-1, true);
			} else {
				core::conformation::ResidueOP newres( core::conformation::ResidueFactory::create_residue(*restypes_pose[j] ));
				pose.conformation().safely_append_polymer_residue_after_seqpos(*newres, j-1, true);
			}
		}
	}
	for ( core::Size i=1; i<=solutionsets_.size(); i++ ) {
		LoopPartialSolutionStore loop = solutionsets_[i];
		for ( int j=1; j<=(int)pose.fold_tree().num_jump(); j++ ) {
			if ( pose.fold_tree().cutpoint_by_jump(j) == (int)loop.get_cutpoint() ) {
				core::kinematics::FoldTree ft_closed = pose.fold_tree();
				ft_closed.delete_jump_and_intervening_cutpoint(j);
				pose.fold_tree(ft_closed);
				break;
			}
		}
	}
	//pose.dump_pdb("afterfill.pdb");
	TR << "post fill fold tree " << pose.fold_tree() << std::endl;

}
void
FragmentExtension::cluster_fragments( utility::vector1<core::fragment::FragSetOP> & fragments, Real fragfilter ){


	Size nfragsets = fragments.size();
	for ( int i=1; i<=(int)nfragsets; ++i ) {
		core::fragment::FragSetOP newfrags = fragments[i]->empty_clone();
		for ( core::fragment::ConstFrameIterator j = fragments[i]->begin(); j != fragments[i]->end(); ++j ) {
			core::fragment::Frame frame = **j;
			core::fragment::FrameOP filteredframe = frame.clone();
			for ( core::Size i_frag=1; i_frag<=frame.nr_frags(); i_frag++ ) {
				core::fragment::FragDataCOP oldfrag = frame.fragment_ptr(i_frag);
				bool addfrag=true;
				core::Size storedfragcount = filteredframe->nr_frags();
				for ( core::Size j_frag=1; j_frag<=storedfragcount && addfrag; j_frag++ ) {
					core::fragment::FragDataCOP storedfrag = filteredframe->fragment_ptr(j_frag);
					core::Real error = 0;
					core::Size ntors = 0;
					for ( core::Size i_res=1; i_res<=oldfrag->size(); ++i_res ) {
						core::fragment::BBTorsionSRFD const & frag_i =
							dynamic_cast< core::fragment::BBTorsionSRFD const &> ( *(oldfrag->get_residue(i_res)) );
						core::fragment::BBTorsionSRFD const & frag_j =
							dynamic_cast< core::fragment::BBTorsionSRFD const &> ( *(storedfrag->get_residue(i_res)) );
						for ( core::Size i_tors=1; i_tors<=frag_i.nbb(); ++i_tors ) {
							core::Real err_i = (frag_i.torsion(i_tors) - frag_j.torsion(i_tors));
							error += err_i*err_i;
							ntors++;
						}
					}
					core::Real RMS = std::sqrt(error/ntors);
					if ( RMS <= fragfilter ) addfrag=false;
				}

				if ( addfrag ) filteredframe->add_fragment(oldfrag);
			}
			newfrags->add(filteredframe);
		}
		fragments[i]=newfrags;
	}
}

}
}
