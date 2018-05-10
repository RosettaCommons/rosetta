// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/GAOptimizer.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#include <protocols/ligand_docking/GALigandDock/GAOptimizer.hh>
#include <protocols/ligand_docking/GALigandDock/LigandConformer.hh>

#include <core/scoring/rms_util.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/lkball/LK_BallEnergy.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <chrono>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

static basic::Tracer TR( "protocols.ligand_docking.GALigandDock.GAOptimizer" );

GAOptimizer::GAOptimizer( GridScorerOP grid ) {
	scorefxn_ = grid;
}

GAOptimizer::~GAOptimizer(){}

// main GA optimization loop
void
GAOptimizer::run( LigandConformers & genes ) {
	show_status( genes, "input pool (scored with soft_rep)" );

	LigandConformers genes_new;
	for ( core::Size i=1; i<=protocol_.size(); ++i ) {
		GADockStageParams const &stage_i = protocol_[i];

		// im not sure about this behavior or not.
		// without this statement, repeats 0 means do a min (only) and rerank at this smoothing
		//   but there is a pretty big overhead in doing so
		// for now, we will take iter=0 to mean "do nothing"
		if ( stage_i.repeats == 0 ) continue;

		TR << "=====================================================================================" << std::endl;
		TR << "  Stage " << i << " : " << std::endl;
		TR << "    iterations = " << stage_i.repeats << " : " << std::endl;
		TR << "    pool size = " << stage_i.pool << " : " << std::endl;
		TR << "    smoothing = " << stage_i.smoothing << " : " << std::endl;
		TR << "    rms threshold = " << stage_i.rmsthreshold << " : " << std::endl;
		TR << "    P(mutation) = " << stage_i.pmut << " : " << std::endl;
		TR << "    minimizer steps = " << stage_i.maxiter << " : " << std::endl;
		TR << "    pack cycles scaling = " << stage_i.packcycles << " : " << std::endl;
		TR << "    fa_rep ramping schedule : ";
		for ( core::Real r : stage_i.ramp_schedule ) {
			TR << r << " ";
		}
		TR << std::endl;
		TR << "=====================================================================================" << std::endl;

		// set up scorefunction for this stage
		bool changed = scorefxn_->set_smoothing( stage_i.smoothing );
		scorefxn_->set_maxiter_minimize( stage_i.maxiter );
		scorefxn_->set_packer_cycles( stage_i.packcycles );

		if ( i==1 || changed ) {
			initialize_rotamer_set_and_scores( genes[1] );
			optimize_generation( genes, stage_i.ramp_schedule ); // optimize with new smoothing factor
			show_status( genes, "rescored pool at start of stage "+utility::to_string( i ) );

			std::chrono::duration<double> pack_time, min_time;
			scorefxn_->report_and_reset_timers( pack_time, min_time );
			if ( ! basic::options::option[ basic::options::OptionKeys::run::no_prof_info_in_silentout ]() ) {
				TR << "Stage " << i << " pre-refine pack/min time = "
					<< pack_time.count() << "s / " << min_time.count() << "s" << std::endl;
			}
		}

		for ( core::Size j = 1; j <= stage_i.repeats; ++j ) {
			update_tags( genes );
			next_generation( genes, genes_new, stage_i.pool, stage_i.pmut );
			optimize_generation( genes_new, stage_i.ramp_schedule );
			show_status( genes_new, "generated in stage "+std::to_string(i)+" iter "+std::to_string(j) );

			update_pool( genes, genes_new, stage_i.pool, stage_i.rmsthreshold );
			show_status( genes, "pool after stage "+std::to_string(i)+" iter "+std::to_string(j) );

			std::chrono::duration<double> pack_time, min_time;
			scorefxn_->report_and_reset_timers( pack_time, min_time );
			if ( ! basic::options::option[ basic::options::OptionKeys::run::no_prof_info_in_silentout ]() ) {
				TR << "Stage " << i << " iter " << j << " pack/min time = "
					<< pack_time.count() << "s / " << min_time.count() << "s" << std::endl;
			}
		}
	}

	max_rot_cumulative_prob_ = 0.99;
	rot_energy_cutoff_ = 1000;
	altcrossover_= false;
}


void
GAOptimizer::show_status(
	LigandConformers & genes,
	std::string comment,
	bool calculate_native_rmsd,
	bool verbose
) {
	using namespace ObjexxFCL::format;
	std::string statuslog("\nReporting score for "+comment+"\n");

	//utility::vector1< core::Real > rmsds( genes.size(), 0.0 );
	core::Real Emin( genes[1].score() );

	for ( core::Size i = 1; i <= genes.size(); ++i ) {
		LigandConformer &gene = genes[i];

		statuslog += "I/score: " + I(3,int(i)) +" " + F(8,3,gene.score());
		std::pair< core::Real, core::Real > internal_d;
		if ( nativegene_.defined() && calculate_native_rmsd ) {

			core::Real nativermsd = distance_slow( nativegene_, gene );
			internal_d = distance_internal( nativegene_, gene );
			statuslog += " " + F(8,3,nativermsd);
			genes[i].rms( nativermsd );
		} else {
			genes[i].rms( 0.0 );
			internal_d = distance_internal( genes[1], gene );
		}

		statuslog += " " + F(8,3,internal_d.first) + " " + F(8,3,internal_d.second);

		if ( genes[i].score() < Emin ) Emin = genes[i].score();

		statuslog += " "+genes[i].generation_tag()+"\n";
		//statuslog += " "+genes[i].to_string()+"\n";
	}
	TR << statuslog << std::endl;

	if ( verbose ) {
		core::scoring::ScoreFunctionOP sfxn = scorefxn_->get_sfxn();

		TR << "DETAIL: \nSCORE:   score  ";
		sfxn->show_line_headers( TR );
		TR << "    rmsd interfaceSC description\n";

		for ( core::Size i = 1; i <= genes.size(); ++i ) {
			// skip if too high energy; hard-coded cutoff
			if ( genes[i].score() - Emin > 20.0 ) continue;

			TR << "SCORE: ";
			core::pose::PoseOP pose( new core::pose::Pose() );
			genes[i].to_pose( pose );
			sfxn->score( *pose );
			sfxn->show_line( TR, *pose );

			std::string gene_index( comment+"_"+I(4,4,i) );
			TR << " " << F(8,3,genes[i].rms()) << " " << F(8,3,0.0) //F(8,3,score_holo-score_apo)
				<< " " << gene_index << "\n";
		}
		TR << std::endl;
	}
}

void
GAOptimizer::optimize_generation( LigandConformers & genes, utility::vector1<core::Real> const &ramping ) {
	for ( core::Size i = 1; i <= genes.size(); ++i ) {
		LigandConformer gene_prv( genes[i] );
		core::Real score_premin = scorefxn_->score( genes[i] );
		core::Real score( 0.0 );

		score = scorefxn_->optimize( genes[i], ramping, rotamer_data_, rotamer_energies_ );
		TR.Debug << "Score before/after min: " << score_premin << " -> " << score
			<< " (wrep=" << scorefxn_->get_w_rep() << ")" << std::endl;

		genes[i].score( score );
	}
}

void
GAOptimizer::update_tags( LigandConformers &genes ) const {
	core::Size npool = genes.size();
	for ( core::Size ii = 1; ii <= npool; ++ii ) {
		genes[ii].set_generation_tag("parent."+std::to_string(ii));
	}
}

void
GAOptimizer::next_generation(
	LigandConformers const & genes,
	LigandConformers & genes_new,
	core::Size npoolout,
	core::Real pmut
) {
	genes_new.clear();
	core::Size npoolin = genes.size();
	std::string move("");
	std::string tag(""), l_sampling("");

	for ( core::Size iparent = 1; iparent <= npoolout; ++iparent ) {
		iparent = (iparent-1)%npoolin+1;
		LigandConformer newgene( genes[iparent] );

		if ( ( numeric::random::rg().uniform() <= pmut ) || (npoolin == 1) ) {
			newgene = mutate( newgene );
			move = "mutate";
			tag = "mut."+std::to_string(iparent)+" ["+newgene.generation_tag()+"]";
		} else {
			core::Size ipartner = iparent;
			while ( ipartner == iparent ) {
				ipartner = numeric::random::rg().random_range(1, npoolin);
			}

			if ( altcrossover_ ) {
				newgene = crossover_ft( newgene, genes[ipartner] );
			} else {
				newgene = crossover( newgene, genes[ipartner] );
			}
			move = "crossover";
			tag = "cross."+std::to_string(iparent)+"."+std::to_string(ipartner)+" ["+newgene.generation_tag()+"]";
		}

		newgene.set_generation_tag( tag );
		genes_new.push_back( newgene );
		l_sampling += " "+tag;
	}

	TR << "Sampling summary: " << l_sampling << std::endl;
}

void
GAOptimizer::update_pool(
	LigandConformers & genes,
	LigandConformers & genes_new,
	core::Size npool,
	core::Real rmscut
) {
	runtime_assert( genes.size() + genes_new.size() > npool );

	std::sort(genes_new.begin(), genes_new.end(),
		[&](LigandConformer const &lig_i, LigandConformer const &lig_j){ return lig_i.score() < lig_j.score(); } );

	// concatenate everything and sort again
	for ( core::Size ii = 1; ii <= genes.size(); ++ii ) {
		genes_new.push_back( genes[ii]);
	}
	genes.clear();

	std::sort(genes_new.begin(), genes_new.end(),
		[&](LigandConformer const &lig_i, LigandConformer const &lig_j){ return lig_i.score() < lig_j.score(); } );

	utility::vector1< bool > selected( genes_new.size(), false );

	for ( core::Size ii = 1; ii <= genes_new.size(); ++ii ) {
		LigandConformer &gene_new = genes_new[ii];

		bool is_similar( false );
		for ( core::Size jj = 1; jj <= genes.size() && !is_similar; ++jj ) {
			core::Real d = distance_fast( gene_new, genes[jj] );
			is_similar = ( d < rmscut );
		}
		if ( !is_similar ) {
			genes.push_back( gene_new );
			selected[ii] = true;
		}

		if ( genes.size() == npool ) break;
	}

	// fill in missing (in case similar ones filtered out too many times)
	// this shouldn't happen often...
	for ( core::Size ii = 1; ii <= genes_new.size() && genes.size()<npool; ++ii ) {
		if ( !selected[ii] ) {
			genes.push_back( genes_new[ii] );
			selected[ii] = true;
		}
	}
}


// initializes the rotamer sets and precomputes 1b and 2b energies.
void
GAOptimizer::initialize_rotamer_set_and_scores(
	LigandConformer lig
) {
	using namespace core::pack::dunbrack;
	using namespace ObjexxFCL::format;
	using namespace core::chemical;

	core::chemical::ResidueTypeSetCOP restype_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	utility::vector1< core::Size > const &movingscs = lig.moving_scs();
	core::Size nSCs = movingscs.size();

	core::pose::PoseOP pose( new core::pose::Pose() );
	lig.to_pose( pose );

	// if this is called a second time (e.g., due to scfxn changing) just replace existing
	bool first_call = (rotamer_data_.size() == 0);
	rotamer_data_.clear();
	rotamer_energies_.clear();

	rotamer_data_.resize( nSCs );

	core::Real rotprob_epsilon = 0.0;

	core::Size totalrotamers = 0;
	for ( core::Size isc = 1; isc <= movingscs.size(); ++isc ) {
		core::Size const resid( movingscs[isc] );

		// design would modify this logic
		utility::vector1< std::string > allowed_types;
		std::string resname = pose->residue( resid ).name();
		allowed_types.push_back( resname );
		if ( resname == "HIS" ) {
			allowed_types.push_back( "HIS_D" );
		} else if ( resname == "HIS_D" ) {
			allowed_types.push_back( "HIS" );
		}

		PlaceableRotamers res_rotamers;
		PlaceableRotamer res_rotamer_i;

		for ( core::Size iallowed = 1; iallowed <= allowed_types.size(); ++iallowed ) {
			// NOTE! this idealizes the input residue!
			core::conformation::Residue res_i ( restype_set->name_map( allowed_types[iallowed] ) , true );
			pose->replace_residue( resid, res_i, true );

			core::conformation::Residue const &rsd = pose->residue( resid );
			core::chemical::AA aa( rsd.aa() );
			core::Size nchi = rsd.nchi();

			// get rotlib
			core::pack::rotamers::SingleResidueRotamerLibraryCOP residue_rotamer_library
				= core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() );
			SingleResidueDunbrackLibraryCOP
				dun_rotlib( utility::pointer::dynamic_pointer_cast< SingleResidueDunbrackLibrary const >( residue_rotamer_library ) );
			utility::fixedsizearray1< core::Real, 5 > bbs; bbs[1] = rsd.mainchain_torsion(1); bbs[2] = rsd.mainchain_torsion(2);
			utility::vector1< DunbrackRotamerSampleData > rotamers = dun_rotlib->get_all_rotamer_samples( bbs );

			// pack info into rotamer data
			res_rotamer_i.resid = resid;
			res_rotamer_i.restype = pose->residue( resid ).type_ptr();
			res_rotamer_i.chis.resize( nchi , 0.0 );

			// generate all rotamers
			for ( core::Size irot = 1; irot <= rotamers.size(); ++irot ) {
				DunbrackRotamerSampleData const &rotamer( rotamers[irot] );

				// opt-H for Ser/Thr/Tyr
				core::Size nProtChiSteps = 1;
				core::Real protChiOffset = 0.0;
				if ( aa==aa_cys || aa==aa_ser || aa==aa_thr ) {
					nProtChiSteps = 3;
					protChiOffset = 60.0;
				} else if ( aa==aa_tyr ) {
					nProtChiSteps = 2;
				}

				core::Real prob_i = rotamer.probability() / nProtChiSteps;
				prob_i = (1.0-rotprob_epsilon) * prob_i + rotprob_epsilon / rotamers.size();
				prob_i /= allowed_types.size();

				res_rotamer_i.prob = prob_i;
				res_rotamer_i.rotno = irot;
				for ( core::Size ichi = 1; ichi <= rotamer.nchi(); ++ichi ) res_rotamer_i.chis[ichi] = rotamer.chi_mean()[ichi];

				for ( core::Size iprotchi = 1; iprotchi <= nProtChiSteps; ++iprotchi ) {
					if ( nProtChiSteps!=1 ) res_rotamer_i.chis[ nchi ] = protChiOffset + (iprotchi-1) * 360.0/nProtChiSteps;

					// build residue in pose, add lkb info, and score against background
					for ( core::Size ichi = 1; ichi <= nchi; ++ichi ) {
						pose->set_chi( ichi, resid, res_rotamer_i.chis[ichi] );
					}
					res_rotamer_i.lkbrinfo =
						core::scoring::lkball::LKB_ResidueInfoOP(new core::scoring::lkball::LKB_ResidueInfo (pose->residue( resid )) );
					res_rotamer_i.score = scorefxn_->get_1b_energy(pose->residue( resid ), res_rotamer_i.lkbrinfo);

					res_rotamers.push_back( res_rotamer_i );
				}
			}
		}

		// filter by: a) cumulative Prob and b) score
		std::sort(res_rotamers.begin(), res_rotamers.end(),
			[&](PlaceableRotamer const &ii, PlaceableRotamer const &jj){ return ii.prob > jj.prob; } );

		PlaceableRotamers res_rotamers_temp = res_rotamers;
		res_rotamers.clear();
		float probsum( 0.0 ), keptprobsum( 0.0 ), maxkeptprob( 1e-6 );
		for ( core::Size irot = 1; irot <= res_rotamers_temp.size(); ++irot ) {
			if ( probsum > max_rot_cumulative_prob_ ) break; // note: this might split proton chi samples

			bool scoreGood = (res_rotamers_temp[irot].score.score(1.0) <= rot_energy_cutoff_);
			bool rotInBox = (res_rotamers_temp[irot].score.penalty_wtd_ <= 1e-6);

			if ( scoreGood && rotInBox ) {
				keptprobsum += res_rotamers_temp[irot].prob;
				maxkeptprob = std::max( res_rotamers_temp[irot].prob, maxkeptprob );
				res_rotamers.push_back( res_rotamers_temp[irot] );
			}
			probsum += res_rotamers_temp[irot].prob;
		}

		// [1] add input sidechain - with max prob
		res_rotamer_i.prob = maxkeptprob; // probability
		res_rotamer_i.rotno = 0; // 0==input
		res_rotamer_i.restype = lig.get_protein_restype( isc );
		res_rotamer_i.chis = lig.get_protein_chis( isc );

		// fd this could be more efficient...
		if ( pose->residue(resid).type().name() != res_rotamer_i.restype->name() ) {
			core::conformation::Residue res_i( res_rotamer_i.restype , true );
			pose->replace_residue( resid, res_i, true );
		}
		for ( core::Size ichi = 1; ichi <= res_rotamer_i.chis.size(); ++ichi ) {
			pose->set_chi( ichi, resid, res_rotamer_i.chis[ichi] );
		}
		res_rotamer_i.lkbrinfo =
			core::scoring::lkball::LKB_ResidueInfoOP(new core::scoring::lkball::LKB_ResidueInfo (pose->residue( resid )) );
		res_rotamer_i.score = scorefxn_->get_1b_energy(pose->residue( resid ), res_rotamer_i.lkbrinfo);
		res_rotamer_i.score.bonus_wtd_ = -favor_native_; // apply bonus to alternate tautomer

		keptprobsum += res_rotamer_i.prob;
		res_rotamers.push_back( res_rotamer_i );

		// special case for HIS: add alternate protonation of input conformation
		if ( resname == "HIS" || resname == "HIS_D" ) {
			std::string alttype = (resname == "HIS") ? "HIS_D" : "HIS";

			if ( pose->residue(resid).type().name() != alttype ) {
				core::conformation::Residue res_i( restype_set->name_map( alttype ) , true );
				pose->replace_residue( resid, res_i, true );
			}
			for ( core::Size ichi = 1; ichi <= res_rotamer_i.chis.size(); ++ichi ) {
				pose->set_chi( ichi, resid, res_rotamer_i.chis[ichi] );
			}
			res_rotamer_i.lkbrinfo =
				core::scoring::lkball::LKB_ResidueInfoOP(new core::scoring::lkball::LKB_ResidueInfo (pose->residue( resid )) );
			res_rotamer_i.score = scorefxn_->get_1b_energy(pose->residue( resid ), res_rotamer_i.lkbrinfo);
			res_rotamer_i.score.bonus_wtd_ = -favor_native_; // apply bonus to alternate tautomer

			res_rotamer_i.restype = pose->residue(resid).type_ptr();

			keptprobsum += res_rotamer_i.prob;
			res_rotamers.push_back( res_rotamer_i );
		}

		// sort by res type
		std::sort(res_rotamers.begin(), res_rotamers.end(),
			[&](PlaceableRotamer const &ii, PlaceableRotamer const &jj){
				return ((ii.restype->name() > jj.restype->name()) || (ii.restype->name() == jj.restype->name() && ii.prob > jj.prob));
			} );
		// [2] add placeholder sidechain (for "include current")  - with max prob ... KEEP THIS LAST!
		res_rotamer_i.prob = maxkeptprob; // probability
		res_rotamer_i.rotno = 0;
		res_rotamer_i.score.reset();
		res_rotamer_i.restype = nullptr;
		keptprobsum += res_rotamer_i.prob;
		res_rotamers.push_back( res_rotamer_i );

		// finally normalize & accum probability and sotr to rotamer_data_
		core::Real probAccum = 0.0;
		for ( core::Size irot = 1; irot <= res_rotamers.size(); ++irot ) {
			res_rotamers[irot].prob /= keptprobsum;
			probAccum += res_rotamers[irot].prob;
			res_rotamers[irot].prob_accum = probAccum;
		}

		rotamer_data_[isc] = res_rotamers;
		totalrotamers += res_rotamers.size();


		// reporting
		if ( first_call ) {
			TR << "-------------------------------------------------------------------" << std::endl;
			TR << "| RESIDUE   i irot prob   aprob  score  ";
			for ( core::Size ichi = 1; ichi <= res_rotamer_i.chis.size(); ++ichi ) {
				TR << " chi"+std::to_string(ichi)+"  ";
			}
			TR << std::endl;
			for ( core::Size irot = 1; irot <= res_rotamers.size(); ++irot ) {
				//if (irot == 4 && res_rotamers.size() > 6) {
				// TR << "| (" << res_rotamers.size()-5 << " rotamers omitted)" << std::endl;
				//}
				//if (irot > 3 && irot <res_rotamers.size()-2) continue;
				std::string restype_name = (res_rotamers[irot].restype == nullptr) ? "   " : res_rotamers[irot].restype->name();
				TR << "| " << I(3,int(resid)) << " " << restype_name
					<< " " << I(3,int(irot)) << " " << I(3,int(res_rotamers[irot].rotno))
					<< " " << F(6,4,res_rotamers[irot].prob)
					<< " " << F(6,4,res_rotamers[irot].prob_accum)
					<< " " << F(6,2,res_rotamers[irot].score.score(1.0));
				for ( core::Size ichi = 1; ichi <= res_rotamers[irot].chis.size(); ++ichi ) {
					TR << " " << F(6,1,res_rotamers[irot].chis[ichi]);
				}
				TR << std::endl;
			}
			TR << "-------------------------------------------------------------------" << std::endl;
		}
	} // isc

	TR << "Built " << totalrotamers << " rotamers with E<" << rot_energy_cutoff_
		<< " and cumulative P<" << max_rot_cumulative_prob_ << std::endl;

	// [3] compute pairwise energies and construct rotamer interaction graph
	for ( core::Size isc = 1; isc <= movingscs.size(); ++isc ) {
		rotamer_energies_.add_residue( rotamer_data_[isc].size() );
	}
	rotamer_energies_.finalize();

	for ( core::Size isc = 1; isc <= movingscs.size(); ++isc ) {
		core::Size resid_i = movingscs[isc];

		for ( core::Size irot = 1; irot < rotamer_data_[isc].size(); ++irot ) {  // last residue is placeholder (to be scored later)

			if ( pose->residue_type(resid_i).name() != rotamer_data_[isc][irot].restype->name() ) {
				core::conformation::Residue res_i ( rotamer_data_[isc][irot].restype , true );
				pose->replace_residue( resid_i, res_i, true );
			}
			for ( core::Size ichi = 1; ichi <= pose->residue_type( resid_i ).nchi(); ++ichi ) {
				pose->set_chi( ichi, resid_i, rotamer_data_[isc][irot].chis[ichi] );
			}

			rotamer_energies_.energy1b( isc, irot ) = rotamer_data_[isc][irot].score; // includes native "bonus"

			for ( core::Size jsc = isc+1; jsc <= movingscs.size(); ++jsc ) {
				core::Size resid_j = movingscs[jsc];

				for ( core::Size jrot = 1; jrot < rotamer_data_[jsc].size(); ++jrot ) {  // last residue is placeholder (to be scored later)
					if ( pose->residue_type(resid_j).name() != rotamer_data_[jsc][jrot].restype->name() ) {
						core::conformation::Residue res_j ( rotamer_data_[jsc][jrot].restype , true );
						pose->replace_residue( resid_j, res_j, true );
					}
					for ( core::Size jchi = 1; jchi <= pose->residue_type( resid_j ).nchi(); ++jchi ) {
						pose->set_chi( jchi, resid_j, rotamer_data_[jsc][jrot].chis[jchi] );
					}

					rotamer_energies_.energy2b( isc, irot, jsc, jrot ) =
						scorefxn_->get_2b_energy(
						pose->residue( resid_i ), rotamer_data_[isc][irot].lkbrinfo,
						pose->residue( resid_j ), rotamer_data_[jsc][jrot].lkbrinfo);
				}
			}
		}
	}
}

}
}
}
