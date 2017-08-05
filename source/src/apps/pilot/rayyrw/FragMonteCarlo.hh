#include <core/types.hh>

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>

#include <utility/vector1.hh>

#include <numeric/random/random.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>
#include <apps/pilot/rayyrw/util.hh>

#ifndef apps_pilot_rayyrw_FragMonteCarlo_hh
#define apps_pilot_rayyrw_FragMonteCarlo_hh

static THREAD_LOCAL basic::Tracer tr("FragMonteCarlo");

using namespace ObjexxFCL::format;

// per fragment info
struct FragID {
	int mer;
	int pos;
	int picker_rank;
	int sph_rank;
};


// options declaration
class FragMonteCarlo;
class FragMonteCarlo {
public:
	// default constructor
	FragMonteCarlo(
		core::Real wt_dens_in,
		core::Real wt_overlap_in,
		core::Real wt_closab_in,
		core::Real wt_clash_in,
		core::Real null_frag_score_in_
	) :
		n_total_frags_(0),
		n_total_rsds_(0),
		wt_dens_(wt_dens_in),
		wt_overlap_(wt_overlap_in),
		wt_closab_(wt_closab_in),
		wt_clash_(wt_clash_in),
		null_frag_score_(null_frag_score_in_) {}


	// 1. load scorefiles into scores_1b_ and scores_2b_
	// 2. assign a null_frag to frag_cands at each position
	// 3. initialize assigned_frags_ with all null_frags
	void load_scorefiles( std::string fragidx_file,
		std::string scorefile_dens,
		std::string scorefile_overlap,
		std::string scorefile_nonoverlap );


	void run( bool verbose=false,
		core::Real sa_start_temp=1000.0,
		core::Real sa_end_temp=1.0,
		core::Size sa_nsteps=25,
		core::Size mc_nsteps=1500 );

	void report_score( core::Real temp=0.0 );

	//void report_results( std::stringstream &outlines,
	//                     int model_ctr=0,
	//                     int runid=0 );
	std::string report_results( int model_ctr=0,
		int runid=0 );

	void updated_covered_rsds( utility::vector1<int> &covered_rsds, FragID const fragidx );

	bool isOverlapping( int i_fragidx, int j_fragidx );

	// with definition
	inline void check_file_exists( std::string file_name ){
		if ( ! utility::file::file_exists( file_name ) ) { utility_exit_with_message( "Unable to open file: " + file_name + '\n' ); }
	}

	inline bool isNullFrag( int fragidx ){ return ( fragidx == 0 ); }

	inline std::string fragid_to_fragfn( FragID fragid ){
		std::string fragfn = "after_rotation_frags." + int2str(fragid.mer) + "." + int2str(fragid.pos) + "." + int2str(fragid.picker_rank) + "." + int2str(fragid.sph_rank) + ".????.pdb";
		return fragfn; }

	inline utility::vector1<int> const & get_candidate_frags( core::Size pos ){ return pos_to_fragcands_[ pos ]; }

	utility::vector1<int> get_assigned_frags(){ return assigned_frags_; }

	void initialize_frag_assignment( core::Size init_type );


private:
	/////////////////////////// variables  ///////////////////////////
	core::Size n_total_frags_;
	core::Size n_total_rsds_;
	core::Real wt_dens_, wt_overlap_, wt_closab_, wt_clash_, null_frag_score_;

	// information storage
	std::map< int, FragID > fragidx_to_fragid_;

	// given one position, what are all the candidates (fragidx)
	std::map< int, utility::vector1<int> > pos_to_fragcands_; // frag_candidates as fragidx

	// the main container to store results of mc sampling
	utility::vector1<int> assigned_frags_; // as fragidx

	utility::vector1< utility::vector1<core::Real> > scores_2b_;
	utility::vector1<core::Real> scores_1b_;
	utility::vector1<core::Real> rmsd_table_;

	/////////////////////////// functions ///////////////////////////
	core::Real cal_frag_score( int candidate_fragidx );
};



void
FragMonteCarlo::
load_scorefiles(
	std::string fragidx_file,
	std::string scorefile_dens,
	std::string scorefile_overlap,
	std::string scorefile_nonoverlap
){
	////////////////////////////////////////////////////////////////////////////////
	// read in frag index file
	// 1. store fragidx_to_fragid_
	// 2. store pos_to_fragcands_
	// 3. get n_total_frags_
	// 4. get n_total_rsds_
	//
	n_total_frags_ = 1;
	n_total_rsds_  = 0;

	check_file_exists( fragidx_file );
	std::ifstream in_fragidx( fragidx_file.c_str() );
	while ( ! in_fragidx.eof() ) {
		FragID fragid;
		int fragidx;
		in_fragidx >> fragidx >> fragid.mer >> fragid.pos >> fragid.picker_rank >> fragid.sph_rank;
		fragidx_to_fragid_[ fragidx ] = fragid;
		pos_to_fragcands_[ fragid.pos ].push_back( fragidx );
		n_total_frags_++;
		n_total_rsds_ = std::max( (int) n_total_rsds_, fragid.pos );
	}

	////////////////////////////////////////////////////////////////////////////////
	// 1. assigned null_frag at each position
	// 2. add null_frag in candidate fragments for each position
	//
	for ( core::Size i=1; i<=n_total_rsds_; ++i ) {
		assigned_frags_.push_back( 0 );  // initialized by assignning null at each position
		pos_to_fragcands_[i].push_back( 0 );
	}

	////////////////////////////////////////////////////////////////////////////////
	// parse score files, and store scores in scores_1b and scores_2b
	//
	scores_1b_.resize( n_total_frags_, 0.0 );
	rmsd_table_.resize( n_total_frags_, 0.0 );
	scores_2b_.resize( n_total_frags_, utility::vector1<core::Real>( n_total_frags_, 0.0 ) );

	// density scores
	check_file_exists( scorefile_dens );
	std::ifstream in_density( scorefile_dens.c_str() );
	tr.Info << "load " << scorefile_dens << " ..." << std::endl;
	while ( ! in_density.eof() ) {
		int i;
		core::Real dens_score, rmsd;
		in_density >> i >> dens_score >> rmsd;
		//tr.Debug << i << " dens_score: " << dens_score << " rmsd: " << rmsd << std::endl;
		scores_1b_[i] = wt_dens_*dens_score;
		rmsd_table_[i] = rmsd;
	}

	// overlap scores
	check_file_exists( scorefile_overlap );
	std::ifstream in_overlap( scorefile_overlap.c_str() );
	tr.Info << "load " << scorefile_overlap << " ..." << std::endl;
	while ( ! in_overlap.eof() ) {
		int i,j;
		core::Real overlap_score;
		in_overlap >> i >> j >> overlap_score;
		if ( isOverlapping(i,j) ) {
			overlap_score = wt_overlap_*overlap_score;
			//tr.Debug << i << " " << fragidx_to_fragid_[i].pos << " " << j << " " << fragidx_to_fragid_[j].pos <<  " " << " overlap_score: " << overlap_score << " rmsd: " << rmsd_table_[i] << " " << rmsd_table_[j] << std::endl;
			scores_2b_[i][j] = overlap_score;
			scores_2b_[j][i] = overlap_score;
		}
	}

	// fill non-overlap scores
	check_file_exists( scorefile_nonoverlap );
	std::ifstream in_nonoverlap( scorefile_nonoverlap.c_str() );
	tr.Info << "load " << scorefile_nonoverlap << " ..." << std::endl;
	while ( ! in_nonoverlap.eof() ) {
		int i,j;
		core::Real closab_score, clash_score;
		core::Real nonoverlap_score=0.0;
		in_nonoverlap >> i >> j >> closab_score >> clash_score;
		if ( ! isOverlapping(i,j) ) {
			nonoverlap_score = wt_closab_*closab_score + wt_clash_*clash_score;
			//tr.Debug << i << " " << fragidx_to_fragid_[i].pos << " " << j << " " << fragidx_to_fragid_[j].pos << " closab_score: " << closab_score << " clash_score: " << clash_score << " nonoverlap_score: " << nonoverlap_score << " rmsd: " << rmsd_table_[i] << " " << rmsd_table_[j] << std::endl;
			scores_2b_[i][j] = nonoverlap_score;
			scores_2b_[j][i] = nonoverlap_score;
		}
	}
}



void
FragMonteCarlo::
initialize_frag_assignment(
	core::Size init_type /*1: null; 2: random; 3: lowrmsd*/
){
	//assign random frag for each pos
	for ( core::Size pos=1; pos<=n_total_rsds_; ++pos ) {
		if ( init_type==1 ) {
			assigned_frags_[pos] = 0;
			tr.Info << "null assigned fragments" << assigned_frags_[pos] << std::endl;

		} else if ( init_type == 2 ) {
			assigned_frags_[pos] = pos_to_fragcands_[pos][ numeric::random::random_range( 1, pos_to_fragcands_[pos].size() ) ];
			tr.Info << "random assigned fragments " << assigned_frags_[pos] << std::endl;

		} else if ( init_type == 3 ) { // primarily used for debugging
			utility::vector1<int> const &frag_cands = get_candidate_frags(pos);
			core::Real best_rmsd=999;

			for ( int j=1, j_end=frag_cands.size(); j<=j_end; ++j ) {
				if ( isNullFrag( frag_cands[j] ) ) {
					continue;
				}

				core::Real rmsd = rmsd_table_[ frag_cands[j] ];
				if ( ( rmsd < best_rmsd ) && ( rmsd <= 3.0 ) ) {
					best_rmsd = rmsd;
					assigned_frags_[pos] = frag_cands[j];
				}
			}

			// Assign null fragment since none of the frag_cands in the position lower than 3.0
			// This step is not necessary in the current algorithms since at the very begining the assigned_frags_ was initialized by null fragment only. One caveat is this function might be call some times other than at the very beginning of simulations.
			if ( best_rmsd == 999 ) {
				assigned_frags_[pos] = 0;
			}

			tr.Info << "lowest_rmsd fragments " << assigned_frags_[pos] << std::endl;

		} else {
			utility_exit_with_message("type has to be 1:null or 2:random or 3:lowrmsd");
		}
	}
}



// Given a candidate_fragidx, return a frag_score of the fragment to the assigned ones,
// 1. null_frag - return null_frag_score
// 2. frag - return frag_score
core::Real
FragMonteCarlo::
cal_frag_score(
	int candidate_fragidx
){
	core::Real thisfrag_score = 0;

	if ( isNullFrag( candidate_fragidx ) ) {
		return null_frag_score_;
	}

	// density_score
	thisfrag_score = scores_1b_[ candidate_fragidx ];
	core::Size cand_pos = fragidx_to_fragid_[candidate_fragidx].pos;
	for ( core::Size pos=1; pos<=n_total_rsds_; ++pos ) {
		int assigned_fragidx = assigned_frags_[pos];
		if ( ( ! isNullFrag( assigned_fragidx ) ) && ( pos != cand_pos ) ) {
			thisfrag_score += scores_2b_[ candidate_fragidx ][ assigned_fragidx ];
		}
	}
	return thisfrag_score;
}


bool
FragMonteCarlo::
isOverlapping(
	int i_fragidx,
	int j_fragidx
){
	FragID i_fragid = fragidx_to_fragid_[ i_fragidx ];
	FragID j_fragid = fragidx_to_fragid_[ j_fragidx ];

	/*assert ( i_fragid.mer == j_fragid.mer );
	if ( std::abs( i_fragid.pos - j_fragid.pos ) <= ( i_fragid.mer - 1 ) ){
	return true;
	} else {
	return false;
	}
	*/
	int i_start = i_fragid.pos;
	int i_end = i_start+i_fragid.mer-1;
	int j_start = j_fragid.pos;
	int j_end = j_start+j_fragid.mer-1;

	if ( i_fragid.mer <= j_fragid.mer ) {
		return ( (i_start>=j_start && i_start<=j_end) || (i_end>=j_start && i_end<=j_end) );
	} else {
		return ( (j_start>=i_start && j_start<=i_end) || (j_end>=i_start && j_end<=i_end) );
	}
}



void
FragMonteCarlo::
report_score(
	core::Real temp /*0.0*/
){
	// score everything
	core::Real score=0.0;
	for ( core::Size i=1; i<=n_total_rsds_; ++i ) {
		if ( isNullFrag( assigned_frags_[i] ) ) {
			tr.Info << "position " << i << ": frag NULL score = " << null_frag_score_ << std::endl;
			score += null_frag_score_;
		} else {
			int thisfrag = assigned_frags_[i];
			core::Real score_i = cal_frag_score( thisfrag );
			FragID this_fragid = fragidx_to_fragid_[thisfrag];
			tr.Info << "position " << i << ": frag " << thisfrag << " score = " << score_i << " rmsd = " << rmsd_table_[thisfrag] << " " << fragid_to_fragfn( this_fragid ) << std::endl;
			score += score_i;
		}
	}
	tr.Info << "temp: " << temp << " total_score: " << score << std::endl;
}



void
FragMonteCarlo::
updated_covered_rsds(
	utility::vector1<int> &covered_rsds,
	FragID const fragid
){
	core::Size range = fragid.mer-1;
	core::Size pos = fragid.pos;

	tr.Debug << "pos: " << int2str(pos);
	for ( core::Size i=pos; i<=pos+range; ++i ) {
		covered_rsds.push_back(i);
		tr.Debug << i;
	}
	tr.Debug << std::endl;

	// remove redundancy
	std::sort( covered_rsds.begin(), covered_rsds.end() );
	covered_rsds.erase( std::unique(covered_rsds.begin(), covered_rsds.end()), covered_rsds.end() );

	tr.Debug << "covered_rsds(" << covered_rsds.size() << "): ";
	for ( core::Size i=1; i<=covered_rsds.size(); ++i ) {
		tr.Debug << int2str(covered_rsds[i]) << " ";
	}
	tr.Debug << std::endl;
}



// this will be the calcualte_state() from
// 1. coverage
// 2.
// 3.
std::string
FragMonteCarlo::
report_results(
	int model_ctr,
	int runid
){
	// score everything
	core::Real score=0.0;
	core::Real rmsd_addup=0.0;
	int n_not_null_frags=0;
	int n_null_frags=0;
	int n_lowrmsd_frags=0;
	int n_highrmsd_frags=0;
	core::Real rmsd_threshold=2.5;

	utility::vector1<int> covered_rsds;
	utility::vector1<int> high_rmsd_pos;

	std::stringstream outlines;
	std::string tag = "run_" + int2str(runid) + "_model_" + int2str(model_ctr) + ":";
	for ( core::Size pos=1; pos<=n_total_rsds_; ++pos ) {
		outlines << LJ( 20, tag );
		if ( isNullFrag( assigned_frags_[pos] ) ) {

			outlines << " rsn: "   << I(5, pos)
				<< ", rmsd: "  << F(6, 3, 0.0)
				<< ", score: " << F(6, 3, null_frag_score_)
				<< ", null"
				<< std::endl;

			n_null_frags++;
			score += null_frag_score_;

		} else {
			int this_fragidx = assigned_frags_[pos];
			core::Real score_i = cal_frag_score( this_fragidx );
			core::Real rmsd = rmsd_table_[this_fragidx];
			rmsd_addup += rmsd;
			FragID this_fragid = fragidx_to_fragid_[this_fragidx];

			outlines << " rsn: "    << I(5, pos)
				<< ", rmsd: "  << F(6, 3, rmsd)
				<< ", score: " << F(6, 3, score_i)
				<< ", " << fragid_to_fragfn( this_fragid )
				<< std::endl;

			if ( rmsd <= rmsd_threshold ) {
				n_lowrmsd_frags++;
				updated_covered_rsds( covered_rsds, this_fragid );
				tr.Info << "n_lowrmsd_frags: " << n_lowrmsd_frags << " n_not_null_frags: " << n_not_null_frags << " lowrmsd_frag: " << this_fragidx << " rmsd: " << rmsd << std::endl;
			} else {
				n_highrmsd_frags++;
				high_rmsd_pos.push_back(pos);
			}

			n_not_null_frags++;
			score += score_i;
		}
	}
	core::Real correct_rate = ((core::Real)n_lowrmsd_frags/(core::Real)n_not_null_frags)*100;
	//core::Size coverage = covered_rsds.size();

	outlines << LJ( 20, tag )
		<< " SCORE: " << F(5, 4, score)
		<< ", correct_rate: " << F(4, 2, correct_rate)
		<< ", n_not_null_frags: " << I(4, n_not_null_frags)
		<< ", n_lowrmsd_frags: " << I(4, n_lowrmsd_frags)
		<< ", n_null_frags: " << I(4, n_null_frags)
		<< ", n_highrmsd_frags: " << I(4, n_highrmsd_frags)
		<< ", rmsd_addup: " << rmsd_addup
		<< std::endl;

	return outlines.str();
}



void
FragMonteCarlo::
run(
	bool verbose, /*false*/
	core::Real sa_start_temp, /*1000.0*/
	core::Real sa_end_temp, /*1.0*/
	core::Size sa_nsteps, /*25*/
	core::Size mc_nsteps /*1500*/
){
	core::Real temp = sa_start_temp;
	core::Real temp_scale = std::pow( sa_end_temp/sa_start_temp, 1.0/((core::Real)sa_nsteps) );

	for ( core::Size temp_ctr=1; temp_ctr<=sa_nsteps; ++temp_ctr ) {
		//tr.Info << "simualated annealing temperature: " << temp << " temp_counter: " << temp_ctr << std::endl;
		// monte carlo
		core::Size counter=1;
		for ( core::Size step=1; step<=mc_nsteps; ++step ) {
			int this_pos = numeric::random::random_range( 1, n_total_rsds_ );

			// reset frag to null at this position, this way
			assigned_frags_[this_pos] = 0;

			// calculate compatibility scores for all the cadidate placements at the given pos
			utility::vector1<int> const &frag_cands = get_candidate_frags( this_pos );
			core::Size n_frag_cands = frag_cands.size();
			utility::vector1<core::Real> frag_scores( n_frag_cands );
			core::Real best_frag_score = 1e30, prob_sum = 0.0;

			for ( core::Size i=1; i<=n_frag_cands; ++i ) {
				frag_scores[i] = cal_frag_score( frag_cands[i] );
				best_frag_score = std::min( frag_scores[i], best_frag_score );
			}

			for ( core::Size i=1; i<=n_frag_cands; ++i ) {
				frag_scores[i] -= best_frag_score;
				frag_scores[i] = std::exp( -std::min( frag_scores[i]/temp, 100.0 ) );
				prob_sum += frag_scores[i];
			}

			// select fragment normalize prob
			core::Real fragpicker = prob_sum*numeric::random::uniform();
			int picked = 0;
			while ( fragpicker>=0 && picked<int(n_frag_cands) ) {
				picked++;
				fragpicker -= frag_scores[picked];
			}
			assert( fragpicker < 0 );
			assigned_frags_[this_pos] = frag_cands[picked];
			counter ++;
		}
		// tr.Info << "monte carlo steps finished: " << counter << std::endl;
		temp *= temp_scale;

		if ( verbose ) {
			report_score( temp );
		}

	}
	report_score( temp );
}
#endif
