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
/// @author Yifan Song
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_loop_grower_LoopGrower_hh
#define INCLUDED_protocols_loop_grower_LoopGrower_hh

#include <protocols/loop_grower/LoopGrower.fwd.hh>

#include <fstream>
#include <iostream>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/Sequence.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/Jump.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/loop_grower/DensSkeleton.hh>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>
#include <numeric/random/WeightedSampler.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>
#include <core/id/SequenceMapping.hh>

#include <queue>



namespace protocols {
namespace loop_grower {

static basic::Tracer TRACER("protocols.loop_grower.LoopGrower");

// 56 bytes
struct ResTorsions {
	ResTorsions(core::Real phi_in, core::Real psi_in, core::Real omega_in, core::Size nchi_in, core::Real chi1_in, core::Real chi2_in, core::Real chi3_in, core::Real chi4_in) {
		phi_=phi_in;
		psi_=psi_in;
		omega_=omega_in;
		nchi_ = nchi_in;
		chi1_=chi1_in; chi2_=chi2_in; chi3_=chi3_in; chi4_=chi4_in;
	}

	ResTorsions(core::pose::Pose const &pose, int i) {
		phi_=pose.phi(i);
		psi_=pose.psi(i);
		omega_=pose.omega(i);
		nchi_ = pose.residue(i).nchi();
		if ( nchi_ > 0 ) chi1_=pose.chi(1,i);
		if ( nchi_ > 1 ) chi2_=pose.chi(2,i);
		if ( nchi_ > 2 ) chi3_=pose.chi(3,i);
		if ( nchi_ > 3 ) chi4_=pose.chi(4,i);
	}

	core::Real phi_, psi_, omega_, chi1_, chi2_, chi3_, chi4_;
	core::Size nchi_;
};

struct SheetPositions {

	SheetPositions(utility::vector1< ResTorsions > residues, numeric::xyzMatrix< core::Real > rotation, numeric::xyzVector< core::Real > translation, core::Size jumpid, core::Size baseres){
		residues_ = residues;
		rotation_ = rotation;
		translation_ = translation;
		jumpid_ = jumpid;
		baseres_ = baseres;
	}

	SheetPositions(core::pose::Pose const &pose, core::Size rangelow, core::Size rangehi, core::Size jumpid, core::Size baseres){
		for ( core::Size i=rangelow; i<=rangehi; i++ ) {
			residues_.push_back( ResTorsions( pose, i ) );
		}
		rotation_ = pose.jump( jumpid ).get_rotation();
		translation_ = pose.jump( jumpid ).get_translation();
		jumpid_ = jumpid;
		baseres_ = baseres;
	}

	utility::vector1< ResTorsions > residues_;
	numeric::xyzMatrix< core::Real > rotation_;
	numeric::xyzVector< core::Real > translation_;
	core::Size jumpid_, baseres_;
};


class LoopPartialSolution {
public:
	LoopPartialSolution (core::pose::Pose const &pose, int reslow, int reshigh, core::Real score) {
		for ( int i=reslow; i<=reshigh; ++i ) {
			residues_.push_back( ResTorsions( pose, i ) );
			utility::vector1< core::Vector > backbone;
			calphas_.push_back( pose.residue(i).atom("CA").xyz() );
			backbone.push_back( pose.residue(i).atom("CA").xyz() );
			backbone.push_back( pose.residue(i).atom("N" ).xyz() );
			backbone.push_back( pose.residue(i).atom("C").xyz() );
			backbone.push_back( pose.residue(i).atom("O").xyz() );
			if ( pose.residue(i).name3()!="GLY" ) backbone.push_back( pose.residue(i).atom("CB").xyz() );
			backbones_.push_back(backbone);
		}
		set_centerofmass();
		score_ = score;
		bonus_score_ = 0;
	}
	LoopPartialSolution(){}

	void store_coordinates(core::pose::Pose const&pose, int reslow, int reshigh, core::Size lower_fasta){
		core::Size k = lower_fasta;
		for ( int i=reslow; i<=reshigh; ++i ) {
			for ( core::Size j=1; j<=pose.residue_type(i).natoms(); j++ ) {
				core::id::AtomID atomid = core::id::AtomID(j,k);
				ids_.push_back( atomid);
				positions_.push_back( pose.xyz(core::id::AtomID(j,i)));
			}
			k++;
		}
	}

	void set_centerofmass();

	void set_coordinates(core::pose::Pose & pose){
		pose.batch_set_xyz(ids_,positions_);
	}

	void store_sheetpositions(core::pose::Pose & pose, core::Size rangelo, core::Size rangehi, core::Size jumpid, core::Size basres){
		sheets_.push_back( SheetPositions( pose, rangelo, rangehi, jumpid, basres ));
	}

	void write_beam_lps(std::ofstream &lpsfile);

	void write_beam(std::ofstream &outbeam);

	void apply( core::pose::Pose &pose, int range1lo, int range1hi, int range2lo, int range2hi ) {
		int ctr = 0;

		calphas_.clear();
		for ( int i=range1lo; i<=range2hi; ++i ) {
			if ( i < range2lo && i> range1hi ) continue;
			ctr++;
			ResTorsions const &res_i = residues_[ctr];
			pose.set_phi( i, res_i.phi_ );
			pose.set_psi( i, res_i.psi_ );
			pose.set_omega( i, res_i.omega_ );
			if ( !pose.is_centroid() && pose.conformation().residue_typeset_mode( true ) != core::chemical::CENTROID_ROT_t ) {
				if ( res_i.nchi_ > 0 ) pose.set_chi( 1, i, res_i.chi1_ );
				if ( res_i.nchi_ > 1 ) pose.set_chi( 2, i, res_i.chi2_ );
				if ( res_i.nchi_ > 2 ) pose.set_chi( 3, i, res_i.chi3_ );
				if ( res_i.nchi_ > 3 ) pose.set_chi( 4, i, res_i.chi4_ );
			}
		}

		if ( ctr != (int)residues_.size() ) {
			TRACER << "ctr and number of residues " << ctr << " " << residues_.size() << std::endl;
		}
		runtime_assert( ctr == (int)residues_.size() );
	}

	void
	apply_sheets(core::pose::Pose & pose);

	void
	add_sheets(core::pose::Pose & pose, core::Size takeoffres, core::Size totalres);

	void set_rms( core::Real RMS ){ RMS_ = RMS; }
	void set_gdt( core::Real GDT ){ GDT_ = GDT; }

	void set_bonus_score( core::Real bonus ){ bonus_score_ = bonus; }

	core::Real get_bonus_score() const { return bonus_score_; }

	core::Real get_rms() const { return RMS_; }

	core::Real get_gdt() const { return GDT_; }

	void apply( core::pose::Pose &pose, int range1lo, int range1hi ) {
		apply(pose,range1lo,range1hi,range1hi,range1hi);
	}

	void set_id( std::string id ){ id_ = id; }


	std::string
	get_id(){ return id_; }

	core::Real size(){ return residues_.size(); }

	void push_back_restorsions(ResTorsions restorsions){ residues_.push_back(restorsions); }
	void push_back_backbone( utility::vector1< core::Vector > backbone ){ backbones_.push_back(backbone); }
	void push_back_sheet( SheetPositions sheet ){ sheets_.push_back(sheet); }

	void set_positions(utility::vector1< numeric::xyzVector<core::Real> > positions) { positions_ = positions; }
	void set_ids(utility::vector1< core::id::AtomID > ids) { ids_ = ids; }

	void set_score( core::Real score ){ score_ = score; }

	void set_history( utility::vector1< std::pair<core::Size,core::Size> > history, std::pair<core::Size,core::Size> newfrag ){
		fragmenthistory_ = history;
		fragmenthistory_.push_back(newfrag);
	}

	void set_history( utility::vector1< std::pair<core::Size,core::Size> > history){
		fragmenthistory_ = history;
	}

	utility::vector1< std::pair<core::Size,core::Size> > get_history(){ return fragmenthistory_; }

	utility::vector1< core::Vector > get_calphas() const { return calphas_; }
	utility::vector1< utility::vector1< core::Vector > > get_backbone() const { return backbones_; }
	core::Real get_score() const { return score_; }

	core::Real partialrms(const LoopPartialSolution&, int fragmelt, int total_lower, bool takeoffonly, int rmswindow ) const;
	core::Real max_calpha_distance( const LoopPartialSolution&, core::Size fragmelt, core::Size total_lower, bool takeoffonly, bool full_loop, bool lower ) const;
	bool operator < (const LoopPartialSolution& lps1) const{ return score_ < lps1.score_; }

private:
	utility::vector1< ResTorsions > residues_;
	utility::vector1< SheetPositions > sheets_;
	utility::vector1 < core::Vector > calphas_;
	utility::vector1< core::id::AtomID > ids_;
	utility::vector1< numeric::xyzVector<core::Real> > positions_;
	utility::vector1 < utility::vector1< core::Vector > > backbones_;
	utility::vector1 < std::pair<core::Size,core::Size> > fragmenthistory_;
	numeric::xyzVector<core::Real> centerofmass_;
	core::Real score_;
	core::Real bonus_score_ = 0, RMS_ = 0, GDT_ = 0;
	std::string id_;
};

struct RMSComparator {
	inline bool operator()(const LoopPartialSolution& lps1, const LoopPartialSolution& lps2)
	{
		return (lps1.get_rms() < lps2.get_rms());
	}
};


class LoopPartialSolutionStore {
public:
	LoopPartialSolutionStore() = default;

	LoopPartialSolutionStore( core::Size N, core::Real rmscutoff, core::Size master_beam_width, core::Real master_beam_cutoff ) { maxelts_=N; rmscutoff_ = rmscutoff; master_beam_width_ = master_beam_width;
		master_beam_cutoff_ = master_beam_cutoff; }


	void
	setfilterparams( core::Size fragmelt, core::Size rmswindow, core::Size parallelcount, core::Real beamscorecutoff, bool dump_errors, bool dump_beam, bool writebeams, bool clustercheck, bool fafilter,
		bool samplesheets, bool filterprevious, bool checksymm, bool asymmdump, bool dumpfinalbeam ){
		fragmelt_ = fragmelt;
		rmswindow_ = rmswindow;
		parallelcount_ = parallelcount;
		beamscorecutoff_ = beamscorecutoff;
		dump_errors_ = dump_errors;
		dump_beam_ = dump_beam;
		writebeams_ = writebeams;
		clustercheck_ = clustercheck;
		fafilter_ = fafilter;
		samplesheets_ = samplesheets;
		filterprevious_ = filterprevious;
		checksymm_ = checksymm;
		asymmdump_ = asymmdump;
		dumpfinalbeam_ = dumpfinalbeam;
	}

	void
	store(LoopPartialSolution const &new_entry){
		solutions_.push_back(new_entry);
	}
	void set_filterprevious( bool filterprevious ){ filterprevious_ = filterprevious; }

	void
	filter( core::pose::Pose& pose, core::Size maxbeam, core::Size total_lower, core::Size lower_res, core::Size upper_res, core::Size cycle);

	void
	diversityfilter(core::Size maxbeamsize, core::Size total_lower);

	void
	parametercheck_filter(core::pose::Pose& pose, core::Size fragmelt, core::Size total_lower, core::Size lower_res, core::Size upper_res, core::Size rmswindow, bool dump_beam,
		core::Size& totalreqbeam, core::Size& totalreqmaxbeam, core::Real& totalreqscorecut, bool writebeams, core::Size cycle, core::Size parallelcount );

	void
	zscoretransform();

	void
	cluster_check(LoopPartialSolution nearnativelps, core::pose::Pose& pose, core::Size total_lower, core::Size lower_res, core::Size upper_res );

	void
	push( core::pose::Pose& pose, LoopPartialSolution const &new_entry, core::Size fragmelt, core::Real beamscorecutoff, core::Size total_lower, bool lower, bool dump_errors, core::Real RMS, core::Size lower_res, core::Size upper_res );

	void
	report_rms_and_scores(core::Size rangelo, core::Size rangehi);

	void
	sort(){ std::sort(solutions_.begin(), solutions_.end()); }

	void
	sortrms(){ std::sort(solutions_.begin(), solutions_.end(), RMSComparator()); }

	bool
	filterprevious( LoopPartialSolution lps, core::Size total_lower, core::Size fragmelt, core::Real rmscutoff, core::Size masterbeamwidth);

	void
	skeleton_filter(core::pose::Pose& pose, DensSkeleton& skeleton, core::Size start_res, core::Size stop_res, core::Size lower_term, core::Size res_gap);

	//LoopPartialSolution
	//pop () {
	// LoopPartialSolution lps = solutions_.pop();
	// return lps;
	//}

	LoopPartialSolution
	operator [] ( core::Size solutions_index ) {
		LoopPartialSolution lps = solutions_[solutions_index];
		return lps;
	}


	void
	clear() {
		solutions_.clear();
	}

	void set_fastas(core::Size lower, core::Size upper){
		lower_fasta_ = lower;
		upper_fasta_ = upper;
	}
	void set_poses(core::Size lower, core::Size upper){
		lower_pose_ = lower;
		upper_pose_ = upper;
	}
	void set_cutpoint(core::Size cutpoint){
		cutpoint_ = cutpoint;
	}

	core::Size get_lower_fasta(){ return lower_fasta_;}
	core::Size get_upper_fasta(){ return upper_fasta_;}
	core::Size get_lower_pose(){ return lower_pose_;}
	core::Size get_upper_pose(){ return upper_pose_;}
	core::Size get_cutpoint(){ return cutpoint_;}

	utility::vector1<LoopPartialSolution> get_solutions(){ return solutions_; }
	utility::vector1<LoopPartialSolution> get_filteronly_solutions(){ return filteronly_solutions_; }
	void store_filteronly_solutions( utility::vector1<LoopPartialSolution> filteronly){ filteronly_solutions_ = filteronly; }

	bool
	is_empty() {
		return (solutions_.empty());
	}
	core::Size size() {
		return solutions_.size();
	}


	void writelpsstore(core::Size loopid, core::Size jobid){
		std::ofstream lpsfile;
		std::string filename = "lpsfile_"+utility::to_string(loopid)+"."+utility::to_string(jobid)+".txt";
		lpsfile.open( filename.c_str() );
		lpsfile << lower_fasta_ << " " << upper_fasta_ << " " << lower_pose_ << " " << upper_pose_ << " " << cutpoint_ << std::endl;
		for ( core::Size i=1; i<=solutions_.size(); i++ ) {
			LoopPartialSolution beam = solutions_[i];
			beam.write_beam_lps(lpsfile);
		}
		lpsfile.close();
	}

private:
	core::Real rmscutoff_ = 0.0; // or 1.2 ??
	core::Real master_beam_cutoff_ = 3.0;
	core::Real beamscorecutoff_ = 0.85;
	core::Size maxelts_ = 1; // ??
	core::Size master_beam_width_ = 1;
	core::Size upper_fasta_ = 0; // ??
	core::Size lower_fasta_ = 0; // ??
	core::Size lower_pose_ = 0; // ??
	core::Size upper_pose_ = 0; // ??
	core::Size cutpoint_ = 0; // ??
	core::Size fragmelt_ = 1;
	core::Size rmswindow_ = 1;
	core::Size parallelcount_ = 0;
	core::Size rank_ = 0; // ??
	utility::vector1<LoopPartialSolution> solutions_, old_solutions_, filteronly_solutions_;
	bool dump_errors_ = false;
	bool dump_beam_ = false;
	bool writebeams_ = false;
	bool clustercheck_ = false;
	bool fafilter_ = true;
	bool famin_ = false;
	bool samplesheets_ = true;
	bool filterprevious_ = false;
	bool checksymm_ = false;
	bool asymmdump_ = true;
	bool dumpfinalbeam_ = false;
};

class LoopGrower: public protocols::moves::Mover {
public:
	LoopGrower() = default; // Default values for member variables set at declaration.

	LoopGrower(core::sequence::SequenceOP seq, protocols::loops::Loop loop, int resstart, int resstop,  utility::vector1<core::fragment::FragSetOP> const &fragments ) :
		seq_(seq), loop_(loop), resstart_(resstart), resstop_(resstop), fragments_(fragments) {

		// Should these be the general defaults?
		pack_min_cycles_ = 2;
		fafilter_pmcycles_ = 2;
		rmscutoff_  = 1.2;

		sf_ = core::scoring::get_score_function();
		cen_sf_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth");
		cenrot_sf_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_cenrot_relax");
		update_fragment_library_pointers();
		native_ = nullptr;
	}

	// run the protocol
	void apply(core::pose::Pose & pose);

	// one growcycle direction
	core::Real single_grow( core::pose::Pose& pose, core::pose::Pose& pose_cen, LoopPartialSolutionStore& solutionset, const core::chemical::ResidueTypeCOPs& restypes_pose,
		const core::chemical::ResidueTypeCOPs& restypes_pos_cen, core::Size lower_pose, core::Size upper_pose, core::Size upper_term, int lower_fasta, int upper_fasta, core::Size torsionrangelo,
		core::Size torsionrangehi, bool insert_lower, core::Size initial_melt_left, bool is_cterm, bool is_nterm, core::Size cycle, int n_to_insert );

	void update_to_stored( core::pose::Pose& growpose, core::pose::Pose& growpose_cen, const core::chemical::ResidueTypeCOPs& restypes_pose,
		const core::chemical::ResidueTypeCOPs& restypes_pose_cen, int & lower_pose, int & upper_pose, int & lower_fasta, int & upper_fasta,
		core::Size newreslow, core::Size newreshi, bool is_cterm, bool is_nterm);

	void update_and_writelps(LoopPartialSolutionStore & solutionset, core::pose::Pose & fa_pose, core::pose::Pose & pose_cen, core::chemical::ResidueTypeCOPs & restypes_pose,
		core::chemical::ResidueTypeCOPs & restypes_pose_cen, int lower_pose, int upper_pose, bool is_nterm, bool is_cterm, core::Size fasta_range_low, core::Size fasta_range_hi, core::Size pose_range_low,
		core::Size pose_range_hi, int torsionrangelo, int torsionrangehi, int tgt_jump, bool update_pose);

	//get the n to n+3 hbonds so we can penalize overwound helices
	core::Real nton3_hbond_score(core::pose::Pose& pose);

	//Get the hbond energy of bonds formed by the residues including and between lower and upper
	core::Real get_resrange_hbond_energy(core::pose::Pose& pose, core::Size lower, core::Size upper);

	core::Real modifiedscore(core::pose::Pose& fapose, core::pose::Pose& cen_pose, core::Size rangelo, core::Size rangehi);
	core::Real modifieddensity(core::pose::Pose& pose, core::Size rangelo, core::Size rangehi, core::Real density_weight, core::Size & includesheets);
	core::Real sheetscore(core::pose::Pose& fapose, core::pose::Pose& cen_pose, core::Size rangelo, core::Size rangehi);

	//Returns a boolean regarding whether the phi and psi meet beta rama space.
	bool is_beta(core::Real phi, core::Real psi);

	//Rescore all the structures in the solution set
	void rescoresolutionset( LoopPartialSolutionStore& solutionset, core::pose::Pose& fa_pose, core::pose::Pose& cen_pose, core::Size torsionrangelo, core::Size torsionrangehi);

	//Check the 5 residues, N or C term of the cutpoint, for all structures in the solutionset against the stored coordinates. If at least one has 4 matches remove all those that don't.
	void coordinate_filter(LoopPartialSolutionStore& solutionset, core::pose::Pose pose, bool lower, core::Size lower_fasta, core::Size upper_fasta, core::Size lower_pose, core::Size rangehi);

	void addnativesolution(LoopPartialSolutionStore& solutionset, core::pose::Pose& fa_pose, core::pose::Pose& cen_pose, int natlowerstart, int natlowerstop, int natupperstart, int natupperstop,
		int lower_pose, core::Size torsionrangelo, core::Size torsionrangehi);

	void
	store_sheets(core::pose::Pose & pose, LoopPartialSolution & lps);

	// one "step" of refinement
	// cen_min + repack + fa_min
	void refine_cycle( core::pose::Pose & pose, core::pose::Pose & cen_pose, int loop_start, int loop_end, bool finalrefinement, int cycle, int beam, int fragment, core::Size is_lower);

	void full_atom_beam( LoopPartialSolutionStore& solutionset, core::pose::Pose & fa_pose, core::pose::Pose & cen_pose, core::Size lower_pos, core::Size upper_pos);

	void set_scorefunction ( core::scoring::ScoreFunctionOP sf ) {
		sf_ = sf;
		//core::scoring::ScoreFunctionOP sf_soft_ = sf_->clone();
		//sf_soft_->set_weight( core::scoring::fa_rep, 0.02 );
	}
	void set_cen_scorefunction ( core::scoring::ScoreFunctionOP sf ) { cen_sf_ = sf; }
	void set_cenrot_scorefunction( core::scoring::ScoreFunctionOP sf) { cenrot_sf_ = sf; }
	void set_beamwidth ( core::Size beamwidth ) { beamwidth_ = beamwidth; }
	void set_fragtrials ( core::Size fragtrials ) { fragtrials_ = fragtrials; }
	void set_chainbreak ( core::Real chainbreak ) { chainbreak_ = chainbreak; }
	void set_continuous_weight ( core::Real continuous_weight ) { continuous_weight_ = continuous_weight; }
	void set_rmscutoff ( core::Real rmscutoff ) { rmscutoff_ = rmscutoff; }
	void set_master_beam_cutoff( core::Real master_beam_cutoff ) { master_beam_cutoff_ = master_beam_cutoff; }
	void set_sheetbonus( core::Real sheetbonus ) { sheetbonus_ = sheetbonus; }
	void set_sheet_tolerance( core::Real sheet_tolerance ) { sheet_tolerance_ = sheet_tolerance; }
	void set_sc_scale( core::Real sc_scale ) { sc_scale_ = sc_scale; }
	void set_windowdensweight( core::Real windowdensweight ){ windowdensweight_ = windowdensweight; }
	void set_master_beam_width( core::Size master_beam_width ) { master_beam_width_ = master_beam_width; }
	void set_rmswindow( core::Size rmswindow ) { rmswindow_ = rmswindow; }
	void set_steps( core::Size steps ) { steps_ = steps; }
	void set_pcount( core::Size parallelcount ) { parallelcount_ = parallelcount; }
	void set_sheetcriteria( core::Size sheetcriteria ) { sheetcriteria_ = sheetcriteria; }
	void set_loopnumber( core::Size loopnumber ) { loopnumber_ = loopnumber; }
	void set_beamscorecutoff ( core::Real beamscorecutoff ) { beamscorecutoff_ = beamscorecutoff; }
	void set_debug ( bool debug ) { debug_ = debug; }
	void set_fragmelt( core::Size fragmelt ) { fragmelt_ = fragmelt; }
	void set_minmelt( core::Size minmelt ) { minmelt_ = minmelt; }
	void set_pack_min_cycles( core::Size pack_min_cycles ) { pack_min_cycles_ = pack_min_cycles; }
	void set_dumpbeam ( bool dumpbeam ) {dumpbeam_ = dumpbeam; }
	void set_dumpfinalbeam ( bool dumpfinalbeam ) {dumpfinalbeam_ = dumpfinalbeam; }
	void set_dumprms( bool dumprms ) {dumprms_ = dumprms; }
	void set_dumperrors( bool dumperrors ) { dumperrors_ = dumperrors;}
	void set_direction ( core::Size direction ) { direction_ = direction; }
	void set_minimize ( bool minimize ) { minimize_ = minimize; }
	void set_native ( core::pose::PoseOP native ) { native_ = native; }
	void set_nativegrow ( bool nativegrow ) { nativegrow_ = nativegrow; }
	void set_greedy( bool greedy ){ greedy_ = greedy; }
	void set_parametercheck( bool parametercheck ) { parametercheck_ = parametercheck; }
	void set_cenrot( bool cenrot ){ cenrot_ = cenrot; }
	void set_writebeams( bool writebeams ) { writebeams_ = writebeams; }
	void set_readbeams( bool readbeams ) { readbeams_ = readbeams; }
	void set_clustercheck( bool clustercheck ) { clustercheck_ = clustercheck; }
	void set_rescorebeams( bool rescorebeams ) { rescorebeams_ = rescorebeams; }
	void set_writelps( bool writelps ) { writelps_ = writelps; }
	void set_fafilter( bool fafilter ) { fafilter_ = fafilter; }
	void set_famin( bool famin ) { famin_ = famin; }
	void set_samplesheets( bool samplesheets ) { samplesheets_ = samplesheets; }
	void set_trackfragments( bool trackfragments ) { trackfragments_ = trackfragments; }
	void set_filterprevious( bool filterprevious ) { filterprevious_ = filterprevious; }
	void set_rephasemap( bool rephasemap ) { rephasemap_ = rephasemap; }
	void set_checksymm( bool checksymm ) { checksymm_ = checksymm; }
	void set_continuous_sheets( bool continuous_sheets ) { continuous_sheets_ = continuous_sheets; }
	void set_auto_stop( bool auto_stop ) { auto_stop_ = auto_stop; }
	void set_storedbeams( std::string storedbeams ) { storedbeams_ = storedbeams; }
	void set_filterbeams( std::string filterbeams ) { filterbeams_ = filterbeams; }
	void set_coordfile( std::string coordfile ) { coordfile_ = coordfile; }
	void set_skeletonfile( std::string skeletonfile ) { skeleton_file_ = skeletonfile; }

	std::string get_name() const { return "LoopGrower"; }

	void update_fragment_library_pointers();

	void add_fragment_csts( core::pose::Pose &pose, core::Size startfasta, core::Size endfasta, core::Size natstoplow, core::Size natstarthi, core::Size startlower );

	void add_user_csts( core::pose::Pose &pose );

	//write everything in the solutionset to the disk. Used in parallelization
	void write_to_disk( LoopPartialSolutionStore solutionset, core::Size step, core::Size added_lower, core::Size added_upper, bool filteronly, bool lower );

	void read_from_disk( LoopPartialSolutionStore & solutionset, int & cycle, bool & lower, bool filterbeams);

	void read_coordfile();

	void fafilter( LoopPartialSolutionStore &solutionset, core::pose::Pose &fapose, core::pose::Pose &cenpose, core::Size total_lower, core::Size torsionrangelo, core::Size torsionrangehi,
		core::Size cycle, core::Size lower_fasta, core::Size upper_fasta, core::Size lower_pose);

	core::Real GDThatonative(core::pose::Pose const &pose, int natlow, int nathi, int natstoplow, int natstarthi, int startlow, int stoplow, int starthi, int stophi);
	core::Real RMStonative(core::pose::Pose const &pose, int natlow, int nathi, int natstoplow, int natstarthi, int startlow, int stoplow, int starthi, int stophi);
	void atomGDTha(core::Vector atom1, core::Vector atom2, core::Size &GDTha);
	core::Size check_coordinates(core::pose::Pose& pose, core::Size lower_pose, core::Size upper_pose, core::Size lower_fasta, core::Real radius);

	//check each residues electron density in the range residue range to determine if there is a dropoff indicating the loop grower has run out of density and should stop
	bool
	check_auto_stop(core::pose::Pose & pose, core::Size lower_res, core::Size upper_res);


	LoopPartialSolutionStore getSolutionSet(){ return solutionset_;}

	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;

private:
	core::Size beamwidth_ = 1;
	core::Size fragtrials_ = 10;
	core::Size fragmelt_ = 1;
	core::Size minmelt_ = 1;
	core::Size lowest_ranked_native_ = 0; // arbitrary init
	core::Size storelow_ = 0;
	core::Size storehi_ = 0;
	core::scoring::ScoreFunctionOP sf_, cen_sf_, cenrot_sf_; //sf_soft,

	bool debug_ = false;
	core::sequence::SequenceOP seq_;
	protocols::loops::Loop loop_;
	int resstart_ = 0;
	int resstop_ = 0;
	int upper_fasta_ = 0; // arbitrary init
	int lower_fasta_ = 0; // arbitrary init
	int rmsrangelo_ = 0; // arbitrary init
	int rmsrangehi_ = 0; // arbitrary init
	core::Real chainbreak_ = 1.5;
	core::Real continuous_weight_ = 0.3;
	core::Real rmscutoff_ = 0.0; // or 1.2 ??
	core::Real beamscorecutoff_ = 0.85;
	core::Real master_beam_cutoff_ = 3.0;
	core::Real windowdensweight_ = 30;
	core::Real startingscore_ = 0.0;
	core::Real censtartingscore_ = 0.0;
	core::Real sheetbonus_ = 0.5;
	core::Real sheet_tolerance_ = 0.7;
	core::Real sc_scale_ = 1;
	utility::vector1<core::fragment::FragSetOP> fragments_;
	utility::vector1<boost::unordered_map<core::Size, core::fragment::Frame> > library_;
	core::Size pack_min_cycles_ = 1; // or 2??
	core::Size fafilter_pmcycles_ = 1; // or 2??
	core::Size direction_ = 0;
	core::Size master_beam_width_ = 1;
	core::Size rmswindow_ = 1;
	core::Size steps_ = 1e4;
	core::Size parallelcount_ = 0;
	core::Size sheetcriteria_ = 2;
	core::Size loopnumber_ = 0;
	core::Size sheetsize_ = 4;
	core::Size total_residues_ = 0;
	core::Size insert_pose_ = 0; // arbitrary init
	core::Size maxfrag_ = 0; // arbitrary init
	core::Size numjumps_ = 0;
	bool dumpbeam_ = false;
	bool dumpfinalbeam_ = false;
	bool dumprms_ = false;
	bool dumperrors_ = false;
	bool minimize_ = true;
	bool nativegrow_ = false;
	bool greedy_ = true;
	bool parametercheck_ = false;
	bool cenrot_ = false;
	bool writebeams_ = false;
	bool readbeams_ = false;
	bool clustercheck_ = false;
	bool rescorebeams_ = false;
	bool writelps_ = false;
	bool fafilter_ = true;
	bool famin_ = false;
	bool samplesheets_ = true;
	bool trackfragments_ = false;
	bool cenrotfilter_ = false;
	bool filterprevious_ = false;
	bool rephasemap_ = false;
	bool checksymm_ = false;
	bool continuous_sheets_ = true;
	bool auto_stop_ = false;
	bool asymmdump_ = true;
	std::string storedbeams_, filterbeams_, coordfile_, skeleton_file_;
	std::map< std::pair< core::Size, core::Size >, utility::vector1< numeric::xyzVector< core::Real > > > scoringcoords_;

	core::pose::PoseOP native_;
	core::pose::Pose fa_seqpose_;
	core::pose::Pose cen_seqpose_;
	DensSkeleton skeleton_;

	LoopPartialSolutionStore solutionset_;

};
} // hybridize
} // protocols

#endif
