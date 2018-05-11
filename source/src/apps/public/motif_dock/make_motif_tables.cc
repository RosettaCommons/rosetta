// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// headers

#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/motif/util.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/mh.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/motif/reference_frames.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <devel/init.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/xyzVector.hh>
#include <protocols/idealize/IdealizeMover.hh>
// #include <protocols/sicdock/Assay.hh>
#include <protocols/sic_dock/SICFast.hh>
#include <protocols/sic_dock/util.hh>
#include <protocols/sic_dock/read_biounit.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/constraint_movers/AddConstraintsToCurrentConformationMover.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/file/file_sys_util.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>
#include <numeric/HomogeneousTransform.hh>

#include <boost/foreach.hpp>


static basic::Tracer TR( "motif_hash_util" );

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;
typedef numeric::xyzTransform<core::Real> Xform;
typedef utility::vector1<core::id::AtomID> AIDs;
using std::make_pair;
using core::chemical::AA;
using numeric::HomogeneousTransform;
using core::id::AtomID;
using basic::options::option;
using namespace basic::options::OptionKeys;
using core::pose::Pose;
using core::Real;
using core::scoring::ScoreFunctionOP;
using core::Size;
using numeric::max;
using numeric::min;
using numeric::random::gaussian;
using numeric::random::uniform;
using numeric::rotation_matrix_degrees;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using namespace ObjexxFCL::format;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using utility::io::izstream;
using utility::io::ozstream;
using utility::file_basename;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_file;
using numeric::geometry::hashing::Real3;
using numeric::geometry::hashing::Real6;
using core::pose::xyzStripeHashPoseCOP;

using core::scoring::motif::ResPairMotifs;
using core::scoring::motif::load_motifs;
using core::scoring::motif::get_nbrs;
using core::scoring::motif::get_residue_pair_rt6;
using core::scoring::motif::get_sasa;
using core::scoring::motif::MotifHash;
using core::scoring::motif::ResPairMotif;
using core::scoring::motif::ResPairMotifMetaBinner;
using core::scoring::motif::ResPairMotifsMap;
using core::scoring::motif::tag_from_pdb_fname;
using core::scoring::motif::write_motifs_binary;
using core::scoring::motif::XformScore;
using core::scoring::motif::XformScoreMap;
using core::scoring::motif::SC_SC;
using core::scoring::motif::SC_BB;
using core::scoring::motif::SC_PH;
using core::scoring::motif::SC_PO;
using core::scoring::motif::BB_BB;
using core::scoring::motif::BB_PH;
using core::scoring::motif::BB_PO;
using core::scoring::motif::PH_PO;
using core::scoring::motif::RPM_Type_NONE;
using core::scoring::motif::MOTIF_HASH_CART_SIZE;
using protocols::sic_dock::KMGT;


inline Size sicdock_max_threads(){
	return 1;
}
inline Size sicdock_num_threads(){
	return 1;
}
inline Size sicdock_thread_num(){
	return 1;
}


inline void print_motifs(std::ostream & out){
	using namespace basic::options::OptionKeys;
	ResPairMotifs motifs;
	load_motifs( option[mh::print_motifs](), motifs );
	for ( ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i ) {
		out << *i << endl;
	}
}

inline void merge_motifs(){
	using namespace basic::options::OptionKeys;
	if ( ! basic::options::option[basic::options::OptionKeys::mh::motif_out_file].user() ) {
		utility_exit_with_message("must sepcify name for merged file -motif_out_file");
	}
	vector1<string> const & fnames ( option[mh::merge_motifs]() );
	ResPairMotifs motifs; load_motifs( fnames, motifs );

	ResPairMotifMetaBinner binner;
	ResPairMotifsMap mmap;

	for ( Size i = 1; i <= motifs.size(); ++i ) {
		ResPairMotif const & m(motifs[i]);
		ResPairMotifMetaBinner::Key k = binner.motif_bin_hash(m);
		mmap[k].push_back(m);
	}

	TR << "write to disk" << endl;
	BOOST_FOREACH ( ResPairMotifsMap::value_type const & p,mmap ) {
		std::string label = binner.hash_to_labal(p.first);
		write_motifs_binary(option[mh::motif_out_file]()+(label.size()?"_":"")+label,p.second);
	}
	TR << "merge_motifs done" << endl;

}

template<class T> T sqr(T const & x){ return x*x; }

inline void harvest_scores(){
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::motif;
	using core::scoring::motif::aa_trustworthiness;

	TR << "harvest_scores: read motifs" << endl;
	vector1<string> const & fnames( option[mh::harvest_scores]() );
	ResPairMotifs motifs; load_motifs( fnames, motifs );

	RPM_Type type = core::scoring::motif::RPM_Type_NONE;
	if ( option[mh::filter::motif_type].user() && option[mh::filter::motif_type]()=="BB_BB" ) type = BB_BB;

	TR << "computing counts with requested bin sizes" << endl;
	ResPairMotifMetaBinner binner;
	utility::vector1<XformScoreMap> xscoremap(sicdock_max_threads());

	bool agg_with_max = option[mh::harvest::agg_with_max]();
	Real const HASH_CART_RESL     = option[mh::harvest::hash_cart_resl]();
	Real const HASH_ANGL_RESL     = option[mh::harvest::hash_angle_resl]();
	Real const SMOOTHING_FAC      = option[mh::harvest::smoothing_factor]();
	Real const EXP_NORM           = 1.0 / sqrt(3.14159) / SMOOTHING_FAC / HASH_CART_RESL;

	Real const global_multiplier = option[mh::harvest::multiplier]();

	cout << "HASH_CART_RESL" << HASH_CART_RESL << endl;
	cout << "HASH_ANGL_RESL" << HASH_ANGL_RESL << endl;
	cout << "SMOOTHING_FAC" << SMOOTHING_FAC << endl;
	cout << "EXP_NORM" << EXP_NORM << endl;

	Real const maxbins = std::numeric_limits<XformScore::Key>::max();

	for ( double d2 = 0.0; d2 < 3.0; d2 += 0.5 ) {
		cout << "WT " << d2 << " " << exp(-d2/SMOOTHING_FAC/SMOOTHING_FAC) * EXP_NORM << endl;
	}


	for ( int i = 1; i <= (int)motifs.size(); ++i ) {
		if ( i%(motifs.size()/10)==0 ) cout << (float)i/(float)motifs.size()*100 << "%" << endl;
		ResPairMotif const & m(motifs[i]);
		if ( type==RPM_Type_NONE ) type=m.type();
		if ( m.type()!=type ) utility_exit_with_message("harvest_scores: mismatched types");
		Real6 const rt6 = m.rt();

		ResPairMotifMetaBinner::Key key1 = binner.motif_bin_hash(m);
		assert( binner.hash_to_labal(key1) == binner.motif_bin_label(m) );
		// cout << binner.hash_to_labal(key1) << " " << binner.motif_bin_label(m) << endl;

		// find correct score to contribute to, create if necessary
		XformScoreMap::iterator xscoreiter = xscoremap[sicdock_thread_num()].find(key1);
		if ( xscoreiter==xscoremap[sicdock_thread_num()].end() ) {
			xscoremap[sicdock_thread_num()][key1] = core::scoring::motif::XformScoreOP(new XformScore(option[mh::harvest::hash_cart_resl](),option[mh::harvest::hash_angle_resl]()));
			xscoreiter = xscoremap[sicdock_thread_num()].find(key1);
			Real numbins = 1.0;
			for ( int i = 1; i <= 6; ++i ) numbins *= (Real)xscoreiter->second->hasher_.dimsizes()[i];
			if ( numbins > maxbins ) {
				cout << "NUM BINS " << F(20,0,numbins) << endl;
				cout << "MAX BINS " << F(20,0,maxbins) << endl;
				utility_exit_with_message("too many bins!");
			}
		}
		XformScore & xsfx( *xscoreiter->second );


		// get weighted score
		float motif_score = global_multiplier;
		{
			// XformScore::Key k = xsfx.bin_index(rt6);
			if ( basic::options::option[basic::options::OptionKeys::mh::harvest::weight_by_energy]() ) {
				motif_score *= std::exp(-m.score() * aa_trustworthiness(m.aa1())* aa_trustworthiness(m.aa2()) );
				// cout << motif_score <<" "<< m.score() <<" "<< m.aa1()<<" "<< aa_trustworthiness(m.aa1()) <<" "<<m.aa2() <<" "<< aa_trustworthiness(m.aa2()) << endl;
			}
		}
		// static int count = 0; if(++count > 30) utility_exit_with_message("arstoine");

		using numeric::geometry::hashing::Bin6D;
		// Bin6D const half = hasher_.halfbin6(rt6);
		Bin6D bin6 = xsfx.hasher_.bin6(rt6);
		XformScore::Key key = xsfx.hasher_.bin_index(bin6);


		// if( agg_with_max ){

		xsfx.max_score(key,-m.score());

		// } else { // loop over neighbors and smooth
		// utility_exit_with_message("this code must be checked");
		utility::fixedsizearray1<Size,6> beg,end;
		for ( numeric::Size i = 1; i <= 3; ++i ) {
			beg[i] = max(                 0             ,(int)bin6[i]-1);
			end[i] = min((int)xsfx.hasher_.dimsizes()[i]-1,(int)bin6[i]+1)+1;
		}
		// 4 & 5 wrap for euler angles
		for ( numeric::Size i = 4; i <= 5; ++i ) {
			beg[i] = bin6[i] > 0                          ? bin6[i]-1 : xsfx.hasher_.dimsizes()[i]-1 ;
			end[i] = bin6[i] < xsfx.hasher_.dimsizes()[i]-1 ? bin6[i]+1 : 0;
			end[i] = end[i] < xsfx.hasher_.dimsizes()[i]-1 ? end[i]+1 : 0;
		}
		beg[6] = max(                 0             ,(int)bin6[6]-1);
		end[6] = min((int)xsfx.hasher_.dimsizes()[6]-1,(int)bin6[6]+1)+1;
		//  if( half[i] && (i==4||i==5||end[i] < xsfx.hasher_.dimsizes()[i]-1)) ++end[i];
		//  if(!half[i] && (i==4||i==5||beg[i] >                0       )) --beg[i];
		// }
		// cout <<I(5,xsfx.hasher_.dimsizes()[4])<<I(5,xsfx.hasher_.dimsizes()[5])<<I(5,xsfx.hasher_.dimsizes()[6])<<I(5,xsfx.hasher_.dimsizes()[1])<<I(5,xsfx.hasher_.dimsizes()[2])<<I(5,xsfx.hasher_.dimsizes()[3])<<endl;
		// cout <<I(5, beg[4])<<I(5, beg[5])<<I(5, beg[6])<<I(5, beg[1])<<I(5, beg[2])<<I(5, beg[3])<<endl;
		// cout <<I(5,bin6[4])<<I(5,bin6[5])<<I(5,bin6[6])<<I(5,bin6[1])<<I(5,bin6[2])<<I(5,bin6[3])<<endl;
		// cout <<I(5, end[4])<<I(5, end[5])<<I(5, end[6])<<I(5, end[1])<<I(5, end[2])<<I(5, end[3])<<endl;
		// cout << endl;


		utility::fixedsizearray1<Size,6> ob6;
		// int count_this_rt=0;
		for ( ob6[1] = beg[1]; ob6[1] != end[1]; ++ob6[1] ) {
			for ( ob6[2] = beg[2]; ob6[2] != end[2]; ++ob6[2] ) {
				for ( ob6[3] = beg[3]; ob6[3] != end[3]; ++ob6[3] ) {
					for ( ob6[4] = beg[4]; ob6[4] != end[4]; ob6[4] = (ob6[4]==xsfx.hasher_.dimsizes()[4]-1) ? 0 : ob6[4]+1 ) {
						for ( ob6[5] = beg[5]; ob6[5] != end[5]; ob6[5] = (ob6[5]==xsfx.hasher_.dimsizes()[5]-1) ? 0 : ob6[5]+1 ) {
							for ( ob6[6] = beg[6]; ob6[6] != end[6]; ++ob6[6] ) {
								// cout <<"ITER "<<I(5,ob6[4])<<I(5,ob6[5])<<I(5,ob6[6])<<I(5,ob6[1])<<I(5,ob6[2])<<I(5,ob6[3])<<endl;
								XformScore::Key other_key = xsfx.hasher_.bin_index(ob6);


								// // get target cells
								// // stupid thing about the SixDofHasher, must do lookups concentrically... shell 5 is the max... speed up by lowering it
								// vector1<XformScore::Key> surrounding_cells;
								// for(Size lookuprad = 0; lookuprad <= HASH_LOOKUP_RADIUS; ++lookuprad){
								//  std::vector<uint64_t> shell_tmp = xsfx.hasher_.radial_bin_index(lookuprad,rt6);
								//  // cout << "RADIAL_BIN_LOOKUP " << lookuprad << " " << shell_tmp.size() << endl;

								//   // loop over all cells close to this motif
								//   BOOST_FOREACH(XformScore::Key const & other_key, shell_tmp ){
								//   // runtime_assert_msg(testkeys.find(other_key)==testkeys.end(),"key already seen!");
								//  // testkeys.insert(other_key);

								Real6 const other_rt6 = xsfx.hasher_.bin_center_point(xsfx.hasher_.bin_from_index(other_key));
								float rt6_dis2 = 0;
								for ( int ii=1; ii<=3; ++ii ) rt6_dis2 += sqr( (rt6[ii]-other_rt6[ii])/HASH_CART_RESL );
								for ( int ii=4; ii<=6; ++ii ) rt6_dis2 += sqr( (rt6[ii]-other_rt6[ii])/HASH_ANGL_RESL );
								Real const wt = exp(-rt6_dis2/SMOOTHING_FAC/SMOOTHING_FAC) * EXP_NORM;
								if ( wt < 0.2 ) continue;
								// cout << "WT " << lookuprad << " " << rt6_dis2 << " " << wt << endl;
								Real const sc = wt*motif_score;
								if ( sc < option[mh::harvest::min_bin_val]()/30.0 ) continue;
								if ( agg_with_max ) xsfx.max_score(other_key,-m.score());
								else             xsfx.add_score(other_key,wt*motif_score);
								// }

							}
						}
					}
				}
			}
		}
		// }
	}


	XformScoreMap & xscoremap_joint(xscoremap.front());

	if ( option[mh::harvest::min_bin_val].user() ) {
		TR << "prune" << endl;
		BOOST_FOREACH ( XformScoreMap::value_type& v,xscoremap_joint ) v.second->prune_small_bins(option[mh::harvest::min_bin_val]());
	}

	TR << xscoremap_joint;

	TR << "write to disk" << endl;
	string PATH = option[mh::motif_out_file]();
	utility::file::create_directory_recursive(PATH);
	BOOST_FOREACH ( XformScoreMap::value_type const & p,xscoremap_joint ) {
		std::string label = binner.hash_to_labal(p.first);
		p.second->write_binary(PATH+"/"+PATH+(label.size()?"_":"")+label+".xh.bin.gz");
	}
	if ( option[mh::harvest::dump]() ) write_motifs_binary(PATH+"/"+PATH+".rpm.bin.gz",motifs);
	TR << PATH << endl;
	TR << "harvest_scores done" << endl;
}

inline void print_scores(){
	using namespace basic::options::OptionKeys;
	//XformScore *xh = NULL;
	core::scoring::motif::XformScoreOP xh( new XformScore(option[mh::harvest::hash_cart_resl](), option[mh::harvest::hash_angle_resl]()));
	XformScore::read_binary(xh,option[mh::print_scores]());
	TR << *xh << endl;
	xh->print_scores(cout);
}

inline bool harvest_motifs_one(
	ResPairMotifs & motifs,
	string const & tag0,
	Pose const & pose,
	core::scoring::dssp::Dssp & dssp,
	//core::scoring::etable::Etable const & etable,
	Size const & ir,
	core::scoring::EnergyEdge const & edge,
	vector1<Real> const & /*occupancy*/,
	vector1<Real> const & bfactors,
	vector1<Real> const & /*rsd_sasa*/,
	vector1<Real> const & nbrs,
	int & nbad, int &/* nbadbfac*/, int & nbadocc, int & /*nbadiface*/, int & /*nbaddist*/
){
	using namespace basic::options::OptionKeys;
	using core::scoring::hbonds::HBond;
	using namespace core::pose::motif;
	using core::scoring::motif::RPM_Type;

	Size const jr( edge.get_second_node_ind() );

	if ( !pose.residue(ir).is_protein() ) { ++nbad; ++nbadocc ; return false; }
	if ( !pose.residue(jr).is_protein() ) { ++nbad; ++nbadocc ; return false; }

	bool const is_disulf = core::conformation::is_disulfide_bond(pose.conformation(),ir,jr);
	bool const is_sspaired = dssp.in_paired_strands(ir,jr);

	Real const fa_atr = edge[core::scoring::fa_atr];
	Real const fa_rep = edge[core::scoring::fa_rep];
	Real const fa_sol = edge[core::scoring::fa_sol];
	//Real const fa_intra_atr_xover4 = edge[core::scoring::fa_intra_atr_xover4];
	//Real const fa_intra_rep_xover4 = edge[core::scoring::fa_intra_rep_xover4];
	Real const fa_intra_sol_xover4 = edge[core::scoring::fa_intra_sol_xover4];
	Real const fa_intra_rep = edge[core::scoring::fa_intra_rep];
	Real const fa_elec = edge[core::scoring::fa_elec];
	//Real const fa_intra_elec = edge[core::scoring::fa_intra_elec];

	//Real const lk_ball = edge[core::scoring::lk_ball];
	Real const lk_ball_wtd = edge[core::scoring::lk_ball_wtd];
	//Real const lk_ball_iso = edge[core::scoring::lk_ball_iso];
	//Real const lk_ball_bridge = edge[core::scoring::lk_ball_bridge];
	//Real const lk_ball_bridge_uncpl = edge[core::scoring::lk_ball_bridge_uncpl];

	Real const pro_close = edge[core::scoring::pro_close];

	Real const hb_sr_bb = edge[core::scoring::hbond_sr_bb];
	Real const hb_lr_bb = edge[core::scoring::hbond_lr_bb];
	Real const hb_bb_sc = edge[core::scoring::hbond_bb_sc];
	Real const hb_sc = edge[core::scoring::hbond_sc];

	//core::scoring::EnergyMap const & i_map = energies.residue_total_energies(ir);
	//core::scoring::EnergyMap const & j_map = energies.residue_total_energies(jr);
	//energies.show();

	//Real dslf_fa13 = 0.0;
	//if(is_disulf){
	// dslf_fa13 = i_map[core::scoring::dslf_fa13] + j_map[core::scoring::dslf_fa13];
	//}
	//Real const rama = i_map[core::scoring::rama] + j_map[core::scoring::rama];
	//Real const rama_prepro = i_map[core::scoring::rama_prepro] + j_map[core::scoring::rama_prepro];
	//Real const omega = i_map[core::scoring::omega] + j_map[core::scoring::omega];
	//Real const p_aa_pp = i_map[core::scoring::p_aa_pp] + j_map[core::scoring::p_aa_pp];
	//Real const fa_dun = i_map[core::scoring::fa_dun] + j_map[core::scoring::fa_dun];
	//Real const fa_dun_rot = i_map[core::scoring::fa_dun_rot] + j_map[core::scoring::fa_dun_rot];
	//Real const fa_dun_dev = i_map[core::scoring::fa_dun_dev] + j_map[core::scoring::fa_dun_dev];
	//Real const fa_dun_semi = i_map[core::scoring::fa_dun_semi] + j_map[core::scoring::fa_dun_semi];
	//Real const yhh_planarity = i_map[core::scoring::yhh_planarity] + j_map[core::scoring::yhh_planarity];
	//Real const hxl_tors = i_map[core::scoring::hxl_tors] + j_map[core::scoring::hxl_tors];
	//Real const ref = i_map[core::scoring::ref] + j_map[core::scoring::ref];

	string tag = tag0;
	if ( is_disulf ) tag = tag+"dsf";


	Real const nb1 = nbrs[ir];
	Real const nb2 = nbrs[jr];
	Real const & bfir(bfactors[ir]);
	Real const & bfjr(bfactors[jr]);


	Real const full_score = 1.0*fa_atr + 0.55*fa_rep + 1.0*fa_sol + 1.0*fa_intra_sol_xover4 + 0.005*fa_intra_rep + 1.0*lk_ball_wtd + 1.0*fa_elec + 1.25*pro_close + 1.0*hb_sr_bb + 1.0*hb_lr_bb + 1.0*hb_bb_sc + 1.0*hb_sc;



	//////////////////////////////////// actual motifs ////////////////////////////////////////////
	Real E_TH = -0.001;

	Xform const bbr1 = core::pose::motif::get_backbone_reference_frame(pose,ir);
	Xform const bbr2 = core::pose::motif::get_backbone_reference_frame(pose,jr);
	if ( bbr1.bad() || bbr2.bad() ) {
		cout << "bad bbr " << ir << " " << jr << endl;
		return false;
	}

	static string DBGDIR = option[out::file::o]()+"/";
	utility::file::create_directory_recursive(DBGDIR);

	static bool header=true;


	////////////////////////////////// 'standard' motifs + sc frame motifs ///////////////////////////////


	ResPairMotif bbbb( tag, pose, (~bbr1*bbr2).rt6(), ir, jr, nb1, nb2, full_score, 0, 0, 0, 0, 0, bfir, bfjr, is_sspaired, BB_BB );
	if ( is_disulf || ( full_score < E_TH ) ) {
		motifs.push_back(bbbb);
		if ( option[mh::harvest::dump]() ) {
			if ( header ) { bbbb.print_header(cout); header=false; }
			cout << bbbb << endl;
			Real rms = bbbb.dump_aligned_motif(DBGDIR+bbbb.tag()+".pdb",pose,ir,jr);
			if ( rms > 0.01 ) utility_exit_with_message("bad rms "+string_of(rms)+" "+DBGDIR+bbbb.tag()+".pdb");
		}
	}


	return true;
}

inline void harvest_motifs(){
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace ObjexxFCL::format;

	core::chemical::AtomTypeSetCAP atypeset = core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
	core::scoring::etable::EtableOptions eopt;
	core::scoring::etable::Etable etable(atypeset,eopt);

	bool allonefile = false;
	std::string motif_out_file = basic::options::option[basic::options::OptionKeys::mh::motif_out_file]()+".rpm.bin.gz";
	utility::io::ozstream allout;
	if ( basic::options::option[basic::options::OptionKeys::mh::motif_out_file].user() ) {
		if ( utility::file::file_exists(motif_out_file) ) utility_exit_with_message("already done!");
		allonefile = true;
		allout.open(motif_out_file,std::ios::out|std::ios::binary);
		if ( !allout.good() ) {
			allout.close();
			utility_exit_with_message("error opening "+motif_out_file );
		}
	}

	vector1<std::string> const & fnames( option[mh::harvest_motifs]() );
	for ( int ifile = 1; ifile <= (int)fnames.size(); ++ifile ) {
		ResPairMotifs motifs;
		std::string fname = fnames[ifile];
		std::string const tag = tag_from_pdb_fname(fname);
		if ( ""==tag ) {
			TR.Error << "WARNING: skipping non-pdb fname: " << fname << endl;
			continue;
		}

		if ( !allonefile && utility::file::file_exists(option[out::file::o]()+"/"+tag.substr(1,2)+"/"+tag+".rpm.bin.gz") ) continue;

		Pose pose; vector1<Real> bfactors,occupancy;
		utility::vector1<int> pdbres;
		std::map<int,char> pdbchain;
		int nresmodel1;
		if ( !protocols::sic_dock::read_biounit(fname,pose,bfactors,occupancy,pdbres,pdbchain,nresmodel1,option[mh::harvest::max_res](),false) ) {
			TR.Error << "FAIL TO READ " << fname << endl;
			continue;
		}
		// core::chemical::ResidueType const & gly_type( (pose.residue(1).residue_type_set().name_map("GLY") ));
		// Pose pose_no_ir(pose), pose_no_jr(pose), pose_no_ir_jr(pose);
		ScoreFunctionOP sf = core::scoring::get_score_function();
		core::scoring::methods::EnergyMethodOptions myopt = sf->energy_method_options();
		myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
		sf->set_energy_method_options(myopt);


		core::scoring::dssp::Dssp dssp(pose);
		dssp.insert_edge_ss_into_pose(pose);
		sf->score(pose);
		core::scoring::hbonds::HBondSet hbset;
		core::scoring::hbonds::fill_hbond_set(pose,false,hbset);
		Energies    const & energies     ( pose.energies() );
		EnergyGraph const & energy_graph ( energies.energy_graph() );
		vector1<Real> rsd_sasa = get_sasa(pose,1.7);
		vector1<Real> nbrs = get_nbrs(pose);

		// int ssbdy = option[mh::harvest_motifs_min_hh_ends]();
		int min_boundary = 0;//ssbdy;

		// scan pairs
		int nbad=0,nbadbfac=0,nbadocc=0,nbadiface=0,nbaddist=0,nbadsselem=0;
		for ( int ir = 1+min_boundary; ir <= nresmodel1-min_boundary; ++ir ) {
			if ( !pose.residue(ir).is_protein() ) { ++nbad; continue; }
			for ( utility::graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node(ir)->const_upper_edge_list_begin(),
					irue = energy_graph.get_node(ir)->const_upper_edge_list_end();
					iru != irue; ++iru
					) {
				EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
				harvest_motifs_one(
					motifs,
					tag,
					pose,
					dssp,
					//etable,
					ir,
					edge,
					occupancy,
					bfactors,
					rsd_sasa,
					nbrs,
					nbad,
					nbadbfac,
					nbadocc,
					nbadiface,
					nbaddist
				);
			}

		}


		// do I/O
		if ( option[mh::filter::filter_harvest]() ) filter_motifs(motifs);
		TR << fname <<' '<< motifs.size() <<" nbad: "<< I(5,nbad) <<" nbadbfac: "<< I(5,nbadbfac) <<" nbadocc: "<< I(5,nbadocc)
			<<" nbadiface: "<< I(5,nbadiface) <<" nbaddist: "<< I(5,nbaddist) <<" nbadsselem: "<<I(5,nbadsselem)<< endl;
		if ( allonefile ) {
			TR << "appending " << motifs.size() << " motifs to " << motif_out_file+".rpm.bin.gz" << endl;
			if ( !write_motifs_binary(allout,motifs) ) utility_exit_with_message("error writing to file "+motif_out_file+".rpm.bin.gz");
		} else {
			std::string outfile = option[out::file::o]()+"/"/*+tag.substr(1,2)+"/"*/+tag+".rpm.bin.gz";
			TR << "dumping " << motifs.size() << " motifs to " << outfile << endl;
			write_motifs_binary(outfile,motifs);
			if ( !utility::file::file_exists(outfile) ) utility_exit_with_message("couldn't make file "+outfile);
		}

	}
	// DONE: ;
	if ( allonefile ) {
		if ( !allout.good() ) {
			allout.close();
			utility_exit_with_message("error after writing to "+motif_out_file+".rpm.bin.gz" );
		}
		allout.close();
	}
	TR << "DONE HARVESTING" << endl;
}

int main(int argc, char *argv[]) {
	try{

		devel::init(argc,argv);

		using namespace basic::options::OptionKeys;
		using namespace core::scoring;
		using namespace core::scoring::motif;

		if ( 64!=sizeof(ResPairMotif) ) utility_exit_with_message("ResPairMotif should be 64 bytes in size!!!!!!");

		if ( option[mh::dump_input_pdb].user() ) {
			TR << "dumping input_pdbs from  -mh:dump_input_pdb" << endl;
			vector1<string> const & fnames( option[mh::dump_input_pdb]() );
			for ( int ifile = 1; ifile <= (int)fnames.size(); ++ifile ) {
				Pose pose; vector1<Real> bfactors,occupancy;
				utility::vector1<int> pdbres;
				std::map<int,char> pdbchain;
				int nresmodel1;
				if ( !protocols::sic_dock::read_biounit(fnames[ifile],pose,bfactors,occupancy,pdbres,pdbchain,nresmodel1,option[mh::harvest::max_res](),true) ) {
					TR.Error << "FAIL TO READ " << fnames[ifile] << endl;
					continue;
				}
			}
		}

		if ( option[mh::harvest_motifs].user() ) {
			TR << "harvesting motifs in files following -mh:harvest_motifs" << endl;
			harvest_motifs();
		}
		if ( option[mh::print_motifs].user() ) {
			TR << "reading motifs in files following -mh:print_motifs" << endl;
			ResPairMotif::print_header(cout);
			print_motifs(cout);
		}
		if ( option[mh::merge_motifs].user() ) {
			TR << "reading motifs in files following -mh:merge_motifs" << endl;
			merge_motifs();
		}
		if ( option[mh::harvest_scores].user() ) {
			TR << "reading motif counts in files following -mh:harvest_scores" << endl;
			harvest_scores();
		}
		if ( option[mh::print_scores].user() ) {
			TR << "reading motifs in files following -mh:print_scores" << endl;
			print_scores();
		}



	} catch ( utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		std::exit( 1 );
	}
	return 0;
}


