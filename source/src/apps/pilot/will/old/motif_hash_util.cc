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
	#include <protocols/simple_moves/MinMover.hh>
	#include <protocols/simple_moves/AddConstraintsToCurrentConformationMover.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/fixedsizearray1.hh>
	#include <utility/file/file_sys_util.hh>
	#include <numeric/geometry/hashing/SixDHasher.hh>
	#include <numeric/HomogeneousTransform.hh>

	#include <apps/pilot/will/will_util.ihh>

	// #include <boost/archive/text_oarchive.hpp>
	// #include <boost/archive/text_iarchive.hpp>
	// #include <boost/archive/binary_oarchive.hpp>
	// #include <boost/archive/binary_iarchive.hpp>

	// // #include <boost/serialization/detail/stack_constructor.hpp>
	// #include <boost/serialization/hash_map.hpp>
	// #include <boost/serialization/hash_collections_save_imp.hpp>
	// #include <boost/serialization/hash_collections_load_imp.hpp>

static THREAD_LOCAL basic::Tracer TR( "motif_hash_util" );

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


#ifdef USE_OPENMP
	#include <omp.h>
	#endif
	Size sicdock_max_threads(){
		#ifdef USE_OPENMP
			return omp_get_max_threads();
		#else
			return 1;
		#endif
	 }
	Size sicdock_num_threads(){
		#ifdef USE_OPENMP
			return omp_get_num_threads();
		#else
			return 1;
		#endif
	 }
	Size sicdock_thread_num(){
		#ifdef USE_OPENMP
			return omp_get_thread_num() + 1;
		#else
			return 1;
		#endif
	 }


char dssp_reduced(char const & c){
	if( c == 'H' || c == 'G' || c == 'I' ) {
		return 'H';
    } else if(c == 'B' || c == 'E') {
	 	return 'E';
	} else return 'L';
 }

void print_motifs(std::ostream & out){
	using namespace basic::options::OptionKeys;
	ResPairMotifs motifs;
	load_motifs( option[mh::print_motifs](), motifs );
	for(ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i){
		out << *i << endl;
	 }
 }

void remove_duplicates(){
	using namespace basic::options::OptionKeys;
	ResPairMotifs motifs;
	load_motifs( option[mh::remove_duplicates](), motifs );

	std::set<ResPairMotif> s;
	for ( ResPairMotif const & rpm : motifs ) s.insert(rpm);

	cout << motifs.size() << " reduced to " << s.size() << endl;

	ResPairMotifs motifs2;
	motifs2.insert(motifs2.begin(),s.begin(),s.end());

	write_motifs_binary(option[mh::motif_out_file]()+".rpm.bin.gz",motifs2);

 }

void merge_motifs(){
	using namespace basic::options::OptionKeys;
	if( ! basic::options::option[basic::options::OptionKeys::mh::motif_out_file].user() ){
		utility_exit_with_message("must sepcify name for merged file -motif_out_file");
	 }
	vector1<string> const & fnames ( option[mh::merge_motifs]() );
	ResPairMotifs motifs; load_motifs( fnames, motifs );
	// TR << motifs.size() << endl;
	// motifs.filter_structurally_identical_motifs();
	// TR << motifs.size() << endl;

	ResPairMotifMetaBinner binner;
	ResPairMotifsMap mmap;

	for(Size i = 1; i <= motifs.size(); ++i){
		ResPairMotif const & m(motifs[i]);
		ResPairMotifMetaBinner::Key k = binner.motif_bin_hash(m);
		mmap[k].push_back(m);
	 }

	TR << "write to disk" << endl;
	for (ResPairMotifsMap::value_type const & p : mmap){
		std::string label = binner.hash_to_labal(p.first);
		write_motifs_binary(option[mh::motif_out_file]()+(label.size()?"_":"")+label,p.second);
	 }
	TR << "merge_motifs done" << endl;

 }

void merge_scores(){
	if( ! basic::options::option[basic::options::OptionKeys::mh::motif_out_file].user() ){
		utility_exit_with_message("must sepcify name for merged file -motif_out_file");
	 }
	vector1<string> const & fnames ( option[mh::merge_scores]() );
	string          const & outfile( option[mh::motif_out_file ]()+".xh.bin.gz" );

	XformScore *xs;
	XformScore::read_binary(xs,fnames);
	xs->prune_small_bins(option[mh::harvest::min_bin_val]());
	cout << *xs << endl;
	xs->write_binary(outfile);
 }

void dump_motif_pdbs(){
	using namespace basic::options::OptionKeys;

	TR << "dump_motif_pdbs: read motifs" << endl;
	vector1<string> const & fnames( option[mh::dump_motif_pdbs]() );
	ResPairMotifs motifs; load_motifs( fnames, motifs );

	TR << "dump_motif_pdbs: hash motifs" << endl;
	MotifHash mh;
	for(ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i){
		// TR << "ADD_KEY " << mh.bin_index(*i) << " " <<  i->rt() << endl;
		mh.add_motif(*i);
	 }

	TR << "dump_motif_pdbs: find largest bucket / check counts" << endl;
	uint64_t maxnum=0,mxkey;
	vector1<uint64_t> keys;
	Real mincount = 9e9;
	if( option[mh::harvest::min_bin_val].user() ) mincount = option[mh::harvest::min_bin_val]();
	for (uint64_t key : mh.keys() ){
		uint64_t n = mh.count_motifs(key);
		if( n >= mincount ){
			keys.push_back(key);
		}
		if( n >= maxnum ){
			maxnum = n;
			mxkey = key;
		}
	 }
	if(keys.size()==0) keys.push_back(mxkey);
	for (uint64_t key : keys) TR << "key " << key << endl;
	TR << "dump_motif_pdbs: max motifs in bucket, dumping that one bucket: " << maxnum << " " << (Real)maxnum/(Real)motifs.size() << endl;
	string tag = string(option[mh::motif_out_file]());
	if(tag=="") tag = "NOTAG";
	ResPairMotif::print_header(TR);
	for(Size k = 1; k <= keys.size(); ++k){
		uint64_t key = keys[k];
		MotifHash::MotifMap::const_iterator l = mh.motif_umap_.equal_range(key).first;
		MotifHash::MotifMap::const_iterator e = mh.motif_umap_.equal_range(key).second;
		string const fn(tag+"_"+string_of(key)+".pdb.gz");
		utility::io::ozstream out(fn);
		Size count(0);
		for(; l != e; ++l,++count){
			if(count > (Size)option[mh::dump::limit_per_pair]()) continue;
			ResPairMotif const & sm(l->second);
			string mdl = string_of(sm.ss1())    +
			             string_of(sm.ss2())    + "_"+
			             string_of(sm.aa1())    +
			             string_of(sm.aa2())    + "_"+
			             sm.pdb()               + "_"+
			             string_of(sm.resi1() ) + "_"+
			             string_of(sm.resi2() ) ;
			// TR << fn << " " << sm << endl;
			l->second.dump_pdb(out,Xform::identity(),mdl);
		}
		out.close();
		TR << "dumped " << count << " motifs to " << fn << endl;
	 }
 }

template<class T> T sqr(T const & x){ return x*x; }

void harvest_scores(){
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::motif;
	using core::scoring::motif::aa_trustworthiness;

	TR << "harvest_scores: read motifs" << endl;
	vector1<string> const & fnames( option[mh::harvest_scores]() );
	ResPairMotifs motifs; load_motifs( fnames, motifs );

	RPM_Type type = core::scoring::motif::RPM_Type_NONE;
	if( option[mh::filter::motif_type].user() && option[mh::filter::motif_type]()=="SC_SC" ) type = SC_SC;
	if( option[mh::filter::motif_type].user() && option[mh::filter::motif_type]()=="SC_BB" ) type = SC_BB;
	if( option[mh::filter::motif_type].user() && option[mh::filter::motif_type]()=="SC_PH" ) type = SC_PH;
	if( option[mh::filter::motif_type].user() && option[mh::filter::motif_type]()=="SC_PO" ) type = SC_PO;
	if( option[mh::filter::motif_type].user() && option[mh::filter::motif_type]()=="BB_BB" ) type = BB_BB;
	if( option[mh::filter::motif_type].user() && option[mh::filter::motif_type]()=="BB_PH" ) type = BB_PH;
	if( option[mh::filter::motif_type].user() && option[mh::filter::motif_type]()=="BB_PO" ) type = BB_PO;
	if( option[mh::filter::motif_type].user() && option[mh::filter::motif_type]()=="PH_PO" ) type = PH_PO;

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

	for(double d2 = 0.0; d2 < 3.0; d2 += 0.5){
		cout << "WT " << d2 << " " << exp(-d2/SMOOTHING_FAC/SMOOTHING_FAC) * EXP_NORM << endl;
	}

	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for(int i = 1; i <= (int)motifs.size(); ++i){
		if(i%(motifs.size()/10)==0) cout << (float)i/(float)motifs.size()*100 << "%" << endl;
		ResPairMotif const & m(motifs[i]);
		if(type==RPM_Type_NONE) type=m.type();
		if(m.type()!=type) utility_exit_with_message("harvest_scores: mismatched types");
		Real6 const rt6 = m.rt();

		ResPairMotifMetaBinner::Key key1 = binner.motif_bin_hash(m);
		assert( binner.hash_to_labal(key1) == binner.motif_bin_label(m) );
		// cout << binner.hash_to_labal(key1) << " " << binner.motif_bin_label(m) << endl;

 		// find correct score to contribute to, create if necessary
		XformScoreMap::iterator xscoreiter = xscoremap[sicdock_thread_num()].find(key1);
		if(xscoreiter==xscoremap[sicdock_thread_num()].end()){
			xscoremap[sicdock_thread_num()][key1] = new XformScore(option[mh::harvest::hash_cart_resl](),option[mh::harvest::hash_angle_resl]());
			xscoreiter = xscoremap[sicdock_thread_num()].find(key1);
			Real numbins = 1.0;
			for(int i = 1; i <= 6; ++i) numbins *= (Real)xscoreiter->second->hasher_.dimsizes()[i];
			if( numbins > maxbins ){
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
			if(basic::options::option[basic::options::OptionKeys::mh::harvest::weight_by_energy]()){
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
			for(numeric::Size i = 1; i <= 3; ++i){
				beg[i] = max(                 1             ,(int)bin6[i]-1);
				end[i] = min((int)xsfx.hasher_.dimsizes()[i],(int)bin6[i]+1)+1;
			}
			// 4 & 5 wrap for euler angles
			for(numeric::Size i = 4; i <= 5; ++i){
				beg[i] = bin6[i] > 1                          ? bin6[i]-1 : xsfx.hasher_.dimsizes()[i] ;
				end[i] = bin6[i] < xsfx.hasher_.dimsizes()[i] ? bin6[i]+1 : 1;
				end[i] = end[i] < xsfx.hasher_.dimsizes()[i] ? end[i]+1 : 1;
			}
			beg[6] = max(                 1             ,(int)bin6[6]-1);
			end[6] = min((int)xsfx.hasher_.dimsizes()[6],(int)bin6[6]+1)+1;
			// 	if( half[i] && (i==4||i==5||end[i] < xsfx.hasher_.dimsizes()[i]-1)) ++end[i];
			// 	if(!half[i] && (i==4||i==5||beg[i] >                0       )) --beg[i];
			// }
			// cout <<I(5,xsfx.hasher_.dimsizes()[4])<<I(5,xsfx.hasher_.dimsizes()[5])<<I(5,xsfx.hasher_.dimsizes()[6])<<I(5,xsfx.hasher_.dimsizes()[1])<<I(5,xsfx.hasher_.dimsizes()[2])<<I(5,xsfx.hasher_.dimsizes()[3])<<endl;
			// cout <<I(5, beg[4])<<I(5, beg[5])<<I(5, beg[6])<<I(5, beg[1])<<I(5, beg[2])<<I(5, beg[3])<<endl;
			// cout <<I(5,bin6[4])<<I(5,bin6[5])<<I(5,bin6[6])<<I(5,bin6[1])<<I(5,bin6[2])<<I(5,bin6[3])<<endl;
			// cout <<I(5, end[4])<<I(5, end[5])<<I(5, end[6])<<I(5, end[1])<<I(5, end[2])<<I(5, end[3])<<endl;
			// cout << endl;

			utility::fixedsizearray1<Size,6> ob6;
			// int count_this_rt=0;
			for(ob6[1] = beg[1]; ob6[1] != end[1]; ++ob6[1])
			for(ob6[2] = beg[2]; ob6[2] != end[2]; ++ob6[2])
			for(ob6[3] = beg[3]; ob6[3] != end[3]; ++ob6[3])
			for(ob6[4] = beg[4]; ob6[4] != end[4]; ob6[4] = (ob6[4]==xsfx.hasher_.dimsizes()[4]) ? 1 : ob6[4]+1)
			for(ob6[5] = beg[5]; ob6[5] != end[5]; ob6[5] = (ob6[5]==xsfx.hasher_.dimsizes()[5]) ? 1 : ob6[5]+1)
			for(ob6[6] = beg[6]; ob6[6] != end[6]; ++ob6[6]){
				// cout <<"ITER "<<I(5,ob6[4])<<I(5,ob6[5])<<I(5,ob6[6])<<I(5,ob6[1])<<I(5,ob6[2])<<I(5,ob6[3])<<endl;
				XformScore::Key other_key = xsfx.hasher_.bin_index(ob6);


				// // get target cells
				// // stupid thing about the SixDofHasher, must do lookups concentrically... shell 5 is the max... speed up by lowering it
				// vector1<XformScore::Key> surrounding_cells;
				// for(Size lookuprad = 0; lookuprad <= HASH_LOOKUP_RADIUS; ++lookuprad){
				// 	std::vector<uint64_t> shell_tmp = xsfx.hasher_.radial_bin_index(lookuprad,rt6);
				// 	// cout << "RADIAL_BIN_LOOKUP " << lookuprad << " " << shell_tmp.size() << endl;

				// 	 // loop over all cells close to this motif
				// 	 for (XformScore::Key const & other_key :  shell_tmp ){
				//  	// runtime_assert_msg(testkeys.find(other_key)==testkeys.end(),"key already seen!");
				// 	// testkeys.insert(other_key);

						Real6 const other_rt6 = xsfx.hasher_.bin_center_point(xsfx.hasher_.bin_from_index(other_key));
						float rt6_dis2 = 0;
						for(int ii=1; ii<=3; ++ii) rt6_dis2 += sqr( (rt6[ii]-other_rt6[ii])/HASH_CART_RESL );
						for(int ii=4; ii<=6; ++ii) rt6_dis2 += sqr( (rt6[ii]-other_rt6[ii])/HASH_ANGL_RESL );
						Real const wt = exp(-rt6_dis2/SMOOTHING_FAC/SMOOTHING_FAC) * EXP_NORM;
						if(wt < 0.2) continue;
						// cout << "WT " << lookuprad << " " << rt6_dis2 << " " << wt << endl;
						Real const sc = wt*motif_score;
						if( sc < option[mh::harvest::min_bin_val]()/30.0 ) continue;
						if(agg_with_max) xsfx.max_score(other_key,-m.score());
						else             xsfx.add_score(other_key,wt*motif_score);
				// }

			}
		// }
	 }

	#ifdef USE_OPENMP
	int mergecount=0;
	XformScoreMap xscoremap_joint;
	for (XformScoreMap & m : xscoremap){
		TR << "merge threaded data " << ++mergecount << " / " << xscoremap.size() << endl;
		for (XformScoreMap::value_type const & v : m){
			if(xscoremap_joint.find(v.first)==xscoremap_joint.end())
				xscoremap_joint[v.first] = new XformScore(option[mh::harvest::hash_cart_resl](),option[mh::harvest::hash_angle_resl]());
			if(agg_with_max) xscoremap_joint[v.first]->aggregate_max(*v.second);
			else             xscoremap_joint[v.first]->aggregate_add(*v.second);
		}
		m.clear();
	}
	#else
	XformScoreMap & xscoremap_joint(xscoremap.front());
	#endif

	if(option[mh::harvest::min_bin_val].user()){
		TR << "prune" << endl;
		for (XformScoreMap::value_type& v : xscoremap_joint) v.second->prune_small_bins(option[mh::harvest::min_bin_val]());
	 }

	TR << xscoremap_joint;

	TR << "write to disk" << endl;
	string PATH = option[mh::motif_out_file]();
	utility::file::create_directory_recursive(PATH);
	for (XformScoreMap::value_type const & p : xscoremap_joint){
		std::string label = binner.hash_to_labal(p.first);
		p.second->write_binary(PATH+"/"+PATH+(label.size()?"_":"")+label+".xh.bin.gz");
	 }
	if(option[mh::harvest::dump]()) write_motifs_binary(PATH+"/"+PATH+".rpm.bin.gz",motifs);
	TR << PATH << endl;
	TR << "harvest_scores done" << endl;
 }

void print_scores(){
	using namespace basic::options::OptionKeys;
	XformScore *xh = NULL;
	XformScore::read_binary(xh,option[mh::print_scores]());
	TR << *xh << endl;
	xh->print_scores(cout);
 }

void dump_matching_motifs(){
	using namespace basic::options::OptionKeys;
	using core::scoring::motif::MotifHashCAP;
	using core::scoring::motif::MotifHits;
	using core::scoring::motif::MotifHashManager;
	using core::scoring::motif::ResPairMotifQuery;

	vector1<string> const & pdbfiles  ( option[mh::dump_matching_motifs]() );
	vector1<string> pdbfiles2;
	if(option[in::file::s].user()) pdbfiles2 = option[in::file::s]();

	for(int ifile = 1; ifile <= (int)pdbfiles.size(); ++ifile ){
		string fname = pdbfiles[ifile];
		std::string outfile = utility::file_basename(fname) + "_motifs.pdb.gz";
		core::pose::PoseOP pose = core::import_pose::pose_from_file(fname, core::import_pose::PDB_file);
		if(option[basic::options::OptionKeys::symmetry::symmetry_definition].user())
			core::pose::symmetry::make_symmetric_pose(*pose);
		core::pose::PoseOP pose2 = pose;
		if(pdbfiles2.size()){
			cout << "2nd file " << pdbfiles2[ifile] << endl;
			pose2 = core::import_pose::pose_from_file(pdbfiles2[ifile], core::import_pose::PDB_file);
		}
		ResPairMotifQuery opt(*pose,*pose2);
		cout << opt << endl;
		MotifHits hits;
		MotifHashManager::get_instance()->get_matching_motifs(opt,hits);
		if( hits.size()==0 ){
			TR << "no hits for " << fname << endl;
		} else {
			TR << "dump " << hits.size() << " (before pruning?) hits to " << outfile << endl;
			utility::io::ozstream pdbout( outfile );
			pdbout << "MODEL MAIN" << endl;
			pose->dump_pdb(pdbout);
			if(pose2!=pose) pose2->dump_pdb(pdbout);
			pdbout << "ENDMDL" << endl;
			hits.dump_motifs_pdb(pdbout);
			pdbout.close();
		}

	 }
 }

Real get_etable_score_for_atoms(
	core::scoring::etable::Etable const & etable,
	core::pose::Pose const & pose1, utility::vector1<core::id::AtomID> const & aids1,
	core::pose::Pose const & /*pose2*/, utility::vector1<core::id::AtomID> const & aids2
 ){
	using namespace core::id;
	Real sc = 0.0;
	for (AtomID const & aid1 : aids1){
		for (AtomID const & aid2 : aids2){
			Real lj_atrE,lj_repE,fa_solE,d2;
			etable.analytic_etable_evaluation(
				pose1.residue(aid1.rsd()).atom(aid1.atomno()),
				pose1.residue(aid2.rsd()).atom(aid2.atomno()),
				lj_atrE, lj_repE, fa_solE, d2);
			sc += 1.0*lj_atrE + 0.2*lj_repE + 0.8*fa_solE;
			// cout << d2 << endl;
		}
	}
	return sc;
 }

bool harvest_motifs_one(
	ResPairMotifs & motifs,
	string const & tag0,
	Pose const & pose,
	// Pose const & pose_no_ir,
	// Pose const & pose_no_jr,
	// Pose const & pose_no_ir_jr,
	core::scoring::dssp::Dssp & dssp,
	core::scoring::hbonds::HBondSet const & hbset,
	core::scoring::etable::Etable const & etable,
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
	// if( ir==1 || jr == pose.size() ) return false;
	// Size const ichain = pose.chain(ir);
	// Size const jchain = pose.chain(jr);

	// if( abs((int)ir-(int)jr) < 10 && ichain==jchain ) { ++nbad; ++nbaddist; return false; } // no same or adjcent res

	if(!pose.residue(ir).is_protein()) { ++nbad; ++nbadocc ; return false; }
	if(!pose.residue(jr).is_protein()) { ++nbad; ++nbadocc ; return false; }
	// if( bfactors [ir] <= 0.00 || bfactors [ir] > option[mh::filter::coorderr]() ) { ++nbad; ++nbadbfac; return false; }
	// if( bfactors [jr] <= 0.00 || bfactors [jr] > option[mh::filter::coorderr]() ) { ++nbad; ++nbadbfac; return false; }
	// if( occupancy[ir] <= 0.99 || occupancy[ir] >      1.01                      ) { ++nbad; ++nbadocc ; return false; }
	// if( occupancy[jr] <= 0.99 || occupancy[jr] >      1.01                      ) { ++nbad; ++nbadocc ; return false; }

	// if( option[mh::filter::maxdist2]() < pose.residue(ir).xyz("CA").distance_squared(pose.residue(jr).xyz("CA")) ) { ++nbad; ++nbaddist; return false; }

	bool const is_disulf = core::conformation::is_disulfide_bond(pose.conformation(),ir,jr);
	bool const is_sspaired = dssp.in_paired_strands(ir,jr) ;
	Real const hb_sc = edge[core::scoring::hbond_sc];

	string tag = tag0;
	if(is_disulf) tag = tag+"dsf";

	Real hbsc_nbb_ir=0, hbsc_cbb_ir=0, hbsc_nbb_jr=0, hbsc_cbb_jr=0, hbbb_i_to_j=0, hbbb_j_to_i=0;
	for(Size ihb = 1; ihb <= hbset.nhbonds(); ++ihb){
		HBond const & hb(hbset.hbond(ihb));
		if( hb.don_res()==ir && hb.acc_res()==jr ){
			// hbonds_i_don_j.push_back(&hb);
			if( hb.don_hatm_is_protein_backbone() && hb.acc_atm_is_protein_backbone() ) hbbb_i_to_j += hb.energy();
			else if( hb.don_hatm_is_protein_backbone()                                ) hbsc_nbb_ir += hb.energy();
			else if( hb.acc_atm_is_protein_backbone()                                 ) hbsc_cbb_jr += hb.energy();
		}
		if( hb.don_res()==jr && hb.acc_res()==ir ){
			// hbonds_j_don_i.push_back(&hb);
			if( hb.don_hatm_is_protein_backbone() && hb.acc_atm_is_protein_backbone() ) hbbb_j_to_i += hb.energy();
			else if( hb.don_hatm_is_protein_backbone()                                ) hbsc_nbb_jr += hb.energy();
			else if( hb.acc_atm_is_protein_backbone()                                 ) hbsc_cbb_ir += hb.energy();
		}
	}

	Real const nb1 = nbrs[ir];
	Real const nb2 = nbrs[jr];
	Real const & bfir(bfactors[ir]);
	Real const & bfjr(bfactors[jr]);
	// cout << I(4,ir) << " " << I(4,jr) << " " << F(7,3,atrss) << " " << F(7,3,fa_atr_0) << " " << F(7,3,fa_atr_1) << " " << F(7,3,fa_atr_2) << " " << F(7,3,fa_atr_3) << endl;
	// if( atrss < -0.3 && fa_atr_1 > -0.001 && fa_atr_2 > -0.001 ){
	// 	pose.dump_pdb("test.pdb");
	// 	utility_exit_with_message("test faatr bb corr");
	// }

	//////////////////////////////////// actual motifs ////////////////////////////////////////////
		Real E_TH = -0.001;

		Xform const bbr1 = core::pose::motif::get_backbone_reference_frame(pose,ir);
		Xform const bbr2 = core::pose::motif::get_backbone_reference_frame(pose,jr);
		if( bbr1.bad() || bbr2.bad() ) {
			cout << "bad bbr " << ir << " " << jr << endl;
			return false;
		}

		static string DBGDIR = option[out::file::o]()+"/";
		utility::file::create_directory_recursive(DBGDIR);

		static bool header=true;

		///////////////////////////////////// bb h-bond motifs ///////////////////////////
		for(int ij_or_ji=0; ij_or_ji <= 1; ++ij_or_ji ){
			Real const & hbe = ij_or_ji? hbbb_i_to_j : hbbb_j_to_i;
			if(hbe > -0.001) continue;
			Xform const pbr1 = core::pose::motif::get_peptide_bond_reference_frame(pose,ir, ij_or_ji);
			Xform const pbr2 = core::pose::motif::get_peptide_bond_reference_frame(pose,jr,!ij_or_ji);
			if( pbr1.bad() || pbr2.bad() ){	cout << "bad pbr " << ir << " " << jr << endl; continue; }
			AIDs const aid1 = core::pose::motif::get_peptide_bond_reference_frame_atomids(pose,ir, ij_or_ji);
			AIDs const aid2 = core::pose::motif::get_peptide_bond_reference_frame_atomids(pose,jr,!ij_or_ji);
			Real const scorebb = get_etable_score_for_atoms(etable,pose,aid1,pose,aid2);
			ResPairMotif motif(tag, pose, (~pbr1*pbr2).rt6(), ir, jr, nb1, nb2, 0, 0, scorebb, 0, 0, hbe, bfactors[ir], bfactors[jr], is_sspaired, PH_PO );
			if( !ij_or_ji ) motif.reverse_in_place_unsafe();
			// if( motif.score() < E_TH  ){
				motifs.push_back(motif);
				if(option[mh::harvest::dump]()){
					if(header){ motif.print_header(cout); header=false; }
					cout << motif << endl;
					Real rms = motif.dump_aligned_motif(DBGDIR+motif.tag()+".pdb",pose,ij_or_ji?ir:jr,ij_or_ji?jr:ir);
					if(rms > 0.01) utility_exit_with_message("bad rms "+string_of(rms)+" "+DBGDIR+motif.tag()+".pdb");
				}
			// }
		}

		///////////////////////////////////// bb /sc h-bond motifs ///////////////////////////
		for(Size iscbb_t = 0; iscbb_t < 4; ++iscbb_t){
			Real hbbbsc_e;
			Xform frm1,frm2;
			AIDs aids1,aids2;
			RPM_Type t; bool rev=false;
			switch(iscbb_t){
				case 0: t=BB_PH; frm2=bbr2; hbbbsc_e=hbsc_nbb_ir; rev=true;
					frm1 = get_nterminal_peptide_bond_reference_frame(pose,ir);
					aids1 = core::pose::motif::get_nterminal_peptide_bond_reference_frame_atomids(pose,ir);
					aids2 = core::pose::motif::get_backbone_reference_frame_atomids_with_downstream(pose,jr);
					break;
				case 1: t=BB_PO; frm2=bbr2; hbbbsc_e=hbsc_cbb_ir; rev=true;
					frm1 = get_cterminal_peptide_bond_reference_frame(pose,ir);
					aids1 = core::pose::motif::get_cterminal_peptide_bond_reference_frame_atomids(pose,ir);
					aids2 = core::pose::motif::get_backbone_reference_frame_atomids_with_downstream(pose,jr);
					break;
				case 2:	 t=BB_PH; frm1=bbr1; hbbbsc_e=hbsc_nbb_jr;
					frm2 = get_nterminal_peptide_bond_reference_frame(pose,jr);
					aids1 = core::pose::motif::get_backbone_reference_frame_atomids_with_downstream(pose,ir);
					aids2 = core::pose::motif::get_nterminal_peptide_bond_reference_frame_atomids(pose,jr);
					break;
				case 3:	 t=BB_PO; frm1=bbr1; hbbbsc_e=hbsc_cbb_jr;
					frm2 =  get_cterminal_peptide_bond_reference_frame(pose,jr);
					aids1 = core::pose::motif::get_backbone_reference_frame_atomids_with_downstream(pose,ir);
					aids2 = core::pose::motif::get_cterminal_peptide_bond_reference_frame_atomids(pose,jr);
					break;
				default: utility_exit_with_message("madness!");
			}
			if( hbbbsc_e > -0.001 ) continue;
			if( frm1.bad() ){ cout << "bad frame " << iscbb_t << " " << t << " " << ir << " " << pose.residue(ir).name() << endl; continue; }
			if( frm2.bad() ){ cout << "bad frame " << iscbb_t << " " << t << " " << jr << " " << pose.residue(jr).name() << endl; continue; }
			Real const score = get_etable_score_for_atoms(etable,pose,aids1,pose,aids2);
			ResPairMotif motif(tag, pose, (~frm1*frm2).rt6(), ir, jr, nb1, nb2, 0, score, 0, 0, hbbbsc_e, 0, bfir, bfjr, is_sspaired, t );
			if(rev) motif.reverse_in_place_unsafe();
			// if( motif.score() < E_TH ){
				motifs.push_back(motif);
				if(option[mh::harvest::dump]()){
					if(header){ motif.print_header(cout); header=false; }
					cout << motif << endl;
					Real rms = motif.dump_aligned_motif(DBGDIR+motif.tag()+".pdb",pose,rev?jr:ir,rev?ir:jr);
					if(rms > 0.01) utility_exit_with_message("bad rms "+string_of(rms)+" "+DBGDIR+motif.tag()+".pdb");
				}
			// }
		}

	 	////////////////////////////////// 'standard' motifs + sc frame motifs ///////////////////////////////

		Xform const scr1 = core::pose::motif::get_sidechain_reference_frame(pose,ir);
		Xform const scr2 = core::pose::motif::get_sidechain_reference_frame(pose,jr);

		AIDs bb1aids,bb2aids,sc1aids,sc2aids;
		for(Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia){
			if(pose.residue(ir).atom_is_backbone(ia)) bb1aids.push_back(AtomID(ia,ir));
			else                                      sc1aids.push_back(AtomID(ia,ir));
		}
		for(Size ia = 1; ia <= pose.residue(jr).natoms(); ++ia){
			if(pose.residue(jr).atom_is_backbone(ia)) bb2aids.push_back(AtomID(ia,jr));
			else                                      sc2aids.push_back(AtomID(ia,jr));
		}
		AIDs sr1aids = get_sidechain_reference_frame_atomids(pose,ir);
		AIDs sr2aids = get_sidechain_reference_frame_atomids(pose,jr);
		Real const et_scsc = get_etable_score_for_atoms(etable,pose,sc1aids,pose,sc2aids);
		Real const et_scbb = get_etable_score_for_atoms(etable,pose,sc1aids,pose,bb2aids);
		Real const et_bbsc = get_etable_score_for_atoms(etable,pose,bb1aids,pose,sc2aids);
		Real const et_bbbb = get_etable_score_for_atoms(etable,pose,bb1aids,pose,bb2aids);
		Real const et_srsr = get_etable_score_for_atoms(etable,pose,sr1aids,pose,sr2aids);
		Real const et_srbb = get_etable_score_for_atoms(etable,pose,sr1aids,pose,bb2aids);
		Real const et_bbsr = get_etable_score_for_atoms(etable,pose,bb1aids,pose,sr2aids);
		Real const et_srsc = get_etable_score_for_atoms(etable,pose,sr1aids,pose,sc2aids);
		Real const et_scsr = get_etable_score_for_atoms(etable,pose,sc1aids,pose,sr2aids);
		// cout << et_scsc << " " << et_scbb << " " << et_bbsc << " " << et_bbbb << endl;


		ResPairMotif bbbb( tag, pose, (~bbr1*bbr2).rt6(), ir, jr, nb1, nb2, et_scsc, et_bbsc+et_scbb, et_bbbb, hb_sc, 0, 0, bfir, bfjr, is_sspaired, BB_BB );
		if( is_disulf || (/*bbbb.score() < E_TH &&*/ ( et_scsc < E_TH || hb_sc < E_TH ) ) ){
			motifs.push_back(bbbb);
			if(option[mh::harvest::dump]()){
				if(header){ bbbb.print_header(cout); header=false; }
				cout << bbbb << endl;
				Real rms = bbbb.dump_aligned_motif(DBGDIR+bbbb.tag()+".pdb",pose,ir,jr);
				if(rms > 0.01) utility_exit_with_message("bad rms "+string_of(rms)+" "+DBGDIR+bbbb.tag()+".pdb");
			}
		}
		if( scr1.bad() || scr2.bad() ) { cout << "bad sc frame " << ir << ' ' << jr << endl; return false; }
		ResPairMotif bbsc( tag, pose, (~bbr1*scr2).rt6(), ir, jr, nb1, nb2, et_scsr, et_bbsr, 0, hb_sc, 0, 0, bfir, bfjr, is_sspaired, SC_BB );
		if(/*bbsc.score() < E_TH &&*/ ( et_scsr < E_TH || hb_sc < E_TH ) ){
			bbsc.reverse_in_place_unsafe();
			motifs.push_back(bbsc);
			if(option[mh::harvest::dump]()){
				if(header){ bbsc.print_header(cout); header=false; }
				cout << bbsc << endl;
				Real rms = bbsc.dump_aligned_motif(DBGDIR+bbsc.tag()+".pdb",pose,jr,ir);
				if(rms > 0.01) utility_exit_with_message("bad rms "+string_of(rms)+" "+DBGDIR+bbsc.tag()+".pdb");
			}
		}
		ResPairMotif scbb( tag, pose, (~scr1*bbr2).rt6(), ir, jr, nb1, nb2, et_srsc, et_srbb, 0, hb_sc, 0, 0, bfir, bfjr, is_sspaired, SC_BB );
		if(/*scbb.score() < E_TH &&*/ ( et_srsc < E_TH || hb_sc < E_TH ) ){
			motifs.push_back(scbb);
			if(option[mh::harvest::dump]()){
				if(header){ scbb.print_header(cout); header=false; }
				cout << scbb << endl;
				Real rms = scbb.dump_aligned_motif(DBGDIR+scbb.tag()+".pdb",pose,ir,jr);
				if(rms > 0.01) utility_exit_with_message("bad rms "+string_of(rms)+" "+DBGDIR+scbb.tag()+".pdb");
			}
		}
		ResPairMotif scsc( tag, pose, (~scr1*scr2).rt6(), ir, jr, nb1, nb2, et_srsr, 0, 0, hb_sc, 0, 0, bfir, bfjr, is_sspaired, SC_SC );
		if(/*scsc.score() < E_TH  &&*/ ( et_srsr < E_TH || hb_sc < E_TH ) ){
			motifs.push_back(scsc);
			if(option[mh::harvest::dump]()){
				if(header){ scsc.print_header(cout); header=false; }
				cout << scsc << endl;
				Real rms = scsc.dump_aligned_motif(DBGDIR+scsc.tag()+".pdb",pose,ir,jr);
				if(rms > 0.01) utility_exit_with_message("bad rms "+string_of(rms)+" "+DBGDIR+scsc.tag()+".pdb");
			}
		}

	return true;
 }

bool harvest_motifs_frag(
	ResPairMotifs & motifs,
	string const & tag0,
	Pose const & pose,
	core::scoring::dssp::Dssp & dssp,
	Size const & ir,
	Size const & jr,
	vector1<Real> const & /*occupancy*/,
	vector1<Real> const & bfactors,
	vector1<Real> const & /*rsd_sasa*/,
	vector1<Real> const & /*nbrs*/
 ){
	using namespace basic::options::OptionKeys;
	using core::scoring::hbonds::HBond;
	using namespace core::pose::motif;
	using core::scoring::motif::RPM_Type;
	if(pose.chain(ir)!=pose.chain(jr)) return false;
	bool const is_sspaired = dssp.in_paired_strands(ir,jr) ;
	std::string const tag = tag0+"frg";
	// Real const nb1 = nbrs[ir];
	// Real const nb2 = nbrs[jr];
	Real const & bfir(bfactors[ir]);
	Real const & bfjr(bfactors[jr]);
	Xform const bbr1 = core::pose::motif::get_backbone_reference_frame(pose,ir);
	Xform const bbr2 = core::pose::motif::get_backbone_reference_frame(pose,jr);
	if( bbr1.bad() || bbr2.bad() ) {
		cout << "bad bbr " << ir << " " << jr << endl;
		return false;
	}
	bool samess=true;
	for(Size i = ir; i <= jr; ++i){
		samess &= (dssp_reduced(pose.secstruct(ir))==dssp_reduced(pose.secstruct(i)));
	}
	Real6 rt6 = (~bbr1*bbr2).rt6();
	if( rt6[1] < -MOTIF_HASH_CART_SIZE || MOTIF_HASH_CART_SIZE < rt6[1] ) return false;
	if( rt6[2] < -MOTIF_HASH_CART_SIZE || MOTIF_HASH_CART_SIZE < rt6[2] ) return false;
	if( rt6[3] < -MOTIF_HASH_CART_SIZE || MOTIF_HASH_CART_SIZE < rt6[3] ) return false;
	ResPairMotif frag( tag, pose, rt6, ir, jr, jr-ir+1, samess?1.0:0.0, 0, 0, 0, 0, 0, 0, bfir, bfjr, is_sspaired, BB_BB );
	motifs.push_back(frag);
	if(option[mh::harvest::dump]()) cout << frag << endl;
	return true;
 }

Real harvest_time_refine(Pose & pose, core::scoring::ScoreFunction const & sf, std::string const & fname){
	using namespace protocols;
	using namespace core;
	using namespace core::scoring;
	// for(Size ir = 1; ir <= pose.size(); ++ir){
	// 	for(Size ia = 5; ia <= pose.residue(ir).nheavyatoms(); ++ia){
	// 		cout << pose.residue(ir).atom_name(ia) << endl;
	// 		cout << pose.dof(core::id::DOF_ID(AtomID(ia,ir),core::id::PHI  )) << " "
	// 		     << pose.dof(core::id::DOF_ID(AtomID(ia,ir),core::id::THETA)) << " "
	// 		     << pose.dof(core::id::DOF_ID(AtomID(ia,ir),core::id::D    )) << endl;
	// 		cout << pose.residue_type(ir).icoor(ia) << endl;
	// 	}
	// }
	Pose orig(pose);
	protocols::simple_moves::AddConstraintsToCurrentConformationMover tether;
	tether.bb_only() = false;
	tether.CA_only() = false;
	tether.use_distance_cst() = false;
	tether.coord_dev  () = 0.8;
	tether.bound_width() = 0.0;
	tether.apply(pose);
	ScoreFunctionOP sfcart = sf.clone();
	sfcart->set_weight(coordinate_constraint,1.0);
	// sfcart->set_weight(pro_close,10.0);
	// sfcart->set_weight(fa_rep,0.8);
	core::kinematics::MoveMapOP mmap = new core::kinematics::MoveMap;
	mmap->set_chi(true);
	mmap->set_bb(true);
	mmap->set_jump(true);
	// pose.dump_pdb("test_0.00.pdb");
	for(Real cart_wt = 0.1; cart_wt < 13.0; cart_wt *= 2.0){
		cout << "cart_min " << cart_wt << endl;
		sfcart->set_weight(core::scoring::cart_bonded,cart_wt);
		simple_moves::MinMoverOP minMover = new simple_moves::MinMover( mmap, sfcart, "lbfgs_armijo", 0.01, true );
		minMover->cartesian(true);
		minMover->apply(pose);
		sfcart->show(pose);
		// pose.dump_pdb("test_"+F(4,1,cart_wt)+".pdb");
	}
	protocols::idealize::IdealizeMover().apply(pose);
	{
		sfcart->set_weight(core::scoring::cart_bonded,0.0);
		simple_moves::MinMoverOP minMover = new simple_moves::MinMover( mmap, sfcart, "lbfgs_armijo", 0.01, true );
		minMover->cartesian(false);
		minMover->apply(pose);
		pose.dump_pdb(option[out::file::o]()+"/"+utility::file_basename(fname)+"_ideal.pdb");
	}

	cout << "ORIG:" << endl;
	sf.show(orig);
	cout << "IDEALIZED:" << endl;
	sf.show(pose);
	Real rms = core::scoring::CA_rmsd(pose,orig);
	Real allrms = core::scoring::all_atom_rmsd(pose,orig);
	cout << "HARVEST_REFINE_RMS " << rms << " " << allrms << endl;

	// for(Size ir = 1; ir <= pose.size(); ++ir){
	// 	for(Size ia = 5; ia <= pose.residue(ir).nheavyatoms(); ++ia){
	// 		cout << pose.residue(ir).atom_name(ia) << endl;
	// 		cout << pose.dof(core::id::DOF_ID(AtomID(ia,ir),core::id::PHI  )) << " "
	// 		     << pose.dof(core::id::DOF_ID(AtomID(ia,ir),core::id::THETA)) << " "
	// 		     << pose.dof(core::id::DOF_ID(AtomID(ia,ir),core::id::D    )) << endl;
	// 		cout << pose.residue_type(ir).icoor(ia) << endl;
	// 	}
	// }

	return allrms;
 }

void harvest_motifs(){
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace ObjexxFCL::format;

	core::chemical::AtomTypeSetCAP atypeset = core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
	core::scoring::etable::EtableOptions eopt;
	core::scoring::etable::Etable etable(atypeset,eopt);

	bool allonefile = false;
	std::string motif_out_file = basic::options::option[basic::options::OptionKeys::mh::motif_out_file]()+".rpm.bin.gz";
	utility::io::ozstream allout;
	if(basic::options::option[basic::options::OptionKeys::mh::motif_out_file].user()){
		if(utility::file::file_exists(motif_out_file)) utility_exit_with_message("already done!");
		allonefile = true;
		allout.open(motif_out_file,std::ios::out|std::ios::binary);
		if(!allout.good()){
			allout.close();
			utility_exit_with_message("error opening "+motif_out_file );
		}
	}

	vector1<std::string> const & fnames( option[mh::harvest_motifs]() );
	for(int ifile = 1; ifile <= (int)fnames.size(); ++ifile ){
		ResPairMotifs motifs;
		std::string fname = fnames[ifile];
		std::string const tag = tag_from_pdb_fname(fname);
		if(""==tag){
			TR.Error << "skipping non-pdb fname: " << fname << endl;
			continue;
		}

		if(!allonefile && utility::file::file_exists(option[out::file::o]()+"/"+tag.substr(1,2)+"/"+tag+".rpm.bin.gz")) continue;

		Pose pose; vector1<Real> bfactors,occupancy;
		utility::vector1<int> pdbres;
		std::map<int,char> pdbchain;
		int nresmodel1(0);
		if( !protocols::sic_dock::read_biounit(fname,pose,bfactors,occupancy,pdbres,pdbchain,nresmodel1,option[mh::harvest::max_res](),false) ){
			TR.Error << "FAIL TO READ " << fname << endl;
			continue;
		}
		// core::chemical::ResidueType const & gly_type( (pose.residue(1).residue_type_set().name_map("GLY") ));
		// Pose pose_no_ir(pose), pose_no_jr(pose), pose_no_ir_jr(pose);
		ScoreFunctionOP sf = core::scoring::get_score_function();
		core::scoring::methods::EnergyMethodOptions myopt = sf->energy_method_options();
		myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
		sf->set_energy_method_options(myopt);

		if( (int)option[mh::path::biounit_ideal]().size() >= ifile ){
			std::string fcache = option[mh::path::biounit_ideal]()[ifile];
			runtime_assert( utility::file::file_exists(fcache) );
			cout << "using cached ideal coods from " << fcache << endl;
			Pose tmp;
			core::import_pose::pose_from_file(tmp,fcache, core::import_pose::PDB_file);
			if( tmp.sequence() != pose.sequence() ) utility_exit_with_message("cached ideal coordts, sequence doesn't match");
			pose = tmp;
		} else if( option[mh::harvest::idealize]() ) {
			Real rms = harvest_time_refine(pose,*sf,fname);
			if(option[mh::harvest::max_rmsd]() < rms) {
				cout << "RMS_FAIL " << fname << endl;
				continue;
			}
		}
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
		for(int ir = 1+min_boundary; ir <= nresmodel1-min_boundary; ++ir){
			if(!pose.residue(ir).is_protein()) { ++nbad; continue; }
			for ( utility::graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node(ir)->const_upper_edge_list_begin(),
					irue = energy_graph.get_node(ir)->const_upper_edge_list_end();
					iru != irue; ++iru
			){
				EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
				harvest_motifs_one(
					motifs,
					tag,
					pose,
					dssp,
					hbset,
					etable,
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
			// for(Size jr=ir+1; jr <= min(nresmodel1-min_boundary,ir+10); ++jr){
			// 	harvest_motifs_frag(
			// 		motifs,
			// 		tag,
			// 		pose,
			// 		dssp,
			// 		ir,
			// 		jr,
			// 		occupancy,
			// 		bfactors,
			// 		rsd_sasa,
			// 		nbrs
			// 	);
			// }
		}


		// do I/O
		if(option[mh::filter::filter_harvest]()) filter_motifs(motifs);
		TR << fname <<' '<< motifs.size() <<" nbad: "<< I(5,nbad) <<" nbadbfac: "<< I(5,nbadbfac) <<" nbadocc: "<< I(5,nbadocc)
		 <<" nbadiface: "<< I(5,nbadiface) <<" nbaddist: "<< I(5,nbaddist) <<" nbadsselem: "<<I(5,nbadsselem)<< endl;
		if(allonefile){
			TR << "appending " << motifs.size() << " motifs to " << motif_out_file+".rpm.bin.gz" << endl;
			if(!write_motifs_binary(allout,motifs)) utility_exit_with_message("error writing to file "+motif_out_file+".rpm.bin.gz");
		} else {
			std::string outfile = option[out::file::o]()+"/"+tag.substr(1,2)+"/"+tag+".rpm.bin.gz";
			TR << "dumping " << motifs.size() << " motifs to " << outfile << endl;
			write_motifs_binary(outfile,motifs);
			if(!utility::file::file_exists(outfile)) utility_exit_with_message("couldn't make file "+outfile);
		}

	}
	// DONE: ;
	if(allonefile){
		if(!allout.good()){
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

	if(64!=sizeof(ResPairMotif)) utility_exit_with_message("ResPairMotif should be 64 bytes in size!!!!!!");

	// tmp_test_euler();
	// utility_exit_with_message("EULER");

	if(option[mh::dump_input_pdb].user()) {
		TR << "dumping input_pdbs from  -mh:dump_input_pdb" << endl;
		vector1<string> const & fnames( option[mh::dump_input_pdb]() );
		for(int ifile = 1; ifile <= (int)fnames.size(); ++ifile ){
			Pose pose; vector1<Real> bfactors,occupancy;
			utility::vector1<int> pdbres;
			std::map<int,char> pdbchain;
			int nresmodel1;
			if( !protocols::sic_dock::read_biounit(fnames[ifile],pose,bfactors,occupancy,pdbres,pdbchain,nresmodel1,option[mh::harvest::max_res](),true) ){
				TR.Error << "FAIL TO READ " << fnames[ifile] << endl;
				continue;
			}
		}
	}
	// if(option[mh::print_motifs_boost].user()){
	// 	TR << "reading motifs in files following -mh:print_motifs_boost" << endl;
	// 	print_motifs_boost(TR);
	// }
	if(option[mh::harvest_motifs].user()){
		TR << "harvesting motifs in files following -mh:harvest_motifs" << endl;
		harvest_motifs();
	}
	if(option[mh::print_motifs].user()){
		TR << "reading motifs in files following -mh:print_motifs" << endl;
		ResPairMotif::print_header(cout);
		print_motifs(cout);
	}
	if(option[mh::merge_motifs].user()){
		TR << "reading motifs in files following -mh:merge_motifs" << endl;
		merge_motifs();
	}
	if(option[mh::merge_scores].user()){
		TR << "reading scores in files following -mh:merge_scores" << endl;
		merge_scores();
	}
	if(option[mh::dump_motif_pdbs].user()){
		TR << "reading motifs in files following -mh:dump_motif_pdbs" << endl;
		dump_motif_pdbs();
	}
	if(option[mh::harvest_scores].user()){
		TR << "reading motif counts in files following -mh:harvest_scores" << endl;
		harvest_scores();
	}
	if(option[mh::print_scores].user()){
		TR << "reading motifs in files following -mh:print_scores" << endl;
		print_scores();
	}
	if(option[mh::remove_duplicates].user()){
		TR << "reading motifs in files following -mh:remove_duplicates" << endl;
		remove_duplicates();
	}
	// if(option[mh::score_pdbs].user()){
	// 	if(!option[mh::xform_score_data].user()) utility_exit_with_message("specify -mh:score_data <datafile>.xh.bin.gz");
	// 	TR << "reading pdbs to score following -mh:score_data " << option[mh::xform_score_data]() << endl;
	// 	TR << "reading scoring data from file " << option[mh::xform_score_data]() << " following -mh:score_data" << endl;
	// 	score_with_counts();
	// }
	// if(option[mh::sequence_recovery].user()){
	// 	if(!option[mh::path::motifs].user()) utility_exit_with_message("specify -mh:input_motifs <datafile>.xh.bin.gz");
	// 	TR << "reading pdbs to score following -mh:input_motifs " << option[mh::path::motifs]() << endl;
	// 	TR << "reading scoring data from file " << option[mh::path::motifs]() << " following -mh:input_motifs" << endl;
	// 	sequence_recovery();
	// }
	if(option[mh::dump_matching_motifs].user()){
		// if(!option[mh::path::motifs].user()) utility_exit_with_message("specify -mh:input_motifs <datafile>.xh.bin.gz");
		// TR << "reading pdbs to score following -mh:input_motifs " << option[mh::path::motifs]() << endl;
		// TR << "reading scoring data from file " << option[mh::path::motifs]() << " following -mh:input_motifs" << endl;
		dump_matching_motifs();
	}

	// {
	// 	utility::io::ozstream ozs("test.mfh.gz");
	// 	boost::archive::binary_oarchive out_archive(ozs);
	// 	out_archive << mh;
	// }

	// utility::io::ozstream ofs1("test1.txt"); ofs1 << mh << endl; ofs1.close();
	// MotifHash mh2;
	// utility::io::izstream ifs1("test1.txt"); ifs1 >> mh2        ; ifs1.close();
	// utility::io::ozstream ofs2("test2.txt"); ofs2 << mh2 << endl; ofs2.close();

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}


