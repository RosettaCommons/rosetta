// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// includes
	#ifndef INCLUDED_core_scoring_motif_motif_hash_stuff_hh
	#define INCLUDED_core_scoring_motif_motif_hash_stuff_hh

	#include <core/scoring/motif/motif_hash_stuff.fwd.hh>

	#include <core/pose/Pose.fwd.hh>
	#include <core/pose/util.hh>
	#include <numeric/geometry/hashing/xyzStripeHash.fwd.hh>
	#include <numeric/xyzVector.hh>
	#include <utility/fixedsizearray1.hh>
	#include <numeric/geometry/hashing/SixDHasher.hh>
	#include <numeric/xyzTransform.hh>
	#include <numeric/HomogeneousTransform.hh>
	#include <boost/unordered_map.hpp>

namespace core {
namespace scoring {
namespace motif {

	extern const numeric::Real MOTIF_HASH_CART_SIZE;
	extern const uint64_t MOTIF_VERSION;

//types

	using core::id::AtomID;
	using core::pose::Pose;
	using core::pose::PoseCOP;
	using core::pose::PoseCOP;
	using core::Real;
	using core::scoring::ScoreFunctionOP;
	using core::Size;
	using std::string;
	using utility::vector1;
	using numeric::geometry::hashing::Real3;
	using numeric::geometry::hashing::Real6;
	// using numeric::geometry::hashing::xyzStripeHash;
	// using numeric::geometry::hashing::xyzStripeHashCOP;
	// using numeric::geometry::hashing::xyzStripeHashCOP;

	class ResPairMotif;
	class MotifHit;
	class MotifHits;
	class MotifRotamerSetOperation;
	typedef utility::pointer::shared_ptr<MotifRotamerSetOperation      > MotifRotamerSetOperationOP;
	typedef utility::pointer::shared_ptr<MotifRotamerSetOperation const> MotifRotamerSetOperationCOP;
	typedef utility::vector1<Real> Reals;
	typedef utility::vector1<Size> Sizes;
	typedef utility::vector1<int>  Ints;
	typedef utility::vector1<float> Floats;
	typedef utility::vector1<bool> Bools;
	typedef numeric::xyzVector<core::Real> Vec;
	typedef numeric::xyzMatrix<core::Real> Mat;
	typedef numeric::xyzTransform<core::Real> Xform;

	static const numeric::Xforms XsI = numeric::Xforms(1,numeric::Xform::identity());


// XformScore
	class XformScore  : public utility::pointer::ReferenceCount {
	 public:
		typedef uint32_t Key;
		typedef float    Score;
		typedef numeric:: geometry::hashing::SixDCoordinateBinner SixDCoordinateBinner;
		typedef numeric::geometry::hashing::bin_index_hasher bin_index_hasher;
		SixDCoordinateBinner hasher_;
		typedef boost::unordered_map< Key, Score, bin_index_hasher > ScoreMap;
		ScoreMap scores_;
		double cart_size_,cart_resl_,angle_resl_;
	 public:
		XformScore(Real cart_resl, Real ang_resl);
		// Key key_of(Real6 const & xform);
		void set_score(Key const &, Score const & count);
		void add_score(Key const &, Score const & count);
		void max_score(Key const &, Score const & count);
		Score score_of_bin(Real6 const & xform) const;
		Score score_of_bin(Xform const & xform) const;
		Score score_of_bin(Key const & key) const;
		Real  score_of_bin_normalized(Real6 const & xform) const;
		bool write_binary(std::ostream & out) const;
		bool write_binary(string const & fname) const;
		static bool read_binary(XformScoreOP & xs, std::istream & in, bool clearme=true, Real const & wt=1.0);
		static bool read_binary(XformScoreOP & xs, string const & fname, bool clearme=true,Real const & wt=1.0);
		static bool read_binary(XformScoreOP & xs, utility::vector1<std::string> const & fnames,Real const & wt=1.0);
		Size num_bins() const;
		Size num_possible_bins() const;
		Score tot_score() const;
		void print_scores(std::ostream & out, std::string const tag="XH ") const;
		void print_scores_norm(std::ostream & out, std::string const tag="XH ") const;
		void clear();
		Real cart_size() const { return cart_size_; }
		Key bin_index(Real6 const & rt) const;
		void print_quantiles(std::ostream & out, int num) const;
		void prune_small_bins(Score thresh);
		void sanity_check() const;
		friend std::ostream & operator << (std::ostream & out, XformScore const & x);
		void aggregate_add(XformScore const & other);
		void aggregate_max(XformScore const & other);
	private:
		XformScore();
	};
	typedef std::map<std::string,XformScore> XformScoreSringMap;
	std::ostream & operator<<(std::ostream & out, XformScore const & x);
	std::ostream & operator<<(std::ostream & out, XformScoreSringMap const & x);


	struct RPM_FilterStats {
		Size F_aa1,                P_aa1;
		Size F_aa2,                P_aa2;
		Size F_coorderr,           P_coorderr;
		Size F_uniformfrag,        P_uniformfrag;
		Size F_dssp1,              P_dssp1;
		Size F_dssp2,              P_dssp2;
		Size F_score,              P_score;
		Size F_faatr,              P_faatr;
		Size F_faatr_or_hb,        P_faatr_or_hb;
		Size F_faatr_or_hbbb,      P_faatr_or_hbbb;
		Size F_hb_bb,              P_hb_bb;
		Size F_hb_bb_sc,           P_hb_bb_sc;
		Size F_hb_sc,              P_hb_sc;
		Size F_max_seqsep,         P_max_seqsep;
		Size F_maxdist2,           P_maxdist2;
		Size F_mindist2,           P_mindist2;
		Size F_type_SC_SC,         P_type_SC_SC;
		Size F_type_SC_BB,         P_type_SC_BB;
		Size F_type_SC_PH,         P_type_SC_PH;
		Size F_type_SC_PO,         P_type_SC_PO;
		Size F_type_BB_BB,         P_type_BB_BB;
		Size F_type_BB_PH,         P_type_BB_PH;
		Size F_type_BB_PO,         P_type_BB_PO;
		Size F_type_PH_PO,         P_type_PH_PO;
		Size F_no_hb_bb,           P_no_hb_bb;
		Size F_nodisulf,           P_nodisulf;
		Size F_noloops,            P_noloops;
		Size F_not_restype,        P_not_restype;
		Size F_not_restype_one,    P_not_restype_one;
		Size F_occupancy,          P_occupancy;
		Size F_oneloop,            P_oneloop;
		Size F_pdb,                P_pdb;
		Size F_lig,                P_lig;
		Size F_restype1,           P_restype1;
		Size F_restype2,           P_restype2;
		Size F_restype,            P_restype;
		Size F_restype_one,        P_restype_one;
		Size F_seqsep,             P_seqsep;
		Size F_ss1,                P_ss1;
		Size F_ss2,                P_ss2;
		friend std::ostream & operator << (std::ostream & out, RPM_FilterStats const & s);
		RPM_FilterStats();
	};



class ResPairMotif {
 private:
	/*  4 */ char ss1_,ss2_,aa1_,aa2_;
	/* 12 */ utility::fixedsizearray1<uint16_t,6> rt6_;
	/*  8 */ utility::fixedsizearray1<uint8_t,4> chi1_,chi2_;
			 uint8_t fa_atr_;
			 uint8_t fa_atr_sc_bb_;
			 uint8_t fa_atr_bb_;
			 uint8_t hb_sc_;
			 uint8_t hb_bb_sc_;
			 uint8_t hb_bb_;
			 uint8_t nbrs1_;
			 uint8_t nbrs2_;
	/*	1 */ uint8_t count_;
	/*  1 */ uint8_t type_;
	/*  2 */ uint8_t misc1_;
	/*  6 */ uint16_t lg1x_,lg1y_,lg1z_; // dummy for lig / metal / water
	/*  6 */ uint16_t lg2x_,lg2y_,lg2z_; // euler ang for lig?
	/*  2 */ uint8_t  bfac1_,bfac2_;
	/*  2 */ char chain1_,chain2_;
	/*  4 */ uint16_t resi1_,resi2_;
	/*  8 */ utility::fixedsizearray1<unsigned char,8> pdb_; // simplify this // last3 are lig resn
 public:
		ResPairMotif();
		ResPairMotif(
			string const & tag,
			Pose const & pose,
			Real6 const & rt6,
			Size const & resi1,
			Size const & resi2,
			Real const & nbrs1,
			Real const & nbrs2,
			Real const & fa_atr,
			Real const & fa_atr_sc_bb,
			Real const & fa_atr_bb,
			Real const & hb_sc,
			Real const & hb_bb_sc,
			Real const & hb_bb,
			Real const & bfac1,
			Real const & bfac2,
			bool is_sspaired,
			RPM_Type _type
		);
		bool is_reversible() const                                      { return type1()==type2(); }
		ResPairMotif & reverse_in_place();
		ResPairMotif & reverse_in_place_unsafe();
		ResPairMotif   reversed() const;
		void reset();
		Xform xform() const;
		void  xform(Xform const & _xform );
		Real6 rt() const;
		void  rt(Real6 const & rt    );
		Real dist2       () const;
		Real bfac1       () const;
		Real bfac2       () const;
		Real fa_atr      () const;
		Real fa_atr_sc_bb() const;
		Real fa_atr_bb   () const;
		Real hb_sc       () const;
		Real hb_bb_sc    () const;
		Real hb_bb       () const;
		Real nbrs1       () const;
		Real nbrs2       () const;
		Real chi11() const;
		Real chi12() const;
		Real chi13() const;
		Real chi14() const;
		Real chi21() const;
		Real chi22() const;
		Real chi23() const;
		Real chi24() const;
		Real chi1(Size const & ichi) const;
		Real chi2(Size const & ichi) const;
		std::string pdb()  const                                       {
            std::string pdb_r(8,' ');
            for(Size ii=0; ii<8; ++ii){
                pdb_r[ii]=pdb_[ii+1];
            }
            return(pdb_r);   }
		bool check_pdb_code(std::string const & qpdb) const            { return pdb_[1]==qpdb[0] && pdb_[2]==qpdb[1] && pdb_[3]==qpdb[2] && pdb_[4]==qpdb[3] && (qpdb.size()==4||pdb_[5]==qpdb[4]); }
		bool check_lig_code(std::string const & qlig) const            { return pdb_[6]==qlig[0] && pdb_[7]==qlig[1] && pdb_[8]==qlig[2]; }
		void bfac1       (Real const & bfac1   );
		void bfac2       (Real const & bfac2   );
		void fa_atr      (Real const & fa_atr  );
		void fa_atr_sc_bb(Real const & fa_atr_sc_bb  );
		void fa_atr_bb   (Real const & fa_atr_bb  );
		void hb_sc       (Real const & hb_sc   );
		void hb_bb_sc    (Real const & hb_bb_sc);
		void hb_bb       (Real const & hb_bb   );
		void nbrs1       (Real const & nbrs    );
		void nbrs2       (Real const & nbrs    );
		char   aa1() const;
		char   aa2() const;
		char   ss1() const;
		char   ss2() const;
		char dssp1() const;
		char dssp2() const;
		Real score() const;
		RPM_Type type () const;
		RM_Type  type1() const;
		RM_Type  type2() const;
		void     type(RPM_Type const & type)                           { type_ = (char)type; }
		// char sscode() const;
		inline Size resi1() const                                      { return resi1_; }
		inline Size resi2() const                                      { return resi2_; }
		inline Size count() const                                      { return count_; }
		inline void count(Size const & count)                                { count_ = count; }
		inline void addcount()                                         { if(count_<255) ++count_; }
		std::string tag() const;
		bool filter(RPM_FilterStats *stat=NULL) const;

		/*TODO;*/
		void dump_pdb(std::ostream & out, numeric::xyzTransform<Real> const & xform=Xform::identity(), std::string tag="1") const;
		void dump_pdb(string const & fname="", string const & tag="1") const;
		Real dump_aligned_motif(std::ostream     & out, Pose const & pose1, Size const & ir, Pose const & pose2, Size const & jr, Size & aotmno, int const & num=1, numeric::Xforms const & xforms=XsI) const;
		Real dump_aligned_motif(std::ostream     & out, Pose const & pose , Size const & ir              , Size const & jr, Size & aotmno, int const & num=1, numeric::Xforms const & xforms=XsI) const;
		Real dump_aligned_motif(std::string const & fn, Pose const & pose1, Size const & ir, Pose const & pose2, Size const & jr, Size & aotmno, int const & num=1, numeric::Xforms const & xforms=XsI) const;
		Real dump_aligned_motif(std::string const & fn, Pose const & pose , Size const & ir              , Size const & jr, Size & aotmno, int const & num=1, numeric::Xforms const & xforms=XsI) const;
		Real dump_aligned_motif(std::ostream     & out, Pose const & pose1, Size const & ir, Pose const & pose2, Size const & jr, int const & num=1, numeric::Xforms const & xforms=XsI) const;
		Real dump_aligned_motif(std::ostream     & out, Pose const & pose , Size const & ir              , Size const & jr, int const & num=1, numeric::Xforms const & xforms=XsI) const;
		Real dump_aligned_motif(std::string const & fn, Pose const & pose1, Size const & ir, Pose const & pose2, Size const & jr, int const & num=1, numeric::Xforms const & xforms=XsI) const;
		Real dump_aligned_motif(std::string const & fn, Pose const & pose , Size const & ir              , Size const & jr, int const & num=1, numeric::Xforms const & xforms=XsI) const;
		void fill_pose_with_motif(Pose & pose, int const & ir=1, int const & jr=2) const;
		core::pose::Pose get_pose(int const & ir=1, int const & jr=2) const;
		// io
		friend std::ostream & operator << (std::ostream & out, ResPairMotif const & x);
		// friend std::istream & operator >> (std::istream & in , ResPairMotif & x);
		// friend class boost::serialization::access;
	 //    template<class Archive> void serialize(Archive & ar, const unsigned int version);
		static void print_header(std::ostream & out);
 		bool operator < (ResPairMotif const & other) const;
		bool operator ==(ResPairMotif const & other) const;

	};
	std::ostream & operator<<(std::ostream & out, RM_Type const & x);
	std::ostream & operator<<(std::ostream & out, RPM_Type const & x);
	std::ostream & operator << (std::ostream & out, ResPairMotif const & x);
	std::istream & operator >> (std::istream & in , ResPairMotif       & x);
	typedef utility::vector1<ResPairMotif const*> ResPairMotifPtrs;

class ResPairMotifMetaBinner {
 public:
	typedef uint64_t Key;
	ResPairMotifMetaBinner();
	virtual ~ResPairMotifMetaBinner(){}
	std::string motif_bin_label(ResPairMotif const & r) const;
	Key         motif_bin_hash (ResPairMotif const & r) const;
	std::string hash_to_labal(Key const & k) const;
	Key	bin0_of_real(Real val, Reals const & breaks) const;
	Reals const sep_lj_,  sep_hb_,  sep_nbrs_,  sep_bfac_,  sep_dist_;
	bool  const sep_aa_,  sep_aa1_, sep_aa2_,  sep_ss_,  sep_dssp_;
 };
 typedef boost::unordered_map<ResPairMotifMetaBinner::Key,XformScoreOP> XformScoreMap;
 std::ostream & operator<<(std::ostream & out, XformScoreMap const & x);

class ResPairMotifs : public utility::vector1<ResPairMotif> {
 public:
	void filter_structurally_identical_motifs();
	void filter_structurally_similar_motifs(Real cart_size, Real cart_resl, Real ang_resl);
	void add_reverse_motifs();
 };
 typedef std::map<std::string,ResPairMotifs> ResPairMotifsStringMap;
 typedef boost::unordered_map<ResPairMotifMetaBinner::Key,ResPairMotifs> ResPairMotifsMap;

bool read_motifs_binary (               string const & fname , ResPairMotifs       & motifs);
bool read_motifs_binary (  std::istream        & in    , ResPairMotifs       & motifs);
bool read_motifs_binary (vector1<string> const & fnames, ResPairMotifs       & motifs);
bool write_motifs_binary(  std::ostream        & out   , ResPairMotifs const & motifs);
bool write_motifs_binary(        string const &        fname , ResPairMotifs const & motifs);

void filter_motifs(	ResPairMotifs const & motifs_in, ResPairMotifs & motifs_out );
void filter_motifs( ResPairMotifs & motifs );

class MotifHit {
  public:
	core::pose::PoseCOP pose1_,pose2_;
	core::pose::PoseOP mpose_;
	int residue1,residue2,index;
	float score1,score2,rms,chi_rmsd_;
	ResPairMotif motif;
	MotifHit():residue1(0),residue2(0),score1(0),score2(0),rms(9e9),chi_rmsd_(9e9) {}
	MotifHit(PoseCOP _pose1, PoseCOP _pose2, int res1, int res2, ResPairMotif const & sm)
	 : pose1_(_pose1),pose2_(_pose2),mpose_(core::pose::PoseOP( new Pose )),residue1(res1),residue2(res2),score1(0),score2(0),rms(9e9),chi_rmsd_(9e9),motif(sm) {}
	bool operator< (MotifHit const & other) const;
	bool operator==(MotifHit const & other) const;
	core::pose::Pose       & mpose()       { return *mpose_; }
	core::pose::Pose const & mpose() const { return *mpose_; }
	core::pose::Pose const & pose1() const { return *pose1_; }
	core::pose::Pose const & pose2() const { return *pose2_; }
	core::pose::PoseOP  mposeptr()       { return mpose_; }
	core::pose::PoseCOP mposeptr() const { return mpose_; }
	core::pose::PoseCOP pose1ptr() const { return pose1_; }
	core::pose::PoseCOP pose2ptr() const { return pose2_; }
	core::Real chi_rmsd() const { return chi_rmsd_; }
 };
 std::ostream & operator<<(std::ostream & out, MotifHit  const & h);
 class MotifHitStats : public std::map<std::string,Real> {
  };

class MotifHits : public utility::vector1<MotifHit> {
  public:
	mutable int rescover,paircover;
	mutable Real total_score;
	typedef boost::unordered_multimap<uint32_t,uint32_t> ResMap;
	ResMap resmap1_,resmap2_;
	void compute_metrics() const;
	void filter_redundant();
	void sort_by_rms_and_energy();
	void dump_motifs_pdb( std::ostream & out, int & count,   Size & atomno, numeric::Xforms const & xforms=XsI ) const;
	void dump_motifs_pdb( string const & fname    ,                Size & atomno, numeric::Xforms const & xforms=XsI ) const;
	void dump_motifs_pdb( string const & fname, Pose const & pose, Size & atomno, numeric::Xforms const & xforms=XsI ) const;
	void dump_motifs_pdb( std::ostream & out, int & count,   numeric::Xforms const & xforms=XsI ) const { Size tmp=0; dump_motifs_pdb(out,count,tmp,xforms); }
	void dump_motifs_pdb( string const & fname    ,                numeric::Xforms const & xforms=XsI ) const { Size tmp=0; dump_motifs_pdb(fname,tmp,xforms); }
	void dump_motifs_pdb( string const & fname, Pose const & pose, numeric::Xforms const & xforms=XsI ) const { Size tmp=0; dump_motifs_pdb(fname,pose,tmp,xforms); }
	void dump_motifs_pdb( std::ostream & out, numeric::Xforms const & xforms=XsI ) const { int d(0); Size tmp=0; dump_motifs_pdb(out,d,tmp,xforms); }
//	int print_motifs( std::ostream & out) const; commented out because failing in the python build
	int stat_motifs( MotifHitStats & stat) const;
	string get_resfile(bool restrict=false) const;
	string get_resfile(bool restrict, std::set<core::Size> & resi_in_resfile) const;
	void dump_resfile(std::string fn) const;
	void push_back(MotifHit const & h);
	void push_back_raw(MotifHit const & h);
	Size num_hits1(Size ir) const;
	Size num_hits2(Size ir) const;
	Size num_hits(Size ir) const;
	ResPairMotifPtrs hits1(Size ir) const;
	ResPairMotifPtrs hits2(Size ir) const;
	// utility::vector1<core::conformation::ResidueOP> rotamers(Size const & ir) { return rotamers_[ir]; }
	friend std::ostream & operator<<(std::ostream & out, MotifHits const & h);
  };

enum MotifMatchOverlapReq {
	NO_OVERLAP=1,
	SMALL_OVERLAP, //gly ala
	MATCH_ROTAMER,
	MATCH_AA,
	MATCH_HYDROPHOBICITY,
	MATCH_ANY
 };
class ResPairMotifQuery {
	friend class MotifHash;
	PoseCOP pose1_,pose2_;
	numeric::geometry::hashing::xyzStripeHash const *inter_clash1_;
	numeric::geometry::hashing::xyzStripeHash const *inter_clash2_;
	Xform inter_clash1_xform_,inter_clash2_xform_;
	MotifMatchOverlapReq overlap1_,overlap2_;
	bool auto_clash_,interface_only_,match_ss1_,match_aa1_,match_ss2_,match_aa2_,clash_check_;
	Real match_radius_,max_ca_dis2_,match_chi_rmsd_;
	Bools useres1_,useres2_,useres1_ph_,useres2_ph_,useres1_po_,useres2_po_;
	void init(Pose const & pose1, Pose const & pose2);
	ResPairMotifQuery(){}
 public:
	ResPairMotifQuery( Pose const & pose1                     ) : clash_check_(true) { init(pose1,pose1); }
	ResPairMotifQuery( Pose const & pose1, Pose const & pose2 ) : clash_check_(true) { init(pose1,pose2); }
	PoseCOP pose1() const { return pose1_; }
	PoseCOP pose2() const { return pose2_; }
	// numeric::geometry::hashing::xyzStripeHash const * intra_clash1() const { return intra_clash1_; }
	// numeric::geometry::hashing::xyzStripeHash const * intra_clash2() const { return intra_clash2_; }
	numeric::geometry::hashing::xyzStripeHash const * inter_clash1() const { return inter_clash1_; }
	numeric::geometry::hashing::xyzStripeHash const * inter_clash2() const { return inter_clash2_; }
	// numeric::geometry::hashing::xyzStripeHash       * intra_clash1()       { return intra_clash1_; }
	// numeric::geometry::hashing::xyzStripeHash       * intra_clash2()       { return intra_clash2_; }
	numeric::geometry::hashing::xyzStripeHash const * & inter_clash1()       { return inter_clash1_; }
	numeric::geometry::hashing::xyzStripeHash const * & inter_clash2()       { return inter_clash2_; }
	Xform & inter_clash1_xform() { return inter_clash1_xform_; }
	Xform & inter_clash2_xform() { return inter_clash2_xform_; }
	Xform const & inter_clash1_xform() const { return inter_clash1_xform_; }
	Xform const & inter_clash2_xform() const { return inter_clash2_xform_; }
	bool  const & auto_clash() const { return auto_clash_; }
	bool        & auto_clash()       { return auto_clash_; }
	bool  const & interface_only() const { return interface_only_; }
	bool        & interface_only()       { return interface_only_; }
	Real  const & match_radius() const { return match_radius_; }
	Real        & match_radius()       { return match_radius_; }
	Real  const & match_chi_rmsd() const { return match_chi_rmsd_; }
	Real        & match_chi_rmsd()       { return match_chi_rmsd_; }
	Real  const & max_ca_dis2() const { return max_ca_dis2_; }
	Real        & max_ca_dis2()       { return max_ca_dis2_; }
	Bools const & useres1() const { return useres1_; }
	Bools       & useres1()       { return useres1_; }
	Bools const & useres2() const { return useres2_; }
	Bools       & useres2()       { return useres2_; }
	Bools const & useres1_ph() const { return useres1_ph_; }
	Bools       & useres1_ph()       { return useres1_ph_; }
	Bools const & useres2_ph() const { return useres2_ph_; }
	Bools       & useres2_ph()       { return useres2_ph_; }
	Bools const & useres1_po() const { return useres1_po_; }
	Bools       & useres1_po()       { return useres1_po_; }
	Bools const & useres2_po() const { return useres2_po_; }
	Bools       & useres2_po()       { return useres2_po_; }
	bool  const & match_ss1() const { return match_ss1_; }
	bool        & match_ss1()       { return match_ss1_; }
	bool  const & match_ss2() const { return match_ss2_; }
	bool        & match_ss2()       { return match_ss2_; }
	bool  const & match_aa1() const { return match_aa1_; }
	bool        & match_aa1()       { return match_aa1_; }
	bool  const & match_aa2() const { return match_aa2_; }
	bool        & match_aa2()       { return match_aa2_; }
	bool  const & clash_check() const { return clash_check_; }
	bool        & clash_check()       { return clash_check_; }
	MotifMatchOverlapReq const & overlap1() const { return overlap1_; }
	MotifMatchOverlapReq       & overlap1()       { return overlap1_; }
	MotifMatchOverlapReq const & overlap2() const { return overlap2_; }
	MotifMatchOverlapReq       & overlap2()       { return overlap2_; }
	friend std::ostream & operator<<(std::ostream & out, ResPairMotifQuery const & opt);
	void copy_opts_swapped(ResPairMotifQuery const & other);
 };
class MotifHash : public utility::pointer::ReferenceCount {
	public:
		typedef ResPairMotif Motif;
		typedef boost::uint64_t Key;
		typedef numeric::geometry::hashing::SixDCoordinateBinner SixDCoordinateBinner;
		typedef numeric::geometry::hashing::bin_index_hasher bin_index_hasher;
		typedef boost::unordered_multimap< Key, ResPairMotif , bin_index_hasher > MotifMap;
		typedef boost::unordered_set<Key> KeySet;
		typedef numeric::geometry::hashing::xyzStripeHashCOP ClashCheckerPtr;
		double const cart_size_,cart_resl_,angle_resl_;
		SixDCoordinateBinner hasher_;
		MotifMap motif_umap_;
		KeySet key_set_;
		RPM_Type type_;
	public:
		MotifHash();
		MotifHash(Real cart_resl, Real ang_resl);
		MotifHash(ResPairMotifs const & motifs);
		void sanity_check() const;
		SixDCoordinateBinner const & hasher() const { return hasher_; };
		SixDCoordinateBinner       & hasher()       { return hasher_; };
		RPM_Type type () const;
		RM_Type  type1() const;
		RM_Type  type2() const;
		bool check_bounds(Real6 const & rt6) const;
		void add_motif(Motif const & d);
		void add_motif(Motif const & d, Key const & key);
		void find_motifs(Key const & k, ResPairMotifs & results ) const;
		void find_motifs(Real6 const & rt, ResPairMotifs & results ) const;
		Size count_motifs(Real6 const & rt) const ;
		Size count_motifs(Key const & k) const;
		void find_motifs_with_radius(Real6 const & rt, Real radius, vector1<Motif> & results) const;
		Size num_motifs() { return motif_umap_.size(); }
		// Undefined, commenting out to fix PyRosetta build  Key key_of(Motif const & d) const;
		Real cart_size() const { return cart_size_; }
		Key bin_index(Real6 const & rt) const;
		Key bin_index(Motif const & m) const;
		// io
		int get_matching_motifs( ResPairMotifQuery const & opt, MotifHits & hits ) const;
		int get_matching_motifs( ResPairMotifQuery const & opt, MotifHits & hits, MotifHits & newhits ) const;
		MotifHits get_matching_motifs( ResPairMotifQuery const & opt ) const;
		friend std::ostream & operator << (std::ostream & out, MotifHash const & x);
		// friend std::istream & operator >> (std::istream & in , MotifHash & x);
	    // friend class boost::serialization::access;
	    // template<class Archive> void serialize(Archive & ar, const unsigned int version);
	    // KeySet::const_iterator begin() { return key_set_.begin(); }
	    // KeySet::const_iterator end  () { return key_set_.end  (); }
	    KeySet const & keys() const { return key_set_; }
	};
	std::ostream & operator << (std::ostream & out, MotifHash const & x);
	std::istream & operator >> (std::istream & in , MotifHash       & x);
	typedef std::map<std::string,MotifHashOP> MotifHashStringMap;


class MotifHashManager {
 public:
	typedef ResPairMotifMetaBinner::Key Key;
	~MotifHashManager();

	static MotifHashManager * get_instance();

	bool have_motifs() const { return motif_set_names_.size(); }
	MotifHashCOP get_motif_hash_BB_BB();
	MotifHashCOP get_motif_hash_SC_BB();
	MotifHashCOP get_motif_hash_SC_SC();
	MotifHashCOP get_motif_hash_BB_PH();
	MotifHashCOP get_motif_hash_BB_PO();
	void add_motif_set_name(std::string const & motifset);
	bool have_motif_set_named(std::string const & motifset) const;
	std::set<std::string> motif_set_names() const { return motif_set_names_; }
	MotifHashCOP get_motif_hash_by_fname(std::string const & fname);
	MotifHashStringMap const & get_motif_hash_by_fname();
	XformScoreCOP get_xform_score_SC_BB(char aa1);
	XformScoreCOP get_xform_score_SC_SC(char aa1, char aa2);
	XformScoreCOP get_xform_score_BB_BB(char ss1, char ss2, char aa1=' ', char aa2=' ');
	XformScoreCOP get_xform_score_frags();
	XformScoreCOP get_xform_score_BB_PH();
	XformScoreCOP get_xform_score_BB_PO();
	XformScoreCOP get_xform_score_PH_PO();

	int get_matching_motifs( ResPairMotifQuery const & opt, MotifHits & hits ) const;

 private:
	MotifHashManager();
	void init();
	static MotifHashManager * instance_;
	MotifHashOP motif_hash_SC_BB_;
	MotifHashOP motif_hash_SC_SC_;
	MotifHashOP motif_hash_BB_BB_;
	MotifHashOP motif_hash_BB_PH_;
	MotifHashOP motif_hash_BB_PO_;
	MotifHashStringMap motifs_by_fname_;
	XformScoreOP xform_score_BB_PH_;
	XformScoreOP xform_score_BB_PO_;
	XformScoreOP xform_score_PH_PO_;
	XformScoreMap xform_scores_BB_BB_;
	XformScoreOP xform_score_frags_;
	XformScoreMap xform_scores_SC_BB_;
	XformScoreMap xform_scores_SC_SC_;
	std::set<std::string> motif_set_names_;
	Key key_mask_BB_BB_/*,key_mask_SC_BB_,key_mask_frags_*/;
	bool done_loading_;
	friend void preload_motif_data(MotifHashManager & mman);

 };

 	void preload_motif_data(MotifHashManager & mman);

	void load_motifs(vector1<string> const & fnames, ResPairMotifs & motifs, ResPairMotifsStringMap * map=NULL);
	void load_motifs(string const & fname, ResPairMotifs & motifs, ResPairMotifsStringMap * map=NULL);





}
}
}

#endif












