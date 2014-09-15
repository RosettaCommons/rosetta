// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// includes
	#ifndef INCLUDED_protocols_motif_hash_motif_hash_stuff_hh
	#define INCLUDED_protocols_motif_hash_motif_hash_stuff_hh

	#include <protocols/motif_hash/motif_hash_stuff.fwd.hh>
	#include <ObjexxFCL/FArray2D.hh>
	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>
	#include <basic/Tracer.hh>
	#include <basic/database/open.hh>
	#include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/option_macros.hh>
	#include <core/chemical/AtomType.hh>
	#include <core/chemical/ChemicalManager.hh>
	#include <core/conformation/symmetry/util.hh>
	#include <core/import_pose/import_pose.hh>
	#include <core/io/silent/SilentFileData.hh>
	#include <core/pose/PDBInfo.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/util.hh>
	#include <core/pack/task/PackerTask.hh>
	#include <core/scoring/Energies.hh>
	#include <core/scoring/EnergyGraph.hh>
	#include <core/scoring/ScoreFunctionFactory.hh>
	#include <core/scoring/ScoreTypeManager.hh>
	#include <core/scoring/dssp/Dssp.hh>
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
	#include <protocols/sic_dock/RigidScore.hh>
	#include <protocols/sic_dock/SICFast.hh>
	#include <protocols/sic_dock/util.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/fixedsizearray1.hh>
	#include <core/pack/rotamer_set/RotamerSetOperation.hh>
	#include <numeric/geometry/hashing/SixDHasher.hh>
	#include <numeric/HomogeneousTransform.hh>

	#include <boost/unordered_map.hpp>

namespace protocols {
namespace motif_hash {


//types
	typedef numeric::xyzVector<core::Real> Vec;
	typedef numeric::xyzMatrix<core::Real> Mat;
	typedef utility::vector1<bool> Bools;

	using core::id::AtomID;
	using basic::options::option;
	using core::pose::Pose;
	using core::pose::PoseCOP;
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
	using ObjexxFCL::format::F;
	using ObjexxFCL::format::I;
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
	using core::import_pose::pose_from_pdb;
	using numeric::geometry::hashing::Real3;
	using numeric::geometry::hashing::Real6;

	class ResPairMotif;
	// KAB - commenting out class definition to remove warning:
	// 'MotifHit' defined as a struct here but previously declared as a class [-Wmismatched-tags]
	// class MotifHit;
	class MotifHits;
	class MotifRotamerSetOperation;
	typedef utility::pointer::owning_ptr<MotifRotamerSetOperation      > MotifRotamerSetOperationOP;
	typedef utility::pointer::owning_ptr<MotifRotamerSetOperation const> MotifRotamerSetOperationCOP;


	Real align_motif_pose_super(
		Pose & pose,
		Pose const & paln1, Size const & ir,
		Pose const & paln2, Size const & jr
	);
	Real align_motif_pose_break(
		Pose & pose,
		Pose const & paln1, Size const & ir,
		Pose const & paln2, Size const & jr
	);


// utility
	// char ss_pair_to_code(char const & ss1, char const & ss2);
	static inline float fasterlog (float x){
		union { float f; uint32_t i; }
		vx = { x };
		float y = vx.i;
		y *= 8.2629582881927490e-8f;
		return y - 87.989971088f;
	}

// XformScore
	class XformScore  : public utility::pointer::ReferenceCount {
	 public:
		typedef uint64_t Key;
		typedef double Count;
		typedef numeric:: geometry::hashing::SixDCoordinateBinner SixDCoordinateBinner;
		typedef numeric::geometry::hashing::bin_index_hasher bin_index_hasher;
		SixDCoordinateBinner hasher_;
		typedef boost::unordered_map< Key, Count, bin_index_hasher > CountMap;
		CountMap total_count_;
		double cart_size_,cart_resl_,angle_resl_;
		Size num_samples_;
	public:
		XformScore();
		void add_xform_count(Real6 const & xform);
		// Key key_of(Real6 const & xform);
		void set_score(Key const &, Count const & count);
		void add_score(Key const &, Count const & count);
		Count score_of_bin(Real6 const & xform) const;
		Count score_of_bin(Key const & key) const;
		Real  score_of_bin_normalized(Real6 const & xform) const;
		bool write_binary(std::ostream & out) const;
		bool write_binary(std::string const & fname) const;
		bool read_binary(std::istream & in, bool clearme=true);
		bool read_binary(std::string const & fname, bool clearme=true);
		bool read_binary(utility::vector1<std::string> const & fnames);
		Size num_bins() const;
		Size num_samples() const;
		void print_scores(std::ostream & out, std::string const tag="XH ") const;
		void print_scores_norm(std::ostream & out, std::string const tag="XH ") const;
		void clear();
		Real cart_size() const { return cart_size_; }
		Key bin_index(Real6 const & rt) const;
		void print_quantiles(std::ostream & out, int num) const;
		friend std::ostream & operator << (std::ostream & out, XformScore const & x);
	};
	std::ostream & operator<<(std::ostream & out, XformScore const & x);

class ResPairMotifs;

// simple motif
class ResPairMotif {
 private:
	/*  4 */ char ss1_,ss2_,aa1_,aa2_;
	/* 12 */ utility::fixedsizearray1<uint16_t,6> rt6_;
	/*  8 */ utility::fixedsizearray1<uint8_t,4> chi1_,chi2_;
	/*  6 */ uint8_t fa_atr_,hb_sc_,hb_bb_sc_,hb_bb_,sasa_,nbrs_;
	/*  6 */ uint8_t misc1_,misc2_,misc3_,misc4_,misc5_,count_;
	/*  6 */ uint16_t lg1x_,lg1y_,lg1z_; // dummy for lig / metal / water
	/*  6 */ uint16_t lg2x_,lg2y_,lg2z_; // euler ang for lig?
	/*  2 */ uint8_t  bfac1_,bfac2_;
	/*  2 */ char chain1_,chain2_;
 public:
	/*  4 */ uint16_t resi1_,resi2_;
	/*  8 */ utility::fixedsizearray1<unsigned char,8> pdb_; // simplify this // last3 are lig resn
 public:
		ResPairMotif();
		ResPairMotif(
			std::string const & tag,
			core::pose::Pose const & pose,
			Real6 const & xform,
			Size const & resi1,
			Size const & resi2,
			Real const & sasa,
			Real const & nbrs,
			Real const & fa_atr,
			Real const & hb_sc,
			Real const & hb_bb_sc,
			Real const & hb_bb,
			Real const & bfac1,
			Real const & bfac2,
			bool is_sspaired
		);
		void reverse_in_place();
		void reset();
		Real6 rt() const;
		void  rt(Real6 const & rt);
		numeric::HomogeneousTransform<Real> ht() const;
		void ht(numeric::HomogeneousTransform<Real> const & h);
		Real dist2   () const;
		Real bfac1   () const;
		Real bfac2   () const;
		Real sasa    () const;
		Real fa_atr  () const;
		Real hb_sc   () const;
		Real hb_bb_sc() const;
		Real hb_bb   () const;
		Real nbrs    () const;
		Real chi11() const;
		Real chi12() const;
		Real chi13() const;
		Real chi14() const;
		Real chi21() const;
		Real chi22() const;
		Real chi23() const;
		Real chi24() const;
		void bfac1   (Real const & bfac1   );
		void bfac2   (Real const & bfac2   );
		void sasa    (Real const & sasa    );
		void fa_atr  (Real const & fa_atr  );
		void hb_sc   (Real const & hb_sc   );
		void hb_bb_sc(Real const & hb_bb_sc);
		void hb_bb   (Real const & hb_bb   );
		void nbrs    (Real const & nbrs    );
		char aa1() const;
		char aa2() const;
		char ss1() const;
		char ss2() const;
		char dssp1() const;
		char dssp2() const;
		// char sscode() const;
		Size count() const { return count_; }
		void count(Size const & count) { count_ = count; }
		void addcount(){ if(count_<255) ++count_; }

		bool filter() const;

		/*TODO;*/
		void dump_pdb(std::string const & fname, std::string const & tag="1") const;
		void dump_pdb(std::ostream & out, numeric::xyzTransform<Real> const & xform, std::string const & tag="1") const;
		Real dump_aligned_motif(std::ostream & out, Pose const & pose1, Size const & ir, Pose const & pose2, Size const & jr, int const & num=1) const;
		Real dump_aligned_motif_break(std::ostream & out, Pose const & pose1, Size const & ir, Pose const & pose2, Size const & jr, int const & num=1) const;
		Real dump_aligned_motif_super(std::ostream & out, Pose const & pose1, Size const & ir, Pose const & pose2, Size const & jr, int const & num=1) const;
		Real dump_aligned_motif(std::ostream & out, Pose const & pose , Size const & ir                    , Size const & jr, int const & num=1) const;
		void fill_pose_with_motif(Pose & pose, int ir=1, int jr=2) const;
		// io
		friend std::ostream & operator << (std::ostream & out, ResPairMotif const & x);
		// friend std::istream & operator >> (std::istream & in , ResPairMotif & x);
		// friend class boost::serialization::access;
	 //    template<class Archive> void serialize(Archive & ar, const unsigned int version);
		static void print_header(std::ostream & out);
	};

	std::ostream & operator << (std::ostream & out, ResPairMotif const & x);
	std::istream & operator >> (std::istream & in , ResPairMotif       & x);

	bool read_motifs_binary (        string  const & fname , ResPairMotifs       & motifs);
	bool read_motifs_binary (  std::istream        & in    , ResPairMotifs       & motifs);
	bool write_motifs_binary(  std::ostream        & out   , ResPairMotifs const & motifs);
	bool write_motifs_binary(        string  const & fname , ResPairMotifs const & motifs);
	bool read_motifs_binary (vector1<string> const & fnames, ResPairMotifs       & motifs);
	// Undefined, commenting out to fix PyRosetta build  void read_motifs_boost  (vector1<string> const & fnames, ResPairMotifs       & motifs);
	// void read_motifs_boost  (        string  const & fname , ResPairMotifs       & motifs);
	// void write_motifs_boost (        string  const & fname , ResPairMotifs const & motifs);

	void filter_motifs(
		ResPairMotifs const & motifs_in,
		ResPairMotifs       & motifs_out
	);

class ResPairMotifs : public utility::vector1<ResPairMotif> {
public:
	void filter_structurally_identical_motifs();
	void filter_structurally_similar_motifs(Real cart_size, Real cart_resl, Real ang_resl);
};

// xfrag
	class Xfres {
		/*  2 */ char aa_,ss_;
		/*  6 */ uint16_t phi_,psi_,omg_;
		/*  4 */ utility::fixedsizearray1<uint8_t,4> chi_;
		public:
		Xfres() {}
		Xfres(char aa, char ss, Real phi, Real psi, Real omg, Real chi1=-180.0, Real chi2=-180.0, Real chi3=-180.0, Real chi4=-180.0);
		char aa  () const { return aa_; }
		char ss  () const { return ss_; }
		Real phi () const;
		Real psi () const;
		Real omg () const;
		Real chi1() const;
		Real chi2() const;
		Real chi3() const;
		Real chi4() const;
		void aa  (char aa  ) { aa_ = aa; }
		void ss  (char ss  ) { ss_ = ss; }
		void phi (Real _phi );
		void psi (Real _psi );
		void omg (Real _omg );
		void chi1(Real chi1);
		void chi2(Real chi2);
		void chi3(Real chi3);
		void chi4(Real chi4);
		void place_sidechain_in_pose(Pose & pose, int ir) const;
		friend std::ostream & operator << (std::ostream & out, Xfres const & x);
	};
	std::ostream & operator << (std::ostream & out, Xfres const & x);
	std::istream & operator >> (std::istream & in , Xfres       & x);
	bool read_xfres_binary (        string  const & fname , vector1<Xfres>       & Xfres);
	bool read_xfres_binary (  std::istream        & in    , vector1<Xfres>       & Xfres);
	bool write_xfres_binary(  std::ostream        & out   , vector1<Xfres> const & Xfres);
	bool write_xfres_binary(        string  const & fname , vector1<Xfres> const & Xfres);
	bool read_xfres_binary (vector1<string> const & fnames, vector1<Xfres>       & Xfres);

	class Xfrag {
		/* 12 */ utility::fixedsizearray1<uint16_t,6> rt6_;
		/*  4 */ uint32_t position_;
		char sscomp_;
		/*  4 */ uint8_t  size_,bfac_,burial_;
		/*  4 */ uint8_t  ex1_,ex2_,ex3_,ex4_;
		public:
		Xfrag(){}
		Xfrag(Real6 _rt, Size _position, uint8_t _size, char _sscomp, Real _bfac, Real _burial);
		Real6   rt6     () const;
		Size    position() const;
		Size    size    () const;
		char    sscomp  () const;
		Real    bfac    () const;
		Real    burial  () const;
		Real    ex1     () const;
		Real    ex2     () const;
		Real    ex3     () const;
		Real    ex4     () const;
		void    rt6     (Real6 const & _rt6);
		void    position(Size    _position);
		void    size    (Size    _size    );
		void    sscomp  (char    _sscomp  );
		void    bfac    (Real    _bfac    );
		void    burial  (Real    _burial  );
		void    ex1     (Real    _ex1     );
		void    ex2     (Real    _ex2     );
		void    ex3     (Real    _ex3     );
		void    ex4     (Real    _ex4     );
		void insert(Pose & pose,  vector1<Xfres> const & xfres, int lowres, int highres) const;
	};
	std::ostream & operator << (std::ostream & out, Xfrag const & x);
	std::istream & operator >> (std::istream & in , Xfrag       & x);
	bool read_xfrag_binary (        string  const & fname , vector1<Xfrag>       & xfrag);
	bool read_xfrag_binary (  std::istream        & in    , vector1<Xfrag>       & xfrag);
	bool write_xfrag_binary(  std::ostream        & out   , vector1<Xfrag> const & xfrag);
	bool write_xfrag_binary(        string  const & fname , vector1<Xfrag> const & xfrag);
	bool read_xfrag_binary (vector1<string> const & fnames, vector1<Xfrag>       & xfrag);


	bool read_xfrag_binary (        string  const & fname , vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres);
	bool read_xfrag_binary (  std::istream        & in    , vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres);
	bool write_xfrag_binary(  std::ostream        & out   , vector1<Xfrag> const & xfrag, vector1<Xfres> const & xfres);
	bool write_xfrag_binary(        string  const & fname , vector1<Xfrag> const & xfrag, vector1<Xfres> const & xfres);
	bool read_xfrag_binary (vector1<string> const & fnames, vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres);

	class XfragSet : public utility::pointer::ReferenceCount {
		typedef boost::uint64_t Key;
		typedef numeric::geometry::hashing::SixDCoordinateBinner SixDCoordinateBinner;
		typedef numeric::geometry::hashing::bin_index_hasher bin_index_hasher;
		typedef boost::unordered_multimap< Key, Xfrag, bin_index_hasher > XfragMap;

		double const cart_resl_, angle_resl_;
		// KAB - below variable commented out (-Wunused-private-field) on 2014-09-11
		// double const cart_size_;

		SixDCoordinateBinner hasher_;
		XfragMap xfragmap_;
		utility::vector1<Xfres> xfres_;
		public:
		XfragSet(core::Real cartsize=10.0, core::Real cartresl=1.0, core::Real angleresl=15.0);
		// Undefined, commenting out to fix PyRosetta build  void write_binary(std::string fname) const;
		// Undefined, commenting out to fix PyRosetta build  void read_binary(std::string fname);
		// Undefined, commenting out to fix PyRosetta build  void add_frags(utility::vector1<Xfrag> const & frags, utility::vector1<Xfres> const & res);
		// Undefined, commenting out to fix PyRosetta build  void lookup_frags(Real6 const & rt6, utility::vector1<Xfrag> & frags) const;
	};


// MotifHash

	class MotifHash : public utility::pointer::ReferenceCount {
	public:
		typedef ResPairMotif Motif;
		typedef boost::uint64_t Key;
		typedef numeric::geometry::hashing::SixDCoordinateBinner SixDCoordinateBinner;
		typedef numeric::geometry::hashing::bin_index_hasher bin_index_hasher;
		typedef boost::unordered_multimap< Key, ResPairMotif , bin_index_hasher > MotifMap;
		typedef boost::unordered_map    < Key, float , bin_index_hasher > ScoreMap;
		typedef boost::unordered_set<Key> KeySet;
		// typedef std::tr1::unordered_multimap< Key, ResPairMotif , bin_index_hasher > MotifMap;
		// typedef std::tr1::unordered_map    < Key,       float , bin_index_hasher > ScoreMap;
		typedef protocols::sic_dock::xyzStripeHashPoseCAP ClashChecker;
		double const cart_size_,cart_resl_,angle_resl_;
		SixDCoordinateBinner hasher_;
		MotifMap motif_hash_;
		ScoreMap score_hash_;
		KeySet key_set_;
	public:
		MotifHash();
		MotifHash(Real cart_size, Real cart_resl, Real ang_resl);
		MotifHash(ResPairMotifs const & motifs);
		SixDCoordinateBinner const & hasher() const { return hasher_; };
		SixDCoordinateBinner       & hasher()       { return hasher_; };
		void add_motif(Motif const & d);
		void add_motif(Motif const & d, Key const & key);
		void add_motif_shift(Motif const & d, Real6 const & shift);
		void set_score(Key const & k, float const & score);
		void find_motifs(Key const & k, ResPairMotifs & results ) const;
		void find_motifs(Real6 const & rt, ResPairMotifs & results ) const;
		int count_motifs(Real6 const & rt) const ;
		int count_motifs(Key const & k) const;
		void find_motifs_with_radius(Real6 const & rt, Real radius, vector1<Motif> & results) const;
		float motif_score(Real6 const & rt) const;
		void generate_scoring_hash();
		Size num_motifs() { return motif_hash_.size(); }
		// Undefined, commenting out to fix PyRosetta build  Key key_of(Motif const & d) const;
		Real cart_size() const { return cart_size_; }
		Key bin_index(Real6 const & rt) const;
		Key bin_index(Motif const & m) const;
		// io
		int dump_matching_motifs( Pose const & pose1, Pose const & pose2, Real radius, std::ostream & out, int & count, ClashChecker clash_check=NULL, bool print=false ) const;
		int dump_matching_motifs( Pose const & pose1, Pose const & pose2,              std::ostream & out, int & count, ClashChecker clash_check=NULL, bool print=false ) const;
		int print_matching_motifs( Pose const & pose1, Pose const & pose2, std::ostream & out, ClashChecker clash_check=NULL ) const;
		int stat_matching_motifs( Pose const & pose1, Pose const & pose2, std::map<std::string,Real> & stats, ClashChecker clash_check=NULL, Real radius=0.0 ) const;

		int get_matching_motifs( Pose const & pose1, Pose const & pose2, MotifHits & hits, ClashChecker clash_check=NULL, Real radius=0.0 ) const;
		int
		get_matching_motifs(
			Pose  const & pose1,
			Pose  const & pose2,
			Bools const & useres1,
			Bools const & useres2,
			MotifHits & hits,
			ClashChecker clash_check=NULL,
			Real radius=0.0
		) const;
		int get_matching_motifs( Pose const & pose, core::pack::task::PackerTaskCOP ptask, MotifHits & hits, ClashChecker clash_check=NULL, Real radius=0.0 ) const;

		MotifRotamerSetOperationCOP get_matching_motifs_rotsetop( Pose const & pose, core::pack::task::PackerTaskCOP ptask=NULL, Real radius=0.6, ClashChecker clash_check=NULL ) const;

		friend std::ostream & operator << (std::ostream & out, MotifHash const & x);
		// friend std::istream & operator >> (std::istream & in , MotifHash & x);
	    // friend class boost::serialization::access;
	    // template<class Archive> void serialize(Archive & ar, const unsigned int version);
	    KeySet::const_iterator begin() { return key_set_.begin(); }
	    KeySet::const_iterator end  () { return key_set_.end  (); }
	    KeySet const & key_set() const { return key_set_; }
	};
	std::ostream & operator << (std::ostream & out, MotifHash const & x);
	std::istream & operator >> (std::istream & in , MotifHash       & x);

class MotifHashManager {
public:
	static MotifHashManager * get_instance();
	MotifHashCAP motif_hash_from_cli();
	XformScoreCAP xform_score_from_cli();
	XformScoreCAP xform_score_sspair_from_cli();
	XformScoreCAP xform_score_he_from_cli();
	XformScoreCAP xform_score_eh_from_cli();
	XformScoreCAP xform_score_ee_from_cli();
	XformScoreCAP xform_score_hh_from_cli();
private:
	MotifHashManager():cli_motif_hash_(NULL),cli_xform_score_(NULL),cli_xform_score_ee_(NULL),cli_xform_score_eh_(NULL),cli_xform_score_he_(NULL),cli_xform_score_hh_(NULL),cli_xform_score_sspair_(NULL) {}
	static MotifHashManager * instance_;
	MotifHash * cli_motif_hash_;
	XformScore * cli_xform_score_,*cli_xform_score_ee_,*cli_xform_score_eh_,*cli_xform_score_he_,*cli_xform_score_hh_,*cli_xform_score_sspair_;
};

struct MotifHit {
	int residue1,residue2;
	Real score;
	ResPairMotif motif;
	MotifHit():residue1(0),residue2(0),score(0){}
	MotifHit(int res1, int res2, Real sc, ResPairMotif const & sm) : residue1(res1),residue2(res2),score(sc),motif(sm) {}
	bool operator< (MotifHit const & other) const;
	bool operator==(MotifHit const & other) const;
};
std::ostream & operator<<(std::ostream & out, MotifHit  const & h);

class MotifHits : public utility::vector1<MotifHit> {
public:
	mutable int rescover,paircover;
	mutable Real total_score;
	void compute_metrics() const;
	void filter_redundant();
};
std::ostream & operator<<(std::ostream & out, MotifHits const & h);


// "protocol" stuff
void merge_motifs();
void dump_motif_pdbs();
void harvest_scores();
void print_scores();
	// Undefined, commenting out to fix PyRosetta build  void sequence_recovery();
void dump_matching_motifs();
void harvest_motifs();
void print_motifs(std::ostream & out);


class MotifRotamerSetOperation : public core::pack::rotamer_set::RotamerSetOperation {
 public:
	MotifRotamerSetOperation( Pose const & refpose, MotifHits const & motifs );

	virtual ~MotifRotamerSetOperation(){}

	core::pack::rotamer_set::RotamerSetOperationOP clone() const { return new MotifRotamerSetOperation(*this); }

	void alter_rotamer_set(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::pack::task::PackerTask const & ptask,
		core::graph::GraphCOP packer_neighbor_graph,
		core::pack::rotamer_set::RotamerSet & rotamer_set
	);

	virtual	Real increase_packer_residue_radius(
		core::pose::Pose const & pose,
		core::pack::task::PackerTaskCOP the_task,
		core::Size residue_in
	) const;

	void dump_pdb(std::string const & fname) const;

 private:
	MotifHits motifs_;
	utility::vector1<utility::vector1<PoseCOP> > res1_poses_;
	utility::vector1<utility::vector1<PoseCOP> > res2_poses_;


 };



}
}


#endif
