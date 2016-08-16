// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/motif/util.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/mh.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/pymol_chains.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/motif/reference_frames.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/xyzStripeHashPose.hh>
#include <core/io/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/dssp/StrandPairing.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/xyzVector.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/string_constants.hh>
#include <numeric/xyzTransform.hh>

#include <numeric/geometry/hashing/SixDHasher.hh>

#include <boost/unordered_set.hpp>

#include <boost/foreach.hpp>
#include <bitset>

#ifdef USE_OPENMP
#include <omp.h>
#ifndef _WIN32
#include <pthread.h>
#endif
#endif


#define MAX_UINT16 65535
#define MAX_UINT8    255


namespace core {
namespace scoring {
namespace motif {

extern const numeric::Real MOTIF_HASH_CART_SIZE = 16.0;
extern const uint64_t      MOTIF_VERSION        = 2;

using numeric::Xforms;


/************************************************* types ************************************************/
static basic::Tracer TR("core.scoring.motif");

using core::pose::PoseCoordPickMode_BB;
using core::pose::PoseCoordPickMode_N_C_O;

typedef utility::fixedsizearray1<float,20> float20;
typedef utility::fixedsizearray1<float,9> float9;
using std::make_pair;
using core::chemical::AA;
using core::id::AtomID;
using basic::options::option;
namespace mh = basic::options::OptionKeys::mh;
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
using std::ostream;
using std::string;
using utility::io::izstream;
using utility::io::ozstream;
using utility::file_basename;
using utility::vector1;
using std::endl;
//using core::import_pose::pose_from_file;
using numeric::geometry::hashing::Real3;
using numeric::geometry::hashing::Real6;
using core::pose::xyzStripeHashPoseCOP;

/************************************************ helpers ***********************************************/
static inline double sqr(double const & x) { return x*x; }
static inline float  sqr(float  const & x) { return x*x; }


static inline Real uint16_to_real(uint16_t const & val, Real const & lb, Real const & ub){
	return min(ub,max(lb,((Real)val)/((Real)MAX_UINT16) * (ub-lb) + lb));
}
static inline uint16_t real_to_uint16(Real const & val, Real const & lb, Real const & ub){
	return (uint16_t)( min(1.0,max(0.0,(val-lb)/(ub-lb))) * (Real)MAX_UINT16 );
}

static inline Real uint8_to_real(uint8_t const & val, Real const & lb, Real const & ub){
	return min(ub,max(lb,((Real)val)/((Real)MAX_UINT8) * (ub-lb) + lb));
}
static inline uint8_t real_to_uint8(Real const & val, Real const & lb, Real const & ub){
	return (uint16_t)( min(1.0,max(0.0,(val-lb)/(ub-lb))) * (Real)MAX_UINT8 );
}

/********************************************** ResPairMotif ********************************************/
void ResPairMotif::reset(){
	pdb_=0;
	rt6_=0;
	resi1_=0;
	resi2_=0;
	chain1_=0;
	chain2_=0;
	ss1_=0;
	ss2_=0;
	aa1_=0;
	aa2_=0;
	chi1_[1]=0;
	chi1_[2]=0;
	chi1_[3]=0;
	chi1_[4]=0;
	chi2_[1]=0;
	chi2_[2]=0;
	chi2_[3]=0;
	chi2_[4]=0;
	fa_atr_=0;
	fa_atr_sc_bb_=0;
	fa_atr_bb_=0;
	hb_sc_=0;
	hb_bb_sc_=0;
	hb_bb_=0;
	nbrs1_=0;
	nbrs2_=0;
	count_=1;
	type_=0;
	misc1_=0;
	bfac1_=0;
	bfac2_=0;
	lg1x_=0;
	lg1y_=0;
	lg1z_=0;
	lg2x_=0;
	lg2y_=0;
	lg2z_=0;
}
ResPairMotif::ResPairMotif(){
	reset();
}
ResPairMotif::ResPairMotif(
	string const & tag,
	Pose const & _pose,
	Real6 const & _xform,
	Size const & _resi1,
	Size const & _resi2,
	Real const & _nbrs1,
	Real const & _nbrs2,
	Real const & _fa_atr,
	Real const & _fa_atr_sc_bb,
	Real const & _fa_atr_bb,
	Real const & _hb_sc,
	Real const & _hb_bb_sc,
	Real const & _hb_bb,
	Real const & _bfac1,
	Real const & _bfac2,
	bool is_sspaired,
	RPM_Type _type
){
	reset();
	if ( tag.size()!=5 && tag.size()!=8 ) utility_exit_with_message("bad tag for ResPairMotif not nchar 5 or 8");
	pdb_[1]=tag[0]; pdb_[2]=tag[1]; pdb_[3]=tag[2]; pdb_[4]=tag[3]; pdb_[5]=tag[4];
	if ( tag.size()==8 ) {
		pdb_[6]=tag[5]; pdb_[7]=tag[6]; pdb_[8]=tag[7];
	} else {
		pdb_[6]='_'; pdb_[7]='_'; pdb_[8]='_';
	}

	this->rt(_xform);

	resi1_  = _resi1;
	resi2_  = _resi2;
	if ( is_sspaired ) {
		ss1_ = 'P';
		ss2_ = 'P';
	} else {
		ss1_ = _pose.secstruct(_resi1);
		ss2_ = _pose.secstruct(_resi2);
	}
	aa1_ = _pose.residue(_resi1).name1();
	aa2_ = _pose.residue(_resi2).name1();
	chain1_ = basic::get_pymol_chain(_pose.chain(_resi1));
	chain2_ = basic::get_pymol_chain(_pose.chain(_resi2));

	this->bfac1   (_bfac1);
	this->bfac2   (_bfac2);
	this->nbrs1   (_nbrs1);
	this->nbrs2   (_nbrs2);
	this->fa_atr  (_fa_atr);
	this->fa_atr_sc_bb(_fa_atr_sc_bb);
	this->fa_atr_bb(_fa_atr_bb);
	this->hb_sc   (_hb_sc);
	this->hb_bb_sc(_hb_bb_sc);
	this->hb_bb   (_hb_bb);
	this->type    (_type);

	for ( Size i = 1; i <= _pose.residue(_resi1).nchi(); ++i ) {
		Real chi = _pose.residue(_resi1).chi(i);
		// if(i==flip_chi_ir){
		//  chi = chi + 180.0;
		//  if(chi >= 180.0) chi -= 360.0;
		//  cout << "ResPairMotif c'tor: flipped CHI " << i << endl;
		// }
		chi1_[i] = real_to_uint8(chi,-180.0,180.0);
	}
	for ( Size i = 1; i <= _pose.residue(_resi2).nchi(); ++i ) {
		Real chi = _pose.residue(_resi2).chi(i);
		// if(i==flip_chi_jr){
		//  chi = chi + 180.0;
		//  if(chi >= 180.0) chi -= 360.0;
		//  cout << "ResPairMotif c'tor: flipped CHI " << i << endl;
		// }
		chi2_[i] = real_to_uint8(chi,-180.0,180.0);
	}
}
char ResPairMotif::aa1  () const { return aa1_; }
char ResPairMotif::aa2  () const { return aa2_; }
char ResPairMotif::ss1  () const { if ( ss1_=='H'||ss1_=='I'||ss1_=='G' ) return 'H'; else if ( ss1_=='P'||ss1_=='E'||ss1_=='U'||ss1_=='B' ) return 'E'; else if ( ss1_=='S'||ss1_=='T'||ss1_==' '||ss1_=='L' ) return 'L'; else utility_exit_with_message("unknown SS "+string_of(ss1_)); }
char ResPairMotif::ss2  () const { if ( ss2_=='H'||ss2_=='I'||ss2_=='G' ) return 'H'; else if ( ss2_=='P'||ss2_=='E'||ss2_=='U'||ss2_=='B' ) return 'E'; else if ( ss2_=='S'||ss2_=='T'||ss2_==' '||ss2_=='L' ) return 'L'; else utility_exit_with_message("unknown SS "+string_of(ss2_)); }
char ResPairMotif::dssp1() const { return ss1_; }
char ResPairMotif::dssp2() const { return ss2_; }
bool ResPairMotif::operator < (ResPairMotif const & other) const {
	return 0 <  memcmp(this,&other,sizeof(ResPairMotif));
}
bool ResPairMotif::operator ==(ResPairMotif const & other) const {
	return 0 == memcmp(this,&other,sizeof(ResPairMotif));
}

ResPairMotif & ResPairMotif::reverse_in_place_unsafe(){
	xform( xform().inverse() );
	std::swap(resi1_  ,resi2_);
	std::swap(chain1_ ,chain2_);
	std::swap(ss1_    ,ss2_);
	std::swap(aa1_    ,aa2_);
	std::swap(chi1_   ,chi2_);
	std::swap(bfac1_  ,bfac2_);
	std::swap(nbrs1_  ,nbrs2_);
	// lg1x_,lg1y_,lg1z_; // dummy for lig / metal / water ???????????
	// lg2x_,lg2y_,lg2z_; // euler ang for lig? ???????????
	return *this;
}
ResPairMotif & ResPairMotif::reverse_in_place(){
	runtime_assert(type1()==type2());
	return reverse_in_place_unsafe();
}
ResPairMotif ResPairMotif::reversed() const {
	ResPairMotif motif(*this);
	motif.reverse_in_place();
	return motif;
}
Real6 ResPairMotif::rt() const {
	Real6 rt6;
	for ( Size i = 1; i <= 3; ++i ) rt6[i] = uint16_to_real(rt6_[i],-MOTIF_HASH_CART_SIZE,  MOTIF_HASH_CART_SIZE );
	for ( Size i = 4; i <= 5; ++i ) rt6[i] = uint16_to_real(rt6_[i],  0.0, 360.0 );
	for ( Size i = 6; i <= 6; ++i ) rt6[i] = uint16_to_real(rt6_[i],  0.0, 180.0 );
	return rt6;
}
void  ResPairMotif::rt(Real6 const & rt_in) {
	Real6 rt6 = rt_in;
	rt6[4] = fmod(rt_in[4],360.0);
	rt6[5] = fmod(rt_in[5],360.0);
	rt6[6] = fmod(rt_in[6],360.0);
	rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
	rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
	rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
	for ( Size i = 1; i <= 3; ++i ) rt6_[i] = real_to_uint16(rt_in[i],-MOTIF_HASH_CART_SIZE,  MOTIF_HASH_CART_SIZE );
	for ( Size i = 4; i <= 5; ++i ) rt6_[i] = real_to_uint16(rt_in[i],  0.0, 360.0 );
	for ( Size i = 6; i <= 6; ++i ) rt6_[i] = real_to_uint16(rt_in[i],  0.0, 180.0 );
}
Xform ResPairMotif::xform() const{
	Real6 rt6 = rt();
	Xform x;
	x.from_euler_angles_deg(numeric::xyzVector<Real>(rt6[4],rt6[5],rt6[6]));
	x.t = (Vec(rt6[1],rt6[2],rt6[3]));
	return x;
}
void ResPairMotif::xform(Xform const & x){
	numeric::xyzVector < Real > euler_angles =  x.euler_angles_deg();
	Real6 rt6;
	rt6[1] = x.px();
	rt6[2] = x.py();
	rt6[3] = x.pz();
	rt6[4] = fmod(euler_angles.x(),360.0);
	rt6[5] = fmod(euler_angles.y(),360.0);
	rt6[6] = fmod(euler_angles.z(),360.0);
	rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
	rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
	rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
	if ( rt6[6] > 180.0 ) utility_exit_with_message("bad rt6[6]");
	rt(rt6);
}

Real ResPairMotif::dist2() const {
	Real sqdis = 0.0;
	for ( Size i = 1; i <= 3; ++i ) {
		Real d = uint16_to_real(rt6_[i],-MOTIF_HASH_CART_SIZE,MOTIF_HASH_CART_SIZE);
		sqdis += d*d;
	}
	return sqdis;
}

Real ResPairMotif::       bfac1() const { return uint8_to_real( bfac1_        ,   0.0,   2.55  ); }
Real ResPairMotif::       bfac2() const { return uint8_to_real( bfac2_        ,   0.0,   2.55  ); }
Real ResPairMotif::       nbrs1() const { return uint8_to_real( nbrs1_        ,   0.0,  63.75  ); }
Real ResPairMotif::       nbrs2() const { return uint8_to_real( nbrs2_        ,   0.0,  63.75  ); }
Real ResPairMotif::fa_atr      () const { return uint8_to_real( fa_atr_       ,  -8.0,   0.0   ); }
Real ResPairMotif::fa_atr_sc_bb() const { return uint8_to_real( fa_atr_sc_bb_ ,  -8.0,   0.0   ); }
Real ResPairMotif::fa_atr_bb   () const { return uint8_to_real( fa_atr_bb_    ,  -8.0,   0.0   ); }
Real ResPairMotif::       hb_sc() const { return uint8_to_real( hb_sc_        ,  -8.0,   0.0   ); }
Real ResPairMotif::    hb_bb_sc() const { return uint8_to_real( hb_bb_sc_     ,  -8.0,   0.0   ); }
Real ResPairMotif::       hb_bb() const { return uint8_to_real( hb_bb_        ,  -8.0,   0.0   ); }
Real ResPairMotif::       chi11() const { return uint8_to_real( chi1_[1]      ,-180.0, 180.0   ); }
Real ResPairMotif::       chi12() const { return uint8_to_real( chi1_[2]      ,-180.0, 180.0   ); }
Real ResPairMotif::       chi13() const { return uint8_to_real( chi1_[3]      ,-180.0, 180.0   ); }
Real ResPairMotif::       chi14() const { return uint8_to_real( chi1_[4]      ,-180.0, 180.0   ); }
Real ResPairMotif::       chi21() const { return uint8_to_real( chi2_[1]      ,-180.0, 180.0   ); }
Real ResPairMotif::       chi22() const { return uint8_to_real( chi2_[2]      ,-180.0, 180.0   ); }
Real ResPairMotif::       chi23() const { return uint8_to_real( chi2_[3]      ,-180.0, 180.0   ); }
Real ResPairMotif::       chi24() const { return uint8_to_real( chi2_[4]      ,-180.0, 180.0   ); }

void ResPairMotif::bfac1       (Real const & _bfac1       ){ bfac1_        = real_to_uint8( _bfac1        ,  0.0,   2.55  ); }
void ResPairMotif::bfac2       (Real const & _bfac2       ){ bfac2_        = real_to_uint8( _bfac2        ,  0.0,   2.55  ); }
void ResPairMotif::nbrs1       (Real const & _nbrs1       ){ nbrs1_        = real_to_uint8( _nbrs1        ,  0.0,  63.75  ); }
void ResPairMotif::nbrs2       (Real const & _nbrs2       ){ nbrs2_        = real_to_uint8( _nbrs2        ,  0.0,  63.75  ); }
void ResPairMotif::fa_atr      (Real const & _fa_atr      ){ fa_atr_       = real_to_uint8( _fa_atr       , -8.0,   0.0  ); }
void ResPairMotif::fa_atr_sc_bb(Real const & _fa_atr_sc_bb){ fa_atr_sc_bb_ = real_to_uint8( _fa_atr_sc_bb , -8.0,   0.0  ); }
void ResPairMotif::fa_atr_bb   (Real const & _fa_atr_bb   ){ fa_atr_bb_    = real_to_uint8( _fa_atr_bb    , -8.0,   0.0  ); }
void ResPairMotif::hb_sc       (Real const & _hb_sc       ){ hb_sc_        = real_to_uint8( _hb_sc        , -8.0,   0.0  ); }
void ResPairMotif::hb_bb_sc    (Real const & _hb_bb_sc    ){ hb_bb_sc_     = real_to_uint8( _hb_bb_sc     , -8.0,   0.0  ); }
void ResPairMotif::hb_bb       (Real const & _hb_bb       ){ hb_bb_        = real_to_uint8( _hb_bb        , -8.0,   0.0  ); }

Real ResPairMotif::chi1(Size const & ichi) const {
	switch(ichi){
	case 1 : return chi11();
	case 2 : return chi12();
	case 3 : return chi13();
	case 4 : return chi14();
	default : utility_exit_with_message("unknown CHI");
	}
	return 0;
}
Real ResPairMotif::chi2(Size const & ichi) const {
	switch(ichi){
	case 1 : return chi21();
	case 2 : return chi22();
	case 3 : return chi23();
	case 4 : return chi24();
	default : utility_exit_with_message("unknown CHI");
	}
	return 0;
}

Real ResPairMotif::score() const{
	// this needs an audit, along with the rest of the metrics, not that there are types....
	switch(type()){
	case SC_SC :
		return 0.75*fa_atr() + 0.00*fa_atr_sc_bb() + 0.00*fa_atr_bb() + 1.00*hb_sc() + 0.00*hb_bb_sc() + 0.00*hb_bb();
	case SC_BB :
		return 0.75*fa_atr() + 0.50*fa_atr_sc_bb() + 0.00*fa_atr_bb() + 1.00*hb_sc() + 0.00*hb_bb_sc() + 0.00*hb_bb();
	case BB_BB :
		return 0.75*fa_atr() + 0.50*fa_atr_sc_bb() + 0.25*fa_atr_bb() + 1.00*hb_sc() + 0.00*hb_bb_sc() + 0.00*hb_bb();
	case SC_PH:
	case SC_PO :
		return 0.00*fa_atr() + 0.75*fa_atr_sc_bb() + 0.00*fa_atr_bb() + 0.00*hb_sc() + 1.00*hb_bb_sc() + 0.00*hb_bb();
	case BB_PH:
	case BB_PO :
		return 0.00*fa_atr() + 0.00*fa_atr_sc_bb() + 0.00*fa_atr_bb() + 0.00*hb_sc() + 1.00*hb_bb_sc() + 0.00*hb_bb();
	case PH_PO :
		return 0.00*fa_atr() + 0.00*fa_atr_sc_bb() + 0.75*fa_atr_bb() + 0.00*hb_sc() + 0.00*hb_bb_sc() + 1.00*hb_bb();
	default : utility_exit_with_message("atartart");
	}
}
RPM_Type ResPairMotif::type() const {
	runtime_assert(type_!=(uint8_t)RPM_Type_NONE);
	return (RPM_Type)type_;
}
RM_Type ResPairMotif::type1() const {
	return rpm_type1(type());
}
RM_Type ResPairMotif::type2() const {
	return rpm_type2(type());
}


RPM_FilterStats::RPM_FilterStats(){
	memset(this,0,sizeof(RPM_FilterStats));
}


std::ostream & operator << (std::ostream & out, RPM_FilterStats const & s){
	out << "RPM_FilterStats     FAIL        PASS" << endl;
	if ( s.F_type_SC_SC     ) out << "type_SC_SC     " << I(9,s.F_type_SC_SC     ) <<" / "<< I(9,s.P_type_SC_SC      ) << endl;
	if ( s.F_type_SC_BB     ) out << "type_SC_BB     " << I(9,s.F_type_SC_BB     ) <<" / "<< I(9,s.P_type_SC_BB      ) << endl;
	if ( s.F_type_SC_PH     ) out << "type_SC_PH     " << I(9,s.F_type_SC_PH     ) <<" / "<< I(9,s.P_type_SC_PH      ) << endl;
	if ( s.F_type_SC_PO     ) out << "type_SC_PO     " << I(9,s.F_type_SC_PO     ) <<" / "<< I(9,s.P_type_SC_PO      ) << endl;
	if ( s.F_type_BB_BB     ) out << "type_BB_BB     " << I(9,s.F_type_BB_BB     ) <<" / "<< I(9,s.P_type_BB_BB      ) << endl;
	if ( s.F_type_BB_PH     ) out << "type_BB_PH     " << I(9,s.F_type_BB_PH     ) <<" / "<< I(9,s.P_type_BB_PH      ) << endl;
	if ( s.F_type_BB_PO     ) out << "type_BB_PO     " << I(9,s.F_type_BB_PO     ) <<" / "<< I(9,s.P_type_BB_PO      ) << endl;
	if ( s.F_type_PH_PO     ) out << "type_PH_PO     " << I(9,s.F_type_PH_PO     ) <<" / "<< I(9,s.P_type_PH_PO      ) << endl;
	if ( s.F_pdb            ) out << "pdb            " << I(9,s.F_pdb            ) <<" / "<< I(9,s.P_pdb             ) << endl;
	if ( s.F_lig            ) out << "lig            " << I(9,s.F_lig            ) <<" / "<< I(9,s.P_lig             ) << endl;
	if ( s.F_seqsep         ) out << "seqsep         " << I(9,s.F_seqsep         ) <<" / "<< I(9,s.P_seqsep          ) << endl;
	if ( s.F_max_seqsep     ) out << "max_seqsep     " << I(9,s.F_max_seqsep     ) <<" / "<< I(9,s.P_max_seqsep      ) << endl;
	if ( s.F_no_hb_bb       ) out << "no_hb_bb       " << I(9,s.F_no_hb_bb       ) <<" / "<< I(9,s.P_no_hb_bb        ) << endl;
	if ( s.F_ss1            ) out << "ss1            " << I(9,s.F_ss1            ) <<" / "<< I(9,s.P_ss1             ) << endl;
	if ( s.F_ss2            ) out << "ss2            " << I(9,s.F_ss2            ) <<" / "<< I(9,s.P_ss2             ) << endl;
	if ( s.F_dssp1          ) out << "dssp1          " << I(9,s.F_dssp1          ) <<" / "<< I(9,s.P_dssp1           ) << endl;
	if ( s.F_dssp2          ) out << "dssp2          " << I(9,s.F_dssp2          ) <<" / "<< I(9,s.P_dssp2           ) << endl;
	if ( s.F_aa1            ) out << "aa1            " << I(9,s.F_aa1            ) <<" / "<< I(9,s.P_aa1             ) << endl;
	if ( s.F_aa2            ) out << "aa2            " << I(9,s.F_aa2            ) <<" / "<< I(9,s.P_aa2             ) << endl;
	if ( s.F_score          ) out << "score          " << I(9,s.F_score          ) <<" / "<< I(9,s.P_score           ) << endl;
	if ( s.F_faatr          ) out << "faatr          " << I(9,s.F_faatr          ) <<" / "<< I(9,s.P_faatr           ) << endl;
	if ( s.F_hb_sc          ) out << "hb_sc          " << I(9,s.F_hb_sc          ) <<" / "<< I(9,s.P_hb_sc           ) << endl;
	if ( s.F_hb_bb_sc       ) out << "hb_bb_sc       " << I(9,s.F_hb_bb_sc       ) <<" / "<< I(9,s.P_hb_bb_sc        ) << endl;
	if ( s.F_hb_bb          ) out << "hb_bb          " << I(9,s.F_hb_bb          ) <<" / "<< I(9,s.P_hb_bb           ) << endl;
	if ( s.F_coorderr       ) out << "coorderr       " << I(9,s.F_coorderr       ) <<" / "<< I(9,s.P_coorderr        ) << endl;
	if ( s.F_uniformfrag    ) out << "uniformfrag    " << I(9,s.F_uniformfrag    ) <<" / "<< I(9,s.P_uniformfrag     ) << endl;
	if ( s.F_mindist2       ) out << "mindist2       " << I(9,s.F_mindist2       ) <<" / "<< I(9,s.P_mindist2        ) << endl;
	if ( s.F_maxdist2       ) out << "maxdist2       " << I(9,s.F_maxdist2       ) <<" / "<< I(9,s.P_maxdist2        ) << endl;
	if ( s.F_faatr_or_hbbb  ) out << "faatr_or_hbbb  " << I(9,s.F_faatr_or_hbbb  ) <<" / "<< I(9,s.P_faatr_or_hbbb   ) << endl;
	if ( s.F_faatr_or_hb    ) out << "faatr_or_hb    " << I(9,s.F_faatr_or_hb    ) <<" / "<< I(9,s.P_faatr_or_hb     ) << endl;
	if ( s.F_noloops        ) out << "noloops        " << I(9,s.F_noloops        ) <<" / "<< I(9,s.P_noloops         ) << endl;
	if ( s.F_oneloop        ) out << "oneloop        " << I(9,s.F_oneloop        ) <<" / "<< I(9,s.P_oneloop         ) << endl;
	if ( s.F_nodisulf       ) out << "nodisulf       " << I(9,s.F_nodisulf       ) <<" / "<< I(9,s.P_nodisulf        ) << endl;
	if ( s.F_restype1       ) out << "restype1       " << I(9,s.F_restype1       ) <<" / "<< I(9,s.P_restype1        ) << endl;
	if ( s.F_restype2       ) out << "restype2       " << I(9,s.F_restype2       ) <<" / "<< I(9,s.P_restype2        ) << endl;
	if ( s.F_restype        ) out << "restype        " << I(9,s.F_restype        ) <<" / "<< I(9,s.P_restype         ) << endl;
	if ( s.F_not_restype    ) out << "not_restype    " << I(9,s.F_not_restype    ) <<" / "<< I(9,s.P_not_restype     ) << endl;
	if ( s.F_restype_one    ) out << "restype_one    " << I(9,s.F_restype_one    ) <<" / "<< I(9,s.P_restype_one     ) << endl;
	//  if(s.F_occupancy      ) out << "occupancy      " << I(9,s.F_occupancy      ) <<" / "<< I(9,s.P_occupancy       ) << endl;
	if ( s.F_not_restype_one ) out << "not_restype_one" << I(9,s.F_not_restype_one) <<" / "<< I(9,s.P_not_restype_one ) << endl;
	return out;
}

bool ResPairMotif::filter(RPM_FilterStats *s) const {
	using std::string;
	namespace f = basic::options::OptionKeys::mh::filter;
	basic::options::OptionCollection & o(basic::options::option);
	Real trust = aa_trustworthiness(aa1()) * aa_trustworthiness(aa2());
	// do this once per run... shouldn't slow things down
	static bool io_check = true;
	Real const & CSIZE(MOTIF_HASH_CART_SIZE);
	if(io_check) {
	if ( o[f::motif_type].user() ) {
		if (
				"SC_SC"!=o[f::motif_type]() &&
				"SC_BB"!=o[f::motif_type]() &&
				"SC_PH"!=o[f::motif_type]() &&
				"SC_PO"!=o[f::motif_type]() &&
				"BB_BB"!=o[f::motif_type]() &&
				"BB_PH"!=o[f::motif_type]() &&
				"BB_PO"!=o[f::motif_type]() &&
				"PH_PO"!=o[f::motif_type]() && true
				) utility_exit_with_message("unknown motif_type request "+o[f::motif_type]()+" allowed types: SC_SC SC_BB SC_PH SC_PO BB_BB BB_PH BB_PO PH_PO");
		// disable this 'autocorrect' for now
		// // put in canonical order so only one check required
		// if( o[f::motif_type]()!="SB" ) o[f::motif_type]("BS");
		// if( o[f::motif_type]()!="PB" ) o[f::motif_type]("BP");
		// if( o[f::motif_type]()!="PS" ) o[f::motif_type]("SP");
	}
	if ( o[f::pdb].user() && ( o[f::pdb]().size() < 4 ||  o[f::pdb]().size() > 5 ) ) utility_exit_with_message("bad pdb code "+o[f::pdb]()); // FIXME
	io_check = false;
	}
	if ( o[f::motif_type].user() ) {
		switch(type()){
		case SC_SC : if ( o[f::motif_type]()!="SC_SC" ) { if ( s ) { ++s->F_type_SC_SC; } return false; } else if ( s ) { ++s->P_type_SC_SC; } break;
		case SC_BB : if ( o[f::motif_type]()!="SC_BB" ) { if ( s ) { ++s->F_type_SC_BB; } return false; } else if ( s ) { ++s->P_type_SC_BB; } break;
		case SC_PH : if ( o[f::motif_type]()!="SC_PH" ) { if ( s ) { ++s->F_type_SC_PH; } return false; } else if ( s ) { ++s->P_type_SC_PH; } break;
		case SC_PO : if ( o[f::motif_type]()!="SC_PO" ) { if ( s ) { ++s->F_type_SC_PO; } return false; } else if ( s ) { ++s->P_type_SC_PO; } break;
		case BB_BB : if ( o[f::motif_type]()!="BB_BB" ) { if ( s ) { ++s->F_type_BB_BB; } return false; } else if ( s ) { ++s->P_type_BB_BB; } break;
		case BB_PH : if ( o[f::motif_type]()!="BB_PH" ) { if ( s ) { ++s->F_type_BB_PH; } return false; } else if ( s ) { ++s->P_type_BB_PH; } break;
		case BB_PO : if ( o[f::motif_type]()!="BB_PO" ) { if ( s ) { ++s->F_type_BB_PO; } return false; } else if ( s ) { ++s->P_type_BB_PO; } break;
		case PH_PO : if ( o[f::motif_type]()!="PH_PO" ) { if ( s ) { ++s->F_type_PH_PO; } return false; } else if ( s ) { ++s->P_type_PH_PO; } break;
		default : utility_exit_with_message("unknown ResPairMotiftype");
		}
	}
	if ( o[f::pdb       ].user() && !this->check_pdb_code(o[f::pdb]()) )                                                                  { if ( s ) { ++s->F_pdb;           } return false; } else if ( s ) { ++s->P_pdb;            }
	if ( o[f::lig       ].user() && !this->check_lig_code(o[f::lig]()) )                                                                  { if ( s ) { ++s->F_lig;           } return false; } else if ( s ) { ++s->P_lig;            }
	if ( o[f::seqsep    ].user() &&  o[f::seqsep]() >= std::abs(resi1_-resi2_) && chain1_==chain2_ )                                      { if ( s ) { ++s->F_seqsep;        } return false; } else if ( s ) { ++s->P_seqsep;         }
	if ( o[f::max_seqsep].user() && (o[f::seqsep]() <  std::abs(resi1_-resi2_) || chain1_!=chain2_ ) )                                     { if ( s ) { ++s->F_max_seqsep;    } return false; } else if ( s ) { ++s->P_max_seqsep;     }
	if ( o[f::no_hb_bb  ].user() && o[f::no_hb_bb ]() && hb_bb() != 0.0 )                                                                  { if ( s ) { ++s->F_no_hb_bb;      } return false; } else if ( s ) { ++s->P_no_hb_bb;       }
	if ( o[f::no_hb_bb  ].user() && o[f::no_hb_bb ]() && (ss1_=='P'||ss2_=='P') )                                                         { if ( s ) { ++s->F_no_hb_bb;      } return false; } else if ( s ) { ++s->P_no_hb_bb;       }
	if ( o[f::ss1       ].user() && o[f::ss1      ]()[0] != ss1()    )                                                                    { if ( s ) { ++s->F_ss1;           } return false; } else if ( s ) { ++s->P_ss1;            }
	if ( o[f::ss2       ].user() && o[f::ss2      ]()[0] != ss2()    )                                                                    { if ( s ) { ++s->F_ss2;           } return false; } else if ( s ) { ++s->P_ss2;            }
	if ( o[f::dssp1     ].user() && o[f::dssp1    ]()[0] != dssp1()  )                                                                    { if ( s ) { ++s->F_dssp1;         } return false; } else if ( s ) { ++s->P_dssp1;          }
	if ( o[f::dssp2     ].user() && o[f::dssp2    ]()[0] != dssp2()  )                                                                    { if ( s ) { ++s->F_dssp2;         } return false; } else if ( s ) { ++s->P_dssp2;          }
	if ( o[f::aa1       ].user() && o[f::aa1      ]()[0] != aa1_     )                                                                    { if ( s ) { ++s->F_aa1;           } return false; } else if ( s ) { ++s->P_aa1;            }
	if ( o[f::aa2       ].user() && o[f::aa2      ]()[0] != aa2_     )                                                                    { if ( s ) { ++s->F_aa2;           } return false; } else if ( s ) { ++s->P_aa2;            }
	if ( o[f::score     ].user() && o[f::score    ]() <  trust*score()     )                                                              { if ( s ) { ++s->F_score;         } return false; } else if ( s ) { ++s->P_score;          }
	if ( o[f::faatr     ].user() && o[f::faatr    ]() <  trust*fa_atr()    )                                                              { if ( s ) { ++s->F_faatr;         } return false; } else if ( s ) { ++s->P_faatr;          }
	if ( o[f::hb_sc     ].user() && o[f::hb_sc    ]() <  trust*hb_sc()     )                                                              { if ( s ) { ++s->F_hb_sc;         } return false; } else if ( s ) { ++s->P_hb_sc;          }
	if ( o[f::hb_bb_sc  ].user() && o[f::hb_bb_sc ]() <  trust*hb_bb_sc()  )                                                              { if ( s ) { ++s->F_hb_bb_sc;      } return false; } else if ( s ) { ++s->P_hb_bb_sc;       }
	if ( o[f::hb_bb     ].user() && o[f::hb_bb    ]() <  trust*hb_bb()     )                                                              { if ( s ) { ++s->F_hb_bb;         } return false; } else if ( s ) { ++s->P_hb_bb;          }
	if ( o[f::coorderr  ].user() && (o[f::coorderr ]() < bfac1() || o[f::coorderr ]() < bfac2() ) )                                       { if ( s ) { ++s->F_coorderr;      } return false; } else if ( s ) { ++s->P_coorderr;       }
	if ( o[f::uniformfrag]() && ( ss1_!=ss2_ || nbrs2()!=1 || !this->check_lig_code("frg") ) )                                            { if ( s ) { ++s->F_uniformfrag;   } return false; } else if ( s ) { ++s->P_uniformfrag;    }
	Real x=uint16_to_real(rt6_[1],-CSIZE,CSIZE), y=uint16_to_real(rt6_[2],-CSIZE,CSIZE), z=uint16_to_real(rt6_[3],-CSIZE,CSIZE);
	Real dis2 = x*x+y*y+z*z;
	if ( o[f::mindist2 ].user() && o[f::mindist2 ]() > dis2 )                                                                             { if ( s ) { ++s->F_mindist2;      } return false; } else if ( s ) { ++s->P_mindist2;       }
	if ( o[f::maxdist2 ].user() && o[f::maxdist2 ]() < dis2 )                                                                             { if ( s ) { ++s->F_maxdist2;      } return false; } else if ( s ) { ++s->P_maxdist2;       }

	if ( o[f::faatr_or_hbbb].user()       &&
			o[f::faatr_or_hbbb] < trust*fa_atr()   &&
			o[f::faatr_or_hbbb] < trust*hb_sc()    &&
			o[f::faatr_or_hbbb] < trust*hb_bb_sc() &&
			o[f::faatr_or_hbbb] < trust*hb_bb()     )                                                                                          { if ( s ) { ++s->F_faatr_or_hbbb; } return false; } else if ( s ) { ++s->P_faatr_or_hbbb;  }

	if ( o[f::faatr_or_hb].user()       &&
			o[f::faatr_or_hb] < trust*fa_atr()   &&
			o[f::faatr_or_hb] < trust*hb_sc()     )                                                                                            { if ( s ) { ++s->F_faatr_or_hb;    } return false; } else if ( s ) { ++s->P_faatr_or_hb;    }
	// if( o[mh::lg1   ].user() && o[mh::lg1   ]() != lg1x_     )                                                                        { if(s) { ++s->F_FULLME;         } return false; } else if(s){ ++s->P_FULLME;         }
	// if( o[mh::lg2   ].user() && o[mh::lg2   ]() != lg2x_     )                                                                        { if(s) { ++s->F_FULLME;         } return false; } else if(s){ ++s->P_FULLME;         }

	if ( o[f::noloops ]() && (ss1()=='L' || ss2() =='L')  )                                                                               { if ( s ) { ++s->F_noloops;        } return false; } else if ( s ) { ++s->P_noloops;        }
	if ( o[f::oneloop ]() && (ss1()=='L' && ss2() =='L')  )                                                                               { if ( s ) { ++s->F_oneloop;        } return false; } else if ( s ) { ++s->P_oneloop;        }
	if ( o[f::nodisulf]() && (aa1_=='C' && aa2_=='C')  )                                                                                  { if ( s ) { ++s->F_nodisulf;       } return false; } else if ( s ) { ++s->P_nodisulf;       }

	size_t const N(string::npos);
	if ( o[f::    restype1   ].user() &&  o[f::    restype1   ]().find(aa1_)==N                                           )               { if ( s ) { ++s->F_restype1;       } return false; } else if ( s ) { ++s->P_restype1;       }
	if ( o[f::    restype2   ].user() &&  o[f::    restype2   ]().find(aa2_)==N                                           )               { if ( s ) { ++s->F_restype2;       } return false; } else if ( s ) { ++s->P_restype2;       }
	if ( o[f::    restype    ].user() && (o[f::    restype    ]().find(aa1_)==N || o[f::    restype    ]().find(aa2_)==N) )               { if ( s ) { ++s->F_restype;        } return false; } else if ( s ) { ++s->P_restype;        }
	if ( o[f::    restype_one].user() && (o[f::    restype_one]().find(aa1_)==N && o[f::    restype_one]().find(aa2_)==N) )               { if ( s ) { ++s->F_restype_one;    } return false; } else if ( s ) { ++s->P_restype_one;    }
	if ( o[f::not_restype    ].user() && (o[f::not_restype    ]().find(aa1_)!=N && o[f::not_restype    ]().find(aa2_)!=N) )               { if ( s ) { ++s->F_not_restype;    } return false; } else if ( s ) { ++s->P_not_restype;    }
	if ( o[f::not_restype_one].user() && (o[f::not_restype_one]().find(aa1_)!=N || o[f::not_restype_one]().find(aa2_)!=N) )               { if ( s ) { ++s->F_not_restype_one;} return false; } else if ( s ) { ++s->P_not_restype_one;}
	return true;
}
Real ResPairMotif::dump_aligned_motif( ostream & out, Pose const & paln1, Size const & ir, Pose const & paln2, Size const & jr, Size & /* atomno */, int const & tag, Xforms const & xforms ) const {
	using namespace core::id;
	using namespace core::pose::motif;
	Pose pose;
	fill_pose_with_motif(pose,tag,tag);
	// Real const motif_align_rms = align_motif_pose_super(pose,paln1,ir,paln2,jr);
	Real const motif_align_rms = align_motif_pose(pose,paln1,ir,paln2,jr,type());
	core::id::AtomID_Mask const mask = get_motif_atom_mask(pose,type(),/*with_Hpol=*/true);
	utility::vector1<Size> resnums(pose.n_residue(),tag);
	BOOST_FOREACH ( Xform const & x,xforms ) {
		xform_pose(pose,x);
		//core::io::pdb::dump_pdb(pose,out,mask,atomno,string_of(tag),'~'/*,resnums*/); JAB XRW
		core::io::pdb::dump_pdb(pose, out, mask, string_of(tag));
		xform_pose(pose,~x);
	}
	return motif_align_rms;
}
Real ResPairMotif::dump_aligned_motif( std::string const & fn, Pose const & paln1, Size const & ir, Pose const & paln2, Size const & jr, Size & atomno, int const & tag, Xforms const & xforms ) const {
	utility::io::ozstream out(fn);
	Real const v = dump_aligned_motif(out,paln1,ir,paln2,jr,atomno,tag,xforms);
	out.close();
	return v;
}
Real ResPairMotif::dump_aligned_motif( ostream & out, Pose const & paln, Size const & ir, Size const & jr, Size & atomno, int const & tag, Xforms const & xforms ) const {
	return dump_aligned_motif(out,paln,ir,paln,jr,atomno,tag,xforms);
}
Real ResPairMotif::dump_aligned_motif( std::string const & fn, Pose const & paln, Size const & ir, Size const & jr, Size & atomno, int const & tag, Xforms const & xforms ) const {
	return dump_aligned_motif(fn,paln,ir,paln,jr,atomno,tag,xforms);
}
Real ResPairMotif::dump_aligned_motif(std::ostream     & out, Pose const & pose1, Size const & ir, Pose const & pose2, Size const & jr, int const & num, numeric::Xforms const & xforms) const{ Size tmp=0; return dump_aligned_motif(out,pose1,ir,pose2,jr,tmp,num,xforms); }
Real ResPairMotif::dump_aligned_motif(std::ostream     & out, Pose const & pose , Size const & ir              , Size const & jr, int const & num, numeric::Xforms const & xforms) const{ Size tmp=0; return dump_aligned_motif(out,pose ,ir,      jr,tmp,num,xforms); }
Real ResPairMotif::dump_aligned_motif(std::string const & fn, Pose const & pose1, Size const & ir, Pose const & pose2, Size const & jr, int const & num, numeric::Xforms const & xforms) const{ Size tmp=0; return dump_aligned_motif( fn,pose1,ir,pose2,jr,tmp,num,xforms); }
Real ResPairMotif::dump_aligned_motif(std::string const & fn, Pose const & pose , Size const & ir              , Size const & jr, int const & num, numeric::Xforms const & xforms) const{ Size tmp=0; return dump_aligned_motif( fn,pose ,ir,      jr,tmp,num,xforms); }

void ResPairMotif::fill_pose_with_motif( Pose & pose, int const & ir, int const & jr ) const {
	pose.clear();
	string const & name1( core::chemical::name_from_aa(core::chemical::aa_from_oneletter_code(aa1())) );
	string const & name2( core::chemical::name_from_aa(core::chemical::aa_from_oneletter_code(aa2())) );
	core::chemical::ResidueType const & rest1 = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard")->name_map(name1);
	core::chemical::ResidueType const & rest2 = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard")->name_map(name2);
	core::conformation::ResidueOP res1 = core::conformation::ResidueFactory::create_residue(rest1);
	core::conformation::ResidueOP res2 = core::conformation::ResidueFactory::create_residue(rest2);

	pose.append_residue_by_jump(*res1,1);
	pose.append_residue_by_jump(*res2,1);
	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( pose, true ) );
	pdb_info->number(1,ir);
	pdb_info->number(2,jr);
	for ( Size i=1; i<=pose.n_residue(); ++i ) pdb_info->chain (i,'~');
	pose.pdb_info( pdb_info );
	for ( Size i=1; i<=pose.residue(1).nchi(); ++i ) { pose.set_chi(i,1,uint8_to_real(chi1_[i],-180.0,180.0)); }
	for ( Size i=1; i<=pose.residue(2).nchi(); ++i ) { pose.set_chi(i,2,uint8_to_real(chi2_[i],-180.0,180.0)); }
	set_residue_pair_xform(xform(),pose,1,2,type());
}
core::pose::Pose ResPairMotif::get_pose(int const & ir, int const & jr) const {
	core::pose::Pose pose;
	fill_pose_with_motif(pose,ir,jr);
	return pose;
}
void ResPairMotif::dump_pdb(
	string const & fname,
	string const & tag
) const {
	utility::io::ozstream out(fname!=""?fname:this->tag()+".pdb");
	this->dump_pdb(out,Xform(),tag);
	out.close();
}
void ResPairMotif::dump_pdb( ostream & out, Xform const & x, string tag ) const {
	Pose pose;
	fill_pose_with_motif(pose);

	core::id::AtomID_Mask const mask = get_motif_atom_mask(pose,type(),/*with_Hpol=*/true);


	Xform xl;
	// xl.from_four_points( pose.residue(1).xyz(2),pose.residue(1).xyz(1),pose.residue(1).xyz(2),pose.residue(1).xyz(3) );
	// xform_pose(pose, x * ~xl );
	if ( x != Xform::identity() ) utility_exit_with_message("not implemeted");
	// cerr << "dump_pdb" << endl;
	if ( tag=="" ) tag = this->tag();
	out << "MODEL " << tag << endl;
	//Size atomno=0;
	//core::io::pdb::dump_pdb(pose,out,mask,atomno,string_of(tag),'~');
	core::io::pdb::dump_pdb(pose, out, mask, string_of(tag));

	out << "ENDMDL" << endl;
}
void ResPairMotif::print_header(ostream & out){
	out << "RPM  type pdbid_lg  pir  pjr AA SS CC bfac1 bfac2 nbr1 nbr2 e cnt   atrsc atrscbb   atrbb hb_scsc hb_bbsc hb_bbbb rt     dx     dy     dz    ex    ey    ez chi chi1 chi2 chi3 chi4 chi5 chi6 chi7 chi8" << endl;
}
std::string ResPairMotif::tag() const {
	return pdb()+"_"+string_of(resi1_)+"_"+string_of(resi2_)+"_"+aa1()+aa2()+"_"+string_of(type() );
}

ostream & operator<<(ostream & out, ResPairMotif const & x){
	Real6 const rt = x.rt();
	out << "RPM " << x.type() << " ";
	out << x.pdb() << " ";
	out << I(4,x.resi1_) <<" ";
	out << I(4,x.resi2_) <<" ";
	out << x.aa1_ <<"";
	out << x.aa2_ <<" ";
	out << x.ss1_ <<"";
	out << x.ss2_ <<" ";
	out << x.chain1_ <<"";
	out << x.chain2_ <<" ";
	out << F(5,3,x.bfac1())  <<" ";
	out << F(5,3,x.bfac2())  <<" ";
	out << F(4,1,x.nbrs1())  <<" ";
	out << F(4,1,x.nbrs2())  <<" ";
	out << "e ";
	out << I(3,x.count())   <<" ";
	out << F(7,3,x.fa_atr())   <<" ";
	out << F(7,3,x.fa_atr_sc_bb())   <<" ";
	out << F(7,3,x.fa_atr_bb())   <<" ";
	out << F(7,3,x.hb_sc())    <<" ";
	out << F(7,3,x.hb_bb_sc()) <<" ";
	out << F(7,3,x.hb_bb())    <<" ";
	out << "rt " << F(6,2,rt[1])<<' '<<F(6,2,rt[2])<<' '<<F(6,2,rt[3])<<' '<<F(5,1,rt[4])<<' '<<F(5,1,rt[5])<<' '<<F(5,1,rt[6]) <<" ";
	out << "chi "
		<< I(4,(int)uint8_to_real(x.chi1_[1],-180.0,180.0))<<' '
		<< I(4,(int)uint8_to_real(x.chi1_[2],-180.0,180.0))<<' '
		<< I(4,(int)uint8_to_real(x.chi1_[3],-180.0,180.0))<<' '
		<< I(4,(int)uint8_to_real(x.chi1_[4],-180.0,180.0))<<' '
		<< I(4,(int)uint8_to_real(x.chi2_[1],-180.0,180.0))<<' '
		<< I(4,(int)uint8_to_real(x.chi2_[2],-180.0,180.0))<<' '
		<< I(4,(int)uint8_to_real(x.chi2_[3],-180.0,180.0))<<' '
		<< I(4,(int)uint8_to_real(x.chi2_[4],-180.0,180.0))   ;
	return out;
}
ostream & operator<<(ostream & out, RM_Type const & x){
	switch(x){
	case RM_Type_NONE : out << "RM_Type_NONE"; break;
	case RM_SC : out << "RM_SC"; break;
	case RM_BB : out << "RM_BB"; break;
	case RM_PH : out << "RM_PH"; break;
	case RM_PO : out << "RM_PO"; break;
	default : utility_exit_with_message("unknown RM_Type "+I(2,x));
	}
	return out;
}
ostream & operator<<(ostream & out, RPM_Type const & x){
	switch(x){
	case RPM_Type_NONE : out << "RPM_Type_NONE"; break;
	case SC_SC : out << "SC_SC"; break;
	case SC_BB : out << "SC_BB"; break;
	case SC_PH : out << "SC_PH"; break;
	case SC_PO : out << "SC_PO"; break;
	case BB_BB : out << "BB_BB"; break;
	case BB_PH : out << "BB_PH"; break;
	case BB_PO : out << "BB_PO"; break;
	case PH_PO : out << "PH_PO"; break;
	default : utility_exit_with_message("unknown RPM_Type "+I(2,x));
	}
	return out;
}

bool write_motifs_binary(ostream & out, ResPairMotifs const & motifs){
	uint64_t const v = MOTIF_VERSION;
	uint64_t const n = motifs.size();
	out.write((char*)&v,sizeof(uint64_t));
	out.write((char*)&n,sizeof(uint64_t));
	for ( ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i ) {
		ResPairMotif const & sm( *i );
		out.write((char*)&sm,sizeof(ResPairMotif));
		if ( !out.good() ) return false;
	}
	return true;
}
bool write_motifs_binary(string const & fname_in, ResPairMotifs const & motifs){
	string fname = fname_in;
	if ( fname.size()<12 || ".rpm.bin.gz"!=fname.substr(fname.size()-11,11) ) fname += ".rpm.bin.gz";
	utility::io::ozstream out(fname,std::ios::out|std::ios::binary);
	if ( !out.good() ) {
		std::cerr << "write_motifs_binary(fname): problem opening motif output file " << fname << endl;
		return false;
	}
	if ( !write_motifs_binary(out,motifs) ) {
		std::cerr << "write_motifs_binary(fname): problem while writing motif file " << fname << endl;
		return false;
	}
	if ( !out.good() ) {
		std::cerr << "write_motifs_binary(fname): problem after writing motif file " << fname << endl;
		out.close();
		return false;
	}
	out.close();
	TR << "wrote motif to " << fname << endl;
	return true;
}
bool read_motifs_binary(std::istream & in, ResPairMotifs & motifs){
	while ( true ) {
		uint64_t n0=motifs.size();
		uint64_t n=123456789;
		uint64_t v=123456789;
		in.read((char*)&v,sizeof(uint64_t));
		if ( !in.good() ) break;
		if ( v != MOTIF_VERSION ) utility_exit_with_message("bad version .rpm.bin.gz: "+string_of(v)+", should be "+string_of(MOTIF_VERSION));
		in.read((char*)&n,sizeof(uint64_t));
		motifs.resize(motifs.size()+n);
		for ( Size i = 1; i <= n; ++i ) {
			in.read((char*)&motifs[n0+i],sizeof(ResPairMotif));
		}
	}
	return true;
}
bool read_motifs_binary(string const & fname, ResPairMotifs & motifs){
	if ( ".rpm.bin.gz"!=fname.substr(fname.size()-11,11) ) {
		std::cerr << "read_motifs_binary(fname): bad extension on input file " << fname << endl;
		return false;
	}
	if ( !utility::file::file_exists(fname) ) utility_exit_with_message("non-existent motif file "+fname);
	utility::io::izstream in(fname,std::ios::in|std::ios::binary);
	if ( !in.good() ) {
		std::cerr << "read_motifs_binary(fname): problem opening motif input file " << fname << endl;
		return false;
	}
	if ( !read_motifs_binary(in,motifs) ) {
		std::cerr << "read_motifs_binary(fname): problem while reading motif file " << fname << endl;
		return false;
	}
	// if(!in.good()){
	//  std::cerr << "problem after reading motif file " << fname << endl;
	//  in.close();
	//  return false;
	// }
	// in.close();
	return true;
}
bool read_motifs_binary(vector1<string> const & fnames, ResPairMotifs & motifs){
	TR << "Nfiles: " << fnames.size() << endl;
	for ( vector1<string>::const_iterator i = fnames.begin(); i != fnames.end(); ++i ) {
		// if( ++count %50 == 0 ){ TR << "... read " << count << endl; }
		if ( fnames.size()<10 ) TR << "reading binary ResPairMotifs " << *i << endl;
		else TR << ".";
		if ( !read_motifs_binary(*i,motifs) ) {
			std::cerr << "read_motifs_binary(fnames): error reading file "+*i << endl;
			if ( !option[basic::options::OptionKeys::mh::ignore_io_errors]() ) utility_exit_with_message("");
		}
	}
	TR << endl;
	return true;
}
void load_motifs(vector1<string> const & fnames, ResPairMotifs & motifs, ResPairMotifsStringMap * map){
	motifs.clear();
	BOOST_FOREACH ( std::string const & fname,fnames ) {
		ResPairMotifs raw_motifs,filt_motifs;
		if ( !read_motifs_binary(fname,raw_motifs) ) utility_exit_with_message("error reading motif files");
		TR << "read " << ((Real)raw_motifs.size())/1000000.0 << " million binary motifs from " << fname << endl;
		if ( option[mh::gen_reverse_motifs_on_load]() ) {
			raw_motifs.add_reverse_motifs();
		}
		if ( option[basic::options::OptionKeys::mh::filter::filter_io]() ) {
			filter_motifs(raw_motifs,filt_motifs);
			TR << "filter: keep " << ((Real)filt_motifs.size())/1000000.0 << " million binary motifs from " << fname << endl;
		} else {
			filt_motifs = raw_motifs;
		}
		if ( option[basic::options::OptionKeys::mh::merge_similar_motifs].user() ) {
			Real cartsize = option[basic::options::OptionKeys::mh::merge_similar_motifs]()[1];
			Real cartresl = option[basic::options::OptionKeys::mh::merge_similar_motifs]()[2];
			Real anglresl = option[basic::options::OptionKeys::mh::merge_similar_motifs]()[3];
			filt_motifs.filter_structurally_similar_motifs(cartsize,cartresl,anglresl);
			TR << "merge:  keep " << ((Real)filt_motifs.size())/1000000.0 << " million binary motifs from " << fname << endl;
		}
		motifs.insert(motifs.begin(),filt_motifs.begin(),filt_motifs.end());
		if ( map ) map->insert(make_pair(fname,filt_motifs));
	}
}
void load_motifs(string const & fname, ResPairMotifs & motifs, ResPairMotifsStringMap * map){
	load_motifs(utility::vector1<string>(1,fname),motifs,map);
}

struct MyHash {
	uint64_t operator()(ResPairMotif const & rpm) const {
		uint64_t const *x = (uint64_t const *)(&rpm);
		return x[0] ^ x[1] ^ x[2] ^ x[3] ^ x[4] ^ x[5];
	}
};
struct MyPred {
	bool operator()(ResPairMotif const & a, ResPairMotif const & b) const {
		return 0 == memcmp(&a,&b,48);
	}
};
typedef boost::unordered_set<ResPairMotif,MyHash,MyPred> MotifSet;
void ResPairMotifs::filter_structurally_identical_motifs(){
	utility_exit_with_message("not impl");
	MotifSet motifset;
	BOOST_FOREACH ( ResPairMotif const & rpm, *this ) {
		ResPairMotif m(rpm);
		if ( m.is_reversible() && m.resi1()>m.resi2() ) m.reverse_in_place();
		motifset.insert(m);
	}
	clear();
	reserve(motifset.size());
	BOOST_FOREACH ( ResPairMotif const & rpm, motifset ) push_back(rpm);
}
void ResPairMotifs::filter_structurally_similar_motifs(Real /*cart_size*/, Real /*cart_resl*/, Real /*ang_resl*/){
	utility_exit_with_message("not impl");
	// for(int iofst = 0; iofst < 30; ++iofst){
	//  Real6 shift(0);
	//  if(iofst>0){ shift[1]=uniform(); shift[2]=uniform(); shift[3]=uniform(); shift[4]=uniform(); shift[5]=uniform(); shift[6]=0; }
	//  MotifHash mh(cart_size + (iofst? cart_resl/2.0 : 0.0)  ,cart_resl,ang_resl);
	//  ResPairMotifs all_keepers;
	//  BOOST_FOREACH(ResPairMotif const & rpm, *this) mh.add_motif_shift(rpm,shift);
	//  for(MotifHash::KeySet::const_iterator i = mh.begin(); i != mh.end(); ++i){
	//   ResPairMotifs this_bin,keepers;
	//   mh.find_motifs(*i,this_bin);
	//   BOOST_FOREACH(ResPairMotif const & m, this_bin){
	//    bool seenit = false;
	//    BOOST_FOREACH(ResPairMotif & n, keepers){
	//     bool seenthis = ( m.aa1()==n.aa1() && m.aa2()==n.aa2() && m.ss1()==n.ss1() && m.ss2()==n.ss2() );
	//     if( min(std::fabs(m.chi11()-n.chi11()),360.0-std::fabs(m.chi11()-n.chi11())) > ang_resl*1.5 ) seenthis = false;
	//     if( min(std::fabs(m.chi21()-n.chi21()),360.0-std::fabs(m.chi21()-n.chi21())) > ang_resl*1.5 ) seenthis = false;
	//     if( min(std::fabs(m.chi12()-n.chi12()),360.0-std::fabs(m.chi12()-n.chi12())) > ang_resl*2.0 ) seenthis = false;
	//     if( min(std::fabs(m.chi22()-n.chi22()),360.0-std::fabs(m.chi22()-n.chi22())) > ang_resl*2.0 ) seenthis = false;
	//     if( min(std::fabs(m.chi13()-n.chi13()),360.0-std::fabs(m.chi13()-n.chi13())) > ang_resl*2.5 ) seenthis = false;
	//     if( min(std::fabs(m.chi23()-n.chi23()),360.0-std::fabs(m.chi23()-n.chi23())) > ang_resl*2.5 ) seenthis = false;
	//     if( min(std::fabs(m.chi14()-n.chi14()),360.0-std::fabs(m.chi14()-n.chi14())) > ang_resl*3.0 ) seenthis = false;
	//     if( min(std::fabs(m.chi24()-n.chi24()),360.0-std::fabs(m.chi24()-n.chi24())) > ang_resl*3.0 ) seenthis = false;
	//     if(seenthis){
	//      seenit = true;
	//      n.addcount();
	//      break;
	//     }
	//    }
	//    if(!seenit) keepers.push_back(m);
	//   }
	//   all_keepers.insert(all_keepers.end(),keepers.begin(),keepers.end());
	//  }
	//  *this = all_keepers;
	//  cout << iofst << " " << size() << endl;
	//  }
}
void ResPairMotifs::add_reverse_motifs(){
	this->reserve(2*this->size());
	ResPairMotifs rev;
	rev.reserve(this->size());
	BOOST_FOREACH ( ResPairMotif const & m,*this ) {
		if ( m.type1()==m.type2() && !m.check_lig_code("frg") ) {
			rev.push_back( m.reversed() );
		}
	}
	this->insert(this->end(),rev.begin(),rev.end());
	TR << "ResPairMotifs::added " << rev.size() << " reverse motifs" << endl;
}


ResPairMotifMetaBinner::ResPairMotifMetaBinner():
	sep_lj_  (option[basic::options::OptionKeys::mh::harvest::sep_lj  ]()),
	sep_hb_  (option[basic::options::OptionKeys::mh::harvest::sep_hb  ]()),
	sep_nbrs_(option[basic::options::OptionKeys::mh::harvest::sep_nbrs]()),
	sep_bfac_(option[basic::options::OptionKeys::mh::harvest::sep_bfac]()),
	sep_dist_(option[basic::options::OptionKeys::mh::harvest::sep_dist]()),
	sep_aa_  (option[basic::options::OptionKeys::mh::harvest::sep_aa  ]()),
	sep_aa1_ (option[basic::options::OptionKeys::mh::harvest::sep_aa1 ]()),
	sep_aa2_ (option[basic::options::OptionKeys::mh::harvest::sep_aa2 ]()),
	sep_ss_  (option[basic::options::OptionKeys::mh::harvest::sep_ss  ]()),
	sep_dssp_(option[basic::options::OptionKeys::mh::harvest::sep_dssp]())
{
	if ( sep_dssp_||sep_lj_.size()||sep_hb_.size()||sep_nbrs_.size()||sep_bfac_.size()||sep_dist_.size() ) {
		utility_exit_with_message("ResPairMotifMetaBinner not implemented yet: sep_dssp sep_lj sep_hb sep_nbrs sep_bfac sep_dist_");
	}
}
std::string
ResPairMotifMetaBinner::motif_bin_label(ResPairMotif const & r) const {
	string lab="";
	if     ( sep_aa1_ ) lab += "_" + string_of(r.aa1())          ;
	else if ( sep_aa2_ ) lab += "_" + string_of(r.aa2())          ;
	else if ( sep_aa_  ) lab += "_" + string_of(r.aa1()) + r.aa2();
	if     ( sep_ss_  ) lab += "_" + string_of(r.ss1()) + r.ss2();
	if ( lab.size() ) lab = lab.substr(1);
	return lab;
}
ResPairMotifMetaBinner::Key ResPairMotifMetaBinner::motif_bin_hash(ResPairMotif const & r) const {
	Key k=0;
	if     ( sep_aa1_ ) k ^= ((Key)r.aa1())<< 0              ;
	else if ( sep_aa2_ ) k ^=                      ((Key)r.aa2())<< 8;
	else if ( sep_aa_  ) k ^= ((Key)r.aa1())<< 0 ^ ((Key)r.aa2())<< 8;
	if     ( sep_ss_  ) k ^= ((Key)r.ss1())<<16 ^ ((Key)r.ss2())<<24;
	return k;
}
std::string ResPairMotifMetaBinner::hash_to_labal(ResPairMotifMetaBinner::Key const & k) const{
	string lab="";
	if     ( sep_aa1_ ) lab += "_" + string_of((char)(k>> 0))          ;
	else if ( sep_aa2_ ) lab += "_" + string_of(                 (char)(k>> 8));
	else if ( sep_aa_  ) lab += "_" + string_of((char)(k>> 0)) + (char)(k>> 8);
	if     ( sep_ss_  ) lab += "_" + string_of((char)(k>>16)) + (char)(k>>24);
	if ( lab.size() ) lab = lab.substr(1);
	return lab;
}

ResPairMotifMetaBinner::Key
ResPairMotifMetaBinner::bin0_of_real(Real val, Reals const & breaks) const {
	Key k=0;
	BOOST_FOREACH ( Real const & b, breaks ) {
		if ( val < b ) break;
		++k;
	}
	return k;
}
std::ostream & operator<<(std::ostream & out, XformScoreMap const & x){
	out << "      LAB XformScore  hcsize  hcresl  haresl   Npossible     tot_sc   num_bins  scPb    qt1   qt2   qt3   qt4   qt5   qt6   qt7   qt8   qt9  qt10  qt11  qt12  qt13  qt14  qt15  qt16  qt17  qt18  qt19  qt20" << endl;
	ResPairMotifMetaBinner binner;
	BOOST_FOREACH ( XformScoreMap::value_type const & p,x ) {
		out << RJ(9,binner.hash_to_labal(p.first)) << " " << *p.second << endl;
	}
	return out;
}

bool MotifHit::operator<(MotifHit const & other) const{
	// !!!!!!!!!!!!!!!!!!! HERE BE DRAGONS !!!!!!!!!!!!!!!!!!!!!!
	return memcmp(this,&other,sizeof(MotifHit)-16) < 0;

}
bool MotifHit::operator==(MotifHit const & other) const{
	// !!!!!!!!!!!!!!!!!!! HERE BE DRAGONS !!!!!!!!!!!!!!!!!!!!!!
	return memcmp(this,&other,sizeof(MotifHit)-16) == 0;
}

void MotifHits::push_back(MotifHit const & h){
	utility::vector1<MotifHit>::push_back(h);
	back().index = size();
	resmap1_.insert(make_pair((uint32_t)back().residue1,(uint32_t)size()));
	resmap2_.insert(make_pair((uint32_t)back().residue2,(uint32_t)size()));
}
void MotifHits::push_back_raw(MotifHit const & h){
	utility::vector1<MotifHit>::push_back(h);
}
void MotifHits::filter_redundant(){
	std::set<MotifHit> hitset(begin(),end());
	clear();
	std::copy(hitset.begin(),hitset.end(),std::back_insert_iterator<MotifHits>(*this));
}
struct RMSandEnergyCMP { bool operator()(MotifHit const & a, MotifHit const & b) const {
	return (sqr(a.rms))/std::fabs(a.motif.score())  <  (sqr(b.rms))/std::fabs(b.motif.score());
}};
void MotifHits::sort_by_rms_and_energy(){
	std::sort(begin(),end(),RMSandEnergyCMP());
}
void MotifHits::compute_metrics() const {
	std::set<int> resset;
	std::set<std::pair<int,int> > pairset;
	total_score = 0.0;
	for ( MotifHits::const_iterator i = begin(); i != end(); ++i ) {
		resset.insert(i->residue1);
		resset.insert(i->residue2);
		pairset.insert(std::make_pair(i->residue1,i->residue2));
		total_score += i->motif.score();
	}
	rescover = resset.size();
	paircover = pairset.size();
}
void MotifHits::dump_motifs_pdb( string const & fn, Size & atomno, Xforms const & xforms ) const {
	utility::io::ozstream out(fn);
	int tmp=0;
	dump_motifs_pdb(out,tmp,atomno,xforms);
	out.close();
}
void MotifHits::dump_motifs_pdb( ostream & out, int & count, Size & atomno, Xforms const & xforms ) const {
	// BOOST_FOREACH(MotifHit const & hit,*this){
	//  hit.motif.dump_aligned_motif(out,*hit.pose1,hit.residue1,*hit.pose2,hit.residue2,++count);
	// }
	// return;
	typedef boost::unordered_map<uint64_t,MotifHits> HitResMap;
	HitResMap byres;
	BOOST_FOREACH ( MotifHit const & hit,*this ) {
		byres[((uint64_t)(hit.residue1))    ].push_back_raw(hit);
		byres[((uint64_t)(hit.residue2))<<32].push_back_raw(hit);
	}
	boost::unordered_set<int> dumpedit;
	BOOST_FOREACH ( HitResMap::value_type & v, byres ) {
		v.second.sort_by_rms_and_energy();
		int rescount=0;
		BOOST_FOREACH ( MotifHit const & hit,v.second ) {
			if ( ++rescount > option[mh::dump::max_per_res]() ) if ( hit.motif.type()!=BB_PH && hit.motif.type()!=BB_PO ) continue;
			if ( hit.rms    > option[mh::dump::max_rms    ]() ) continue;
			if ( dumpedit.find(hit.index)!=dumpedit.end() ) continue;
			hit.motif.dump_aligned_motif(out,hit.pose1(),hit.residue1,hit.pose2(),hit.residue2,atomno,++count,xforms);
			dumpedit.insert(hit.index);

			// Size ir(i->residue1);
			// Size jr(i->residue2);
			// ResPairMotif const & m(i->motif);
			// core::pack::rotamer_set::RotamerSetOP rotset1,rotset2;
			// {
			//  Pose pose(pose1);
			//  core::scoring::ScoreFunction dummy_sfxn;
			//  dummy_sfxn( pose );
			//  core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
			//  dummy_task->initialize_from_command_line();
			//  dummy_task->nonconst_residue_task( ir ).restrict_to_repacking();
			//  dummy_task->nonconst_residue_task( ir ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
			//  dummy_task->nonconst_residue_task( ir ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
			//  core::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
			//  core::pack::rotamer_set::RotamerSetFactory rsf;
			//  rotset1 = rsf.create_rotamer_set( pose.residue( ir ) );
			//  rotset1->set_resid( ir );
			//  rotset1->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
			// }
			// {
			//  Pose pose(pose2);
			//  core::scoring::ScoreFunction dummy_sfxn;
			//  dummy_sfxn( pose );
			//  core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
			//  dummy_task->initialize_from_command_line();
			//  dummy_task->nonconst_residue_task( jr ).restrict_to_repacking();
			//  dummy_task->nonconst_residue_task( jr ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
			//  dummy_task->nonconst_residue_task( jr ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
			//  core::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
			//  core::pack::rotamer_set::RotamerSetFactory rsf;
			//  rotset2 = rsf.create_rotamer_set( pose.residue( jr ) );
			//  rotset2->set_resid( jr );
			//  rotset2->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
			// }
		}

	}
}
void MotifHits::dump_motifs_pdb( string const & fname, Pose const & pose, Size & atomno, numeric::Xforms const & xforms ) const {
	if ( xforms.size()!=1 || xforms.front()!=Xform::identity() ) utility_exit_with_message("not implemented");
	utility::io::ozstream pdbout(fname);
	pdbout << "MODEL MAIN" << endl;
	pose.dump_pdb(pdbout); // TODO maintain atomno!
	pdbout << "ENDMDL" << endl;
	int tmp=0;
	this->dump_motifs_pdb(pdbout,tmp,atomno,xforms);
	pdbout.close();
}

int MotifHits::stat_motifs( MotifHitStats & stats) const {
	using std::map;
	using std::pair;
	using std::make_pair;
	Real hh=0,he=0,hl=0,ee=0,el=0,ll=0,pp=0,totscore=0,totsqstscore=0,count=0;
	Reals scores;
	map<pair<int,int>,Real> pairscores;
	map<int,Real> resscores;
	////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	MotifHits const & hits(*this);
	BOOST_FOREACH ( MotifHit const & hit,hits ) {
		Size ir = hit.residue1;
		Size jr = hit.residue2;
		Pose const & pose1(hit.pose1());
		Pose const & pose2(hit.pose2());
		char aa1 = pose1.residue(ir).name1(), ss1=pose1.secstruct(ir);
		char aa2 = pose2.residue(jr).name1(), ss2=pose2.secstruct(jr);
		// Real6 rt = get_residue_pair_rt6(pose1,ir,pose2,jr);
		Real score=1.0;
		ResPairMotif const & sm(hit.motif);
		Real err = 0;//rt6_rt6_bb_dis2(rt,sm.rt());
		score *= exp( -2.0 * err );
		if ( sm.aa1()==aa1 && sm.aa2()==aa2 ) score *= 1.00;
		else                                 score *= 0.50;
		if ( sm.ss1()==ss1 && sm.ss2()==ss2 ) score *= 1.00;
		else                                 score *= 0.30;
		totscore += score;
		totsqstscore += sqrt(score);
		scores.push_back(score);
		pairscores[make_pair(ir,jr)] += score;
		resscores [ir        ] += score;
		resscores [jr+1000000] += score;
		count += 1.0;
		if ( (sm.dssp1()=='P' && sm.dssp2()=='P') ) pp += 1.0;
		else if ( (sm.ss1()=='H' && sm.ss2()=='H') /*|| (sm.ss1()=='H' && sm.ss2()=='H')*/ ) hh += 1.0;
		else if ( (sm.ss1()=='H' && sm.ss2()=='E') || (sm.ss1()=='E' && sm.ss2()=='H') ) he += 1.0;
		else if ( (sm.ss1()=='H' && sm.ss2()=='L') || (sm.ss1()=='L' && sm.ss2()=='H') ) hl += 1.0;
		else if ( (sm.ss1()=='E' && sm.ss2()=='E') /*|| (sm.ss1()=='E' && sm.ss2()=='E')*/ ) ee += 1.0;
		else if ( (sm.ss1()=='E' && sm.ss2()=='L') || (sm.ss1()=='L' && sm.ss2()=='E') ) el += 1.0;
		else if ( (sm.ss1()=='L' && sm.ss2()=='L') /*|| (sm.ss1()=='L' && sm.ss2()=='L')*/ ) ll += 1.0;
		else std::cerr << "BAD SS IN MOTIF " << sm.ss1() << " " << sm.ss2() << endl;
	}

	if ( stats.find("M_NUM"    )==stats.end() ) stats["M_NUM"    ] = 0;
	// if( stats.find("M_xpsc"   )==stats.end()) stats["M_xpsc"   ] = 0;
	// if( stats.find("M_xpscrt" )==stats.end()) stats["M_xpscrt" ] = 0;
	// if( stats.find("M_PP"     )==stats.end()) stats["M_PP"     ] = 0;
	if ( stats.find("M_HH"     )==stats.end() ) stats["M_HH"     ] = 0;
	if ( stats.find("M_HE"     )==stats.end() ) stats["M_HE"     ] = 0;
	// if( stats.find("M_HL"     )==stats.end()) stats["M_HL"     ] = 0;
	if ( stats.find("M_EE"     )==stats.end() ) stats["M_EE"     ] = 0;
	// if( stats.find("M_EL"     )==stats.end()) stats["M_EL"     ] = 0;
	// if( stats.find("M_LL"     )==stats.end()) stats["M_LL"     ] = 0;

	stats["M_NUM"     ] += count;
	// stats["M_xpsc"    ] += totscore;
	// stats["M_xpscrt" ] += totsqstscore;
	// stats["M_PP" ] += pp;
	stats["M_HH"      ] += hh;
	stats["M_HE"      ] += he;
	// stats["M_HL"      ] += hl;
	stats["M_EE"      ] += ee;
	// stats["M_EL"      ] += el;
	// stats["M_LL"      ] += ll;

	// cout <<endl<< "stat_matching_motifs " << count << " " << hits.size() << " " << stats["M_NUM"     ] << endl<<endl;

	return hits.size();
}
string MotifHits::get_resfile(bool restrict, std::set<Size> & resi_in_resfile) const {
	string resfile = restrict ? "NATRO\nstart\n" : "PIKAA ACDEFGHIKLMNPQRSTVWY\nstart\n";
	// string resfile = restrict ? "NATRO\nstart\n" : "PIKAA AEFIKLMQRVWY\nstart\n";
	// string resfile = restrict ? "NATRO\nstart\n" : "PIKAA AILMV\nstart\n";
	if ( size()==0 ) return resfile;
	MotifHits const & hits(*this);
	PoseCOP pose1op ((at(1).pose1ptr()));
	PoseCOP pose2op ((at(1).pose2ptr()));
	Pose const &  pose1   ((at(1).pose1()));
	Pose const &  pose2   ((at(1).pose2()));
	Size nres1=pose1.n_residue(), nres2=pose2.n_residue();
	if ( core::pose::symmetry::is_symmetric(pose1) ) nres1 = core::pose::symmetry::symmetry_info(pose1)->num_independent_residues();
	if ( core::pose::symmetry::is_symmetric(pose2) ) nres2 = core::pose::symmetry::symmetry_info(pose2)->num_independent_residues();
	bool sameseq = pose1.annotated_sequence()==pose2.annotated_sequence();
	vector1<std::map<char,float> > resdat1(nres1),resdat2(nres2);
	BOOST_FOREACH ( MotifHit const & hit,hits ) {
		PoseCOP thispose1op((hit.pose1ptr()));
		PoseCOP thispose2op((hit.pose2ptr()));
		Pose const &  thispose1  ((hit.pose1()));
		Pose const &  thispose2  ((hit.pose2()));
		bool thissameseq = thispose1.annotated_sequence()==thispose2.annotated_sequence();
		if ( thispose1op != pose1op || thispose2op != pose2op ) {
			cout << "FIXME mismatched MotifHit poses!  " << pose1.n_residue() << " " << pose2.n_residue() << " " << thispose1.n_residue() << " " << thispose2.n_residue() << endl;
			if ( thispose1op           != thispose1op       ) utility_exit_with_message("fixme");
			if ( thispose1.n_residue() != pose1.n_residue() ) utility_exit_with_message("fixme");
		}
		// ResPairMotif const & sm(hit.motif);
		// cout << "get_resfile motif: " << hit << endl;
		Size const ir = hit.residue1, irsym = (ir-1)%nres1+1;
		Size const jr = hit.residue2, jrsym = (jr-1)%nres2+1;
		char const paa1 = thispose1.residue(ir).name1(), pss1=thispose1.secstruct(ir);
		char const paa2 = thispose2.residue(jr).name1(), pss2=thispose2.secstruct(jr);
		char const aa1 = hit.motif.aa1(), ss1=hit.motif.ss1();
		char const aa2 = hit.motif.aa2(), ss2=hit.motif.ss2();
		// cout << "aa " << aa1 << " " << aa2 << " ss " << ss1 << " "  << ss2 << endl;
		// Real6 rt = get_residue_pair_rt6(thispose1,ir,thispose2,jr);
		Real score=1.0;
		Real err = 0;//rt6_rt6_bb_dis2(rt,sm.rt());
		score *= exp( -2.0 * err );
		if ( paa1==aa1 && paa2==aa2 ) score *= 1.00;
		else                         score *= 0.50;
		if ( pss1==ss1 && pss2==ss2 ) score *= 1.00;
		else                         score *= 0.30;
		if ( score >= option[mh::dump::resfile_min_pair_score]() ) {
			if ( resdat1[irsym].find(aa1)==resdat1[irsym].end() ) resdat1[irsym][aa1]=0;
			resdat1[irsym][aa1] += score;
			std::map<char,float> & tmp((thissameseq?resdat1:resdat2)[jrsym]);
			if ( tmp.find(aa1)==tmp.end() ) tmp[aa1]=0;
			tmp[aa2] += score;
		}
	}
	typedef std::pair<char,float> PCF;
	for ( Size ir = 1; ir <= nres1; ++ir ) {
		char const chain1 = utility::UPPERCASE_ALPHANUMERICS[(pose1.chain(ir)-1)%36];
		string aas;
		BOOST_FOREACH ( PCF const & p,resdat1[ir] ) {
			if ( p.second >= option[mh::dump::resfile_min_tot_score]() ) {
				aas += string_of(p.first);
			}
		}
		if ( aas.size()>0 ) {
			if ( aas[0]!='A' ) aas = 'A'+aas;
			resfile += I(4,ir) + " "+chain1+" PIKAA " + aas + "\n";
			resi_in_resfile.insert(ir);
		}
	}
	if ( !sameseq ) {
		for ( Size jr = 1; jr <= nres2; ++jr ) {
			char const chain2 = utility::UPPERCASE_ALPHANUMERICS[(pose1.conformation().num_chains()+pose2.chain(jr)-1)%36];
			string aas;
			BOOST_FOREACH ( PCF const & p,resdat2[jr] ) {
				if ( p.second >= option[mh::dump::resfile_min_tot_score]() ) {
					aas += string_of(p.first);
				}
			}
			if ( aas.size()>0 ) {
				if ( aas[0]!='A' ) aas = 'A'+aas;
				resfile += I(4,jr+nres1) + " "+chain2+" PIKAA " + aas + "\n";
				resi_in_resfile.insert(jr);
			}
		}
	}
	return resfile;
}
string MotifHits::get_resfile(bool restrict) const {
	std::set<Size> resi_in_resfile_dummy;
	return get_resfile(restrict,resi_in_resfile_dummy);
}
void MotifHits::dump_resfile(std::string fn) const {
	utility::io::ozstream out(fn);
	out << get_resfile();
	out.close();
}
ostream & operator<<(ostream & out, MotifHit  const & h){
	out << I(4,h.residue1) << " " << I(4,h.residue2) << " " << F(7,3,h.score1) << " " << F(7,3,h.score2) << " " << h.motif;
	return out;
}
ostream & operator<<(ostream & out, MotifHits const & h){
	h.compute_metrics();
	out << "MotifHits" << " score: " << h.total_score  << ", num: " << h.size() << ", rescover: " << h.rescover << ", paircover: " << h.paircover;
	for ( MotifHits::const_iterator i = h.begin(); i != h.end(); ++i ) {
		out << endl << "    HIT " << *i;
	}
	return out;
}
Size MotifHits::num_hits1(Size ir) const { return resmap1_.count((uint32_t)ir); }
Size MotifHits::num_hits2(Size ir) const { return resmap2_.count((uint32_t)ir); }
Size MotifHits::num_hits (Size ir) const { return resmap1_.count((uint32_t)ir) + resmap2_.count((uint32_t)ir); }
ResPairMotifPtrs MotifHits::hits1(Size ir) const {
	ResPairMotifPtrs ret;
	BOOST_FOREACH ( ResMap::value_type const & v, resmap1_.equal_range((uint32_t)ir) ) {
		ret.push_back(&at(v.second).motif);
	}
	return ret;
}
ResPairMotifPtrs MotifHits::hits2(Size ir) const {
	ResPairMotifPtrs ret;
	BOOST_FOREACH ( ResMap::value_type const & v, resmap2_.equal_range((uint32_t)ir) ) {
		ret.push_back(&at(v.second).motif);
	}
	return ret;
}

/************************************************ MotifHash *********************************************/

MotifHash::MotifHash():
	cart_size_ (MOTIF_HASH_CART_SIZE),
	cart_resl_ ( basic::options::option[basic::options::OptionKeys::mh::harvest::hash_cart_resl]() ),
	angle_resl_( basic::options::option[basic::options::OptionKeys::mh::harvest::hash_angle_resl]() ),
	hasher_(
	numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> >(Vec(-cart_size_),Vec( cart_size_)),
	utility::fixedsizearray1<Size,3>(0),
	get_bins(cart_resl_,angle_resl_)
	),
	type_(RPM_Type_NONE)
{
	sanity_check();
}
MotifHash::MotifHash(Real cart_resl, Real ang_resl):
	cart_size_ (MOTIF_HASH_CART_SIZE),
	cart_resl_ (cart_resl ),
	angle_resl_( ang_resl),
	hasher_(
	numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> >(Vec(-cart_size_),Vec( cart_size_)),
	utility::fixedsizearray1<Size,3>(0),
	get_bins(cart_resl_,angle_resl_)
	),
	type_(RPM_Type_NONE)
{
	sanity_check();
}
MotifHash::MotifHash( ResPairMotifs const & motifs ):
	cart_size_ (MOTIF_HASH_CART_SIZE),
	cart_resl_ ( basic::options::option[basic::options::OptionKeys::mh::harvest::hash_cart_resl]() ),
	angle_resl_( basic::options::option[basic::options::OptionKeys::mh::harvest::hash_angle_resl]() ),
	hasher_(
	numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> >(Vec(-cart_size_),Vec( cart_size_)),
	utility::fixedsizearray1<Size,3>(0),
	get_bins(cart_resl_,angle_resl_)
	),
	type_(RPM_Type_NONE)
{
	for ( ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i ) {
		add_motif(*i);
	}
	sanity_check();
}
RPM_Type MotifHash::type() const {
	runtime_assert(type_!=(uint8_t)RPM_Type_NONE);
	return (RPM_Type)type_;
}
RM_Type MotifHash::type1() const {
	return rpm_type1(type());
}
RM_Type MotifHash::type2() const {
	return rpm_type2(type());
}

void MotifHash::sanity_check() const {
	runtime_assert_msg( cart_resl_ > 0.001 && angle_resl_ > 0.001,"cart & angle resl must be > 0.001");
	runtime_assert_msg( cart_resl_ < 999.0 && angle_resl_ < 999.0,"cart & angle resl must be < 999.0");
	runtime_assert_msg( std::fabs( fmod(180.0/angle_resl_,1.0) ) < 0.000000001, "angle resolution should evenly divide 180.0");
	runtime_assert_msg( cart_size_ == MOTIF_HASH_CART_SIZE, "cart_size != MOTIF_HASH_CART_SIZE");
}

MotifHash::Key MotifHash::bin_index(Real6 const & rt) const {
	Real6 rt6 = rt;
	rt6[4] = fmod(rt6[4],360.0);
	rt6[5] = fmod(rt6[5],360.0);
	rt6[6] = fmod(rt6[6],360.0);
	rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
	rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
	rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
	if ( !hasher_.contains(rt6) ) {
		return std::numeric_limits<Key>::max();
		// AMW: cppcheck flags these statements as following a return
		/*cout << rt6 << endl;
		utility_exit_with_message("out of bounds");*/
	}
	return hasher_.bin_index(rt6);
}
MotifHash::Key MotifHash::bin_index(Motif const & m) const {
	return bin_index(m.rt());
}
void MotifHash::add_motif(Motif const & d, Key const & key){
	if ( type_==RPM_Type_NONE ) type_ = d.type();
	motif_umap_.insert( std::make_pair( key , d ) );
	key_set_.insert(key);
	runtime_assert(d.type()==type_);
	// Motif tmp(d);         // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! should I do this???
	// tmp.reverse_in_place();
	// motif_umap_.insert( std::make_pair( bin_index(tmp.rt()) , tmp ) );
}
void MotifHash::add_motif(Motif const & d){
	Real6 rt( d.rt() );
	if ( std::fabs(rt[1]) > cart_size_ ) return;
	if ( std::fabs(rt[2]) > cart_size_ ) return;
	if ( std::fabs(rt[3]) > cart_size_ ) return;
	Key key = bin_index(rt);
	add_motif(d,key);
}

void MotifHash::find_motifs(Key const & k, ResPairMotifs & results ) const {
	std::pair<MotifMap::const_iterator,MotifMap::const_iterator> range = motif_umap_.equal_range(k);
	for ( MotifMap::const_iterator it = range.first; it != range.second; ++it ) {
		results.push_back( it->second );
	}
}
void MotifHash::find_motifs(Real6 const & rt6, ResPairMotifs & results ) const {
	find_motifs(bin_index(rt6),results);
}
Size MotifHash::count_motifs(Key const & k ) const {
	return motif_umap_.count(k);
}
Size MotifHash::count_motifs(Real6 const & rt6) const {
	return count_motifs(bin_index(rt6));
}
bool MotifHash::check_bounds(Real6 const & rt6) const {
	return std::fabs(rt6[1]) < cart_size_ &&
		std::fabs(rt6[2]) < cart_size_ &&
		std::fabs(rt6[3]) < cart_size_ &&
		0.0 <= rt6[4] && rt6[4] <  360.0 &&
		0.0 <= rt6[5] && rt6[5] <  360.0 &&
		0.0 <= rt6[6] && rt6[6] <= 180.0;
}
void MotifHash::find_motifs_with_radius(Real6 const & rt6, Real radius, vector1<Motif> & results) const {
	using numeric::geometry::hashing::Bin6D;
	if ( !check_bounds(rt6) ) return;
	Real const ang2dis = cart_resl_/angle_resl_;
	Xform const x1(rt6);
	Bin6D const half = hasher_.halfbin6(rt6);
	Bin6D bin6 = hasher_.bin6(rt6);
	utility::fixedsizearray1<int,6> lb,ub;
	for ( numeric::Size i = 1; i <= 6; ++i ) {
		lb[i] = bin6[i];
		ub[i] = bin6[i];
		if ( half[i] && (i==4||i==5||ub[i] < (int)hasher_.dimsizes()[i]-1) ) ++ub[i];
		if ( !half[i] && (i==4||i==5||lb[i] >                0            ) ) --lb[i];
	}
	// cout <<I(5,hasher_.dimsizes()[4])<<I(5,hasher_.dimsizes()[5])<<I(5,hasher_.dimsizes()[6])<<I(5,hasher_.dimsizes()[1])<<I(5,hasher_.dimsizes()[2])<<I(5,hasher_.dimsizes()[3])<<endl;
	// cout <<I(5,lb[4])<<I(5,lb[5])<<I(5,lb[6])<<I(5,lb[1])<<I(5,lb[2])<<I(5,lb[3])<<endl;
	// cout <<I(5,bin6[4])<<I(5,bin6[5])<<I(5,bin6[6])<<I(5,bin6[1])<<I(5,bin6[2])<<I(5,bin6[3])<<endl;
	// cout <<I(5,ub[4])<<I(5,ub[5])<<I(5,ub[6])<<I(5,ub[1])<<I(5,ub[2])<<I(5,ub[3])<<endl;
	// cout << endl;

	int count_this_rt=0;
	for ( int i4 = lb[4]; i4 <= ub[4]; ++i4 ) {
		if     ( i4 == (int)hasher_.dimsizes()[4] ) bin6[4] = 0;
		else if ( i4 <  0                    ) bin6[4] = hasher_.dimsizes()[4]-1;
		else                                 bin6[4] = i4;
		for ( int i5 = lb[5]; i5 <= ub[5]; ++i5 ) {
			if     ( i5 == (int)hasher_.dimsizes()[5] ) bin6[5] = 0;
			else if ( i5 <  0                    ) bin6[5] = hasher_.dimsizes()[5]-1;
			else                                 bin6[5] = i5;
			for ( bin6[6] = lb[6]; (int)bin6[6] <= ub[6]; ++bin6[6] ) {
				for ( bin6[1] = lb[1]; (int)bin6[1] <= ub[1]; ++bin6[1] ) {
					for ( bin6[2] = lb[2]; (int)bin6[2] <= ub[2]; ++bin6[2] ) {
						for ( bin6[3] = lb[3]; (int)bin6[3] <= ub[3]; ++bin6[3] ) {

							// cout <<I(5,bin6[4])<<I(5,bin6[5])<<I(5,bin6[6])<<I(5,bin6[1])<<I(5,bin6[2])<<I(5,bin6[3])<<endl;

							Key key = hasher_.bin_index(bin6);
							MotifHash::MotifMap::const_iterator i = motif_umap_.equal_range(key).first;
							MotifHash::MotifMap::const_iterator e = motif_umap_.equal_range(key).second;
							for ( ; i != e; ++i ) {
								ResPairMotif const & sm(i->second);
								Real6 const ot6 = sm.rt();
								Xform const x2(ot6);
								Real const dis2 = x1.t.distance_squared(x2.t);
								Mat const r = x1.R * x2.R.transposed();
								Real const cos_theta = (r.xx()+r.yy()+r.zz()-1.0)/2.0;
								// Real sin2_theta = 1.0-cos_theta*cos_theta;
								Real const theta = degrees(acos(cos_theta));
								Real const angdis2 = ang2dis*ang2dis*theta*theta;
								if ( dis2+angdis2 <= radius*radius ) {
									results.push_back(sm);
									if ( ++count_this_rt>=option[mh::dump::limit_per_pair]() ) return;
								}
							}

						}
					}
				}
			}
		}
	}
	// cout << endl;
	// if(bin6[4]==0) utility_exit_with_message("arst");
	// if(bin6[5]==0) utility_exit_with_message("arst");
}

ostream & operator << (ostream & out, MotifHash const & x){
	numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> > const & bb(x.hasher_.bounding_box() );
	Real3 const & eo( x.hasher_.euler_offsets() );
	Real6 const & bw( x.hasher_.bin_widths() );
	out << bb.lower() << ' ' << bb.upper() << endl;
	out << eo[1] << ' ' << eo[2] << ' ' << eo[3] << endl;
	out << bw[1] << ' ' << bw[2] << ' ' << bw[3] << ' ' << bw[4] << ' ' << bw[5] << ' ' << bw[6] << endl;
	out << x.motif_umap_.size() << endl;
	out << x.motif_umap_.begin()->second << endl;
	// for(vector1<MotifHash::Motif>::const_iterator i = x.motifs_.begin(); i != x.motifs_.end(); ++i){
	//  out << *i << endl;
	// }
	return out;
}

void ResPairMotifQuery::init( Pose const & pose1_in, Pose const & pose2_in){
	core::pose::PoseOP pose1op( new Pose(pose1_in) );
	bool dssp1=true; BOOST_FOREACH ( char ss,pose1_in.secstruct() ) dssp1 &= ss=='L'; // dssp iff all L
	if ( dssp1 ) core::scoring::dssp::Dssp(*pose1op).insert_ss_into_pose(*pose1op,false);
	pose1_ = pose1op;
	pose2_ = NULL;
	if ( &pose1_in != &pose2_in ) {
		core::pose::PoseOP pose2op( new Pose(pose2_in) );
		bool dssp2=true; BOOST_FOREACH ( char ss,pose2_in.secstruct() ) dssp2 &= ss=='L'; // dssp iff all L
		if ( dssp2 ) core::scoring::dssp::Dssp(*pose2op).insert_ss_into_pose(*pose2op,false);
		pose2_ = pose2op;
	}

	inter_clash1_ = inter_clash2_ = NULL;
	auto_clash_ = false;
	match_ss1_ = option[mh::match::ss1].user() ? option[mh::match::ss1]() : option[mh::match::ss]();
	match_ss2_ = option[mh::match::ss2].user() ? option[mh::match::ss2]() : option[mh::match::ss]();
	match_aa1_ = option[mh::match::aa1].user() ? option[mh::match::aa1]() : option[mh::match::aa]();
	match_aa2_ = option[mh::match::aa2].user() ? option[mh::match::aa2]() : option[mh::match::aa]();
	match_chi_rmsd_ = 9e9;
	match_radius_ = option[mh::motif_match_radius]();
	interface_only_ = option[mh::match::interface_only]();
	overlap1_ = MATCH_ANY;
	overlap2_ = MATCH_ANY;
	max_ca_dis2_ = 2.0*MOTIF_HASH_CART_SIZE*MOTIF_HASH_CART_SIZE;
	if ( useres1_.size()==0 ) useres1_.resize(pose1_->n_residue(),true);
	if ( !pose2_ ) {
		pose2_ = pose1_;
		useres2_.resize(pose1_->n_residue(),false);
		if ( core::pose::symmetry::is_symmetric(*pose1_) ) {
			using namespace core::conformation::symmetry;
			using namespace core::pose::symmetry;
			SymmetryInfoCOP syminfo = symmetry_info(*pose1_);
			for ( Size ir = 1; ir <= pose1_->n_residue(); ++ir ) {
				useres1_[ir] = true;
				useres2_[ir] = syminfo->chi_is_independent(ir);
			}
		} else { // no symmetry
			useres2_ = useres1_;
		}
	} else {
		runtime_assert( !core::pose::symmetry::is_symmetric(*pose1_) );
	}
	if ( useres2_.size()==0 ) useres2_.resize(pose2_->n_residue(),true);
	useres1_ph_.resize(pose1_->n_residue(),true);
	useres2_ph_.resize(pose2_->n_residue(),true);
	useres1_po_.resize(pose1_->n_residue(),true);
	useres2_po_.resize(pose2_->n_residue(),true);

	runtime_assert(useres1_.size()==pose1_->n_residue());
	runtime_assert(useres2_.size()==pose2_->n_residue());
}

std::ostream & operator<<(std::ostream & out, ResPairMotifQuery const & opt){
	out << "ResPairMotifQuery(";
	out << " same_pose=" << (opt.pose1_==opt.pose2_);
	out << " auto_clash=" << opt.auto_clash_;
	out << " iface_only=" << opt.interface_only_;
	out << " radius=" << opt.match_radius_;
	out << ")";
	// out << "    useres1 " << opt.useres1_ << endl;
	// out << "    useres2 " << opt.useres2_ << endl;
	return out;
}

void ResPairMotifQuery::copy_opts_swapped(ResPairMotifQuery const & other){
	inter_clash1_ = other.inter_clash2_;
	inter_clash2_ = other.inter_clash1_;
	inter_clash1_xform_ = other.inter_clash2_xform_;
	inter_clash2_xform_ = other.inter_clash1_xform_;
	this->overlap1() = other.overlap2();
	this->overlap2() = other.overlap1();
	this->match_aa1() = other.match_aa2();
	this->match_aa2() = other.match_aa1();
	this->match_ss1() = other.match_ss2();
	this->match_ss2() = other.match_ss1();
	this->useres1() = other.useres2();
	this->useres2() = other.useres1();
	this->useres1_ph() = other.useres2_ph();
	this->useres2_ph() = other.useres1_ph();
	this->useres1_po() = other.useres2_po();
	this->useres2_po() = other.useres1_po();
	this->max_ca_dis2() = other.max_ca_dis2();
	this->match_chi_rmsd() = other.match_chi_rmsd();
	this->match_radius() = other.match_radius();
	this->interface_only() = other.interface_only();
	this->auto_clash() = other.auto_clash();
	this->clash_check() = other.clash_check();
}

MotifHits MotifHash::get_matching_motifs( ResPairMotifQuery const & opt ) const {
	MotifHits h;
	this->get_matching_motifs(opt,h);
	return h;
}

int MotifHash::get_matching_motifs( ResPairMotifQuery const & opt, MotifHits & hits ) const{
	MotifHits dummy;
	return get_matching_motifs(opt,hits,dummy);
}

int MotifHash::get_matching_motifs(ResPairMotifQuery const & opt, MotifHits & hits, MotifHits & newhits ) const {
	// std::cout << "PING /Users/sheffler/rosetta/scheme/source/src/core/scoring/motif/motif_hash_stuff.cc/motif_hash_stuff.cc:1678" << std::endl;  runtime_assert(motif_umap_.size()>0);
	runtime_assert(type_!=RPM_Type_NONE);
	using namespace core::pose::symmetry;
	TR.Info << opt << endl;
	newhits.clear();
	PoseCOP pose1cop(opt.pose1_);
	PoseCOP pose2cop(opt.pose2_);
	Pose const & pose1(*pose1cop);
	Pose const & pose2(*pose2cop);
	bool samepose = (pose1cop==pose2cop);
	Size res_begin1 = 1, res_end1 = pose1.total_residue();
	Size res_begin2 = 1, res_end2 = pose2.total_residue();
	Bools const &useres1(opt.useres1_), &useres2(opt.useres2_);
	// cout << opt << endl;

	bool const symmetry = core::pose::symmetry::is_symmetric(pose1);
	bool const multicomponent = symmetry && core::pose::symmetry::is_multicomponent(pose1);
	core::conformation::symmetry::SymmetryInfoCOP sym_info;
	Sizes intra_subs1, intra_subs2;
	if ( symmetry ) {
		runtime_assert( samepose );
		sym_info = symmetry_info(pose1);
		res_end1 = sym_info->num_total_residues_without_pseudo();
		res_end2 = sym_info->num_independent_residues();
		if ( multicomponent ) {
			if ( sym_dof_names(pose1).size() != 2 ) {
				cout << "WARNING: more than two symdofname, you better know what you're doing..." << endl;
			}
			// BOOST_FOREACH(std::string const & s,sym_dof_names(pose1)) cout << s << endl;
			intra_subs1 = get_jump_name_to_subunits(pose1,sym_dof_names(pose1)[1]);
			intra_subs2 = get_jump_name_to_subunits(pose1,sym_dof_names(pose1)[2]);
			// BOOST_FOREACH(Size i,intra_subs1) cout << "intra_subs1 " << i << endl;
			// BOOST_FOREACH(Size i,intra_subs2) cout << "intra_subs2 " << i << endl;
			// cout << opt.interface_only() << endl;
		}
	}

	xyzStripeHashPoseCOP ccheck1bb32=0,ccheck2bb32=0;
	if ( opt.auto_clash_ ) {
		ccheck1bb32 = xyzStripeHashPoseCOP( new xyzStripeHashPose(pose1,PoseCoordPickMode_BB,2.8) );
		if ( !samepose ) ccheck2bb32 = xyzStripeHashPoseCOP( new xyzStripeHashPose(pose2,PoseCoordPickMode_BB,2.8) );
		else ccheck2bb32 = ccheck1bb32;
	}

	for ( Size ir = res_begin1; ir <= res_end1; ++ir ) {

		if ( !useres1[ir] ) { /*<<"not useres1"<<endl;*/ continue; }
		if ( hits.num_hits1(ir)>0 && opt.overlap1()==NO_OVERLAP ) { /*<<"overlap1"<<endl;*/ continue; }
		char aa1 = pose1.residue(ir).name1(), ss1=pose1.secstruct(ir);
		// cout << ir << " " << aa1 << " " << ss1 << endl;
		for ( Size jr = res_begin2; jr <= res_end2; ++jr ) {
			runtime_assert( !symmetry || sym_info->subunit_index(jr)==1 );
			if ( !useres2[jr] ) { /*<<"not useres2"<<endl;*/ continue; };
			if ( hits.num_hits2(jr)>0 && opt.overlap2()==NO_OVERLAP ) { /*<<"overlap2"<<endl;*/ continue; }

			bool is_interface = true;
			if ( multicomponent && get_component_of_residue(pose1,ir) == get_component_of_residue(pose1,jr) ) {
				char comp = get_component_of_residue(pose1,jr);
				Sizes const & intra_subs(comp=='A'?intra_subs1:intra_subs2);
				if ( find(intra_subs.begin(), intra_subs.end(), sym_info->subunit_index(ir)) != intra_subs.end() ) { /*cout<<"intra sub"<<endl;*/ is_interface=false; };
			} else if (  symmetry && !multicomponent && sym_info->subunit_index(ir)==sym_info->subunit_index(jr) && pose1.chain(ir)==pose1.chain(jr) ) { /*<<"subunit index"<<endl;*/ is_interface=false; }
			else if ( !symmetry && samepose        && pose1.chain(ir)            ==pose1.chain(jr)             ) { /*<<"no sym, same chain"<<endl;*/ is_interface=false; };

			if ( !is_interface && opt.interface_only() ) continue;

			char aa2 = pose2.residue(jr).name1(), ss2=pose2.secstruct(jr);
			float dsqCA = pose1.xyz(AtomID(2,ir)).distance_squared(pose2.xyz(AtomID(2,jr)));
			if ( dsqCA > opt.max_ca_dis2() ) { /*<<"too far"<<endl;*/ continue; };

			// cout << ir << " "
			//      << jr << " "
			//      << get_component_of_residue(pose1,ir) << " "
			//      << get_component_of_residue(pose1,jr) << " "
			//      << sym_info->subunit_index(ir) << " "
			//        << sym_info->subunit_index(jr) << " "
			//        << endl;

			ResPairMotifs motifs; {
				Real6 const rt = get_residue_pair_rt6(pose1,ir,pose2,jr,type_);
				if ( rt[1]!=rt[1] ) { /*cout << "bad rt6" << endl;*/ continue; }
				find_motifs_with_radius(rt,opt.match_radius(),motifs);
			}
			// cout << "MOTIFS " << motifs.size() << endl;
			for ( ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i ) {

				///////////////////////////////////// hack intended to ignore strand pairs ////////////////////////////////////////
				if ( i->ss1()=='E' && i->ss2()=='E' && dsqCA < 10025.0 ) continue;
				///////////////////////////////////////////////////////////////////////////////


				MotifHit h(pose1cop,pose2cop,ir,jr,*i);
				if (  opt.match_ss1() && i->type1()!=RM_SC && i->type1()!=RM_PH && i->type1()!=RM_PO && ss1!=i->ss1() ) { /* << "ss mismatch" << endl;*/ continue; };
				if (  opt.match_ss1() && i->type2()!=RM_SC && i->type2()!=RM_PH && i->type2()!=RM_PO && ss2!=i->ss2() ) { /* << "ss mismatch" << endl;*/ continue; };
				if ( (opt.match_aa1() || i->type1()==RM_SC) && i->aa1() != aa1 && i->type1()!=RM_PH && i->type1()!=RM_PO ) { /* << "aa mismatch"<<endl;*/ continue; } // SC motif require aa match
				if ( (opt.match_aa2() || i->type2()==RM_SC) && i->aa2() != aa2 && i->type2()!=RM_PH && i->type2()!=RM_PO ) { /* << "aa mismatch"<<endl;*/ continue; }// SC motif require aa match

				i->fill_pose_with_motif(h.mpose());
				h.rms = align_motif_pose_super(h.mpose(),pose1,ir,pose2,jr,i->type());
				if ( h.rms > option[mh::dump::max_rms]() ) continue;
				// compute chi rmsd if same residues, weight chi1 highest, etc...
				if ( i->aa1() == aa1 && i->aa2() == aa2 ) {
					Real chirms = 0.0, totchi = 0.0;
					if ( i->type1()==RM_BB ) {
						for ( Size ichi = 1; ichi <= pose1.residue(ir).nchi(); ++ichi ) {
							chirms += sqr(angle_distance( pose1.residue(ir).chi(ichi), i->chi1(ichi) )) / pose1.residue(ir).nchi();
							totchi += 1.0/pose1.residue(ir).nchi();
						}
					}
					if ( i->type2()==RM_BB ) {
						for ( Size jchi = 1; jchi <= pose2.residue(jr).nchi(); ++jchi ) {
							chirms += sqr(angle_distance( pose2.residue(jr).chi(jchi), i->chi2(jchi) )) / pose2.residue(jr).nchi();
							totchi += 1.0/pose2.residue(jr).nchi();
						}
					}
					if ( totchi == 0.0 ) chirms = 0.0;
					else              chirms = sqrt(chirms/totchi);
					if ( chirms > opt.match_chi_rmsd() ) continue;

					h.chi_rmsd_ = chirms;
					// chirms = 0.0, totchi = 0.0;
					// if(i->type1()==RM_BB){
					//  for(Size ichi = 1; ichi <= pose1.residue(ir).nchi(); ++ichi ){
					//   chirms += sqr(angle_distance( pose1.residue(ir).chi(ichi), i->chi1(ichi) ));
					//   totchi += 1.0;
					//   cout << ir << " " << jr << " " << ichi << " " << pose1.residue(ir).chi(ichi) << " " << i->chi1(ichi) << " " << angle_distance( pose1.residue(ir).chi(ichi), i->chi1(ichi) ) << endl;
					//  }
					// }
					// if(i->type2()==RM_BB){
					//  for(Size jchi = 1; jchi <= pose2.residue(jr).nchi(); ++jchi ){
					//   chirms += sqr(angle_distance( pose2.residue(jr).chi(jchi), i->chi2(jchi) ));
					//   totchi += 1.0;
					//   cout << ir << " " << jr << " " << jchi << " " << pose2.residue(jr).chi(jchi) << " " << i->chi2(jchi) << " " << angle_distance( pose2.residue(jr).chi(jchi), i->chi2(jchi) ) << endl;
					//  }
					// }
					// cout << chirms << " " << totchi << " " << sqrt(chirms/totchi) << " " << opt.match_chi_rmsd() << endl;
				}

				bool motif_clashes = false;
				if ( opt.clash_check() ) {

					Size nst1er=6,nst2er=6;
					if ( ccheck1bb32 && ccheck2bb32 ) {
						if ( !samepose ) {
							for ( Size ia = nst1er; ia <= h.mpose().residue(1).nheavyatoms(); ++ia ) { if ( ccheck1bb32->clash_not_resid(h.mpose().xyz(AtomID(ia,1)),ir) ) { motif_clashes = true; }} if ( motif_clashes ) continue;
							for ( Size ia = nst1er; ia <= h.mpose().residue(1).nheavyatoms(); ++ia ) { if ( ccheck2bb32->clash          (h.mpose().xyz(AtomID(ia,1))   ) ) { motif_clashes = true; }} if ( motif_clashes ) continue;
							for ( Size ia = nst2er; ia <= h.mpose().residue(2).nheavyatoms(); ++ia ) { if ( ccheck1bb32->clash          (h.mpose().xyz(AtomID(ia,2))   ) ) { motif_clashes = true; }} if ( motif_clashes ) continue;
							for ( Size ia = nst2er; ia <= h.mpose().residue(2).nheavyatoms(); ++ia ) { if ( ccheck2bb32->clash_not_resid(h.mpose().xyz(AtomID(ia,2)),jr) ) { motif_clashes = true; }} if ( motif_clashes ) continue;
						} else {
							for ( Size ia = nst1er; ia <= h.mpose().residue(1).nheavyatoms(); ++ia ) { if ( ccheck1bb32->clash_not_resid(h.mpose().xyz(AtomID(ia,1)),ir,jr) ) { motif_clashes = true; }} if ( motif_clashes ) continue;
							for ( Size ia = nst2er; ia <= h.mpose().residue(2).nheavyatoms(); ++ia ) { if ( ccheck1bb32->clash_not_resid(h.mpose().xyz(AtomID(ia,2)),ir,jr) ) { motif_clashes = true; }} if ( motif_clashes ) continue;
						}
					}
					// cout << "MOTIF CLASH CHECK " << is_interface << " " << i->type() << " " << opt.inter_clash2() << " " << opt.inter_clash1() << endl;

					// THIS IS BROKEN SOMEHOW...
					// if( is_interface && i->type1()==RM_BB && opt.inter_clash2() != 0 ){
					//  // cout << "CLASH_CHECK " << *i << endl;
					//  // cout << "clash check inter2 " << endl;
					//  for(Size ia = 1; ia <= h.mpose().residue(1).nheavyatoms(); ++ia){
					//   if(opt.inter_clash2()->clash( opt.inter_clash2_xform() * h.mpose().xyz(AtomID(ia,1)) )) motif_clashes = true;
					//   if(motif_clashes) break;
					//  }
					//  // if(motif_clashes) continue;
					// }
					// if( is_interface && i->type2()==RM_BB && opt.inter_clash1() != 0 ){
					//  // cout << "CLASH_CHECK " << *i << endl;
					//  // cout << "clash check inter1 " << endl;
					//  for(Size ia = 1; ia <= h.mpose().residue(2).nheavyatoms(); ++ia){
					//   if(opt.inter_clash1()->clash(  opt.inter_clash1_xform() * h.mpose().xyz(AtomID(ia,2)) )) motif_clashes = true;
					//   if(motif_clashes) break;
					//  }
					//  // if(motif_clashes) continue;
					// }
				}
				if ( motif_clashes ) {
					// cout << "fail" << endl;
					// ++nbadmotifs; // FIX ABOVE CONTINUES
					continue;
				}
				// cout << h << endl;
				// h.mpose().dump_pdb("test1.pdb");
				// pose1.dump_pdb("pose1.pdb");
				// pose2.dump_pdb("pose2.pdb");
				// utility_exit_with_message("foo");

				newhits.push_back(h);
			}
		}
	}

	// here do intra if one in pair is at interface

	// if(newhits.size()==0) TR << "no hits, maybe you want -mh:pack:interface_only false?" << endl;
	BOOST_FOREACH ( MotifHit const & h,newhits ) hits.push_back(h);
	return newhits.size();
}

void filter_motifs(
	ResPairMotifs const & motifs_in,
	ResPairMotifs       & motifs_out
){
	RPM_FilterStats s;
	for ( ResPairMotifs::const_iterator i = motifs_in.begin(); i != motifs_in.end(); ++i ) {
		if ( i->filter(&s) ) motifs_out.push_back(*i);
	}
	TR << s;
}
void filter_motifs(
	ResPairMotifs & motifs
){
	ResPairMotifs tmp;
	filter_motifs(motifs,tmp);
	motifs = tmp;
}

/********************************************** XformScore **********************************************/
// this is a dummy
XformScore::XformScore():
	hasher_(numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> >(
	Vec(-MOTIF_HASH_CART_SIZE), Vec( MOTIF_HASH_CART_SIZE)
	),utility::fixedsizearray1<Size,3>(0),get_bins( -1, -1 )), cart_resl_(-1),angle_resl_(-1)
{
	sanity_check();
}

XformScore::XformScore(Real cart_resl, Real ang_resl):
	hasher_(
	numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> >(
	Vec(-MOTIF_HASH_CART_SIZE),
	Vec( MOTIF_HASH_CART_SIZE)
	),
	utility::fixedsizearray1<Size,3>(0),
	get_bins( cart_resl, ang_resl )
	),
	cart_size_(MOTIF_HASH_CART_SIZE),
	cart_resl_(cart_resl),
	angle_resl_(ang_resl)
{
	sanity_check();
}
void XformScore::sanity_check() const {
	runtime_assert_msg( cart_resl_ > 0.001 && angle_resl_ > 0.001,"cart & angle resl must be > 0.001");
	runtime_assert_msg( cart_resl_ < 999.0 && angle_resl_ < 999.0,"cart & angle resl must be < 999.0");
	runtime_assert_msg( std::fabs( fmod(180.0/angle_resl_,1.0) ) < 0.000000001, "angle resolution should evenly divide 180.0");
	runtime_assert_msg( cart_size_ == MOTIF_HASH_CART_SIZE, "cart_size != MOTIF_HASH_CART_SIZE");
}
void XformScore::clear() {
	scores_.clear();
}
Size XformScore::num_bins() const { return scores_.size(); }
XformScore::Score XformScore::tot_score() const {
	Score tot=0;
	BOOST_FOREACH ( ScoreMap::value_type const & v, scores_ ) {
		tot += v.second;
	}
	return tot;
}
Size XformScore::num_possible_bins() const {
	numeric::geometry::hashing::Size6 d = hasher_.dimsizes();
	return d[1]*d[2]*d[3]*d[4]*d[5]*d[6];
}

void XformScore::prune_small_bins(XformScore::Score thresh){
	ScoreMap n;
	BOOST_FOREACH ( ScoreMap::value_type const & v, scores_ ) if ( v.second>=thresh ) n.insert(v);
	scores_ = n;
}
void XformScore::print_quantiles(ostream & out, int num) const {
	if ( scores_.size()==0 ) {
		out << "EMPTY"<<endl;
		return;
	}
	Score mx = -9e9;
	utility::vector1<Score> vals;
	for ( ScoreMap::const_iterator i=scores_.begin(); i != scores_.end(); ++i ) {
		vals.push_back(i->second);
		mx = max(i->second,mx);
	}
	int W = 5;
	std::sort(vals.begin(),vals.end());
	for ( int i = 1; i <= num; ++i ) {
		Real idx = min((Real)1.0,max((Real)0.0,(Real)(i-1)/(Real)(num-1))) * (Real)(vals.size()) + 1.0;
		Size idx_size = min(vals.size(),max((Size)1,(Size)idx));
		out <<" " << F(W,1,vals[idx_size]);
	}
}
XformScore::Key XformScore::bin_index(Real6 const & rt) const {
	Real6 rt6 = rt;
	rt6[4] = rt6[4]< 0.0   ?  rt6[4]+360.0 : rt6[4];
	rt6[5] = rt6[5]< 0.0   ?  rt6[5]+360.0 : rt6[5];
	rt6[6] = rt6[6]< 0.0   ? -rt6[6]       : rt6[6];
	rt6[4] = rt6[4]>=360.0 ?  rt6[4]-360.0 : rt6[4];
	rt6[5] = rt6[5]>=360.0 ?  rt6[5]-360.0 : rt6[5];
	rt6[6] = rt6[6]>=360.0 ?  180.0-rt6[6] : rt6[6];
	return hasher_.bin_index(rt6);
}
void XformScore::set_score(Key const & key, Score const & count){
	ScoreMap::iterator it = scores_.find(key);
	if ( it==scores_.end() ) scores_.insert( std::make_pair( key, count  ) );
	else it->second = count;
}
void XformScore::add_score(Key const & key, Score const & count){
	ScoreMap::iterator it = scores_.find(key);
	if ( it==scores_.end() ) scores_.insert( std::make_pair( key, count  ) );
	else it->second += count;
}
void XformScore::max_score(Key const & key, Score const & count){
	ScoreMap::iterator it = scores_.find(key);
	if ( it==scores_.end() ) scores_.insert( std::make_pair( key, count  ) );
	else it->second = max(it->second,count);
}

void XformScore::aggregate_add(XformScore const & other){
	BOOST_FOREACH ( ScoreMap::value_type const & v,other.scores_ ) {
		add_score(v.first,v.second);
	}
}

void XformScore::aggregate_max(XformScore const & other){
	BOOST_FOREACH ( ScoreMap::value_type const & v,other.scores_ ) {
		max_score(v.first,v.second);
	}
}

XformScore::Score XformScore::score_of_bin(Real6 const & xform) const {
	Key key = bin_index(xform);
	ScoreMap::const_iterator i = scores_.find(key);
	if ( i != scores_.end() ) {
		return i->second;
	} else {
		return 0;
	}
}
XformScore::Score XformScore::score_of_bin(Xform const & xform) const {
	if ( xform.bad() ) return 0;
	return score_of_bin(xform.rt6());
}
XformScore::Score XformScore::score_of_bin(Key const & key) const {
	ScoreMap::const_iterator i = scores_.find(key);
	if ( i != scores_.end() ) {
		return i->second;
	} else {
		return 0;
	}
}
Real XformScore::score_of_bin_normalized(Real6 const & xform) const {
	Key key = bin_index(xform);
	Real6 center = hasher_.bin_center_point(hasher_.bin_from_index(key));
	Real count = score_of_bin(key);
	Real freq = count / sin(radians(center[6]));
	return freq;
}
bool XformScore::write_binary(ostream & out) const {
	numeric::Real  cart_size =  cart_size_;
	numeric::Real  cart_resl =  cart_resl_;
	numeric::Real angle_resl = angle_resl_;
	double   totsc   = tot_score();
	uint64_t ncounts = num_bins();
	uint64_t version = MOTIF_VERSION;
	out.write((char*)&   version,sizeof(uint64_t));
	out.write((char*)&   ncounts,sizeof(uint64_t));
	out.write((char*)&     totsc,sizeof(uint64_t));
	out.write((char*)& cart_size,sizeof(numeric::Real));
	out.write((char*)& cart_resl,sizeof(numeric::Real));
	out.write((char*)&angle_resl,sizeof(numeric::Real));
	for ( ScoreMap::const_iterator i = scores_.begin(); i != scores_.end(); ++i ) {
		Key key = i->first;
		Score count = i->second;
		out.write((char*)&key  ,sizeof(Key));
		out.write((char*)&count,sizeof(Score));
	}
	return true;
}
bool XformScore::write_binary(string const & fname) const {
	utility::io::ozstream out(fname,std::ios::out|std::ios::binary);
	if ( !out.good() ) {
		std::cerr << "XformScore::write_binary(fname): problem opening XformScore output file " << fname << endl;
		return false;
	}
	if ( !write_binary(out) ) {
		std::cerr << "XformScore::write_binary(fname): problem while writing XformScore file " << fname << endl;
		out.close();
		return false;
	}
	out.close();
	return true;
}
bool XformScore::read_binary(XformScoreOP & xs, std::istream & in, bool clearme, Real const & wt){
	numeric::Real  cart_size=0;
	numeric::Real  cart_resl=0;
	numeric::Real angle_resl=0;
	uint64_t ncounts,version;
	double totsc;
	in.read((char*)&   version,sizeof(uint64_t));
	if ( version != MOTIF_VERSION ) utility_exit_with_message("bad version .xh.bin.gz, should be "+string_of(MOTIF_VERSION));
	in.read((char*)&   ncounts,sizeof(uint64_t));
	in.read((char*)&     totsc,sizeof(uint64_t));
	in.read((char*)& cart_size,sizeof(numeric::Real));
	in.read((char*)& cart_resl,sizeof(numeric::Real));
	in.read((char*)&angle_resl,sizeof(numeric::Real));
	if ( !xs ) xs = XformScoreOP(new XformScore(cart_resl,angle_resl));
	else if ( clearme ) xs->clear();
	if ( cart_size != xs->cart_size_ || cart_resl!= xs->cart_resl_ || angle_resl != xs->angle_resl_ ) {
		cout << "from file " << cart_size  << " " << cart_resl  << " " << angle_resl  << endl;
		cout << "from this " << xs->cart_size_ << " " << xs->cart_resl_ << " " << xs->angle_resl_ << endl;
		utility_exit_with_message("hash param mismatch!");
	}
	for ( Size i = 1; i <= ncounts; ++i ) {
		Key key;
		Score score;
		in.read((char*)&key  ,sizeof(Key));
		in.read((char*)&score,sizeof(Score));
		xs->add_score(key,score*(Score)wt);
	}
	return true;
}
bool XformScore::read_binary(XformScoreOP & xs, string const & fname, bool clearme, Real const & wt){
	TR << "XformScore::read_binary(fname): " << fname << endl;
	if ( !utility::file::file_exists(fname) ) {
		TR.Error << "non-existent .xh.bin.gz file "+fname << endl;
		return false;
	}
	utility::io::izstream in(fname,std::ios::in|std::ios::binary);
	if ( !in.good() ) {
		std::cerr << "XformScore::read_binary(fname): problem opening XformScore input file " << fname << endl;
		return false;
	}
	if ( !read_binary(xs,in,clearme,wt) ) {
		std::cerr << "XformScore::read_binary(fname): problem while writing XformScore file " << fname << endl;
		in.close();
		return false;
	}
	in.close();
	return true;
}
bool XformScore::read_binary(XformScoreOP & xs, vector1<string> const & fnames, Real const & wt){
	if ( xs ) xs->clear();
	Size count = 0;
	TR << "Nfiles: " << fnames.size() << endl;
	for ( vector1<string>::const_iterator i = fnames.begin(); i != fnames.end(); ++i ) {
		if ( ++count %50 == 0 ) { TR << "... read " << count << endl; }
		if ( !read_binary(xs,*i,false,wt) ) {
			std::cerr << "XformScore::read_binary(fnames): error reading file "+*i << endl;
			if ( !option[basic::options::OptionKeys::mh::ignore_io_errors]() ) utility_exit_with_message("XformScore::read_binary(fnames): see tracer output");
		}
	}
	TR << endl;
	return true;
}

void XformScore::print_scores(ostream & out, std::string const tag) const {
	for ( ScoreMap::const_iterator i = scores_.begin(); i != scores_.end(); ++i ) {
		Key key = i->first;
		// Real6 center = hasher_.bin_center_point(hasher_.bin_from_index(key));
		Score count = i->second;
		out << tag << " " << key << " " << count << " " << hasher_.bin_center_point(hasher_.bin_from_index(key)) << endl;
	}
}

void XformScore::print_scores_norm(ostream & out, std::string const tag) const {
	for ( ScoreMap::const_iterator i = scores_.begin(); i != scores_.end(); ++i ) {
		Key key = i->first;
		Real6 center = hasher_.bin_center_point(hasher_.bin_from_index(key));
		Real count = i->second;
		Real freq = count / sin(radians(center[6]));
		out << tag << " " << key << " " << freq << endl;
	}
}

ostream & operator<<(ostream & out, XformScore const & xh){
	out << "XformScore "  << F(7,3,xh.cart_size_) << " " << F(7,3,xh.cart_resl_) << " " << F(7,3,xh.angle_resl_) << " " << F(10,3,(Real)xh.num_possible_bins()/1000000.0)<<"M";
	out << " " << F(9,3,xh.tot_score()/1000000.0) << "M " <<  F(9,3,(Real)xh.num_bins()/1000000.0) << "M " << F(5,1,xh.tot_score()/(Real)xh.num_bins()) << " ";
	xh.print_quantiles(out,20);
	return out;
}
ostream & operator<<(ostream & out, XformScoreSringMap const & x){
	BOOST_FOREACH ( XformScoreSringMap::value_type const & p,x ) {
		out << p.first << " " << p.second << endl;
	}
	return out;
}

/******************************************** MotifHashManager ******************************************/
MotifHashManager * MotifHashManager::instance_( 0 );

MotifHashManager * MotifHashManager::get_instance() {
	if ( instance_ == 0 ) {
		if ( option[basic::options::OptionKeys::mh::harvest::sep_aa1 ]() &&
				option[basic::options::OptionKeys::mh::harvest::sep_aa2 ]() ) {
			utility_exit_with_message("use -mh:harvest:sep_aa");
		}
#ifdef USE_OPENMP
		// omp_set_lock(&mman_lock);
		#pragma omp critical
		if ( instance_ == 0 ) {
			cout << "MotifHashManager initializing global instance..."; cout.flush();
			instance_ = new MotifHashManager();
			instance_->init();
			cout << "DONE!" << endl;
		}
		// omp_unset_lock(&mman_lock);
#else
		instance_ = new MotifHashManager();
		instance_->init();
#endif

	}
	return instance_;
}

void fill_xform_score_from_file(XformScoreOP & xs,string const & datfile, Real const & wt=1.0){
	// else utility_exit_with_message("already populated: "+ datfile);
	if ( !XformScore::read_binary(xs,datfile,false,wt) ) {
		TR << "WARNING: fail to read: " << datfile << endl;
		if ( utility::file::file_exists(datfile) ) {
			TR<<"file exists, error rading"<<endl;
		} else {
			TR<<"file does not exist " << datfile << endl;
		}
		xs = NULL;
	}
}

#ifdef USE_OPENMP
static pthread_t __preload_motifs_thread__;
#endif

void* preload_motif_data_pthread_wrapper(void* ptr)
{
	TR << "preload_motif_data_pthread_wrapper" << std::endl;
#ifdef USE_OPENMP
	#pragma omp critical
	cout << "WARNING: pre-loading big-ass motif file(s) in a separate thread" << endl;
#endif
	MotifHashManager & mman( *((MotifHashManager*)ptr) );

	preload_motif_data(mman);

	return NULL;
}

void preload_motif_data(MotifHashManager & mman){
	TR << "preload_motif_data" << std::endl;
	ResPairMotifsStringMap tmp_map;
	if ( option[basic::options::OptionKeys::mh::path::motifs_BB_PH]().size() > 0 ) {
		ResPairMotifs motifs;
		load_motifs( option[basic::options::OptionKeys::mh::path::motifs_BB_PH](), motifs, &tmp_map );
		MotifHashManager* mman = MotifHashManager::get_instance();
		mman->motif_hash_BB_PH_ =  MotifHashOP(new MotifHash(motifs));
		// /*kill me!*/mman->motif_hash_BB_PH_->hasher().tree_init(5);
		BOOST_FOREACH ( string const & fn,option[basic::options::OptionKeys::mh::path::motifs_BB_PH]() ) cout << "  ...loaded " << fn << endl;
	}
	if ( option[basic::options::OptionKeys::mh::path::motifs_BB_PO]().size() > 0 ) {
		ResPairMotifs motifs;
		load_motifs( option[basic::options::OptionKeys::mh::path::motifs_BB_PO](), motifs, &tmp_map );
		MotifHashManager* mman = MotifHashManager::get_instance();
		mman->motif_hash_BB_PO_ = MotifHashOP(new MotifHash(motifs));
		// /*kill me!*/mman->motif_hash_BB_PO_->hasher().tree_init(5);
		BOOST_FOREACH ( string const & fn,option[basic::options::OptionKeys::mh::path::motifs_BB_PO]() ) cout << "  ...loaded " << fn << endl;
	}
	if ( option[basic::options::OptionKeys::mh::path::motifs_BB_BB]().size() > 0 ) {
		ResPairMotifs motifs;
		load_motifs( option[basic::options::OptionKeys::mh::path::motifs_BB_BB](), motifs, &tmp_map );
		MotifHashManager* mman = MotifHashManager::get_instance();
		mman->motif_hash_BB_BB_ = MotifHashOP(new MotifHash(motifs));
		// /*kill me!*/mman->motif_hash_BB_BB_->hasher().tree_init(5);
		BOOST_FOREACH ( string const & fn,option[basic::options::OptionKeys::mh::path::motifs_BB_BB]() ) cout << "  ...loaded " << fn << endl;
	}
	if ( option[basic::options::OptionKeys::mh::path::motifs]().size() > 0 ) {
		ResPairMotifs motifs;
		load_motifs( option[basic::options::OptionKeys::mh::path::motifs](), motifs, &tmp_map );
		// MotifHashManager* mman = MotifHashManager::get_instance();
		cout << "WARNING: using -motifs is depricated" << endl;
		BOOST_FOREACH ( string const & fn,option[basic::options::OptionKeys::mh::path::motifs]() ) cout << "  ...loaded " << fn << endl;
	}
	if ( option[basic::options::OptionKeys::mh::path::motifs_SC_BB]().size() > 0 ) {
		ResPairMotifs motifs;
		load_motifs( option[basic::options::OptionKeys::mh::path::motifs_SC_BB](), motifs, &tmp_map );
		MotifHashManager* mman = MotifHashManager::get_instance();
		mman->motif_hash_SC_BB_ = MotifHashOP(new MotifHash(motifs));
		// /*kill me!*/mman->motif_hash_SC_BB_->hasher().tree_init(5);
		BOOST_FOREACH ( string const & fn,option[basic::options::OptionKeys::mh::path::motifs_SC_BB]() ) cout << "  ...loaded " << fn << endl;
	}
	BOOST_FOREACH ( ResPairMotifsStringMap::value_type const & v, tmp_map ) {
		string tag = utility::file_basename(v.first);
		if ( tag.size()>=11 && tag.substr(tag.size()-11)==".rpm.bin.gz" ) tag = tag.substr(0,tag.size()-11);
		TR << "adding to motifs_by_fname " << tag << endl;
		if ( mman.motifs_by_fname_.find(tag)!=mman.motifs_by_fname_.end() ) utility_exit_with_message("duplicate motif tag: "+tag);
		mman.motifs_by_fname_[tag] = utility::pointer::shared_ptr<class core::scoring::motif::MotifHash>( new MotifHash(v.second) );
		// /*kill me!*/mman.motifs_by_fname_[tag]->hasher().tree_init(5);
	}
	mman.done_loading_ = true;
}

MotifHashManager::MotifHashManager() :
	motif_hash_SC_BB_(NULL),
	motif_hash_BB_BB_(NULL),
	motif_hash_BB_PH_(NULL),
	motif_hash_BB_PO_(NULL),
	xform_score_BB_PH_(NULL),
	xform_score_BB_PO_(NULL),
	xform_score_PH_PO_(NULL),
	xform_score_frags_(NULL),
	done_loading_(false)
{}
void MotifHashManager::init(){
	TR << "init MotifHashManager" << endl;
	if ( option[mh::path::scores_BB_BB].user() ) {
		TR << "reading xform_score_data_BB_BB" << endl;
		key_mask_BB_BB_ = 0;
		if ( option[mh::score::use_ss1]() ) key_mask_BB_BB_ ^= 255u;
		if ( option[mh::score::use_ss2]() ) key_mask_BB_BB_ ^= 255u<<8;
		if ( option[mh::score::use_aa1]() ) key_mask_BB_BB_ ^= 255u<<16;
		if ( option[mh::score::use_aa2]() ) key_mask_BB_BB_ ^= 255u<<24;
		// cout << std::bitset<32>(key_mask_BB_BB_) << endl;
		const char * AA1 = option[mh::score::use_aa1]()? "ACDEFGHIKLMNPQRSTVWY" : " ";
		const char * AA2 = option[mh::score::use_aa2]()? "ACDEFGHIKLMNPQRSTVWY" : " ";
		const char * SS1 = option[mh::score::use_ss1]()? "EHL" : " ";
		const char * SS2 = option[mh::score::use_ss2]()? "EHL" : " ";
		bool allnull = true;
		BOOST_FOREACH ( char const & ss1,SS1 ) {
			BOOST_FOREACH ( char const & ss2,SS2 ) {
				BOOST_FOREACH ( char const & aa1,AA1 ) {
					BOOST_FOREACH ( char const & aa2,AA2 ) {
						Key k = key_mask_BB_BB_ & ((Key)ss1 | ((Key)ss2<<8) | ((Key)aa1<<16) | ((Key)aa2<<24));
						xform_scores_BB_BB_[k] = NULL;
						// cout << ss1 << " " << ss2 << " " << aa1 << " " << aa2 << " " << k << endl;
						string ss1s = ss1!=' ' ? "_"+string_of(ss1) : "";
						string ss2s = ss2!=' ' ?     string_of(ss2) : "";
						string aa1s = aa1!=' ' ?     string_of(aa1) : "";
						string aa2s = aa2!=' ' ?     string_of(aa2) : "";
						BOOST_FOREACH ( string const & fname_base, option[mh::path::scores_BB_BB]() ) {
							string fname = fname_base + aa1s+aa2s+ss1s+ss2s+".xh.bin.gz";
							Real wt = 1.0;
							// if( (ss1=='H' && ss2=='E') || (ss1=='E' && ss2=='H') ){ wt =  5.5; /*cout << "FIXME SS wt" << endl;*/ }
							// if( (ss1=='E' && ss2=='E')                           ){ wt = 10.2; /*cout << "FIXME SS wt" << endl;*/ }
							fill_xform_score_from_file(xform_scores_BB_BB_[k],fname,wt);
							allnull &= (NULL==xform_scores_BB_BB_[k]);
							// if(xform_scores_BB_BB_[k]) cout << "read " << fname << " " << k << endl;
							// else                 cout << "fail " << fname << endl;
						}
					}
				}
			}
		}
		if ( allnull ) {
			BOOST_FOREACH ( string s, option[mh::path::scores_BB_BB]() ) cout << "mh:path:scores_BB_BB " << s << endl;
			utility_exit_with_message("no score files for mh::path::scores_BB_BB");
		}
	}
	if ( option[mh::path::scores_frags].user() ) {
		TR << "reading xform_score_data_frags" << endl;
		BOOST_FOREACH ( string const & fname, option[mh::path::scores_frags]() ) {
			fill_xform_score_from_file(xform_score_frags_,fname,1.0);
		}
		// key_mask_frags_ = 0;
		// if(option[mh::score::use_ss1]()) key_mask_frags_ ^= 255u;
		// if(option[mh::score::use_ss2]()) key_mask_frags_ ^= 255u<<8;
		// const char * SS1 = option[mh::score::use_ss1]()? "EHL" : " ";
		// const char * SS2 = option[mh::score::use_ss2]()? "EHL" : " ";
		// bool allnull = true;
		// BOOST_FOREACH(char const & ss1,SS1){
		// BOOST_FOREACH(char const & ss2,SS2){
		//  Key k = key_mask_frags_ & ((Key)ss1 | ((Key)ss2<<8));
		//  xform_scores_frags_[k] = new XformScore;
		//  // cout << ss1 << " " << ss2 << " " << aa1 << " " << aa2 << " " << k << endl;
		//  string ss1s = ss1!=' ' ? "_"+string_of(ss1) : "";
		//  string ss2s = ss2!=' ' ?     string_of(ss2) : "";
		//  BOOST_FOREACH( string const & fname_base, option[mh::path::scores_frags]() ){
		//   string fname = fname_base + ss1s+ss2s+".xh.bin.gz";
		//   fill_xform_score_from_file(xform_scores_frags_[k],fname,1.0);
		//   // if(xform_scores_frags_[k]) cout << "read " << fname << " " << k << endl;
		//   // else                 cout << "fail " << fname << endl;
		//   allnull &= (NULL==xform_scores_frags_[k]);
		//  }
		// }}
		// if(allnull) utility_exit_with_message("no score files for option[mh::path::scores_frags");
	}
	if ( option[mh::path::scores_SC_SC].user() ) {
		TR << "reading xform_score_data_SC_SC" << endl;
		const char * AA = "ACDEFGHIKLMNPQRSTVWY";
		bool allnull = true;
		BOOST_FOREACH ( char const & aa1,AA ) {
			BOOST_FOREACH ( char const & aa2,AA ) {
				Key const k = (((Key)aa1<<0) | ((Key)aa2<<8));
				xform_scores_SC_SC_[k] = NULL;
				BOOST_FOREACH ( string const & fname,option[mh::path::scores_SC_SC]() ) {
					fill_xform_score_from_file(xform_scores_SC_SC_[k],fname+"_"+aa1+aa2+".xh.bin.gz",1.0);
					allnull &= (NULL==xform_scores_SC_SC_[k]);
				}
			}
		}
		if ( allnull ) {
			BOOST_FOREACH ( string s, option[mh::path::scores_SC_SC]() ) cout << "mh:path:scores_SC_SC " << s << endl;
			utility_exit_with_message("no score files for mh::path::scores_SC_SC");
		}
	}
	if ( option[mh::path::scores_SC_BB].user() ) {
		TR << "reading xform_score_data_SC_BB" << endl;
		const char * AA = "ACDEFGHIKLMNPQRSTVWY";
		bool allnull = true;
		BOOST_FOREACH ( char const & aa,AA ) {
			xform_scores_SC_BB_[aa] = NULL;
			BOOST_FOREACH ( string const & fname,option[mh::path::scores_SC_BB]() ) {
				fill_xform_score_from_file(xform_scores_SC_BB_[aa],fname+"_"+aa+".xh.bin.gz",1.0);
				allnull &= (NULL==xform_scores_SC_BB_[aa]);
			}
		}
		if ( allnull ) {
			BOOST_FOREACH ( string s, option[mh::path::scores_SC_BB]() ) cout << "mh:path:scores_SC_BB " << s << endl;
			utility_exit_with_message("no score files for mh::path::scores_SC_BB");
		}
	}
	if ( option[mh::path::scores_BB_PH].user() ) {
		TR << "reading xform_score_data_BB_PH" << endl;
		xform_score_BB_PH_ = NULL;
		BOOST_FOREACH ( string const & fname,option[mh::path::scores_BB_PH]() ) {
			fill_xform_score_from_file(xform_score_BB_PH_,fname,1.0);
		}
	}
	if ( option[mh::path::scores_BB_PO].user() ) {
		TR << "reading xform_score_data_BB_PO" << endl;
		xform_score_BB_PO_ = NULL;
		BOOST_FOREACH ( string const & fname,option[mh::path::scores_BB_PO]() ) {
			fill_xform_score_from_file(xform_score_BB_PO_,fname,1.0);
		}
	}
	if ( option[mh::path::scores_PH_PO].user() ) {
		TR << "reading xform_score_data_PH_PO" << endl;
		xform_score_PH_PO_ = NULL;
		BOOST_FOREACH ( string const & fname,option[mh::path::scores_PH_PO]() ) {
			fill_xform_score_from_file(xform_score_PH_PO_,fname);
		}
	}
	BOOST_FOREACH ( std::string const & s,option[mh::path::motifs      ]() ) add_motif_set_name(utility::file_basename(s));
	BOOST_FOREACH ( std::string const & s,option[mh::path::motifs_BB_BB]() ) add_motif_set_name(utility::file_basename(s));
	BOOST_FOREACH ( std::string const & s,option[mh::path::motifs_SC_BB]() ) add_motif_set_name(utility::file_basename(s));
	BOOST_FOREACH ( std::string const & s,option[mh::path::motifs_BB_PH]() ) add_motif_set_name(utility::file_basename(s));
	BOOST_FOREACH ( std::string const & s,option[mh::path::motifs_BB_PO]() ) add_motif_set_name(utility::file_basename(s));
	// BOOST_FOREACH(std::string const & s,motif_set_names_){
	//  string key(utility::file_basename(s));
	//  if(key.size()>=11 && key.substr(key.size()-11)==".rpm.bin.gz") key = key.substr(0,key.size()-11);
	//  add_motif_set_name(key);
	// }
#ifdef USE_OPENMP
	//#pragma omp critical
	{
		static bool launched = false;
		if ( !launched ) {
			launched = true;
			pthread_create( &__preload_motifs_thread__, NULL, preload_motif_data_pthread_wrapper, this );
		}
		// pthread_join(__preload_motifs_thread__,NULL); // test...
	}
#else
	static bool launched = false;
	if ( !launched ) preload_motif_data(*this);
	launched = true;
#endif

}

void MotifHashManager::add_motif_set_name(std::string const & motifset){
	TR << "add_motif_set_name " << motifset << std::endl;
	string s = motifset;
	if ( s.size() > 11 && s.substr(s.size()-11) == ".rpm.bin.gz" ) {
		s = s.substr(0,s.size()-11);
	}
	motif_set_names_.insert(s);
}

MotifHashManager::~MotifHashManager(){
	/* transitioned from naked pointers to OP for Luki's code transition. This code should be automatically deleted
	BOOST_FOREACH(XformScoreMap::value_type & v,xform_scores_BB_BB_) delete v.second;
	BOOST_FOREACH(XformScoreMap::value_type & v,xform_scores_SC_BB_) delete v.second;
	if(motif_hash_SC_BB_) delete motif_hash_SC_BB_;
	if(motif_hash_BB_BB_) delete motif_hash_BB_BB_;
	if(motif_hash_BB_PH_) delete motif_hash_BB_PH_;
	if(motif_hash_BB_PO_) delete motif_hash_BB_PO_;
	if(xform_score_BB_PH_) delete xform_score_BB_PH_;
	if(xform_score_BB_PO_) delete xform_score_BB_PO_;
	if(xform_score_PH_PO_) delete xform_score_PH_PO_;
	*/
}

MotifHashCOP MotifHashManager::get_motif_hash_BB_BB(){
	if ( !done_loading_ ) {
		if ( option[basic::options::OptionKeys::mh::path::motifs_BB_BB]().size() > 0 ) {
#ifdef USE_OPENMP
			cout << "fancy threaded preload of BB_BB Motifs not done in time.... joining..." << endl;
			pthread_join(__preload_motifs_thread__,NULL);
#else
			utility_exit_with_message("BB_BB motif data should have been loaded at this point!");
#endif
		} else return NULL;
	}
	return motif_hash_BB_BB_;
}
MotifHashCOP MotifHashManager::get_motif_hash_SC_BB(){
	if ( !done_loading_ ) {
		if ( option[basic::options::OptionKeys::mh::path::motifs_SC_BB]().size() > 0 ) {
#ifdef USE_OPENMP
			cout << "fancy threaded preload of SC_BB Motifs not done in time.... joining..." << endl;
			pthread_join(__preload_motifs_thread__,NULL);
#else
			utility_exit_with_message("SC_BB motif data should have been loaded at this point!");
#endif
		} else return NULL;
	}
	return motif_hash_SC_BB_;
}
MotifHashCOP MotifHashManager::get_motif_hash_SC_SC(){
	if ( !done_loading_ ) {
		if ( option[basic::options::OptionKeys::mh::path::motifs_SC_SC]().size() > 0 ) {
#ifdef USE_OPENMP
			cout << "fancy threaded preload of SC_SC Motifs not done in time.... joining..." << endl;
			pthread_join(__preload_motifs_thread__,NULL);
#else
			utility_exit_with_message("SC_SC motif data should have been loaded at this point!");
#endif
		} else return NULL;
	}
	return motif_hash_SC_SC_;
}
MotifHashCOP MotifHashManager::get_motif_hash_BB_PH(){
	if ( !done_loading_ ) {
		if ( option[basic::options::OptionKeys::mh::path::motifs_BB_PH]().size() > 0 ) {
#ifdef USE_OPENMP
			cout << "fancy threaded preload of BB_PH Motifs not done in time.... joining..." << endl;
			pthread_join(__preload_motifs_thread__,NULL);
#else
			utility_exit_with_message("BB_PH motif data should have been loaded at this point!");
#endif
		} else return NULL;
	}
	return motif_hash_BB_PH_;
}
MotifHashCOP MotifHashManager::get_motif_hash_BB_PO(){
	if ( !done_loading_ ) {
		if ( option[basic::options::OptionKeys::mh::path::motifs_BB_PO]().size() > 0 ) {
#ifdef USE_OPENMP
			cout << "fancy threaded preload of BB_PO Motifs not done in time.... joining..." << endl;
			pthread_join(__preload_motifs_thread__,NULL);
#else
			utility_exit_with_message("BB_PO motif data should have been loaded at this point!");
#endif
		} else return NULL;
	}
	return motif_hash_BB_PO_;
}
MotifHashStringMap const & MotifHashManager::get_motif_hash_by_fname(){
	if ( !done_loading_ ) { // proxy...
		if ( option[basic::options::OptionKeys::mh::path::motifs      ]().size() > 0 ||
				option[basic::options::OptionKeys::mh::path::motifs_BB_BB]().size() > 0 ||
				option[basic::options::OptionKeys::mh::path::motifs_SC_BB]().size() > 0 ||
				option[basic::options::OptionKeys::mh::path::motifs_BB_PH]().size() > 0 ||
				option[basic::options::OptionKeys::mh::path::motifs_BB_PO]().size() > 0
				) {
#ifdef USE_OPENMP
			#pragma omp critical
			{
				cout << "fancy threaded preload of Motifs not done in time.... joining..." << endl;
				pthread_join(__preload_motifs_thread__,NULL);
			}
#else
			utility_exit_with_message("motif data should have been loaded at this point!");
#endif
		}
	}
	return motifs_by_fname_;
}
MotifHashCOP MotifHashManager::get_motif_hash_by_fname(std::string const & fname){
	string key(utility::file_basename(fname));
	if ( key.size()>=11 && key.substr(key.size()-11)==".rpm.bin.gz" ) key = key.substr(0,key.size()-11);
	get_motif_hash_by_fname(); // for sync
	if ( motifs_by_fname_.find(key) == motifs_by_fname_.end() ) {
		for ( MotifHashStringMap::const_iterator i = motifs_by_fname_.begin(); i!=motifs_by_fname_.end(); ++i ) {
			cout << "in motifs_by_fname_ " << i->first << endl;
		}
		utility_exit_with_message("unknown motif file name "+key);
	}
	return motifs_by_fname_[key];
}
bool MotifHashManager::have_motif_set_named(std::string const & fname) const {
	string key(utility::file_basename(fname));
	if ( key.size()>=11 && key.substr(key.size()-11)==".rpm.bin.gz" ) key = key.substr(0,key.size()-11);
	return motif_set_names_.find(key)!=motif_set_names_.end();
}

XformScoreCOP MotifHashManager::get_xform_score_PH_PO(){
	return xform_score_PH_PO_;
}

XformScoreCOP MotifHashManager::get_xform_score_BB_BB(char ss1, char ss2, char aa1, char aa2){
	Key const k = key_mask_BB_BB_ & ((Key)ss1 | ((Key)ss2<<8) | ((Key)aa1<<16) | ((Key)aa2<<24));
	// cout << "lookup " << k << " " << ss1 <<' '<< ss2 <<' '<< aa1 <<' '<< aa2 << endl;
	XformScoreMap::const_iterator i = xform_scores_BB_BB_.find(k);
	if ( i==xform_scores_BB_BB_.end() ) return NULL;
	return XformScoreCOP(i->second);
}
XformScoreCOP MotifHashManager::get_xform_score_frags(){
	// Key k = key_mask_frags_ & ((Key)ss1 | ((Key)ss2<<8));
	// XformScoreMap::const_iterator i = xform_scores_frags_.find(k);
	// if(i==xform_scores_frags_.end()) return NULL;
	// return i->second;
	return XformScoreCOP(xform_score_frags_);
}
XformScoreCOP MotifHashManager::get_xform_score_SC_BB(char aa1){
	Key const k = aa1;
	XformScoreMap::const_iterator i = xform_scores_SC_BB_.find(k);
	if ( i==xform_scores_SC_BB_.end() ) return NULL;
	return XformScoreCOP(i->second);
}
XformScoreCOP MotifHashManager::get_xform_score_SC_SC(char aa1, char aa2){
	Key const k = (((Key)aa1<<0) | ((Key)aa2<<8));
	XformScoreMap::const_iterator i = xform_scores_SC_SC_.find(k);
	if ( i==xform_scores_SC_SC_.end() ) return NULL;
	return XformScoreCOP(i->second);
}
XformScoreCOP MotifHashManager::get_xform_score_BB_PH(){
	return XformScoreCOP(xform_score_BB_PH_);
}
XformScoreCOP MotifHashManager::get_xform_score_BB_PO(){
	return xform_score_BB_PO_;
}

int MotifHashManager::get_matching_motifs( ResPairMotifQuery const & opt, MotifHits & hits ) const {
	if ( motif_set_names_.size()==0 ) return 0;
	ResPairMotifQuery * optswap = NULL;
	if ( opt.pose1()!=opt.pose2() ) {
		optswap = new ResPairMotifQuery(*opt.pose2(),*opt.pose1());
		optswap->copy_opts_swapped(opt);
	}
	MotifHits newhits;
	int count=0;
	BOOST_FOREACH ( std::string const & mset, motif_set_names_/*option[mh::path::motifs]()*/ ) {
		TR << "check file " << mset << endl;
		MotifHashCOP mh = MotifHashManager::get_instance()->get_motif_hash_by_fname(mset);
		int n = mh->get_matching_motifs(opt,hits,newhits);
		// cout << "MotifHashManager::get_matching_motifs " << hits.size() << " " << newhits.size() << endl;
		count += n;
		// AMW: if you want to re-enable the below tracer, reenable this conditional!
		// but as it was, the inside was actually doing nothing because n
		// would be reassigned after being increased
		//if( opt.pose1()!=opt.pose2() && mh->type1()!=mh->type2() ){
		// cout << "==================== swapping motif pose order ==================" << endl;
		//n += mh->get_matching_motifs(*optswap,hits,newhits);
		//}
		// TR << n << " new hits from mset: " << mset << endl;
		// newhits.dump_motifs_pdb(utility::file_basename(fname)+"_motifs"+string_of(++count)+".pdb.gz");
	}
	if ( optswap ) delete optswap;
	return count;
}

/********************************************** // packing **********************************************/

//  MotifRotamerSetOperationCOP
//  MotifHash::get_matching_motifs_rotsetop(
//  Pose const & pose,
//  core::pack::task::PackerTaskCOP ptask,
//  Real radius,
//  ClashCheckerPtr clash_check
//  ) const {
//  if(!ptask) ptask = core::pack::task::TaskFactory::create_packer_task(pose);
//  // cout << *ptask << endl;
//  MotifHits hits;
//  get_matching_motifs(pose,ptask,hits,clash_check,radius);
//  return new MotifRotamerSetOperation(pose,hits);
//  }

//  MotifRotamerSetOperation::MotifRotamerSetOperation(
//  Pose const & refpose,
//  MotifHits const & motifs
//  ):
//  motifs_(motifs)
//  {
//  BOOST_FOREACH(MotifHit const & h, motifs_){
//   Pose *pose = new Pose;
//   h.motif.fill_pose_with_motif(*pose);
//   align_motif_pose_super(*pose,refpose,h.residue1,refpose,h.residue2);
//   if(h.residue1>(int)res1_poses_.size()) res1_poses_.resize(h.residue1);
//   if(h.residue2>(int)res2_poses_.size()) res2_poses_.resize(h.residue2);
//   res1_poses_[h.residue1].push_back(pose);
//   res2_poses_[h.residue2].push_back(pose);
//  }
//  }

//  void
//  MotifRotamerSetOperation::alter_rotamer_set(
//  Pose const & /*pose*/,
//  core::scoring::ScoreFunction const & /*sfxn*/,
//  core::pack::task::PackerTask const & /*ptask*/,
//  core::graph::GraphCOP packer_neighbor_graph,
//  core::pack::rotamer_set::RotamerSet & rotamer_set
//  ){
//  BOOST_FOREACH(PoseCOP pp, res1_poses_[rotamer_set.resid()]){
//   rotamer_set.add_rotamer(pp->residue(1));
//  }
//  BOOST_FOREACH(PoseCOP pp, res2_poses_[rotamer_set.resid()]){
//   rotamer_set.add_rotamer(pp->residue(2));
//  }
//  }

//  void
//  MotifRotamerSetOperation::dump_pdb(string const & fname) const {
//  utility::io::ozstream out(fname);
//  Size ano = 0;
//  Size nres = res1_poses_.size();
//  for(Size ir = 1; ir <= nres; ++ir){
//   BOOST_FOREACH(PoseCOP pose, res1_poses_[ir]){
//    out << "MODEL" << endl;
//    core::io::pdb::dump_pdb_residue(pose->residue(1),ano,out);
//    core::io::pdb::dump_pdb_residue(pose->residue(2),ano,out);
//    out << "ENDMDL" << endl;
//   }
//   BOOST_FOREACH(PoseCOP pose, res2_poses_[ir]){
//    out << "MODEL" << endl;
//    core::io::pdb::dump_pdb_residue(pose->residue(1),ano,out);
//    core::io::pdb::dump_pdb_residue(pose->residue(2),ano,out);
//    out << "ENDMDL" << endl;
//   }
//  }
//  out.close();
//  }

//  Real
//  MotifRotamerSetOperation::increase_packer_residue_radius(
//  Pose const & /*pose*/,
//  core::pack::task::PackerTaskCOP /*the_task*/,
//  core::Size /*residue_in*/
//  ) const {
//  return 0.0;
//  }

}
}
}

