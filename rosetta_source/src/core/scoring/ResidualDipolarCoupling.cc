// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ResidualDipolarCoupling.cc
/// @brief  Uses NMR RDC for scoring
/// @author Oliver Lange
/// @author Srivatsan Raman
/// @author Nikolas Sgourakis


//Unit headers
#include <core/scoring/ResidualDipolarCoupling.hh>

// Package headers

// Project headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>

// AUTO-REMOVED #include <basic/options/after_opts.hh>
#include <basic/options/option.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/Tracer.hh>

//#include <protocols/loops/Loops.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

//C++ headers
#include <iostream>
#include <fstream>
#include <string>



/// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rdc.OptionKeys.gen.hh>

static basic::Tracer tr("core.scoring.ResidualDipolarCoupling");

namespace core {
namespace scoring {

inline Real sqr(Real x) {
	return x * x;
}
//////////////////////////////////////////////////////
//@brief reads in RDC data file
//////////////////////////////////////////////////////
extern void store_RDC_in_pose(ResidualDipolarCouplingOP rdc_info,
		core::pose::Pose& pose) {
	// ////using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA;
	pose.data().set(core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA, rdc_info);
}

extern ResidualDipolarCouplingCOP retrieve_RDC_from_pose(
		core::pose::Pose const& pose) {
	// ////using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA;
	if (pose.data().has(core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA)) {
		return static_cast<ResidualDipolarCoupling const *> (pose.data().get_const_ptr(
				core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA)());
	};
	return NULL;
}

extern ResidualDipolarCouplingOP retrieve_RDC_from_pose(core::pose::Pose& pose) {
	// ////using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA;
	if (pose.data().has(core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA)) {
		return static_cast<ResidualDipolarCoupling*> (pose.data().get_ptr(
				core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA)());
	};
	return NULL;
}

void ResidualDipolarCoupling::show(std::ostream& out) const {
	Size ct=0;
	for (RDC_lines::const_iterator it=All_RDC_lines_.begin();it!=All_RDC_lines_.end();it++){
				out << "RDC "<<++ct << "     ";
		out << (*it) << std::endl;
	}
}

void RDC::show(std::ostream& out) const {
	using namespace ObjexxFCL::fmt;
	out << RJ(4, res1_) << RJ(3, atom1_) << " " << RJ(4, res2_)
			<< RJ(3, atom2_) << " " << RJ(5, Jdipolar_);
}

std::ostream& operator<<(std::ostream& out, RDC const& rdc) {
	rdc.show(out);
	return out;
}

std::ostream& operator<<(std::ostream& out, ResidualDipolarCoupling const& rdc) {
	rdc.show(out);
	return out;
}

ResidualDipolarCoupling::ResidualDipolarCoupling(ResidualDipolarCoupling const& other) :
	basic::datacache::CacheableData(other) {
	All_RDC_lines_ = other.All_RDC_lines_;
	preprocess_data();
	reserve_buffers();
}


//explicit assignment operator to initialize buffers
ResidualDipolarCoupling&
ResidualDipolarCoupling::operator=(ResidualDipolarCoupling const & other) {
	basic::datacache::CacheableData::operator=(other);
	All_RDC_lines_ = other.All_RDC_lines_;
	release_buffers();
	preprocess_data();
	reserve_buffers();
	return *this;
}

void ResidualDipolarCoupling::read_RDC_file( Size expid, std::string const& filename ) {

	std::string line;
	utility::io::izstream infile(filename.c_str());
	if ( !infile.good() ) {
		throw( utility::excn::EXCN_FileNotFound( filename ) );
	}
	//	std::cout << "Reading RDC file " << filename << std::endl;
	tr.Info << "Reading RDC file " << filename << std::endl;
	Size lines_previously_read( All_RDC_lines_.size() );
	while (getline(infile, line)) {
		std::istringstream line_stream(line);
		std::string atom1, atom2;
		Size res1, res2;
		Real Jdipolar;
		line_stream >> res1 >> atom1 >> res2 >> atom2 >> Jdipolar;

		if ( atom1 == "HN" ) atom1 = "H"; //take care of typical NMR community notation
		if ( atom2 == "HN" ) atom2 = "H";
		if ( res1 == 1 && atom1 == "H" ) atom1 == "H1"; //or should it be ignored ?
		if ( res2 == 1 && atom2 == "H" ) atom2 == "H1";
		if ( line_stream.fail() ) {
			tr.Error << "couldn't read line " << line << " in rdc-file " << filename << std::endl;
			throw( utility::excn::EXCN_BadInput(" invalid line "+line+" in rdc-file "+filename));
		}

		Real weight(1.0);
		line_stream >> weight;
		if (line_stream.fail()) {
			tr.Debug << " set weight for RDC " << res1 << " to 1.0 ";
			weight = 1.0;
		}
		//			if ( core.size() == 0 || ( core.has( res1 ) && core.has( res2 )  ) ) {
		tr.Debug << "... added to dataset!";
		All_RDC_lines_.push_back(RDC(res1, atom1, res2, atom2, Jdipolar,
				weight, expid - 1 /*C-style counting*/));
		//			}
		tr.Debug << std::endl;
	}
	if ( All_RDC_lines_.size() == lines_previously_read ) {
		tr.Error << "file empty ? " << std::endl;
		throw( utility::excn::EXCN_BadInput(" no valid RDCs found in file " + filename + "\n" ));
	}
	lines_previously_read = All_RDC_lines_.size();


}

void ResidualDipolarCoupling::read_RDC_file() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !option[ OptionKeys::in::file::rdc ].user() ) {
		tr.Warning << "no RDC file specified" << std::endl;
		return;
	}

	nex_ = 0;
	//	protocols::loops::Loops core;
	//	if ( option[ OptionKeys::rdc::select_residues_file ].user() ) {
	//		core.read_loop_file( option[ OptionKeys::rdc::select_residues_file ](), false, "RIGID" );
	//	}
	for ( Size expid = 1; expid <= option[ OptionKeys::in::file::rdc ]().size(); expid++ ) {
		nex_ = expid;
		std::string filename( option[ OptionKeys::in::file::rdc ]()[expid] );
		read_RDC_file( expid, filename);
	}



	//extra weight file?
	if (option[OptionKeys::rdc::weights].user()) {
		std::string filename(option[OptionKeys::rdc::weights]().name());
		std::ifstream infile(filename.c_str());
		std::string line;
		while (getline(infile, line)) {
			std::istringstream line_stream(line);
			Real weight;
			Size res1;
			line_stream >> res1 >> weight;
			if (line_stream.fail()) {
				tr.Error << "[Error] reading rdc-weight-file " << filename
						<< std::endl;
				throw(utility::excn::EXCN_BadInput(" invalid line " + line
						+ " in rdc-weight-file " + filename));
			}
			for (RDC_lines::iterator it = All_RDC_lines_.begin(); it
					!= All_RDC_lines_.end(); ++it) {
				if (it->res1() == res1) {
					it->weight(weight);
					tr.Info << "set tensor-weight for RDCs to " << weight
							<< " at residue (res1) " << res1 << std::endl;
				}
			}
		}
	}

	release_buffers();
	preprocess_data();
	reserve_buffers();
}

ResidualDipolarCoupling::~ResidualDipolarCoupling() {
	release_buffers();
}

void ResidualDipolarCoupling::reserve_buffers() {
	Size nex(nex_);
	Size nrows(nrows_);
	if (nex_ < 1)
		nex = 1; //keep always 1 around so that pointers are never empty
	if (nrows_ < 1)
		nrows = 1;
	tr.Trace << "reserve buffers for nex: " << nex << " and nrows " << nrows
			<< std::endl;
	D_ = new rvec5[nrows];
	rhs_ = new rvec5[nex];
	T_ = new Tensor5[nex];
	S_ = new Tensor[nex];
	EV_ = new rvec[nex];
	EIG_ = new Tensor[nex];
	SD_ = new Tensor[nex];
	FA_ =new core::Real[nex];
	trace_=new core::Real[nex];
	maxz_=new core::Real[nex];
}

void ResidualDipolarCoupling::release_buffers() {
	delete[] D_;
	delete[] rhs_;
	delete[] T_;
	delete[] S_;
	delete[] SD_;
	delete[] EV_;
	delete[] FA_;
	delete[] trace_;
	delete[] maxz_;
	delete[] EIG_;
}

//initialize local buffers ( S, T, etc. )
//call whenever All_RDC_lines_ is changed
void ResidualDipolarCoupling::preprocess_data() {

	//find highest expid
	nex_ = 0;
	utility::vector1<core::scoring::RDC>::const_iterator it;
	for (it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it) {
		if (it->expid() >= nex_)
			nex_ = it->expid() + 1; //C-counting
	}

	nrows_ = All_RDC_lines_.size();
}

std::string element_string(std::string atom) {
	if (atom == "HN" || atom == "H" || atom == "HA")
		return "H";
	if (atom == "C" || atom == "CA")
		return "C";
	if (atom == "N")
		return "N";
	throw(utility::excn::EXCN_BadInput("unknown atom for RDC: " + atom));
	return ""; //to make compile happy.
}
//////////////////////////////////////////////////////
//@brief returns type of RDC data N-H, HA-CA etc
//////////////////////////////////////////////////////
RDC::RDC_TYPE RDC::get_RDC_data_type(std::string const & atom1,
		std::string const & atom2) {
	std::string elem1(element_string(atom1));
	std::string elem2(element_string(atom2));

	RDC_TYPE RDC_type;
	if ((elem1 == "N" && elem2 == "H") || (elem1 == "H" && elem2 == "N"))
		RDC_type = RDC_TYPE_NH;
	else if ((elem1 == "C" && elem2 == "H") || (elem1 == "H" && elem2 == "C"))
		RDC_type = RDC_TYPE_CH;
	else if ((elem1 == "C" && elem2 == "N") || (elem1 == "N" && elem2 == "C"))
		RDC_type = RDC_TYPE_NC;
	else if ((elem1 == "C" && elem2 == "C"))
		RDC_type = RDC_TYPE_CC;
	else
		throw(utility::excn::EXCN_BadInput(
				"unknown combination of atoms for RDC " + atom1 + " " + atom2));
	return RDC_type;
}

typedef core::Real Tensor[3][3];
typedef core::Real Tensor5[5][5];
typedef core::Real rvec[3];
int m_inv_gen(Tensor5 m, int n, Tensor5 minv);
void jacobi(Real a[5][5], Real d[], Real v[5][5], int *nrot);
void jacobi3(Real a[3][3], Real d[], Real v[3][3], int *nrot);

//void sort_rvec(rvec rv);
#define XX 0
#define YY 1
#define ZZ 2
typedef core::Real matrix[3][3];
typedef core::Real rvec[3];
typedef core::Real rvec5[5];
inline Real iprod(const rvec a, const rvec b) {
	return (a[XX] * b[XX] + a[YY] * b[YY] + a[ZZ] * b[ZZ]);
}
inline void mvmul(matrix a, const rvec src, rvec dest) {
	dest[XX] = a[XX][XX] * src[XX] + a[XX][YY] * src[YY] + a[XX][ZZ] * src[ZZ];
	dest[YY] = a[YY][XX] * src[XX] + a[YY][YY] * src[YY] + a[YY][ZZ] * src[ZZ];
	dest[ZZ] = a[ZZ][XX] * src[XX] + a[ZZ][YY] * src[YY] + a[ZZ][ZZ] * src[ZZ];
}


inline int compare_by_abs(const void * a, const void * b)
{
	//	int a_i = (int*) a;
//	int b_i = (const int*) b;
//	std::cout << std::abs(*(int*)a) << ' '<< std::abs(*(int*)b) << ' ' <<std::abs(*(core::Real*)a) - std::abs(*(core::Real*)b) <<std::endl;
//	std::cout << "james_hack: " << a_i << ' ' << b_i << std::endl;
//	int temp = std::abs(*(int*)a) - std::abs(*(int*)b) ;
//	std::cout <<temp <<std::endl;

	if ( *(core::Real*)a ==  *(core::Real*)b ) {
		return 0;
	} else if ( std::abs( *(core::Real*)a) - std::abs( *(core::Real*)b) > 0){
		return 1;
	} else {return -1;}

}


Real ResidualDipolarCoupling::iterate_tensor_weights(
		core::pose::Pose const& pose, core::Real sigma2, core::Real tolerance,
		bool reset) {
	utility::vector1<core::scoring::RDC>::iterator it;
	if (reset) {
		for (it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it) {
			it->weight(1.0);
		}
	}

	Real wsum_old = 100;
	Real wsum(1.0);
	Real invn(1.0 / All_RDC_lines_.size());
	tr.Debug << "wRDC: iter  wsum   wRDC" << std::endl;
	Size ct(0);
	Real score(999);
	while (std::abs(wsum_old - wsum) > tolerance) {
		score = compute_dipscore(pose);
		wsum_old = wsum;
		wsum = 0.0;
		Real wMSD = 0.0;
		for (it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it) {
			Real dev = it->Jcomputed() - it->Jdipolar();
			it->weight(exp(-dev * dev / sigma2));
			wsum += it->weight() * invn;
			wMSD += it->weight() * dev * dev;
		}
		tr.Debug << "wRMSD: " << ++ct << " " << wsum << " " << sqrt(wMSD)
				<< std::endl;
		if (ct > 200) {
			tr.Warning << "no convergence in wRDC aka iterate_tensor_weights"
					<< std::endl;
			return score;
			break;
		}
	}
	return score;
}

Real ResidualDipolarCoupling::compute_dipscore(core::pose::Pose const& pose) {
	if ( nex_ == 0 || nrows_ == 0 ) return 0;

	//taken from gromacs/gmxlib/orires.c
	//ref: B. Hess and RM Scheek, Journal of Magnetic Resonance 164 ( 2003 ) 19-27
	utility::vector1<core::scoring::RDC>::const_iterator it;
	bool const
			correct_NH(
					basic::options::option[basic::options::OptionKeys::rdc::correct_NH_length]);
	bool const
			bReduced(
					basic::options::option[basic::options::OptionKeys::rdc::reduced_couplings]);
	Size nrow(0);
	for (it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it) {
		if ( it->res1() > pose.total_residue() || it->res2() > pose.total_residue() ) {
			if ( tr.Debug.visible() ) tr.Debug << "non-existing residue, ignore RDC" << std::endl;
			continue;
		}

		//check for cutpoints!!!
		kinematics::FoldTree const& ft(pose.fold_tree());
		if ((ft.is_cutpoint(std::min((int) it->res1(), (int) it->res2())))
				&& it->res1() != it->res2()) {
			if ( tr.Trace.visible() ) tr.Trace << "cutpoint: ignore RDC " << *it << std::endl;
			continue;
		}


		++nrow;
		numeric::xyzVector<Real> r(
				pose.residue(it->res1()).atom(it->atom1()).xyz()
						- pose.residue(it->res2()).atom(it->atom2()).xyz());

		core::Real r2 = r.norm_squared();
		if (it->type() == RDC::RDC_TYPE_NH && correct_NH)
			r2 = 1.04 * 1.04;
		core::Real invr = 1.0 / sqrt(r2);

		core::Real pfac = it->Dconst() * invr * invr;
		bool bCSA(false);// hook up for later... to compute chemical shift anisotropy
		if (!bCSA) {
			pfac *= invr * invr * invr;
		}
		Size const d(nrow - 1);
		D_[d][0] = 3* pfac * (2* r [0] * r[0] + r[1] * r[1] - r2);
		D_[d][1] = 3* pfac * (2* r [0] * r[1]);
		D_[d][2] = 3* pfac * (2* r [0] * r[2]);
		D_[d][3] = 3* pfac * (2* r [1] * r[1] + r[0] * r[0] - r2);
		D_[d][4] = 3* pfac * (2* r [1] * r[2]);
		// 			Real umn_x = umn.x()/it->fixed_dist();
		// 			Real umn_y = umn.y()/it->fixed_dist();
		// 			Real umn_z = umn.z()/it->fixed_dist();
	} //cycle over atoms

	/* Calculate the order tensor S for each experiment via optimization */
	for (core::Size ex = 0; ex < nex_; ex++) {
		for (Size i = 0; i < 5; i++) {
			rhs_[ex][i] = 0;
			for (Size j = 0; j <= i; j++) {
				T_[ex][i][j] = 0;
			}
		}
	}

	for (core::Size d = 0; d < nrow; d++) {

		//type   = forceatoms[fa];
		core::Size ex = All_RDC_lines_[d + 1].expid(); //only one experiment now
		core::Real weight = All_RDC_lines_[d + 1].weight(); //force constant
		core::Real obs = All_RDC_lines_[d + 1].Jdipolar();
		/* Calculate the vector rhs and half the matrix T for the 5 equations */
		for (Size i = 0; i < 5; i++) {
			rhs_[ex][i] += D_[d][i] * obs * weight;
			for (Size j = 0; j <= i; j++)
				T_[ex][i][j] += D_[d][i] * D_[d][j] * weight;
		}
	}

	/* Now we have all the data we can calculate S */
	runtime_assert( nex_ < 200 );
	Real Smax[200];
	for (Size ex = 0; ex < nex_; ex++) {
		/* Correct corrfac and copy one half of T to the other half */
		//		std::cout << " T_[ex][j][i]: ";
		for (Size i = 0; i < 5; i++) {
			for (Size j = 0; j < i; j++) {
				T_[ex][j][i] = T_[ex][i][j];
				//			std::cout << " " << T_[ex][j][i] ;
			}
		}
		//		std::cout << std::endl;
		try {
			m_inv_gen(T_[ex], 5, T_[ex]);
//			std::cout << "performing SVD decomposition" << std::endl;
		} catch (utility::excn::EXCN_BadInput &excn) {
			if (tr.Debug) {
				pose.dump_pdb("failed_jacobi.pdb");
			}
			throw excn;
		}
		/* Calculate the orientation tensor S for this experiment */
		S_[ex][0][0] = 0;
		S_[ex][0][1] = 0;
		S_[ex][0][2] = 0;
		S_[ex][1][1] = 0;
		S_[ex][1][2] = 0;
		for (Size i = 0; i < 5; i++) {
			S_[ex][0][0] += 1.5 * T_[ex][0][i] * rhs_[ex][i];
			S_[ex][0][1] += 1.5 * T_[ex][1][i] * rhs_[ex][i];
			S_[ex][0][2] += 1.5 * T_[ex][2][i] * rhs_[ex][i];
			S_[ex][1][1] += 1.5 * T_[ex][3][i] * rhs_[ex][i];
			S_[ex][1][2] += 1.5 * T_[ex][4][i] * rhs_[ex][i];
		}
		S_[ex][1][0] = S_[ex][0][1];
		S_[ex][2][0] = S_[ex][0][2];
		S_[ex][2][1] = S_[ex][1][2];
		S_[ex][2][2] = -S_[ex][0][0] - S_[ex][1][1];
		Smax[ex] = sqrt(sqr(S_[ex][0][0]) + sqr(S_[ex][0][1]) + sqr(
				S_[ex][0][2]) + sqr(S_[ex][1][1]) + sqr(S_[ex][1][2]));


		//	std::cout << "AL.TENSOR (molecular frame): " << S_[ex][0][0] << ' ' <<  S_[ex][0][1]  << ' ' << 	S_[ex][0][2] << ' ' << S_[ex][1][1] << ' ' << S_[ex][1][2] << std::endl;



		{
			using namespace basic::options;

			if ( option[ OptionKeys::rdc::fix_normAzz ].user() ) {
				if ( option[ OptionKeys::rdc::fix_normAzz ]().size() != nex_ ) {
					utility_exit_with_message("fix_normAzz must have one value for each alignment medium !");
				}
				Smax[ ex ] = option[ OptionKeys::rdc::fix_normAzz ]()[ ex+1 ];
			}
		}
		if ( tr.Debug.visible() ) tr.Debug << "Smax( " << ex << " ): " << Smax[ex] << std::endl;

		//		std::cout << 'AL.TENSOR elements : ' << S_[ex][0][0]) << ' ' <<  S_[ex][0][1]  << ' ' << 	S_[ex][0][2] << ' '
		//				<< S_[ex][1][1] << ' ' << S_[ex][1][2] << std::endl;



	} // for ( ex = 0 .. nex )



	Real wsv2 = 0;
	Real sw = 0;
	Real two_thr = 2.0 / 3.0;
	Real vtot = 0;
	Real Q = 0;
	Real Qnorm = 0;
	core::Real const Rohl2Hess(2.5);

	Size irow(0);
	for (utility::vector1<core::scoring::RDC>::iterator it =
			All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it) {

		//check for cutpoints!!!
		kinematics::FoldTree const& ft(pose.fold_tree());
		if ((ft.is_cutpoint(std::min((int) it->res1(), (int) it->res2())))
				&& it->res1() != it->res2()) {
			if ( tr.Trace.visible() ) tr.Trace << "cutpoint: ignore RDC " << *it << std::endl;
			continue;
		}

		++irow;
		Size const d(irow - 1);

		//for( Size d = 0; d<nrow; d++ ) {
		Size ex = it->expid(); //exp_id const 1 ... fix later

		Real computed_coupling = it->Jdipolar_computed_ = two_thr
				* (S_[ex][0][0] * D_[d][0] + S_[ex][0][1] * D_[d][1]
						+ S_[ex][0][2] * D_[d][2] + S_[ex][1][1] * D_[d][3]
						+ S_[ex][1][2] * D_[d][4]);

		//		pfac  = fc*ip[type].orires.c*invr2;
		//for(i=0; i<power; i++)
		//    pfac *= invr;
		Size const power(3); //this will be 0 for CSA see above
		RDC& rdc = *it;
		numeric::xyzVector<Real> r(
				pose.residue(rdc.res1()).atom(rdc.atom1()).xyz()
						- pose.residue(rdc.res2()).atom(rdc.atom2()).xyz());
		core::Real r2 = r.norm_squared();
		core::Real invr = 1.0 / sqrt(r2);
		core::Real invr2 = sqr(invr);
		core::Real pfac = rdc.Dconst() * invr2 * invr2 * invr * Rohl2Hess
				/ (Smax[ex] * Smax[ex]) / All_RDC_lines_.size();

		core::Real const pfac_NH = 6.088 / 1.05 * Rohl2Hess / (Smax[ex]
				* Smax[ex]) / All_RDC_lines_.size();
		Real obs = it->Jdipolar();
		Real dev = computed_coupling - obs;
		if (bReduced) {
			if ( tr.Trace.visible() ) tr.Trace << "reducing coupling for " << rdc << " dev: " << dev
																			 << " pfac: " << pfac << " pfac_NH " << pfac_NH;
			dev *= pfac_NH / pfac;
			obs *= pfac_NH / pfac;
			if ( tr.Trace.visible() ) tr.Trace << " new dev: " << dev << std::endl;
		}

		rvec Sr;
		rvec rgmx;
		rgmx[0] = r[0];
		rgmx[1] = r[1];
		rgmx[2] = r[2];
		mvmul(S_[ex], rgmx, Sr);

		if (bReduced)
			pfac = pfac_NH;
		for (Size i = 0; i < 3; i++) {
			rdc.fij_[i] = -pfac * dev * (4* Sr [i] - 2* (2 + power ) * invr2 *
			iprod (Sr ,rgmx) * rgmx[ i]);
		}

		Real weight = it->weight(); //force constant
		//	std::cout << "WEIGHT " << weight <<std::endl;
		vtot += 0.5*sqr( dev )/( Smax[ ex ] * Smax[ ex ] ); //*weight if we want that
		//      vtot += sqrt( dev );
		wsv2 += weight*sqr(dev);
		sw += weight;
		Q += sqr( dev );
		Qnorm += sqr( obs );
	}
	R_ = sqrt( Q/Qnorm/2 );
	rmsd_ = sqrt(wsv2/sw);
	//	Real rohl_score = wsv2 / Azz; //Azz axial component of tensor?

	//  /* Rotate the S matrices back, so we get the correct grad(tr(S D)) */
	//   for(ex=0; ex<od->nex; ex++) {
	//     tmmul(R,S_[ex],TMP);
	//     mmul(TMP,R,S_[ex]);
	//   }
	//this is a factor that brings the score roughly back in the range where it was
	//with the old implementation: i.e., norm by Azz*Azz and deviation of "reduced" couplings instead of the real ones.
	//std::cout  << "DEV_TOT " << vtot <<std::endl;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::fmt;
	if ( option[ OptionKeys::rdc::print_rdc_values ].user() ) {
		try{
			compute_tensor_stats();
		} catch (utility::excn::EXCN_BadInput &excn) {
			if ( tr.Debug.visible() ) {
				pose.dump_pdb("failed_jacobi.pdb");
			}
			throw excn;
		}
		std::string filename( option[ OptionKeys::rdc::print_rdc_values ]() );
		utility::io::ozstream out;

		out.open_append( filename ) ;
		using namespace core::pose::datacache;
		Size const width( 8 );
		Size const width_large(6);
		// mjo comment out precision because it is not used and causes a warning.
		//Size const precision( 2 );
		std::string tag( core::pose::tag_from_pose(pose) );
		for (Size ex = 0; ex < nex_; ex++) {
			Real Smax = sqrt(sqr(S_[ex][0][0]) + sqr(S_[ex][0][1]) + sqr(S_[ex][0][2]) + sqr(S_[ex][1][1]) + sqr(S_[ex][1][2]));
			out << A( width_large, "TAG ")   << A( width, tag ) << A( width, "EXP       ")   << I( width, ex ) << I( width, Smax )
					<< I( width, vtot )<< I( width, Rohl2Hess )<< std::endl;
			show_tensor_stats( out, ex );
			show_tensor_matrix( out, ex );
			show_rdc_values( out, ex );
		}
		out << "//" <<std::endl;
		out.close();
		//std::cout << "RDC values " << obs << ' ' << computed_coupling <<std::endl;
	}

	return Rohl2Hess*vtot/All_RDC_lines_.size();
}

void ResidualDipolarCoupling::show_tensor_stats( std::ostream& out, Size ex ) const {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::fmt;

	Real Da = EV_[ex][0]/2;
	Real Dr = (EV_[ex][1]-EV_[ex][2])/3;
// 	Real Aa = Da*2;
// 	Real Ar = Dr*2;
	Size const width( 10 );
	Size const precision( 2 );
	Real rhombicity = Dr/Da;
	out << A( width, "Da" ) << " " << A( width, "Dr" ) << A( width, "rhombicity" ) << std::endl;
	out << F( width, precision, Da ) << F( width, precision, Dr ) << F( width, precision, rhombicity ) << std::endl;
}

void ResidualDipolarCoupling::show_tensor_matrix( std::ostream& out, Size ex ) const {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::fmt;

	Size const width( 8 );
	//mjo comment out width_large because it is unused and causes a warning
	//Size const width_large(6);
	Size const precision( 2 );
	out << A( width, "AL.EVAL   ")   << F(width, precision, EV_[ex][0]) <<  F(width, precision, EV_[ex][1])  <<  F(width, precision, EV_[ex][2])  << std::endl
			<< A( width, "AL.EVEC XX")   << F(width, precision, EIG_[ex][0][0]) <<  F(width, precision, EIG_[ex][1][0])  <<  F(width, precision, EIG_[ex][2][0])  << std::endl
			<< A( width, "AL.EVEC YY")   << F(width, precision, EIG_[ex][0][1]) <<  F(width, precision, EIG_[ex][1][1])  <<  F(width, precision, EIG_[ex][2][1])  << std::endl
			<< A( width, "AL.EVEC ZZ")   << F(width, precision, EIG_[ex][0][2]) <<  F(width, precision, EIG_[ex][1][2])  <<  F(width, precision, EIG_[ex][2][2])  << std::endl;
}

void ResidualDipolarCoupling::show_rdc_values( std::ostream& out, Size ex ) const {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::fmt;

	Size const width( 8 );
	//mjo comment out width_large because it is unused and causes a warning
	//Size const width_large(6);
	Size const precision( 2 );
	utility::vector1<core::scoring::RDC>::const_iterator it;
	for (it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it) {
		if(it->expid() == ex){
			Size count( it->res1() );
			std::string type;
			if ( it->type() == RDC::RDC_TYPE_NH ) {type ="NH";}
			if ( it->type() == RDC::RDC_TYPE_CH ) {type ="CH";}
			if ( it->type() == RDC::RDC_TYPE_NC ) {type ="NC";}
			if ( it->type() == RDC::RDC_TYPE_CC ) {type ="CC";}
			Real rdc_exp( it->Jdipolar() );
			Real rdc_computed ( it->Jcomputed());
			out   << A( width, "RDC" )
						<< A( width, type )
						<< I( width, count )
						<< F( width, precision, rdc_exp )
						<< F( width, precision, rdc_computed )
						<< std::endl;
		}
	}
}

void ResidualDipolarCoupling::compute_tensor_stats() {
	for (Size ex = 0; ex < nex_; ex++) {
		// This section performs diagonlization of the recently optimized alignment tensor (NGS)

    SD_[ex][0][0] =  S_[ex][0][0];
		SD_[ex][0][1] =  S_[ex][0][1];
		SD_[ex][0][2] =  S_[ex][0][2];
		SD_[ex][1][1] =  S_[ex][1][1];
		SD_[ex][1][2] =  S_[ex][1][2];
		SD_[ex][1][0] =  S_[ex][1][0];
		SD_[ex][2][0] =  S_[ex][2][0];
		SD_[ex][2][1] =  S_[ex][2][1];
		SD_[ex][2][2] =  S_[ex][2][2];


		//	std::cout << "AL.TENSOR (molecular frame): " << SD_[ex][0][0] << ' ' <<  SD_[ex][0][1]  << ' ' << 	SD_[ex][0][2] << ' ' << SD_[ex][1][1] << ' ' << SD_[ex][1][2] << std::endl;


  	Tensor v2;
		int nrot2;
		rvec ev2;

		jacobi3(SD_[ex],ev2,v2,&nrot2);

		EIG_[ex][0][0] =  v2[0][0];
		EIG_[ex][0][1] =  v2[0][1];
		EIG_[ex][0][2] =  v2[0][2];
		EIG_[ex][1][1] =  v2[1][1];
		EIG_[ex][1][2] =  v2[1][2];
		EIG_[ex][1][0] =  v2[1][0];
		EIG_[ex][2][0] =  v2[2][0];
		EIG_[ex][2][1] =  v2[2][1];
		EIG_[ex][2][2] =  v2[2][2];


		EV_[ex][0] = ev2[2];
		EV_[ex][1] = ev2[1];
		EV_[ex][2] = ev2[0];
		qsort(ev2, 3 ,sizeof(core::Real),compare_by_abs);
		//		std::cout << "AL_TENSOR eigenvalues : " << EV_[ex][0] << ' ' << EV_[ex][1]  << ' ' << 	EV_[ex][2]  << ' ' << nrot2 << std::endl;

		/*std::cout << "AL_TENSOR eigenvectors  : " << std::endl;
			std::cout << v2[0][0] << "\t" << v2[1][0]  << "\t" << 	v2[2][0]  << std::endl;
			std::cout << v2[0][1] << "\t" << v2[1][1]  << "\t" << 	v2[2][1]  << std::endl;
			std::cout << v2[0][2] << "\t" << v2[1][2]  << "\t" << 	v2[2][2]  << std::endl;*/
		trace_[ex] = v2[0][0] + v2[1][1] + v2[2][2];
		rvec tempvec;
		tempvec[0] =  v2[0][2];
		tempvec[1] =  v2[1][2];
		tempvec[2] =  v2[2][2];
		qsort(tempvec, 3, sizeof(core::Real),compare_by_abs);
		maxz_[ex] = tempvec[0];
		//	Real ev_trace =  (ev2[2] + ev2[1] +ev2[0]) / 3;
		//		FA_[ex] = 4.0;
		//			FA_[ex] = sqrt (0.5) *  sqrt( (ev2[0] - ev_trace)*(ev2[0] - ev_trace) +  (ev2[1] - ev_trace)*(ev2[1] - ev_trace) +  (ev2[2] - ev_trace)*(ev2[2] - ev_trace) ) / sqrt (ev2[0] * ev2[0] + ev2[1] * ev2[1] + ev2[2] * ev2[2]);
		//	std::cout << "Fractional Anisotropy " << FA_[ex] << std::endl;

		//	std::cout << "AL.TENSOR (molecular frame after): " << S_[ex][0][0] << ' ' <<  S_[ex][0][1]  << ' ' << 	S_[ex][0][2] << ' ' << S_[ex][1][1] << ' ' << S_[ex][1][2] << std::endl;

		//	delete[] v2;
		//	delete[] ev2;

	} //for ( ex = 0 .. nex )
}

int m_inv_gen(Tensor5 m,int n,Tensor5 minv)
{
	Tensor5 md,v;
	rvec5 eig;
	Real tol,s;
	int nzero,i,j,k,nrot;
	//   md = new Real*[n];
				// 	v = new Real*[n];
				//   for(i=0; i<n; i++) {
				//     md[i] = new Real[n];
				//     v[i] = new Real[n];
				// 	}
				// 	eig = new Real[ n ];
				for(i=0; i<n; i++)
				for(j=0; j<n; j++)
				md[i][j] = m[i][j];

				tol = 0;
				for(i=0; i<n; i++)
				tol += fabs(md[i][i]);
				tol = 1e-6*tol/n;

				jacobi(md,eig,v,&nrot);

				nzero = 0;
				for(i=0; i<n; i++)
				if (fabs(eig[i]) < tol) {
					eig[i] = 0;
					nzero++;
				} else
				eig[i] = 1.0/eig[i];

				for(i=0; i<n; i++)
				for(j=0; j<n; j++) {
					s = 0;
					for(k=0; k<n; k++)
					s += eig[k]*v[i][k]*v[j][k];
					minv[i][j] = s;
				}

				//   delete[] eig;
				//   for(i=0; i<n; i++) {
				//     delete[] v[i];
				//     delete[] md[i];
				// 	}

				//   delete[] v;
				// 	delete[] md;

				return nzero;
			}

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
  a[k][l]=h+s*(g-h*tau);

					void jacobi(Real a[5][5],Real d[],Real v[5][5],int *nrot)
					{
						int j,i;
						int iq,ip;
						Real tresh,theta,tau,t,sm,s,h,g,c;
						Real b[5];
						Real z[5];
						int const n( 5 );
						for (ip=0; ip<n; ip++) {
							for (iq=0; iq<n; iq++) v[ip][iq]=0.0;
							v[ip][ip]=1.0;
						}
						for (ip=0; ip<n;ip++) {
							b[ip]=d[ip]=a[ip][ip];
							z[ip]=0.0;
						}
						*nrot=0;
						for (i=1; i<=50; i++) {
							sm=0.0;
							for (ip=0; ip<n-1; ip++) {
								for (iq=ip+1; iq<n; iq++)
								sm += fabs(a[ip][iq]);
							}
							if (sm == 0.0) {
								return;
							}
							if (i < 4) //first 3 iterations
							tresh=0.2*sm/(n*n);
							else
							tresh=0.0;
							for (ip=0; ip<n-1; ip++) {
								for (iq=ip+1; iq<n; iq++) {
									g=100.0*fabs(a[ip][iq]);
									if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
					&& fabs(d[iq])+g == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if (fabs(h)+g == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0; j<ip; j++) {
						ROTATE(a,j,ip,j,iq)
							}
					for (j=ip+1; j<iq; j++) {
						ROTATE(a,ip,j,j,iq)
							}
					for (j=iq+1; j<n; j++) {
						ROTATE(a,ip,j,iq,j)
							}
					for (j=0; j<n; j++) {
						ROTATE(v,j,ip,j,iq)
							}
					++(*nrot);
				}
      }
    }
    for (ip=0; ip<n; ip++) {
      b[ip] +=  z[ip];
      d[ip]  =  b[ip];
      z[ip]  =  0.0;
    }
  }
	//probably different type of Exception is better suited
	throw( utility::excn::EXCN_BadInput(" too many iterations in Jacobi when compute RDC tensor") );
}


		void jacobi3(Real a[3][3],Real d[],Real v[3][3],int *nrot)
					{
						int j,i;
						int iq,ip;
						Real tresh,theta,tau,t,sm,s,h,g,c;
						Real b[3];
						Real z[3];
						int const n( 3 );
						for (ip=0; ip<n; ip++) {
							for (iq=0; iq<n; iq++) v[ip][iq]=0.0;
							v[ip][ip]=1.0;
						}
						for (ip=0; ip<n;ip++) {
							b[ip]=d[ip]=a[ip][ip];
							z[ip]=0.0;
						}
						*nrot=0;
						for (i=1; i<=50; i++) {
							sm=0.0;
							for (ip=0; ip<n-1; ip++) {
								for (iq=ip+1; iq<n; iq++)
								sm += fabs(a[ip][iq]);
							}
							if (sm == 0.0) {
								return;
							}
							if (i < 4) //first 3 iterations
							tresh=0.2*sm/(n*n);
							else
							tresh=0.0;
							for (ip=0; ip<n-1; ip++) {
								for (iq=ip+1; iq<n; iq++) {
									g=100.0*fabs(a[ip][iq]);
									if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
					&& fabs(d[iq])+g == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if (fabs(h)+g == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0; j<ip; j++) {
						ROTATE(a,j,ip,j,iq)
							}
					for (j=ip+1; j<iq; j++) {
						ROTATE(a,ip,j,j,iq)
							}
					for (j=iq+1; j<n; j++) {
						ROTATE(a,ip,j,iq,j)
							}
					for (j=0; j<n; j++) {
						ROTATE(v,j,ip,j,iq)
							}
					++(*nrot);
				}
      }
    }
    for (ip=0; ip<n; ip++) {
      b[ip] +=  z[ip];
      d[ip]  =  b[ip];
      z[ip]  =  0.0;
    }
  }
	//probably different type of Exception is better suited
	throw( utility::excn::EXCN_BadInput(" too many iterations in Jacobi when compute RDC tensor") );
}


/*void sort_rvec(rvec rv){

 int round(0);
 int winner(0);
 rvec rv_s;
 rvec indeces = {1.0,1.0,1.0};
 while (round < 2){

	 for (int i = 0; i < 3; ++i) {
		 for(int j = i; j < 3; ++j){
			 if( abs(rv[j]) > abs(rv[i]) && indeces[i] == 1.0 && indeces[j] ==1.0 ) { winner = j; }
		}

		//	std::cout << array[i] <<std::endl;
    }
 rv_s[round] = rv[winner];
 indeces[winner] = 0.0;
	++round;
	winner = 0;
 }
 for(int i=0; i <3; ++i){if(indeces[i])rv_s[2]=rv[i];}

 for(int i=0; i <3; ++i){rv[i]=rv_s[i];}

 }*/
} //namespace Scoring
} //namespace core
