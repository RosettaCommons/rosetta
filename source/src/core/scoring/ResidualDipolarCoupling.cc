// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ResidualDipolarCoupling.cc
/// @brief  Uses NMR RDC for scoring
/// @author Oliver Lange
/// @author Srivatsan Raman
/// @author Nikolas Sgourakis
/// @author Lei Shi (nls and modifications to original code)

//Unit headers
#include <core/scoring/ResidualDipolarCoupling.hh>

// Package headers

// Project headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>

#include <basic/options/option.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/Tracer.hh>

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

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rdc.OptionKeys.gen.hh>

#include <numeric/nls/lmmin.hh>

#include <utility/vector1.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/NumericTraits.hh>

static basic::Tracer tr( "core.scoring.ResidualDipolarCoupling" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector0.srlz.hh>
#include <utility/vector1.srlz.hh>
#include <utility/fixedsizearray0.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

#define NLS 0
#define NLSDA 1
#define NLSR 2
#define NLSDAR 3

inline Real sqr(Real x) {
	return x * x;
}

Real const ResidualDipolarCoupling::COMMON_DENOMINATOR = 36.5089/1.041/1.041/1.041;


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
	if ( pose.data().has(core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA) ) {
		return utility::pointer::static_pointer_cast< core::scoring::ResidualDipolarCoupling const > ( pose.data().get_const_ptr(
			core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA) );
	};
	return nullptr;
}

extern ResidualDipolarCouplingOP retrieve_RDC_from_pose(core::pose::Pose& pose) {
	// ////using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA;
	if ( pose.data().has(core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA) ) {
		return utility::pointer::static_pointer_cast< core::scoring::ResidualDipolarCoupling > ( pose.data().get_ptr(
			core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA) );
	};
	return nullptr;
}

void ResidualDipolarCoupling::show(std::ostream& out) const {
	Size ct=0;
	for ( auto const & All_RDC_line : All_RDC_lines_ ) {
		out << "RDC "<<++ct << "     ";
		out << All_RDC_line << std::endl;
	}
}

void RDC::show(std::ostream& out) const {
	using namespace ObjexxFCL::format;
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
	basic::datacache::CacheableData(other)
{
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
		throw CREATE_EXCEPTION(utility::excn::FileNotFound,  filename ) ;
	}
	// std::cout << "Reading RDC file " << filename << std::endl;
	tr.Info << "Reading RDC file " << filename << std::endl;
	Size lines_previously_read( All_RDC_lines_.size() );
	while ( getline(infile, line) ) {
		std::istringstream line_stream(line);
		std::string atom1, atom2;
		Size res1, res2;
		Real Jdipolar;
		line_stream >> res1 >> atom1 >> res2 >> atom2 >> Jdipolar;

		if ( atom1 == "HN" ) atom1 = "H"; //take care of typical NMR community notation
		if ( atom2 == "HN" ) atom2 = "H";

		if ( line_stream.fail() ) {
			tr.Error << "couldn't read line " << line << " in rdc-file " << filename << std::endl;
			throw CREATE_EXCEPTION(utility::excn::BadInput, " invalid line "+line+" in rdc-file "+filename);
		}

		if ( res1 < 1 || res2 < 1 ) {
			tr.Error << "negative residue number in line " << line << " in rdc-file " << filename << std::endl;
			throw CREATE_EXCEPTION(utility::excn::BadInput, " invalid line "+line+" in rdc_file " + filename );
		}

		Real weight(1.0);
		line_stream >> weight;
		if ( line_stream.fail() ) {
			tr.Debug << " set weight for RDC " << res1 << " to 1.0 ";
			weight = 1.0;
		}
		tr.Debug << "... added to dataset!";
		All_RDC_lines_.push_back(RDC(res1, atom1, res2, atom2, Jdipolar,
			weight, expid - 1 /*C-style counting*/));
		tr.Debug << std::endl;
	}
	if ( All_RDC_lines_.size() == lines_previously_read ) {
		tr.Error << "file empty ? " << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput, " no valid RDCs found in file " + filename + "\n" );
	}
	//lines_previously_read = All_RDC_lines_.size(); //unused
}

void ResidualDipolarCoupling::read_RDC_file() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !option[ OptionKeys::in::file::rdc ].user() ) {
		tr.Warning << "no RDC file specified" << std::endl;
		return;
	}

	nex_ = 0;
	for ( Size expid = 1; expid <= option[ OptionKeys::in::file::rdc ]().size(); expid++ ) {
		nex_ = expid;
		std::string filename( option[ OptionKeys::in::file::rdc ]()[expid] );
		read_RDC_file( expid, filename);
	}

	//extra weight file?
	if ( option[OptionKeys::rdc::weights].user() ) {
		std::string filename(option[OptionKeys::rdc::weights]().name());
		std::ifstream infile(filename.c_str());
		std::string line;
		while ( getline(infile, line) ) {
			std::istringstream line_stream(line);
			Real weight;
			Size res1;
			line_stream >> res1 >> weight;
			if ( line_stream.fail() ) {
				tr.Error << "reading rdc-weight-file " << filename
					<< std::endl;
				throw CREATE_EXCEPTION(utility::excn::BadInput, " invalid line " + line + " in rdc-weight-file " + filename);
			}
			for ( auto & All_RDC_line : All_RDC_lines_ ) {
				if ( All_RDC_line.res1() == res1 ) {
					All_RDC_line.weight(weight);
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
	if ( nex_ < 1 ) {
		nex = 1; //keep always 1 around so that pointers are never empty
	}
	if ( nrows_ < 1 ) {
		nrows = 1;
	}
	tr.Trace << "reserve buffers for nex: " << nex << " and nrows " << nrows << std::endl;
	D_.resize( nrows );        // D_ = new rvec5[nrows];
	rhs_.resize( nex );        // rhs_ = new rvec5[nex];
	T_.resize( nex );          // T_ = new Tensor5[nex];
	S_.resize( nex );          // S_ = new Tensor[nex];
	EV_.resize( nex );         // EV_ = new rvec[nex];
	EIG_.resize( nex );        // EIG_ = new Tensor[nex];
	SD_.resize( nex );         // SD_ = new Tensor[nex];
	FA_.resize( nex );         // FA_ =new core::Real[nex];
	trace_.resize( nex );      // trace_=new core::Real[nex];
	maxz_.resize( nex );       // maxz_=new core::Real[nex];
	r0_.resize( nrows );       // r0_=new core::Real[nrows];
	r1_.resize( nrows );       // r1_=new core::Real[nrows];
	r2_.resize( nrows );       // r2_=new core::Real[nrows];
	exprdc_.resize( nrows );   // exprdc_=new core::Real[nrows];
	rdcconst_.resize( nrows ); // rdcconst_=new core::Real[nrows];
	rdcweight_.resize( nrows );// rdcweight_=new core::Real[nrows];
	lenex_.resize( nex + 1 );  // lenex_=new core::Size[nex+1];
}

void ResidualDipolarCoupling::release_buffers() {
	D_.clear();        // delete[] D_;
	rhs_.clear();    // delete[] rhs_;
	T_.clear();     // delete[] T_;
	S_.clear();     // delete[] S_;
	SD_.clear();    // delete[] SD_;
	EV_.clear();    // delete[] EV_;
	FA_.clear();    // delete[] FA_;
	trace_.clear();   // delete[] trace_;
	maxz_.clear();   // delete[] maxz_;
	EIG_.clear();    // delete[] EIG_;
	r0_.clear();    // delete[] r0_;
	r1_.clear();    // delete[] r1_;
	r2_.clear();    // delete[] r2_;
	exprdc_.clear();  // delete[] exprdc_;
	rdcconst_.clear(); // delete[] rdcconst_;
	rdcweight_.clear();// delete[] rdcweight_;
	lenex_.clear();   // delete[] lenex_;
}

//initialize local buffers ( S, T, etc. )
//call whenever All_RDC_lines_ is changed
void ResidualDipolarCoupling::preprocess_data() {

	//find highest expid
	nex_ = 0;
	utility::vector1<core::scoring::RDC>::const_iterator it;
	for ( it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it ) {
		if ( it->expid() >= nex_ ) {
			nex_ = it->expid() + 1; //C-counting
		}
	}

	nrows_ = All_RDC_lines_.size();
}

std::string element_string(std::string atom) {
	if ( atom == "HN" || atom == "H" || atom == "HA" || atom == "1H" ) {
		return "H";
	}
	if ( atom == "C" || atom == "CA" ) {
		return "C";
	}
	if ( atom == "N" ) {
		return "N";
	}
	throw CREATE_EXCEPTION(utility::excn::BadInput, "unknown atom for RDC: " + atom);
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
	if ( (elem1 == "N" && elem2 == "H") || (elem1 == "H" && elem2 == "N") ) {
		RDC_type = RDC_TYPE_NH;
	} else if ( (elem1 == "CA" && elem2 == "HA") || (elem1 == "HA" && elem2 == "CA") ) {
		RDC_type = RDC_TYPE_CH;
	} else if ( (elem1 == "C" && elem2 == "N") || (elem1 == "N" && elem2 == "C") ) {
		RDC_type = RDC_TYPE_NC;
	} else if ( (elem1 == "C" && elem2 == "C") ) {
		RDC_type = RDC_TYPE_CC;
	} else if ( (elem1 == "C" && elem2 == "H") || (elem1 == "H" && elem2 == "C") ) {
		RDC_type = RDC_TYPE_CHN;
	} else if ( (elem1 == "N" && elem2 == "CA") || (elem1 == "CA" && elem2 == "N") ) {
		RDC_type = RDC_TYPE_NCA;
	} else {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "unknown combination of atoms for RDC " + atom1 + " " + atom2);
	}
	return RDC_type;
}

int m_inv_gen( ResidualDipolarCoupling::Tensor5 const & m, int n, ResidualDipolarCoupling::Tensor5 & minv );
void jacobi( ResidualDipolarCoupling::Tensor5 & a, ResidualDipolarCoupling::rvec5 & d, ResidualDipolarCoupling::Tensor5 & v, int & nrot );
void jacobi3( ResidualDipolarCoupling::Tensor & a, ResidualDipolarCoupling::rvec  & d, ResidualDipolarCoupling::Tensor  & v, int & nrot );

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
	// int a_i = (int*) a;
	// int b_i = (const int*) b;
	// std::cout << std::abs(*(int*)a) << ' '<< std::abs(*(int*)b) << ' ' <<std::abs(*(core::Real*)a) - std::abs(*(core::Real*)b) <<std::endl;
	// std::cout << "james_hack: " << a_i << ' ' << b_i << std::endl;
	// int temp = std::abs(*(int*)a) - std::abs(*(int*)b) ;
	// std::cout <<temp <<std::endl;

	if ( *(core::Real*)a ==  *(core::Real*)b ) {
		return 0;
	} else if ( std::abs( *(core::Real*)a) - std::abs( *(core::Real*)b) > 0 ) {
		return 1;
	} else { return -1; }
}


Real ResidualDipolarCoupling::iterate_tensor_weights(
	core::pose::Pose const& pose, core::Real sigma2, core::Real tolerance,
	bool reset) {
	utility::vector1<core::scoring::RDC>::iterator it;
	if ( reset ) {
		for ( it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it ) {
			it->weight(1.0);
		}
	}

	Real wsum_old = 100;
	Real wsum(1.0);
	Real invn(1.0 / All_RDC_lines_.size());
	tr.Debug << "wRDC: iter  wsum   wRDC" << std::endl;
	Size ct(0);
	Real score(999);
	while ( std::abs(wsum_old - wsum) > tolerance ) {
		score = compute_dipscore(pose);
		wsum_old = wsum;
		wsum = 0.0;
		Real wMSD = 0.0;
		for ( it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it ) {
			Real dev = it->Jcomputed() - it->Jdipolar();
			it->weight(exp(-dev * dev / sigma2));
			wsum += it->weight() * invn;
			wMSD += it->weight() * dev * dev;
		}
		tr.Debug << "wRMSD: " << ++ct << " " << wsum << " " << sqrt(wMSD)
			<< std::endl;
		if ( ct > 200 ) {
			tr.Warning << "no convergence in wRDC aka iterate_tensor_weights"
				<< std::endl;
			return score;
			//break;
		}
	}
	return score;
}

Real ResidualDipolarCoupling::compute_dipscore(core::pose::Pose const& pose) {
	if ( nex_ == 0 || nrows_ == 0 ) return 0;

	//taken from gromacs/gmxlib/orires.c
	//ref: B. Hess and RM Scheek, Journal of Magnetic Resonance 164 ( 2003 ) 19-27
	utility::vector1<core::scoring::RDC>::const_iterator it;
	bool const correct_NH( basic::options::option[basic::options::OptionKeys::rdc::correct_NH_length]);
	bool const bReduced( basic::options::option[basic::options::OptionKeys::rdc::reduced_couplings]);
	Size nrow(0);

	/* Calculate the order tensor S for each experiment via optimization */
	for ( core::Size ex = 0; ex < nex_; ex++ ) {
		for ( Size i = 0; i < 5; i++ ) {
			rhs_[ex][i] = 0;
			for ( Size j = 0; j <= i; j++ ) {
				T_[ex][i][j] = 0;
			}
		}
	}

	for ( it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it ) {
		if ( it->res1() > pose.size() || it->res2() > pose.size() ) {
			tr.Warning << "non-existing residue, ignore RDC" << std::endl;
			continue;
		}

		if ( it->res1() == 1 || it->res2() == 1  ) {
			tr.Warning << "residue 1, always problematic with 1H vs H, just ignore RDC" << std::endl;
			continue;
		}

		//check for cutpoints!!!
		kinematics::FoldTree const& ft(pose.fold_tree());
		if ( (ft.is_cutpoint(std::min((int) it->res1(), (int) it->res2())))
				&& it->res1() != it->res2() ) {
			tr.Warning << "cutpoint: ignore RDC " << *it << std::endl;
			continue;
		}


		++nrow;
		numeric::xyzVector<Real> r( pose.residue(it->res1()).atom(it->atom1()).xyz() - pose.residue(it->res2()).atom(it->atom2()).xyz());

		core::Real r2 = r.norm_squared();
		core::Real invr = 1.0 / sqrt(r2);
		if ( it->type() == RDC::RDC_TYPE_NH && correct_NH ) {
			r.normalize(1.041);
			r2 = 1.041 * 1.041;
			invr = 1.0 / sqrt(r2);
		}

		core::Real pfac = it->Dconst();
		bool bCSA(false);// hook up for later... to compute chemical shift anisotropy
		if ( !bCSA ) {
			pfac *= invr * invr * invr;
		}
		Size const d(nrow - 1);

		//put a length there should change it
		D_[d][0] = pfac * (2* r [0] * r[0] + r[1] * r[1] - r.norm_squared());
		D_[d][1] = pfac * (2* r [0] * r[1]);
		D_[d][2] = pfac * (2* r [0] * r[2]);
		D_[d][3] = pfac * (2* r [1] * r[1] + r[0] * r[0] - r.norm_squared());
		D_[d][4] = pfac * (2* r [1] * r[2]);

		core::Size ex = it->expid(); //only one experiment now
		core::Real weight = it->weight(); //force constant
		core::Real obs = it->Jdipolar();
		/* Calculate the vector rhs and half the matrix T for the 5 equations */
		for ( Size i = 0; i < 5; i++ ) {
			rhs_[ex][i] += D_[d][i] * obs * weight;
			for ( Size j = 0; j <= i; j++ ) {
				T_[ex][i][j] += D_[d][i] * D_[d][j] * weight;
			}
		}
	} //cycle over atoms

	/* Now we have all the data we can calculate S */
	utility::vector1<Real> Smax(nex_);
	for ( Size ex = 0; ex < nex_; ex++ ) {
		/* Correct corrfac and copy one half of T to the other half */
		//  std::cout << " T_[ex][j][i]: ";
		for ( Size i = 0; i < 5; i++ ) {
			for ( Size j = 0; j < i; j++ ) {
				T_[ex][j][i] = T_[ex][i][j];
				//   std::cout << " " << T_[ex][j][i] ;
			}
		}
		//  std::cout << std::endl;
		try {
			m_inv_gen(T_[ex], 5, T_[ex]);
			//   std::cout << "performing SVD decomposition" << std::endl;
		} catch (utility::excn::BadInput &excn) {
			if ( tr.Debug.visible() ) {
				pose.dump_pdb("failed_jacobi.pdb");
			}
			// changed from throw excn--you do not want to copy the exception you catch
			throw;
		}

		/* Calculate the orientation tensor S for this experiment */
		S_[ex][0][0] = 0;
		S_[ex][0][1] = 0;
		S_[ex][0][2] = 0;
		S_[ex][1][1] = 0;
		S_[ex][1][2] = 0;
		for ( Size i = 0; i < 5; i++ ) {
			S_[ex][0][0] += T_[ex][0][i] * rhs_[ex][i];
			S_[ex][0][1] += T_[ex][1][i] * rhs_[ex][i];
			S_[ex][0][2] += T_[ex][2][i] * rhs_[ex][i];
			S_[ex][1][1] += T_[ex][3][i] * rhs_[ex][i];
			S_[ex][1][2] += T_[ex][4][i] * rhs_[ex][i];
		}
		S_[ex][1][0] = S_[ex][0][1];
		S_[ex][2][0] = S_[ex][0][2];
		S_[ex][2][1] = S_[ex][1][2];
		S_[ex][2][2] = -S_[ex][0][0] - S_[ex][1][1];

		//Always use one, otherwise it is overwritten.
		Smax[ex] = 1;

		// std::cout << "AL.TENSOR (molecular frame): " << S_[ex][0][0] << ' ' <<  S_[ex][0][1]  << ' ' <<  S_[ex][0][2] << ' ' << S_[ex][1][1] << ' ' << S_[ex][1][2] << std::endl;

		{
			using namespace basic::options;

			if ( option[ OptionKeys::rdc::fix_normAzz ].user() ) {
				if ( option[ OptionKeys::rdc::fix_normAzz ]().size() != nex_ ) {
					utility_exit_with_message("fix_normAzz must have one value for each alignment medium !");
				}
				Smax[ ex ] = option[ OptionKeys::rdc::fix_normAzz ]()[ ex+1 ];
			}
		}

		tr.Debug << "Smax( " << ex << " ): " << Smax[ex] << std::endl;

		//  std::cout << 'AL.TENSOR elements : ' << S_[ex][0][0]) << ' ' <<  S_[ex][0][1]  << ' ' <<  S_[ex][0][2] << ' '
		//    << S_[ex][1][1] << ' ' << S_[ex][1][2] << std::endl;
	} // for ( ex = 0 .. nex )

	//Always diagonalize the matrix to find out the alignment
	try {
		compute_tensor_stats();
	} catch (utility::excn::BadInput &excn) {
		if ( tr.Debug.visible() ) {
			pose.dump_pdb("failed_jacobi.pdb");
		}
		throw;
	}

	//Compute the fitting stats
	Real wsv2 = 0;
	Real sw = 0;
	Real vtot = 0;
	Real Q = 0;
	Real Qnorm = 0;

	Size irow(0);
	Size ex(0);
	for ( auto & All_RDC_line : All_RDC_lines_ ) {

		//check for cutpoints!!!
		kinematics::FoldTree const& ft(pose.fold_tree());
		if ( All_RDC_line.res1() > pose.size() || All_RDC_line.res2() > pose.size() ) {
			tr.Warning << "non-existing residue, ignore RDC" << All_RDC_line << std::endl;
			continue;
		}
		if ( (ft.is_cutpoint(std::min((int) All_RDC_line.res1(), (int) All_RDC_line.res2())))
				&& All_RDC_line.res1() != All_RDC_line.res2() ) {
			tr.Warning << "cutpoint: ignore RDC " << All_RDC_line << std::endl;
			continue;
		}
		if ( All_RDC_line.res1() == 1 || All_RDC_line.res2() == 1  ) {
			tr.Warning << "residue 1, always problematic with 1H vs H, just ignore RDC" << std::endl;
			continue;
		}

		++irow;
		Size const d(irow - 1);

		ex = All_RDC_line.expid();

		Real computed_coupling = All_RDC_line.Jdipolar_computed_ = S_[ex][0][0] * D_[d][0] + S_[ex][0][1] * D_[d][1] + S_[ex][0][2] * D_[d][2] + S_[ex][1][1] * D_[d][3] + S_[ex][1][2] * D_[d][4];

		RDC& rdc = All_RDC_line;
		numeric::xyzVector<Real> r( pose.residue(rdc.res1()).atom(rdc.atom1()).xyz() - pose.residue(rdc.res2()).atom(rdc.atom2()).xyz());
		core::Real r2 = r.norm_squared();
		core::Real invr = 1.0 / sqrt(r2);
		core::Real invr2 = sqr(invr);
		Real obs = All_RDC_line.Jdipolar();
		Real dev = computed_coupling - obs;
		Real weight = All_RDC_line.weight()*Smax[ex]; //force constant

		//compute derivatives
		//prefactor used in derivative calculations
		core::Real pfac = weight* rdc.Dconst() * invr2 * invr;
		core::Real const pfac_NH = weight * COMMON_DENOMINATOR;

		if ( bReduced ) {
			if ( tr.Trace.visible() ) {
				tr.Trace << "reducing coupling for " << rdc << " dev: " << dev
					<< " pfac: " << pfac << " pfac_NH " << pfac_NH;
			}
			pfac = pfac_NH;
			dev *= pfac_NH / pfac;
			obs *= pfac_NH / pfac;
			if ( tr.Trace.visible() ) tr.Trace << " new dev: " << dev << std::endl;
		}

		rdc.fij_[0]= -dev * pfac * (S_[ex][0][0]*2*r.x()+S_[ex][0][1]*2*r.y()+S_[ex][0][2]*2*r.z()+S_[ex][1][1]*0+S_[ex][1][2]*0);
		rdc.fij_[1]= -dev * pfac * (S_[ex][0][0]*0+S_[ex][0][1]*2*r.x()+S_[ex][0][2]*0+S_[ex][1][1]*2*r.y()+S_[ex][1][2]*2*r.z());
		rdc.fij_[2]= -dev * pfac * (-S_[ex][0][0]*2*r.z()+S_[ex][0][1]*0+S_[ex][0][2]*2*r.x()-S_[ex][1][1]*2*r.z()+S_[ex][1][2]*2*r.y());

		// std::cout << "WEIGHT " << weight <<std::endl;
		//COMMON_DENOMINATOR is the Dcnst_NH/1.041^3
		vtot += 0.5*sqr( dev )*weight;
		wsv2 += weight*sqr(dev);
		sw += weight;
		Q += sqr( dev );
		Qnorm += sqr( obs );

	}
	//compute size for each experiment and normalize energy by size
	R_ = sqrt( Q/Qnorm/2 );
	rmsd_ = sqrt(wsv2/sw);

	//this is a factor that brings the score roughly back in the range where it was
	//with the old implementation: i.e., norm by Azz*Azz and deviation of "reduced" couplings instead of the real ones.
	//std::cout  << "DEV_TOT " << vtot <<std::endl;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;
	if ( option[ OptionKeys::rdc::print_rdc_values ].user() ) {
		std::string filename( option[ OptionKeys::rdc::print_rdc_values ]() );
		utility::io::ozstream out;

		out.open_append( filename ) ;
		using namespace core::pose::datacache;
		// mjo comment out precision because it is not used and causes a warning.
		//std::string tag( core::pose::tag_from_pose(pose) );
		for ( Size ex = 0; ex < nex_; ex++ ) {
			//Real Smax = sqrt(sqr(S_[ex][0][0]) + sqr(S_[ex][0][1]) + sqr(S_[ex][0][2]) + sqr(S_[ex][1][1]) + sqr(S_[ex][1][2]));
			//out << A( width_large, "TAG ")   << A( width, tag ) << A( width, "EXP       ")   << I( width, ex ) << I( width, Smax )
			//  << I( width, vtot )<< I( width, Rohl2Hess )<< std::endl;
			show_tensor_stats( out, ex );
			show_tensor_matrix( out, ex );
			show_rdc_values( out, ex );
		}
		out << "//" <<std::endl;
		out.close();
		//std::cout << "RDC values " << obs << ' ' << computed_coupling <<std::endl;
	}

	return vtot;
}//end of compute_dipscore

//added by LS Aug 2011
//Auxillary functions for compute_dipscore_nls
double frdc( double r0, double r1, double r2, double rdcconst, const double *par)
{
	double Ax=par[0];
	double Ay=par[1];
	double Az=-par[0]-par[1];
	//use radius instead
	double a=par[2];
	double b=par[3];
	double c=par[4];
	double rdcx = cos(b)*cos(c)*r0+(-cos(a)*sin(c)+sin(a)*sin(b)*cos(c))*r1+(sin(a)*sin(c)+cos(a)*sin(b)*cos(c))*r2;
	double rdcy = cos(b)*sin(c)*r0+(cos(a)*cos(c)+sin(a)*sin(b)*sin(c))*r1+(-sin(a)*cos(c)+cos(a)*sin(b)*sin(c))*r2;
	double rdcz = -sin(b)*r0+sin(a)*cos(b)*r1+cos(a)*cos(b)*r2;
	return rdcconst*(rdcx*rdcx*Ax+rdcy*rdcy*Ay+rdcz*rdcz*Az);
}//frdc

double frdcDa( double r0, double r1, double r2, double rdcconst, double const tensorDa, const double *par)
{
	double Ax=(3.0*par[0]/2.0-1.0)*tensorDa/(36.5089/1.041/1.041/1.041);
	double Ay=-(3.0*par[0]/2.0+1.0)*tensorDa/(36.5089/1.041/1.041/1.041);
	double Az=2.0*tensorDa/(36.5089/1.041/1.041/1.041);
	//use radius instead
	double a=par[1];
	double b=par[2];
	double c=par[3];
	double rdcx = cos(b)*cos(c)*r0+(-cos(a)*sin(c)+sin(a)*sin(b)*cos(c))*r1+(sin(a)*sin(c)+cos(a)*sin(b)*cos(c))*r2;
	double rdcy = cos(b)*sin(c)*r0+(cos(a)*cos(c)+sin(a)*sin(b)*sin(c))*r1+(-sin(a)*cos(c)+cos(a)*sin(b)*sin(c))*r2;
	double rdcz = -sin(b)*r0+sin(a)*cos(b)*r1+cos(a)*cos(b)*r2;
	return rdcconst*(rdcx*rdcx*Ax+rdcy*rdcy*Ay+rdcz*rdcz*Az);
}//frdcDa

double frdcR( double r0, double r1, double r2, double rdcconst, double const tensorR, const double *par)
{
	double Ax=(3.0*tensorR/2.0-1.0)*par[0]/(36.5089/1.041/1.041/1.041);
	double Ay=-(3.0*tensorR/2.0+1.0)*par[0]/(36.5089/1.041/1.041/1.041);
	double Az=2.0*par[0]/(36.5089/1.041/1.041/1.041);
	//use radius instead
	double a=par[1];
	double b=par[2];
	double c=par[3];
	double rdcx = cos(b)*cos(c)*r0+(-cos(a)*sin(c)+sin(a)*sin(b)*cos(c))*r1+(sin(a)*sin(c)+cos(a)*sin(b)*cos(c))*r2;
	double rdcy = cos(b)*sin(c)*r0+(cos(a)*cos(c)+sin(a)*sin(b)*sin(c))*r1+(-sin(a)*cos(c)+cos(a)*sin(b)*sin(c))*r2;
	double rdcz = -sin(b)*r0+sin(a)*cos(b)*r1+cos(a)*cos(b)*r2;
	return rdcconst*(rdcx*rdcx*Ax+rdcy*rdcy*Ay+rdcz*rdcz*Az);
}//frdcR

double frdcDaR( double r0, double r1, double r2, double rdcconst, double const tensorDa, double const tensorR, const double *par)
{
	double Ax=(3.0*tensorR/2.0-1.0)*tensorDa/(36.5089/1.041/1.041/1.041);
	double Ay=-(3.0*tensorR/2.0+1.0)*tensorDa/(36.5089/1.041/1.041/1.041);
	double Az=2.0*tensorDa/(36.5089/1.041/1.041/1.041);
	//use radius instead
	double a=par[0];
	double b=par[1];
	double c=par[2];
	double rdcx = cos(b)*cos(c)*r0+(-cos(a)*sin(c)+sin(a)*sin(b)*cos(c))*r1+(sin(a)*sin(c)+cos(a)*sin(b)*cos(c))*r2;
	double rdcy = cos(b)*sin(c)*r0+(cos(a)*cos(c)+sin(a)*sin(b)*sin(c))*r1+(-sin(a)*cos(c)+cos(a)*sin(b)*sin(c))*r2;
	double rdcz = -sin(b)*r0+sin(a)*cos(b)*r1+cos(a)*cos(b)*r2;
	return rdcconst*(rdcx*rdcx*Ax+rdcy*rdcy*Ay+rdcz*rdcz*Az);
}//frdcDaR

class data_struct {
public:

	data_struct( double* _r0, double* _r1, double* _r2, double* _rdc, double* _rdcconst, double* _rdcweight, double _tensorDa, double _tensorR, int _type_of_computation) :
		r0( _r0 ),
		r1( _r1 ),
		r2( _r2 ),
		rdc( _rdc ),
		rdcconst( _rdcconst ),
		rdcweight( _rdcweight ),
		tensorDa( _tensorDa ),
		tensorR( _tensorR ),
		type_of_computation( _type_of_computation )
	{ }

	double *r0, *r1, *r2;
	double *rdc;
	double *rdcconst;
	double *rdcweight;
	double const tensorDa;
	double const tensorR;
	int type_of_computation;

};

//Evaluaterdc function required by lmmin
void evaluaterdc(const double *par, int m_dat, const void *data, double *fvec, int * ) {
	// amw: we could make this more concise but we are being very careful to be super explicit
	// that we are not using tensorR in cases where it is its default for example
	data_struct *mydata;
	mydata = const_cast< data_struct*>(static_cast< const data_struct*>( data ) );
	int i;
	for ( i = 0; i < m_dat; i++ ) {
		fvec[i] = mydata->rdc[i];
		if ( mydata->type_of_computation == NLS ) {
			fvec[i] -= frdc( mydata->r0[i], mydata->r1[i],mydata->r2[i],mydata->rdcconst[i], par );
		} else if ( mydata->type_of_computation == NLSR ) {
			fvec[i] -= frdcR( mydata->r0[i], mydata->r1[i],mydata->r2[i],mydata->rdcconst[i], mydata->tensorR, par );
		} else if ( mydata->type_of_computation == NLSDA ) {
			fvec[i] -= frdcDa( mydata->r0[i], mydata->r1[i],mydata->r2[i],mydata->rdcconst[i], mydata->tensorDa, par );
		} else if ( mydata->type_of_computation == NLSDAR ) {
			fvec[i] -= frdcDaR( mydata->r0[i], mydata->r1[i],mydata->r2[i],mydata->rdcconst[i], mydata->tensorDa, mydata->tensorR, par );
		}
		fvec[i] *= sqrt(mydata->rdcweight[i]);
	}
}//evaluaterdc

Real ResidualDipolarCoupling::compute_dipscore_nls(core::pose::Pose const& pose) {
	// empty vectors so this legacy function is linked to the One True Function
	utility::vector1<Real> tensorDa;
	utility::vector1<Real> tensorR;
	return compute_dipscore_nls( pose, tensorDa, tensorR );
}

Real ResidualDipolarCoupling::compute_dipscore_nlsDa(core::pose::Pose const& pose, utility::vector1<Real> const & tensorDa) {
	// empty vector so this legacy function is linked to the One True Function
	utility::vector1<Real> tensorR;
	return compute_dipscore_nls( pose, tensorDa, tensorR );
}

Real ResidualDipolarCoupling::compute_dipscore_nlsR(core::pose::Pose const& pose, utility::vector1<Real> const & tensorR) {
	// empty vector so this legacy function is linked to the One True Function
	utility::vector1<Real> tensorDa;
	return compute_dipscore_nls( pose, tensorDa, tensorR );
}

Real ResidualDipolarCoupling::compute_dipscore_nlsDaR(core::pose::Pose const& pose, utility::vector1<Real> const & tensorDa, utility::vector1<Real> const & tensorR) {
	return compute_dipscore_nls( pose, tensorDa, tensorR );
}


Real ResidualDipolarCoupling::compute_dipscore_nls(
	core::pose::Pose const& pose,
	utility::vector1<Real> const & tensorDa,
	utility::vector1<Real> const & tensorR
) {

	// Really explicitly separate the four cases
	Size type_of_computation = 0;
	if ( tensorDa.size() == 0 ) {
		if ( tensorR.size() == 0 ) {
			type_of_computation = NLS;
		} else {
			type_of_computation = NLSR;
		}
	} else {
		if ( tensorR.size() == 0 ) {
			type_of_computation = NLSDA;
		} else {
			type_of_computation = NLSDAR;
		}
	}

	if ( nex_ == 0 || nrows_ == 0 ) return 0;

	//non-linear square fitting of RDC data
	utility::vector1<core::scoring::RDC>::const_iterator it;
	bool const correct_NH( basic::options::option[basic::options::OptionKeys::rdc::correct_NH_length]);
	core::Size nrow(0);
	core::Real obs(0.0);

	//initialize the cnt
	for ( core::Size id = 0; id < nex_+1; id++ ) {
		lenex_[id]=0;
	}

	core::Real scale_to_NH = COMMON_DENOMINATOR;

	for ( it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it ) {
		if ( it->res1() > pose.size() || it->res2() > pose.size() ) {
			if ( tr.Debug.visible() ) tr.Debug << "non-existing residue, ignore RDC" << std::endl;
			continue;
		}

		//check for cutpoints!!!
		kinematics::FoldTree const& ft(pose.fold_tree());
		if ( (ft.is_cutpoint(std::min((int) it->res1(), (int) it->res2()))) && it->res1() != it->res2() ) {
			if ( tr.Trace.visible() ) tr.Trace << "cutpoint: ignore RDC " << *it << std::endl;
			continue;
		}

		++nrow;

		numeric::xyzVector<Real> r( pose.residue(it->res1()).atom(it->atom1()).xyz() - pose.residue(it->res2()).atom(it->atom2()).xyz());
		core::Real r2 = r.norm_squared();
		core::Real invr = 1.0 / sqrt(r2);
		if ( correct_NH )  {
			do_correct_NH( it, r, r2, invr );
		}

		core::Real pfac = it->Dconst();
		bool bCSA(false);// hook up for later... to compute chemical shift anisotropy
		if ( !bCSA ) {
			pfac *= invr * invr * invr;
		}

		//check the -1 if it is correct
		core::Size id = All_RDC_lines_[nrow].expid();
		obs = All_RDC_lines_[nrow].Jdipolar()*(scale_to_NH)/pfac;

		r0_[nrow-1] = r.normalized().x();
		r1_[nrow-1] = r.normalized().y();
		r2_[nrow-1] = r.normalized().z();
		rdcconst_[nrow-1]=scale_to_NH;
		exprdc_[nrow-1]=obs;
		rdcweight_[nrow-1]=it->weight();
		lenex_[id+1]=lenex_[id+1]+1;
		//tr.Trace << "weight: " << it->weight() << std::endl;
	} //cycle over atoms


	//parameters
	int n_par = 3; // number of parameters in model function frdcR
	int nrepeat = basic::options::option[ basic::options::OptionKeys::rdc::nlsrepeat ](); // number of repeat lmfit

	std::vector<double> parbest(n_par*nex_);
	std::vector<double> par(n_par*nex_);

	int i,j;
	Size prelen=0;

	//optional weighting provided by user
	utility::vector1<Real> Smax(nex_);

	//experimental data array
	for ( Size ex = 0; ex < nex_; ex++ ) {

		//compute the length of previous exps
		prelen=0;
		for ( Size cnt=0; cnt<=ex; cnt++ ) {
			prelen+=lenex_[cnt];
		}

		//perform lmfit on each exp
		double p1 = -99999;
		double p2 = -99999;
		if ( type_of_computation == NLSR ) {
			p2 = tensorR[ex+1];
		} else if ( type_of_computation == NLSDA ) {
			p1 = tensorDa[ex+1];
		} else if ( type_of_computation == NLSDAR ) {
			p1 = tensorDa[ex+1];
			p2 = tensorR[ex+1];
		}

		data_struct data( &r0_[prelen], &r1_[prelen], &r2_[prelen], &exprdc_[prelen], &rdcconst_[prelen], &rdcweight_[prelen], p1, p2, type_of_computation );

		//definition of auxiliary parameters
		numeric::nls::lm_status_struct status;

		Smax[ex] = 1;
		{
			using namespace basic::options;
			if ( option[ OptionKeys::rdc::fix_normAzz ].user() ) {
				if ( option[ OptionKeys::rdc::fix_normAzz ]().size() != nex_ ) {
					utility_exit_with_message("fix_normAzz must have one value for each alignment medium !");
				}
				Smax[ ex ] = option[ OptionKeys::rdc::fix_normAzz ]()[ ex+1 ];
			}
		}

		double bestnorm = 1e12;

		for ( j = 0; j < nrepeat; j++ ) {
			//random starting value
			par[ex*n_par+0]=2.0*numeric::NumericTraits<Real>::pi()*numeric::random::uniform();//alpha
			par[ex*n_par+1]=2.0*numeric::NumericTraits<Real>::pi()*numeric::random::uniform();//beta
			par[ex*n_par+2]=2.0*numeric::NumericTraits<Real>::pi()*numeric::random::uniform();//gamma
			//call lmmin
			numeric::nls::lmmin( n_par, &par[ex*n_par], lenex_[ex+1], (const void*) &data, evaluaterdc, &status,numeric::nls::lm_printout_std);
			if ( tr.Trace.visible() ) {
				tr.Trace << std::endl;
				tr.Trace << "Iteration: " << j << " status.fnorm: " << status.fnorm << " bestnorm: "<< bestnorm << std::endl;
			}
			//save to best fitting parameter
			if ( status.fnorm < bestnorm ) {
				tr.Trace << "status.fnorm: " << status.fnorm << " replaced by bestnorm: " << bestnorm << std::endl;
				bestnorm=status.fnorm;
				for ( i=0; i<n_par ; i++ ) {
					parbest[ex*n_par+i]=par[ex*n_par+i];
				}
			}

		}//repeat lmfit five times with different starting parameters

		//copy back to par
		for ( i=0; i<n_par ; i++ ) {
			par[ex*n_par+i]=parbest[ex*n_par+i];
		}

		double Ax = 0;
		double Ay = 0;
		//sort the right order Ax<Ay and calculate Da and R
		if ( type_of_computation == NLS ) {
			if ( parbest[ex*n_par+0]>0 ) {
				if ( parbest[ex*n_par+1]>0 ) {
					if ( parbest[ex*n_par+0]<parbest[ex*n_par+1] ) {
						Ax=parbest[ex*n_par+0];
						Ay=parbest[ex*n_par+1];
					} else {
						Ax=parbest[ex*n_par+1];
						Ay=parbest[ex*n_par+0];
					}
				} else {
					if ( parbest[ex*n_par+0]<-parbest[ex*n_par+1] ) {
						Ax=std::min(parbest[ex*n_par+0],-parbest[ex*n_par+1]-parbest[ex*n_par+0]);
						Ay=std::max(parbest[ex*n_par+0],-parbest[ex*n_par+1]-parbest[ex*n_par+0]);
					} else {
						Ax=std::max(parbest[ex*n_par+1],-parbest[ex*n_par+1]-parbest[ex*n_par+0]);
						Ay=std::min(parbest[ex*n_par+1],-parbest[ex*n_par+1]-parbest[ex*n_par+0]);
					}
				}
			} else {
				if ( parbest[ex*n_par+1]>0 ) {
					if ( -parbest[ex*n_par+0]<parbest[ex*n_par+1] ) {
						Ax=std::max(parbest[ex*n_par+0],-parbest[ex*n_par+1]-parbest[ex*n_par+0]);
						Ay=std::min(parbest[ex*n_par+0],-parbest[ex*n_par+1]-parbest[ex*n_par+0]);
					} else {
						Ax=std::min(parbest[ex*n_par+1],-parbest[ex*n_par+1]-parbest[ex*n_par+0]);
						Ay=std::max(parbest[ex*n_par+1],-parbest[ex*n_par+1]-parbest[ex*n_par+0]);
					}
				} else {
					if ( parbest[ex*n_par+0]<parbest[ex*n_par+1] ) {
						Ax=parbest[ex*n_par+1];
						Ay=parbest[ex*n_par+0];
					} else {
						Ax=parbest[ex*n_par+0];
						Ay=parbest[ex*n_par+1];
					}
				}
			}
			//store in the temporarily in parbest
			parbest[ex*n_par+0]=1.0/2.0*(-Ax-Ay);
			parbest[ex*n_par+1]=2.0/3.0*(Ay-Ax)/(Ax+Ay);
		}

		//debug
		tr.Trace << std::endl;
		tr.Trace << "ex: " << ex << std::endl;
		if ( type_of_computation == NLSDA || type_of_computation == NLSDAR ) {
			tr.Trace << " tensorDa["<<ex<<"]: "<<tensorDa[ex+1]<< std::endl;
		}
		if ( type_of_computation == NLSR || type_of_computation == NLSDAR ) {
			tr.Trace << " tensorR["<<ex<<"]: "<<tensorR[ex+1]<< std::endl;
		}
		tr.Trace << "Ax: " << (3.0*tensorR[ex+1]/2.0-1.0)*tensorDa[ex+1]/COMMON_DENOMINATOR<< std::endl;
		tr.Trace << "Ay: " << -(3.0*tensorR[ex+1]/2.0+1.0)*tensorDa[ex+1]/COMMON_DENOMINATOR<< std::endl;
		tr.Trace << "alpha: " << par[ex*n_par+0]<< std::endl;
		tr.Trace << "beta: " << par[ex*n_par+1]<< std::endl;
		tr.Trace << "gamma:" << par[ex*n_par+2]<< std::endl;
		if ( type_of_computation == NLS ) {
			tr.Trace << "Da:" << 1.0/2.0*(-Ax-Ay)*COMMON_DENOMINATOR<< std::endl;
			tr.Trace << "R:" << 2.0/3.0*(Ay-Ax)/(Ax+Ay)<< std::endl;
		}
		tr.Trace << "norm:" << bestnorm<<std::endl;
		if ( type_of_computation == NLS ) {
			tr.Trace << "Pales Da:" << 3.0/4.0*(-Ax-Ay)<< std::endl;
			tr.Trace << "Pales Dr:" << 1.0/2.0*(Ax-Ay)<< std::endl;
		}
	}//end of loop through all exps

	//compute scores and other stats for all exps
	Real wsv2 = 0;
	Real sw = 0;
	Real vtot = 0;
	Real Q = 0;
	Real Qex = 0;
	Real Qnorm = 0;

	Size irow(0);
	for ( auto it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it ) {

		Size ex = it->expid();

		//compute the length of previous exps
		prelen = 0;
		for ( Size cnt = 0; cnt <= ex; cnt++ ) {
			prelen += lenex_[ cnt ];
		}

		// amw: another place where I need to split behavior
		Real computed_coupling;
		if ( type_of_computation == NLS ) {
			computed_coupling = frdc(r0_[prelen+irow], r1_[prelen+irow], r2_[prelen+irow], rdcconst_[prelen+irow], &par[ex*n_par]);
		} else if ( type_of_computation == NLSDA ) {
			computed_coupling = frdcDa(r0_[prelen+irow], r1_[prelen+irow], r2_[prelen+irow], rdcconst_[prelen+irow], tensorDa[ex+1], &par[ex*n_par]);
		} else if ( type_of_computation == NLSR ) {
			computed_coupling = frdcR(r0_[prelen+irow], r1_[prelen+irow], r2_[prelen+irow], rdcconst_[prelen+irow], tensorR[ex+1], &par[ex*n_par]);
		} else {
			computed_coupling = frdcDaR(r0_[prelen+irow], r1_[prelen+irow], r2_[prelen+irow], rdcconst_[prelen+irow], tensorDa[ex+1], tensorR[ex+1], &par[ex*n_par]);
		}

		// std::cout << "WEIGHT " << weight <<std::endl;
		//normalized by the tensor values and number of experiments

		//compute derivatives
		RDC& rdc = *it;
		numeric::xyzVector<Real> r( pose.residue(rdc.res1()).atom(rdc.atom1()).xyz() - pose.residue(rdc.res2()).atom(rdc.atom2()).xyz());
		numeric::xyzVector<Real> r_copy( pose.residue(rdc.res1()).atom(rdc.atom1()).xyz() - pose.residue(rdc.res2()).atom(rdc.atom2()).xyz());
		core::Real r2 = r.norm_squared();
		core::Real invr = 1.0 / sqrt(r2);

		if ( correct_NH )  {
			do_correct_NH( it, r, r2, invr );
		}
		core::Real invr2 = sqr(invr);

		//prefactor used in derivative calculations
		Real weight = it->weight()*Smax[ex]; //force constant
		core::Real pfac = scale_to_NH*(1/r.length())*(1/r_copy.length())*weight;
		Real obs = rdc.Jdipolar()*(scale_to_NH)/(rdc.Dconst() * invr2 * invr);
		Real dev = computed_coupling - obs;
		it->Jdipolar_computed_ = computed_coupling/((scale_to_NH)/(rdc.Dconst() * invr2 * invr));

		core::Real Axx, Ayy, Azz, a, b, c;
		//parameters after fitting
		if ( type_of_computation == NLS ) {
			Axx= par[ex*n_par+0];
			Ayy= par[ex*n_par+1];
			Azz= -par[ex*n_par+0]-par[ex*n_par+1];
			a=par[ex*n_par+2];
			b=par[ex*n_par+3];
			c=par[ex*n_par+4];
		} else if ( type_of_computation == NLSR ) {
			Axx=(3.0*tensorR[ex+1]/2.0-1.0)*par[ex*n_par+0]/COMMON_DENOMINATOR;
			Ayy=-(3.0*tensorR[ex+1]/2.0+1.0)*par[ex*n_par+0]/COMMON_DENOMINATOR;
			Azz=2.0*par[ex*n_par+0]/COMMON_DENOMINATOR;
			a=par[ex*n_par+1];
			b=par[ex*n_par+2];
			c=par[ex*n_par+3];

		} else if ( type_of_computation == NLSDA ) {
			Axx=(3.0*par[ex*n_par+0]/2.0-1.0)*tensorDa[ex+1]/ COMMON_DENOMINATOR;
			Ayy=-(3.0*par[ex*n_par+0]/2.0+1.0)*tensorDa[ex+1]/ COMMON_DENOMINATOR;
			Azz=2.0*tensorDa[ex+1]/ COMMON_DENOMINATOR;
			a=par[ex*n_par+1];
			b=par[ex*n_par+2];
			c=par[ex*n_par+3];
		} else {
			Axx=(3.0*tensorR[ex+1]/2.0-1.0)*tensorDa[ex+1]/ COMMON_DENOMINATOR;
			Ayy=-(3.0*tensorR[ex+1]/2.0+1.0)*tensorDa[ex+1]/ COMMON_DENOMINATOR;
			Azz=2.0*tensorDa[ex+1]/ COMMON_DENOMINATOR;
			a=par[ex*n_par+0];
			b=par[ex*n_par+1];
			c=par[ex*n_par+2];
		}
		//rotation matrix
		core::Real r00=cos(b)*cos(c);
		core::Real r01=-cos(a)*sin(c)+sin(a)*sin(b)*cos(c);
		core::Real r02=sin(a)*sin(c)+cos(a)*sin(b)*cos(c);
		core::Real r10=cos(b)*sin(c);
		core::Real r11=cos(a)*cos(c)+sin(a)*sin(b)*sin(c);
		core::Real r12=-sin(a)*cos(c)+cos(a)*sin(b)*sin(c);
		core::Real r20=-sin(b);
		core::Real r21=sin(a)*cos(b);
		core::Real r22=cos(a)*cos(b);
		//projects in the PAF
		core::Real rx=r00*r.x()+r01*r.y()+r02*r.z();
		core::Real ry=r10*r.x()+r11*r.y()+r12*r.z();
		core::Real rz=r20*r.x()+r21*r.y()+r22*r.z();
		//derivatives
		rdc.fij_[0] = - dev * pfac * 2 * (Axx*rx*r00+Ayy*ry*r10+Azz*rz*r20) ;
		rdc.fij_[1] = - dev * pfac * 2 * (Axx*rx*r01+Ayy*ry*r11+Azz*rz*r21);
		rdc.fij_[2] = - dev * pfac * 2 * (Axx*rx*r02+Ayy*ry*r12+Azz*rz*r22) ;

		//compute energy
		vtot += 0.5*sqr( dev )*weight;
		wsv2 += weight*sqr(dev);
		sw += weight;
		Q += sqr( dev );
		Qex +=sqr( dev );
		Qnorm += sqr( obs );

		//printout Qbax_ for each experiment
		if ( irow==lenex_[ex+1]-1 ) {
			if ( tr.Trace.visible() ) {
				Real qbax;
				if ( type_of_computation == NLS ) {
					qbax = sqrt(Qex/lenex_[ex+1])/sqrt(sqr(parbest[ex*n_par+0]*COMMON_DENOMINATOR)*(4+3*sqr(parbest[ex*n_par+1]))/5);
				} else if ( type_of_computation == NLSR ) {
					qbax = sqrt(Qex/lenex_[ex+1])/sqrt(sqr(par[ex*n_par+0])*(4+3*sqr(tensorR[ex+1]))/5);
				} else if ( type_of_computation == NLSDA ) {
					qbax = sqrt(Qex/lenex_[ex+1])/sqrt(sqr(tensorDa[ex+1])*(4+3*sqr(par[ex*n_par+0]))/5);
				} else {
					qbax = sqrt(Qex/lenex_[ex+1])/sqrt(sqr(tensorDa[ex+1])*(4+3*sqr(tensorR[ex+1]))/5);
				}
				tr.Trace << "ex: " << ex << " Qbax_: " << qbax << std::endl; //JACS 2003 125(30) 9179-9191 Table 2 lagend
				tr.Trace << "ex: " << ex << " rms_dc: " << sqrt(Qex/lenex_[ex+1]) << std::endl;
			}
			Qex=0;
		}

		//increament the array
		if ( irow < lenex_[ ex + 1 ] - 1 ) {
			irow++;
		} else {
			irow=0;
		}
	}

	R_ = sqrt( Q/Qnorm/2 );
	rmsd_ = sqrt(wsv2/sw);

	tr.Trace << "R_: " << R_ << std::endl;
	tr.Trace << "rmsd_: " << rmsd_ << std::endl;
	tr.Trace << "All_RDC_lines_.size(): " << All_RDC_lines_.size() << std::endl;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;

	//print the rdc values
	if ( option[ OptionKeys::rdc::print_rdc_values ].user() ) {
		std::string filename( option[ OptionKeys::rdc::print_rdc_values ]() );
		utility::io::ozstream out;

		out.open_append( filename ) ;
		using namespace core::pose::datacache;
		for ( Size ex = 0; ex < nex_; ex++ ) {
			show_rdc_values( out, ex );
		}
		out << "//" <<std::endl;
		out.close();
	}

	return vtot;
}

void
ResidualDipolarCoupling::do_correct_NH(
	utility::vector1<core::scoring::RDC>::const_iterator it,
	numeric::xyzVector<Real> & r,
	core::Real & r2,
	core::Real & invr
) {
	if ( it->type() == RDC::RDC_TYPE_NH && std::abs((int) it->res1()-(int) it->res2())==0 ) {
		r.normalize(1.041);
		r2 = r.norm_squared();
		invr = 1.0 / sqrt(r2);
	} else if ( it->type() == RDC::RDC_TYPE_NC && std::abs((int) it->res1()-(int) it->res2())==1 ) {
		r.normalize(1.329);
		r2 = r.norm_squared();
		invr = 1.0 / sqrt(r2);
	} else if ( it->type() == RDC::RDC_TYPE_CH && std::abs((int) it->res1()-(int) it->res2())==0 ) {
		r.normalize(1.107);
		r2 = r.norm_squared();
		invr = 1.0 / sqrt(r2);
	} else if ( it->type() == RDC::RDC_TYPE_CC && std::abs((int) it->res1()-(int) it->res2())==0 ) {
		r.normalize(1.525);
		r2 = r.norm_squared();
		invr = 1.0 / sqrt(r2);
	} else if ( it->type() == RDC::RDC_TYPE_CHN && std::abs((int) it->res1()-(int) it->res2())==1 ) {
		r.normalize(2.085);
		r2 = r.norm_squared();
		invr = 1.0 / sqrt(r2);
	} else if ( it->type() == RDC::RDC_TYPE_NCA && std::abs((int) it->res1()-(int) it->res2())==0 ) {
		r.normalize(1.458);
		r2 = r.norm_squared();
		invr = 1.0 / sqrt(r2);
	} else {
		tr.Error << "unreognized type or residue sequence separation does not allow using correct_NH" << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput, "unreognized type or residue sequence separation does not allow using correct_NH ");
	}
}

void ResidualDipolarCoupling::show_tensor_stats( std::ostream& out, Size ex ) const {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;

	Real Aa = EV_[ex][0]/2;
	Real Ar = (EV_[ex][2]-EV_[ex][1])/3;

	Size const width( 10 );
	Size const precision( 2 );
	Real rhombicity = Ar/Aa;
	out << A( width, "Ev[ex][0]" ) << " " << A( width, "Ev[ex][1]" ) << A( width, "Ev[ex][2]" ) << std::endl;
	out << F( width, precision, EV_[ex][0] ) << F( width, precision, EV_[ex][1] ) << F( width, precision, EV_[ex][2]) << std::endl;
	out << A( width, "Aa" ) << " " << A( width, "Ar" ) << A( width, "rhombicity" ) << std::endl;
	out << F( width, precision, Aa ) << F( width, precision, Ar ) << F( width, precision, rhombicity ) << std::endl;
}

void ResidualDipolarCoupling::show_tensor_stats_nls( std::ostream& out, Size, const double *par) const {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;

	Real Ax = par[0];
	Real Ay = par[1];
	Real Az = -par[0]-par[1];
	Size const width( 10 );
	Size const precision( 2 );
	out << A( width, "Ax" ) << " " << A( width, "Ay" ) << A( width, "Az" ) << std::endl;
	out << F( width, precision, Ax ) << F( width, precision, Ay ) << F( width, precision, Az) << std::endl;
}

void ResidualDipolarCoupling::show_tensor_matrix( std::ostream& out, Size ex ) const {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;

	Size const width( 8 );
	//mjo comment out width_large because it is unused and causes a warning
	Size const precision( 2 );
	out << A( width, "AL.EVAL   ")   << F(width, precision, EV_[ex][0]) <<  F(width, precision, EV_[ex][1])  <<  F(width, precision, EV_[ex][2])  << std::endl
		<< A( width, "AL.EVEC XX")   << F(width, precision, EIG_[ex][0][0]) <<  F(width, precision, EIG_[ex][1][0])  <<  F(width, precision, EIG_[ex][2][0])  << std::endl
		<< A( width, "AL.EVEC YY")   << F(width, precision, EIG_[ex][0][1]) <<  F(width, precision, EIG_[ex][1][1])  <<  F(width, precision, EIG_[ex][2][1])  << std::endl
		<< A( width, "AL.EVEC ZZ")   << F(width, precision, EIG_[ex][0][2]) <<  F(width, precision, EIG_[ex][1][2])  <<  F(width, precision, EIG_[ex][2][2])  << std::endl;
}

void ResidualDipolarCoupling::show_rdc_values( std::ostream& out, Size ex ) const {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;

	Size const width( 8 );
	//mjo comment out width_large because it is unused and causes a warning
	Size const precision( 3 );
	utility::vector1<core::scoring::RDC>::const_iterator it;
	for ( it = All_RDC_lines_.begin(); it != All_RDC_lines_.end(); ++it ) {
		if ( it->expid() != ex ) continue;

		Size count( it->res1() );
		std::string type;
		if ( it->type() == RDC::RDC_TYPE_NH ) { type ="NH";}
		if ( it->type() == RDC::RDC_TYPE_CH ) { type ="CH";}
		if ( it->type() == RDC::RDC_TYPE_NC ) { type ="NC";}
		if ( it->type() == RDC::RDC_TYPE_CC ) { type ="CC";}
		if ( it->type() == RDC::RDC_TYPE_CHN ) { type ="CHN";}
		if ( it->type() == RDC::RDC_TYPE_NCA ) { type ="NCA";}
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

void ResidualDipolarCoupling::compute_tensor_stats() {
	for ( Size ex = 0; ex < nex_; ex++ ) {
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

		// std::cout << "AL.TENSOR (molecular frame): " << SD_[ex][0][0] << ' ' <<  SD_[ex][0][1]  << ' ' <<  SD_[ex][0][2] << ' ' << SD_[ex][1][1] << ' ' << SD_[ex][1][2] << std::endl;

		Tensor v2;
		int nrot2;
		rvec ev2;

		jacobi3(SD_[ex],ev2,v2,nrot2);

		EIG_[ex][0][0] =  v2[0][0];
		EIG_[ex][0][1] =  v2[0][1];
		EIG_[ex][0][2] =  v2[0][2];
		EIG_[ex][1][1] =  v2[1][1];
		EIG_[ex][1][2] =  v2[1][2];
		EIG_[ex][1][0] =  v2[1][0];
		EIG_[ex][2][0] =  v2[2][0];
		EIG_[ex][2][1] =  v2[2][1];
		EIG_[ex][2][2] =  v2[2][2];

		qsort(&ev2[0], 3 ,sizeof(core::Real),compare_by_abs);

		EV_[ex][0] = ev2[2];
		EV_[ex][1] = ev2[1];
		EV_[ex][2] = ev2[0];

		//  std::cout << "AL_TENSOR eigenvalues : " << EV_[ex][0] << ' ' << EV_[ex][1]  << ' ' <<  EV_[ex][2]  << ' ' << nrot2 << std::endl;

		/*std::cout << "AL_TENSOR eigenvectors  : " << std::endl;
		std::cout << v2[0][0] << "\t" << v2[1][0]  << "\t" <<  v2[2][0]  << std::endl;
		std::cout << v2[0][1] << "\t" << v2[1][1]  << "\t" <<  v2[2][1]  << std::endl;
		std::cout << v2[0][2] << "\t" << v2[1][2]  << "\t" <<  v2[2][2]  << std::endl;*/
		trace_[ex] = v2[0][0] + v2[1][1] + v2[2][2];
		rvec tempvec;
		tempvec[0] =  v2[0][2];
		tempvec[1] =  v2[1][2];
		tempvec[2] =  v2[2][2];
		qsort(&tempvec[0], 3, sizeof(core::Real),compare_by_abs);
		maxz_[ex] = tempvec[0];
		// std::cout << "Fractional Anisotropy " << FA_[ex] << std::endl;

		// std::cout << "AL.TENSOR (molecular frame after): " << S_[ex][0][0] << ' ' <<  S_[ex][0][1]  << ' ' <<  S_[ex][0][2] << ' ' << S_[ex][1][1] << ' ' << S_[ex][1][2] << std::endl;

	} //for ( ex = 0 .. nex )
}

int m_inv_gen( ResidualDipolarCoupling::Tensor5 const & m, int n, ResidualDipolarCoupling::Tensor5 & minv )
{
	ResidualDipolarCoupling::Tensor5 md,v;
	ResidualDipolarCoupling::rvec5 eig;
	Real tol,s;
	int nzero,i,j,k,nrot;
	for ( i=0; i<n; i++ ) {
		for ( j=0; j<n; j++ ) {
			md[i][j] = m[i][j];
		}
	}

	tol = 0;
	for ( i=0; i<n; i++ ) {
		tol += fabs(md[i][i]);
	}
	tol = 1e-6*tol/n;

	jacobi(md,eig,v,nrot);

	nzero = 0;
	for ( i=0; i<n; i++ ) {
		if ( fabs(eig[i]) < tol ) {
			eig[i] = 0;
			nzero++;
		} else {
			eig[i] = 1.0/eig[i];
		}
	}

	for ( i=0; i<n; i++ ) {
		for ( j=0; j<n; j++ ) {
			s = 0;
			for ( k=0; k<n; k++ ) {
				s += eig[k]*v[i][k]*v[j][k];
			}
			minv[i][j] = s;
		}
	}

	return nzero;
}

////// ITS REALLY ANNOYING THAT YOU ARE USING THIS MACRO.
////// MACRO USE IS A VIOLATION OF THE CODING CONVENTIONS
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
  a[k][l]=h+s*(g-h*tau)

void jacobi( ResidualDipolarCoupling::Tensor5 & a, ResidualDipolarCoupling::rvec5 & d, ResidualDipolarCoupling::Tensor5 & v, int & nrot )
{
	int j,i;
	int iq,ip;
	Real tresh,theta,tau,t,sm,s,h,g,c;
	Real b[5];
	Real z[5];
	int const n( 5 );
	for ( ip=0; ip<n; ip++ ) {
		for ( iq=0; iq<n; iq++ ) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for ( ip=0; ip<n; ip++ ) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for ( i=1; i<=50; i++ ) {
		sm=0.0;
		for ( ip=0; ip<n-1; ip++ ) {
			for ( iq=ip+1; iq<n; iq++ ) {
				sm += fabs(a[ip][iq]);
			}
		}
		if ( sm == 0.0 ) {
			return;
		}
		if ( i < 4 ) { //first 3 iterations
			tresh=0.2*sm/(n*n);
		} else {
			tresh=0.0;
		}
		for ( ip=0; ip<n-1; ip++ ) {
			for ( iq=ip+1; iq<n; iq++ ) {
				g=100.0*fabs(a[ip][iq]);
				if ( i > 4 && fabs(d[ip])+g == fabs(d[ip])
						&& fabs(d[iq])+g == fabs(d[iq]) ) {
					a[ip][iq]=0.0;
				} else if ( fabs(a[ip][iq]) > tresh ) {
					h=d[iq]-d[ip];
					if ( fabs(h)+g == fabs(h) ) {
						t=(a[ip][iq])/h;
					} else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if ( theta < 0.0 ) t = -t;
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
					for ( j=0; j<ip; j++ ) {
						ROTATE(a,j,ip,j,iq);
					}
					for ( j=ip+1; j<iq; j++ ) {
						ROTATE(a,ip,j,j,iq);
					}
					for ( j=iq+1; j<n; j++ ) {
						ROTATE(a,ip,j,iq,j);
					}
					for ( j=0; j<n; j++ ) {
						ROTATE(v,j,ip,j,iq);
					}
					++nrot;
				}
			}
		}
		for ( ip=0; ip<n; ip++ ) {
			b[ip] +=  z[ip];
			d[ip]  =  b[ip];
			z[ip]  =  0.0;
		}
	}
	//probably different type of Exception is better suited
	throw CREATE_EXCEPTION(utility::excn::BadInput, " too many iterations in Jacobi when compute RDC tensor");
}

void jacobi3(
	ResidualDipolarCoupling::Tensor & a,
	ResidualDipolarCoupling::rvec & d,
	ResidualDipolarCoupling::Tensor & v,
	int & nrot
) {
	int j,i;
	int iq,ip;
	Real tresh,theta,tau,t,sm,s,h,g,c;
	Real b[3];
	Real z[3];
	int const n( 3 );
	for ( ip=0; ip<n; ip++ ) {
		for ( iq=0; iq<n; iq++ ) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for ( ip=0; ip<n; ip++ ) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for ( i=1; i<=50; i++ ) {
		sm=0.0;
		for ( ip=0; ip<n-1; ip++ ) {
			for ( iq=ip+1; iq<n; iq++ ) {
				sm += fabs(a[ip][iq]);
			}
		}
		if ( sm == 0.0 ) {
			return;
		}
		if ( i < 4 ) { //first 3 iterations
			tresh=0.2*sm/(n*n);
		} else {
			tresh=0.0;
		}
		for ( ip=0; ip<n-1; ip++ ) {
			for ( iq=ip+1; iq<n; iq++ ) {
				g=100.0*fabs(a[ip][iq]);
				if ( i > 4 && fabs(d[ip])+g == fabs(d[ip])
						&& fabs(d[iq])+g == fabs(d[iq]) ) {
					a[ip][iq]=0.0;
				} else if ( fabs(a[ip][iq]) > tresh ) {
					h=d[iq]-d[ip];
					if ( fabs(h)+g == fabs(h) ) {
						t=(a[ip][iq])/h;
					} else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if ( theta < 0.0 ) t = -t;
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
					for ( j=0; j<ip; j++ ) {
						ROTATE(a,j,ip,j,iq);
					}
					for ( j=ip+1; j<iq; j++ ) {
						ROTATE(a,ip,j,j,iq);
					}
					for ( j=iq+1; j<n; j++ ) {
						ROTATE(a,ip,j,iq,j);
					}
					for ( j=0; j<n; j++ ) {
						ROTATE(v,j,ip,j,iq);
					}
					++nrot;
				}
			}
		}
		for ( ip=0; ip<n; ip++ ) {
			b[ip] +=  z[ip];
			d[ip]  =  b[ip];
			z[ip]  =  0.0;
		}
	}
	//probably different type of Exception is better suited
	throw CREATE_EXCEPTION(utility::excn::BadInput, " too many iterations in Jacobi when compute RDC tensor");
}

} //namespace Scoring
} //namespace core



#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::RDC::save( Archive & arc ) const {
	arc( CEREAL_NVP( type_ ) ); // enum core::scoring::RDC::RDC_TYPE
	arc( CEREAL_NVP( res1_ ) ); // Size
	arc( CEREAL_NVP( res2_ ) ); // Size
	arc( CEREAL_NVP( Jdipolar_ ) ); // Real
	arc( CEREAL_NVP( Reduced_Jdipolar_ ) ); // Real
	arc( CEREAL_NVP( weight_ ) ); // Real
	arc( CEREAL_NVP( Jdipolar_computed_ ) ); // Real
	arc( CEREAL_NVP( fij_ ) ); // core::Vector
	arc( CEREAL_NVP( expid_ ) ); // Size
	arc( CEREAL_NVP( atom1_ ) ); // std::string
	arc( CEREAL_NVP( atom2_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::RDC::load( Archive & arc ) {
	arc( type_ ); // enum core::scoring::RDC::RDC_TYPE
	arc( res1_ ); // Size
	arc( res2_ ); // Size
	arc( Jdipolar_ ); // Real
	arc( Reduced_Jdipolar_ ); // Real
	arc( weight_ ); // Real
	arc( Jdipolar_computed_ ); // Real
	arc( fij_ ); // core::Vector
	arc( expid_ ); // Size
	arc( atom1_ ); // std::string
	arc( atom2_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::RDC );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::ResidualDipolarCoupling::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( All_RDC_lines_ ) ); // RDC_lines
	arc( CEREAL_NVP( EV_ ) ); // utility::vector0<rvec>
	arc( CEREAL_NVP( D_ ) ); // utility::vector0<rvec5>
	arc( CEREAL_NVP( rhs_ ) ); // utility::vector0<rvec5>
	arc( CEREAL_NVP( S_ ) ); // utility::vector0<Tensor>
	arc( CEREAL_NVP( T_ ) ); // utility::vector0<Tensor5>
	arc( CEREAL_NVP( nex_ ) ); // core::Size
	arc( CEREAL_NVP( nrows_ ) ); // core::Size
	arc( CEREAL_NVP( R_ ) ); // core::Real
	arc( CEREAL_NVP( rmsd_ ) ); // core::Real
	arc( CEREAL_NVP( SD_ ) ); // utility::vector0<Tensor>
	arc( CEREAL_NVP( EIG_ ) ); // utility::vector0<Tensor>
	arc( CEREAL_NVP( FA_ ) ); // utility::vector0<core::Real>
	arc( CEREAL_NVP( trace_ ) ); // utility::vector0<core::Real>
	arc( CEREAL_NVP( maxz_ ) ); // utility::vector0<core::Real>
	arc( CEREAL_NVP( r0_ ) ); // utility::vector0<core::Real>
	arc( CEREAL_NVP( r1_ ) ); // utility::vector0<core::Real>
	arc( CEREAL_NVP( r2_ ) ); // utility::vector0<core::Real>
	arc( CEREAL_NVP( exprdc_ ) ); // utility::vector0<core::Real>
	arc( CEREAL_NVP( rdcconst_ ) ); // utility::vector0<core::Real>
	arc( CEREAL_NVP( rdcweight_ ) ); // utility::vector0<core::Real>
	arc( CEREAL_NVP( lenex_ ) ); // utility::vector0<core::Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::ResidualDipolarCoupling::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( All_RDC_lines_ ); // RDC_lines
	arc( EV_ ); // utility::vector0<rvec>
	arc( D_ ); // utility::vector0<rvec5>
	arc( rhs_ ); // utility::vector0<rvec5>
	arc( S_ ); // utility::vector0<Tensor>
	arc( T_ ); // utility::vector0<Tensor5>
	arc( nex_ ); // core::Size
	arc( nrows_ ); // core::Size
	arc( R_ ); // core::Real
	arc( rmsd_ ); // core::Real
	arc( SD_ ); // utility::vector0<Tensor>
	arc( EIG_ ); // utility::vector0<Tensor>
	arc( FA_ ); // utility::vector0<core::Real>
	arc( trace_ ); // utility::vector0<core::Real>
	arc( maxz_ ); // utility::vector0<core::Real>
	arc( r0_ ); // utility::vector0<core::Real>
	arc( r1_ ); // utility::vector0<core::Real>
	arc( r2_ ); // utility::vector0<core::Real>
	arc( exprdc_ ); // utility::vector0<core::Real>
	arc( rdcconst_ ); // utility::vector0<core::Real>
	arc( rdcweight_ ); // utility::vector0<core::Real>
	arc( lenex_ ); // utility::vector0<core::Size>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::ResidualDipolarCoupling );
CEREAL_REGISTER_TYPE( core::scoring::ResidualDipolarCoupling )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_ResidualDipolarCoupling )
#endif // SERIALIZATION
