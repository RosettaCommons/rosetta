// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @Author Frank Dimaio


#include <protocols/cryst/spacegroup.hh>
#include <protocols/cryst/refinable_lattice.hh>
#include <apps/pilot/frank/cryst_reporters.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/Remarks.hh>
#include <core/io/CrystInfo.hh>
#include <core/pose/motif/reference_frames.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/types.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/motif/util.hh>
#include <core/scoring/motif/xfrags.hh>
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/symmetry/SymmetricEnergies.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>


#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/cryst/refinable_lattice.hh>
#include <protocols/cryst/cryst_rms.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/flxbb/LayerDesignOperation.hh>
#include <protocols/protein_interface_design/design_utils.hh>
#include <protocols/cryst/util.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/constants.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/UniformRotationSampler.hh>
#include <numeric/random/random_permutation.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

// option includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <fstream>
#include <iostream>
#include <math.h>

#include <sstream>
#include <string>
#include <queue>


using namespace basic;
using namespace core;
using namespace core::pose;
using namespace core::io::pdb;
using namespace core::io;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::cryst;

static THREAD_LOCAL basic::Tracer TR("cryst.design");

OPT_1GRP_KEY(String, crystdock, mode)
OPT_1GRP_KEY(String, crystdock, spacegroup)
OPT_1GRP_KEY(Real, crystdock, A)
OPT_1GRP_KEY(Real, crystdock, B)
OPT_1GRP_KEY(Real, crystdock, C)
OPT_1GRP_KEY(Real, crystdock, alpha)
OPT_1GRP_KEY(Real, crystdock, beta)
OPT_1GRP_KEY(Real, crystdock, gamma)
OPT_1GRP_KEY(Real, crystdock, maxclash)
OPT_1GRP_KEY(Real, crystdock, mininterface)
OPT_1GRP_KEY(Real, crystdock, mininterfacesum)
OPT_1GRP_KEY(Real, crystdock, trans_step)
OPT_1GRP_KEY(Real, crystdock, rot_step)
OPT_1GRP_KEY(Real, crystdock, cb_radius)
OPT_1GRP_KEY(Integer, crystdock, nmodels)
OPT_1GRP_KEY(Integer, crystdock, rotnum)
OPT_1GRP_KEY(Integer, crystdock, symnum)
OPT_1GRP_KEY(Boolean, crystdock, ssonly)
OPT_1GRP_KEY(Boolean, crystdock, debug)
OPT_1GRP_KEY(Boolean, crystdock, debug_exact)
OPT_1GRP_KEY(Boolean, crystdock, eval_native)
OPT_1GRP_KEY(Real, crystdock, n_clashdist)
OPT_1GRP_KEY(Real, crystdock, ca_clashdist)
OPT_1GRP_KEY(Real, crystdock, c_clashdist)
OPT_1GRP_KEY(Real, crystdock, o_clashdist)
OPT_1GRP_KEY(Real, crystdock, cb_clashdist)
OPT_1GRP_KEY(Real, crystdock, sigwidth)
OPT_1GRP_KEY(Real, crystdock, interfacedist)
OPT_1GRP_KEY(Real, crystdock, interface_sigwidth)
OPT_1GRP_KEY(Real, crystdock, cluster_cutoff)
OPT_1GRP_KEY(Real, crystdock, motif_initial_cut)
OPT_1GRP_KEY(Real, crystdock, motif_final_cut)
OPT_1GRP_KEY(Real, crystdock, energy_cut)
OPT_1GRP_KEY(Real, crystdock, favor_native_bonus)
OPT_1GRP_KEY(Real, crystdock, K)
OPT_1GRP_KEY(Integer, crystdock, dock_ncycle)
OPT_1GRP_KEY(Real, crystdock, dock_vdw_cut)
OPT_1GRP_KEY(Integer, crystdock, ncyc)
OPT_1GRP_KEY(Integer, crystdock, nmut)
OPT_1GRP_KEY(IntegerVector, crystdock, fixedres)
OPT_1GRP_KEY(IntegerVector, crystdock, hits_to_dock)
OPT_1GRP_KEY(StringVector, crystdock, adjust_ref_weights)
OPT_1GRP_KEY(Boolean, crystdock, nofilt)
OPT_1GRP_KEY(Boolean, crystdock, nocys)

OPT_1GRP_KEY(String, crystdock, hits_in)
OPT_1GRP_KEY(String, crystdock, hits_out)


///////////////////////////////////////
///////////////////////////////////////

#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.295779513082323


////////////////////////////////////////////////
inline double sign(double x) {return (x >= 0.0) ? 1.0 : -1.0;}
inline double norm4(double a, double b, double c, double d) {return sqrt(a * a + b * b + c * c + d * d);}

////////////////////////////////////////////////
void
adjust_ref_weights( core::scoring::ScoreFunction & scorefxn ) {
	if ( !option[ crystdock::adjust_ref_weights ].user() ) return;

	core::scoring::methods::EnergyMethodOptions options( scorefxn.energy_method_options() );
	utility::vector1< core::Real > ref_weights( options.method_weights( core::scoring::ref ) );

	utility::vector1< std::string > const l( option[ crystdock::adjust_ref_weights ]() );
	runtime_assert( l.size()%2 == 0 );
	core::Size const nwts( l.size()/2 );

	for ( Size ii=0; ii< nwts; ++ii ) {
		core::chemical::AA const aa( core::chemical::aa_from_oneletter_code( l[2*ii+1 ][0] ) );
		Real const adjustment( float_of( l[2*ii+2] ) );
		TR << "Updating reference energy for " << aa << " from " << ref_weights[aa] << " to " << ref_weights[aa]+adjustment << std::endl;
		ref_weights[ aa ] += adjustment;
	}
	options.set_method_weights( core::scoring::ref, ref_weights );
	scorefxn.set_energy_method_options( options );
}

/// backward interpolation (unrotated->rotated)
///    used to generate grids
template <class S>
core::Real interp_linear(
	ObjexxFCL::FArray3D< S > const & data ,
	numeric::xyzVector< core::Real > const & idxX ) {

	numeric::xyzVector<int> pt000, pt111;
	numeric::xyzVector< core::Real > fpart,neg_fpart;
	numeric::xyzVector<int> srcgrid(data.u1(),data.u2(),data.u3());

	pt000[0] = (int)(std::floor(idxX[0]+1e-6));
	pt000[1] = (int)(std::floor(idxX[1]+1e-6));
	pt000[2] = (int)(std::floor(idxX[2]+1e-6));

	// interpolation coeffs
	fpart[0] = idxX[0]-pt000[0]; neg_fpart[0] = 1-fpart[0];
	fpart[1] = idxX[1]-pt000[1]; neg_fpart[1] = 1-fpart[1];
	fpart[2] = idxX[2]-pt000[2]; neg_fpart[2] = 1-fpart[2];
	S retval = (S)0.0;

	// bound check
	if ( pt000[0] < -srcgrid[0]/2 || pt000[0] > srcgrid[0]/2 ) return 0.0;
	if ( pt000[1] < -srcgrid[1]/2 || pt000[1] > srcgrid[1]/2 ) return 0.0;
	if ( pt000[2] < -srcgrid[2]/2 || pt000[2] > srcgrid[2]/2 ) return 0.0;
	if ( pt000[0] < 0 ) pt000[0] += srcgrid[0];
	if ( pt000[1] < 0 ) pt000[1] += srcgrid[1];
	if ( pt000[2] < 0 ) pt000[2] += srcgrid[2];

	pt111[0] = (pt000[0]+1); if ( pt111[0]>=srcgrid[0] ) pt111[0]=0;
	pt111[1] = (pt000[1]+1); if ( pt111[1]>=srcgrid[1] ) pt111[1]=0;
	pt111[2] = (pt000[2]+1); if ( pt111[2]>=srcgrid[2] ) pt111[2]=0;

	retval += neg_fpart[0]*neg_fpart[1]*neg_fpart[2] * data(pt000[0]+1,pt000[1]+1,pt000[2]+1);
	retval += neg_fpart[0]*neg_fpart[1]*    fpart[2] * data(pt000[0]+1,pt000[1]+1,pt111[2]+1);
	retval += neg_fpart[0]*    fpart[1]*neg_fpart[2] * data(pt000[0]+1,pt111[1]+1,pt000[2]+1);
	retval += neg_fpart[0]*    fpart[1]*    fpart[2] * data(pt000[0]+1,pt111[1]+1,pt111[2]+1);
	retval +=     fpart[0]*neg_fpart[1]*neg_fpart[2] * data(pt111[0]+1,pt000[1]+1,pt000[2]+1);
	retval +=     fpart[0]*neg_fpart[1]*    fpart[2] * data(pt111[0]+1,pt000[1]+1,pt111[2]+1);
	retval +=     fpart[0]*    fpart[1]*neg_fpart[2] * data(pt111[0]+1,pt111[1]+1,pt000[2]+1);
	retval +=     fpart[0]*    fpart[1]*    fpart[2] * data(pt111[0]+1,pt111[1]+1,pt111[2]+1);
	return (core::Real)retval;
}


// rotation stuff
void
euler2rot( core::Real a, core::Real b, core::Real g, numeric::xyzMatrix<Real> &R) {
	R.xx() = -sin(DEG2RAD*a)*cos(DEG2RAD*b)*sin(DEG2RAD*g) + cos(DEG2RAD*a)*cos(DEG2RAD*g);
	R.xy() =  cos(DEG2RAD*a)*cos(DEG2RAD*b)*sin(DEG2RAD*g) + sin(DEG2RAD*a)*cos(DEG2RAD*g);
	R.xz() =  sin(DEG2RAD*b)*sin(DEG2RAD*g);
	R.yx() = -sin(DEG2RAD*a)*cos(DEG2RAD*b)*cos(DEG2RAD*g) - cos(DEG2RAD*a)*sin(DEG2RAD*g);
	R.yy() =  cos(DEG2RAD*a)*cos(DEG2RAD*b)*cos(DEG2RAD*g) - sin(DEG2RAD*a)*sin(DEG2RAD*g);
	R.yz() =  sin(DEG2RAD*b)*cos(DEG2RAD*g);
	R.zx() =  sin(DEG2RAD*a)*sin(DEG2RAD*b);
	R.zy() = -cos(DEG2RAD*a)*sin(DEG2RAD*b);
	R.zz() =  cos(DEG2RAD*b);
}

numeric::xyzVector<Real>
center_pose_at_origin( Pose & pose ) {
	numeric::xyzVector<Real> com(0,0,0);
	int count=0;
	for ( int i=1; i<=(int)pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		com += pose.residue(i).atom(2).xyz();
		count++;
	}
	com /= count;

	for ( int i=1; i<=(int)pose.size(); ++i ) {
		for ( int j=1; j<=(int)pose.residue(i).natoms(); ++j ) {
			pose.set_xyz(id::AtomID(j,i), pose.residue(i).atom(j).xyz()-com );
		}
	}
	return com;
}

void
add_crystinfo_to_pose( Pose & pose, Spacegroup const &sg ) {
	CrystInfo ci;
	ci.A(sg.A()); ci.B(sg.B()); ci.C(sg.C()); ci.alpha(sg.alpha()); ci.beta(sg.beta()); ci.gamma(sg.gamma());
	ci.spacegroup(sg.pdbname());
	pose.pdb_info()->set_crystinfo(ci);
}


///////////////////////////////////////////////////////////////////////////////

// fpd this class is used for storing hits over all rotations
//    same top-N priority_queue but much larger store
struct SpacegroupHit {
	SpacegroupHit( std::string line ) {
		std::istringstream iss(line);
		std::string sgname;
		core::Real A,B,C,alpha,beta,gamma;
		iss >> tag >> sgname >> A >> B >> C >> alpha >> beta >> gamma
			>> score >> R.xx() >> R.xy() >> R.xz() >> R.yx() >> R.yy() >> R.yz() >> R.zx() >> R.zy() >> R.zz()
			>> T[0] >> T[1] >> T[2];
		sg = protocols::cryst::Spacegroup(sgname);
		sg.set_parameters(A,B,C,alpha,beta,gamma);
	}

	std::string tag;
	protocols::cryst::Spacegroup sg;
	Real score;
	Real x,y,z;
	numeric::xyzMatrix<Real> R;
	numeric::xyzVector<Real> T;

	std::string
	to_string() {
		std::ostringstream oss;

		oss << tag << " " << sg.name() << " "
			<< sg.A() << " " << sg.B() << " " << sg.C() << " " << sg.alpha() << " " << sg.beta() << " " << sg.gamma() << " "
			<< score << " " << R.xx() << " " << R.xy() << " " << R.xz() << " " << R.yx() << " " << R.yy() << " " << R.yz() << " " << R.zx() << " " << R.zy() << " " << R.zz() << " "
			<< T[0] << " " << T[1] << " " << T[2];

		return oss.str();
	}
};

// fpd this class is used for storing hits over all rotations
//    same top-N priority_queue but much larger store
struct InterfaceHit {
	InterfaceHit( Real score_in, Real x_in, Real y_in, Real z_in, Size rot_index_in, utility::vector1<SingleInterface> iinfo_in ) {
		scoresum=score=score_in;
		x=x_in; y=y_in; z=z_in;
		rot_index=rot_index_in;
		iinfo = iinfo_in;
	}
	InterfaceHit( Real score_in, Real scoresum_in, Real x_in, Real y_in, Real z_in, Size rot_index_in, utility::vector1<SingleInterface> iinfo_in ) {
		score=score_in;
		scoresum=scoresum_in;
		x=x_in; y=y_in; z=z_in;
		rot_index=rot_index_in;
		iinfo = iinfo_in;
	}

	Real score,scoresum;
	Real x,y,z;
	Size rot_index;
	utility::vector1<SingleInterface> iinfo;

	std::string
	to_string(numeric::UniformRotationSampler const &urs) {
		numeric::xyzMatrix<Real> R;
		urs.get(rot_index, R);
		numeric::xyzVector<Real> T (x, y, z);

		std::ostringstream oss;
		oss << score << " " << R.xx() << " " << R.xy() << " " << R.xz() << " " << R.yx() << " " << R.yy() << " " << R.yz() << " " << R.zx() << " " << R.zy() << " " << R.zz() << " "
			<< T[0] << " " << T[1] << " " << T[2];
		return oss.str();
	}
};

class InterfaceHitComparitor {
public:
	bool operator()(InterfaceHit& t1, InterfaceHit& t2) {
		return (t1.score > t2.score);
	}
};

class InterfaceHitDatabase  {
private:
	core::Size N_;
	std::priority_queue<InterfaceHit, std::vector<InterfaceHit>, InterfaceHitComparitor> queue_;

public:
	InterfaceHitDatabase(int N_in) { N_ = N_in; }

	void add_interface( InterfaceHit interface ) {
		if ( queue_.size() <N_ ) {
			queue_.push( interface );
		} else if ( interface.score > queue_.top().score ) {
			queue_.pop();
			queue_.push( interface );
		}
	}

	InterfaceHit pop() {
		InterfaceHit retval = queue_.top();
		queue_.pop();
		return retval;
	}

	InterfaceHit top() {
		InterfaceHit retval = queue_.top();
		return retval;
	}

	Size size() { return queue_.size(); }
};



///////////////////////////////////////////////////////////////////////////////
class CrystRMS : public protocols::moves::Mover {
public:
	CrystRMS() { }

	virtual std::string get_name() const { return "CrystRMS"; }

	void apply( Pose & pose) {
		if ( !native ) {
			native = core::pose::PoseOP( new core::pose::Pose() );
			core::import_pose::pose_from_file( *native, option[in::file::native]() , core::import_pose::PDB_file);
		}

		core::Real rms = protocols::cryst::crystRMS( *native, pose ); // will symmetrize once
		setPoseExtraScore( pose, "cryst_rms", rms );
		TR << "RMS = " << rms << std::endl;
	}

	core::pose::PoseOP native;
};


///////////////////////////////////////////////////////////////////////////////
class CrystDesign : public protocols::moves::Mover {
public:
	CrystDesign(bool revertonly) {
		revertonly_ = revertonly;

		DDG_CUT_ = -10.0;
		WEAKE_CUT_ = -5.0;
		WFY_MUT_CUT_ = 1, M_MUT_CUT_ = 1;
		UNSATS_SOFT_CUT_ = 2.5;
		UNSATS_CUT_ = 1.1;
		SC_CUT_ = 0.4;
		MUT_CUT_ = option[crystdock::nmut]();
		K_ = option[crystdock::K]();
		NCYC_ = option[crystdock::ncyc]();
		NOFILT_ = option[crystdock::nofilt]();
	}

	virtual std::string get_name() const { return "CrystDesign"; }

	void
	filter_and_report( core::pose::Pose & pose, core::scoring::ScoreFunctionOP sf, bool &pass, bool &need_greedy, bool runfilter, bool report ) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::chemical;
		using namespace core::conformation::symmetry;

		bool filter = (!NOFILT_ && runfilter);

		pass = true;
		need_greedy=false;

		//
		SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
		//core::Size nsubunits = symm_info->subunits();
		core::Size nres_asu = symm_info->num_independent_residues();

		// check 1: min interface E
		core::Real total_score, total_interfaceE, best_score_weakint;
		get_interface_energies( setup_, pose, sf, total_score, total_interfaceE, best_score_weakint);

		if ( filter && best_score_weakint > WEAKE_CUT_ ) {   // HARD CODED FILTER
			TR << "Fail on weakE! " << best_score_weakint << " v " << WEAKE_CUT_ << std::endl;
			pass = false;
			return;
		}

		// ddg + unsats + sasa + sc
		utility::vector1<core::Real> temp;
		core::Real easymm, esymm;
		core::Real ddgv = ddg( pose, sf, temp, esymm, easymm, NCYC_, 3 );

		if ( filter &&  ddgv > DDG_CUT_  ) {   // HARD CODED FILTER
			TR << "Fail on ddg! " << ddgv << " v " << DDG_CUT_ << std::endl;
			pass = false;
			return;
		}

		// # mutations
		core::Size nmut=0, nFmut=0, nYmut=0, nWmut=0, nMmut=0;
		std::string before_seq = native_.sequence();
		std::string after_seq = pose.sequence();
		before_seq = before_seq.substr( 0, nres_asu );
		after_seq = after_seq.substr( 0, nres_asu );
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			if ( before_seq[i] != after_seq[i] ) {
				nmut++;
				if ( after_seq[i] == 'W' ) nWmut++;
				if ( after_seq[i] == 'F' ) nFmut++;
				if ( after_seq[i] == 'Y' ) nYmut++;
				if ( after_seq[i] == 'M' ) nMmut++;
			}
		}

		if ( filter &&  ( (nmut > MUT_CUT_) || (nWmut+nFmut+nYmut > WFY_MUT_CUT_) || (nMmut>M_MUT_CUT_) ) ) {   // HARD CODED FILTER
			TR << "Fail res identity!" << std::endl;
			TR << "nmuts */W/F/Y/M " << nmut << "/" << nWmut << "/" << nFmut << "/" << nYmut << "/" << nMmut << std::endl;
			pass = false;
			return;
		}

		utility::vector1< bool > scs_buried;
		core::Real sc_tot, sc_weak, sasa_tot, sasa_weak;
		get_sc( setup_, pose, sc_tot, sc_weak, sasa_tot, sasa_weak);

		core::Real unsat10 = unsatisfied_buried_polars( pose, sf, 1);
		core::Real unsat12 = unsatisfied_buried_polars( pose, sf, 1.2);
		core::Real unsat14 = unsatisfied_buried_polars( pose, sf, 1.4);
		core::Real total_sasa, weakest_sasa, total_nonpolar_sasa, weakest_nonpolar_sasa;
		get_sasa( setup_, pose, total_sasa, weakest_sasa, total_nonpolar_sasa, weakest_nonpolar_sasa);

		// HARD CODED FILTER
		if ( filter && sc_tot < SC_CUT_ ) {
			TR << "Fail on sc! " << sc_tot << " v " << SC_CUT_ << std::endl;
			pass = false;
			return;
		}

		if ( filter && unsat10 > UNSATS_SOFT_CUT_ ) {
			TR << "Fail on unsat10! " << unsat10 << " v " << UNSATS_SOFT_CUT_ << std::endl;
			pass = false;
			return;
		}

		// >>> greedy opt unsats if we are close <<<
		if ( filter && unsat10 > UNSATS_CUT_ ) {
			TR << "High unsat10 (" << unsat10 << ")" << std::endl;
			pass = false;
			need_greedy = true;
			return;
		}


		TR << "Passed all filters!" << std::endl;
		TR << "   ENERGIES:   totE "+utility::to_string(total_score)+"   intE "+utility::to_string(total_interfaceE)
			+"   weakE "+utility::to_string(best_score_weakint);
		TR << "   SC:   sc "+utility::to_string(sc_tot)+"   sc_weak "+utility::to_string(sc_weak)+
			"   sasa "+utility::to_string(sasa_tot)+"   sasa_weak "+utility::to_string(sasa_weak);
		TR << "   MUTATIONS:   nmuts "+utility::to_string(nmut)+
			"  W/F/Y/M "+utility::to_string(nWmut)+" / "+utility::to_string(nFmut)+
			" / "+utility::to_string(nYmut)+" / "+utility::to_string(nMmut);
		TR << "   DDG:   ddg "+utility::to_string(ddgv)+
			"   esymm "+utility::to_string(esymm)+"   easymm "+utility::to_string(easymm);
		TR << "   UNSATS:   unsat10 "+utility::to_string(unsat10)+
			"   unsat12 "+utility::to_string(unsat12)+"   unsat14 "+utility::to_string(unsat14);
		TR << "   SASA:   total_sasa "+utility::to_string(total_sasa)+"   weak_sasa "+utility::to_string(weakest_sasa)+
			"   total_np_sasa "+utility::to_string(total_nonpolar_sasa)+"   weak_np_sasa "+utility::to_string(weakest_nonpolar_sasa);
		TR << std::endl;

		if ( !report ) {
			return;
		}

		// report in PDB header
		core::io::RemarkInfo remark;
		remark.num = 1; remark.value = "ENERGIES:   totE "+utility::to_string(total_score)+"   intE "+utility::to_string(total_interfaceE)
			+"   weakE "+utility::to_string(best_score_weakint);
		pose.pdb_info()->remarks().push_back( remark );
		remark.num = 1; remark.value =
			"SC:   sc "+utility::to_string(sc_tot)+"   sc_weak "+utility::to_string(sc_weak)+
			"   sasa "+utility::to_string(sasa_tot)+"   sasa_weak "+utility::to_string(sasa_weak);
		pose.pdb_info()->remarks().push_back( remark );
		remark.num = 1; remark.value =
			"MUTATIONS:   nmuts "+utility::to_string(nmut)+
			"  W/F/Y/M "+utility::to_string(nWmut)+" / "+utility::to_string(nFmut)+
			" / "+utility::to_string(nYmut)+" / "+utility::to_string(nMmut);
		pose.pdb_info()->remarks().push_back( remark );
		remark.num = 1; remark.value =
			"DDG:   ddg "+utility::to_string(ddgv)+
			"   esymm "+utility::to_string(esymm)+"   easymm "+utility::to_string(easymm);
		pose.pdb_info()->remarks().push_back( remark );
		remark.num = 1; remark.value ="UNSATS:   unsat10 "+utility::to_string(unsat10)+
			"   unsat12 "+utility::to_string(unsat12)+"   unsat14 "+utility::to_string(unsat14);
		pose.pdb_info()->remarks().push_back( remark );
		remark.num = 1; remark.value =
			"SASA:   total_sasa "+utility::to_string(total_sasa)+"   weak_sasa "+utility::to_string(weakest_sasa)+
			"   total_np_sasa "+utility::to_string(total_nonpolar_sasa)+"   weak_np_sasa "+utility::to_string(weakest_nonpolar_sasa);
		pose.pdb_info()->remarks().push_back( remark );

	}

	void design_cycle( Pose & pose, core::scoring::ScoreFunctionOP sf, core::scoring::ScoreFunctionOP sfsoft, utility::vector1<bool> interface ) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::chemical;
		using namespace core::conformation::symmetry;

		utility::vector1<int> fixedres = option[crystdock::fixedres];

		///
		/// SETUP
		///
		SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		SymmetryInfoCOP symm_info = ( SymmConf.Symmetry_Info() );
		core::Size nres_asu = symm_info->num_independent_residues();

		// disallow PG mutations
		utility::vector1<bool> allowed_aas( num_canonical_aas, true );
		allowed_aas[ aa_pro ] = false;
		allowed_aas[ aa_gly ] = false;
		if ( option[crystdock::nocys]() ) { allowed_aas[ aa_cys ] = false; }

		core::pack::task::PackerTaskOP designtask = core::pack::task::TaskFactory::create_packer_task( pose );

		designtask->restrict_to_residues(interface);

		designtask->initialize_from_command_line();  // ex1 ex2 from command line
		if ( option[ packing::resfile ].user() ) {
			parse_resfile( pose, *designtask);
		}
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			designtask->nonconst_residue_task( i ).restrict_absent_canonical_aas( allowed_aas );
			designtask->nonconst_residue_task( i ).or_include_current( true );
		}

		// dont repack fixed residues
		for ( core::Size i=1; i<=fixedres.size(); ++i ) {
			designtask->nonconst_residue_task( fixedres[i] ).prevent_repacking(  );
		}

		kinematics::MoveMapOP mm (new kinematics::MoveMap());
		mm->set_bb  ( false ); mm->set_chi ( false ); mm->set_jump( true );
		for ( int i=1; i<=(int)nres_asu; ++i ) mm->set_chi(i, interface[i]);
		core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
		options.nblist_auto_update( false );
		core::optimization::symmetry::SymAtomTreeMinimizer minimizer;

		///
		/// DESIGN
		///
		// 1 soft design+min+hard pack+min
		{
			protocols::simple_moves::symmetry::SymPackRotamersMoverOP design_soft (
				new protocols::simple_moves::symmetry::SymPackRotamersMover( sfsoft, designtask ) );
			design_soft->apply( pose );
			minimizer.run( pose, *mm, *sf, options );

			core::pack::task::PackerTaskOP repacktask = core::pack::task::TaskFactory::create_packer_task( pose );
			repacktask->restrict_to_residues(interface);
			repacktask->restrict_to_repacking();
			repacktask->initialize_from_command_line();  // ex1 ex2 from command line
			for ( core::Size i=1; i<=nres_asu; ++i ) repacktask->nonconst_residue_task( i ).or_include_current( true );
			for ( core::Size i=1; i<=fixedres.size(); ++i ) repacktask->nonconst_residue_task( fixedres[i] ).prevent_repacking(  );

			protocols::simple_moves::symmetry::SymPackRotamersMoverOP repack_hard (
				new protocols::simple_moves::symmetry::SymPackRotamersMover( sf, repacktask ) );
			repack_hard->apply( pose );
			minimizer.run( pose, *mm, *sf, options );
		}

		// 2 perturb docking
		protocols::cryst::DockLatticeMover docking;
		docking.set_trans_mag(0.2);
		docking.set_rot_mag(1.0);
		docking.set_ncycles(option[crystdock::dock_ncycle]());
		docking.set_fullatom(true);
		docking.set_scorefunction(sf);
		docking.set_temp(0.5); // low temp

		if ( option[crystdock::dock_ncycle]() > 0 ) {
			docking.apply(pose);
		}

		// 3 hard design+min+hard pack+min
		{
			protocols::simple_moves::symmetry::SymPackRotamersMoverOP design_hard (
				new protocols::simple_moves::symmetry::SymPackRotamersMover( sf, designtask ) );
			design_hard->apply( pose );
			minimizer.run( pose, *mm, *sf, options );

			core::pack::task::PackerTaskOP repacktask = core::pack::task::TaskFactory::create_packer_task( pose );
			repacktask->restrict_to_residues(interface);
			repacktask->restrict_to_repacking();
			repacktask->initialize_from_command_line();  // ex1 ex2 from command line
			for ( core::Size i=1; i<=nres_asu; ++i ) repacktask->nonconst_residue_task( i ).or_include_current( true );
			for ( core::Size i=1; i<=fixedres.size(); ++i ) repacktask->nonconst_residue_task( fixedres[i] ).prevent_repacking(  );

			protocols::simple_moves::symmetry::SymPackRotamersMoverOP repack_hard (
				new protocols::simple_moves::symmetry::SymPackRotamersMover( sf, repacktask ) );
			repack_hard->apply( pose );
			minimizer.run( pose, *mm, *sf, options );
		}
	}

	//
	void
	get_interface_residues( Pose & pose, core::scoring::ScoreFunctionOP sf, utility::vector1< bool > &interface) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::conformation::symmetry;

		// interface parameters
		//   if (angle_degrees > K*exp(b*distance_angstrom) we have an interface interaction
		//   reduce K to increase # detected interface residues
		//      max dist = (1/b)*ln(180/k)
		//   at K= 10, maxdist = 10.32
		//         12             9.67
		//         14             9.12
		//         16             8.64
		core::Real K=K_;
		core::Real b=0.28;

		utility::vector1<int> fixedres = option[crystdock::fixedres];

		// symminfo
		SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
		//core::Size nsubunits = symm_info->subunits();
		core::Size nres_asu = symm_info->num_independent_residues();

		utility::vector1< Size > protein_residues;
		for ( Size i=1, i_end = nres_asu; i<= i_end; ++i ) {
			if ( pose.residue( i ).is_protein() ) {
				protein_residues.push_back( i );
			}
		}

		interface.clear();
		interface.resize( pose.size(), false );

		// initialize energy graph
		(*sf)(pose);
		Energies & energies( pose.energies() );
		EnergyGraph & energy_graph( energies.energy_graph() );

		// stage 1
		//   find interface residues
		for ( Size i=1, i_end = nres_asu; i<= i_end; ++i ) {
			conformation::Residue const & rsd1( pose.residue( i ) );

			if ( std::find(fixedres.begin(), fixedres.end(), i)!=fixedres.end() ) continue;

			// ignore native gly / pro
			core::chemical::AA aa_src = pose.residue(i).aa();
			if ( aa_src == core::chemical::aa_pro || aa_src == core::chemical::aa_gly ) continue;

			for ( utility::graph::Graph::EdgeListIter
					iru  = energy_graph.get_node(i)->edge_list_begin(),
					irue = energy_graph.get_node(i)->edge_list_end();
					iru != irue; ++iru ) {
				EnergyEdge & edge( static_cast< EnergyEdge & > (**iru) );

				Size const j = edge.get_other_ind( i );
				conformation::Residue const & rsd2( pose.residue( j ) );

				if ( i==j ) continue;  // don't think this is necessary
				if ( !rsd1.is_protein() || !rsd2.is_protein() ) continue;
				if ( symm_info->bb_is_independent(rsd2.seqpos()) ) continue;  // only care about interactions across the symmetric interface

				// if CB-CB distance < 8A and CA-CB-CB angles are >75 deg then design
				core::Real dist, angle1, angle2;
				if ( rsd2.aa() == core::chemical::aa_gly ) {
					dist = (rsd1.atom("CB").xyz() - rsd2.atom("CA").xyz()).length() + 1.0;
					angle1 = numeric::angle_degrees(rsd1.atom("CA").xyz(), rsd1.atom("CB").xyz(), rsd2.atom("CA").xyz() ) ;
					angle2 = 180.0;
				} else {
					dist = (rsd1.atom("CB").xyz() - rsd2.atom("CB").xyz()).length();
					angle1 = numeric::angle_degrees(rsd1.atom("CA").xyz(), rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz() ) ;
					angle2 = numeric::angle_degrees(rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz(), rsd2.atom("CA").xyz() ) ;
				}
				core::Real angle_tgt = K*exp(b*dist);

				if ( angle_tgt < 180 && angle1 > angle_tgt && angle2 > angle_tgt ) {
					interface[i] = true;
				}
			}
		}
	}


	//
	void
	get_neighbor_residues(
		Pose & pose,
		core::Size res_i,
		utility::vector1< bool > &neighbor)
	{
		using namespace core;
		using namespace core::scoring;
		using namespace core::conformation;
		using namespace core::conformation::symmetry;

		// interface parameters
		//   same logic as above
		core::Real K=option[crystdock::K]();
		core::Real b=0.28;

		// symminfo
		SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
		//core::Size nsubunits = symm_info->subunits();
		core::Size nres_asu = symm_info->num_independent_residues();

		// make pose polyA
		utility::vector1< Size > protein_residues;
		for ( Size i=1, i_end = nres_asu; i<= i_end; ++i ) {
			if ( pose.residue( i ).is_protein() ) protein_residues.push_back( i );
		}

		neighbor.clear();
		neighbor.resize( pose.size(), false );


		conformation::Residue const & rsd1( pose.residue( res_i ) );
		for ( Size j=1, j_end = pose.size(); j<= j_end; ++j ) {
			conformation::Residue const & rsd2( pose.residue( j ) );

			if ( res_i==j ) continue;
			if ( !rsd1.is_protein() || !rsd2.is_protein() ) continue;

			core::Real dist, angle1, angle2;
			if ( rsd1.aa() == core::chemical::aa_gly && rsd2.aa() == core::chemical::aa_gly ) {
				dist = (rsd1.atom("CA").xyz() - rsd2.atom("CA").xyz()).length() + 2.0;
				angle1 = 180.0;
				angle2 = 180.0;
			} else if ( rsd1.aa() == core::chemical::aa_gly ) {
				dist = (rsd1.atom("CB").xyz() - rsd2.atom("CA").xyz()).length() + 1.0;
				angle1 = 180.0;
				angle2 = numeric::angle_degrees(rsd1.atom("CA").xyz(), rsd2.atom("CB").xyz(), rsd2.atom("CA").xyz() ) ;
			} else if ( rsd2.aa() == core::chemical::aa_gly ) {
				dist = (rsd1.atom("CB").xyz() - rsd2.atom("CA").xyz()).length() + 1.0;
				angle1 = numeric::angle_degrees(rsd1.atom("CA").xyz(), rsd1.atom("CB").xyz(), rsd2.atom("CA").xyz() ) ;
				angle2 = 180.0;
			} else {
				dist = (rsd1.atom("CB").xyz() - rsd2.atom("CB").xyz()).length();
				angle1 = numeric::angle_degrees(rsd1.atom("CA").xyz(), rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz() ) ;
				angle2 = numeric::angle_degrees(rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz(), rsd2.atom("CA").xyz() ) ;
			}

			core::Real angle_tgt = K*exp(b*dist);

			if ( angle_tgt < 180 && angle1 > angle_tgt && angle2 > angle_tgt ) {
				core::Size j_master = j;
				if ( !symm_info->bb_is_independent(j) ) j_master=symm_info->bb_follows(j);
				neighbor[j_master]=true;
			}
		}
	}


	void do_reversion( Pose & pose, core::scoring::ScoreFunctionOP sf, core::scoring::ScoreFunctionOP sfsoft ) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::conformation;
		using namespace core::conformation::symmetry;

		//
		SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
		//core::Size nsubunits = symm_info->subunits();
		core::Size nres_asu = symm_info->num_independent_residues();

		// get baseline stats
		core::Real total_score, interfaceE, weakintE;
		core::Real easymm,esymm,ddgv,sc_tot,sc_weak,sasa_tot,sasa_weak,unsat10;
		core::Real total_sasa, weakest_sasa, total_nonpolar_sasa, weakest_nonpolar_sasa;
		utility::vector1<core::Real> temp;

		get_interface_energies( setup_, pose, sf, total_score, interfaceE, weakintE);
		ddgv = ddg( pose, sf, temp, esymm, easymm, option[crystdock::ncyc](), 1 ); // fast...
		get_sc( setup_, pose, sc_tot, sc_weak, sasa_tot, sasa_weak);
		unsat10 = unsatisfied_buried_polars( pose, sf, 1);
		get_sasa( setup_, pose, total_sasa, weakest_sasa, total_nonpolar_sasa, weakest_nonpolar_sasa);

		core::Real interfaceE_baseline=interfaceE;
		core::Real weakintE_baseline=weakintE;
		core::Real ddgv_baseline=ddgv;
		core::Real unsat10_baseline=unsat10;
		core::Real sc_tot_baseline=sc_tot;
		core::Real sc_weak_baseline=sc_weak;

		// get mutations
		utility::vector1< core::Size > res_to_mutate;
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			if ( pose.residue(i).aa() != native_.residue(i).aa() ) res_to_mutate.push_back(i);
		}

		numeric::random::random_permutation( res_to_mutate, numeric::random::rg() );

		for ( int i=1; i<=(int)res_to_mutate.size(); ++i ) {
			core::pose::Pose pose_curr = pose;

			// mutate
			pose.replace_residue(res_to_mutate[i], native_.residue(res_to_mutate[i]), true );

			// setup
			utility::vector1< bool > neighbor_i;
			get_neighbor_residues( pose, res_to_mutate[i],neighbor_i);
			neighbor_i[ res_to_mutate[i] ] = true;
			core::pack::task::PackerTaskOP repacktask = core::pack::task::TaskFactory::create_packer_task( pose );
			repacktask->restrict_to_residues(neighbor_i);
			repacktask->restrict_to_repacking();
			repacktask->initialize_from_command_line();  // ex1 ex2 from command line

			kinematics::MoveMapOP mm (new kinematics::MoveMap());
			mm->set_bb  ( false ); mm->set_chi ( false ); mm->set_jump( true );
			for ( int j=1; j<=(int)nres_asu; ++j ) mm->set_chi(j, neighbor_i[j]);
			core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
			options.nblist_auto_update( false );
			core::optimization::symmetry::SymAtomTreeMinimizer minimizer;

			// repack interface + rbmin
			protocols::simple_moves::symmetry::SymPackRotamersMoverOP repack_soft (
				new protocols::simple_moves::symmetry::SymPackRotamersMover( sfsoft, repacktask ) );
			repack_soft->apply( pose );
			minimizer.run( pose, *mm, *sf, options );

			// repack interface + rbmin
			protocols::simple_moves::symmetry::SymPackRotamersMoverOP repack (
				new protocols::simple_moves::symmetry::SymPackRotamersMover( sf, repacktask ) );
			repack->apply( pose );
			minimizer.run( pose, *mm, *sf, options );

			// eval
			get_interface_energies( setup_, pose, sf, total_score, interfaceE, weakintE);
			ddgv = ddg( pose, sf, temp, esymm, easymm, option[crystdock::ncyc](), 1 ); // fast...
			get_sc( setup_, pose, sc_tot, sc_weak, sasa_tot, sasa_weak);
			unsat10 = unsatisfied_buried_polars( pose, sf, 1);
			get_sasa( setup_, pose, total_sasa, weakest_sasa, total_nonpolar_sasa, weakest_nonpolar_sasa);

			//
			if (
					1.0  + interfaceE_baseline > interfaceE &&
					1.0  + ddgv_baseline > ddgv &&
					0.5  + weakintE_baseline > weakintE &&
					0.25 + unsat10_baseline > unsat10 &&  // dangerous?
					-0.1 + sc_tot_baseline < sc_tot
					) {
				TR << "ACCEPT reversion at " << res_to_mutate[i] << " with"
					<< " intE/weakE = " << interfaceE << " / " << weakintE
					<< " ddg = " << ddgv
					<< " unsat10 = " << unsat10
					<< " sc/scweak = " << sc_tot << " / " << sc_weak << std::endl;

				// update
				interfaceE_baseline = interfaceE;
				ddgv_baseline = ddgv;
				weakintE_baseline = weakintE;
				unsat10_baseline = unsat10;
				sc_tot_baseline = sc_tot;
			} else {
				TR << "REJECT reversion at " << res_to_mutate[i] << " with"
					<< " intE/weakE = " << interfaceE << " / " << weakintE
					<< " ddg = " << ddgv
					<< " unsat10 = " << unsat10
					<< " sc/scweak = " << sc_tot << " / " << sc_weak << std::endl;
				TR << "                WAS"
					<< " intE/weakE = " << interfaceE_baseline << " / " << weakintE_baseline
					<< " ddg = " << ddgv_baseline
					<< " unsat10 = " << unsat10_baseline
					<< " sc/scweak = " << sc_tot_baseline << " / " << sc_weak_baseline << std::endl;
				pose = pose_curr;
			}
		}
	}

	void greedy_revert_unsats( Pose & pose, core::scoring::ScoreFunctionOP sf, core::scoring::ScoreFunctionOP sfsoft ) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::conformation;
		using namespace core::conformation::symmetry;

		//
		SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
		//core::Size nsubunits = symm_info->subunits();
		core::Size nres_asu = symm_info->num_independent_residues();

		// get baseline stats
		utility::vector1< bool > sc_unsats;
		core::Real total_score, interfaceE, weakintE;

		get_interface_energies( setup_, pose, sf, total_score, interfaceE, weakintE);
		core::Real unsat10 = unsatisfied_buried_polars( pose, sf, 1, sc_unsats);

		core::Real interfaceE_baseline=interfaceE;
		core::Real unsat10_baseline=unsat10;

		// get mutations
		utility::vector1< core::Size > res_to_mutate;
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			if ( sc_unsats[i] ) res_to_mutate.push_back(i);
		}

		numeric::random::random_permutation( res_to_mutate, numeric::random::rg() );

		for ( int i=1; i<=(int)res_to_mutate.size(); ++i ) {
			core::pose::Pose pose_curr = pose;

			for ( int j=1; j<=20; ++j ) {
				if ( ((core::chemical::AA)j)==core::chemical::aa_pro || ((core::chemical::AA)j)==core::chemical::aa_gly ) continue;

				// mutate
				core::chemical::ResidueTypeCOP restype =
					core::chemical::ChemicalManager::get_instance()->residue_type_set(
					core::chemical::FA_STANDARD )->get_representative_type_aa((core::chemical::AA)j);
				core::conformation::ResidueOP rsdnew = ResidueFactory::create_residue( *(restype) );
				pose.replace_residue(res_to_mutate[i], *rsdnew, true );

				// setup
				utility::vector1< bool > neighbor_i;
				get_neighbor_residues( pose, res_to_mutate[i],neighbor_i);
				neighbor_i[ res_to_mutate[i] ] = true;
				core::pack::task::PackerTaskOP repacktask = core::pack::task::TaskFactory::create_packer_task( pose );
				repacktask->restrict_to_residues(neighbor_i);
				repacktask->restrict_to_repacking();
				repacktask->initialize_from_command_line();  // ex1 ex2 from command line

				kinematics::MoveMapOP mm (new kinematics::MoveMap());
				mm->set_bb  ( false ); mm->set_chi ( false ); mm->set_jump( true );
				for ( int k=1; k<=(int)nres_asu; ++k ) mm->set_chi(k, neighbor_i[j]);
				core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
				options.nblist_auto_update( false );
				core::optimization::symmetry::SymAtomTreeMinimizer minimizer;

				// repack interface + rbmin
				protocols::simple_moves::symmetry::SymPackRotamersMoverOP repack_soft (
					new protocols::simple_moves::symmetry::SymPackRotamersMover( sfsoft, repacktask ) );
				repack_soft->apply( pose );
				minimizer.run( pose, *mm, *sf, options );

				// repack interface + rbmin
				protocols::simple_moves::symmetry::SymPackRotamersMoverOP repack (
					new protocols::simple_moves::symmetry::SymPackRotamersMover( sf, repacktask ) );
				repack->apply( pose );
				minimizer.run( pose, *mm, *sf, options );

				// eval
				get_interface_energies( setup_, pose, sf, total_score, interfaceE, weakintE);
				unsat10 = unsatisfied_buried_polars( pose, sf, 1);

				if ( unsat10 < unsat10_baseline ) {
					// accept
					TR << "ACCEPT mutation to " << ((core::chemical::AA)j) <<  " at " << res_to_mutate[i] << " with"
						<< " unsat10 = " << unsat10
						<< " intE/weakE = " << interfaceE << " / " << weakintE << std::endl;
					interfaceE_baseline = interfaceE;
					unsat10_baseline = unsat10;
				} else if ( unsat10 == unsat10_baseline && interfaceE < interfaceE_baseline ) {
					// accept
					TR << "WEAK ACCEPT mutation to " << ((core::chemical::AA)j) <<  " at " << res_to_mutate[i] << " with"
						<< " unsat10 = " << unsat10
						<< " intE/weakE = " << interfaceE << " / " << weakintE << std::endl;
					interfaceE_baseline = interfaceE;
					unsat10_baseline = unsat10;
				} else {
					TR << "REJECT mutation to " << j <<  " at " << res_to_mutate[i] << " with"
						<< " unsat10 = " << unsat10
						<< " intE/weakE = " << interfaceE << " / " << weakintE << std::endl;
					TR << "                WAS"
						<< " unsat10 = " << unsat10_baseline
						<< " intE/weakE = " << interfaceE_baseline << " / na "  << std::endl;
					pose = pose_curr;
				}
			}
		}
	}

	void apply( Pose & pose) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::chemical;
		using namespace core::conformation::symmetry;

		native_ = pose;

		//// //////////////// ////
		////                  ////
		////      SETUP       ////
		////                  ////
		//// //////////////// ////
		// 1: mutate to native
		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_file( native_, option[in::file::native]() , core::import_pose::PDB_file);
			runtime_assert( native_.size() == pose.size() );

			if ( !revertonly_ ) {
				for ( core::Size i=1; i<pose.size(); ++i ) {
					pose.replace_residue(i, native_.residue(i), true );
				}
			}
		}

		runtime_assert( ! core::pose::symmetry::is_symmetric( pose ) );
		setup_.contact_dist(12);
		setup_.apply( pose );

		//
		SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
		//core::Size nsubunits = symm_info->subunits();
		core::Size nres_asu = symm_info->num_independent_residues();

		core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
		sf->set_weight( core::scoring::res_type_constraint , 1.0 );
		core::scoring::ScoreFunctionOP sfsoft = core::scoring::get_score_function();
		sfsoft->set_weight( core::scoring::fa_rep , 0.02 );
		sfsoft->set_weight( core::scoring::res_type_constraint , 1.0 );

		adjust_ref_weights( *sfsoft );
		adjust_ref_weights( *sf );

		for ( core::Size i=1; i<= nres_asu;  ++i ) {
			pose.add_constraint(
				core::scoring::constraints::ConstraintOP(
				new core::scoring::constraints::ResidueTypeConstraint( native_, i,  option[ crystdock::favor_native_bonus]()) ) );
		}

		// 2: find interface residues
		utility::vector1<bool> interface;
		get_interface_residues( pose, sf, interface);
		int ndes=0;
		for ( int i=1; i<=(int)interface.size(); ++i ) {
			if ( interface[i] ) {
				ndes++;
			}
		}

		bool pass, need_greedy;
		if ( !revertonly_ ) {
			//// //////////////// ////
			////                  ////
			////      DESIGN      ////
			////                  ////
			//// //////////////// ////
			design_cycle( pose, sf, sfsoft, interface );

			filter_and_report( pose, sf, pass, need_greedy, true, false );

			//// //////////////// ////
			////                  ////
			////      GREEDY      ////
			////                  ////
			//// //////////////// ////
			if ( !pass && need_greedy ) {
				TR << "Running greedy..." << std::endl;
				greedy_revert_unsats( pose, sf, sfsoft );
				TR << "Greedy opt complete.  Rerun filters:" << std::endl;
				filter_and_report( pose, sf, pass, need_greedy, true, false );
			}

			if ( !pass ) {
				TR << "FAILED!" << std::endl;
				set_last_move_status( protocols::moves::FAIL_RETRY );
				return;
			}
		}

		//// //////////////// ////
		////                  ////
		////      REVERT      ////
		////                  ////
		//// //////////////// ////
		TR << "All filters passed.  Running reversion." << std::endl;
		if ( option[ in::file::native ].user() ) {
			do_reversion(pose, sf, sfsoft );
		}

		// finally, report stats
		filter_and_report( pose, sf, pass, need_greedy, false, true );

		// crystinfo stuff
		protocols::cryst::UpdateCrystInfo update_cryst1;
		update_cryst1.apply(pose);
	}

private:
	protocols::cryst::MakeLatticeMover setup_;

	// filters
	core::Real DDG_CUT_,WEAKE_CUT_;
	core::Real UNSATS_SOFT_CUT_,UNSATS_CUT_,SC_CUT_;
	core::Size WFY_MUT_CUT_,M_MUT_CUT_, MUT_CUT_;

	// interface definition
	core::Real K_;

	// other options
	core::Size NCYC_;
	bool NOFILT_;

	// native
	core::pose::Pose native_;

	bool revertonly_;  // run in revert mode?
};


///////////////////////////////////////////////////////////////////////////////
class CrystRelax : public protocols::moves::Mover {
public:
	CrystRelax(bool relax=true) {
		passed_filts=false;
		relax_=relax;
	}

	virtual std::string get_name() const { return "CrystRelax"; }

	void apply( Pose & pose) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::conformation::symmetry;

		runtime_assert( ! core::pose::symmetry::is_symmetric( pose ) );
		protocols::cryst::MakeLatticeMover setup;
		setup.contact_dist(12);
		setup.apply( pose );

		core::kinematics::MoveMapOP movemap (new core::kinematics::MoveMap);
		movemap->set_bb(false); movemap->set_chi(true); movemap->set_jump(true);
		//movemap->set( core::id::THETA, true );

		core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
		//adjust_ref_weights_from_command_line( *sf );
		protocols::relax::FastRelax relax_prot( sf, option[ crystdock::ncyc ]() );
		relax_prot.min_type("lbfgs_armijo_nonmonotone");
		relax_prot.set_movemap( movemap );
		if ( relax_ ) relax_prot.apply(pose);

		core::Real total_score, total_interfaceE, best_score_weakint;
		get_interface_energies( setup, pose, sf, total_score, total_interfaceE, best_score_weakint);
		core::io::RemarkInfo remark;
		remark.num = 1; remark.value = "score(relax) = "+utility::to_string(total_score)+" | "+utility::to_string(total_interfaceE)
			+" | "+utility::to_string(best_score_weakint);
		pose.pdb_info()->remarks().push_back( remark );

		utility::vector1<core::Real> temp;
		core::Real easymm, esymm;
		core::Real ddgv = ddg( pose, sf, temp, esymm, easymm, option[crystdock::ncyc](), 3 );
		remark.num = 1; remark.value =
			"DDG:   ddg "+utility::to_string(ddgv)+
			"   esymm "+utility::to_string(esymm)+
			"   easymm "+utility::to_string(easymm);
		pose.pdb_info()->remarks().push_back( remark );

		passed_filts =
			(best_score_weakint < option[ crystdock::energy_cut ]());

		protocols::cryst::UpdateCrystInfo update_cryst1;
		update_cryst1.apply(pose);
	}

	bool passed_filts, relax_;
};

///////////////////////////////////////////////////////////////////////////////
class CrystCluster : public protocols::moves::Mover {
private:
public:
	CrystCluster() {}

	Pose
	transform_pdb( Pose const &pose, SpacegroupHit sh ) {
		Pose posecopy = pose;
		posecopy.apply_transform_Rx_plus_v( sh.R,sh.T );
		add_crystinfo_to_pose( posecopy, sh.sg );
		return posecopy;
	}

	void apply( Pose & pose) {
		numeric::xyzVector<Real> native_shift = center_pose_at_origin( pose );

		std::string infile = option[ crystdock::hits_in ];
		std::string outfile = option[ crystdock::hits_out ];
		core::Real cluster_cutoff = option[crystdock::cluster_cutoff]();

		if ( cluster_cutoff == 0 ) {
			TR << "Specify a non-zero cluster cutoff to continue!" << std::endl;
			exit(1);
		}

		utility::vector1< SpacegroupHit > in_hits, out_hits;

		std::ifstream in( infile.c_str() );
		while ( !in.eof() ) {
			std::string line;
			std::getline( in, line );
			if ( line.length() < 2 ) continue;
			in_hits.push_back( SpacegroupHit(line) );
		}

		TR << "Read " << in_hits.size() << " placements." << std::endl;

		for ( core::Size i=1; i<=in_hits.size(); ++i ) {
			TR << i << "/" << in_hits.size() << std::endl;
			core::pose::Pose pose1 = transform_pdb( pose, in_hits[i] );
			core::pose::Pose pose1c = transform_pdb( pose, in_hits[i] );
			bool keep = true;
			for ( core::Size j=1; j<=out_hits.size() && keep; ++j ) {
				core::pose::Pose pose2 = transform_pdb( pose, out_hits[j] );
				core::Real rms = protocols::cryst::crystRMSfast( pose1, pose2 );
				if ( rms < cluster_cutoff+3 ) {
					rms = protocols::cryst::crystRMS( pose1, pose2 );
				}
				if ( rms < cluster_cutoff ) keep = false;
			}

			if ( keep ) out_hits.push_back( in_hits[i] );
		}

		TR << "Have " << out_hits.size() << " after clustering." << std::endl;

		std::ofstream outs( outfile.c_str() );
		for ( core::Size i=1; i<=out_hits.size(); ++i ) {
			outs << out_hits[i].to_string() << std::endl;
		}
	}

	virtual std::string get_name() const { return "CrystSlideDock"; }
};

///////////////////////////////////////////////////////////////////////////////
class CrystSlideDock : public protocols::moves::Mover {
private:
	core::Real motif_init_cut_, motif_final_cut_, vdw_cut_;
	core::Real mininterface_;
	core::Real interfacedist_;
	Size dock_ncyc_;

public:
	CrystSlideDock() {
		interfacedist_ = option[crystdock::interfacedist];
		mininterface_ = option[ crystdock::mininterface ]();
		vdw_cut_ = option[crystdock::dock_vdw_cut]();
		motif_init_cut_ = option[crystdock::motif_initial_cut]();
		motif_final_cut_ = option[crystdock::motif_final_cut]();
		dock_ncyc_ = option[crystdock::dock_ncycle]();
	}

	Pose
	transform_pdb( Pose const &pose, SpacegroupHit sh ) {
		Pose posecopy = pose;
		posecopy.apply_transform_Rx_plus_v( sh.R,sh.T );
		add_crystinfo_to_pose( posecopy, sh.sg );
		return posecopy;
	}

	void
	get_docked_interface_score(
		Pose & pose,
		Real &score_vdw,
		Real &score_motif,
		utility::vector1<Real> &score_perinterface,
		Real &ncontacts,
		utility::vector1<Real> &ncontacts_perinterface

	) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::conformation::symmetry;

		core::scoring::motif::MotifHashManager & mman(*core::scoring::motif::MotifHashManager::get_instance());
		SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
		core::Size nsubunits = symm_info->subunits();

		core::scoring::symmetry::SymmetricScoreFunction sfvdw;
		sfvdw.set_weight( core::scoring::vdw, 1.0 );

		score_vdw=sfvdw(pose);  // even if wt is 0 do this to init energy graph

		Energies & energies( pose.energies() );
		EnergyGraph & energy_graph( energies.energy_graph() );

		core::Real all_bb_bb=0;

		score_motif = 0.0;
		score_perinterface.clear();
		score_perinterface.resize(nsubunits,0.0);

		ncontacts = 0.0;
		ncontacts_perinterface.clear();
		ncontacts_perinterface.resize(nsubunits,0.0);

		for ( Size i=1, i_end = pose.size(); i<= i_end; ++i ) {
			conformation::Residue const & rsd1( pose.residue( i ) );
			if ( !symm_info->bb_is_independent(rsd1.seqpos()) ) continue;  // not sure if this is necessary

			for ( utility::graph::Graph::EdgeListIter
					iru  = energy_graph.get_node(i)->upper_edge_list_begin(),
					irue = energy_graph.get_node(i)->upper_edge_list_end();
					iru != irue; ++iru ) {
				EnergyEdge & edge( static_cast< EnergyEdge & > (**iru) );

				Size const j( edge.get_second_node_ind() );
				conformation::Residue const & rsd2( pose.residue( j ) );

				if ( i==j ) continue;
				if ( !rsd1.is_protein() || !rsd2.is_protein() ) continue;
				if ( symm_info->bb_is_independent(rsd2.seqpos()) ) continue;

				// get reference frames
				numeric::xyzTransform<Real> const xbb1 = core::pose::motif::get_backbone_reference_frame(pose,i);
				numeric::xyzTransform<Real> const xbb2 = core::pose::motif::get_backbone_reference_frame(pose,j);

				// get score tables
				char const &ss1(pose.secstruct(i)), &ss2(pose.secstruct(j)), &aa1(rsd1.name1()), &aa2(rsd2.name1());

				core::Real wt=1.0;
				//if (ss1 == 'L' || ss2 == 'L') wt = 0.5; // downweight non-SS interactions by half

				core::scoring::motif::XformScoreCOP xs_bb_bb_fwd = mman.get_xform_score_BB_BB(ss1,ss2,aa1,aa2);
				core::scoring::motif::XformScoreCOP xs_bb_bb_rev = mman.get_xform_score_BB_BB(ss2,ss1,aa2,aa1);
				if ( xs_bb_bb_fwd || xs_bb_bb_rev ) {
					core::Real dist = xbb1.t.distance_squared(xbb2.t);
					if ( dist < ((interfacedist_)*(interfacedist_)) ) {
						core::Real this_bb_bb = 0;
						if ( xs_bb_bb_fwd ) this_bb_bb -= xs_bb_bb_fwd->score_of_bin( xbb1.inverse() * xbb2 );
						if ( xs_bb_bb_rev ) this_bb_bb -= xs_bb_bb_rev->score_of_bin( xbb2.inverse() * xbb1 );

						// fade to zero at long range (from maxdis/2 to maxdis)
						dist = std::sqrt(dist);
						core::Real halfdist = interfacedist_/2.0;
						core::Real fade1 = std::max(dist-halfdist,0.0)/halfdist;
						core::Real fade2 = (1-fade1*fade1);
						this_bb_bb *= fade2*fade2;

						// SS fade (if specified)
						this_bb_bb *= wt;

						all_bb_bb += this_bb_bb;
						score_perinterface[symm_info->subunit_index( j )] += this_bb_bb;

						// add a contact
						ncontacts += 1.0;
						ncontacts_perinterface[symm_info->subunit_index( j )] += 1.0;
					}
				}


			}
		}

		//score_motif = all_bb_bb;

		// return max interface
		score_motif=0.0;
		for ( int i=1; i<=(int)score_perinterface.size(); ++i ) {
			score_motif = std::min(score_motif, score_perinterface[i]);
		}
	}


	// get the score per interface
	core::Real
	get_weakest_score(
		Pose & pose,
		utility::vector1<core::kinematics::RT> const &rts,
		protocols::cryst::MakeLatticeMover &setup,  // should be const
		utility::vector1<Real> const &perint_scores,
		bool get_lowest ) {
		utility::vector1<SingleInterface> motif_interfaces;
		core::Real scale = get_lowest? -1.0 : 1.0;
		for ( int i=2; i<=(int)perint_scores.size(); ++i ) {
			numeric::xyzMatrix<core::Real> Ri;
			numeric::xyzVector<core::Real> Ti;
			setup.getRT( i, Ri, Ti );
			motif_interfaces.push_back( SingleInterface( Ri, Ti, scale*perint_scores[i] ) );
		}
		return scale*get_interface_score (pose, motif_interfaces, rts );
	}


	void apply( Pose & pose) {
		numeric::xyzVector<Real> native_shift = center_pose_at_origin( pose );

		std::string infile = option[ crystdock::hits_in ];
		utility::vector1< SpacegroupHit > in_hits;

		std::ifstream in( infile.c_str() );
		while ( !in.eof() ) {
			std::string line;
			std::getline( in, line );
			if ( line.length() < 2 ) continue;
			in_hits.push_back( SpacegroupHit(line) );
		}

		TR << "Read " << in_hits.size() << " placements." << std::endl;

		protocols::cryst::MakeLatticeMover setup;
		setup.contact_dist(16); // probably big enough for most applications

		protocols::cryst::DockLatticeMover docking;
		docking.set_trans_mag(0.2);
		docking.set_rot_mag(1.0);

		// save SCs
		protocols::simple_moves::ReturnSidechainMoverOP restore_sc (new protocols::simple_moves::ReturnSidechainMover( pose ));

		protocols::cryst::UpdateCrystInfo update_cryst1;
		protocols::simple_moves::SwitchResidueTypeSetMover to_cen("centroid");
		to_cen.apply(pose);

		// monomer bump score
		core::scoring::ScoreFunctionOP sf_vdw ( new core::scoring::symmetry::SymmetricScoreFunction() );
		sf_vdw->set_weight( core::scoring::vdw, 1.0 );
		core::Real monomer_bump=0;
		monomer_bump = 2*((*sf_vdw)( pose )); // * SymmConf.Symmetry_Info()->subunits();

		core::Size I_MIN=1, I_MAX=in_hits.size();

		if ( option[ crystdock::hits_to_dock ].user() ) {
			runtime_assert( option[ crystdock::hits_to_dock ]().size() >=2 );
			I_MIN = std::max(1,option[ crystdock::hits_to_dock ]()[1]);
			I_MAX = std::min(option[ crystdock::hits_to_dock ]()[2], (int)in_hits.size());
		}

		for ( core::Size i=I_MIN; i<=I_MAX; ++i ) {
			SpacegroupHit sh = in_hits[i];
			TR << i << ":  score = " << sh.score << std::endl;

			utility::vector1<core::kinematics::RT> const &rts = sh.sg.symmops();

			// Treat tags as file names so that we put the number before the extension.
			std::string base_name = protocols::jd2::current_input_tag();
			utility::vector1< std::string > temp_out_names= utility::split( base_name );
			utility::file::FileName out_name = utility::file::combine_names( temp_out_names );
			base_name = out_name.base();
			std::string outname = base_name+option[ out::suffix ]()+"_"+right_string_of( i, 8, '0' )+".pdb";

			core::pose::Pose pose_xform = transform_pdb( pose, sh );
			setup.apply( pose_xform );

			// get SS ... motif score uses it
			{ core::scoring::dssp::Dssp dssp( pose_xform ); dssp.insert_ss_into_pose( pose_xform ); }

			core::Real best_score=999;
			core::pose::Pose best_pose = pose_xform;
			docking.init( pose_xform );

			// various scores used in this stage
			core::Real score_motif=0, score_vdw=0, ncontacts=0;
			utility::vector1<Real> score_perinterface, ncontacts_per_interface;
			core::Real score_weakinterface, ncontacts_weakinterface;
			core::Real total_score;

			get_docked_interface_score( pose_xform, score_vdw, score_motif, score_perinterface, ncontacts, ncontacts_per_interface );
			score_weakinterface = get_weakest_score( pose_xform, rts, setup, score_perinterface, true );
			TR << "score before slide " << score_vdw << " / " << score_motif << " / " << score_weakinterface << std::endl;

			//docking.slide_lattice( pose_xform );
			docking.min_lattice( pose_xform );
			get_docked_interface_score( pose_xform, score_vdw, score_motif, score_perinterface, ncontacts, ncontacts_per_interface );

			TR << "score after slide " << score_vdw << " / " << score_motif << " / " << score_weakinterface << std::endl;

			// FILTER 0 -- slide failure
			if ( score_vdw-monomer_bump > vdw_cut_ ) {
				TR << "Fail! Slide move did not remove clashes." << std::endl;
				continue;
			}

			// if the slide moves the structure significantly, then we may be clashing with non-contacting subunits
			//   check that here
			update_cryst1.apply( pose_xform );
			setup.apply( pose_xform );

			// get SS again
			{ core::scoring::dssp::Dssp dssp( pose_xform ); dssp.insert_ss_into_pose( pose_xform ); }

			get_docked_interface_score( pose_xform, score_vdw, score_motif, score_perinterface, ncontacts, ncontacts_per_interface );
			score_weakinterface = get_weakest_score( pose_xform, rts, setup, score_perinterface, true );
			ncontacts_weakinterface = get_weakest_score( pose_xform, rts, setup, ncontacts_per_interface, false );

			// recheck contacts
			//if (ncontacts_weakinterface < mininterface_) {
			// TR << "Fail! Interface too small (" << ncontacts_weakinterface << ") after slide move!" << std::endl;
			// continue;
			//}

			// FILTER 1
			if ( score_vdw-monomer_bump > vdw_cut_ ) {
				TR << "Fail!  Regeneration after slide move introduced clashes." << std::endl;
				continue;
			}

			total_score = score_weakinterface;

			// FILTER 2
			if ( total_score > motif_init_cut_ ) {
				TR << "Fail!  Initial motif score: " << total_score << " > " << motif_init_cut_ << std::endl;
				continue;
			}

			// initialize for docking
			best_pose = pose_xform;
			best_score = total_score;

			// MAIN DOCKING LOOP
			// we now have to reinitialize docking mover since the symmetry was regenerated
			docking.init( pose_xform );

			for ( core::Size j=1; j<=dock_ncyc_; ++j ) {
				docking.perturb_trial( pose_xform );

				get_docked_interface_score( pose_xform, score_vdw, score_motif, score_perinterface, ncontacts, ncontacts_per_interface );
				score_weakinterface = get_weakest_score( pose_xform, rts, setup, score_perinterface, true );
				ncontacts_weakinterface = get_weakest_score( pose_xform, rts, setup, ncontacts_per_interface, false );

				///
				total_score = score_weakinterface;
				///

				if ( total_score<best_score && score_vdw-monomer_bump < vdw_cut_ ) {
					TR << "New best " << total_score << std::endl;
					best_score = total_score;
					best_pose = pose_xform;
				} else {
					pose_xform = best_pose;
				}
			}


			// FINAL REPORT
			get_docked_interface_score( best_pose, score_vdw, score_motif, score_perinterface, ncontacts, ncontacts_per_interface );
			score_weakinterface = get_weakest_score( pose_xform, rts, setup, score_perinterface, true );
			ncontacts_weakinterface = get_weakest_score( pose_xform, rts, setup, ncontacts_per_interface, false );

			///
			total_score = score_weakinterface;
			///

			if ( total_score > motif_final_cut_ ) {
				TR << "Failed final motif energy filter: " << total_score << " > " << motif_final_cut_ << std::endl;
				continue;
			}

			core::io::RemarkInfo remark;
			remark.num = 1;
			remark.value = "score = "
				+ utility::to_string(ncontacts) + " | "
				+ utility::to_string(ncontacts_weakinterface) + " | "
				+ utility::to_string(score_motif) + " | "
				+ utility::to_string(score_weakinterface) + " | "
				+ utility::to_string(total_score);
			best_pose.pdb_info()->remarks().push_back( remark );
			remark.num = 2;

			for ( int k=2; k<=(int)score_perinterface.size(); ++k ) {
				core::io::RemarkInfo remark;
				remark.value = "score [" + utility::to_string(k) + "] = "
					+ utility::to_string(score_perinterface[k]) + " | "
					+ utility::to_string(ncontacts_per_interface[k]);
				best_pose.pdb_info()->remarks().push_back( remark );
			}

			update_cryst1.apply(best_pose);

			if ( option[in::file::native].user() ) {
				CrystRMS eval;
				eval.apply(best_pose);
			}

			best_pose.dump_pdb( outname );
		}
	}

	virtual std::string get_name() const { return "CrystSlideDock"; }
};

///////////////////////////////////////////////////////////////////////////////
class CrystFFTDock : public protocols::moves::Mover {
private:
	core::Real ca_clashdist_, cb_clashdist_, n_clashdist_, c_clashdist_, o_clashdist_;
	core::Real interfacedist_;
	core::Real cluster_cutoff_;
	numeric::xyzVector<int> grid_, oversamplegrid_;
	numeric::xyzMatrix<Real> i2c_, c2i_;
	protocols::cryst::Spacegroup sg_;

	// parameters from options
	core::Real maxclash_, mininterface_, mininterface_sum_, trans_step_, rot_step_;
	Size nmodels_, rotnum_;
	bool ss_only_, eval_native_;
	bool debug_, debug_exact_;

public:
	CrystFFTDock() {
		// to do: make these options
		n_clashdist_  =  1.40;   // ros VDW=1.75
		ca_clashdist_ =  2.00;   // ros VDW=2.0
		c_clashdist_  =  2.00;   // ros VDW=2.0
		o_clashdist_  =  1.30;   // ros VDW=1.55
		cb_clashdist_ =  1.50;   // ros VDW=2.0    .. make this a bit smaller to encourage better packing

		//
		interfacedist_ = option[crystdock::interfacedist];

		// filters
		maxclash_ = option[ crystdock::maxclash ]();
		mininterface_ = option[ crystdock::mininterface ]();
		mininterface_sum_ = option[ crystdock::mininterfacesum ]();

		// options: scoring
		ss_only_ = option[ crystdock::ssonly ]();

		// options: sampling
		cluster_cutoff_ = option[crystdock::cluster_cutoff]();
		trans_step_ = option[ crystdock::trans_step ]();
		rot_step_ = option[ crystdock::rot_step ]();
		nmodels_ = option[ crystdock::nmodels ]();

		// debugging options
		rotnum_ = option[ crystdock::rotnum ]();
		debug_ = option[ crystdock::debug ]();
		debug_exact_ = option[ crystdock::debug_exact ]();
		eval_native_ = option[ crystdock::eval_native ]();
	}

	virtual std::string get_name() const { return "CrystFFTDock"; }

	// write density grids in MRC format for debugging
	void
	writeMRC(FArray3D<Real> density, std::string mapfilename, bool is_oversampled=false, bool fftshift=false);

	// build occupancy and interaction masks from pose
	// sample in bounding box around molecule (rather than unit cell)
	void setup_maps( Pose & pose, FArray3D<Real> &rho_ca, FArray3D<Real> &rho_cb, Real trans_step);

	// resample maps subject to rotation
	// get self clashes and self rotations
	core::Real resample_maps_and_get_self(
		FArray3D<Real> const &rho_ca, FArray3D<Real> const &rho_cb,
		numeric::xyzMatrix<Real> R, protocols::cryst::Spacegroup const &sg,
		FArray3D<Real> &r_rho_ca, FArray3D<Real> &r_rho_cb,
		utility::vector1<SingleInterface> &p1_interface_map );

	// aplpy xform and dump PDB
	void dump_transformed_pdb( Pose pose, InterfaceHit ih, numeric::UniformRotationSampler const &urs, std::string outname );


	Pose transform_pdb( Pose const &pose, InterfaceHit ih, numeric::UniformRotationSampler const &urs );

	// nearest-neighbor interpolation subject to grid-space transform
	void transform_map(
		FArray3D<Real> const &rho,
		numeric::xyzMatrix<Real> S, numeric::xyzVector<Real> T,
		FArray3D<Real> &Srho);

	// same as previous function, assumes offset of 0
	void transform_map_offset0(
		FArray3D<Real> const &rho,
		numeric::xyzMatrix<Real> S,
		FArray3D<Real> &Srho);

	// get maximum radius from a pose
	Real get_radius_of_pose( Pose & pose );

	//Calculate transformed distance
	Real get_transform_distance (InterfaceHit ih_vec, InterfaceHit ih_vec_clustered , numeric::UniformRotationSampler const &urs, Real radius);

	// wrapper does a series of 2d ffts
	//    this could be an option to the low-level fft code to possibly save time
	void
	fft2dslice( FArray3D<Real> const &rho, FArray3D< std::complex<Real> > &Frho, int axisSlice );

	// wrapper does a series of 2d iffts
	void
	ifft2dslice( FArray3D< std::complex<Real> > const &Frho, FArray3D<Real> &rho, int axisSlice );

	// project scores along some axis, not necessarily orthogonal to unit cell
	void
	project_along_axis( FArray3D<Real> &rho, numeric::xyzVector<Real> axis );

	// the main convolution operation used for both CA and CB maps
	void
	do_convolution( FArray3D<Real> const &rho, FArray3D<Real> const &Srho, numeric::xyzMatrix<Real> S, FArray3D<Real> &conv_out);

	// get the score per interface
	void
	get_interfaces_allatom(
		Pose pose,
		utility::vector1<core::kinematics::RT> const &rts,
		numeric::xyzMatrix<Real> R,
		numeric::xyzVector<Real> xyz_grid,
		utility::vector1<SingleInterface> &allInterfaces );

	// get the exact clash score (debugging/eval_native only)
	core::Real
	get_clash_score_exact(
		numeric::xyzVector<int> xyz_grid,
		numeric::xyzMatrix<Real> R,
		numeric::xyzVector<Real> T,
		FArray3D<Real> const &r_rho_ca);

	void apply( Pose & pose);
};

// write density grids in MRC --  debugging only for now
void
CrystFFTDock::writeMRC(FArray3D<Real> density, std::string mapfilename, bool is_oversampled /*=false*/, bool fftshift /*=false*/) {
	const int CCP4HDSIZE = 1024;  // size of CCP4/MRC header
	std::fstream outx( (mapfilename).c_str() , std::ios::binary | std::ios::out );

	float buff_f, buff_vf[3];
	int buff_i, buff_vi[3], symBytes = 0;

	if ( !outx ) {
		TR << "Error opening " << mapfilename << " for writing." << std::endl;
		return;
	}

	// extent
	buff_vi[0] = density.u1(); buff_vi[1] = density.u2(); buff_vi[2] = density.u3();
	outx.write(reinterpret_cast <char*>(buff_vi), sizeof(int)*3);

	// mode
	buff_i = 2;
	outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));

	// origin
	int ori_int[3] = {0,0,0};
	outx.write(reinterpret_cast <char*>(ori_int), sizeof(int)*3);

	// grid
	int grid[3] = {density.u1(),density.u2(),density.u3()};
	outx.write(reinterpret_cast <char*>(&grid[0]), sizeof(int)*3);

	// cell params
	Real Aeff=sg_.A(), Beff=sg_.B(), Ceff=sg_.C();
	if ( is_oversampled ) {
		Aeff *= ((Real)oversamplegrid_[0]) / ((Real)grid_[0]);
		Beff *= ((Real)oversamplegrid_[1]) / ((Real)grid_[1]);
		Ceff *= ((Real)oversamplegrid_[2]) / ((Real)grid_[2]);
	}

	float cellDimensions[3] = {(float)Aeff,(float)Beff,(float)Ceff};
	float cellAngles[3] = {(float)sg_.alpha(),(float)sg_.beta(),(float)sg_.gamma()};
	outx.write(reinterpret_cast <char*>(&cellDimensions), sizeof(float)*3);
	outx.write(reinterpret_cast <char*>(&cellAngles), sizeof(float)*3);

	// crs2xyz
	buff_vi[0] = 1; buff_vi[1] = 2; buff_vi[2] = 3;
	outx.write(reinterpret_cast <char*>(buff_vi), sizeof(int)*3);

	// min, max, mean dens
	buff_vf[0] = -100.0; buff_vf[1] = 100.0; buff_vf[2] = 0.0;
	outx.write(reinterpret_cast <char*>(buff_vf), sizeof(float)*3);

	// 4 bytes junk
	buff_i = 0;
	outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));

	// symmops
	outx.write(reinterpret_cast <char*>(&symBytes), sizeof(int));

	// 104 bytes junk
	buff_i = 0;
	for ( int i=0; i<25; ++i ) {
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
	}

	// alt origin (MRC)
	float ori_float[3]={0,0,0};
	outx.write(reinterpret_cast <char*>(ori_float), sizeof(float)*3);

	// Write "MAP" at byte 208, indicating a CCP4 file.
	char buff_s[80]; strcpy(buff_s, "MAP DD");
	outx.write(reinterpret_cast <char*>(buff_s), 8);

	// fill remainder of head with junk
	int nJunkWords = (CCP4HDSIZE - 216) /4;
	buff_i = 0;
	for ( int i=0; i<nJunkWords; ++i ) {
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
	}

	// data
	if ( fftshift ) {
		int coord[3], rcoord[3];
		for ( coord[2] = 1; coord[2] <= density.u3(); coord[2]++ ) {
			for ( coord[1] = 1; coord[1] <= density.u2(); coord[1]++ ) {
				for ( coord[0] = 1; coord[0] <= density.u1(); coord[0]++ ) {
					rcoord[0] = pos_mod( coord[0]-1+density.u1()/2, density.u1())+1;
					rcoord[1] = pos_mod( coord[1]-1+density.u2()/2, density.u2())+1;
					rcoord[2] = pos_mod( coord[2]-1+density.u3()/2, density.u3())+1;
					buff_f = (float) density(rcoord[0],rcoord[1],rcoord[2]);
					outx.write(reinterpret_cast <char*>(&buff_f), sizeof(float));
				}
			}
		}
	} else {
		int coord[3];
		for ( coord[2] = 1; coord[2] <= density.u3(); coord[2]++ ) {
			for ( coord[1] = 1; coord[1] <= density.u2(); coord[1]++ ) {
				for ( coord[0] = 1; coord[0] <= density.u1(); coord[0]++ ) {
					buff_f = (float) density(coord[0],coord[1],coord[2]);
					outx.write(reinterpret_cast <char*>(&buff_f), sizeof(float));
				}
			}
		}
	}
}


// build occupancy and interaction masks from pose
// since we will be sampling rotations from this map, we sample over a larger volume than the unit cell
//  ASSUMES POSE IS CENTERED
void
CrystFFTDock::setup_maps( Pose & pose, FArray3D<Real> &rho_ca, FArray3D<Real> &rho_cb, Real trans_step) {
	Real ATOM_MASK_PADDING = 2.0;
	Real UNIT_CELL_PADDING = 6.0;  // IN GRID POINTS!
	Real sigwidth=option[crystdock::sigwidth];
	Real interface_sigwidth=option[crystdock::interface_sigwidth];
	Size minmult = sg_.minmult();

	// find true grid
	grid_[0] = minmult*(Size)std::floor(sg_.A()/(minmult*trans_step) + 0.5);
	grid_[1] = minmult*(Size)std::floor(sg_.B()/(minmult*trans_step) + 0.5);
	grid_[2] = minmult*(Size)std::floor(sg_.C()/(minmult*trans_step) + 0.5);

	numeric::xyzMatrix<Real> i2f, f2i;
	f2i = numeric::xyzMatrix<Real>::rows( (Real)grid_[0],0,0,  0,(Real)grid_[1],0,  0,0,(Real)grid_[2] );
	i2f = numeric::xyzMatrix<Real>::rows( 1.0/((Real)grid_[0]),0,0,  0,1.0/((Real)grid_[1]),0,  0,0,1.0/((Real)grid_[2]) );

	i2c_ = sg_.f2c()*i2f;
	c2i_ = f2i*sg_.c2f();

	// find oversampled grid
	oversamplegrid_ = numeric::xyzVector<int>(0,0,0);
	for ( int i=1 ; i<=(int)pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		for ( int j=1; j<=4; ++j ) {
			numeric::xyzVector< core::Real> xyz_j = pose.residue(i).atom(j).xyz();
			numeric::xyzVector< core::Real> atm_idx = c2i_*xyz_j;
			for ( int k=0; k<3; ++k ) {
				oversamplegrid_[k] = std::max(oversamplegrid_[k], 2*(int)std::floor( (std::abs(atm_idx[k])+UNIT_CELL_PADDING) )+1 );
			}
		}
	}
	TR << "Base grid = [" <<  oversamplegrid_[0] << " , " <<  oversamplegrid_[1]<< " , " << oversamplegrid_[2] << "]" << std::endl;

	rho_ca.dimension( oversamplegrid_[0], oversamplegrid_[1], oversamplegrid_[2] ); rho_ca = 1;
	rho_cb.dimension( oversamplegrid_[0], oversamplegrid_[1], oversamplegrid_[2] ); rho_cb = 1;

	// loop over bb heavyatoms
	for ( int i=1 ; i<=(int)pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;

		for ( int j=1; j<=4; ++j ) {
			core::Real clashdist=0.0;
			if ( j==1 ) clashdist = n_clashdist_;
			if ( j==2 ) clashdist = ca_clashdist_;
			if ( j==3 ) clashdist = c_clashdist_;
			if ( j==4 ) clashdist = o_clashdist_;

			numeric::xyzVector< core::Real> xyz_j = pose.residue(i).atom(j).xyz();
			numeric::xyzVector< core::Real> atm_idx = c2i_*xyz_j;
			numeric::xyzVector< core::Real> atm_j, del_ij;

			for ( int z=1; z<=oversamplegrid_[2]; ++z ) {
				atm_j[2] = z-1;
				del_ij[2] = min_mod(atm_idx[2] - atm_j[2], (Real)oversamplegrid_[2]);
				del_ij[0] = del_ij[1] = 0.0;
				if ( (i2c_*del_ij).length_squared() > (clashdist+ATOM_MASK_PADDING)*(clashdist+ATOM_MASK_PADDING) ) continue;
				for ( int y=1; y<=oversamplegrid_[1]; ++y ) {
					atm_j[1] = y-1;
					del_ij[1] = min_mod(atm_idx[1] - atm_j[1], (Real)oversamplegrid_[1]);
					del_ij[0] = 0.0;
					if ( (i2c_*del_ij).length_squared() > (clashdist+ATOM_MASK_PADDING)*(clashdist+ATOM_MASK_PADDING) ) continue;
					for ( int x=1; x<=oversamplegrid_[0]; ++x ) {
						atm_j[0] = x-1;
						del_ij[0] = min_mod(atm_idx[0] - atm_j[0], (Real)oversamplegrid_[0]);
						numeric::xyzVector< core::Real > cart_del_ij = i2c_*del_ij;
						core::Real d2 = (cart_del_ij).length_squared();
						if ( d2 <= (clashdist+ATOM_MASK_PADDING)*(clashdist+ATOM_MASK_PADDING) ) {
							core::Real doff = sqrt(d2) - clashdist;
							core::Real sig = 1 / ( 1 + exp ( -sigwidth*doff ) );   //
							rho_ca(x,y,z) *= sig;
						}
					}
				}
			}
		}
	}

	// loop over CBs
	for ( int i=1 ; i<=(int)pose.size(); ++i ) {
		if ( pose.residue(i).aa() == core::chemical::aa_gly ) continue;
		numeric::xyzVector< core::Real> CB = pose.residue(i).atom(5).xyz();
		numeric::xyzVector< core::Real> atm_idx = c2i_*CB;
		numeric::xyzVector< core::Real> atm_j, del_ij;

		for ( int z=1; z<=oversamplegrid_[2]; ++z ) {
			atm_j[2] = z-1;
			del_ij[2] = min_mod(atm_idx[2] - atm_j[2], (Real)oversamplegrid_[2]);
			del_ij[0] = del_ij[1] = 0.0;
			if ( (i2c_*del_ij).length_squared() > (interfacedist_+ATOM_MASK_PADDING)*(interfacedist_+ATOM_MASK_PADDING) ) continue;
			for ( int y=1; y<=oversamplegrid_[1]; ++y ) {
				atm_j[1] = y-1;
				del_ij[1] = min_mod(atm_idx[1] - atm_j[1], (Real)oversamplegrid_[1]);
				del_ij[0] = 0.0;
				if ( (i2c_*del_ij).length_squared() > (interfacedist_+ATOM_MASK_PADDING)*(interfacedist_+ATOM_MASK_PADDING) ) continue;
				for ( int x=1; x<=oversamplegrid_[0]; ++x ) {
					atm_j[0] = x-1;
					del_ij[0] = min_mod(atm_idx[0] - atm_j[0], (Real)oversamplegrid_[0]);
					numeric::xyzVector< core::Real > cart_del_ij = i2c_*del_ij;
					core::Real d2 = (cart_del_ij).length_squared();
					if ( d2 <= (interfacedist_+ATOM_MASK_PADDING)*(interfacedist_+ATOM_MASK_PADDING) ) {
						core::Real doff = sqrt(d2) - cb_clashdist_;
						core::Real sig = 1 / ( 1 + exp ( -sigwidth*doff ) );   // '6' gives sigmoid dropoff
						rho_ca(x,y,z) *= sig;

						if ( !ss_only_ || pose.secstruct(i)!='L' ) {
							doff = sqrt(d2) - interfacedist_;
							sig = 1 / ( 1 + exp ( -interface_sigwidth*doff ) );
							rho_cb(x,y,z) *= sig;
						}
					}
				}
			}
		}
	}

	// factor CA mask out of CB mask; invert both
	for ( int i=0 ; i<oversamplegrid_[0]*oversamplegrid_[1]*oversamplegrid_[2]; ++i ) {
		rho_cb[i] = (1-rho_cb[i])*rho_ca[i];
		rho_ca[i] = (1-rho_ca[i]);
	}
}


// resample maps subject to rotation
// get self clashes and self rotations
core::Real
CrystFFTDock::resample_maps_and_get_self(
	FArray3D<Real> const &rho_ca, FArray3D<Real> const &rho_cb,
	numeric::xyzMatrix<Real> R, protocols::cryst::Spacegroup const &sg,
	FArray3D<Real> &r_rho_ca, FArray3D<Real> &r_rho_cb,
	utility::vector1<SingleInterface> &p1_interface_map )
{
	r_rho_ca.dimension( grid_[0], grid_[1], grid_[2] ); r_rho_ca=0;
	r_rho_cb.dimension( grid_[0], grid_[1], grid_[2] ); r_rho_cb=0;

	FArray3D<Real> r_rho_ca_base = r_rho_ca;
	FArray3D<Real> r_rho_cb_base = r_rho_cb;

	numeric::xyzMatrix<Real> Ri = numeric::inverse(R);
	numeric::xyzMatrix<Real> Rigridspace = c2i_*Ri*i2c_;
	numeric::xyzMatrix<Real> Rgridspace = c2i_*R*i2c_;

	// rotate "oversmaple grid and find the boundaries
	Real xmax=0, ymax=0, zmax=0;
	{
		numeric::xyzVector<Real> boundbox;
		boundbox = Rgridspace * numeric::xyzVector<Real>(0, 0, oversamplegrid_[2]);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(0, oversamplegrid_[1], 0);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(oversamplegrid_[0], 0, 0);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(0, oversamplegrid_[1], oversamplegrid_[2]);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(oversamplegrid_[0], 0, oversamplegrid_[2]);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(oversamplegrid_[0], oversamplegrid_[1], 0);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(oversamplegrid_[0], oversamplegrid_[1], oversamplegrid_[2]);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
	}

	int AMAX = (int)std::ceil( 0.5*(xmax/grid_[0]-1) );
	int BMAX = (int)std::ceil( 0.5*(ymax/grid_[1]-1) );
	int CMAX = (int)std::ceil( 0.5*(zmax/grid_[2]-1) );

	// calculate base transformation
	for ( int z=1; z<=grid_[2]; ++z ) {
		for ( int y=1; y<=grid_[1]; ++y ) {
			for ( int x=1; x<=grid_[0]; ++x ) {
				int cx=x-1,cy=y-1,cz=z-1;
				if ( cx>grid_[0]/2 ) cx -= grid_[0];
				if ( cy>grid_[1]/2 ) cy -= grid_[1];
				if ( cz>grid_[2]/2 ) cz -= grid_[2];

				numeric::xyzVector<Real> rx = Rigridspace*numeric::xyzVector<Real>(cx,cy,cz);

				r_rho_ca_base(x,y,z) = interp_linear( rho_ca, rx );
				r_rho_cb_base(x,y,z) = interp_linear( rho_cb, rx );
			}
		}
	}

	r_rho_ca = r_rho_ca_base;
	r_rho_cb = r_rho_cb_base;

	Real ca_overlap = 0;
	Real Npoints = grid_[0] * grid_[1] * grid_[2];
	Real voxel_volume = sg.volume() / Npoints;

	// offset transformations
	for ( int a=-AMAX; a<=AMAX; ++a ) {
		for ( int b=-BMAX; b<=BMAX; ++b ) {
			for ( int c=-CMAX; c<=CMAX; ++c ) {
				if ( a==0 && b==0 && c==0 ) continue;

				Real cb_overlap_abc = 0;

				for ( int z=1; z<=grid_[2]; ++z ) {
					for ( int y=1; y<=grid_[1]; ++y ) {
						for ( int x=1; x<=grid_[0]; ++x ) {
							int cx=x-1,cy=y-1,cz=z-1;
							if ( cx>grid_[0]/2 ) cx -= grid_[0];
							if ( cy>grid_[1]/2 ) cy -= grid_[1];
							if ( cz>grid_[2]/2 ) cz -= grid_[2];

							cx += a*grid_[0]; cy += b*grid_[1]; cz += c*grid_[2];

							numeric::xyzVector<Real> rx = Rigridspace*numeric::xyzVector<Real>(cx,cy,cz);

							Real rho_ca_rx = interp_linear( rho_ca, rx );
							Real rho_cb_rx = interp_linear( rho_cb, rx );

							// compute exact overlap
							ca_overlap     += rho_ca_rx*r_rho_ca_base(x,y,z);
							cb_overlap_abc += rho_cb_rx*r_rho_cb_base(x,y,z);

							// add to unit cell (will be used if this is non-overlapping)
							r_rho_ca(x,y,z) += rho_ca_rx;
							r_rho_cb(x,y,z) += rho_cb_rx;
						}
					}
				}

				// add interface if large enough
				cb_overlap_abc *= voxel_volume;
				if ( cb_overlap_abc > mininterface_ ) {
					p1_interface_map.push_back(
						SingleInterface(
						numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1),
						numeric::xyzVector<core::Real>(a,b,c),
						cb_overlap_abc
						)
					);
				}
			}
		}
	}

	return ca_overlap;
}


Real
CrystFFTDock::get_radius_of_pose( Pose & pose ) {
	Real radius = pose.residue(1).atom(2).xyz().length();
	for ( int i=2; i<=(int)pose.size(); ++i ) {
		Real r=pose.residue(i).atom(2).xyz().length();
		if ( radius < r ) radius=r;
	}
	return radius;
}

Real
CrystFFTDock::get_transform_distance (InterfaceHit ih_vec, InterfaceHit ih_vec_clustered , numeric::UniformRotationSampler const &urs, Real radius){
	Real x1=ih_vec.x, y1=ih_vec.y, z1=ih_vec.z, x2=ih_vec_clustered.x, y2=ih_vec_clustered.y, z2=ih_vec_clustered.z;
	Size rot1=ih_vec.rot_index, rot2=ih_vec_clustered.rot_index;

	numeric::xyzMatrix<Real> R1, R2;
	urs.get(rot1,R1);
	urs.get(rot2,R2);
	numeric::xyzVector< core::Real > vec1;
	numeric::xyzVector< core::Real > vec2;
	vec1[0]=x1;vec1[1]=y1;vec1[2]=z1;
	vec2[0]=x2;vec2[1]=y2;vec2[2]=z2;
	Real dist=(vec1-vec2).length();
	R2=numeric::inverse(R2);
	Real angle=DEG2RAD*numeric::urs_R2ang(R1*R2);
	Real transform_dist=angle*radius+dist;

	return transform_dist;
}


Pose
CrystFFTDock::transform_pdb( Pose const &pose, InterfaceHit ih, numeric::UniformRotationSampler const &urs ) {
	numeric::xyzMatrix<Real> R;
	urs.get( ih.rot_index, R );
	numeric::xyzVector<Real> T (ih.x, ih.y, ih.z);

	Pose posecopy = pose;
	posecopy.apply_transform_Rx_plus_v( R,T );
	return posecopy;
}

void
CrystFFTDock::dump_transformed_pdb( Pose pose, InterfaceHit ih, numeric::UniformRotationSampler const &urs, std::string outname ) {
	numeric::xyzMatrix<Real> R;
	urs.get( ih.rot_index, R );
	numeric::xyzVector<Real> T (ih.x, ih.y, ih.z);

	// add score to header
	core::io::RemarkInfo remark;
	std::ostringstream oss;
	oss << "  score = " << ih.score;
	remark.num = 1; remark.value = oss.str();
	pose.pdb_info()->remarks().push_back( remark );

	pose.apply_transform_Rx_plus_v( R,T );
	pose.dump_pdb( outname );
}

// nearest-neighbor interpolation subject to grid-space transform
void
CrystFFTDock::transform_map(
	FArray3D<Real> const &rho,
	numeric::xyzMatrix<Real> S, numeric::xyzVector<Real> T,
	FArray3D<Real> &Srho) {
	Srho.dimension( rho.u1(), rho.u2(), rho.u3() );
	for ( int z=1; z<=rho.u3(); ++z ) {
		for ( int y=1; y<=rho.u2(); ++y ) {
			for ( int x=1; x<=rho.u1(); ++x ) {
				int cx=x-1,cy=y-1,cz=z-1;
				int rx = 1 + pos_mod( (int)std::floor( (S.xx()*cx) + (S.xy()*cy) + (S.xz()*cz) + (T[0]*rho.u1()) + 0.5 ) , rho.u1());
				int ry = 1 + pos_mod( (int)std::floor( (S.yx()*cx) + (S.yy()*cy) + (S.yz()*cz) + (T[1]*rho.u2()) + 0.5 ) , rho.u2());
				int rz = 1 + pos_mod( (int)std::floor( (S.zx()*cx) + (S.zy()*cy) + (S.zz()*cz) + (T[2]*rho.u3()) + 0.5 ) , rho.u3());
				Real rho_sx = rho(rx,ry,rz);
				Srho(x,y,z) = rho_sx;
			}
		}
	}
}

// same as previous function, applies an offset of 0
void
CrystFFTDock::transform_map_offset0(
	FArray3D<Real> const &rho,
	numeric::xyzMatrix<Real> S,
	FArray3D<Real> &Srho) {
	Srho.dimension( rho.u1(), rho.u2(), rho.u3() );
	for ( int z=1; z<=rho.u3(); ++z ) {
		for ( int y=1; y<=rho.u2(); ++y ) {
			for ( int x=1; x<=rho.u1(); ++x ) {
				int cx=x-1,cy=y-1,cz=z-1;
				int rx = 1 + pos_mod( (int)std::floor( (S.xx()*cx) + (S.xy()*cy) + (S.xz()*cz) ) , rho.u1());
				int ry = 1 + pos_mod( (int)std::floor( (S.yx()*cx) + (S.yy()*cy) + (S.yz()*cz) ) , rho.u2());
				int rz = 1 + pos_mod( (int)std::floor( (S.zx()*cx) + (S.zy()*cy) + (S.zz()*cz) ) , rho.u3());
				Real rho_sx = rho(rx,ry,rz);
				Srho(x,y,z) = rho_sx;
			}
		}
	}
}

// wrapper does a series of 2d ffts
//    this could be an option to the low-level fft code to possibly save time
void
CrystFFTDock::fft2dslice( FArray3D<Real> const &rho, FArray3D< std::complex<Real> > &Frho, int axisSlice ) {
	FArray2D<Real> rhoSlice;
	FArray2D< std::complex<Real> > FrhoSlice;
	int xi = rho.u1(), yi = rho.u2(), zi = rho.u3();
	Frho.dimension(xi,yi,zi);

	if ( axisSlice == 1 ) {
		rhoSlice.dimension(yi,zi);
		for ( int ii=1; ii<=xi; ii++ ) {
			for ( int jj=1; jj<=yi; jj++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					rhoSlice(jj,kk) = rho(ii,jj,kk);
				}
			}
			numeric::fourier::fft2(rhoSlice, FrhoSlice);
			for ( int jj=1; jj<=yi; jj++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					Frho(ii,jj,kk) = FrhoSlice(jj,kk);
				}
			}
		}
	} else if ( axisSlice == 2 ) {
		rhoSlice.dimension(xi,zi);
		for ( int jj=1; jj<=yi; jj++ ) {
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					rhoSlice(ii,kk) = rho(ii,jj,kk);
				}
			}
			numeric::fourier::fft2(rhoSlice, FrhoSlice);
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					Frho(ii,jj,kk) = FrhoSlice(ii,kk);
				}
			}
		}
	} else if ( axisSlice == 3 ) {
		rhoSlice.dimension(xi,yi);
		for ( int kk=1; kk<=zi; kk++ ) {
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int jj=1; jj<=yi; jj++ ) {
					rhoSlice(ii,jj) = rho(ii,jj,kk);
				}
			}
			numeric::fourier::fft2(rhoSlice, FrhoSlice);
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int jj=1; jj<=yi; jj++ ) {
					Frho(ii,jj,kk) = FrhoSlice(ii,jj);
				}
			}
		}
	} else {
		utility_exit_with_message( "ERROR! Bad axis specified!");
	}
}

// wrapper does a series of 2d iffts
void
CrystFFTDock::ifft2dslice( FArray3D< std::complex<Real> > const &Frho, FArray3D<Real> &rho, int axisSlice ) {
	FArray2D<Real> rhoSlice;
	FArray2D< std::complex<Real> > FrhoSlice;
	int xi = Frho.u1(), yi = Frho.u2(), zi = Frho.u3();
	rho.dimension(xi,yi,zi);
	rho = 0;

	if ( axisSlice == 1 ) {
		FrhoSlice.dimension(yi,zi);
		for ( int ii=1; ii<=xi; ii++ ) {
			for ( int jj=1; jj<=yi; jj++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					FrhoSlice(jj,kk) = Frho(ii,jj,kk);
				}
			}
			numeric::fourier::ifft2(FrhoSlice, rhoSlice);
			for ( int jj=1; jj<=yi; jj++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					rho(ii,jj,kk) = rhoSlice(jj,kk);
				}
			}
		}
	} else if ( axisSlice == 2 ) {
		FrhoSlice.dimension(xi,zi);
		for ( int jj=1; jj<=yi; jj++ ) {
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					FrhoSlice(ii,kk) = Frho(ii,jj,kk);
				}
			}
			numeric::fourier::ifft2(FrhoSlice, rhoSlice);
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					rho(ii,jj,kk) = rhoSlice(ii,kk);
				}
			}
		}
	} else if ( axisSlice == 3 ) {
		FrhoSlice.dimension(xi,yi);
		for ( int kk=1; kk<=zi; kk++ ) {
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int jj=1; jj<=yi; jj++ ) {
					FrhoSlice(ii,jj) = Frho(ii,jj,kk);
				}
			}
			numeric::fourier::ifft2(FrhoSlice, rhoSlice);
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int jj=1; jj<=yi; jj++ ) {
					rho(ii,jj,kk) = rhoSlice(ii,jj);
				}
			}
		}
	} else {
		utility_exit_with_message( "ERROR! Bad axis specified!");
	}
}

void
CrystFFTDock::project_along_axis( FArray3D<Real> &rho, numeric::xyzVector<Real> axis ) {
	FArray3D<Real> rho_input = rho;

	Real scalefact = std::max(std::fabs(axis[0]), std::max(std::fabs(axis[1]),std::fabs(axis[2])));
	axis /= scalefact;

	int ngrid = (std::fabs(axis[0]>0.999))? grid_[0] : ((std::fabs(axis[1])>0.999)? grid_[1] : grid_[2]);

	// for pp = 1:D-1; mfft = mfft+circshift( mfft_orig, [pp*slide(1,1),pp*slide(2,2),pp*slide(3,3)]);
	for ( int P=1; P<ngrid; ++P ) {  // one less than full cycle
		for ( int z=1; z<=grid_[2]; ++z ) {
			for ( int y=1; y<=grid_[1]; ++y ) {
				for ( int x=1; x<=grid_[0]; ++x ) {
					int xt = pos_mod((int)floor( x+P*axis[0]-0.5 ), grid_[0]) + 1;
					int yt = pos_mod((int)floor( y+P*axis[1]-0.5 ), grid_[1]) + 1;
					int zt = pos_mod((int)floor( z+P*axis[2]-0.5 ), grid_[2]) + 1;
					rho(x,y,z) += rho_input(xt,yt,zt);
				}
			}
		}
	}
}

// the main convolution operation used for both CA and CB maps
void
CrystFFTDock::do_convolution( FArray3D<Real> const &rho, FArray3D<Real> const &Srho, numeric::xyzMatrix<Real> S, FArray3D<Real> &conv_out) {
	FArray3D<Real> working_s, working_trans, working_strans;
	FArray3D< std::complex<Real> > Fworking_trans, Fworking_strans;
	Size Npoints = rho.u1()*rho.u2()*rho.u3();

	// find out what plane we are slicing along
	numeric::urs_Quat Sinv_Q(S);

	// many different possibilities explain the convoluted logic below
	//    - this may return (1,0,0), (sqrt(1/2),sqrt(1/2),0) or (1/2,1/2,1/2) for X, XY, and XYZ symm axes
	//    - also in some spacegps( e.g. P3112[5] ) it may be (1,1/4,0) -- since it's not a rotation but a skew
	//    - in others ( e.g. P6122[8] ) it may also be (1,1/4,0) but we need to pre-skew
	//    - in these cases, we want to process as 'x' rotation with the skew handled in the projection
	int slice_x = (Sinv_Q.x_>=0.5) ? 1 : ((Sinv_Q.x_<=-0.5) ? -1:0);
	int slice_y = (Sinv_Q.y_>=0.5) ? 1 : ((Sinv_Q.y_<=-0.5) ? -1:0);
	int slice_z = (Sinv_Q.z_>=0.5) ? 1 : ((Sinv_Q.z_<=-0.5) ? -1:0);

	numeric::xyzVector<Real> rotaxis(
		std::sqrt(std::fabs(Sinv_Q.x_)),
		std::sqrt(std::fabs(Sinv_Q.y_)),
		std::sqrt(std::fabs(Sinv_Q.z_)));
	if ( Sinv_Q.x_<0 ) rotaxis[0] *= -1;
	if ( Sinv_Q.y_<0 ) rotaxis[1] *= -1;
	if ( Sinv_Q.z_<0 ) rotaxis[2] *= -1;

	// do we need to do fft convolution at all?
	if ( slice_x==0 && slice_y==0 && slice_z==0 ) {

		// translation only, just take dot product

		if ( debug_||debug_exact_ ) {
			TR << "no slice"  << std::endl;
			TR << "                : " << Sinv_Q.x_ << " " << Sinv_Q.y_ << " " << Sinv_Q.z_ << std::endl;
		}
		core::Real dotProd=0;
		for ( int i=0; i<(int)Npoints; ++i ) dotProd += rho[i]*Srho[i];
		conv_out.dimension( rho.u1(),rho.u2(),rho.u3() );
		conv_out = dotProd;
	} else {

		// do the fft convolution

		// Rinv maps xyz->pij (our reindexed coordinates with p is perpendicular to the symmaxis)
		// R maps pij->xyz
		numeric::xyzMatrix<Real> Rinv;
		int slice_axis = 0;
		bool R_is_identity=false;
		if ( slice_x != 0 ) {
			slice_axis = 1;
			Rinv.row_x( Vector( slice_x, slice_y, slice_z ) );
			R_is_identity = (slice_y==0 && slice_z==0);
		} else {
			Rinv.row_x( Vector( 1, 0, 0 ) );
		}

		if ( slice_x == 0 && slice_y != 0 ) {
			slice_axis = 2;
			Rinv.row_y( Vector( slice_x, slice_y, slice_z ) );
			R_is_identity = (slice_z==0);
		} else {
			Rinv.row_y( Vector( 0, 1, 0 ) );
		}
		if ( slice_x == 0 && slice_y == 0 && slice_z != 0 ) {
			slice_axis = 3;
			Rinv.row_z( Vector( slice_x, slice_y, slice_z ) );
			R_is_identity = true;
		} else {
			Rinv.row_z( Vector( 0, 0, 1 ) );
		}

		// oddball case (+ rotations of it)
		if ( S == numeric::xyzMatrix<Real>::rows(1,-1,0,   0,-1,0,   0,0,-1) ) {
			Rinv = numeric::xyzMatrix<Real>::rows( 2,-1,0,  1,0,0,  0,0,1 );
			R_is_identity = false;
			rotaxis = Vector(1,0,0);
		}
		if ( S == numeric::xyzMatrix<Real>::rows(-1,0,0, -1,1,0, 0,0,-1) ) {
			Rinv = numeric::xyzMatrix<Real>::rows( 0,1,0, -1,2,0, 0,0,1 );
			R_is_identity = false;
			rotaxis = Vector(0,1,0);
		}
		if ( S == numeric::xyzMatrix<Real>::rows(-1,0,0, 0,-1,0, 0,-1,1) ) {
			Rinv = numeric::xyzMatrix<Real>::rows(1,0,0, 0,0,1, 0,-1,2);
			R_is_identity = false;
			rotaxis = Vector(0,0,1);
		}
		if ( S == numeric::xyzMatrix<Real>::rows(1,0,-1, 0,-1,0, 0,0,-1) ) {
			Rinv = numeric::xyzMatrix<Real>::rows( 2,0,-1, 0,1,0, 1,0,0 );
			R_is_identity = false;
			rotaxis = Vector(1,0,0);
		}
		if ( S == numeric::xyzMatrix<Real>::rows(0,0,-1, -1,0,0, 0,1,-1) ) {
			Rinv = numeric::xyzMatrix<Real>::rows( 0,1,0, 1,0,0, 0,2,-1 );
			R_is_identity = false;
			rotaxis = Vector(0,0,1);
		}
		if ( S == numeric::xyzMatrix<Real>::rows(0,-1,0, -1,0,1, -1,0,0) ) {
			Rinv = numeric::xyzMatrix<Real>::rows( 0,1,0, -1,0,2, 0,0,1 );
			R_is_identity = false;
			rotaxis = Vector(0,1,0);
		}

		numeric::xyzMatrix<Real> R = numeric::inverse(Rinv);

		// Q is our symmetric transformation in our new coordinates
		numeric::xyzMatrix<Real> Q = Rinv*S*R;

		// Qstar is the FFT transformation
		//   Q-I on diagonal elements not on the slice axis
		//   partially apply the elements in the column of the slice axis (only applied to symm-sampled map)
		//        to compensate for point group sliding as we move up the symm axis in resampled groups
		numeric::xyzMatrix<Real> QstarF, QstarG;
		if ( slice_axis == 1 ) {
			QstarF = numeric::xyzMatrix<Real>::rows( 1,0,0,   Q.yx(),Q.yy()-1,Q.yz(),   Q.zx(),Q.zy(),Q.zz()-1 );
			QstarG = numeric::xyzMatrix<Real>::rows( 1,0,0,        0,Q.yy()-1,Q.yz(),        0,Q.zy(),Q.zz()-1 );
		} else if ( slice_axis == 2 ) {
			QstarF = numeric::xyzMatrix<Real>::rows( Q.xx()-1,  Q.xy(),Q.xz(), 0,1,0, Q.zx(),  Q.zy(),Q.zz()-1 );
			QstarG = numeric::xyzMatrix<Real>::rows( Q.xx()-1,       0,Q.xz(), 0,1,0, Q.zx(),       0,Q.zz()-1 );
		} else { //slice_axis == 3
			QstarF = numeric::xyzMatrix<Real>::rows( Q.xx()-1,Q.xy(),  Q.xz(), Q.yx(),Q.yy()-1,  Q.yz(), 0,0,1 );
			QstarG = numeric::xyzMatrix<Real>::rows( Q.xx()-1,Q.xy(),       0, Q.yx(),Q.yy()-1,       0, 0,0,1 );
		}

		if ( debug_||debug_exact_ ) {
			TR << "  slice on axis :    " << slice_axis << std::endl;
			TR << "slice along axis: " << slice_x << " " << slice_y << " " << slice_z << std::endl;
			TR << "                : " << Sinv_Q.x_ << " " << Sinv_Q.y_ << " " << Sinv_Q.z_ << std::endl;
			TR << "                : " << rotaxis[0] << " " << rotaxis[1] << " " << rotaxis[2] << std::endl;
			TR << "         S = [" << S.xx() << "," <<  S.xy() << "," <<  S.xz() << ";  "
				<< S.yx() << "," <<  S.yy() << "," <<  S.yz() << ";  "
				<< S.zx() << "," <<  S.zy() << "," <<  S.zz() << "]" << std::endl;
			if ( !R_is_identity ) {
				TR << "      Rinv = [" << Rinv.xx() << "," << Rinv.xy() << "," << Rinv.xz() << ";  "
					<< Rinv.yx() << "," << Rinv.yy() << "," << Rinv.yz() << ";  "
					<< Rinv.zx() << "," << Rinv.zy() << "," << Rinv.zz() << "]" << std::endl;
				TR << "         R = [" << R.xx() << "," <<  R.xy() << "," <<  R.xz() << ";  "
					<< R.yx() << "," <<  R.yy() << "," <<  R.yz() << ";  "
					<< R.zx() << "," <<  R.zy() << "," <<  R.zz() << "]" << std::endl;
			}
			TR << "         Q = [" << Q.xx() << "," << Q.xy() << "," << Q.xz() << ";  "
				<< Q.yx() << "," << Q.yy() << "," << Q.yz() << ";  "
				<< Q.zx() << "," << Q.zy() << "," << Q.zz() << "]" << std::endl;
			TR << "        Qf = [" << QstarF.xx() << "," <<  QstarF.xy() << "," <<  QstarF.xz() << ";  "
				<< QstarF.yx() << "," <<  QstarF.yy() << "," <<  QstarF.yz() << ";  "
				<< QstarF.zx() << "," <<  QstarF.zy() << "," <<  QstarF.zz() << "]" << std::endl;
			TR << "        Qg = [" << QstarG.xx() << "," <<  QstarG.xy() << "," <<  QstarG.xz() << ";  "
				<< QstarG.yx() << "," <<  QstarG.yy() << "," <<  QstarG.yz() << ";  "
				<< QstarG.zx() << "," <<  QstarG.zy() << "," <<  QstarG.zz() << "]" << std::endl;
		}

		// transform
		transform_map_offset0( rho, R, working_trans);
		transform_map_offset0( Srho, R, working_strans);

		// transform in each plane
		// fft in each plane
		FArray3D<Real> rworking_trans, rworking_strans;
		transform_map_offset0( working_trans, QstarF, rworking_trans);
		transform_map_offset0( working_strans, QstarG, rworking_strans);

		fft2dslice( rworking_trans, Fworking_trans, slice_axis );
		fft2dslice( rworking_strans, Fworking_strans, slice_axis );
		for ( int i=0; i<(int)Npoints; ++i ) Fworking_trans[i] *= std::conj(Fworking_strans[i]);


		if ( R_is_identity ) {
			ifft2dslice(Fworking_trans, conv_out, slice_axis);
		} else {
			ifft2dslice(Fworking_trans, working_trans, slice_axis);
			transform_map_offset0( working_trans, Rinv, conv_out);
		}
		project_along_axis( conv_out, rotaxis );
	}
}

// get the score per interface
void
CrystFFTDock::get_interfaces_allatom(
	Pose pose, // make a copy
	utility::vector1<core::kinematics::RT> const &/*rts*/,
	numeric::xyzMatrix<Real> R,
	numeric::xyzVector<Real> xyz_grid,
	utility::vector1<SingleInterface> &allInterfaces )
{
	pose.apply_transform_Rx_plus_v( R, xyz_grid );

	//////////////////////
	////
	//// recenter
	////
	//////////////////////
	Size nres = pose.size();

	Vector com(0,0,0);
	for ( Size i=1; i<= nres; ++i ) com += pose.residue(i).xyz(2);
	com /= nres;

	Real mindis2(1e6);
	for ( Size i=1; i<= nres; ++i ) {
		Real const dis2( com.distance_squared(  pose.residue(i).xyz(2) ) );
		if ( dis2 < mindis2 ) {
			mindis2 = dis2;
		}
	}
	Size nsymm = sg_.nsymmops();
	Size bestxform=0;
	Vector bestoffset(0,0,0);
	mindis2=1e6;
	com = sg_.c2f()*com;
	for ( Size i=1; i<=nsymm; ++i ) {
		Vector foffset = sg_.symmop(i).get_rotation()*com + sg_.symmop(i).get_translation(), rfoffset;
		rfoffset[0] = min_mod( foffset[0], 1.0 );
		rfoffset[1] = min_mod( foffset[1], 1.0 );
		rfoffset[2] = min_mod( foffset[2], 1.0 );
		Real dist = (sg_.f2c()*rfoffset).length_squared();
		if ( dist<mindis2 ) {
			mindis2=dist;
			bestxform=i;
			bestoffset = foffset - rfoffset;
		}
	}
	numeric::xyzMatrix<Real> Rx = sg_.f2c()*sg_.symmop(bestxform).get_rotation()*sg_.c2f();
	numeric::xyzVector<Real> Tx = sg_.f2c()*(sg_.symmop(bestxform).get_translation() - bestoffset);
	pose.apply_transform_Rx_plus_v( Rx,Tx );

	//////////////////////
	////
	//// find contacts
	////
	//////////////////////
	Real contact_dist=interfacedist_;
	Real radius = 0;
	utility::vector1<Vector> monomer_cbs;
	for ( Size i=1; i<= pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		Size atm=(pose.residue(i).aa() == core::chemical::aa_gly)?2:5;
		Vector cb_i = pose.residue(i).xyz(atm);
		monomer_cbs.push_back(cb_i);
		radius = std::max( (cb_i).length_squared() , radius );
	}
	Size nres_monomer = monomer_cbs.size();
	radius = sqrt(radius);

	for ( int s=1; s<=(int)sg_.nsymmops(); ++s ) {
		numeric::xyzMatrix<Real> R_i = sg_.symmop(s).get_rotation();

		for ( int i=-1; i<=1; ++i ) {
			for ( int j=-1; j<=1; ++j ) {
				for ( int k=-1; k<=1; ++k ) {
					if ( s==1 && i==0 && j==0 && k==0 ) continue;

					numeric::xyzVector<Real> T_i = sg_.symmop(s).get_translation() + numeric::xyzVector<Real>(i,j,k);

					// pass 1 check vrt-vrt dist to throw out very distant things
					Real disVRT = T_i.length();
					if ( disVRT>contact_dist+2*radius ) continue;

					// pass 2 check cb-cb dists
					numeric::xyzMatrix<core::Real> R_i_realspace =  sg_.f2c()*R_i*sg_.c2f();
					Size contact=0;
					for ( Size jj=1; jj<= nres_monomer; ++jj ) {
						Vector Xi = R_i_realspace*monomer_cbs[jj] + sg_.f2c()*T_i;
						for ( Size kk=1; kk<= nres_monomer; ++kk ) {
							if ( (Xi-monomer_cbs[kk]).length_squared() < contact_dist*contact_dist ) contact++;
						}
					}

					if ( contact>mininterface_ ) {
						// make minipose & score
						core::pose::Pose poseCopy = pose;
						Real motif_score = (Real)contact;
						allInterfaces.push_back( SingleInterface( R_i, T_i, motif_score ) );
					}
				}
			}
		}
	}
}

core::Real
CrystFFTDock::get_clash_score_exact(
	numeric::xyzVector<int> xyz_grid,
	numeric::xyzMatrix<Real> R,
	numeric::xyzVector<Real> T,
	FArray3D<Real> const &r_rho_ca
) {
	FArray3D<Real> shift_rho_ca, s_shift_rho_ca;
	shift_rho_ca.dimension( grid_[0] , grid_[1] , grid_[2] );
	Size Npoints = grid_[0]*grid_[1]*grid_[2];

	for ( int z=1; z<=grid_[2]; ++z ) {
		for ( int y=1; y<=grid_[1]; ++y ) {
			for ( int x=1; x<=grid_[0]; ++x ) {
				int cx = pos_mod( x-xyz_grid[0]-1 , grid_[0] ) + 1;
				int cy = pos_mod( y-xyz_grid[1]-1 , grid_[1] ) + 1;
				int cz = pos_mod( z-xyz_grid[2]-1 , grid_[2] ) + 1;
				shift_rho_ca(x,y,z) = r_rho_ca(cx,cy,cz);
			}
		}
	}
	transform_map( shift_rho_ca, R,T, s_shift_rho_ca);

	Real retval = 0;
	for ( int i=0; i<(int)Npoints; ++i ) {
		retval += shift_rho_ca[i] * s_shift_rho_ca[i];
	}
	//  Use it to dump maps

	//  static int dump_count=0;
	//  std::ostringstream oss; oss << "symmmap" << dump_count << ".mrc";
	//  writeMRC( s_shift_rho_ca, oss.str() );
	//
	//  std::ostringstream oss2; oss2 << "map" << dump_count << ".mrc";
	//  writeMRC( shift_rho_ca, oss2.str() );
	//  dump_count++;
	return retval;
}




void
CrystFFTDock::apply( Pose & pose) {
	// set up crystal info
	sg_.set_spacegroup( option[ crystdock::spacegroup ] );
	sg_.set_parameters(
		option[ crystdock::A ],option[ crystdock::B ],option[ crystdock::C ],
		option[ crystdock::alpha ],option[ crystdock::beta ],option[ crystdock::gamma ] );

	Real rotstep = rot_step_;
	Real maxclash = maxclash_;

	// center pose at origin
	numeric::xyzVector<Real> native_shift = center_pose_at_origin( pose );

	// get SS
	core::scoring::dssp::Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	// lookup symmops
	utility::vector1<core::kinematics::RT> const &rts = sg_.symmops();
	protocols::cryst::CheshireCell cc = sg_.cheshire_cell();

	if ( debug_||debug_exact_ ) {
		TR << "search [ " << cc.low[0] <<","<< cc.low[1] <<","<< cc.low[2] << " : "
			<< cc.high[0] <<","<< cc.high[1] <<","<< cc.high[2] << "]" << std::endl;
	}

	// compute rho_ca, rho_cb
	FArray3D<Real> rho_ca, rho_cb, r_rho_ca, r_rho_cb;

	setup_maps( pose, rho_ca, rho_cb, trans_step_);

	Size Npoints = grid_[0]*grid_[1]*grid_[2];
	Real voxel_volume = sg_.volume() / Npoints;

	numeric::xyzVector<Size> ccIndexLow(
		(Size)std::floor(cc.low[0]*grid_[0]+1.5),
		(Size)std::floor(cc.low[1]*grid_[1]+1.5),
		(Size)std::floor(cc.low[2]*grid_[2]+1.5));
	numeric::xyzVector<Size> ccIndexHigh(
		(Size)std::floor(cc.high[0]*grid_[0]+0.5),
		(Size)std::floor(cc.high[1]*grid_[1]+0.5),
		(Size)std::floor(cc.high[2]*grid_[2]+0.5));

	for ( int i=0; i<3; i++ ) {
		if ( ccIndexHigh[i]<ccIndexLow[i] ) {
			ccIndexHigh[i]=ccIndexLow[i];
		}
	}

	if ( debug_||debug_exact_ ) {
		writeMRC( rho_ca, "ca_mask.mrc", true, true );
		writeMRC( rho_cb, "cb_mask.mrc", true, true );
	}

	// space for intermediate results
	FArray3D<Real> working_s, conv_out;

	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
	numeric::xyzMatrix<Real> r_local;

	// the collection of hits
	InterfaceHitDatabase IDB(2*nmodels_);  // oversample 2x to start
	core::Size nnonclashing = 0, nconnected = 0;

	// foreach rotation
	numeric::UniformRotationSampler urs( rotstep );
	for ( core::Size i=1; i<=rts.size(); ++i ) {
		numeric::xyzMatrix<Real> S = rts[i].get_rotation();
		if ( S.xy()==0 && S.xz()==0 && S.yz()==0 && S.yx()==0 && S.zx()==0 && S.zy()==0 && S.xx()==1 && S.yy()==1 && S.zz()==1 ) {
			continue; // identity
		}
		numeric::xyzMatrix<Real> Scart = i2c_*S*c2i_;
		urs.remove_redundant( Scart );
	}

	Size rot_lb = 1, rot_ub = urs.nrots();
	if ( option[ crystdock::rotnum ].user() ) {
		rot_lb = rot_ub = option[ crystdock::rotnum ]();
	}
	Size sym_lb = 2, sym_ub = rts.size();
	if ( option[ crystdock::symnum ].user() ) {
		sym_lb = sym_ub = option[ crystdock::symnum ]();
	}

	// hook into other code for evaluating native interfaces
	if ( eval_native_ ) {
		mininterface_ = 0.0001;  // we want to evaluate all interfaces

		utility::vector1<Size> symmoplist;

		numeric::xyzVector<Real> offset_grid;
		numeric::xyzVector<int> offset_grid_pt;
		utility::vector1<SingleInterface> iinfo;

		offset_grid = c2i_*native_shift;

		// round to nearest grid
		offset_grid_pt[0] = (int) std::floor( offset_grid[0] + 0.5 );
		offset_grid_pt[1] = (int) std::floor( offset_grid[1] + 0.5 );
		offset_grid_pt[2] = (int) std::floor( offset_grid[2] + 0.5 );
		offset_grid = numeric::xyzVector<Real>(offset_grid_pt[0], offset_grid_pt[1], offset_grid_pt[2]);
		native_shift = i2c_*offset_grid;

		// exact ca overlap volume
		core::Real ca_score = 0;
		ca_score += resample_maps_and_get_self( rho_ca, rho_cb, identity, sg_, r_rho_ca, r_rho_cb, iinfo );
		for ( int s=2; s<=(int)rts.size(); ++s ) {
			symmoplist.push_back(s);
			ca_score += get_clash_score_exact( offset_grid_pt, rts[s].get_rotation(), rts[s].get_translation(), r_rho_ca );
		}

		// exact cb contact count
		iinfo.clear();
		get_interfaces_allatom( pose, rts, identity, native_shift, iinfo );
		core::Real cb_score = get_interface_score( pose, iinfo, rts );

		TR << "OVERLAP score = " << ca_score << std::endl;
		TR << "INTERACTION score (cb count)= " << cb_score << std::endl;

		// add this to database with a super-good score
		IDB.add_interface( InterfaceHit( 1e9, native_shift[0],native_shift[1],native_shift[2], 0, iinfo ) );
	}

	for ( int ctr=(int)rot_lb; ctr<=(int)rot_ub ; ++ctr ) {
		TR << "Rotation " << ctr << " of " << urs.nrots() << std::endl;
		urs.get(ctr, r_local);

		// interface_map:           stores the exact transformation and area of all interfaces > minintarea
		// ambiguous_interface_map: stores the index (in 'rts') of all interfaces with a sum > minintarea
		// p1_interface_map:        self interactions, independent of translation
		FArray3D< Real > sum_interface_area;
		FArray3D< utility::vector1<Size> > ambiguous_interface_map;
		utility::vector1<SingleInterface> p1_interface_map;

		ambiguous_interface_map.dimension( grid_[0] , grid_[1] , grid_[2] );
		sum_interface_area.dimension( grid_[0] , grid_[1] , grid_[2] ); sum_interface_area=0;
		Real self_ca = resample_maps_and_get_self( rho_ca, rho_cb, r_local, sg_, r_rho_ca, r_rho_cb, p1_interface_map );

		if ( debug_||debug_exact_ ) {
			std::ostringstream oss1; oss1 << "rot"<<ctr<<".mrc";
			writeMRC( r_rho_ca, oss1.str(), false, true );
			std::ostringstream oss2; oss2 << "rot"<<ctr<<".pdb";
			dump_transformed_pdb( pose, InterfaceHit( 0.,0.,0.,0., ctr, utility::vector1<SingleInterface>() ), urs, oss2.str() );
		}

		if ( self_ca >= maxclash ) {
			TR << "   self clashing!" << std::endl;
			if ( debug_ ) writeMRC( r_rho_ca, "clash.mrc" );
			continue; // next rotation
		}

		// P1 interactions are OK, now compute the rest
		FArray3D<Real> collision_map, ex_collision_map;
		collision_map.dimension( grid_[0] , grid_[1] , grid_[2] );
		collision_map=self_ca;

		Real sum_interface_p1 = 0;
		for ( int i=1; i<=(int)p1_interface_map.size(); ++i ) {
			sum_interface_p1 += p1_interface_map[i].cb_overlap_;
		}
		sum_interface_area = sum_interface_p1;

		if ( debug_exact_ ) {
			ex_collision_map.dimension( grid_[0] , grid_[1] , grid_[2] );
			ex_collision_map=self_ca;
		}

		// do the convolution with each symmop to find configurations that:
		//   (a) are not clashing
		//   (b) _might_ be fully connected in the lattice
		for ( int s=sym_lb; s<=(int)sym_ub; ++s ) {   // s==1 is always identity
			numeric::xyzMatrix<Real> s_i = rts[s].get_rotation();
			numeric::xyzMatrix<Real> s_inv = numeric::inverse(s_i);
			numeric::xyzVector<Real> t_inv = s_inv*(-rts[s].get_translation());

			transform_map( r_rho_ca, s_inv, t_inv, working_s);

			// all the magic is in here
			do_convolution( r_rho_ca, working_s, s_i, conv_out);

			for ( int z=(int)ccIndexLow[2]; z<=(int)ccIndexHigh[2]; ++z ) {
				for ( int y=(int)ccIndexLow[1]; y<=(int)ccIndexHigh[1]; ++y ) {
					for ( int x=(int)ccIndexLow[0]; x<=(int)ccIndexHigh[0]; ++x ) {
						collision_map(x,y,z) += conv_out(x,y,z);
					}
				}
			}

			// debug: exact CA fft
			if (  debug_exact_ ) {
				for ( int sz=(int)ccIndexLow[2]-1; sz<=(int)ccIndexHigh[2]-1; ++sz ) {
					for ( int sy=(int)ccIndexLow[1]-1; sy<=(int)ccIndexHigh[1]-1; ++sy ) {
						for ( int sx=(int)ccIndexLow[0]-1; sx<=(int)ccIndexHigh[0]-1; ++sx ) {
							ex_collision_map(sx+1,sy+1,sz+1) +=
								get_clash_score_exact( numeric::xyzVector<int>(sx,sy,sz), rts[s].get_rotation(), rts[s].get_translation(), r_rho_ca);
						}
					}
				}
			}


			// do CB fft
			transform_map( r_rho_cb, s_inv, t_inv, working_s);
			do_convolution( r_rho_cb, working_s, s_i, conv_out);
			for ( int z=(int)ccIndexLow[2]; z<=(int)ccIndexHigh[2]; ++z ) {
				for ( int y=(int)ccIndexLow[1]; y<=(int)ccIndexHigh[1]; ++y ) {
					for ( int x=(int)ccIndexLow[0]; x<=(int)ccIndexHigh[0]; ++x ) {
						if ( conv_out(x,y,z)*voxel_volume > mininterface_ ) {
							sum_interface_area(x,y,z) += conv_out(x,y,z)*voxel_volume;
							ambiguous_interface_map(x,y,z).push_back( s );
						}
					}
				}
			}
		}

		if ( debug_ || debug_exact_ ) {
			std::ostringstream oss; oss << "collisionmap_"<<ctr<<".mrc";
			FArray3D<Real> collision_map_dump = collision_map;
			for ( int i=0; i<(int)Npoints; ++i ) collision_map_dump[i] = /*maxclash-*/collision_map[i];
			writeMRC( collision_map_dump, oss.str() );

			if ( debug_exact_ ) {
				std::ostringstream oss2; oss2 << "ex_collisionmap_"<<ctr<<".mrc";
				for ( int i=0; i<(int)Npoints; ++i ) collision_map_dump[i] = /*maxclash-*/ex_collision_map[i];
				writeMRC( collision_map_dump, oss2.str() );

				// get correl
				Real x2=0,y2=0,xy=0,x=0,y=0, Ncc=0;
				for ( int sz=(int)ccIndexLow[2]-1; sz<=(int)ccIndexHigh[2]-1; ++sz ) {
					for ( int sy=(int)ccIndexLow[1]-1; sy<=(int)ccIndexHigh[1]-1; ++sy ) {
						for ( int sx=(int)ccIndexLow[0]-1; sx<=(int)ccIndexHigh[0]-1; ++sx ) {
							x2 += collision_map(sx+1,sy+1,sz+1) * collision_map(sx+1,sy+1,sz+1);
							y2 += ex_collision_map(sx+1,sy+1,sz+1) * ex_collision_map(sx+1,sy+1,sz+1);
							xy += collision_map(sx+1,sy+1,sz+1) * ex_collision_map(sx+1,sy+1,sz+1);
							x  += collision_map(sx+1,sy+1,sz+1);
							y  += ex_collision_map(sx+1,sy+1,sz+1);
							Ncc++;
						}
					}
				}
				Real correl = 1, slope = 0;
				Real sx = (Ncc*x2-x*x);
				Real sy = (Ncc*y2-y*y);
				if ( sx*sy > 0 ) {
					correl = (Ncc*xy - x*y) / std::sqrt( (sx) * (sy) );
					slope = correl * sx/sy;
				}
				TR << "correl = " << correl << "   scale = " << slope << std::endl;

				collision_map = ex_collision_map;
			}
		}

		// finally add nonclashing interfaces to the DB
		Real mininterfacesum_filter = std::max( 3*mininterface_, mininterface_sum_);
		for ( int z=(int)ccIndexLow[2]; z<=(int)ccIndexHigh[2]; ++z ) {
			for ( int y=(int)ccIndexLow[1]; y<=(int)ccIndexHigh[1]; ++y ) {
				for ( int x=(int)ccIndexLow[0]; x<=(int)ccIndexHigh[0]; ++x ) {
					if ( collision_map(x,y,z) < maxclash ) {
						nnonclashing++;
					}

					if ( collision_map(x,y,z) < maxclash && sum_interface_area(x,y,z) > mininterfacesum_filter ) {
						// get_interface_score populates iinfo
						//    then computes the weakest connection necessary to construct the lattice
						utility::vector1<SingleInterface> iinfo; // = p1_interface_map;
						numeric::xyzVector<Real> xyz((Real)x-1,(Real)y-1,(Real)z-1);
						xyz = i2c_*xyz;

						get_interfaces_allatom( pose, rts, r_local, xyz, iinfo);

						Real score_xyz_max=0;
						for ( int inum=1; inum<=(int)iinfo.size(); ++inum ) {
							score_xyz_max = std::max(score_xyz_max, iinfo[inum].cb_overlap_);
						}

						Real score_xyz = get_interface_score (pose, iinfo, rts );

						if ( score_xyz > mininterface_ ) {
							nconnected++;
							IDB.add_interface( InterfaceHit( score_xyz, score_xyz_max, xyz[0],xyz[1],xyz[2], ctr, iinfo ) );
						}
					}
				}
			}
		}
		if ( IDB.size()>0 ) {
			TR << IDB.size() << " of " << nnonclashing << " nonclashing and " << nconnected << " connected "
				<< " configurations; min_score = " << IDB.top().score << std::endl;
		} else {
			TR << IDB.size() << " of " << nnonclashing << " nonclashing configurations" << std::endl;
		}

	}


	utility::vector1< InterfaceHit > ih_vec;
	int nhits = IDB.size();
	for ( int i=1; i<=nhits; ++i ) {
		InterfaceHit ih = IDB.pop();
		ih_vec.push_back( ih );
	}
	std::reverse(ih_vec.begin(),ih_vec.end());


	Real pose_radius = get_radius_of_pose( pose );
	Real cluster_cutoff = cluster_cutoff_;
	utility::vector1< InterfaceHit > ih_vec_clustered;

	// don't bother searching if we're not clustering
	if ( cluster_cutoff==0 ) {
		ih_vec_clustered = ih_vec;
	} else {
		for ( int i=1; i<=nhits; ++i ) {
			int nclust = ih_vec_clustered.size();
			bool found_match=false;
			for ( int j=1; j<=nclust && !found_match; ++j ) {
				Real dist_ij = get_transform_distance( ih_vec[i], ih_vec_clustered[j] , urs, pose_radius);
				found_match = (dist_ij <= cluster_cutoff);
			}
			if ( !found_match ) {
				ih_vec_clustered.push_back( ih_vec[i] );
			}
		}
	}

	int nhits_after_cluster = std::min( ih_vec_clustered.size() , nmodels_ );
	TR << "Have " << nhits_after_cluster << " after clustering." << std::endl;

	//////////////////////////////////////
	std::string base_name = protocols::jd2::current_input_tag();
	utility::vector1< std::string > temp_out_names= utility::split( base_name );
	utility::file::FileName out_name = utility::file::combine_names( temp_out_names );
	base_name = out_name.base();
	std::string outname = option[ out::prefix ]+base_name+option[ out::suffix ];

	if ( nhits_after_cluster>0 ) {
		std::ofstream out(outname.c_str());
		for ( int i=1; i<=nhits_after_cluster; ++i ) {
			InterfaceHit ih = ih_vec_clustered[i];
			out << outname << " " << sg_.name() << " "
				<< sg_.A() << " " << sg_.B() << " " << sg_.C() << " " << sg_.alpha() << " " << sg_.beta() << " " << sg_.gamma()
				<< " " << ih.to_string( urs ) << std::endl;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols;
	using namespace protocols::jd2;

	SequenceMoverOP seq( new SequenceMover() );

	std::string strategy = option[ crystdock::mode ]();

	// force some options
	option[ basic::options::OptionKeys::in::preserve_crystinfo ].value(true);
	option[ basic::options::OptionKeys::symmetry::symmetry_definition ].value("dummy");  // get right rot type set for relax

	if ( strategy == "fftdock" ) {
		seq->add_mover( MoverOP(new CrystFFTDock()) );
		option[ out::nooutput ].value(true);
	} else if ( strategy == "cluster" ) {
		seq->add_mover(  MoverOP(new CrystCluster()) );
		option[ out::nooutput ].value(true);
	} else if ( strategy == "slidedock" ) {
		seq->add_mover(  MoverOP(new CrystSlideDock()) );
		option[ out::nooutput ].value(true);
		option[ basic::options::OptionKeys::symmetry::detect_bonds ].value(false); // avoid forming disulfides across the crystal interface
	}  else if ( strategy == "design" ) {
		seq->add_mover(  MoverOP(new CrystDesign(false)) );
	} else if ( strategy == "revert" ) {
		seq->add_mover(  MoverOP(new CrystDesign(true)) );
	} else if ( strategy == "relax" ) {
		seq->add_mover(  MoverOP(new CrystRelax()) );
	} else if ( strategy == "score" ) {
		seq->add_mover(  MoverOP(new CrystRelax(false)) );
	} else if ( strategy == "rms" ) {
		SilentFileJobOutputterOP jobout (new SilentFileJobOutputter);
		jobout->set_write_no_structures();
		jobout->set_write_separate_scorefile(true);
		protocols::jd2::JobDistributor::get_instance()->set_job_outputter( JobDistributorFactory::create_job_outputter( jobout ));
		seq->add_mover(  MoverOP(new CrystRMS()) );
	} else {
		utility_exit_with_message( "Unknown mode:"+ option[ crystdock::mode ]() );
	}

	// main loop
	protocols::jd2::JobDistributor::get_instance()->go( seq );

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
	try {
		NEW_OPT(crystdock::mode, "mode", "dock");
		NEW_OPT(crystdock::spacegroup, "spacegroup (no spaces)", "P23");
		NEW_OPT(crystdock::A, "unit cell A", 50);
		NEW_OPT(crystdock::B, "unit cell B", 50);
		NEW_OPT(crystdock::C, "unit cell C", 50);
		NEW_OPT(crystdock::alpha, "unit cell alpha", 90);
		NEW_OPT(crystdock::beta, "unit cell beta", 90);
		NEW_OPT(crystdock::gamma, "unit cell gamma", 90);

		// FFT MEASURES
		NEW_OPT(crystdock::maxclash, "max allowed clashscore", 20);
		NEW_OPT(crystdock::mininterface, "min allowed interface (#cbs)", 50);
		NEW_OPT(crystdock::mininterfacesum, "min interface volume", 4000);
		NEW_OPT(crystdock::trans_step, "translational stepsize (A)", 1);
		NEW_OPT(crystdock::rot_step, "rotational stepsize (degrees) ([debug] 0 searches input rotation only)", 10);

		NEW_OPT(crystdock::nmodels, "number of models to output (docking)", 1000);
		NEW_OPT(crystdock::motif_initial_cut, "initial cut on motif energy", -20);
		NEW_OPT(crystdock::motif_final_cut, "final cut on motif energy", -100);

		NEW_OPT(crystdock::energy_cut, "cut on weakest interface energy", 1e6);
		NEW_OPT(crystdock::favor_native_bonus, "favor native bonus (design only)", 0.25);
		NEW_OPT(crystdock::ssonly, "limit interface calcs to secstruct only", false);
		NEW_OPT(crystdock::rotnum, "[debug] only run a single rotation", 0);
		NEW_OPT(crystdock::symnum, "[debug] only run a single symmop", 0);
		NEW_OPT(crystdock::debug, "[debug] dump intermediate info", false);
		NEW_OPT(crystdock::debug_exact, "[debug] debug mode with exact (non-FFT) calculations (slow!)", false);
		NEW_OPT(crystdock::eval_native, "[debug] evaluate input structure without docking", false);
		NEW_OPT(crystdock::sigwidth, "sigwidth", 18.00);
		NEW_OPT(crystdock::interfacedist, "interface_dist", 6.0);
		NEW_OPT(crystdock::interface_sigwidth, "interface_sigwidth", 1.00);
		NEW_OPT(crystdock::cluster_cutoff, "cluster_cutoff", 0.0);
		NEW_OPT(crystdock::dock_ncycle, "dock_ncycle", 200);
		NEW_OPT(crystdock::dock_vdw_cut, "dock_vdw", 2.0);
		NEW_OPT(crystdock::ncyc, "ncyc", 1);
		NEW_OPT(crystdock::nmut, "nmut", 12);
		NEW_OPT(crystdock::K, "K", 14.0);
		NEW_OPT(crystdock::fixedres, "fixedpos", utility::vector1<Size>(0));
		NEW_OPT(crystdock::nofilt, "nofilt", false);
		NEW_OPT(crystdock::adjust_ref_weights, "adjust_ref_weights", utility::vector1<std::string>(0));

		NEW_OPT(crystdock::hits_in, "hits_in", "");
		NEW_OPT(crystdock::hits_to_dock, "hits_to_dock", utility::vector1<Size>(0));
		NEW_OPT(crystdock::hits_out, "hits_out", "");

		devel::init( argc, argv );
		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

