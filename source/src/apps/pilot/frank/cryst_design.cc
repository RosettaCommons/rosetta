/// @file
/// @brief


#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/types.hh>
#include <core/pack/task/ResfileReader.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/docking/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/ddG.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>


#include <fstream>
#include <iostream>
#include <math.h>

#include <sstream>
#include <string>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

using std::string;
using ObjexxFCL::string_of;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using namespace core::scoring::packing;
using namespace ObjexxFCL::fmt;
using core::scoring::ScoreFunctionOP;
using core::pose::Pose;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef vector1<Size> Sizes;

static basic::Tracer TR("cryst.design");

OPT_1GRP_KEY(String, crystdes, prefix)
OPT_1GRP_KEY(Real, crystdes, alpha)
OPT_1GRP_KEY(Boolean, crystdes, cen_min)
OPT_1GRP_KEY(Boolean, crystdes, smallres_only)
OPT_1GRP_KEY(Real, crystdes, rot_step)
OPT_1GRP_KEY(Real, crystdes, trans_step)
OPT_1GRP_KEY(Real, crystdes, rot_max)
OPT_1GRP_KEY(Real, crystdes, trans_max)
OPT_1GRP_KEY(Real, crystdes, vdw)
OPT_1GRP_KEY(Real, crystdes, sa1)
OPT_1GRP_KEY(Real, crystdes, sa2)
OPT_1GRP_KEY(Real, crystdes, sa3)
OPT_1GRP_KEY(Real, crystdes, sc1)
OPT_1GRP_KEY(Real, crystdes, sc2)
OPT_1GRP_KEY(Real, crystdes, sc3)
OPT_1GRP_KEY(Real, crystdes, ddg)
OPT_1GRP_KEY(Real, crystdes, air)
OPT_1GRP_KEY(Boolean, crystdes, print_des)

///////////////////////////////////////////////////////////////////////////////

utility::vector1<Real>
sidechain_sasa(Pose const & pose, Real probe_radius) {
	using core::id::AtomID;
	utility::vector1<Real> rsd_sasa(pose.n_residue(),0.0);
	core::id::AtomID_Map<Real> atom_sasa;
	core::id::AtomID_Map<bool> atom_mask;
	core::pose::initialize_atomid_map(atom_sasa,pose,0.0);
	core::pose::initialize_atomid_map(atom_mask,pose,false);
	for (Size i = 1; i <= pose.n_residue(); i++) {
		for (Size j = 1; j <= pose.residue(i).nheavyatoms(); j++) {
			atom_mask[AtomID(j,i)] = true;
		}
	}
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false, atom_mask );
	utility::vector1<Real> sc_sasa(pose.n_residue(),0.0);
	for(Size i = 1; i <= pose.n_residue(); i++) {
		// Use CA as the side chain for Glys
		if(pose.residue(i).name3()=="GLY") sc_sasa[i] += atom_sasa[AtomID(2,i)];		
		for(Size j = 5; j <= pose.residue(i).nheavyatoms(); j++) {
			sc_sasa[i] += atom_sasa[AtomID(j,i)];
		}
	}
	return sc_sasa;
}

void
new_sc(Pose &pose, Real& sa1, Real& sc1, Real& sa2, Real& sc2, Real& sa3, Real& sc3 ) {
	using namespace core;

	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);

	utility::vector1<Real> sas;
	utility::vector1<Real> scs;

	for (Size i=2; i<=symm_info->subunits(); ++i) {
		core::scoring::sc::ShapeComplementarityCalculator scc;
		scc.Init();

		// Figure out which chains touch chain A, and add the residues from those chains
		// into the sc surface objects	
		Size nres_monomer = symm_info->num_independent_residues();
		for (Size ir=1; ir<=nres_monomer; ++ir) {
			scc.AddResidue(0, pose.residue(ir));
		}
	
		bool contact = false;
		Size start = (i-1)*nres_monomer;
		for (Size ir=1; ir<=nres_monomer; ir++) {
			if (pose.energies().residue_total_energies(ir+start)[core::scoring::fa_atr] < 0) {
				contact = true;
				break;
			}
		}
		if (contact) {
			for (Size ir=1; ir<=nres_monomer; ir++) {
				scc.AddResidue(1, pose.residue(ir+start));
			}
			if (scc.Calc()) {
				scs.push_back( scc.GetResults().sc );
				sas.push_back( scc.GetResults().surface[2].trimmedArea );
			}
		}
	}

	// find max uniques
	sa1=sa2=sa3=0;
	sc1=sc2=sc3=0;
	for (Size i=1; i<=sas.size(); ++i) {
		// check if duplicate
		if ( fabs(sa1-sas[i]) < 1e-1 || fabs(sa2-sas[i]) < 1e-1 || fabs(sa3-sas[i]) < 1e-1 ) continue;
		if (sas[i] > sa1) {
			sa3 = sa2; sc3 = sc2;
			sa2 = sa1; sc2 = sc1;
			sa1 = sas[i]; sc1 = scs[i];
		} else if (sas[i] > sa2) {
			sa3 = sa2; sc3 = sc2;
			sa2 = sas[i]; sc2 = scs[i];
		} else if (sas[i] > sa3) {
			sa3 = sas[i]; sc3 = scs[i];
		}
	}
}


void
design(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos, bool small_only) {
	using namespace core;
	using namespace pack;
	using namespace task;
	using namespace conformation;
	using namespace conformation::symmetry;
	using namespace scoring;
	using namespace chemical;

	// Set allowed AAs.
	vector1<bool> allowed_aas(20, true);
	allowed_aas[aa_cys] = false;
	allowed_aas[aa_gly] = false;
	allowed_aas[aa_pro] = false;

	if(small_only == true) {
		allowed_aas[aa_ala] = true;
		allowed_aas[aa_asp] = true;
		allowed_aas[aa_glu] = false;
		allowed_aas[aa_phe] = false;
		allowed_aas[aa_his] = false;
		allowed_aas[aa_ile] = true;
		allowed_aas[aa_lys] = false;
		allowed_aas[aa_leu] = true;
		allowed_aas[aa_met] = false;
		allowed_aas[aa_asn] = true;
		allowed_aas[aa_pro] = false;
		allowed_aas[aa_gln] = false;
		allowed_aas[aa_arg] = false;
		allowed_aas[aa_ser] = true;
		allowed_aas[aa_thr] = true;
		allowed_aas[aa_val] = true;
		allowed_aas[aa_trp] = false;
		allowed_aas[aa_tyr] = false;
	}

	// Get the symmetry info and make the packer task
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	PackerTaskOP task( TaskFactory::create_packer_task( pose ));

	if (basic::options::option[basic::options::OptionKeys::packing::resfile].user()) {
		core::pack::task::parse_resfile( pose, *task);
	}

	// Set which residues can be designed
	for (Size i=1; i<=pose.n_residue(); i++) {
		if (!sym_info->bb_is_independent(i)) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else if (pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY") {
			// Don't mess with Pros or Glys at the interfaces
			task->nonconst_residue_task(i).prevent_repacking();
		} else if (find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else {
			bool temp = allowed_aas[pose.residue(i).aa()];
			allowed_aas[pose.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
			task->nonconst_residue_task(i).or_include_current(true);
			task->nonconst_residue_task(i).initialize_from_command_line();
			allowed_aas[pose.residue(i).aa()] = temp;
		}
	}

	// Actually perform design.
	make_symmetric_PackerTask_by_truncation(pose, task);
	protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
	packer->apply(pose);
}


void
repack(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos) {
  using namespace core;
  using namespace pack;
  using namespace task;
  using namespace conformation;
  using namespace conformation::symmetry;
  using namespace scoring;
  using namespace chemical;

  // Get the symmetry info and make the packer task
  SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
  PackerTaskOP task( TaskFactory::create_packer_task( pose ));

  // Set which residues can be repacked
  for (Size i=1; i<=pose.n_residue(); i++) {
    if (!sym_info->bb_is_independent(i)) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else if (find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else {
			vector1<bool> allowed_aas(20, false);
      allowed_aas[pose.residue(i).aa()] = true;
      task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
      task->nonconst_residue_task(i).initialize_from_command_line();
    }
  }

  // Actually repack.
  make_symmetric_PackerTask_by_truncation(pose, task);
  protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
  packer->apply(pose);
}


void
minimize(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos, bool move_bb, bool move_sc, bool move_rb) {
	// Initialize a MoveMap
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(move_rb);
	movemap->set_bb(false);
	movemap->set_chi(false);
	
	// Set allowable move types at interface positions
	// Currently, only sc moves allowed
	for (utility::vector1<Size>::iterator i = design_pos.begin(); i != design_pos.end(); i++) {
		movemap->set_bb (*i, move_bb);
		movemap->set_chi(*i, move_sc);
	}
	
	// Make MoveMap symmetric, apply it to minimize the pose
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
	m.apply(pose);
}

///////////////////////////////////////////////////////////////////////////////

class RBSymmSampler : public protocols::moves::Mover {
private:
	core::Real ang_step_, ang_max_;
	core::Real alpha_;
	core::Real trans_min_, trans_step_, trans_max_;
	core::Real fav_nat_bonus_;

	core::Real cen_vdw_filter_;
	core::Real sa_filter_1_, sa_filter_2_, sa_filter_3_;
	core::Real sc_filter_1_, sc_filter_2_, sc_filter_3_;
	core::Real ddg_filter_;
	core::Real air_filter_;

	bool cen_min_, smallres_only_;

	core::scoring::ScoreFunctionOP cen_scorefxn_, cenmin_scorefxn_, fa_scorefxn_;

public:
	RBSymmSampler(){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		ang_step_   =  option[crystdes::rot_step]();
		ang_max_    =  option[crystdes::rot_max]();
		trans_step_ =  option[crystdes::trans_step]();
		trans_max_  =  option[crystdes::trans_max]();
		trans_min_  = -trans_max_;

		cen_min_ =  option[crystdes::cen_min]();
		smallres_only_ = option[crystdes::smallres_only]();
		alpha_ = option[crystdes::alpha]();

		// filters
		cen_vdw_filter_ = option[crystdes::vdw]();
		sa_filter_1_    = option[crystdes::sa1]();
		sa_filter_2_    = option[crystdes::sa2]();
		sa_filter_3_    = option[crystdes::sa3]();
		sc_filter_1_    = option[crystdes::sc1]();
		sc_filter_2_    = option[crystdes::sc2]();
		sc_filter_3_    = option[crystdes::sc3]();
		ddg_filter_     = option[crystdes::ddg]();
		air_filter_     = option[crystdes::air]();
		
		cen_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("score0");
		cenmin_scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction();
		cenmin_scorefxn_->set_weight( core::scoring::vdw, 0.5 );
		cenmin_scorefxn_->set_weight( core::scoring::cbeta_smooth, 1.0 );
		cenmin_scorefxn_->set_weight( core::scoring::cenpack_smooth, 1.0 );
		cenmin_scorefxn_->set_weight( core::scoring::cen_env_smooth, 1.0 );
		cenmin_scorefxn_->set_weight( core::scoring::hbond_lr_bb, 1.0 );
		cenmin_scorefxn_->set_weight( core::scoring::hbond_sr_bb, 1.0 );

		fa_scorefxn_  = core::scoring::getScoreFunction();

		// Get the favor_native_residue weight from the command line
		fav_nat_bonus_ = 0.0;
		if (option[matdes::design::fav_nat_bonus]()) {
	      fav_nat_bonus_ = option[matdes::design::fav_nat_bonus]();
		}

		if (option[crystdes::print_des]()) {
			trans_min_ = trans_max_ = ang_max_ = alpha_ = 0;
			cen_vdw_filter_ = 99999;
		}
	}

	virtual std::string get_name() const { return "RBSymmSampler"; }

	// helper func
	//    euler angles (degrees) to rotation matrix
	void euler2rot( core::Real a, core::Real b, core::Real g, numeric::xyzMatrix<core::Real> &R) {
		const core::Real deg2rad = numeric::NumericTraits<Real>::pi()/180.0;
		R.xx() = -sin(deg2rad*a)*cos(deg2rad*b)*sin(deg2rad*g) + cos(deg2rad*a)*cos(deg2rad*g);
		R.xy() =  cos(deg2rad*a)*cos(deg2rad*b)*sin(deg2rad*g) + sin(deg2rad*a)*cos(deg2rad*g);
		R.xz() =  sin(deg2rad*b)*sin(deg2rad*g);
		R.yx() = -sin(deg2rad*a)*cos(deg2rad*b)*cos(deg2rad*g) - cos(deg2rad*a)*sin(deg2rad*g);
		R.yy() =  cos(deg2rad*a)*cos(deg2rad*b)*cos(deg2rad*g) - sin(deg2rad*a)*sin(deg2rad*g);
		R.yz() =  sin(deg2rad*b)*cos(deg2rad*g);
		R.zx() =  sin(deg2rad*a)*sin(deg2rad*b);
		R.zy() = -cos(deg2rad*a)*sin(deg2rad*b);
		R.zz() =  cos(deg2rad*b);
	}

	//
	void apply( Pose & pose) {
		using namespace core;
		using namespace core::conformation;
		using namespace core::conformation::symmetry;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		utility::vector1<id::AtomID> atoms;

		SymmetricConformation const & symm_conf (
		      dynamic_cast<SymmetricConformation const & > ( pose.conformation() ) );

		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
		Size nres_monomer = symm_info->num_independent_residues();

		for (core::Size resid=1; resid<=symm_info->num_total_residues_without_pseudo(); ++resid) {
			if (!symm_info->bb_is_independent(resid)) continue;
			core::conformation::Residue const &rsd_i = pose.residue(resid);
			for (core::Size atmid=1; atmid<=rsd_i.natoms(); ++atmid) {
				atoms.push_back( id::AtomID(atmid,resid) );
			}
		}

		Pose mono;
		core::pose::symmetry::extract_asymmetric_unit(pose, mono);
		utility::vector1<Real> sc_sasa = sidechain_sasa(mono,2.5);

		// steal coords from ASU
		utility::vector1< numeric::xyzVector<core::Real> > xyzs( atoms.size() );
		pose.batch_get_xyz( atoms, xyzs );
		utility::vector1< numeric::xyzVector<core::Real> > xyz_xform=xyzs;
		core::Size natoms = xyzs.size();

		// get CoM
		numeric::xyzVector<core::Real> com(0,0,0);
		for (int i=1; i<=natoms; ++i)
			com += xyzs[i];
		com /= (core::Real)natoms;

		protocols::moves::MoverOP tocentroid( new protocols::simple_moves::SwitchResidueTypeSetMover("centroid") );
		protocols::moves::MoverOP tofa( new protocols::simple_moves::SwitchResidueTypeSetMover("fa_standard") );

		core::io::silent::SilentFileData sfd;	

		// loop over all positions
		core::Real gamma=0;  // fpd no need to sample both alpha & gamma for small values of beta
		core::Real alpha=alpha_;
		for (core::Real beta=0; beta<=ang_max_; beta+=ang_step_)
		for (core::Real x=trans_min_; x<=trans_max_; x+=trans_step_)
		for (core::Real y=trans_min_; y<=trans_max_; y+=trans_step_)
		for (core::Real z=trans_min_; z<=trans_max_; z+=trans_step_) {
			numeric::xyzMatrix<core::Real> R_abg;
			euler2rot(alpha,beta,gamma, R_abg);
			numeric::xyzVector<core::Real> T_xyz(x,y,z);

			// apply transform
			for (int i=1; i<=natoms; ++i)
				xyz_xform[i] = R_abg*(xyzs[i]-com)+T_xyz+com;

			pose.batch_set_xyz( atoms, xyz_xform );

			// CENTROID FILTER
			core::pose::Pose cenpose = pose;
			tocentroid->apply( cenpose );

			// opt centroid rigid-body min
			if (cen_min_) {
				core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
				movemap->set_jump(true);movemap->set_bb(false);movemap->set_chi(false);
				core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
				protocols::simple_moves::symmetry::SymMinMover m( movemap, cenmin_scorefxn_, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
				m.apply(cenpose);
			}

			core::Real cen_score = (*cen_scorefxn_)(cenpose);
			if (cen_score > cen_vdw_filter_) {
				TR << alpha<<"_"<<beta<<"_"<<x<<"_"<<y<<"_"<<z << "    fail centroid filter (" << cen_score << ")" << std::endl;
				continue;
			}

			// steal sidechains
			tofa->apply( cenpose );
			for ( Size ii = 1; ii <= nres_monomer; ++ii ) {
				core::conformation::ResidueOP new_residue( pose.residue( ii ).clone() );
				cenpose.replace_residue ( ii, *new_residue, true );
			}
			pose = cenpose;

			// DESIGN STEP
			// 1 find interface reses
			// Find out which positions are near the inter-subunit interfaces
			// These will be further screened below, then passed to design() and minimize()
			Pose pose_for_design = pose;
			Real const contact_dist = option[matdes::design::contact_dist]();
			Real const contact_dist_sq = contact_dist * contact_dist;
			vector1<bool> indy_resis = symm_info->independent_residues();
			Sizes interface_pos;
			for (Size ir=1; ir<=symm_info->num_total_residues_without_pseudo(); ir++) {
				if (!indy_resis[ir]) continue;
				std::string atom_i = "";
				atom_i = (pose_for_design.residue(ir).name3() == "GLY")? "CA":"CB";
				for (Size jr=1; jr<=symm_info->num_total_residues_without_pseudo(); jr++) {
					if (indy_resis[jr]) continue;
					std::string atom_j = "";
					atom_j = (pose_for_design.residue(jr).name3() == "GLY")? "CA":"CB";
					if (pose_for_design.residue(ir).xyz(atom_i).distance_squared(pose_for_design.residue(jr).xyz(atom_j)) <= contact_dist_sq) {
					  interface_pos.push_back(ir);
					  break;
					}
				}
			}

			// Finally, filter the positions for design based on surface accessibility.
			// We don't want to design positions that are near to the interface but pointed
			// in toward the core of the monomer (e.g., positions on the insides of helices).
			// At the same time, create ResidueTypeConstraints favoring the native residue
			// at each design position if a fav_nat_bonus option is passed.
			Sizes design_pos;
			utility::vector1<core::scoring::constraints::ConstraintOP> favor_native_constraints;
			favor_native_constraints.clear();
			for(Size i = 1; i <= interface_pos.size(); ++i) {
				if( sc_sasa[interface_pos[i]] > 0.0 ) {
					design_pos.push_back(interface_pos[i]);
					if (fav_nat_bonus_ != 0.0) {
						core::scoring::constraints::ConstraintOP resconstraint
						     = new core::scoring::constraints::ResidueTypeConstraint(pose_for_design, interface_pos[i], fav_nat_bonus_);
						favor_native_constraints.push_back(resconstraint);
						//resconstraint->show(TR);
						//TR << std::endl;
					}
				}
			}

			if ( option[crystdes::print_des]() ) {
				std::cout << "design_pos " << design_pos[1];
				for (Size index=2; index<=design_pos.size(); index++) {
					std::cout << "," << design_pos[index];
				}
				std::cout << std::endl;
				return; // output only
			}

			if (fav_nat_bonus_ != 0.0) {
				pose_for_design.add_constraints(favor_native_constraints);
			}

			// Design
			design(pose_for_design, fa_scorefxn_, design_pos, smallres_only_);
			minimize(pose_for_design, fa_scorefxn_, design_pos, false, true, true);
			design(pose_for_design, fa_scorefxn_, design_pos, smallres_only_);

			// Repack and minimize using scorefxn
			ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
			repack(pose_for_design, scorefxn, design_pos);
			minimize(pose_for_design, scorefxn, design_pos, false, true, false);
			scorefxn->score(pose_for_design);

			// Calculate the surface area and surface complementarity for the interface
			Real sa1=0,sc1=0,sa2=0,sc2=0,sa3=0,sc3=0;
			new_sc(pose_for_design, sa1,sc1,sa2,sc2,sa3,sc3);

			// Calculate the ddG of the monomer in the assembled and unassembled states
			protocols::simple_moves::ddG ddG_mover = protocols::simple_moves::ddG(scorefxn, 1, true);
			ddG_mover.calculate(pose_for_design);
			Real ddG = ddG_mover.sum_ddG();

			// Calculate per-residue energies for interface residues
			Real interface_energy = 0;
			core::scoring::EnergyMap em;
			Real avg_interface_energy = 0;
			for (Size index=1; index<=design_pos.size(); index++) {
				interface_energy += pose_for_design.energies().residue_total_energy(design_pos[index]);
				em += pose_for_design.energies().residue_total_energies(design_pos[index]);
			}
			avg_interface_energy = interface_energy / design_pos.size();
			em *= fa_scorefxn_->weights();

			// fa filters
			if (sa1 < sa_filter_1_) {
				TR << alpha<<"_"<<beta<<"_"<<x<<"_"<<y<<"_"<<z << "    fail sa1 filter (" << sa1 << ")" << std::endl;
				continue;
			}
			if (sa2 < sa_filter_2_) {
				TR << alpha<<"_"<<beta<<"_"<<x<<"_"<<y<<"_"<<z << "    fail sa2 filter (" << sa2 << ")" << std::endl;
				continue;
			}
			if (sa3 < sa_filter_3_) {
				TR << alpha<<"_"<<beta<<"_"<<x<<"_"<<y<<"_"<<z << "    fail sa3 filter (" << sa3 << ")" << std::endl;
				continue;
			}
			if (sc1 < sc_filter_1_) {
				TR << alpha<<"_"<<beta<<"_"<<x<<"_"<<y<<"_"<<z << "    fail sc1 filter (" << sc1 << ")" << std::endl;
				continue;
			}
			if (sc2 < sc_filter_2_) {
				TR << alpha<<"_"<<beta<<"_"<<x<<"_"<<y<<"_"<<z << "    fail sc2 filter (" << sc2 << ")" << std::endl;
				continue;
			}
			if (sc3 < sc_filter_3_) {
				TR << alpha<<"_"<<beta<<"_"<<x<<"_"<<y<<"_"<<z << "    fail sc3 filter (" << sc3 << ")" << std::endl;
				continue;
			}
			if (ddG > ddg_filter_) {
				TR << alpha<<"_"<<beta<<"_"<<x<<"_"<<y<<"_"<<z << "    fail ddG filter (" << ddG << ")" << std::endl;
				continue;
			}
			if (avg_interface_energy > air_filter_) {
				TR << alpha<<"_"<<beta<<"_"<<x<<"_"<<y<<"_"<<z << "    fail air filter (" << avg_interface_energy << ")" << std::endl;
				continue;
			}

			// finally ... dump the output struct
			TR << alpha<<"_"<<beta<<"_"<<x<<"_"<<y<<"_"<<z << "    PASSES ALL FILTERS!" << std::endl;


			std::string tag = string_of(numeric::random::uniform()).substr(2,4);
			std::string fn = option[crystdes::prefix]()+"_"+string_of(alpha)+"_"+string_of(beta)+"_"+string_of(x)+"_"+string_of(y)+"_"+string_of(z)+"_"+tag+"_final.pdb.gz";

			// Write the pdb file of the design
			utility::io::ozstream out( option[out::file::o]() + "/" + fn );
			pose_for_design.dump_pdb(out);
			core::io::pdb::extract_scores(pose_for_design,out);
			out.close();

			// Create a scorefile struct, add custom metrics to it
			core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
			ss_out->fill_struct(pose_for_design,fn);
			ss_out->add_energy("ddG", ddG);
			ss_out->add_energy("air_energy", avg_interface_energy);
			ss_out->add_energy("air_fa_atr", em[core::scoring::fa_atr] / design_pos.size());
			ss_out->add_energy("air_fa_rep", em[core::scoring::fa_rep] / design_pos.size());
			ss_out->add_energy("des_pos", design_pos.size());
			ss_out->add_energy("sa1", sa1);
			ss_out->add_energy("sc1", sc1);
			ss_out->add_energy("sa2", sa2);
			ss_out->add_energy("sc2", sc2);
			ss_out->add_energy("sa3", sa3);
			ss_out->add_energy("sc3", sc3);

			// Write the scorefile
			sfd.write_silent_struct( *ss_out, option[out::file::o]() + "/" + option[ out::file::silent ]() );
		}
	}
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols::simple_moves::symmetry;

	SequenceMoverOP seq( new SequenceMover() );
	seq->add_mover( new SetupForSymmetryMover() );
	seq->add_mover( new RBSymmSampler() );

	protocols::jd2::JobDistributor::get_instance()->go( seq );

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
    try {
	NEW_OPT(crystdes::prefix, "out-prefix", "S");
	NEW_OPT(crystdes::cen_min, "cen_min", true);
	NEW_OPT(crystdes::smallres_only, "smallres_only", true);
	NEW_OPT(crystdes::alpha, "alpha", 0.0);
	NEW_OPT(crystdes::rot_step, "rot_step", 1.5);
	NEW_OPT(crystdes::trans_step, "trans_step", 0.75);
	NEW_OPT(crystdes::rot_max, "rot_max", 7.5);
	NEW_OPT(crystdes::trans_max, "trans_max", 1.5);
	NEW_OPT(crystdes::vdw, "vdw", 1);
	NEW_OPT(crystdes::sa1, "sa1", 150.0);
	NEW_OPT(crystdes::sa2, "sa2", 150.0);
	NEW_OPT(crystdes::sa3, "sa3", 0.0);
	NEW_OPT(crystdes::sc1, "sc1", 0.5);
	NEW_OPT(crystdes::sc2, "sc2", 0.5);
	NEW_OPT(crystdes::sc3, "sc3", 0.0);
	NEW_OPT(crystdes::ddg, "ddg", -10);
	NEW_OPT(crystdes::air, "air", -1);
	NEW_OPT(crystdes::print_des, "print_des", false);

	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
                                  }
        return 0;
}
