// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/coiled_coil.cc
/// @brief samples coiled coils

#include <basic/database/open.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/VirtualCoordinate.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/flxbb/FlxbbDesign.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/symmetry/SymDockingInitialPerturbation.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
// #include <devel/init.hh>

#include <apps/pilot/will/mynamespaces.ihh>

typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

static basic::Tracer TR("smhybrid");

static numeric::random::RandomGenerator RG(60542);

string & replace_string(string & s, string const & f, string const & r) {
	size_t i = s.find(f);
	bool replaced = false;
	while( i != s.npos ) {
		// TR << i << " " << s.size() << " " << s.npos << std::endl;
		s = s.replace(i,f.size(),r);
		i = s.find(f,i);
		replaced = true;
	}
	if(!replaced) utility_exit_with_message("no string replaced!!! probably mismatch tags in input and sym file");
	// TR << i << " " << s.size() << " " << s.npos << std::endl;
	return s;
}

core::pose::Pose append_jumpless_pose(core::pose::Pose const & pose1, core::pose::Pose const & pose2) {
	core::pose::Pose pose(pose1);
	pose.append_residue_by_jump(pose2.residue(1),pose.n_residue());
	for(Size i = 2; i <= pose2.n_residue(); ++i ) {
		pose.append_residue_by_bond(pose2.residue(i));
	}
	return pose;
}

string process_ss_str(string in) {
	char prevchar = NULL;
	string outss = "";
	Size i = 0;
	do {
		if('{'==in[i]) {
			++i;
			string beg,end;
			while( i < in.size() && ','!=in[i] ) beg+=in[i++];
			i++; // skip ','
			while( i < in.size() && '}'!=in[i] ) end+=in[i++];
			if( beg.size()==0 || end.size()==0 ) {
				utility_exit_with_message( "can't parse outss string '"+in+"'" );
			}
			Size ibeg = max(0,utility::string2int(beg));
			Size iend = utility::string2int(end);
			Size n = std::ceil(((Real)(iend-ibeg+1))*uniform()) - 1;
			for(Size j = 1; j <= ibeg+n; ++j) outss += prevchar;
			prevchar = NULL;
		} else {
			if( prevchar ) {
				TR << "append to outss " << prevchar << std::endl;
				outss += prevchar;
			}
			prevchar = in[i];
		}
		++i;
	} while( i < in.size() );
	if( prevchar && prevchar != '}' ) outss += prevchar;
	// std::cerr << "!!!!!!!!!!!!!!!!!! '" << in << "' '" << outss << "'" << std::endl;
	return outss;
}

int up_jump_tree(Pose & pose, int rsd) {
	for(int i = 1; i <= (int)pose.fold_tree().num_jump(); ++i) {
		if( pose.fold_tree().downstream_jump_residue(i)==rsd ) {
			return pose.fold_tree().upstream_jump_residue(i);
		}
	}
	return 0;
}


// void addcc(core::pose::Pose & pose, core::id::AtomID aid, Size anchor, core::Real mult = 1.0 ) {
// 	core::id::StubID sid(AtomID(1,anchor),AtomID(2,anchor),AtomID(3,anchor));
// 	core::scoring::constraints::ConstraintOP cc = new core::scoring::constraints::LocalCoordinateConstraint(
// 		aid, sid, pose.xyz(aid), new core::scoring::constraints::HarmonicFunc(0,mult) );
// 	pose.add_constraint(cc);
// }

void add_apc(core::pose::Pose & pose, core::id::AtomID aid1, core::id::AtomID aid2, core::Real mean, core::Real sd,
			 core::scoring::ScoreType st = core::scoring::atom_pair_constraint)
{
	core::scoring::constraints::ConstraintOP cc = new core::scoring::constraints::AtomPairConstraint(
		aid1, aid2, new core::scoring::constraints::HarmonicFunc(mean,sd), st );
	pose.add_constraint(cc);
}

void
add_agc(
	core::pose::Pose & pose,
	core::id::AtomID aid1,
	core::id::AtomID aid2,
	core::id::AtomID aid3,
	core::Real mean,
	core::Real sd
) {
	core::scoring::constraints::ConstraintOP cc = new core::scoring::constraints::AngleConstraint(
		aid1, aid2, aid3, new core::scoring::constraints::HarmonicFunc(mean,sd) );
	pose.add_constraint(cc);
}


void change_floating_sc_geometry(Pose & pose, Size rsd_in, Size nres) {
	Pose before = pose;
	for(Size rsd = rsd_in; rsd <= pose.n_residue(); rsd+=nres) {
		Vec metal;
		for(Size j = 1; j <= pose.fold_tree().num_jump(); ++j) {
			if(pose.fold_tree().downstream_jump_residue(j) == (int)rsd) {
				metal = pose.xyz(AtomID(1,pose.fold_tree().upstream_jump_residue(j)));
				break;
			}
		}
		if( pose.residue(rsd).name3()=="HSC" ) {
			Size oldath=pose.residue(rsd).atom_index("NE2"), newath=pose.residue(rsd).atom_index("ND1"), ce1=pose.residue(rsd).atom_index("CE1");
			string newatm = "ND1";
			if( pose.xyz(AtomID(newath,rsd)).distance(metal) < pose.xyz(AtomID(oldath,rsd)).distance(metal) ) {
				newath=pose.residue(rsd).atom_index("NE2");
				oldath=pose.residue(rsd).atom_index("ND1");
				newatm = "NE2";
			}
			Vec axis = (pose.xyz(AtomID(oldath,rsd))-pose.xyz(AtomID(ce1,rsd))).cross(pose.xyz(AtomID(newath,rsd))-pose.xyz(AtomID(ce1,rsd)));
			Vec newcen = pose.xyz(AtomID(newath,rsd));
			Vec oldcen = pose.xyz(AtomID(oldath,rsd));
			Mat m = rotation_matrix_degrees(axis,142.3);
			for(Size i = 1; i <= pose.residue_type(rsd).natoms(); ++i) {
				// TR << "change_floating_sc_geometry_reverse: moving atom " << i << " " << pose.residue(rsd).atom_name(i) << " res " << rsd << std::endl;
				pose.set_xyz(AtomID(i,rsd),  m*(pose.xyz(AtomID(i,rsd))-newcen)+oldcen );
			}
			FoldTree ft = pose.fold_tree();
			for(Size j = 1; j <= pose.fold_tree().num_jump(); ++j) {
				if((pose.fold_tree().downstream_jump_residue(j)-1)%nres+1 == rsd && pose.fold_tree().upstream_atom(j).size()) {
					ft.set_jump_atoms(j,ft.upstream_atom(j),newatm);
				}
			}
			pose.fold_tree(ft);
		} else if( pose.residue(rsd).name3()=="BPY" ) {
			Vec axis = (pose.xyz(AtomID(pose.residue(rsd).atom_index("NE1"),rsd)) + pose.xyz(AtomID(pose.residue(rsd).atom_index("NN1"),rsd)))/2.0;
			axis = (axis-metal).normalized();
			Mat m = rotation_matrix_degrees(axis.normalized(),180.0);
			Vec oldnn1 = pose.xyz(AtomID(pose.residue(rsd).atom_index("NN1"),rsd));
			Vec oldne1 = pose.xyz(AtomID(pose.residue(rsd).atom_index("NE1"),rsd));
			for(Size i = 1; i <= pose.residue_type(rsd).natoms(); ++i) {
				// TR << "change_floating_sc_geometry_reverse: moving atom " << i << " " << pose.residue(rsd).atom_name(i) << " res " << rsd << std::endl;
				pose.set_xyz(AtomID(i,rsd),  m*(pose.xyz(AtomID(i,rsd))-metal)+metal );
			}
			if( std::fabs(pose.xyz(AtomID(pose.residue(rsd).atom_index("NE1"),rsd)).distance(oldnn1)) > 0.1 ||
			    std::fabs(pose.xyz(AtomID(pose.residue(rsd).atom_index("NN1"),rsd)).distance(oldne1)) > 0.1 )
			{
				before.dump_pdb("change_floating_sc_geometry_before.pdb");
				std::ofstream out("change_floating_sc_geometry_before.pdb",std::ios::app);
				out << "HETATM 9999 ZN    ZN A 999       " << F(5,3,metal.x()) << "   " << F(5,3,metal.y()) << "   " << F(5,3,metal.z()) << "  1.00  0.00              " << std::endl;
				out << "HETATM 9999 PT    PT A 999       " << F(5,3,(metal+axis).x()) << "   " << F(5,3,(metal+axis).y()) << "   " << F(5,3,(metal+axis).z()) << "  1.00  0.00              " << std::endl;
				pose  .dump_pdb("change_floating_sc_geometry_after.pdb");
				std::ofstream out2("change_floating_sc_geometry_after.pdb",std::ios::app);
				out2 << "HETATM 9999 ZN    ZN A 999       " << F(5,3,metal.x()) << "   " << F(5,3,metal.y()) << "   " << F(5,3,metal.z()) << "  1.00  0.00              " << std::endl;
				out2 << "HETATM 9999 PT    PT A 999       " << F(5,3,(metal+axis).x()) << "   " << F(5,3,(metal+axis).y()) << "   " << F(5,3,(metal+axis).z()) << "  1.00  0.00              " << std::endl;
				utility_exit_with_message("ERROR in change_floating_sc_geometry_reverse! BPY");
			}
			/// below only works for 2 around z
			// if(rsd != rsd_in) {
			// 	numeric::xyzVector<Real> xyz1 = pose.xyz(AtomID(1,rsd_in));
			// 	numeric::xyzVector<Real> xyz2 = pose.xyz(AtomID(1,rsd_in+nres));
			// 	if( fabs(xyz1.x()+xyz2.x()) > 0.1 || fabs(xyz1.y()+xyz2.y()) > 0.1 || fabs(xyz1.z()-xyz2.z()) > 0.1 ) {
			// 		before.dump_pdb("change_floating_sc_geometry_before.pdb");
			// 		std::ofstream out("change_floating_sc_geometry_before.pdb",std::ios::app);
			// 		out << "HETATM 9999 ZN    ZN A 999       " << F(5,3,metal.x()) << "   " << F(5,3,metal.y()) << "   " << F(5,3,metal.z()) << "  1.00  0.00              " << std::endl;
			// 		out << "HETATM 9999 PT    PT A 999       " << F(5,3,(metal+axis).x()) << "   " << F(5,3,(metal+axis).y()) << "   " << F(5,3,(metal+axis).z()) << "  1.00  0.00              " << std::endl;
			// 		pose  .dump_pdb("change_floating_sc_geometry_after.pdb");
			// 		std::ofstream out2("change_floating_sc_geometry_after.pdb",std::ios::app);
			// 		out2 << "HETATM 9999 ZN    ZN A 999       " << F(5,3,metal.x()) << "   " << F(5,3,metal.y()) << "   " << F(5,3,metal.z()) << "  1.00  0.00              " << std::endl;
			// 		out2 << "HETATM 9999 PT    PT A 999       " << F(5,3,(metal+axis).x()) << "   " << F(5,3,(metal+axis).y()) << "   " << F(5,3,(metal+axis).z()) << "  1.00  0.00              " << std::endl;
			// 		utility_exit_with_message("ERROR in change_floating_sc_geometry_reverse! BPY");
			// 	}
			// }
		}
	}

}

void change_floating_sc_geometry_reverse(Pose & pose, Size rsd, Size nres) {
	change_floating_sc_geometry(pose,rsd,nres);
}


class AbsFunc : public core::scoring::constraints::Func {
public:
	AbsFunc( Real const x0_in, Real const sd_in ): x0_( x0_in ), sd_( sd_in ){}

	core::scoring::constraints::FuncOP
	clone() const { return new AbsFunc( *this ); }

	Real func( Real const x ) const {
		Real const z = ( x-x0_ )/sd_;
		if(z < 0) return -z;
		else return z;
	}
	Real dfunc( Real const x ) const {
		if(x-x0_ < 0) return -1.0/sd_;
		else return 1.0/sd_;
	}

	void read_data( std::istream & in ){
		in >> x0_ >> sd_;
	}

	void show_definition( std::ostream &out ) const {
		out << "ABS " << x0_ << " " << sd_ << std::endl;
	}

	Real x0() const {
		return x0_;
	}

	Real sd() const {
		return sd_;
	}

	void x0( Real x ) {
		x0_ = x;
	}

	void sd( Real sd ) {
		sd_ = sd;
	}

	// Size
	// show_violations( std::ostream& out, Real x, Size verbose_level, core::Real threshold = 1 ) const;

private:
	Real x0_;
	Real sd_;
};

class RepFunc : public core::scoring::constraints::Func {
public:
	RepFunc( Real const range, Real const height ): range_( range )
	{
		height_ = height / (range_*range_);
	}

	core::scoring::constraints::FuncOP
	clone() const { return new RepFunc( *this ); }

	Real func( Real const x ) const {
		if(x >= range_) return 0;
		return (x-range_)*(x-range_)*height_;
	}
	Real dfunc( Real const x ) const {
		if(x >= range_) return 0;
		return 2*(x-range_)*height_;
	}

	void read_data( std::istream & in ){
		in >> range_ >> height_;
	}

	void show_definition( std::ostream &out ) const {
		out << "ABS " << range_ << " " << height_ << std::endl;
	}

private:
	Real range_,height_;
};

class CoupledFunc : public core::scoring::constraints::Func {
public:
	CoupledFunc( PoseCAP pose, AtomID a1, AtomID a2, AtomID b1, AtomID b2, core::scoring::constraints::FuncOP func )
	: pose_(pose),a1_(a1),a2_(a2),b1_(b1),b2_(b2),func_(func) {}

	core::scoring::constraints::FuncOP
	clone() const { return new CoupledFunc( *this ); }

	Real dist() const {
		return fabs(pose_->xyz(a1_).distance(pose_->xyz(a2_)) - pose_->xyz(b1_).distance(pose_->xyz(b2_)));
	}

	Real func( Real const /*x*/ ) const {
		// TR << "CoupledFunc " << pose_->xyz(a1_).distance(pose_->xyz(a2_)) << " " << pose_->xyz(b1_).distance(pose_->xyz(b2_)) << " " << dist() << std::endl;
		return func_->func(dist());
	}
	Real dfunc( Real const /*x*/ ) const {
		return func_->dfunc(dist());
	}

	void read_data( std::istream & /*in*/ ){
		utility_exit_with_message("CoupledFunc::read_data NOT IMPLEMENTED!");
		// in >> range_ >> height_;
	}

	void show_definition( std::ostream &out ) const {
		out << "Coupled " << a1_ << " " << a2_ << " " << b1_ << " " << b2_ << std::endl;
	}

private:
	PoseCAP pose_;
	AtomID a1_,a2_,b1_,b2_;
	core::scoring::constraints::FuncOP func_;
};




struct PoseWrap : public ReferenceCount {
	string tag_;
	core::pose::Pose pose,orig_pose,validate_reference_;
	Size nsub,nres,primary_subsub,root_atomno_;
	string ss_;
	vector1<core::pose::Pose> poses_;
	vector1<string> attach_atom_,ss_pref_,ss_suff_,attached_scattach_res_name_;
	vector1<Size> attach_rsd_,jmpu_,jmpd_,cuts_,user_cuts_,subsub_fixed_begin,subsub_fixed_end,floating_scs_,attach_as_sc_;
	vector1<Size> floating_scs_subsub_,floating_scs_sub_,design_res_user,fixed_res_user,frag_res_user,virtual_res_user,subsub_ends_,subsub_starts_;
	vector1<Size> attached_scattach_res_,linker_res_,rep_edge_res_;
	std::map<Size,Size> attach_to_region_;
	// vector1<Size> floating_scs_attach_;
	vector1<bool> is_repackable_, is_designable_, chainbreaks_;
	vector1<core::chemical::ResidueTypeCOP> res_types_;
	bool hascst, has_floating_sc;
	core::conformation::symmetry::SymmDataOP symm_data,symm_data_full;
	core::chemical::ResidueTypeSetCAP fa_residue_set_,cen_residue_set_;
	vector1<vector1<Size> > scattach_res_user,cst_sub_files_;
	bool debug;
	char ss(Size i) const { return ss_[i-1]; }
	string ss(Size i, Size l) const { return ss_.substr(i-1,l); }
	void ss(Size i, char c) { ss_[i-1] = c; }
	Size which_subsub(Size i) {
		i = (i-1)%nres+1;
		Size j = 1;
		while( j <= subsub_ends_.size() && i > subsub_ends_[j] ) j++;
		return j;
	}
	Size which_sub(Size i) {
		return( std::ceil(i/nres) );
	}

	PoseWrap(
		string tag_in,
		vector1<core::pose::Pose> const & poses,
		vector1<Size> const & rel_attach,
		vector1<string> const & attach_atom,
		vector1<string> const & ss_pref,
		vector1<string> const & ss_suff,
		vector1<string> const & res_pref,
		vector1<string> const & res_suff,
		vector1<Size>        const & attach_as_sc_in,
		string symm_def_template,
		string symm_def_template_reduced,
		vector1<vector1<Size> > const & design_res_in ,//= vector1<vector1<Size> >(),
		vector1<vector1<Size> > const & rep_edge_res_in ,//= vector1<vector1<Size> >(),
		vector1<vector1<Size> > const & fixed_res_in,// = vector1<vector1<Size> >(),
		vector1<vector1<Size> > const & frag_res_in,// = vector1<vector1<Size> >(),
		vector1<vector1<Size> > const & virtual_res_in,// = vector1<vector1<Size> >(),
		vector1<vector1<Size> > const & scattach_res_in,// = vector1<vector1<Size> >(),
		vector1<vector1<xyzVector<Size> > > const & jumpcut,// = vector1<vector1<xyzVector<Size> > >(),
		vector1<vector1<Size> > const & cst_sub_files,// = vector1<vector1<Size> >()
		vector1<bool>           const & chainbreaks_in
	) : tag_(tag_in),poses_(poses), attach_atom_(attach_atom), ss_pref_(ss_pref), ss_suff_(ss_suff), hascst(false)
	{
		using namespace core;
		using namespace chemical;
		using namespace conformation;
		using namespace pose;

		if(poses.size()!=rel_attach .size()) utility_exit_with_message( "attach_rsd  must be specified for every pose" );
		if(poses.size()!=attach_atom.size()) utility_exit_with_message( "attach_atom must be specified for every pose" );
		if(poses.size()!=ss_pref    .size()) utility_exit_with_message( "ss_pref     must be specified for every pose" );
		if(poses.size()!=ss_suff    .size()) utility_exit_with_message( "ss_suff     must be specified for every pose" );

		scattach_res_user.push_back(vector1<Size>());

		primary_subsub = option[ basic::options::OptionKeys::smhybrid::primary_subsubunit ]();
		debug = option[basic::options::OptionKeys::smhybrid::debug]();

		for( Size i = 1; i <= poses.size(); ++i ) {
			if( attach_as_sc_in[i] && (ss_pref[i].size()>0 || ss_suff[i].size()>0) ) {
				utility_exit_with_message("can't have prefix or suffix on pose attached as side chain!! "+string_of(i)+" pref:  '"+ss_pref[i]+"' suff: '"+ss_suff[i]+"'");
			}
		}

		 fa_residue_set_ = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		cen_residue_set_ = ChemicalManager::get_instance()->residue_type_set( CENTROID    );

		has_floating_sc = false;
		// TR << "BLAH " << rel_attach.size() << " " << ss_pref.size() << " " << ss_suff.size() << std::endl;

		for(Size i = 1, e = poses_.size(); i <= e; ++i) {
			core::pose::Pose tmp(poses_[i]);
			// tmp.dump_pdb("tmp"+string_of(i)+".pdb");
			// set ss
			scoring::dssp::Dssp dssp(tmp);
			dssp.insert_ss_into_pose(tmp);
			// TR << "append pose " << rel_attach[i] << " '" << ss_pref[i] << "' '" << ss_suff[i] << "'" << std::endl;

			ss_ += ss_pref[i];
			scattach_res_user.push_back(vector1<Size>());
			for(Size j = 1; j <=   design_res_in[i].size(); ++j)   design_res_user   .push_back(  design_res_in[i][j]+ss_.size());
			for(Size j = 1; j <=    fixed_res_in[i].size(); ++j)    fixed_res_user   .push_back(   fixed_res_in[i][j]+ss_.size());
			for(Size j = 1; j <=     frag_res_in[i].size(); ++j)     frag_res_user   .push_back(    frag_res_in[i][j]+ss_.size());
			for(Size j = 1; j <=  virtual_res_in[i].size(); ++j)  virtual_res_user   .push_back( virtual_res_in[i][j]+ss_.size());
			for(Size j = 1; j <= scattach_res_in[i].size(); ++j) scattach_res_user[i].push_back(scattach_res_in[i][j]+ss_.size());

			// Size ins_frag_segment_edge = 3;
			bool marked_subsub_fixed = false;
			for(Size j=1;j<=tmp.n_residue();++j) {
				// if( frag_res_in[i].size() ) {
					if(std::find(frag_res_in[i].begin(),frag_res_in[i].end(),j)!=frag_res_in[i].end()) {
						ss_ += tmp.secstruct(j);
					} else {
						ss_ += 'F';
					}
					if(ss_[ss_.size()]=='F' || std::find(fixed_res_in[i].begin(),fixed_res_in[i].end(),j)!=fixed_res_in[i].end()) {
						if(!marked_subsub_fixed) {
							marked_subsub_fixed = true;
							subsub_fixed_begin.push_back(ss_.size());
							subsub_fixed_end  .push_back(ss_.size());
						}
						subsub_fixed_end.back() = ss_.size();
					}
			}
			if( ! marked_subsub_fixed ) {
				subsub_fixed_begin.push_back(0);
				subsub_fixed_end.push_back(0);
			}

			ss_ += ss_suff[i];
			TR << "subunit " << i << " ss: " << ss_ << std::endl;
			// create tmp pose with pref and suff
			if(tmp.residue(              1).is_protein()) remove_lower_terminus_type_from_pose_residue(tmp,1);
			if(tmp.residue(tmp.n_residue()).is_protein()) remove_upper_terminus_type_from_pose_residue(tmp,tmp.n_residue());
			for( Size j = 1; j <= ss_pref[i].size(); ++j ) {
				char c = 'L'; if(res_pref[i].size()==ss_pref[i].size()) c = res_pref[i][ss_pref.size()-i];
				string name3 = name_from_aa(aa_from_oneletter_code(c));
				tmp.prepend_polymer_residue_before_seqpos(*ResidueFactory::create_residue(fa_residue_set_->name_map(name3)),1,true);
				// tmp.set_phi  (1,-60); tmp.set_psi  (1,-45); tmp.set_omega(1,180);
			}
			for( Size j = 1; j <= ss_suff[i].size(); ++j ) {
				char c = 'L'; if(res_suff[i].size()==ss_suff[i].size()) c = res_suff[i][ss_suff.size()-i];
				string name3 = name_from_aa(aa_from_oneletter_code(c));
				tmp.append_residue_by_bond(*ResidueFactory::create_residue(fa_residue_set_->name_map(name3)),true);
				// tmp.set_phi  (tmp.n_residue(),-60); tmp.set_psi  (tmp.n_residue(),-45); tmp.set_omega(tmp.n_residue(),180);
			}
			// append to main pose
			attach_rsd_.push_back(pose.n_residue()+ss_pref[i].size()+rel_attach[i]);
			pose.append_residue_by_jump(tmp.residue(1),pose.n_residue());
			for(Size j = 2; j <= tmp.n_residue(); ++j ) {
				if(! tmp.residue(j).is_lower_terminus() ) {
					pose.append_residue_by_bond(tmp.residue(j));
				} else {
					pose.append_residue_by_jump(tmp.residue(j),pose.n_residue());
				}
			}

			for(Size j = 1; j <= rep_edge_res_in[i].size(); ++j) {
				rep_edge_res_.push_back( rep_edge_res_in[i][j] + pose.n_residue() - tmp.n_residue() + ss_pref[i].size() );
			}

			if(i != 1 ) cuts_.push_back( pose.n_residue() - tmp.n_residue() );
			if(i != 1 ) jmpu_.push_back(attach_rsd_[i-1]);
			if(i != 1 ) jmpd_.push_back(attach_rsd_[i  ]);
			if( jumpcut[i].size() > 0 ) {
				for( Size j = 1; j <= jumpcut[i].size(); ++j ) {
					jmpu_     .push_back( jumpcut[i][j].x() + pose.n_residue() - tmp.n_residue() + ss_pref[i].size() );
					jmpd_     .push_back( jumpcut[i][j].y() + pose.n_residue() - tmp.n_residue() + ss_pref[i].size() );
					cuts_     .push_back( jumpcut[i][j].z() + pose.n_residue() - tmp.n_residue() + ss_pref[i].size() );
					user_cuts_.push_back( jumpcut[i][j].z() + pose.n_residue() - tmp.n_residue() + ss_pref[i].size() );
					TR << "adding jmps/cuts from user input " << jmpu_[jmpu_.size()] << " " << jmpd_[jmpd_.size()] << " " << cuts_[cuts_.size()] << std::endl;
				}
			} else {
				FoldTree ft = poses_[i].fold_tree();
				for( Size j = 1; j <= ft.num_jump(); ++j ) {
					jmpu_     .push_back( ft.  upstream_jump_residue(j) + pose.n_residue() - tmp.n_residue() + ss_pref[i].size() );
					jmpd_     .push_back( ft.downstream_jump_residue(j) + pose.n_residue() - tmp.n_residue() + ss_pref[i].size() );
					cuts_     .push_back( ft.cutpoint_by_jump       (j) + pose.n_residue() - tmp.n_residue() + ss_pref[i].size() );
					TR << "adding jmps/cuts from orig subsub pose " << jmpu_[jmpu_.size()] << " " << jmpd_[jmpd_.size()] << " " << cuts_[cuts_.size()] << std::endl;
				}
			}
			for(Size j = 1; j <= ss_pref[i].size(); ++j) linker_res_.push_back(pose.n_residue()-tmp.n_residue()+j);
			for(Size j = ss_suff[i].size(); j >= 1; --j) linker_res_.push_back(pose.n_residue()-j+1);

			attach_as_sc_.push_back(attach_as_sc_in[i]);
			if(attach_as_sc_in[i] > 0) {
				floating_scs_.push_back(attach_rsd_[i]);
				if(attach_as_sc_in[i] > poses.size()) utility_exit_with_message("subsub "+string_of(attach_as_sc_in[i])+" doesn't exist!");
				floating_scs_subsub_.push_back(attach_as_sc_in[i]);
				has_floating_sc = true;
			}
			// pose.dump_pdb("PoseWrap_"+string_of(i)+".pdb");
			chainbreaks_  .push_back(chainbreaks_in[i]);
			subsub_ends_  .push_back(pose.n_residue());
			subsub_starts_.push_back(pose.n_residue()-tmp.n_residue()+1);
		}
		nres = pose.n_residue();
		// Size lt = 1, ut = pose.n_residue();
		// while( lt <= attach_as_sc_in.size() && attach_as_sc_in[lt] ) lt += poses[lt].n_residue(); // TODO: make more general
		// if( attach_as_sc_in[poses.size()] ) ut -=  poses[poses.size()].n_residue(); // TODO: make more general
		// add_lower_terminus_type_to_pose_residue(pose,lt);
		// add_upper_terminus_type_to_pose_residue(pose,ut);


		// TR << ss_.size() << " " << pose.n_residue() << std::endl;
		assert( ss_.size() == pose.n_residue() );

		if(debug) dump_pdb("init0.pdb");

		TR << "PoweWrap constructor attach:"; for(Size i=1;i<=attach_rsd_.size();++i) TR << " " << attach_rsd_[i]; TR << std::endl;
		TR << "jumps/cuts size " << jmpu_.size() << " " << jmpd_.size() << " " << cuts_.size() << std::endl;
		TR << "PoweWrap constructor jumps/cuts :"; for(Size i=1;i<=cuts_.size();++i) TR << jmpu_[i] << " " << jmpd_[i] << " " << cuts_[i] << ",  "; TR << std::endl;


		//////// add variants
		if( floating_scs_.size() > 1 ) {
			// utility_exit_with_message("only support one floating sc at the moment");
			TR << "WARNING: morethan one floating SC may be trouble!!!" << std::endl;
		}
		for(Size i = 1; i <= user_cuts_.size(); ++i) {
			// if not a floating SC and not a terminus, put cut points
			if( find(floating_scs_.begin(),floating_scs_.end(),user_cuts_[i])   != floating_scs_.end() ) continue;
			if( find(floating_scs_.begin(),floating_scs_.end(),user_cuts_[i]+1) != floating_scs_.end() ) continue;
			if( pose.residue(user_cuts_[i]).is_upper_terminus() ) continue;
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, user_cuts_[i]   );
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, user_cuts_[i]+1 );
		}
		for(Size i = 1; i <= chainbreaks_.size()-1; ++i) {
			if(chainbreaks_[i]) {
				if(pose.residue(subsub_ends_[i]).is_protein() && pose.residue(subsub_ends_[i]+1).is_protein() ) {
					core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, subsub_ends_[i]   );
					core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, subsub_ends_[i]+1 );
				} else {
					utility_exit_with_message("tried to add cutpoint to non-protein residue!");
				}
			} else {
				if(pose.residue(subsub_ends_[i]+1).is_protein()) add_lower_terminus_type_to_pose_residue(pose,subsub_starts_[i+1]);
				if(pose.residue(subsub_ends_[i]  ).is_protein()) add_upper_terminus_type_to_pose_residue(pose,subsub_ends_  [i  ]);
			}
		}
		if( chainbreaks_[chainbreaks_.size()] ) {
			if(pose.residue(pose.n_residue()).is_protein()) core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, pose.n_residue() );
			if(pose.residue(1               ).is_protein()) core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, 1    );
		} else {
			if(pose.residue(1               ).is_protein()) pose::add_lower_terminus_type_to_pose_residue(pose,1);
			if(pose.residue(pose.n_residue()).is_protein()) pose::add_upper_terminus_type_to_pose_residue(pose,pose.n_residue());
		}
		for( Size i = 1; i <= floating_scs_.size(); ++i ) {
			core::pose::add_variant_type_to_pose_residue(pose,"VIRTUAL_BBCB",floating_scs_[i]);
		}
		if(option[basic::options::OptionKeys::smhybrid::virtual_nterm]()) {
			if(debug) pose.dump_pdb("before_vnterm.pdb");
			// remove_lower_terminus_type_from_pose_residue(pose,1);
			if(!pose.residue(1).is_lower_terminus()) utility_exit_with_message("virtual_nterm requested, but res 1 isn't terminus");
			if(pose.residue(1).is_upper_terminus()) pose::remove_upper_terminus_type_from_pose_residue(pose,1);
			core::pose::add_variant_type_to_pose_residue(pose,"VIRTUAL_NTERM",1);
		}

		if(debug) dump_pdb("init1.pdb");

		//////////// setup fold tree
		// TR << "FT " << pose.fold_tree() << std::endl;
		std::map<Size,string> jump_atom;
		for(Size i = 1; i <= attach_rsd_.size(); ++i ) {
			jump_atom[attach_rsd_[i]] = attach_atom_[i];
		}
		ObjexxFCL::FArray2D_int jumps(2,cuts_.size());
		ObjexxFCL::FArray1D_int cuts(cuts_.size());
		for(Size i = 1; i <= cuts_.size(); ++i ) {
			cuts(   i) = cuts_[i];
			jumps(1,i) = jmpu_[i];//attach_rsd_[i];
			jumps(2,i) = jmpd_[i];//attach_rsd_[i+1];
		}
		kinematics::FoldTree ft = pose.fold_tree();
		ft.tree_from_jumps_and_cuts(pose.n_residue(),cuts_.size(),jumps,cuts,attach_rsd_[1]);
		TR << "tree_from_jumps_and_cuts:" << std::endl;
		TR << ft << std::endl;
		ft.reorder(attach_rsd_[primary_subsub]);
		std::vector<Size> seenit;
		for(Size i = 1; i <= ft.num_jump(); ++i ) {
			string up = "N", dn = "N";
			if( jump_atom.find(ft.  upstream_jump_residue(i)) != jump_atom.end() ) up = jump_atom[ft.  upstream_jump_residue(i)];
			if( jump_atom.find(ft.downstream_jump_residue(i)) != jump_atom.end() ) dn = jump_atom[ft.downstream_jump_residue(i)];
			ft.set_jump_atoms(i, up, dn );
			seenit.push_back(ft.downstream_jump_residue(i));
			TR << "removing " << ft.downstream_jump_residue(i) << " from root contention" << std::endl;
		}
		Size count = 0;
		for(std::map<Size,string>::iterator i = jump_atom.begin(); i != jump_atom.end(); ++i) {
			// TR << "POSEWRAP checking jump " << i->first << " " << i->second << std::endl;
			if(std::find(seenit.begin(),seenit.end(),i->first)!=seenit.end()) continue;
			if( jump_atom.find(i->first) == jump_atom.end() ) continue;
			TR << "POSEWRAP setting fold tree root " << i->first << " " << i->second << " " << pose.residue(i->first).atom_index(i->second) << std::endl;
			ft.reorder(i->first);
			// TR << "POSEWRAP setting fold tree root atomno " << i->first << " " << i->second << " " << pose.residue(i->first).atom_index(i->second) << std::endl;
			set_root_atomno( pose.residue(i->first).atom_index(i->second) );
			++count;
			if(count > 1) utility_exit_with_message( "should only be 1 anchor left as root after adding jumps" );
		}
		TR << "POSEWRAP done setting up fold tree" << std::endl;
		pose.fold_tree(ft);
		// TR << "FT " << pose.fold_tree() << std::endl;

		// std::exit(-1);

		TR << "expanding symm def template" << std::endl;
		utility::options::StringVectorOption & symm_file_tag = option[ basic::options::OptionKeys::smhybrid::symm_file_tag ];
		if( primary_subsub > poses.size() ) utility_exit_with_message("invalid primary subsub: " + string_of(primary_subsub) + " atch: " + string_of(attach_rsd_[primary_subsub]));
		if( attach_rsd_[primary_subsub] < 1                ) utility_exit_with_message("invalid primary subsub: " + string_of(primary_subsub) + " atch: " + string_of(attach_rsd_[primary_subsub]));
		if( attach_rsd_[primary_subsub] > pose.n_residue() ) utility_exit_with_message("invalid primary subsub: " + string_of(primary_subsub) + " atch: " + string_of(attach_rsd_[primary_subsub]));
		string pstr = symm_file_tag[primary_subsub];
		Size n_eq_primary = 0;
		for(Size i = 1; i <= symm_file_tag.size(); ++i) {
			if(symm_file_tag[i]==pstr) n_eq_primary++;
			if(symm_file_tag[i]!=pstr) symm_def_template = replace_string(symm_def_template,symm_file_tag[i],string_of(attach_rsd_[i]));
			else                       symm_def_template = replace_string(symm_def_template,symm_file_tag[i]," ");
			if(symm_file_tag[i]!=pstr) symm_def_template_reduced = replace_string(symm_def_template_reduced,symm_file_tag[i],string_of(attach_rsd_[i]));
			else                       symm_def_template_reduced = replace_string(symm_def_template_reduced,symm_file_tag[i]," ");
		}
		if(n_eq_primary!=1) {
			utility_exit_with_message("must be exactly one primary subunit... this is seriously messed up "+primary_subsub);
		}
		TR << "making symm data" << std::endl;

		symm_data = new core::conformation::symmetry::SymmData(pose.n_residue(),pose.fold_tree().num_jump());
		std::istringstream iss(symm_def_template_reduced);
		symm_data->read_symmetry_data_from_stream(iss);

		symm_data_full = new core::conformation::symmetry::SymmData(pose.n_residue(),pose.fold_tree().num_jump());
		std::istringstream iss2(symm_def_template);
		symm_data_full->read_symmetry_data_from_stream(iss2);

		orig_pose = pose;

		if(debug) dump_pdb("init2.pdb");

		// TR << "POSEWRAP switch_to_residue_type_set cen" << std::endl;
		// core::chemical::switch_to_residue_type_set( pose, core::chemical::CENTROID );

		// dump_pdb("init_cen.pdb");
		TR << "POSEWRAP make_symmetric_pose" << " " << std::endl;
		core::conformation::symmetry::make_symmetric_pose( pose, *symm_data );
		TR << "POSEWRAP make_symmetric_pose DONE" << std::endl;
		std::map<Size,core::conformation::symmetry::SymDof> dofs = symmetry::symmetry_info(pose)->get_dofs();
		std::map<Size,core::conformation::symmetry::SymDof>::iterator i;
		for( i = dofs.begin(); i != dofs.end(); ++i ) {
			TR << "DOF " << i->first << " " << i->second << std::endl;
		}


		// if( !option[basic::options::OptionKeys::smhybrid::refine]() ) {
		// 	TR << "POSEWRAP set symm dofs" << std::endl;
		// 	protocols::rigid::RigidBodyDofSeqRandomizeMover symsetup(dofs);
		// 	symsetup.apply(pose);
		// 	// for(Size i = 1; i <= 3; ++i)	change_floating_sc_geometry(pose,floating_scs_[i],nres);
		// 	// for(Size i = 1; i <= 2; ++i)	change_floating_sc_geometry(pose,floating_scs_[i],nres);
		// 	// for(Size i = 1; i <= 1; ++i)	change_floating_sc_geometry(pose,floating_scs_[i],nres);
		// 	for(Size i = 0; i < 30; i++) {
		// 		symsetup.apply(pose);
		// 		pose.dump_pdb("sym_init"+string_of(i)+".pdb");
		// 	}
		// 	// for(Size i = 1; i <= 3; ++i)	change_floating_sc_geometry(pose,floating_scs_[3],nres);
		// 	// // symsetup.apply(pose);
		// 	// pose.dump_pdb("sym_init6.pdb");
		// 	// symsetup.apply(pose);
		// 	// pose.dump_pdb("sym_init7.pdb");
		// 	// symsetup.apply(pose);
		// 	// pose.dump_pdb("sym_init8.pdb");
		// 	// symsetup.apply(pose);
		// 	// pose.dump_pdb("sym_init9.pdb");
		// 	std::exit(-1);
		// }

		nsub = symmetry::symmetry_info(pose)->subunits(); // TODO read from SymmData

		if(debug) dump_pdb("sym_init1.pdb");

		// if( !option[basic::options::OptionKeys::smhybrid::refine]() ) {
			// for(std::map<Size,core::conformation::symmetry::SymDof>::iterator i = dofs.begin(); i != dofs.end(); ++i) {
			// 	if( i->second.allow_dof(1) && i->second.allow_dof(2) && i->second.allow_dof(3) &&
			// 		i->second.allow_dof(4) && i->second.allow_dof(5) && i->second.allow_dof(6) ) continue; // ignore the "RB" 6 dof jumps
			// 	core::kinematics::Jump j(pose.jump(i->first));
			// 	j.set_translation( j.get_translation() + Vec(-40.0,0,0) );
			// 	pose.set_jump(i->first,j);
			// }

			TR << "POSEWRAP set phispsi" << std::endl;
			for(Size i = 1; i <= nres; ++i) {
				if(std::find(linker_res_.begin(),linker_res_.end(),i)==linker_res_.end()) continue;
				// if('F'==ss(i)) continue;
				pose.set_phi(i,-60); pose.set_psi(i,-45); pose.set_omega(i,180);
			}
		// }
		if(debug)pose.dump_pdb("sym_init2.pdb");
		switch_to_cen();
		if(debug) {
			pose.dump_pdb("sym_init_switch_to_cen1.pdb");
			switch_to_fa();
			pose.dump_pdb("sym_init_switch_to_fa.pdb");
			switch_to_cen();
			pose.dump_pdb("sym_init_switch_to_cen2.pdb");
			std::exit(-1);
		}


		// pose.dump_pdb("sym_init2.pdb");

		// pose.dump_pdb("sym_init3.pdb");
		// pose.dump_pdb("sym_init1.pdb");
		// pose.set_chi(1,1,pose.chi(1,1)+10.0);
		// pose.dump_pdb("sym_init2.pdb");
		//
		// std::exit(-1);


		// core::scoring::ScoreFunction sf;
		// sf.set_weight(scoring::atom_pair_constraint,1.0);
		// core::scoring::symmetry::SymmetricScoreFunction ssf(sf);
		// TR << ssf(pose) << std::endl;

		ft = pose.fold_tree();
		TR << ft << std::endl;
		for(Size j = 1; j <= ft.num_jump(); ++j ) {
			if( ft.downstream_jump_residue(j) > (int)nres && ft.upstream_jump_residue(j) > (int)nres ) continue;
			TR << "EDGE " << ft.jump_edge(j) << std::endl;
		}
		// TR << "add_constraints_from_cmdline_to_pose BEFORE NCST: " << pose.constraint_set()->get_all_constraints().size() << std::endl;
		core::scoring::constraints::add_constraints_from_cmdline_to_pose( pose );
		// core::scoring::constraints::ConstraintCOPs cs = pose.constraint_set()->get_all_constraints();
		// for(Size i = 1; i <= cs.size(); i++) {
		// 	TR << "cst " << i << std::endl;
		// 	cs[i]->show(TR);
		// }
		TR << "add_constraints_from_cmdline_to_pose AFTER  NCST: " << pose.constraint_set()->get_all_constraints().size() << std::endl;

		add_floating_sc_csts();

		cst_sub_files_ = cst_sub_files;
		// dump_pdb("cst_test.pdb");
		// std::exit(-1);


		if( basic::options::option[basic::options::OptionKeys::smhybrid::linker_cst]() ) {
			using namespace core::scoring::constraints;
			ConstraintCOPs ac;
			Size beg=1, end=subsub_starts_.size();
			while(beg < subsub_starts_.size() && attach_as_sc_[beg]) beg++;
			while(end > 0 && attach_as_sc_[end]) end--;
			if( beg < end ) {
				ac.push_back(new AtomPairConstraint( AtomID(1,subsub_starts_[beg]), AtomID(3,   nres              ), new AbsFunc(0,1), core::scoring::coordinate_constraint ) );
				ac.push_back(new AtomPairConstraint( AtomID(1,subsub_starts_[end]), AtomID(3,subsub_starts_[end]-1), new AbsFunc(0,1), core::scoring::coordinate_constraint ) );
				pose.add_constraint(new AmbiguousConstraint(ac));
			}
		}
		if( basic::options::option[basic::options::OptionKeys::smhybrid::subsubs_attract]() ) {
			using namespace core::scoring::constraints;
			vector1<Size> real_subsubs;
			for(Size i = 1; i <= attach_as_sc_.size(); ++i) if(!attach_as_sc_[i]) real_subsubs.push_back(i);
			for(Size i = 1; i <= real_subsubs.size(); ++i) {
				for(Size j = 1; j <= real_subsubs.size(); ++j) {
					if(   j > i) {
						Size ires=0,jres=0;
						for(Size jmp = 1; jmp <= pose.fold_tree().num_jump(); ++jmp) {
							if(pose.fold_tree().downstream_jump_residue(jmp)==(int)attach_rsd_[real_subsubs[i]]) ires = pose.fold_tree().upstream_jump_residue(jmp);
							if(pose.fold_tree().downstream_jump_residue(jmp)==(int)attach_rsd_[real_subsubs[j]]) jres = pose.fold_tree().upstream_jump_residue(jmp);
						}
						if( ires * jres == 0 ) utility_exit_with_message("ires or jres is 0!!!");
						pose.add_constraint(new AtomPairConstraint(AtomID(1,ires),AtomID(1,jres),new AbsFunc(0,1),core::scoring::coordinate_constraint ));
						TR << "adding subsub's attract cst sub 1-" << i << " 2-" << j << " " << ires << " " << jres << std::endl;
					}
					if(nsub > 1) {
						Size ires=0,jres=0;
						for(Size jmp = 1; jmp <= pose.fold_tree().num_jump(); ++jmp) {
							if(pose.fold_tree().downstream_jump_residue(jmp)==(int)attach_rsd_[real_subsubs[i]]          ) ires = pose.fold_tree().upstream_jump_residue(jmp);
							if(pose.fold_tree().downstream_jump_residue(jmp)==(int)attach_rsd_[real_subsubs[j]]+(int)nres) jres = pose.fold_tree().upstream_jump_residue(jmp);
						}
						if( ires * jres == 0 ) utility_exit_with_message("ires or jres is 0!!!");
						pose.add_constraint(new AtomPairConstraint(AtomID(1,ires),AtomID(1,jres),new AbsFunc(0,1),core::scoring::coordinate_constraint ));
						TR << "adding subsub's attract cst sub 2-" << i << " 2-" << j << " " << ires << " " << jres << std::endl;
					}
				}
			}
		}
		// pose.add_constraint(new AtomPairConstraint(AtomID(1,attach_rsd_[5]),AtomID(1,attach_rsd_[6]     ),new AbsFunc(0,1) ));
		// pose.add_constraint(new AtomPairConstraint(AtomID(1,attach_rsd_[5]),AtomID(1,attach_rsd_[5]+nres),new AbsFunc(0,1) ));
		// pose.add_constraint(new AtomPairConstraint(AtomID(1,attach_rsd_[5]),AtomID(1,attach_rsd_[6]+nres),new AbsFunc(0,1) ));
		// pose.add_constraint(new AtomPairConstraint(AtomID(1,attach_rsd_[6]),AtomID(1,attach_rsd_[5]+nres),new AbsFunc(0,1) ));
		// pose.add_constraint(new AtomPairConstraint(AtomID(1,attach_rsd_[6]),AtomID(1,attach_rsd_[6]+nres),new AbsFunc(0,1) ));

		if( basic::options::option[basic::options::OptionKeys::smhybrid::pseudosym]() ) {
			if( pose.n_residue() < 2*nres ) {
				utility_exit_with_message("pseudosym requires at least 2 subunits!!!");
			}
			using namespace core::scoring::constraints;
			FuncOP h = new HarmonicFunc(0,3.0);
			Size sp = 10;
			for(Size i = 1; i <=nres; i+=sp) {
				for(Size j = i+sp; j <=nres; j+=sp) {
					Size a1=i, a2=j+nres, b1=i+nres, b2=j;
					TR << "Pseudo-sym cst: " << a1 << " " << a2 << " " << b1 << " " << b2 << std::endl;
					CoupledFunc* cf = new CoupledFunc(&pose,AtomID(2,a1),AtomID(2,a2),AtomID(2,b1),AtomID(2,b2),h);
					pose.add_constraint(new AtomPairConstraint(AtomID(2,a1),AtomID(2,a2),cf));
					pose.add_constraint(new AtomPairConstraint(AtomID(2,b1),AtomID(2,b2),cf));
				}
			}
		}


		if( option[basic::options::OptionKeys::edensity::mapfile].user() ) {
			protocols::electron_density::SetupForDensityScoringMover m;
			m.apply(pose);
		}

		if(debug)pose.dump_pdb("sym_init4.pdb");

		// std::exit(-1);
		// for(Size i=1; i < 10; i+=1) {
		// 	for(Size j=1; j <= pose.fold_tree().num_jump(); ++j) {
		// 		if(pose.fold_tree().downstream_jump_residue(j)==1) {
		// 			core::kinematics::Jump jmp(pose.jump(j));
		// 			jmp.set_rotation(rotation_matrix_degrees(Vec(1,0,0),(gaussian()*10.0))*jmp.get_rotation());
		// 			pose.set_jump(j,jmp);
		// 		}
		// 		if(pose.fold_tree().downstream_jump_residue(j)==2) {
		// 			core::kinematics::Jump jmp(pose.jump(j));
		// 			jmp.set_rotation(rotation_matrix_degrees(Vec(1,0,0),(gaussian()*10.0))*jmp.get_rotation());
		// 			pose.set_jump(j,jmp);
		// 		}
		// 		// if(pose.fold_tree().downstream_jump_residue(j)==3) {
		// 		// 	core::kinematics::Jump jmp(pose.jump(j));
		// 		// 	jmp.set_rotation(rotation_matrix_degrees(Vec(1,0,0),(gaussian()*10.0))*jmp.get_rotation());
		// 		// 	pose.set_jump(j,jmp);
		// 		// }
		// 		// if(pose.fold_tree().downstream_jump_residue(j)==4) {
		// 		// 	core::kinematics::Jump jmp(pose.jump(j));
		// 		// 	jmp.set_rotation(rotation_matrix_degrees(Vec(1,0,0),(gaussian()*10.0))*jmp.get_rotation());
		// 		// 	pose.set_jump(j,jmp);
		// 		// }
		// 	}
		// 	dump_pdb("sym_init_"+string_of(i)+".pdb");
		// }
		// std::exit(-1);
		// TR << "score after adding csts: " << ssf(pose) << std::endl;
		//
		//
		// using namespace core::scoring::constraints;
		// ConstraintCOPs cs = pose.constraint_set()->get_all_constraints();
		// TR << "num cst: " << cs.size() << std::endl;
		// for(Size i = 1; i <= cs.size(); i++) {
		// 	TR << "cst " << i << std::endl;
		// 	cs[i]->show(TR);
		// }
		// // std::exit(-1);
		//
		// // TR << "constraints " <<  std::endl;
		// // pose.constraint_set()->show(TR);
		//
		// pose.set_xyz(AtomID(1,scattach_res_user[1]),pose.xyz(AtomID(1,floating_scs_[1])));
		// pose.set_xyz(AtomID(2,scattach_res_user[1]),pose.xyz(AtomID(2,floating_scs_[1])));
		// pose.set_xyz(AtomID(3,scattach_res_user[1]),pose.xyz(AtomID(3,floating_scs_[1])));
		// pose.set_xyz(AtomID(5,scattach_res_user[1]),pose.xyz(AtomID(5,floating_scs_[1])));
		//
		// for(Real i = -10; i <= 10; i+=1.0 ) {
		// 	pose.set_xyz(AtomID(1,scattach_res_user[1]),pose.xyz(AtomID(1,floating_scs_[1]))+Vec(i,0,0));
		// 	pose.set_xyz(AtomID(2,scattach_res_user[1]),pose.xyz(AtomID(2,floating_scs_[1]))+Vec(i,0,0));
		// 	pose.set_xyz(AtomID(3,scattach_res_user[1]),pose.xyz(AtomID(3,floating_scs_[1]))+Vec(i,0,0));
		// 	pose.set_xyz(AtomID(5,scattach_res_user[1]),pose.xyz(AtomID(5,floating_scs_[1]))+Vec(i,0,0));
		// 	TR << i << " " << ssf(pose) << std::endl;
		// }
		// ssf.show(pose);
		//
		// std::exit(-1);

		// dump_pdb("PoseWrap_symm_0.pdb");
		// for(Size iter = 1; iter < 10; iter++ ) {
		// 	for(std::map<Size,core::conformation::symmetry::SymDof>::iterator i = dofs.begin(); i != dofs.end(); ++i) {
		// 			core::kinematics::Jump j(pose.jump(i->first));
		// 			j.set_translation( j.get_translation() + Vec(gaussian()*4,0,0) );
		// 			j.set_rotation(numeric::rotation_matrix(Vec(1,0,0),gaussian())*j.get_rotation());
		// 			pose.set_jump(i->first,j);
		// 		}
		// 	// for(Size j = 1; j <= 4; ++j) {
		// 	// 	pose.set_chi( j, attach_rsd_[1], pose.chi(j,attach_rsd_[1]) + 10.0);
		// 	// }
		// 	dump_pdb("PoseWrap_symm_"+string_of(iter)+".pdb");
		// }
		//
		// std::exit(-1);

		// pose.dump_pdb("pw_init_done.pdb");

		if (option[ basic::options::OptionKeys::parser::view ]()) {
			protocols::viewer::add_conformation_viewer(pose.conformation(),"smhybrid",1000,1000);
		}


		// pose.dump_pdb("test0.pdb");
		// for(Size i = 1; i < 10; ++i) {
		// 	pose.set_chi(4,1,pose.chi(4,1)+10.0);
		// 	pose.dump_pdb("test"+string_of(i)+".pdb");
		// }
		// std::exit(-1);


	}

	void add_floating_sc_csts() {
		using namespace core::scoring::constraints;

		// TR << "NUM CST BEFORE: " << pose.constraint_set()->get_all_constraints().size() << std::endl;
		ConstraintCOPs cs = pose.constraint_set()->get_all_constraints();
		ConstraintCOPs rem;
		for(Size i = 1; i <= cs.size(); i++) {
			if( cs[i]->type() != "AmbiguousConstraint" ) continue;
			rem.push_back(cs[i]);
		}
		pose.remove_constraints(rem);
		// TR << "NUM CST AFTER: " << pose.constraint_set()->get_all_constraints().size() << std::endl;

		FuncOP afunc1 = new core::scoring::constraints::HarmonicFunc(0.0,1.0/1.0);
		FuncOP afunc2 = new core::scoring::constraints::HarmonicFunc(0.0,1.0/2.0);
		FuncOP afunc3 = new core::scoring::constraints::HarmonicFunc(0.0,1.0/3.0);
		FuncOP afunc6 = new core::scoring::constraints::HarmonicFunc(0.0,1.0/6.0);
		// FuncOP afunc1 = new AbsFunc(0.0,1.0/1.0);
		// FuncOP afunc2 = new AbsFunc(0.0,1.0/2.0);
		// FuncOP afunc3 = new AbsFunc(0.0,1.0/3.0);
		// FuncOP afunc6 = new AbsFunc(0.0,1.0/6.0);
		// hfunc->show(TR);
		// Size last_subsub = 0;
		floating_scs_sub_ = basic::options::option[basic::options::OptionKeys::smhybrid::attach_as_sc_sub]();
		while(floating_scs_sub_.size() < floating_scs_sub_.size()) floating_scs_sub_.push_back(1);
		for(Size i = 1; i <= floating_scs_.size(); ++i) {
			TR << i << " Adding floating scs constraints " << floating_scs_[i] << " " << floating_scs_subsub_[i] << " " << scattach_res_user[floating_scs_subsub_[i]].size() << " " << floating_scs_sub_[i] << std::endl;
			vector1<Size> scattach_res = scattach_res_user[floating_scs_subsub_[i]];
			// if(!last_subsub) last_subsub = floating_scs_subsub_[i];
			// if( floating_scs_subsub_[i] != last_subsub && pose.n_residue() > 2*nres ) { // TODO FIX THIS HACK.. asusmes inter-subunit
			// 	TR << "shifting scattach res to next subunit" << std::endl;
			// 	for(Size j = 1; j <= scattach_res.size(); ++j ) {
			// 		scattach_res[j] = scattach_res[j] + nres;
			// 	}
			// 	floating_scs_sub_.push_back(2);
			// } else {
				// floating_scs_sub_.push_back(scattach_sub[i]);
			// }
			// last_subsub = floating_scs_subsub_[i];
			ConstraintCOPs ac;
			for(Size j = 1; j <= scattach_res.size(); ++j) {
				// TR << "scattach_res " << i << " " << j << " " << scattach_res[j] << std::endl;
				if( pose.residue(scattach_res[j]).name3() == "GLY" ) continue;
				if( pose.residue(scattach_res[j]).name3() == "PRO" ) continue;
				ConstraintCOPs mc;
				FuncOP bb,cb;
				if(pose.residue(floating_scs_[i]).name3()=="BPY") {
					bb = afunc1; cb = afunc3;
				} else {
					bb = afunc2; cb = afunc6;
				}
				mc.push_back(new AtomPairConstraint( AtomID(pose.residue(scattach_res[j]).atom_index( "N"),scattach_res[j]+(floating_scs_sub_[i]-1)*nres), AtomID(pose.residue(floating_scs_[i]).atom_index( "N"),floating_scs_[i]), bb ) );
				mc.push_back(new AtomPairConstraint( AtomID(pose.residue(scattach_res[j]).atom_index("CA"),scattach_res[j]+(floating_scs_sub_[i]-1)*nres), AtomID(pose.residue(floating_scs_[i]).atom_index("CA"),floating_scs_[i]), bb ) );
				mc.push_back(new AtomPairConstraint( AtomID(pose.residue(scattach_res[j]).atom_index( "C"),scattach_res[j]+(floating_scs_sub_[i]-1)*nres), AtomID(pose.residue(floating_scs_[i]).atom_index( "O"),floating_scs_[i]), bb ) );
				mc.push_back(new AtomPairConstraint( AtomID(pose.residue(scattach_res[j]).atom_index("CB"),scattach_res[j]+(floating_scs_sub_[i]-1)*nres), AtomID(pose.residue(floating_scs_[i]).atom_index("CB"),floating_scs_[i]), cb ) );
				ac.push_back(new MultiConstraint(mc));
				TR << "adding cst " << i << " " << j << " " << scattach_res[j]+(floating_scs_sub_[i]-1)*nres << std::endl;
			}
			if(ac.size()) pose.add_constraint(new AmbiguousConstraint(ac));
			else utility_exit_with_message("no floating sidechain attachments speficied!");
		}

		// for(Size i = 1; i <= floating_scs_.size(); ++i) {
		// 	for(Size j = 1; j <= floating_scs_.size(); ++j) {
		// 		for(Size k = 0; k < nsub; ++k) {
		// 			Size r1 = floating_scs_[i];
		// 			Size r2 = floating_scs_[j]+k*nres;
		// 			if(r1 >= r2) continue;
		// 			for(Size a1 = 1; a1 <= pose.residue(r1).nheavyatoms(); ++a1) {
		// 				for(Size a2 = 1; a2 <= pose.residue(r2).nheavyatoms(); ++a2) {
		// 					Real d = 4.0;
		// 					if( a1 < 5 && a2 < 5 ) d = 8.0;
		// 					pose.add_constraint(new core::scoring::constraints::AtomPairConstraint(AtomID(a1,r1),AtomID(a2,r2), new RepFunc(d,20.0) ), core::scoring::coordinate_constraint);
		// 				}
		// 			}
		// 		}
		// 	}
		// }
		// prevents res in same subsub from overlapping
		// if(basic::options::option[basic::options::OptionKeys::smhybrid::floating_scs_rep]() && nsub > 1) {
		for(Size k = 1; k <= nsub; ++k) {
			for(Size i = 1; i <= floating_scs_.size(); ++i) {
				Size s = 1;	if(k == 1) s = i+1;
				for(Size j = s; j <= floating_scs_.size(); ++j) {
					ConstraintCOPs ac;
					TR << "adding fsc CA rep cstres1=" << floating_scs_[i] << " cstres2=" << floating_scs_[j]+(k-1)*nres << " | i=" << i << " j=" << j << " k=" << k << " n=" << nres << std::endl;
					ac.push_back(new core::scoring::constraints::AtomPairConstraint(
					   AtomID(2,floating_scs_[i]),AtomID(2,floating_scs_[j]+(k-1)*nres),
					   new RepFunc(8.0,1000.0), core::scoring::coordinate_constraint ));
					pose.add_constraint(new AmbiguousConstraint(ac));
				}
			}
		}

		// }
		// std::exit(-1);
		// TR << "BPY CB IS: " << pose.residue(1).atom_index("CB") << std::endl;
		// for(Size i = 1; i <= pose.residue(1).natoms(); ++i) {
		// 	TR << "BPY atom " << i << " " << pose.residue(1).atom_name(i) << std::endl;
		// }

	}

	void set_root_atomno(Size i) { root_atomno_ = i; }
	Size set_root_atomno() { return root_atomno_; }

	bool move_chi(Size rsd) {
		for(Size i = 1; i <= attach_rsd_.size(); ++i ) {
			if( rsd == attach_rsd_[i] && attach_atom_[i] != "N" && pose.residue(rsd).is_protein() ) return true;
		}
		return false;
	}


	void validate_set_reference() {
		validate_reference_ = pose;
	}
	bool validate() {
		check_scattach_res();
		return true;
	}

	string tag() {
		string r = tag_;
		for(Size i = 1; i <= attached_scattach_res_.size(); ++i) r += "_"+string_of(attached_scattach_res_[i]);
		return r;
	}
	// core::pose::Pose make_mono_pose( core::pose::Pose const & in ) const {
	// 	// create_subpose(
	// 	// 	Pose const & src,
	// 	// 	vector1< Size > const & positions,
	// 	// 	kinematics::FoldTree const & f,
	// 	// 	Pose & pose
	// 	// );
	//
	// 	ObjexxFCL::FArray2D_int jumps(2,cuts_.size());
	// 	ObjexxFCL::FArray1D_int cuts(cuts_.size());
	// 	for(Size i = 1; i <= cuts_.size(); ++i ) {
	// 		cuts(i) = cuts_[i];
	// 		jumps(1,i) = attach_rsd_[i];
	// 		jumps(2,i) = attach_rsd_[i+1];
	// 	}
	// 	core::kinematics::FoldTree ft = pose.fold_tree();
	// 	ft.tree_from_jumps_and_cuts(nres,cuts_.size(),jumps,cuts,attach_rsd_[1]);
	//
	// 	vector1<Size> pos;
	// 	for(Size i = 1; i <= nres; ++i) pos.push_back(i);
	//
	// 	core::pose::Pose pose;
	// 	core::pose::create_subpose(in,pos,ft,pose);
	// 	return pose;
	// }
	// void make_full_pose(core::pose::Pose & mono_pose) const {
	// 	core::conformation::symmetry::make_symmetric_pose( mono_pose, *symm_data_full );
	// }
	// Size get_region(Size rsd) {
	// 	for(Size i = 1; i<= cuts_.size(); ++i) if( rsd <= cuts_[i] ) return i;
	// 	return cuts_.size()+1;
	// }
	void dump_pdb(string fname, bool fullsymm=false) {
		using namespace ObjexxFCL::format;

		core::pose::Pose out_pose;
		if( fullsymm ) {
		// 	// TR << "init pose nres " << pose.n_residue() << std::endl;
		// 	out_pose = make_mono_pose(pose);
		// 	// TR << "mono pose nres " << out_pose.n_residue() << std::endl;
		// 	make_full_pose(out_pose);
		// 	// TR << "full pose nres " << out_pose.n_residue() << std::endl;
		//
			out_pose = pose;
		} else  {
			out_pose = pose;
		}

		utility::io::ozstream out(fname);

		out << "REMARK   1" << std::endl;
		out << "REMARK   2 " << tag_ << std::endl;
		for(Size i = 1; i <= attached_scattach_res_.size(); ++i) {
			string s = " ABABABABABABABABABABABABABABABABAB";
			out << "REMARK   3 " << attached_scattach_res_[i] << " " << s[which_subsub(attached_scattach_res_[i])] << std::endl;
		}

		if(basic::options::option[basic::options::OptionKeys::smhybrid::add_metal_at_0]()) {
			out << "HETATM 9999 ZN    ZN A 999       0.000   0.000   0.000  1.00  0.00              " << std::endl;
		}

		vector1<xyzVector<Real> > added;
		{
			utility::options::StringVectorOption & add_atom_at_cen = option[ basic::options::OptionKeys::smhybrid::add_atom_at_cen ];
			core::kinematics::FoldTree const & ft( out_pose.fold_tree());
			for( Size i = 1; i <= ft.num_jump(); ++i ) {
				for( Size j = 1; j <= attach_rsd_.size(); ++j ) {
					if( add_atom_at_cen[j] == "_" ) continue;
					if( (ft.downstream_jump_residue(i)-1)%nres+1 == attach_rsd_[j] && ft.downstream_jump_residue(i) <= (int)(nsub*nres) ) {
						Vec o = out_pose.xyz(core::id::AtomID(1,ft.upstream_jump_residue(i)));
						added.push_back(o);
						// TR << "output metal at origin of " << ft.upstream_jump_residue(i) << " " << add_atom_at_cen[j] << " " << o << std::endl;
					}
				}
			}
		}

		numeric::xyzVector<Real> shift;
		// center output on added atom, if only one
		Real mxdis = 0.0;
		for(Size i = 1; i <= added.size(); ++i) {
			for(Size j = i; j <= added.size(); ++j) {
				mxdis = numeric::max(added[i].distance_squared(added[j]),mxdis);
			}
		}
		// if(added.size() > 0) {
		if(mxdis < 0.01 && added.size() > 0) {
			shift = added[1];
		} else {
			Real N = 0;
			for(Size i = 1; i <= pose.n_residue(); ++i) {
				for(Size j = 1; j <= pose.residue_type(i).natoms(); ++j) {
					N += 1.0;
					shift += out_pose.xyz(AtomID(j,i));
				}
			}
			shift /= N;
		}

		for(Size i = 1; i <= pose.n_residue(); ++i) {
			for(Size j = 1; j <= pose.residue_type(i).natoms(); ++j) {
				out_pose.set_xyz( AtomID(j,i), out_pose.xyz(AtomID(j,i)) - shift );
			}
		}

		out_pose.dump_pdb(out);

		{
			utility::options::StringVectorOption & add_atom_at_cen = option[ basic::options::OptionKeys::smhybrid::add_atom_at_cen ];
			core::kinematics::FoldTree const & ft( out_pose.fold_tree());
			for( Size i = 1; i <= ft.num_jump(); ++i ) {
				for( Size j = 1; j <= attach_rsd_.size(); ++j ) {
					if( add_atom_at_cen[j] == "_" ) continue;
					if( (ft.downstream_jump_residue(i)-1)%nres+1 == attach_rsd_[j] && ft.downstream_jump_residue(i) <= (int)(nsub*nres) ) {
						Vec o = out_pose.xyz(core::id::AtomID(1,ft.upstream_jump_residue(i)));
						Size subunit = (ft.downstream_jump_residue(i)-1) / nres;
						char chain = string("ABCDEFGHIJKLMNOPQRSTUVWXYZ")[(subunit-1)%26];
						string l = "HETATM 9999 "+add_atom_at_cen[j]+"    "+add_atom_at_cen[j]+" "+chain+" 999     "+F(7,3,o.x())+" "+F(7,3,o.y())+" "+F(7,3,o.z())+"  1.00  0.00";
						out << l << std::endl;
						// TR << "output metal at origin of " << ft.upstream_jump_residue(i) << " " << add_atom_at_cen[j] << " " << o << std::endl;
					}
				}
			}
		}

		if(basic::options::option[basic::options::OptionKeys::smhybrid::add_cavities]()) {
			core::scoring::packstat::output_packstat_pdb(out_pose,out);
		}

		// std::map< string, core::conformation::symmetry::VirtualCoordinate > const & m( symm_data->get_virtual_coordinates() );
		// for( std::map< string, core::conformation::symmetry::VirtualCoordinate >::const_iterator i = m.begin(); i != m.end(); ++i ) {
		// 	TR << i->first << " " << i->second.get_x() << " " << i->second.get_y() << " " << i->second.get_origin() << std::endl;
		// }

	}

	bool is_floating_sc(Size ir) {
		return std::find(floating_scs_.begin(),floating_scs_.end(),ir)!=floating_scs_.end();
	}

	void align_orig_pose( Size beg, Size end ) {
		using namespace core;
		id::AtomID_Map< id::AtomID > atom_map;
		id::initialize( atom_map, orig_pose, id::BOGUS_ATOM_ID ); // maps every atomid to bogus atom

		for ( Size i=beg; i<=end; ++i ) {
		  id::AtomID const id1( orig_pose.residue(i).atom_index("CA"), i );
		  id::AtomID const id2( pose     .residue(i).atom_index("CA"), i );
		  atom_map[ id1 ] = id2;
		}
		scoring::superimpose_pose( orig_pose, pose, atom_map );

	}

	bool switch_to_fa() {
		TR << "switch to fa" << std::endl;
		if(pose.is_fullatom()) {
			TR << "WARNING: switch to FA when already in fa!" << std::endl;
		} else {

			if(debug) pose.dump_pdb("before_switch_to_fa.pdb");

			if( basic::options::option[basic::options::OptionKeys::smhybrid::centroid_all_val]() ) {
				for(Size i = 1; i <= nres; ++i) {
					res_types_.push_back(&(pose.residue(i).type()));
					if(find(floating_scs_.begin(),floating_scs_.end(),i)==floating_scs_.end()) {
						core::chemical::replace_pose_residue_copying_existing_coordinates( pose, i, *res_types_[i] );
					}
				}
			}

			core::pose::Pose const cen( pose );
			// TR << "switch to fa:" << std::endl;
			core::chemical::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
			// check_scattach_res();
			// dump_pdb("switch_to_fa_1.pdb");
			//get_closest_res_for_sc_attach

			for(Size i = 1; i <= subsub_fixed_begin.size(); ++i) {
				if(subsub_fixed_begin[i]==0 || subsub_fixed_end[i]==0) continue;
				TR << "align orig_pose " << subsub_fixed_begin[i] << " " << subsub_fixed_end[i] << std::endl;
				align_orig_pose( subsub_fixed_begin[i], subsub_fixed_end[i] );
				TR << "set coords " << subsub_fixed_begin[i] << " " << subsub_fixed_end[i] << std::endl;
				for(Size rsd = subsub_fixed_begin[i]; rsd <= subsub_fixed_end[i]; ++rsd) {
					for(Size j = 1; j <= pose.residue(rsd).natoms(); ++j) {
						pose.set_xyz(core::id::AtomID(j,rsd),orig_pose.xyz(core::id::AtomID(j,rsd)));
					}
					// TR << "adding coordinate_constraint " << rsd << std::endl;
					// addcc(pose,AtomID(2,rsd),attach_rsd_[i],1.0);
				}
			}
			// dump_pdb("switch_to_fa_2.pdb");

			//TODO is the following necessary? need to switch to iterate over attach_rsd_?
			core::id::AtomID_Map<bool> missing;
			core::id::initialize(missing,pose,false);
			bool anymissing = false;
			for(Size ath = 1; ath <= attach_rsd_.size(); ++ath) {
				for(Size attach = attach_rsd_[ath]; attach < nsub*nres; attach += nres) {
					for(Size i = 1; i <= pose.residue(attach).natoms(); ++i) {
						string aname = pose.residue(attach).atom_name(i);
						if( cen.residue(attach).type().has(aname) ) {
							// TR << "switch_to_fa: setting coords for res " << attach << " atom " << aname << std::endl;
							pose.set_xyz(core::id::AtomID(i,attach),cen.residue(attach).xyz(aname));
						} else {
							// TR << "switch_to_fa: res " << attach << " missing " << aname << std::endl;
							missing[core::id::AtomID(i,attach)] = true;
							anymissing = true;
						}
					}
				}
			}
			if(anymissing) pose.conformation().fill_missing_atoms( missing );
			check_scattach_res();
		}
		// add_apcs_floating_sc_csts();
		check_scattach_res();
		add_apcs_local(true);
		check_scattach_res();
		if(debug) pose.dump_pdb("after_switch_to_fa.pdb");
		if(!replace_scattach_res()) {
			TR << "switch_to_fa: replace_scattach_res failed!" << std::endl;
			return false;
		}

		for(Size i = 1; i <= rep_edge_res_.size(); ++i) {
			core::core::pose::add_variant_type_to_pose_residue(pose,"SHOVE_STRAND",rep_edge_res_[i]);
			TR << "switch_to_fa setting edge strand repulsive: " << rep_edge_res_[i] << " " << pose.residue(rep_edge_res_[i]).atom_name(1) << std::endl;

		}

		check_scattach_res();
		// for(Size i = 1; i <= pose.residue(1).natoms(); ++i) {
		// 	TR << "BPY atom " << i << " " << pose.residue(1).atom_name(i) << std::endl;
		// }
		// std::exit(-1);
		// dump_pdb("switch_to_fa_3.pdb");

		// for( Size scid = 1; scid <= floating_scs_.size(); ++scid ) {
		// 	Vec f = pose.residue(floating_scs_[scid]).xyz("CA");
		// 	Real best = 9e9;
		// 	Size bsti = 0;
		// 	for( Size i = 1; i <= nres; ++i ) {
		// 		if(i == floating_scs_[scid]) continue;
		// 		Real const d = pose.residue(i).xyz("CA").distance(f);
		// 		if( d < best ) {
		// 			best = d;
		// 			bsti = i;cen_fold
		// 		}
		// 	}
		// 	floating_scs_attach_.push_back(bsti);
		//
		// }
		return true;
	}
	bool check() {
		for(Size i = 1; i <=attached_scattach_res_.size(); ++i) {
			if(pose.residue(attached_scattach_res_[i]).name3() != attached_scattach_res_name_[i]) {
				if(debug) std::cerr << "FA FAIL: " << pose.residue(attached_scattach_res_[i]).name3() << " " << attached_scattach_res_name_[i] << std::endl;
				return false;
			}
		}
		return true;
	}
	void switch_to_cen() {
		if(!pose.is_fullatom()) {
			TR << "WARNING: tried to switch to CEN when already in CEN! doing nothing. " << std::endl;
			return;
		}
		core::pose::Pose const fa( pose );
		// TR << "switch_to_cen switch type set" << std::endl;
		core::chemical::switch_to_residue_type_set( pose, core::chemical::CENTROID );
		// pose.dump_pdb("switch_to_cen0.pdb");
		//TODO is the following necessary? need to switch to iterate over attach_rsd_?
		// TR << "switch_to_cen copy coords" << std::endl;
		if( basic::options::option[basic::options::OptionKeys::smhybrid::centroid_all_val]() ) {
			res_types_.clear();
			for(Size i = 1; i <= nres; ++i) {
				res_types_.push_back(&(pose.residue(i).type()));
				if(find(floating_scs_.begin(),floating_scs_.end(),i)==floating_scs_.end()) continue;
				bool skip = false;
				for(Size j = 1; j <= subsub_starts_.size(); ++j) if(subsub_starts_[j]==i && subsub_ends_[j]==i) skip=true;
				if(skip) continue;
				if(i==1 && option[basic::options::OptionKeys::smhybrid::virtual_nterm]() ) continue;
				core::chemical::replace_pose_residue_copying_existing_coordinates( pose, i, cen_residue_set_->name_map("PHE") );
			}
		}
		// pose.dump_pdb("switch_to_cen1.pdb");
		// for(Size i = 1; i <= subsub_fixed_begin.size(); ++i) {
		// 	if(subsub_fixed_begin[i]==0 || subsub_fixed_end[i]==0) continue;
		// 	for(Size rsd = subsub_fixed_begin[i]; rsd <= subsub_fixed_end[i]; ++rsd) {
		// 		// addcc(pose,AtomID(2,rsd),attach_rsd_[i],1.0);
		// 	}
		// }


		core::id::AtomID_Map<bool> missing;
		core::id::initialize(missing,pose,false);
		bool anymissing = false;
		for(Size ath = 1; ath <= attach_rsd_.size(); ++ath) {
			for(Size attach = attach_rsd_[ath]; attach <= nsub*nres; attach += nres) {
				for(Size i = 1; i <= pose.residue(attach).natoms(); ++i) {
					string aname = pose.residue(attach).atom_name(i);
					// if(aname=="CEN")
					if( fa.residue(attach).type().has(aname) ) {
						// TR << "switch_to_cen: setting coords for res " << attach << " atom " << aname << std::endl;
						pose.set_xyz(core::id::AtomID(i,attach),fa.residue(attach).xyz(aname));
					} else {
						TR << "WARNING: switch_to_cen: res " << attach << " missing " << aname << std::endl;
						missing[core::id::AtomID(i,attach)] = true;
						anymissing = true;
					}
				}
			}
		}

		for(Size i = 1; i <= rep_edge_res_.size(); ++i) {
			core::core::pose::add_variant_type_to_pose_residue(pose,"SHOVE_STRAND",rep_edge_res_[i]);
			TR << "switch_to_cen setting edge strand repulsive: " << rep_edge_res_[i] << " " << pose.residue(rep_edge_res_[i]).atom_name(1) << std::endl;
		}
		// pose.dump_pdb("switch_to_cen2.pdb");
		// TR << "switch_to_cen fill missing" << std::endl;
		if(anymissing) pose.conformation().fill_missing_atoms( missing );
		// TR << "switch_to_cen DONE" << std::endl;
		// pose.dump_pdb("switch_to_cen3.pdb");
		add_apcs_local(false);

	}


	Size get_closest_res_for_sc_attach(Size float_rsd) {
		if(attached_scattach_res_.size()) {
			TR << "WARNING: get_closest_res_for_sc_attach called after floating res attached" << std::endl;
			return 0;
		}
		using namespace core::scoring::constraints;
		ConstraintCOPs cs = pose.constraint_set()->get_all_constraints();

		// TR << pose.constraint_set()->get_all_constraints()[1] << std::endl;
		// TR << pose.constraint_set()->get_all_constraints()[1] << std::endl;
		// TR << pose.constraint_set()->get_all_constraints()[1] << std::endl;
		// std::exit(-1);

		// core::scoring::EnergyMap weights,emap;
		// weights[core::scoring::atom_pair_constraint] = 1;
		// core::scoring::constraints::ResiduePairXYZ xyzfunc(pose.residue(float_rsd),pose.residue(float_rsd));

		core::scoring::ScoreFunctionOP sf = new core::scoring::ScoreFunction;
		sf->set_weight(core::scoring::atom_pair_constraint,1.0);
		sf = new core::scoring::symmetry::SymmetricScoreFunction(*sf);
		(*sf)(pose);

		for(Size i = 1; i <= cs.size(); i++) {
			if( cs[i]->type() != "AmbiguousConstraint" ) continue;
			AmbiguousConstraint const * ac = (AmbiguousConstraint const *)(cs[i]());
			if(ac->active_constraint()->type() != "MultiConstraint") continue;
			MultiConstraint     const * mc = (MultiConstraint     const *)(ac->active_constraint()());
			AtomPairConstraint  const * pc = (AtomPairConstraint  const *)(mc->member_constraints()[1]());
			AtomID paid = pc->atom(1);
			AtomID faid = pc->atom(2);
			TR << "PAID " << paid.rsd() << " FAID " << faid.rsd() << std::endl;
			if(faid.rsd()==float_rsd) {
				return (paid.rsd()-1)%nres+1;
			}
		}
		utility_exit_with_message("couldn't find float sc attach point");
		return 0;


		// core::conformation::symmetry::SymmetricConformation const & SymmConf (
		// 	dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		// core::conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
		//
		// Vec floater( pose.residue(float_rsd).xyz("CA"));
		//
		// // vector1<bool> nearcut(pose.n_residue(),false);
		// // for(Size i = 1; i <= 4; ++i) nearcut[i] = true;
		// // for(Size i = 1; i <= (Size)pose.fold_tree().num_cutpoint(); ++i) {
		// // 	Size cp( pose.fold_tree().cutpoint(i) );
		// // 	if( !symm_info->bb_is_independent(cp) ) continue;
		// // 	for(Size j = numeric::max(((Size)1),cp-4); j <= numeric::min(cp+4,pose.n_residue()); ++j) nearcut[j] = true;
		// // }
		//
		// Size best = 9e9;
		// Size sc_attch_rsd = 0;
		// for( Size i = 1; i <= nres; ++i ) {
		// 	if ( ! pose.residue(i).is_protein() ) continue;
		// 	if( i == float_rsd ) continue;
		// 	bool isattach = false;
		// 	for(Size j = 1; j <= attach_rsd_.size(); ++j) if( attach_rsd_[j]==i && attach_atom_[j]!="N" ) isattach = true;
		// 	if(isattach) continue;
		// 	// if( nearcut[i] ) continue;
		//
		// 	Real dis = (floater.distance( pose.residue(i).xyz(2) )); dis = abs(dis);
		// 	if( dis < best ) {
		// 		// TR << "floating_sc " << i << " " << dis << std::endl;
		// 		sc_attch_rsd = i;
		// 		best = dis;
		// 	}
		// }
		// return sc_attch_rsd;
	}

	void update_designable_packable() {
		is_repackable_.resize(nres);
		is_designable_.resize(nres);
		for(Size i = 1; i <= nres; ++i) {
			is_designable_[i] = true;
			is_repackable_[i] = true;
		}
		for(Size i = 1; i <= nres; ++i) {
			if(ss(i)=='F') {
				is_designable_[i] = false;
				is_repackable_[i] = true;
			}
			if( find(linker_res_.begin(),linker_res_.end(),i) != linker_res_.end() ) {
				is_designable_[i] = basic::options::option[basic::options::OptionKeys::smhybrid::design_linker]();
			}
		}
		for(Size i = 1; i <= design_res_user.size(); ++i) {
			is_designable_[design_res_user[i]] = true;
			is_repackable_[design_res_user[i]] = true;
		}
		for(Size i = 1; i <= fixed_res_user.size(); ++i) {
			is_designable_[fixed_res_user[i]] = false;
			is_repackable_[fixed_res_user[i]] = false;
		}
		for(Size j = 1; j <= floating_scs_.size(); ++j) {
			Size r = get_closest_res_for_sc_attach(floating_scs_[j]);
			if(!r) continue;
			if( pose.is_fullatom() ) {
				core::chemical::replace_pose_residue_copying_existing_coordinates( pose, r, fa_residue_set_->name_map("ALA") );
				core::core::pose::add_variant_type_to_pose_residue(pose,"VIRTUAL_SC",r);
				TR << "setting res " << r << " to alanine" << std::endl;
			}
			is_designable_[r] = false;
			is_repackable_[r] = false;
		}
		if( basic::options::option[basic::options::OptionKeys::smhybrid::restrict_design_to_interface]() ) {
			for(Size i = 1; i <= nres; ++i) {
				// if(!is_designable_[i]) continue;
				Size is = which_sub(i);
				bool iface = false;
				for(Size j = 1; j <= nres*nsub; ++j) {
					if(which_sub(j)==is) continue;
					if( pose.xyz(AtomID(5,i)).distance_squared(pose.xyz(AtomID(2,j))) < 144.0 ) {
						iface = true;
						break;
					}
				}
				if(!iface) is_designable_[i] = false;
		 		else is_designable_[i] = true;
			}
		}
		if( basic::options::option[basic::options::OptionKeys::smhybrid::restrict_design_to_subsub_interface]() ) {
			for(Size i = 1; i <= nres; ++i) {
				if(!is_designable_[i]) continue;
				Size is = which_subsub(i);
				bool iface = false;
				for(Size j = 1; j <= nres*nsub; ++j) {
					if(which_subsub(j)==is) continue;
					if( pose.xyz(AtomID(5,i)).distance_squared(pose.xyz(AtomID(2,j))) < 144.0 ) {
						iface = true;
						break;
					}
				}
				if(!iface) is_designable_[i] = false;
			}
		}
		for(Size j = 1; j <= attach_rsd_.size(); ++j) {
			if(	attach_atom_[j]!="N") {
				is_designable_[attach_rsd_[j]] = false;
				is_repackable_[attach_rsd_[j]] = false;
			}
		}
		if( !basic::options::option[basic::options::OptionKeys::smhybrid::design]() ) {
			for(Size i = 1; i <= nres; ++i) is_designable_[i] = false;
		}
		for(Size i = 1; i <= attached_scattach_res_.size(); ++i) {
			is_designable_[attached_scattach_res_[i]] = false;
			is_repackable_[attached_scattach_res_[i]] = true;
		}
	}
	bool is_rsd_designable(Size i) const { return is_designable_[i]; }
	bool is_rsd_repackable(Size i) const { return is_repackable_[i]; }

	void add_apcs_local(bool fullatom = false) {
		// std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		// return;
		vector1<Size> all_subsub;
		// for(Size ss = 1; ss <= subsub_starts_.size(); ++ss) {
			// Size start=subsub_starts_[ss], end=subsub_ends_[ss];
		for(Size ss = 1; ss <= subsub_fixed_begin.size(); ++ss) {
			Size start=subsub_fixed_begin[ss], end=subsub_fixed_end[ss];
			vector1<Size> subs;
			if(cst_sub_files_.size() >= ss) subs = cst_sub_files_[ss];
			if(subs.size()==0) subs.push_back(1);

			vector1<Size> res;
			for(vector1<Size>::iterator i = subs.begin(); i != subs.end(); ++i) {
				for(Size j = start; j <= end; ++j) {
					res.push_back(j+(*i-1)*nres);
					all_subsub.push_back(j+(*i-1)*nres);
				}
			}

			for(Size i = start; i <= end; ++i) {
				if(std::find(frag_res_user.begin(),frag_res_user.end(),i) == frag_res_user.end()) { // not frag res
					if(!fullatom) continue;
					else TR << "FA adding local apcs to resi " << i << std::endl;
				} else { // is frag res
					if(fullatom) continue;
					else TR << "CEN adding local apcs to resi " << i << std::endl;
				}
				for(vector1<Size>::iterator j = res.begin(); j != res.end(); ++j) {
					if( i+5 > *j ) continue;
					Real d = pose.xyz(AtomID(2,i)).distance(pose.xyz(AtomID(2,*j)));
					if( d < 8.0 ) {
						// TR << "distance resi " << i << " and name CA, resi " << *j << " and name CA" << std::endl;
						Real wt = 2;
						// if(!fullatom) wt = 0.5;
						add_apc(pose,AtomID(2,i),AtomID(2,*j),d,wt,core::scoring::coordinate_constraint);
					}
				}
			}
		}
		if( basic::options::option[basic::options::OptionKeys::smhybrid::inter_subsub_cst]() ) {
			for(Size i = 1; i <= all_subsub.size(); ++i) {
				for(Size j = i+1; j <= all_subsub.size(); ++j) {
					Real d = pose.xyz(AtomID(2,i)).distance(pose.xyz(AtomID(2,j)));
					if( d < 10.0 ) {
						// TR << "distance resi " << i << " and name CA, resi " << *j << " and name CA" << std::endl;
						Real wt = 	4;
						if(fullatom) wt = 0.4;
						add_apc(pose,AtomID(5,i),AtomID(5,j),d,wt,core::scoring::coordinate_constraint);
					}
				}
			}
		}

		// trying to make sheet....
// 		using namespace core::scoring::constraints;
// 		ConstraintCOPs ac;
// 		FuncOP func = new HarmonicFunc(4,0.5);
// 		std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
// 		int a[2][2],b[2][2];
// 		a[0][0] =  30; a[0][1] =  34;
// 		a[1][0] =  81; a[1][1] =  85;
// 		b[0][0] = 135+nres; b[0][1] = 139+nres;
// 		b[1][0] = 172+nres; b[1][1] = 176+nres;
// 		for(int i = 0; i <= 1; ++i) {
// 			for(int j = 0; j <= 1; ++j) {
// 				TR << "SS pair " << a[i][0] << "-" << b[j][0] << " " << a[i][1] << "-" << b[j][1] << " " << std::endl;
// 				TR << "SS pair " << a[i][0] << "-" << b[j][1] << " " << a[i][1] << "-" << b[j][0] << " " << std::endl;
// 				add_apc(pose,AtomID(1,a[i][0]),AtomID(4,b[j][0]),3,1,core::scoring::coordinate_constraint);
// 				add_apc(pose,AtomID(1,a[i][1]),AtomID(4,b[j][1]),3,1,core::scoring::coordinate_constraint);

// 				{
// 				ConstraintCOPs mc;
// 				mc.push_back(new AtomPairConstraint( AtomID(2,a[i][0]), AtomID(2,b[j][0]), func ) );
// 				mc.push_back(new AtomPairConstraint( AtomID(2,a[i][1]), AtomID(2,b[j][1]), func ) );
// 				ac.push_back(new MultiConstraint(mc));
// 				}
// 				{
// 				ConstraintCOPs mc;
// 				mc.push_back(new AtomPairConstraint( AtomID(2,a[i][0]), AtomID(2,b[j][1]), func ) );
// 				mc.push_back(new AtomPairConstraint( AtomID(2,a[i][1]), AtomID(2,b[j][0]), func ) );
// 				ac.push_back(new MultiConstraint(mc));
// 				}

// 			}
// 		}
// 		if(ac.size()) pose.add_constraint(new AmbiguousConstraint(ac));
//


	}

	// void add_fsc_harmonic_cst() {
	// 	using namespace core::scoring::constraints;
	// 	for(Size i = 1; i <= floating_scs_.size(); ++i ) {
	// 		Size r = attached_scattach_res_[i], fsc = floating_scs_[i] + (floating_scs_sub_[i]-1)*nres;
	// 		for(Size j = 1; j <= pose.residue(fsc).nheavyatoms(); ++j) {
	// 			if( pose.residue(r).has(pose.residue(fsc).atom_name(j)) ) {
	// 				Real wt = 1.0;
	// 				if(j < 5) wt = 3.0;
	// 				pose.add_constraint( new AtomPairConstraint(
	// 					AtomID(pose.residue(r).atom_index(pose.residue(fsc).atom_name(j)),r),
	// 					AtomID(j,fsc), new HarmonicFunc(0.0,wt) ) );
	// 			}
	// 		}
	// 	}
	//
	// }

	bool check_scattach_res() {
		vector1<Vec> other_tetra;
		for(Size i = 1; i <= floating_scs_.size(); ++i) {
			if(pose.residue(floating_scs_[i]).name3()=="HSC") {
				Vec cen = pose.xyz(AtomID(1,up_jump_tree(pose,floating_scs_[i])));
				Vec ne2 = pose.residue(floating_scs_[i]).xyz("NE2");
				Vec nd1 = pose.residue(floating_scs_[i]).xyz("ND1");
				Vec ce1 = pose.residue(floating_scs_[i]).xyz("CE1");
				// Vec cd2 = pose.residue(floating_scs_[i]).xyz("CD2");
				// Vec cg  = pose.residue(floating_scs_[i]).xyz("CG" );
				Real dd = cen.distance(nd1), de = cen.distance(ne2);
				bool fail = false;
				if(dd<de) {
					if( fabs(dd-2.0274) > 0.02 ) fail = true;
					// TR << "ANGLE2 " << numeric::angle_radians(cen,nd1,ce1) << std::endl;
					if( fabs(numeric::angle_radians(cen,nd1,ce1)-2.158972) > 0.01 ) {
						fail = true;
						TR << "check_scattach_res FAIL angle cen nd1 ce1 " << i << " " << fabs(numeric::angle_radians(cen,nd1,ce1)-2.21751) << std::endl;
						TR << cen << std::endl;
						TR << nd1 << std::endl;
						TR << ce1 << std::endl;
					}
					for(Size j = 1; j <= other_tetra.size(); j++) {
						// TR << "tetra: " << numeric::angle_radians(nd1,cen,other_tetra[j]) << std::endl;
						if(fabs(numeric::angle_radians(nd1,cen,other_tetra[j])-1.911136) > 0.02) {
							fail = true;
							TR << "check_scattach_res FAIL tetra " << i << " " << j << " " << fabs(numeric::angle_radians(nd1,cen,other_tetra[j])-1.911136) << std::endl;
							TR << nd1 << std::endl;
							TR << cen << std::endl;
							TR << other_tetra[j] << std::endl;
						}
					}
					other_tetra.push_back(nd1);
				} else {
					if( fabs(de-2.0274) > 0.02 ) fail = true;
					// TR << "ANGLE2 " << numeric::angle_radians(cen,ne2,ce1) << std::endl;
					if( fabs(numeric::angle_radians(cen,ne2,ce1)-2.21751) > 0.01 ) {
						fail = true;
						TR << "check_scattach_res FAIL angle cen ne2 ce1 " << i << " " << fabs(numeric::angle_radians(cen,ne2,ce1)-2.21751) << std::endl;
						TR << cen << std::endl;
						TR << ne2 << std::endl;
						TR << ce1 << std::endl;
					}
					for(Size j = 1; j <= other_tetra.size(); j++) {
						// TR << "tetra: " << numeric::angle_radians(ne2,cen,other_tetra[j]) << std::endl;
						if(fabs(numeric::angle_radians(ne2,cen,other_tetra[j])-1.911136) > 0.02) {
							fail = true;
							TR << "check_scattach_res FAIL tetra " << i << " " << j << " " << fabs(numeric::angle_radians(ne2,cen,other_tetra[j])-1.911136) << std::endl;
							TR << ne2 << std::endl;
							TR << cen << std::endl;
							TR << other_tetra[j] << std::endl;
						}
					}
					other_tetra.push_back(ne2);
				}
				if( option[basic::options::OptionKeys::smhybrid::stop_on_bad_geom]() && fail) {
					dump_pdb("ERROR.pdb");
					utility_exit_with_message("check_scattach_res failed!");
				}
			} else {
				if( option[basic::options::OptionKeys::smhybrid::stop_on_bad_geom]()) {
					utility_exit_with_message("don't know how to check geom of res "+pose.residue(floating_scs_[i]).name());
				}
			}
		}
		return true;
	}

	bool replace_scattach_res() {
		TR << "REPLACE SCATTACH RES!!!!!!!!!!!!!!!!!" << std::endl;
		TR << pose.fold_tree() << std::endl;
		TR << pose.residue(1).name() << std::endl;
		vector1<Size> scattach_res;
		for(Size i = 1; i <= floating_scs_.size(); ++i ) {
			Size r = get_closest_res_for_sc_attach(floating_scs_[i]);
			scattach_res.push_back(r);
			TR << "scattach res " << i << " " << r << std::endl;
		}
		FoldTree ft = pose.fold_tree();
		// std::exit(-1);
		// for(Size i = 1; i <= floating_scs_.size(); ++i ) {
		// 	for(Size j = 1; j <= (Size)ft.num_jump(); ++j ) {
		// 		if( ((Size)ft.downstream_jump_residue(j)-1)%nres+1 == floating_scs_[i] ) {
		// 			ft.set_jump_atoms(j,"","");
		// 		}
		// 	}
		// }
		// pose.fold_tree(ft);


		TR << "NCST: " << pose.constraint_set()->get_all_constraints().size() << std::endl;
		using namespace core::scoring::constraints;
		ConstraintCOPs cs = pose.constraint_set()->get_all_constraints();
		ConstraintCOPs rem;
		for(Size i = 1; i <= cs.size(); i++) {
			if( cs[i]->type() != "AmbiguousConstraint" ) continue;
			rem.push_back(cs[i]);
		}
		pose.remove_constraints(rem);

		TR << "NCST: " << pose.constraint_set()->get_all_constraints().size() << std::endl;
		TR << "=========================================================" << std::endl;
		TR << pose.fold_tree() << std::endl;
		check_scattach_res();

		attached_scattach_res_.clear();
		attached_scattach_res_name_.clear();
		for(Size i = 1; i <= floating_scs_.size(); ++i ) {
			Size r = scattach_res[i], fsc = floating_scs_[i] + (floating_scs_sub_[i]-1)*nres;
			r = (r-1) % (nres) + 1;
			if( pose.residue(fsc).name3() == "HSC" ) {
				core::chemical::replace_pose_residue_copying_existing_coordinates( pose, r, fa_residue_set_->name_map("HIS_D") );
				// for(Size j = 1; j <= pose.fold_tree().num_jump(); ++j) {
				// 	if(pose.fold_tree().downstream_jump_residue(j) == (int)fsc) {
				// 		TR << "adding dist cst from HIS N to virt ZN" << std::endl;
				// 		Size upres = pose.fold_tree().upstream_jump_residue(j);
				// 		pose.add_constraint( new AtomPairConstraint(
				// 			AtomID(pose.residue(r).atom_index(pose.fold_tree().downstream_atom(j)),r),
				// 			AtomID(1,upres), new AbsFunc(2.0,0.5) ) );
				// 		break;
				// 	}
				// }
			} else if ( pose.residue(fsc).name3() == "BPY" ) {
				core::chemical::replace_pose_residue_copying_existing_coordinates( pose, r, fa_residue_set_->name_map("BPY") );
			}
			check_scattach_res();
			if(debug) dump_pdb("before_virt.pdb");
			TR << "make floating_sc " << floating_scs_[i] << " virtual" << std::endl;
			vector1<core::conformation::ResidueOP> tmp;
			for(Size l = 0; l < nsub; ++l) {
				tmp.push_back(new core::conformation::Residue(pose.residue(floating_scs_[i]+l*nres)));
			}
			core::core::pose::add_variant_type_to_pose_residue(pose,"VIRTUAL",floating_scs_[i]);
			for(Size l = 0; l < nsub; ++l) {
				for(Size k = 1; k <= tmp[l+1]->natoms(); k++) {
					if(pose.residue(floating_scs_[i]+l*nres).has(tmp[l+1]->atom_name(k))) {
						pose.set_xyz(AtomID(pose.residue(floating_scs_[i]+l*nres).atom_index(tmp[l+1]->atom_name(k)),floating_scs_[i]+l*nres),tmp[l+1]->xyz(k));
					}
				}
			}
			if(debug) dump_pdb("after_virt.pdb");
			check_scattach_res();
			for(Size j = 1; j <= pose.residue(fsc).nheavyatoms(); ++j) {
				if( pose.residue(r).has(pose.residue(fsc).atom_name(j)) ) {
					Real wt = 1.0;
					if(j < 5) wt = 3.0;
					pose.add_constraint( new AtomPairConstraint(
						AtomID(pose.residue(r).atom_index(pose.residue(fsc).atom_name(j)),r),
						AtomID(j,fsc), new AbsFunc(0.0,wt) ) );
				} else {
					TR << "floating sc FA cst missing atom " << pose.residue(fsc).atom_name(j) << std::endl;
				}
			}
			attached_scattach_res_.push_back(r);
			attached_scattach_res_name_.push_back(pose.residue(r).name3());
			check_scattach_res();
		}

		for(Size i = 1; i <= attached_scattach_res_.size(); ++i) {
			for(Size j = i+1; j <= attached_scattach_res_.size(); ++j) {
				if(attached_scattach_res_[i]==attached_scattach_res_[j]) return false;
			}
		}
		// ft = pose.fold_tree();
		// for(Size j = 1; j <= (Size)ft.num_jump(); ++j ) {
		// 	for(Size i = 1; i <= scattach_res.size(); ++i ) {
		// 		if( ((Size)ft.downstream_jump_residue(j)-1)%nres+1 == floating_scs_[i] ) ft.set_jump_atoms(j,"ORIG","ORIG");
		// 	}
		// }
		// pose.fold_tree(ft);

		return true;
	}

	Real rms_to_orig_subsubs() {
		Real rms = 0.0;
		Size n = 0;
		for(Size i = 1; i <= subsub_starts_.size(); ++i) {
			// TR << "align orig_pose " << subsub_starts_[i] << " " << subsub_ends_[i] << std::endl;
			align_orig_pose( subsub_starts_[i], subsub_ends_[i] );
			// TR << "set coords " << subsub_starts_[i] << " " << subsub_ends_[i] << std::endl;
			for(Size rsd = subsub_starts_[i]; rsd <= subsub_ends_[i]; ++rsd) {
				for(Size j = 1; j <= pose.residue(rsd).nheavyatoms(); ++j) {
					if( j > orig_pose.residue(rsd).nheavyatoms() ) continue;
					if( pose.residue(rsd).atom_name(j) != orig_pose.residue(rsd).atom_name(j) ) continue;
					rms += pose.xyz(core::id::AtomID(j,rsd)).distance_squared( orig_pose.xyz(core::id::AtomID(j,rsd)));
					n += 1;
				}
			}
		}
		return std::sqrt(rms/n);
	}

	Real absrms_to_orig_subsubs() {
		Real rms = 0.0;
		Size n = 0;
		for(Size i = 1; i <= subsub_starts_.size(); ++i) {
			// TR << "set coords " << subsub_starts_[i] << " " << subsub_ends_[i] << std::endl;
			for(Size rsd = subsub_starts_[i]; rsd <= subsub_ends_[i]; ++rsd) {
				for(Size j = 1; j <= pose.residue(rsd).nheavyatoms(); ++j) {
					if( j > orig_pose.residue(rsd).nheavyatoms() ) continue;
					if( j > 4 ) continue;
					if( pose.residue(rsd).atom_name(j) != orig_pose.residue(rsd).atom_name(j) ) continue;
					rms += pose.xyz(core::id::AtomID(j,rsd)).distance_squared( orig_pose.xyz(core::id::AtomID(j,rsd)));
					n += 1;
				}
			}
		}
		return std::sqrt(rms/n);
	}


	void center() {
		;
	}
};
typedef owning_ptr<PoseWrap> PoseWrapOP;
typedef access_ptr<PoseWrap> PoseWrapAP;
typedef owning_ptr<const PoseWrap> PoseWrapCOP;
typedef access_ptr<const PoseWrap> PoseWrapCAP;


vector1<Size> read_res_list(string fn) {
	vector1<Size> l;
	if(fn=="_") return l;
	if(fn.size()==1 && fn[0]==(char)0) return l;
	izstream in(fn);
	if(!in.good()) {
		utility_exit_with_message("can't open res list file '"+fn+"'");
	}
	Size r;
	while( in >> r ) l.push_back(r);
	return l;
}

vector1<xyzVector<Size> > read_jumpcut_file(string fn) {
	vector1<numeric::xyzVector<Size> > l;
	if(fn=="") return l;
	if(fn=="_") return l;
	if(fn.size()==1 && fn[0]==(char)0) return l;
	izstream in(fn);
	if(!in.good()) {
		utility_exit_with_message("can't open file "+fn);
	}
	Size x,y,z;
	while(in >> x >> y >> z) {
		l.push_back(xyzVector<Size>(x,y,z));
	}
	return l;
}

string select_string(string in, Size & fp) {
	vector1<string> files(1,"");
	Size i=0;
	do {
		if(in[i]!='|') files.back() += in[i];
		else files.push_back("");
		++i;
	} while(i < in.size());
	if(files.size()==1) return files[1];
	if(fp==0) fp = std::ceil(uniform()*files.size());
	return files[fp];
}

Size to_integer(string s) {
	std::istringstream iss(s);
	Size i;
	iss >> i;
	return i;
}

PoseWrap posewrap_from_command_line(string symm_def_template = "", string symm_def_template_reduced = "") {
	using namespace utility::options;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	basic::options::FileVectorOption & pdbs = option[ in::file::s ];
	vector1<core::pose::Pose> poses;
	// TR << "found poses " << pdbs.size() << std::endl;
	vector1<Size> filepick;
	string tag;
	vector1<string> attach_rsd = option[ smhybrid::attach_rsd ]();
	vector1<string> attach_atom = option[ smhybrid::attach_atom ]();
	vector1<string> ss_pref = option[ smhybrid::add_ss_before ]();
	vector1<string> ss_suff = option[ smhybrid::add_ss_after ]();
	vector1<string> res_pref = option[ smhybrid::add_res_before ]();
	vector1<string> res_suff = option[ smhybrid::add_res_after ]();
	vector1<Size> attach_as_sc;
	vector1<vector1<Size> > design_res,fixed_res,frag_res,scattach_res,virtual_res,rep_edge_res;
	vector1<vector1<xyzVector<Size> > > jumpcut;
	vector1<vector1<Size> > cst_sub_files;
	vector1<bool> chainbreaks;

	if(option[smhybrid::attach_as_sc].user()) attach_as_sc = option[ smhybrid::attach_as_sc ]();
	// convert to vector1s
	vector1<Size> myattach;
	vector1<Size> mysc;
	vector1<string> myattach_atom,mypref,mysuff,myrespref,myressuff;

	for(Size i = 1; i <= pdbs.size(); ++i ) {
		// TR << "reading pose from " << pdbs[i] << std::endl;
		core::pose::Pose tmp;
		Size fp = 0;
		string posefile = select_string(pdbs[i],fp);
		if( attach_as_sc.size() < i || !attach_as_sc[i] ) {
			if(tag.size() > 1) tag += "-";
			tag += utility::file_basename(posefile.substr(0,posefile.size()-4));
		}
		core::io::pdb::pose_from_pdb(tmp,posefile);
		filepick.push_back(fp);
		// TR << pdbs[i] << " " << select_string(pdbs[i]) << std::endl;
		poses.push_back( tmp );
	}

	for(Size i = 1; i <= pdbs.size(); ++i ) {
		if( attach_rsd.size() > 1 && attach_rsd.size() != pdbs.size() ) utility_exit_with_message("attach_rsd wrong size!");
		if( attach_rsd.size() < i ) myattach.push_back(1);
		else myattach.push_back( to_integer(select_string(attach_rsd[i],filepick[i])));
		if( attach_atom.size() > 1 && attach_atom.size() != pdbs.size() ) utility_exit_with_message("attach_atom wrong size!");
		if( attach_atom.size() < i ) myattach_atom.push_back("N");
		else myattach_atom.push_back(select_string(attach_atom[i],filepick[i]));
		if( ss_pref.size() > 1 && ss_pref.size() != pdbs.size() ) utility_exit_with_message("ss_pref wrong size!");
		if( ss_pref.size() < i ) mypref.push_back("");
		else if(ss_pref[i][0]=='_') mypref.push_back("");
		else mypref.push_back( process_ss_str(select_string(ss_pref[i],filepick[i])) );
		if( ss_suff.size() > 1 && ss_suff.size() != pdbs.size() ) utility_exit_with_message("ss_suff wrong size!");
		if( ss_suff.size() < i ) mysuff.push_back("");
		else if(ss_suff[i][0]=='_') mysuff.push_back("");
		else mysuff.push_back( process_ss_str(select_string(ss_suff[i],filepick[i])) );
		if( res_pref.size() > 1 && res_pref.size() != pdbs.size() ) utility_exit_with_message("res_pref wrong size!");
		if( res_pref.size() < i ) myrespref.push_back("");
		else if(res_pref[i][0]=='_') myrespref.push_back("");
		else myrespref.push_back( select_string(res_pref[i],filepick[i]) );
		if( res_suff.size() > 1 && res_suff.size() != pdbs.size() ) utility_exit_with_message("res_suff wrong size!");
		if( res_suff.size() < i ) myressuff.push_back("");
		else if(res_suff[i][0]=='_') myressuff.push_back("");
		else myressuff.push_back( select_string(res_suff[i],filepick[i]) );
		if( attach_as_sc.size() > 1 && attach_as_sc.size() != pdbs.size() ) utility_exit_with_message("attach_as_sc wrong size!");
		if( attach_as_sc.size() < i ) mysc.push_back(0);
		else mysc.push_back(attach_as_sc[i]);
		if( option[smhybrid::design_res_files]().size() > 1 && option[smhybrid::design_res_files]().size() != pdbs.size() ) utility_exit_with_message("option[smhybrid::design_res_files]() wrong size!");
		if( option[smhybrid::design_res_files]().size() < i) design_res.push_back(vector1<Size>());
		else design_res.push_back(read_res_list(select_string(option[smhybrid::design_res_files]()[i],filepick[i])));

		if( option[smhybrid::rep_edge_files]().size() > 1 && option[smhybrid::rep_edge_files]().size() != pdbs.size() ) utility_exit_with_message("option[smhybrid::rep_edge_files]() wrong size!");
		if( option[smhybrid::rep_edge_files]().size() < i) rep_edge_res.push_back(vector1<Size>());
		else rep_edge_res.push_back(read_res_list(select_string(option[smhybrid::rep_edge_files]()[i],filepick[i])));

		if( option[smhybrid::fixed_res_files]().size() > 1 && option[smhybrid::fixed_res_files]().size() != pdbs.size() ) utility_exit_with_message("option[smhybrid::fixed_res_files]() wrong size!");
		if( option[smhybrid::fixed_res_files]().size() < i) fixed_res.push_back(vector1<Size>());
		else fixed_res.push_back(read_res_list(select_string(option[smhybrid::fixed_res_files]()[i],filepick[i])));
		if( option[smhybrid::frag_res_files]().size() > 1 && option[smhybrid::frag_res_files]().size() != pdbs.size() ) utility_exit_with_message("option[smhybrid::frag_res_files]() wrong size!");
		if( option[smhybrid::frag_res_files]().size() < i) frag_res.push_back(vector1<Size>());
		else frag_res.push_back(read_res_list(select_string(option[smhybrid::frag_res_files]()[i],filepick[i])));
		if( option[smhybrid::virtual_res_files]().size() > 1 && option[smhybrid::virtual_res_files]().size() != pdbs.size() ) utility_exit_with_message("option[smhybrid::virtual_res_files]() wrong size!");
		if( option[smhybrid::virtual_res_files]().size() < i) virtual_res.push_back(vector1<Size>());
		else virtual_res.push_back(read_res_list(select_string(option[smhybrid::virtual_res_files]()[i],filepick[i])));
		if( attach_rsd.size() > 1 && attach_rsd.size() != pdbs.size() ) utility_exit_with_message("attach_rsd wrong size!");
		if( option[smhybrid::scattach_res_files]().size() < i) scattach_res.push_back(vector1<Size>());
		else scattach_res.push_back(read_res_list(select_string(option[smhybrid::scattach_res_files]()[i],filepick[i])));
		if( option[smhybrid::jumpcut_files]().size() > 1 && option[smhybrid::jumpcut_files]().size() != pdbs.size() ) utility_exit_with_message("option[smhybrid::jumpcut_files]() wrong size!");
		if( option[smhybrid::jumpcut_files]().size() < i) jumpcut.push_back(vector1<xyzVector<Size> >());
		else jumpcut.push_back(read_jumpcut_file(select_string(option[smhybrid::jumpcut_files]()[i],filepick[i])));
		if( option[smhybrid::cst_sub_files]().size() > 1 && option[smhybrid::cst_sub_files]().size() != pdbs.size() ) utility_exit_with_message("option[smhybrid::cst_sub_files]() wrong size!");
		if( option[smhybrid::cst_sub_files]().size() < i) cst_sub_files.push_back(vector1<Size>());
		else cst_sub_files.push_back(read_res_list(select_string(option[smhybrid::cst_sub_files]()[i],filepick[i])));
		if( option[smhybrid::chainbreaks]().size() > 1 && option[smhybrid::chainbreaks]().size() != pdbs.size() ) utility_exit_with_message("option[smhybrid::chainbreaks]() wrong size!");
		if( option[smhybrid::chainbreaks]().size() < i) chainbreaks.push_back(false);
		else chainbreaks.push_back(/*select_string(*/option[smhybrid::chainbreaks]()[i]/*,filepick[i])*/);

	}
	if( symm_def_template.size()==0 ) {
		std::string fname = option[smhybrid::symm_def_template]();
		izstream in(fname.c_str());
		if(!in.good()) utility_exit_with_message("can't open sym template file: "+fname);
		char buf[999];
		while( in.getline(buf,999) ) {
			symm_def_template += buf;
			symm_def_template += "\n";
		}
		in.close();
	}
	if( symm_def_template_reduced.size()==0 ) {
		std::string fname = option[smhybrid::symm_def_template_reduced]();
		izstream in(fname.c_str());
		if(!in.good()) utility_exit_with_message("can't open sym template file: "+fname);
		char buf[999];
		while( in.getline(buf,999) ) {
			symm_def_template_reduced += buf;
			symm_def_template_reduced += "\n";
		}
		in.close();
	}
	return PoseWrap(tag,poses,myattach,myattach_atom,mypref,mysuff,myrespref,myressuff,mysc,symm_def_template,symm_def_template_reduced,design_res,rep_edge_res,fixed_res,frag_res,virtual_res,scattach_res,jumpcut,cst_sub_files,chainbreaks);
}


void read_fragdata( vector1< core::fragment::FragDataOP > & fds, std::istream & in, bool /*design = false*/ ) {
	using namespace core::fragment;
	Size n,count=0;
	while( in >> n ) {
	 	string pdb;
		char buf[999];
		FragDataOP fd = new FragData;
		for( Size i = 1; i <= n; ++i ) {
			utility::pointer::owning_ptr<SingleResidueFragData> srfd;
			srfd = new BBTorsionSRFD;
			in >> pdb;
			in.getline(buf,999);
			std::istringstream iss(buf);
			iss >> *srfd;
			fd->add_residue(srfd);
		}
		fd->set_valid(true);
		fds.push_back(fd);
		count++;
		// if( count >= nread ) break;
	}
}

std::map<string, vector1<core::fragment::FragDataOP> >
get_frags_map( bool design = false ) {
	using namespace core::fragment;
	TR << "reading frags" << std::endl;
	utility::io::izstream in;
	std::map<string,vector1<FragDataOP> > fds;
	core::io::database::open(in,"sampling/ss_fragfiles/EEE.fragfile"); read_fragdata(fds["EEE"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/EEH.fragfile"); read_fragdata(fds["EEH"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/EEL.fragfile"); read_fragdata(fds["EEL"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/EHH.fragfile"); read_fragdata(fds["EHH"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/ELE.fragfile"); read_fragdata(fds["ELE"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/ELH.fragfile"); read_fragdata(fds["ELH"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/ELL.fragfile"); read_fragdata(fds["ELL"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/HEE.fragfile"); read_fragdata(fds["HEE"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/HEH.fragfile"); read_fragdata(fds["HEH"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/HEL.fragfile"); read_fragdata(fds["HEL"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/HHE.fragfile"); read_fragdata(fds["HHE"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/HHH.fragfile"); read_fragdata(fds["HHH"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/HHL.fragfile"); read_fragdata(fds["HHL"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/HLE.fragfile"); read_fragdata(fds["HLE"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/HLH.fragfile"); read_fragdata(fds["HLH"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/HLL.fragfile"); read_fragdata(fds["HLL"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/LEE.fragfile"); read_fragdata(fds["LEE"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/LEH.fragfile"); read_fragdata(fds["LEH"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/LEL.fragfile"); read_fragdata(fds["LEL"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/LHH.fragfile"); read_fragdata(fds["LHH"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/LLE.fragfile"); read_fragdata(fds["LLE"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/LLH.fragfile"); read_fragdata(fds["LLH"],in,design); in.close();
	core::io::database::open(in,"sampling/ss_fragfiles/LLL.fragfile"); read_fragdata(fds["LLL"],in,design); in.close();
	return fds;
}


core::fragment::FragSetOP make_frag_set(PoseWrap const & pw, std::map<string, vector1<core::fragment::FragDataOP> > fds) {
	using namespace core::fragment;
	FragSetOP frags = new ConstantLengthFragSet();
	int const stop = pw.nres + 1 - 3;
	if((int)1 >= stop) return NULL;
	for( Size i = 1; i <= (Size)stop; ++i ) {
		string ss3 = pw.ss(i,3);
		bool mkframe = true;
		for(Size j = 0; j < ss3.size(); ++j) if(ss3[j]!='H'&&ss3[j]!='E'&&ss3[j]!='L'&&ss3[j]!='*') mkframe = false;
		for(Size j = 1; j <= pw.cuts_.size(); ++j ) {
			Size c = pw.cuts_[j];
			if( i==c ) mkframe=false; if( i==c-1 ) mkframe=false;
		}
		if( !mkframe ) continue;
		FrameOP frame = new Frame(i,3);
		vector1<char> ss0,ss1,ss2;
		if('*'==ss3[0]) { ss0.push_back('H'); ss0.push_back('E'); ss0.push_back('L'); } else ss0.push_back(ss3[0]);
		if('*'==ss3[1]) { ss1.push_back('H'); ss1.push_back('E'); ss1.push_back('L'); } else ss1.push_back(ss3[1]);
		if('*'==ss3[2]) { ss2.push_back('H'); ss2.push_back('E'); ss2.push_back('L'); } else ss2.push_back(ss3[2]);
		for( Size j = 1; j <= ss0.size(); ++j ) {
		for( Size k = 1; k <= ss1.size(); ++k ) {
		for( Size l = 1; l <= ss2.size(); ++l ) {
			string ss=""; ss+=ss0[j]; ss+=ss1[k]; ss+=ss2[l];
			// TR << "adding ss " << ss << " '" << ss0[j] << "' '" << ss1[k] << "' '" << ss2[l] << "'" << std::endl;
			vector1<FragDataOP>::iterator beg = fds[ss].begin();
			vector1<FragDataOP>::iterator end = fds[ss].end();
			for( vector1<FragDataOP>::iterator fi = beg; fi != end; ++fi ) {
				frame->add_fragment(*fi);
			}
		}}}
		frags->add(frame);
		TR << "make frag " << i << ": " << ss3 << std::endl;
	}
	if(frags->size() == 0) return NULL;
	return frags;
}

core::fragment::FragSetOP make_frag_set_9mers(Size nres) {
	using namespace core::fragment;

	vector1<core::fragment::FragDataOP> fds9;
	std::ifstream in("input/frags/loop_helix.9mers");
	read_fragdata(fds9,in,false);

	FragSetOP frags = new ConstantLengthFragSet();
	for( Size i = 1; i <= nres-8; ++i ) {
		FrameOP frame = new Frame(i,9);
		for( vector1<FragDataOP>::iterator fi = fds9.begin(); fi != fds9.end(); ++fi ) {
			frame->add_fragment(*fi);
		}
		frags->add(frame);
	}
	if(frags->size() == 0) return NULL;
	return frags;
}

class ChiMover : public protocols::moves::Mover {
	Size residue_,nchi;
	core::pack::dunbrack::SingleResidueRotamerLibraryCAP lib_;
	core::pack::dunbrack::RotamerLibraryScratchSpace scratch_;
public:
	ChiMover(core::pose::Pose const & pose, Size residue) : residue_(residue) {
		using namespace core::pack::dunbrack;
		using namespace core::chemical;
		string name3 = pose.residue(residue).name3();
		if(name3=="BPY") name3 = "TRP"; //TODO find a better way to get rotamers...
		if(name3=="NPP") name3 = "TYR";
		if(name3=="NPD") name3 = "TYR";
		if(name3=="NPA") name3 = "TYR";
		if(name3=="TIA") name3 = "LYS";
		if(name3=="HSC") name3 = "HIS";
		ResidueTypeSetCAP rs( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		lib_ = core::scoring::ScoringManager::get_instance()->get_RotamerLibrary().get_rsd_library( rs->name_map(name3) );
		assert(lib_);
		nchi = pose.residue(residue_).nchi();
	}
	void apply( core::pose::Pose & pose ) {
		using namespace core::pack::dunbrack;
		ChiVector chis;
		lib_->assign_random_rotamer_with_bias(pose.residue(residue_),scratch_,RG,chis,true);
		for(Size i = 1; i <= chis.size(); ++i)	{
			pose.set_chi(i,residue_,chis[i]);
		}
		for(Size i = chis.size()+1; i <= nchi; ++i ) {
			pose.set_chi(i,residue_,uniform()*360.0);
		}
		// Real a =  5.0*numeric::random::gaussian();
		// Real b = 30.0*numeric::random::gaussian();
		// Real c = 30.0*numeric::random::uniform();
		// Real d = 30.0*numeric::random::uniform();
		// pose.set_chi(1,1, pose.chi(1,1) + a );
		// pose.set_chi(2,1, pose.chi(2,1) + b );
		// pose.set_chi(3,1, pose.chi(3,1) + c );
		// pose.set_chi(4,1, pose.chi(4,1) + d );
		// pose.set_psi(2,pose.psi(2)+10.0*numeric::random::gaussian());
	}
	std::string get_name() const { return "ChiMover"; }
	void
	parse_my_tag(
		utility::Tag::TagCOP const,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const &
	) {}
};

class MySlideMover : public protocols::moves::Mover {
	std::map<Size,core::conformation::symmetry::SymDof> dofs_;
	protocols::moves::MoverOP slide_;
public:
	MySlideMover(core::pose::Pose & pose) {
		dofs_ = core::conformation::symmetry::symmetry_info(pose)->get_dofs();
		slide_ = new protocols::symmetric_docking::SymDockingSlideIntoContact(dofs_);
	}
	void apply( core::pose::Pose & pose ) {
		for(std::map<Size,core::conformation::symmetry::SymDof>::iterator i = dofs_.begin(); i != dofs_.end(); ++i) {
			// TR << "MySlideMover: set jump " << i->first << " translation to 0" << std::endl;
			core::kinematics::Jump j(pose.jump(i->first));
			// TR << j << std::endl;
			j.set_translation(Vec(0,0,0));
			// TR << j << std::endl;
			pose.set_jump(i->first,j);
		}
		if(dofs_.size()) slide_->apply(pose);
	}
	std::string get_name() const { return "MySlideMover"; }
};

class MyTransMover : public protocols::moves::Mover {
	std::map<Size,core::conformation::symmetry::SymDof> dofs_;
public:
	PoseWrapAP pw_;
	Real mag_,postfrac_;
	bool concerted_,debug_;
	MoverOP post_;
	Size concert_sub_;
	vector1<int> symjumps,rbjumps,comjump;
	MyTransMover(PoseWrap & pw, Real mag, bool concerted=false, MoverOP post = NULL, Real postfrac = 0, Size concert_sub=0) : mag_(mag), postfrac_(postfrac), concerted_(concerted), debug_(pw.debug), post_(post), concert_sub_(concert_sub) {
		pw_ = pw;
		core::pose::Pose & pose(pw.pose);
		dofs_ = core::conformation::symmetry::symmetry_info(pose)->get_dofs();
		if(concerted_) {
			for(Size j = 1; j <= pose.fold_tree().num_jump(); ++j) {
				if( pose.fold_tree().downstream_jump_residue(j)==(int)concert_sub_ ) concert_sub_ = pose.fold_tree().upstream_jump_residue(j);
			}
		}
		for(std::map<Size,core::conformation::symmetry::SymDof>::iterator i = dofs_.begin(); i != dofs_.end(); ++i) {
			if(        i->second.allow_dof(1) &&  i->second.allow_dof(2) &&  i->second.allow_dof(3) &&
			          !i->second.allow_dof(4) && !i->second.allow_dof(5) && !i->second.allow_dof(6) ) {
				comjump.push_back(i->first);
			} else if( i->second.allow_dof(1) && !i->second.allow_dof(2) && !i->second.allow_dof(3) &&
			           i->second.allow_dof(4) && !i->second.allow_dof(5) && !i->second.allow_dof(6) ) {
				symjumps.push_back(i->first);
			} else if(/*i->second.allow_dof(1)&&*/i->second.allow_dof(2) && i->second.allow_dof(3) &&
			            i->second.allow_dof(4) && i->second.allow_dof(5) && i->second.allow_dof(6) ) {
				rbjumps.push_back(i->first);
			}
		}
		if(comjump.size() > 1) utility_exit_with_message("more than one COM jump!");
	}
	std::string get_name() const { return "MyTransMover"; }
	void apply( core::pose::Pose & pose ) {
		pw_->validate_set_reference();
		// if(debug_) TR << "MyTransMover " << concerted_ << std::endl;
		if( concerted_ ) {
			// pose.dump_pdb("trans_before.pdb");
			if( (rbjumps.size()==0 || uniform() < 0.6) && symjumps.size() > 0 ) {
				Real magx = gaussian()*mag_;
				for(Size idx = 1; idx <= symjumps.size(); ++idx) {
					core::kinematics::Jump j(pose.jump(symjumps[idx]));
					j.set_translation(j.get_translation()+(magx)*Vec(1,0,0));
					pose.set_jump(symjumps[idx],j);
				}
			} else {
				Real magx = gaussian()*mag_; Real magy = gaussian()*mag_; Real magz = gaussian()*mag_;
				for(Size idx = 1; idx <= rbjumps.size(); ++idx) {
					core::kinematics::Jump j(pose.jump(rbjumps[idx]));
					if(dofs_[rbjumps[idx]].allow_dof(1)) j.set_translation(j.get_translation()+(magx)*Vec(1,0,0));
					if(dofs_[rbjumps[idx]].allow_dof(2)) j.set_translation(j.get_translation()+(magy)*Vec(0,1,0));
					if(dofs_[rbjumps[idx]].allow_dof(3)) j.set_translation(j.get_translation()+(magz)*Vec(0,0,1));
					pose.set_jump(rbjumps[idx],j);
				}
			}
			// pose.dump_pdb("trans_after.pdb");
			// if(uniform() < 0.01) std::exit(-1);
		} else {
			if( (rbjumps.size()==0 || uniform() < 0.6) && symjumps.size() > 0 ) {
				Size idx = std::ceil(uniform()*symjumps.size());
				core::kinematics::Jump j(pose.jump(symjumps[idx]));
				j.set_translation(j.get_translation()+(gaussian()*mag_)*Vec(1,0,0));
				pose.set_jump(symjumps[idx],j);
			} else if(rbjumps.size()>0) { // move rb jump
				Size idx = std::ceil(uniform()*rbjumps.size());
				core::kinematics::Jump j(pose.jump(rbjumps[idx]));
				j.set_translation(j.get_translation()+Vec((gaussian()*mag_),(gaussian()*mag_),(gaussian()*mag_)));
				pose.set_jump(rbjumps[idx],j);
			} else {
				static bool printed(false);
				if(!printed) {
					TR << "MyTransMover WARNING! no jump with translation MOVES!" << std::endl;
					printed=true;
				}
			}
		}

		if(post_ && uniform() < postfrac_ && !concerted_) post_->apply(pose);

		pw_->validate();
		pw_->center();
	}
};

void get_res_downstream_of_jump(FoldTree const & ft, Size j, vector1<Size> & resout) {
	// TR << "get_res_downstream_of_jump " << resout << std::endl;
	bool isleaf = true;
	for(Size i = 1; i <= ft.num_jump(); ++i) {
		if(ft.upstream_jump_residue(i)==ft.downstream_jump_residue(j)) {
			// TR << "get_res_downstream_of_jump " << j << " recursing on " << i << std::endl;
			get_res_downstream_of_jump(ft,i,resout);
			isleaf = false;
		}
	}
	if(isleaf) {
		resout.push_back(ft.downstream_jump_residue(j));
		// TR << "get_res_downstream_of_jump IS LEAF " << j << " " << resout->size() << std::endl;
	}
}

class MyRotMover : public protocols::moves::Mover {
	PoseWrapAP pw_;
	std::map<Size,core::conformation::symmetry::SymDof> dofs_;
	Real mag_,postfrac_;
	bool counterrot_,concerted_,debug_;
	MoverOP post_;
	Size concert_sub_;
	vector1<int> symjumps,rbjumps,rotjumps,comjumps;
	vector1<Size> rbparentres,starts_,stops_;
public:
	std::string get_name() const { return "MyRotMover"; }
	MyRotMover(PoseWrap & pw, Real mag, bool counterrot=false, bool concerted=false, MoverOP post = NULL, Real postfrac=0, Size concert_sub=0) : mag_(mag), postfrac_(postfrac), counterrot_(counterrot), concerted_(concerted), debug_(pw.debug), post_(post), concert_sub_(concert_sub) {
		pw_ = pw;
		dofs_ = core::conformation::symmetry::symmetry_info(pw.pose)->get_dofs();
		if(counterrot_) {
			for(Size j = 1; j <= pw.pose.fold_tree().num_jump(); ++j) {
				if( pw.pose.fold_tree().downstream_jump_residue(j)==(int)concert_sub_ ) concert_sub_ = pw.pose.fold_tree().upstream_jump_residue(j);
			}
		}
		for(std::map<Size,core::conformation::symmetry::SymDof>::iterator i = dofs_.begin(); i != dofs_.end(); ++i) {
			if( i->second.allow_dof(1) && !i->second.allow_dof(2) && !i->second.allow_dof(3) &&
			    i->second.allow_dof(4) && !i->second.allow_dof(5) && !i->second.allow_dof(6) ) {
				symjumps.push_back(i->first);
			} else if( !i->second.allow_dof(1) && !i->second.allow_dof(2) && !i->second.allow_dof(3) &&
			            i->second.allow_dof(4) && !i->second.allow_dof(5) && !i->second.allow_dof(6) ) {
				rotjumps.push_back(i->first);
			} else {
			  	rbjumps.push_back(i->first);
				rbparentres.push_back(pw.pose.fold_tree().downstream_jump_residue(i->first));
				// TR << "ROT RBjump " << i->first << " " << pw.pose.fold_tree().downstream_jump_residue(i->first) << std::endl;
				vector1<Size> downstream;
				get_res_downstream_of_jump(pw.pose.fold_tree(),i->first,downstream);
				for(vector1<Size>::iterator ir = downstream.begin(); ir != downstream.end(); ++ir) {
					// TR << "checking res " << *ir << std::endl;
					if( *ir <= pw.nres*pw.nsub ) {
						starts_.push_back(pw.subsub_starts_[pw.which_subsub(*ir)]);
						stops_ .push_back(pw.subsub_ends_  [pw.which_subsub(*ir)]);
					}
				}
			}
		}
		if(starts_.size() == 0 && rbjumps.size() != 0) {
			utility_exit_with_message("have RB jump with no downstream residues!");
		}

	}
	void apply( core::pose::Pose & pose ) {
		pw_->validate_set_reference();
		// if(debug_) TR << "MyRotMover " << concert_sub_ << " " << concerted_ << " " << counterrot_ << std::endl;
		// TR << "ROT MOVER APPLY" << std::endl;
		if( !concerted_ && !counterrot_ ) {
			Real MAG = (gaussian()*mag_);
			if( ( (symjumps.size()==0 && rbjumps.size()==0) || uniform() < 0.333) && rotjumps.size() > 0  ) {
				// TR << "ROT not concerted rotjumps" << std::endl;
					Size idx = std::ceil(uniform()*rotjumps.size());
					core::kinematics::Jump j(pose.jump(rotjumps[idx]));
					j.set_rotation(x_rotation_matrix_degrees(MAG)*j.get_rotation());
					pose.set_jump(rotjumps[idx],j);
					// TR << "ROT MOVER MOVE 1 " << MAG << std::endl;
				} else if( (rbjumps.size()==0 || uniform() < 0.5) && symjumps.size() > 0 ) {
				// TR << "ROT not concerted symjumps" << std::endl;
				Size idx = std::ceil(uniform()*symjumps.size());
				core::kinematics::Jump j(pose.jump(symjumps[idx]));
				j.set_rotation(x_rotation_matrix_degrees(MAG)*j.get_rotation());
				pose.set_jump(symjumps[idx],j);
				// TR << "ROT MOVER MOVE 2 X " << MAG << std::endl;
			} else {
				// TR << "ROT not concerted rbjumps" << std::endl;
				Size idx = std::ceil(uniform()*rbjumps.size());
				core::kinematics::Jump j(pose.jump(rbjumps[idx]));
				if(dofs_[rbjumps[idx]].allow_dof(4)) j.set_rotation(x_rotation_matrix_degrees(gaussian()*mag_)*j.get_rotation());
				if(dofs_[rbjumps[idx]].allow_dof(5)) j.set_rotation(y_rotation_matrix_degrees(gaussian()*mag_)*j.get_rotation());
				if(dofs_[rbjumps[idx]].allow_dof(6)) j.set_rotation(z_rotation_matrix_degrees(gaussian()*mag_)*j.get_rotation());
				pose.set_jump(rbjumps[idx],j);
			}
		} else if(rbjumps.size() > 0 && concerted_ ) {
			// pose.dump_pdb("rot_before.pdb");c
			// for(Size i = 1; i <= rbparentres.size(); ++i) std::cerr << "CONCERTED rbparentres: " << rbparentres[i] << std::endl;
			// TR << "ROT concerted_ switch=" << concert_sub_ << std::endl;
			Real magx = (gaussian()*mag_); Real magy = (gaussian()*mag_); Real magz = (gaussian()*mag_);
			Vec CEN;
			{
				// pose.dump_pdb("rot_before.pdb");
				Vec tmp(0,0,0); // TODO: set to center of subunit!!!!
				Size n = 0;
				for(Size idx = 1; idx <= starts_.size(); ++idx) {
					// TR << "ROT computing CEN from " << starts_[idx] << " " << stops_[idx] << std::endl;
					for(Size i = starts_[idx]; i <= stops_[idx]; ++i) {
						tmp += pose.xyz(AtomID(2,i));
						n++;
					}
				}
				CEN = Vec(tmp.x()/n,tmp.y()/n,tmp.z()/n); // TODO: FIX THIS HACK!!!!!!!!!!! assumes z,x,y coord frame in sym file!!!!!!!
				// Vec CEN(0,0,0);
				// pose.dump_pdb("rot_before.pdb");
				// std::ofstream out1("rot_before.pdb",std::ios::app);
				// out1 << "HETATM 9999 ZN    ZN A 999      " << LJ(8,F(5,3,CEN.x())) << LJ(8,F(5,3,CEN.y())) << LJ(8,F(5,3,CEN.z())) << "1.00  0.00" << std::endl;
				// std::cerr << "CEN  " << CEN << std::endl;
				// for(Size idx = 1; idx <= rbjumps.size(); ++idx) {
				// 	Vec LCEN = pose.xyz(AtomID(1,rbparentres[idx]));
				// 	std::cerr << "LCEN " << LCEN << std::endl;
				// 	out1 << "HETATM 9999 ZN    ZN A "+string_of(999-idx)+"      " << LJ(8,F(5,3,LCEN.x())) << LJ(8,F(5,3,LCEN.y())) << LJ(8,F(5,3,LCEN.z())) << "1.00  0.00" << std::endl;
				// }
			}

			// numeric::xyzVector<Real> com(0,0,0);
			// if(comjumps.size()) com = pose.jump(comjumps[1]).get_translation();
			for(Size idx = 1; idx <= rbjumps.size(); ++idx) {
				core::kinematics::Jump j(pose.jump(rbjumps[idx]));
				// std::cerr << "JMP: " << idx << " " << j.get_translation() << " " << pose.xyz(AtomID(1,rbparentres[idx])) << std::endl;
				// if( (j.get_translation()+com-pose.xyz(AtomID(1,rbparentres[idx]))).length() > 0.0001 ) {
				// 	std::cerr << "COM " << j.get_translation() << " " << pose.xyz(AtomID(1,rbparentres[idx])) << " " << com << std::endl;
				// 	utility_exit_with_message("coordinate frames do not match!!!");
				// }
				// std::cerr << "ROT: " << j.get_rotation() << std::endl;
				if(dofs_[rbjumps[idx]].allow_dof(1)&&dofs_[rbjumps[idx]].allow_dof(4)) j.set_translation(x_rotation_matrix_degrees(magx)*(j.get_translation()-CEN)+CEN);
				if(dofs_[rbjumps[idx]].allow_dof(2)&&dofs_[rbjumps[idx]].allow_dof(5)) j.set_translation(y_rotation_matrix_degrees(magy)*(j.get_translation()-CEN)+CEN);
				if(dofs_[rbjumps[idx]].allow_dof(3)&&dofs_[rbjumps[idx]].allow_dof(6)) j.set_translation(z_rotation_matrix_degrees(magz)*(j.get_translation()-CEN)+CEN);
				if(dofs_[rbjumps[idx]].allow_dof(1)&&dofs_[rbjumps[idx]].allow_dof(4)) j.set_rotation   (x_rotation_matrix_degrees(magx)*(j.get_rotation()));
				if(dofs_[rbjumps[idx]].allow_dof(2)&&dofs_[rbjumps[idx]].allow_dof(5)) j.set_rotation   (y_rotation_matrix_degrees(magy)*(j.get_rotation()));
				if(dofs_[rbjumps[idx]].allow_dof(3)&&dofs_[rbjumps[idx]].allow_dof(6)) j.set_rotation   (z_rotation_matrix_degrees(magz)*(j.get_rotation()));
				pose.set_jump(rbjumps[idx],j);
			}

			// {
			// 	pose.dump_pdb("rot_after.pdb");
			// 	Vec tmp(0,0,0); // TODO: set to center of subunit!!!!
			// 	Size n = 0;
			// 	for(Size idx = 1; idx <= starts_.size(); ++idx) {
			// 		TR << "ROT computing CEN from " << starts_[idx] << " " << stops_[idx] << std::endl;
			// 		for(Size i = starts_[idx]; i <= stops_[idx]; ++i) {
			// 			tmp += pose.xyz(AtomID(2,i));
			// 			n++;
			// 		}
			// 	}
			// 	Vec CEN = Vec(tmp.x()/n,tmp.y()/n,tmp.z()/n); // TODO: FIX THIS HACK!!!!!!!!!!! assumes z,x,y coord frame in sym file!!!!!!!
			// 	pose.dump_pdb("rot_after.pdb");
			// 	std::ofstream out("rot_after.pdb",std::ios::app);
			// 	out << "HETATM 9999 ZN    ZN A 999      " << LJ(8,F(5,3,CEN.x())) << LJ(8,F(5,3,CEN.y())) << LJ(8,F(5,3,CEN.z())) << "1.00  0.00" << std::endl;
			// 	for(Size idx = 1; idx <= rbjumps.size(); ++idx) {
			// 		Vec LCEN = pose.xyz(AtomID(1,rbparentres[idx]));
			// 		out << "HETATM 9999 ZN    ZN A "+string_of(999-idx)+"      " << LJ(8,F(5,3,LCEN.x())) << LJ(8,F(5,3,LCEN.y())) << LJ(8,F(5,3,LCEN.z())) << "1.00  0.00" << std::endl;
			// 	}
			// }

			// pose.dump_pdb("rot_after.pdb");
			// std::exit(-1);
		} else if(rbjumps.size() > 0 && counterrot_) {
			// TR << "ROT counterrot" << std::endl;
			Real magx = (gaussian()*mag_); Real magy = (gaussian()*mag_); Real magz = (gaussian()*mag_);
			Size idx = std::ceil(uniform()*rbjumps.size());
			core::kinematics::Jump j(pose.jump(rbjumps[idx]));
			if(dofs_[rbjumps[idx]].allow_dof(4)) j.set_rotation(x_rotation_matrix_degrees(magx)*j.get_rotation());
			if(dofs_[rbjumps[idx]].allow_dof(5)) j.set_rotation(y_rotation_matrix_degrees(magy)*j.get_rotation());
			if(dofs_[rbjumps[idx]].allow_dof(6)) j.set_rotation(z_rotation_matrix_degrees(magz)*j.get_rotation());
			pose.set_jump(rbjumps[idx],j);
			// TR << "ROT MOVER MOVE 4 " << magx << " " << magy << " " << magz << std::endl;
			if(rbjumps.size() > 1 && counterrot_) {
				Size idx2 = std::ceil(uniform()*rbjumps.size());
				while(idx2 == idx) idx2 = std::ceil(uniform()*rbjumps.size());
				core::kinematics::Jump j(pose.jump(rbjumps[idx2]));
				if(dofs_[rbjumps[idx2]].allow_dof(4)) j.set_rotation(x_rotation_matrix_degrees(magx)*j.get_rotation());
				if(dofs_[rbjumps[idx2]].allow_dof(5)) j.set_rotation(y_rotation_matrix_degrees(magy)*j.get_rotation());
				if(dofs_[rbjumps[idx2]].allow_dof(6)) j.set_rotation(z_rotation_matrix_degrees(magz)*j.get_rotation());
				pose.set_jump(rbjumps[idx2],j);
				// TR << "ROT MOVER MOVE 5 " << magx << " " << magy << " " << magz << std::endl;
			}
		} else {
			static bool printed(false);
			if(!printed) {
				TR << "MyRotMover WARNING! no jump with allowed rotation MOVES! counterrot_ " << counterrot_ << " rbj: " << rbjumps.size() << " concerted_: " << concerted_ << std::endl;
				printed=true;
			}
		}
		if(post_ && uniform() < postfrac_ && !concerted_) post_->apply(pose);

		pw_->validate();

	}
};


class FloatScRotMover : public protocols::moves::Mover {
	vector1<Size> floating_scs_from_pw_;
	Size nres_from_pw_;
public:
	FloatScRotMover(PoseWrap & pw) : floating_scs_from_pw_(pw.floating_scs_), nres_from_pw_(pw.nres) {}
	void apply( core::pose::Pose & pose ) {
		Size r = floating_scs_from_pw_[ std::ceil(uniform()*floating_scs_from_pw_.size()) ];
		change_floating_sc_geometry(pose,r,nres_from_pw_);
	}
	std::string get_name() const { return "FloatScRotMover"; }
};


void print_movemap(core::kinematics::MoveMap const & movemap) {
	using namespace core::id;
	using namespace core::kinematics;
	TR << "movemap " << std::endl;
	for(std::map< TorsionType, bool >::const_iterator i = movemap.torsion_type_begin(); i != movemap.torsion_type_end(); ++i) {
		TR << "TorsionType " << i->first << " " << i->second << std::endl;
	}
	for(std::map< std::pair< Size, TorsionType >, bool >::const_iterator i = movemap.movemap_torsion_id_begin(); i != movemap.movemap_torsion_id_end(); ++i) {
		TR << "MoveMapTorsionID (" << i->first.first << "," << i->first.second << ") " << i->second << std::endl;
	}
	for(std::map< TorsionID, bool >::const_iterator i = movemap.torsion_id_begin(); i != movemap.torsion_id_end(); ++i) {
		TR << "TorsionID " << i->first << " " << i->second << std::endl;
	}
	for(std::map< DOF_Type, bool >::const_iterator i = movemap.dof_type_begin(); i != movemap.dof_type_end(); ++i) {
		TR << "DOF_Type " << i->first << " " << i->second << std::endl;
	}
	for(std::map< DOF_ID, bool >::const_iterator i = movemap.dof_id_begin(); i != movemap.dof_id_end(); ++i) {
		TR << "DOF_ID " << i->first << " " << i->second << std::endl;
	}
	for(std::map< JumpID, bool >::const_iterator i = movemap.jump_id_begin(); i != movemap.jump_id_end(); ++i) {
		TR << "JumpID " << i->first << " " << i->second << std::endl;
	}

}

void jumps_to_root(Pose & pose, Size rsd, vector1<Size> & jumps) {
	FoldTree const & ft(pose.fold_tree());
	if( (int)rsd == ft.root() ) return;
	for(Size i = 1; i <= ft.num_jump(); ++i) {
		// TR << "jumphunt " << i << " " << ft.downstream_jump_residue(i) << " " << ft.upstream_jump_residue(i) << " " << rsd << std::endl;
		if(ft.downstream_jump_residue(i)==(int)rsd) {
			jumps.push_back(i);
			return jumps_to_root(pose,ft.upstream_jump_residue(i),jumps);
		}
	}
	utility_exit_with_message("jumps_to_root: not a jump residue!");
}

MoverOP
make_float_sc_min_mover(core::pose::Pose & pose, vector1<Size> rsd_, Real tol) {
	vector1<Size> jumps;
	for(Size i = 1; i <= rsd_.size(); ++i) jumps_to_root(pose,rsd_[i],jumps);

	ScoreFunctionOP sf_ = new core::scoring::ScoreFunction;
	sf_ = new core::scoring::symmetry::SymmetricScoreFunction(*sf_);

	// sf_->set_weight(core::scoring::floating_sc,0.01);
	sf_->set_weight(core::scoring::atom_pair_constraint,1.0);
	// sf_->set_weight(core::scoring::vdw,1.0);
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_chi(false);
	movemap->set_bb(false);
	movemap->set_jump(false);
	for(Size i = 1; i <= rsd_.size(); ++i) movemap->set_chi(rsd_[i],true);
	for(Size i = 1; i <= rsd_.size(); ++i) movemap->set_bb(rsd_[i],true);
	// TR << pose.fold_tree() << std::endl;
	for(Size i = 1; i <= jumps.size(); ++i) {
		// TR << "JUMP " << jumps[i] << std::endl;
		movemap->set_jump(jumps[i],true);
	}
	// print_movemap(*movemap);
	core::conformation::symmetry::make_symmetric_movemap( pose, *movemap );
	// print_movemap(*movemap);
	// std::exit(-1);
	MoverOP min_ = new protocols::simple_moves::symmetry::SymMinMover( movemap, sf_, "dfpmin_armijo_nonmonotone", tol, true, false, false );
	return min_;
}

class FloatScRandomChi : public protocols::moves::Mover {
	vector1<Size> rsd_;
public:
	FloatScRandomChi(vector1<Size> rsd) : rsd_(rsd) {}
	void apply(Pose & pose) {
		Size i = rsd_[ std::ceil(uniform()*rsd_.size()) ];
		Size j = rsd_[ std::ceil(uniform()*pose.residue(i).nchi()) ];
		pose.set_chi(j,i,uniform()*360.0);
	}
	std::string get_name() const { return "FloatScRandomChi"; }
};

class FloatScMonteCarlo : public protocols::moves::Mover {
	vector1<Size> rsd_;
	Size N_;
public:
	FloatScMonteCarlo(vector1<Size> rsd, Size N) : rsd_(rsd), N_(N) {}
	void apply(Pose & pose) {
		using namespace protocols::moves;
		core::scoring::symmetry::SymmetricScoreFunction sf;
		sf.set_weight(core::scoring::atom_pair_constraint,1.0);
		Real temp = 1.0;
		MonteCarloOP mc = new MonteCarlo( pose, sf, temp );
		mc->set_autotemp( true, temp );
		mc->set_temperature( temp );
		TR << "floating_scs_monte_carlo start " << sf(pose) << std::endl;
		RepeatMover( new TrialMover( new FloatScRandomChi(rsd_), mc ), N_ ).apply( pose );
		TR << "floating_scs_monte_carlo end " << sf(pose) << std::endl;

	}
	std::string get_name() const { return "FloatScMonteCarlo"; }
};

class ScMinMover : public protocols::moves::Mover {
	MoverOP min_;
	ScoreFunctionOP sf_;
	vector1<Size> rsd_;
	Size nres_;
public:
	ScMinMover(core::pose::Pose & pose, vector1<Size> rsd, Real tol, Size nres) : rsd_(rsd),nres_(nres) {
		min_ = make_float_sc_min_mover(pose,rsd,tol);
		ScoreFunctionOP sf = new core::scoring::ScoreFunction;
		sf_ = new core::scoring::symmetry::SymmetricScoreFunction(*sf);
		// sf_->set_weight(core::scoring::floating_sc,1.0);
		sf_->set_weight(core::scoring::atom_pair_constraint,1.0);
	}
	void apply(core::pose::Pose & pose) {
		// TR << "before ScMinMover " << (*sf_)(pose) << std::endl;
		for(Size i = 1; i<=rsd_.size();++i) {
			Pose tmp = pose;
			change_floating_sc_geometry(pose,rsd_[i],nres_);
			if((*sf_)(tmp) < (*sf_)(pose)) {
				pose = tmp;
			} else {
				TR << "ScMinMover: performed sc swap on res " << i << std::endl;
			}
		}
		// FloatScMonteCarlo(rsd_,100).apply(pose);
		min_->apply(pose);
		// TR << "after ScMinMover  " << (*sf_)(pose) << std::endl;
	}
	std::string get_name() const { return "ScMinMover"; }
};


class BBMover : public protocols::moves::Mover {
	Size nres_;
	Real mag_;
public:
	BBMover(Size nres, Real mag) : nres_(nres),mag_(mag) {}
	void apply(core::pose::Pose & pose) {
		Size i = std::ceil(uniform()*nres_);
		if(uniform()<0.5) pose.set_phi(i,pose.phi(i)+gaussian()*mag_);
		else              pose.set_psi(i,pose.psi(i)+gaussian()*mag_);
	}
	std::string get_name() const { return "BBMover"; }
};

void repack(PoseWrap & pw, ScoreFunctionOP sf, bool fsc_only=false) {
	using namespace core::pack::task;
	core::pose::Pose & pose(pw.pose);
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	task->restrict_to_repacking();
	task->or_include_current(true);
	pw.update_designable_packable();
	for(Size i = 1; i <= pw.nres; ++i) {
		if(!pw.is_rsd_repackable(i)) {
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}
	if(fsc_only) {
		for(Size i = 1; i <= pw.nres; ++i) {
			if(std::find(pw.attached_scattach_res_.begin(),pw.attached_scattach_res_.end(),i)!=pw.attached_scattach_res_.end()) continue;
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}
	for(Size i = 1; i <= pw.attached_scattach_res_.size(); ++i) {
		TR << "scattach res packing mode: " << pw.attached_scattach_res_[i] << std::endl;
		task->nonconst_residue_task(pw.attached_scattach_res_[i]).or_ex1(true);
		task->nonconst_residue_task(pw.attached_scattach_res_[i]).or_ex2(true);
		task->nonconst_residue_task(pw.attached_scattach_res_[i]).or_ex1_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		if(fsc_only) task->nonconst_residue_task(pw.attached_scattach_res_[i]).or_ex2_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		task->nonconst_residue_task(pw.attached_scattach_res_[i]).and_extrachi_cutoff(0);
		task->nonconst_residue_task(pw.attached_scattach_res_[i]).restrict_to_repacking();
	}
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);
}

void design(PoseWrap & pw, ScoreFunctionOP sf, bool do_trp = false) {
	core::pose::Pose & pose(pw.pose);
	using namespace core::pack::task;
	// vector1< bool > aas(20,true);
	// aas[core::chemical::aa_cys] = false;
	// aas[core::chemical::aa_his] = false;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	task->or_include_current(true);
	// task->nonconst_residue_task(pw.attach).prevent_repacking();
	// for(Size i = 1; i <= task->total_residue(); ++i) {
	// 	task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
	// }

	vector1< bool > aas(20,true);
	// aas[core::chemical::aa_his] = false;
	aas[core::chemical::aa_pro] = false;
	aas[core::chemical::aa_cys] = false;
	aas[core::chemical::aa_met] = false;
	aas[core::chemical::aa_lys] = false;
	if(pw.pose.is_fullatom()) {
		if( basic::options::option[basic::options::OptionKeys::smhybrid::design_hydrophobic]() ) {
			for(Size i = 1; i <= aas.size(); ++i) aas[i] = false;
			// aas[core::chemical::aa_ala] = true;
			aas[core::chemical::aa_ile] = true;
			aas[core::chemical::aa_val] = true;
			aas[core::chemical::aa_leu] = true;
			aas[core::chemical::aa_phe] = true;
			aas[core::chemical::aa_tyr] = true;
			aas[core::chemical::aa_met] = true;
			aas[core::chemical::aa_pro] = true;
			if(do_trp) aas[core::chemical::aa_trp] = true;
		}
	} else {
		utility_exit_with_message("design() called on cen pose... use cen_design()?");
	}

	pw.update_designable_packable();
	for(Size i = 1; i <= pw.nres; ++i) {
		if(!pw.is_rsd_repackable(i)) {
			TR << "fixed     " << i << std::endl;
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_ex1_sample_level( core::pack::task::EX_ONE_STDDEV );
			task->nonconst_residue_task(i).or_ex2_sample_level( core::pack::task::EX_ONE_STDDEV );
		} else if(!pw.is_rsd_designable(i)) {
			TR << "repacking " << i << std::endl;
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_ex1_sample_level( core::pack::task::EX_ONE_STDDEV );
			task->nonconst_residue_task(i).or_ex2_sample_level( core::pack::task::EX_ONE_STDDEV );
		} else {
			TR << "designing " << i << std::endl;
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
		}
	}
	// std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	// for(Size i = 1; i <= pw.nres; ++i) {
	// 	if(std::find(pw.attached_scattach_res_.begin(),pw.attached_scattach_res_.end(),i)!=pw.attached_scattach_res_.end()) continue;
	// 	std::cerr << "SLKDFJSKLDJ " << i << std::endl;
	// 	task->nonconst_residue_task(i).prevent_repacking();
	// }
	for(Size i = 1; i <= pw.attached_scattach_res_.size(); ++i) {
		TR << "scattach res packing mode: " << pw.attached_scattach_res_[i] << std::endl;
		task->nonconst_residue_task(pw.attached_scattach_res_[i]).or_ex1(true);
		task->nonconst_residue_task(pw.attached_scattach_res_[i]).or_ex2(true);
		task->nonconst_residue_task(pw.attached_scattach_res_[i]).or_ex1_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		// task->nonconst_residue_task(pw.attached_scattach_res_[i]).or_ex2_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		task->nonconst_residue_task(pw.attached_scattach_res_[i]).and_extrachi_cutoff(0);
		task->nonconst_residue_task(pw.attached_scattach_res_[i]).restrict_to_repacking();
	}
	// steal Nobu's smart rules about what can do where...
	protocols::flxbb::DesignLayerOperationOP op = new protocols::flxbb::DesignLayerOperation(true,true,true);
	op->apply(pose,*task);
	// for(Size i = 1; i <= pw.nres; ++i) {
	// 	task->nonconst_residue_task(i).allow_aa(pose.residue(i).aa());
	// }


	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	TR << *task << std::endl;
	repack.apply(pose);
	pw.check_scattach_res();
}


void cen_design(PoseWrap & pw, ScoreFunctionOP sf) {
	core::pose::Pose & pose(pw.pose);
	using namespace core::pack::task;
	vector1< bool > aas(20,false);
	aas[core::chemical::aa_ile] = true;
	aas[core::chemical::aa_val] = true;
	aas[core::chemical::aa_leu] = true;
	aas[core::chemical::aa_phe] = true;
	// vector1< bool > aas(20,true);
	// aas[core::chemical::aa_cys] = false;
	// aas[core::chemical::aa_met] = false;
	// aas[core::chemical::aa_lys] = false;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	task->or_include_current(true);
	for(Size i = 1; i <= pw.nres; ++i) {
		task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
	}
	for(Size i = 1; i <= pw.floating_scs_.size(); ++i) {
		task->nonconst_residue_task(pw.floating_scs_[i]).prevent_repacking();
	}
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);
	pw.check_scattach_res();
}


void minimize(PoseWrap & pw, ScoreFunctionOP sf, int bb=0, bool jmp=true, bool chi1=false) {
	core::pose::Pose & pose(pw.pose);
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	// core::conformation::symmetry::make_symmetric_movemap(pose,*movemap);
	movemap->set_chi(true);
	if(pw.pose.residue(1).name3()=="NPA") { movemap->set_chi(1,chi1); }//movemap->set(DOF_ID(AtomID(pose.residue(1).atom_index("C11"),1),core::id::PHI),true); }
	if(pw.pose.residue(1).name3()=="NPD") { movemap->set_chi(1,chi1); }//movemap->set(DOF_ID(AtomID(pose.residue(1).atom_index("C11"),1),core::id::PHI),true); }
	if(pw.pose.residue(1).name3()=="NPH") { movemap->set_chi(1,chi1); }//movemap->set(DOF_ID(AtomID(pose.residue(1).atom_index("C11"),1),core::id::PHI),true); }
	movemap->set_bb(false);
	movemap->set_jump(jmp);
	// if(bb == 1) for(Size i = 1; i <= pw.nres; ++i) if(pw.pose.secstruct(i)=='L') movemap->set_bb(i,true);
	if(bb > 0) { for(Size i = 1; i <= pw.nres; ++i) movemap->set_bb(i,true); }
	else       { for(Size i = 1; i <= pw.nres; ++i) movemap->set_bb(i,false); }
	// if(bb >= 2) for(Size i = 1; i <= pw.nres; ++i) movemap->set_bb(i,true);
	if(bb > 0) TR << "minimizing with flexible bb! " << bb << std::endl;


	// set jumps, chi, bb true for floating sc res always
	vector1<Size> jumps;
	for(Size i = 1; i <= pw.floating_scs_.size(); ++i) jumps_to_root(pose,pw.floating_scs_[i],jumps);
	for(Size i = 1; i <= pw.floating_scs_.size(); ++i ) {
		TR << "minbb floating sc res " << pw.floating_scs_[i] << std::endl;
		movemap->set_bb(pw.floating_scs_[i],true);
	}
	for(Size i = 1; i <= jumps.size(); ++i) movemap->set_jump(jumps[i],true);
	// movemap->set_jump(jumps.back(),false);

	// // set bb free for linker always
	// for(Size i = 1; i <= pw.linker_res_.size(); ++i) {
	// 	TR << "minbb linker res " << pw.linker_res_[i] << std::endl;
	// 	movemap->set_bb(pw.linker_res_[i],true);
	// }

	print_movemap(*movemap);
	TR << "///////////////////////////////" << std::endl;
	core::conformation::symmetry::make_symmetric_movemap( pose, *movemap );
	print_movemap(*movemap);
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );

	// TR << "movemap " << std::endl;
	// for(std::map< DOF_ID, bool >::const_iterator i = movemap->dof_id_begin(); i != movemap->dof_id_end(); ++i) {
		// TR << "DOF " << i->first << " " << i->second << std::endl;
	// }
	// std::exit(-1);
	// TR << "minimize before " << (*sf)(pose) << std::endl;

	// if(pw.floating_scs_.size()) {
	// 	FloatScMonteCarlo(pw.floating_scs_,1000).apply(pw.pose);
	// 	for(Size i = 1; i <= pw.floating_scs_.size();++i) {
	// 		Pose tmp = pose;
	// 		change_floating_sc_geometry(pose,pw.floating_scs_[i],pw.nres);
	// 		if((*sf)(tmp) < (*sf)(pose)) {
	// 			pose = tmp;
	// 		} else {
	// 			TR << "ScMinMover: performed sc swap on res " << i << std::endl;
	// 		}
	// 	}
	// }

	m.apply(pose);
	// TR << "minimize after  " << (*sf)(pose) << std::endl;
	pw.check_scattach_res();

}

// void add_sheet_cst(PoseWrap & pw) {
// 	Size sheet_start=999999, sheet_end=0;
// 	for(size_t i = 0; i < pw.nres); ++i) {
// 		if(pw.ss(i)=='E') {
// 			sheet_start = numeric::min(sheet_start,i+1);
// 			sheet_end = numeric::max(sheet_end,i+1);
// 		}
// 	}
// 	if((sheet_end-sheet_start)%2==1){
// 		TR << "even sheet length!!!" << std::endl;
// 		std::exit(-1);
// 	}
// 	using core::id::AtomID;
// 	for(size_t i = sheet_start; i <= sheet_end; ++i) {
// 		Size other = sheet_end+sheet_start-i;
// 		add_apc(pw.pose,AtomID(2,i),AtomID(2,pw.nres*3+other),5.0,1.0);
// 		TR << "apc " << i << " " << pw.nres*3+other << std::endl;
// 	}
//
// }

core::scoring::ScoreFunctionOP
cen_fold(PoseWrap & pw, Size NCYCLES, core::fragment::FragSetOP frags0, core::fragment::FragSetOP frags3 ) {
	using namespace core;
	using namespace scoring;
	using namespace protocols::moves;
	if(pw.pose.is_fullatom()) pw.switch_to_cen();
	core::pose::Pose & pose(pw.pose);

	Real temp   = basic::options::option[basic::options::OptionKeys::smhybrid::temperature]();
	Real rb_mag = basic::options::option[basic::options::OptionKeys::smhybrid::rb_mag]();
	Size csub   = basic::options::option[basic::options::OptionKeys::smhybrid::switch_concert_sub]();
	MoverOP slide = new MySlideMover(pw.pose);
	slide->apply(pw.pose);
	// Size cres = 0; if(csub) cres = pw.attach_rsd_[pw.primary_subsub];
	MoverOP minsc,minscpicky,changesc;
	if(pw.has_floating_sc) {
		minsc = new ScMinMover(pw.pose,pw.floating_scs_,0.05,pw.nres);//make_float_sc_min_mover(pw.pose,pw.floating_scs_[1]);
		minscpicky = new ScMinMover(pw.pose,pw.floating_scs_,0.0001,pw.nres);//make_float_sc_min_mover(pw.pose,pw.floating_scs_[1]);
		changesc = new FloatScRotMover(pw);
	}
	RandomMoverOP mover0 = new RandomMover;
	if( frags0() != NULL && frags0->size() > 0 ) {
		protocols::moves::MoverOP fragins = new protocols::abinitio::ClassicFragmentMover(frags0);
		// TODO FIX HACK
		TR << "doing frag moves" << std::endl;
		mover0->add_mover(fragins,5.0);
	} else {
		TR << "doing random bb moves instead of frags0!!" << std::endl;
		protocols::moves::MoverOP bbmove = new BBMover(pw.nres,30.0);
		mover0->add_mover(bbmove,5.0);
	}
	if( core::conformation::symmetry::symmetry_info(pw.pose)->get_dofs().size() ) {
		// mover0->add_mover(slide,0.005);
		MoverOP rot      = new MyRotMover  (pw,rb_mag*180.0,false,false,minsc,0.15,csub);
		MoverOP rotc     = new MyRotMover  (pw,rb_mag*180.0,false,true ,minsc,0.15,csub);
		MoverOP rotsm    = new MyRotMover  (pw,rb_mag* 20.0,false,false,minsc,0.05,csub);
		MoverOP rotsmc   = new MyRotMover  (pw,rb_mag* 20.0,false,true ,minsc,0.05,csub);
		MoverOP trans    = new MyTransMover(pw,rb_mag*  3.0,      false,minsc,0.15,csub);
		MoverOP transc   = new MyTransMover(pw,rb_mag*  3.0,      true ,minsc,0.15,csub);
		MoverOP transsm  = new MyTransMover(pw,rb_mag*  1.0,      false,minsc,0.05,csub);
		MoverOP transsmc = new MyTransMover(pw,rb_mag*  1.0,      true ,minsc,0.05,csub);
		// if( option[basic::options::OptionKeys::smhybrid::refine]() ) {
		// 	MoverOP rot      = new MyRotMover  (pw,rb_mag*2.0,false,false,minsc,0.15,csub);
		// 	MoverOP rotc     = new MyRotMover  (pw,rb_mag*1.0,false,true ,minsc,0.15,csub);
		// 	MoverOP rotsm    = new MyRotMover  (pw,rb_mag*1.0,false,false,minsc,0.15,csub);
		// 	MoverOP rotsmc   = new MyRotMover  (pw,rb_mag*1.0,false,true ,minsc,0.15,csub);
		// 	MoverOP trans    = new MyTransMover(pw,rb_mag*1.0,      false,minsc,0.15,csub);
		// 	MoverOP transc   = new MyTransMover(pw,rb_mag*1.0,      true ,minsc,0.15,csub);
		// 	MoverOP transsm  = new MyTransMover(pw,rb_mag*0.1,      false,minsc,0.15,csub);
		// 	MoverOP transsmc = new MyTransMover(pw,rb_mag*0.1,      true ,minsc,0.15,csub);
		// }
		if(true) mover0->add_mover(rot     ,0.3);
		if(true) mover0->add_mover(rotsm   ,3.0);
		if(true) mover0->add_mover(trans   ,0.1);
		if(true) mover0->add_mover(transsm ,1.0);
		if(csub) mover0->add_mover(rotc    ,0.1);
		if(csub) mover0->add_mover(rotsmc  ,1.0);
		if(csub) mover0->add_mover(transc  ,0.1);
		if(csub) mover0->add_mover(transsmc,1.0);
		/////////////// if(changesc) mover0->add_mover(changesc,0.8);
		/////////////// if(minsc) mover0->add_mover(minscpicky,0.2);
		// std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		// for(Size i = 1; i <= 40; i++) transsm->apply(pose);
	}
	TR << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	// for(Size i = 1; i <= pw.attach_rsd_.size(); ++i) {
	// 	if( pw.move_chi( pw.attach_rsd_[i] ) ) {
	// 		TR << "add chi mover res " << pw.attach_rsd_[i] << std::endl;
	// 		mover0->add_mover(new ChiMover(pw.pose,pw.attach_rsd_[i]),1.0);
	// 	}
	// }



	scoring::ScoreFunctionOP sf0 = scoring::ScoreFunctionFactory::create_score_function( "score0" );
	scoring::ScoreFunctionOP sf1 = scoring::ScoreFunctionFactory::create_score_function( "score1" );
	scoring::ScoreFunctionOP sf2 = scoring::ScoreFunctionFactory::create_score_function( "score2" );
	scoring::ScoreFunctionOP sf3 = scoring::ScoreFunctionFactory::create_score_function( "score3" );
	scoring::ScoreFunctionOP sf5 = scoring::ScoreFunctionFactory::create_score_function( "score5" );
	sf0 = new scoring::symmetry::SymmetricScoreFunction(*sf0);
	sf1 = new scoring::symmetry::SymmetricScoreFunction(*sf1);
	sf2 = new scoring::symmetry::SymmetricScoreFunction(*sf2);
	sf3 = new scoring::symmetry::SymmetricScoreFunction(*sf3);
	sf5 = new scoring::symmetry::SymmetricScoreFunction(*sf5);

	// sf0->set_weight(scoring::pair   ,0.5);
	// sf0->set_weight(scoring::hs_pair,0.5);
	// sf0->set_weight(scoring::ss_pair,0.5);
	// sf0->set_weight(scoring::rsigma ,0.5);
	// sf0->set_weight(scoring::sheet  ,0.5);
	// sf1->set_weight(scoring::pair   ,1.0);
	// sf1->set_weight(scoring::hs_pair,1.0);
	// sf1->set_weight(scoring::ss_pair,1.0);
	// sf1->set_weight(scoring::rsigma ,1.0);
	// sf1->set_weight(scoring::sheet  ,1.0);
	// sf2->set_weight(scoring::pair   ,1.0);
	// sf2->set_weight(scoring::hs_pair,2.0);
	// sf2->set_weight(scoring::ss_pair,2.0);
	// sf2->set_weight(scoring::rsigma ,2.0);
	// sf2->set_weight(scoring::sheet  ,2.0);
	// sf3->set_weight(scoring::pair   ,2.0);
	// sf3->set_weight(scoring::hs_pair,2.0);
	// sf3->set_weight(scoring::ss_pair,2.0);
	// sf3->set_weight(scoring::rsigma ,2.0);
	// sf3->set_weight(scoring::sheet  ,2.0);
	// sf5->set_weight(scoring::pair   ,1.0);
	// sf5->set_weight(scoring::hs_pair,1.0);
	// sf5->set_weight(scoring::ss_pair,1.0);
	// sf5->set_weight(scoring::rsigma ,1.0);
	// sf5->set_weight(scoring::sheet  ,1.0);
	//
	// sf0->set_weight(scoring::vdw,10.0);
	// sf1->set_weight(scoring::vdw,10.0);
	// sf2->set_weight(scoring::vdw,10.0);
	// sf3->set_weight(scoring::vdw,10.0);
	// sf5->set_weight(scoring::vdw,10.0);
	// sf0->set_weight(scoring::linear_chainbreak, 2.0);
	// sf1->set_weight(scoring::linear_chainbreak,10.0);
	// sf2->set_weight(scoring::chainbreak,10.0);
	// sf3->set_weight(scoring::chainbreak,10.0);
	// sf5->set_weight(scoring::chainbreak, 5.0);
	sf0->set_weight(scoring::rg,5.5);
	sf1->set_weight(scoring::rg,5.5);
	sf2->set_weight(scoring::rg,2.0);
	sf3->set_weight(scoring::rg,1.0);
	sf5->set_weight(scoring::rg,5.5);
	sf0->set_weight(scoring::coordinate_constraint,0.5);
	sf1->set_weight(scoring::coordinate_constraint,1);
	sf2->set_weight(scoring::coordinate_constraint,5);
	sf3->set_weight(scoring::coordinate_constraint,5);
	sf5->set_weight(scoring::coordinate_constraint,1);
	// if(pw.has_floating_sc) {
	sf0->set_weight(scoring::atom_pair_constraint,0.001);
	sf1->set_weight(scoring::atom_pair_constraint,0.01);
	sf2->set_weight(scoring::atom_pair_constraint,0.05);
	sf3->set_weight(scoring::atom_pair_constraint,1.0);
	sf5->set_weight(scoring::atom_pair_constraint,0.5);

	// sf0->set_weight(scoring::spline,0.02);
	// sf1->set_weight(scoring::spline,0.03);
	// sf2->set_weight(scoring::spline,0.04);
	// sf3->set_weight(scoring::spline,0.05);
	// sf5->set_weight(scoring::spline,0.02);
		// sf0->set_weight(scoring::floating_sc,0.05);
		// sf1->set_weight(scoring::floating_sc,0.05);
		// sf2->set_weight(scoring::floating_sc,0.05);
		// sf3->set_weight(scoring::floating_sc,0.05);
		// sf5->set_weight(scoring::floating_sc,0.05);
	// } else {
	// 	; //TR << "not using floating sc" << std::endl;
	// }
	// sf0->set_weight(scoring::coordinate_constraint,1.0);
	// sf1->set_weight(scoring::coordinate_constraint,1.0);
	// sf3->set_weight(scoring::coordinate_constraint,1.0);
	// sf5->set_weight(scoring::coordinate_constraint,1.0);


	if( option[basic::options::OptionKeys::edensity::mapfile].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *sf0 );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *sf1 );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *sf2 );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *sf3 );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *sf5 );
	}

	// add_sheet_cst(pw);
	// addcc(pw.pose,core::id::AtomID(2,5),core::id::AtomID(1,pw.nres*pw.nsub+1),0.1);


	core::pose::Pose init = pose;
	if(pw.debug) pw.dump_pdb("cen_init.pdb");

	bool skip0=false,skip1=false,skip2=false,filter=options::option[basic::options::OptionKeys::smhybrid::filter].user();
	Real filter0=0,filter1=0,filter2=0,filter3=0;
	Size ntries=1;
	if(filter) {
		assert(options::option[basic::options::OptionKeys::smhybrid::filter]().size()==5);
		ntries  = options::option[basic::options::OptionKeys::smhybrid::filter]()[1];
		filter0 = options::option[basic::options::OptionKeys::smhybrid::filter]()[2];
		filter1 = options::option[basic::options::OptionKeys::smhybrid::filter]()[3];
		filter2 = options::option[basic::options::OptionKeys::smhybrid::filter]()[4];
		filter3 = options::option[basic::options::OptionKeys::smhybrid::filter]()[5];
	}
	for(Size TRIES = 1; TRIES <= ntries; TRIES++ ) {

		if(!skip0) {
			pose = init;

			TR << "stage 0" << std::endl;
			MonteCarloOP mc = new MonteCarlo( pose, *sf0, temp );
			mc->set_autotemp( true, temp );
			mc->set_temperature( temp );


			RepeatMover( new TrialMover( mover0, mc ), NCYCLES/5 ).apply( pose );
			mc->reset( pose );
			cen_design(pw,sf0);
			if(pw.debug) sf0->show(TR,pose);
			if(pw.debug) pw.dump_pdb("cen0.pdb");
			if(filter) {
				if( pw.pose.energies().total_energies()[scoring::linear_chainbreak] > 50.0 ||
				 	pw.pose.energies().total_energies()[scoring::atom_pair_constraint] > filter0 ) {
					TR << "linear_chainbreak or atom_pair_constraint FAIL after score0! redoing centroid 0: "
					   << pw.pose.energies().total_energies()[scoring::linear_chainbreak] << " "
					   << pw.pose.energies().total_energies()[scoring::atom_pair_constraint] << std::endl;
					if( pw.pose.energies().total_energies()[scoring::linear_chainbreak]       > 50.0 ) {
						sf0->set_weight(scoring::linear_chainbreak, sf0->get_weight(scoring::linear_chainbreak)*1.2 );
						TR << "upweight linear_chainbreak" << std::endl;
					}
					if( pw.pose.energies().total_energies()[scoring::atom_pair_constraint]       > filter0 ) {
						sf0->set_weight(scoring::atom_pair_constraint, sf0->get_weight(scoring::atom_pair_constraint)*1.2 );
						TR << "upweight atom_pair_constraint" << std::endl;
					}
					skip0 = false;
					continue;
				} else {
					init = pose;
					skip0 = true;
				}
			}
		}

		if(!skip1) {

			TR << "stage 1" << std::endl;
			MonteCarloOP mc = new MonteCarlo( pose, *sf1, 2.0 );
			mc->set_autotemp( true, temp );
			mc->set_temperature( temp );
			RepeatMover( new TrialMover( mover0, mc ), NCYCLES ).apply( pose );
			mc->reset( pose );
			cen_design(pw,sf1);
			if(pw.debug) pw.dump_pdb("cen1.pdb");

			if(filter) {
				if( pw.pose.energies().total_energies()[scoring::linear_chainbreak] > 30.0 ||
				 	pw.pose.energies().total_energies()[scoring::atom_pair_constraint]       > filter1 ) {
					TR << "linear_chainbreak or atom_pair_constraint FAIL after score1! redoing centroid 1 "
					   << pw.pose.energies().total_energies()[scoring::linear_chainbreak] << " "
					   << pw.pose.energies().total_energies()[scoring::atom_pair_constraint] << std::endl;
					if( pw.pose.energies().total_energies()[scoring::linear_chainbreak]       > 30.0 ) {
						sf1->set_weight(scoring::linear_chainbreak, sf1->get_weight(scoring::linear_chainbreak)*1.2 );
						TR << "upweight linear_chainbreak" << std::endl;
					}
					if( pw.pose.energies().total_energies()[scoring::atom_pair_constraint]       > filter1 ) {
						sf1->set_weight(scoring::atom_pair_constraint, sf1->get_weight(scoring::atom_pair_constraint)*1.2 );
						TR << "upweight atom_pair_constraint" << std::endl;
					}
					skip1 = false;
					continue;
				} else {
					skip1 = true;
					// init = pose;
				}
			}
		}

		if(pw.debug) sf1->show(TR,pose);
		// pose.dump_pdb("stage1.pdb");

		// NO SHEETS AT THE MOMENT....
		if(!skip2) {
			for( Size i = 1; i <= 5; ++i ) {
				TR << "stage 2 " << i << std::endl;
				MonteCarloOP mc = new MonteCarlo( pose, *sf2, 2.0 );
				mc->set_autotemp( true, temp );
				mc->set_temperature( temp );
				RepeatMover( new TrialMover( mover0, mc ), NCYCLES/3 ).apply( pose );
				mc->reset( pose );
				cen_design(pw,sf2);
				// if(i>3) design(pw,sf2);
				if(pw.debug) sf2->show(TR,pose);

				mc = new MonteCarlo( pose, *sf5, 2.0 );
				mc->set_autotemp( true, temp );
				mc->set_temperature( temp );
				RepeatMover( new TrialMover( mover0, mc ), NCYCLES/3 ).apply( pose );
				mc->reset( pose );
				cen_design(pw,sf5);
				// if(i>3) design(pw,sf5);
				if(pw.debug) sf5->show(TR,pose);
				if(pw.debug) pw.dump_pdb("cen25.pdb");
			}
			if(filter) {
				if( pw.pose.energies().total_energies()[scoring::chainbreak] > 10.0 ||
				    pw.pose.energies().total_energies()[scoring::atom_pair_constraint]       > filter2 ) {
					TR << "chainbreak or atom_pair_constraint FAIL after score25! redoing centroid 2 "
					   << pw.pose.energies().total_energies()[scoring::chainbreak] << " "
					   << pw.pose.energies().total_energies()[scoring::atom_pair_constraint] << std::endl;
					if( pw.pose.energies().total_energies()[scoring::chainbreak]       > 10.0 ) {
						sf2->set_weight(scoring::chainbreak, sf2->get_weight(scoring::chainbreak)*1.2 );
						sf5->set_weight(scoring::chainbreak, sf5->get_weight(scoring::chainbreak)*1.2 );
						TR << "upweight chainbreak" << std::endl;
					}
					if( pw.pose.energies().total_energies()[scoring::atom_pair_constraint]       > filter2 ) {
						sf2->set_weight(scoring::atom_pair_constraint, sf2->get_weight(scoring::atom_pair_constraint)*1.2 );
						sf5->set_weight(scoring::atom_pair_constraint, sf2->get_weight(scoring::atom_pair_constraint)*1.2 );
						TR << "upweight atom_pair_constraint" << std::endl;
					}
					skip2 = false;
					continue;
				} else {
					skip2 = true;
					init = pose;
				}
			}
		}

		RandomMoverOP mover3 = new RandomMover;
/*
		if(frags0 != NULL && frags3 == NULL) {
			TR << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			for(Size i = 0; i < pw.ss_.size(); ++i) pw.ss_[i] = pw.pose.secstruct(i+1);
			core::fragment::FragSetOP frags3 = make_frag_set(pw,fds);
			protocols::moves::MoverOP fragins3 = new protocols::abinitio::ClassicFragmentMover(frags3);
		} else*/
 if( frags3() != NULL && frags3->size() > 0 ) {
			protocols::moves::MoverOP fragins = new protocols::abinitio::ClassicFragmentMover(frags3);
			mover3->add_mover(fragins,5.0);
		} else {
			TR << "doing random bb moves instead of frags3!!" << std::endl;
			protocols::moves::MoverOP bbmove = new BBMover(pw.nres,5.0);
			mover3->add_mover(bbmove,5.0);
		}
		if( core::conformation::symmetry::symmetry_info(pw.pose)->get_dofs().size() ) {
			// mover3->add_mover(slide,0.05);
			MoverOP rot      = new MyRotMover(  pw,rb_mag*10.0,false,false,minsc,0.5,csub);
			MoverOP rotc     = new MyRotMover(  pw,rb_mag*10.0,false,true ,minsc,0.5,csub);
			MoverOP rotsm    = new MyRotMover(  pw,rb_mag* 2.0,false,false,minsc,0.1,csub);
			MoverOP rotsmc   = new MyRotMover(  pw,rb_mag* 2.0,false,true ,minsc,0.1,csub);
			MoverOP trans    = new MyTransMover(pw,rb_mag* 1.0,false,      minsc,0.5,csub);
			MoverOP transc   = new MyTransMover(pw,rb_mag* 1.0,true ,      minsc,0.5,csub);
			MoverOP transsm  = new MyTransMover(pw,rb_mag* 0.1,false,      minsc,0.1,csub);
			MoverOP transsmc = new MyTransMover(pw,rb_mag* 0.1,true ,      minsc,0.1,csub);
			// if( option[basic::options::OptionKeys::smhybrid::refine]() ) {
			// 	MoverOP rot      = new MyRotMover(  pw,rb_mag*1.0,false,false,minsc,0.2,csub);
			// 	MoverOP rotc     = new MyRotMover(  pw,rb_mag*1.0,false,true ,minsc,0.2,csub);
			// 	MoverOP rotsm    = new MyRotMover(  pw,rb_mag*0.1,false,false,minsc,0.2,csub);
			// 	MoverOP rotsmc   = new MyRotMover(  pw,rb_mag*0.1,false,true ,minsc,0.2,csub);
			// 	MoverOP trans    = new MyTransMover(pw,rb_mag*0.5,      false,minsc,0.2,csub);
			// 	MoverOP transc   = new MyTransMover(pw,rb_mag*0.5,      true ,minsc,0.2,csub);
			// 	MoverOP transsm  = new MyTransMover(pw,rb_mag*0.1,      false,minsc,0.2,csub);
			// 	MoverOP transsmc = new MyTransMover(pw,rb_mag*0.1,      true ,minsc,0.2,csub);
			// }
			if(true) mover3->add_mover(rot     ,0.2);
			if(true) mover3->add_mover(rotsm   ,1.0);
			if(true) mover3->add_mover(trans   ,0.1);
			if(true) mover3->add_mover(transsm ,1.0);
			if(csub) mover3->add_mover(rotc    ,1.2);
			if(csub) mover3->add_mover(rotsmc  ,4.0);
			if(csub) mover3->add_mover(transc  ,0.5);
			if(csub) mover3->add_mover(transsmc,2.5);
			//////////////////////////if(minsc) mover3->add_mover(minscpicky,0.2);
			//////////////////////////if(changesc) mover0->add_mover(changesc,0.8);

		}
		TR << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		// for(Size i = 1; i <= pw.attach_rsd_.size(); ++i) {
		// 	if( pw.move_chi( pw.attach_rsd_[i] ) ) {
		// 		TR << "add chi mover res " << pw.attach_rsd_[i] << std::endl;
		// 		mover3->add_mover(new ChiMover(pw.pose,pw.attach_rsd_[i]),1.0);
		// 	}
		// }

		for( Size i = 1; i <= 10; ++i ) {
			TR << "stage 3" << i << std::endl;
			MonteCarloOP mc = new MonteCarlo( pose, *sf3, 2.0 );
			mc->set_autotemp( true, temp );
			mc->set_temperature( temp );
			RepeatMover( new TrialMover( mover3, mc ), NCYCLES/2 ).apply( pose );
			mc->reset( pose );
			cen_design(pw,sf3);
			if(pw.debug) sf3->show(TR,pose);
			if(pw.debug) pw.dump_pdb("cen3.pdb");
		}
		TR << std::endl;

		if(filter) {
			if( pw.pose.energies().total_energies()[scoring::chainbreak] > 5.0 ||
			 	pw.pose.energies().total_energies()[scoring::atom_pair_constraint]       > filter3 ) {
				TR << "chainbreak or atom_pair_constraint FAIL after score3! redoing centroid 3 "
				   << pw.pose.energies().total_energies()[scoring::chainbreak] << " "
				   << pw.pose.energies().total_energies()[scoring::atom_pair_constraint] << std::endl;
				if( pw.pose.energies().total_energies()[scoring::chainbreak]       > 5.0 ) {
					sf3->set_weight(scoring::chainbreak, sf3->get_weight(scoring::chainbreak)*1.2 );
					TR << "upweight chainbreak" << std::endl;
				}
				if( pw.pose.energies().total_energies()[scoring::atom_pair_constraint]       > filter3 ) {
					sf3->set_weight(scoring::atom_pair_constraint, sf3->get_weight(scoring::atom_pair_constraint)*1.2 );
					TR << "upweight atom_pair_constraint" << std::endl;
				}
				continue;
			}
		}

		break;
	}

	sf3->show(TR,pose);
	if(pw.debug) pw.dump_pdb("cen.pdb");

	return sf3;
}

core::scoring::ScoreFunctionOP
fa_refine_and_design(PoseWrap & pw, Size NCYCLE) {
	using namespace core;
	using namespace scoring;
	using namespace ObjexxFCL::format;
	core::pose::Pose & pose(pw.pose);
	ScoreFunctionOP sf1,sf2,sf3,sf4;
	sf1 = core::scoring::get_score_function();
	sf2 = core::scoring::get_score_function();
	sf3 = core::scoring::get_score_function();
	sf4 = core::scoring::get_score_function();

	sf1->set_weight(fa_rep,0.03);
	sf2->set_weight(fa_rep,0.11);
	sf3->set_weight(fa_rep,0.22);
	sf4->set_weight(fa_rep,0.44);

	sf1->set_weight(fa_atr,1.000);
	sf2->set_weight(fa_atr,1.000);
	sf3->set_weight(fa_atr,1.000);
	sf4->set_weight(fa_atr,1.000);

	sf1->set_weight(fa_dun,0.100);
	sf2->set_weight(fa_dun,0.150);
	sf3->set_weight(fa_dun,0.300);
	sf4->set_weight(fa_dun,0.600);

	sf1->set_weight(omega,0.2);
	sf2->set_weight(omega,0.4);
	sf3->set_weight(omega,0.6);
	sf4->set_weight(omega,1.0);

	sf1->set_weight(fa_sol,0.800);
	sf2->set_weight(fa_sol,0.800);
	sf3->set_weight(fa_sol,0.800);
	sf4->set_weight(fa_sol,0.800);

	// sf1->set_weight(overlap_chainbreak,1.00);
	// sf2->set_weight(overlap_chainbreak,2.00);
	// sf3->set_weight(overlap_chainbreak,3.00);
	// sf4->set_weight(overlap_chainbreak,4.00);

// std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	sf1->set_weight(chainbreak,1.00);
	sf2->set_weight(chainbreak,2.00);
	sf3->set_weight(chainbreak,3.00);
	sf4->set_weight(chainbreak,4.00);

	sf1->set_weight(coordinate_constraint,0.4);
	sf2->set_weight(coordinate_constraint,0.6);
	sf3->set_weight(coordinate_constraint,0.8);
	sf4->set_weight(coordinate_constraint,1.0);

	sf1->set_weight(atom_pair_constraint,8.0);
	sf2->set_weight(atom_pair_constraint,5.0);
	sf3->set_weight(atom_pair_constraint,3.0);
	sf4->set_weight(atom_pair_constraint,2.0);
	if(option[basic::options::OptionKeys::smhybrid::refine].user()) {
		sf1->set_weight(atom_pair_constraint,80.0);
		sf2->set_weight(atom_pair_constraint,50.0);
		sf3->set_weight(atom_pair_constraint,30.0);
		sf4->set_weight(atom_pair_constraint,20.0);
	}

	// sf1->set_weight(rg,1.00);
	// sf2->set_weight(rg,1.00);
	// sf3->set_weight(rg,1.00);
	// sf4->set_weight(rg,1.00);

	// sf1->set_weight(scoring::spline,1.0);
	// sf2->set_weight(scoring::spline,0.5);
	// sf3->set_weight(scoring::spline,0.1);
	// sf4->set_weight(scoring::spline,0.0);

	// sf1->set_weight(scoring::sym_lig,10.0);
	// sf2->set_weight(scoring::sym_lig,10.0);
	// sf3->set_weight(scoring::sym_lig,10.0);
	// sf4->set_weight(scoring::sym_lig,10.0);
	sf1 = new symmetry::SymmetricScoreFunction(*sf1);
	sf2 = new symmetry::SymmetricScoreFunction(*sf2);
	sf3 = new symmetry::SymmetricScoreFunction(*sf3);
	sf4 = new symmetry::SymmetricScoreFunction(*sf4);
	// core::pose::Pose best = pose;

	string tmp = string_of(uniform());
	if(pw.debug) pw.dump_pdb("fa_init"+tmp+".pdb");
	// pw.add_fsc_harmonic_cst();
	sf1->set_weight(atom_pair_constraint,100);
	pw.check_scattach_res();
	repack(pw,sf1,true);
	pw.check_scattach_res();
	// if(pw.debug) pw.dump_pdb("fa0PACK"+tmp+".pdb");
	minimize(pw,sf1,0,false);
	pw.check_scattach_res();
	if(pw.debug) pw.dump_pdb("fa0_"+tmp+".pdb");
	sf1->set_weight(atom_pair_constraint,1);
	// std::exit(-1);
	pw.check_scattach_res();

	int bb = basic::options::option[basic::options::OptionKeys::smhybrid::minbb]();
	int chi1 = -1;
	for(Size i = 1; i <= NCYCLE; ++i ) {
		// for(Size k = 1; k <= 4; ++k) {
		// 	repack(pw,sf1); minimize(pw,sf1);
		// 	repack(pw,sf2); minimize(pw,sf2);
		// 	repack(pw,sf3); minimize(pw,sf3);
		// 	repack(pw,sf4); minimize(pw,sf4);
		// 	TR << "repack/min1234    " << I(2,k) << " " << F(10,3,(*sf4)(pose)) << " " << pose.sequence().substr(0,nres_mono) << std::endl;
		// 	// if( (*sf4)(best) >= (*sf4)(pose) ) best = pose;
		// }
		// // pose = best;
		// for(Size k = 1; k <= 1; ++k) {
			TR << "bb = " << bb << std::endl;
			// std::cerr << "SCORE 1" << std::endl;
			design(pw,sf1); if(!pw.check()) return NULL;
			if(pw.debug) pw.dump_pdb("fa1A_"+tmp+".pdb");
			minimize(pw,sf1,bb-1,false,chi1>0);
			if(pw.debug) sf1->show(std::cerr,pose);
			if(pw.debug) pw.dump_pdb("fa1B_"+tmp+".pdb");
			// minimize(pw,sf1,bb-1,true);
			// if(pw.debug) sf1->show(std::cerr,pose);
			// if(pw.debug) pw.dump_pdb("fa1_"+tmp+".pdb");

			// std::cerr << "SCORE 2" << std::endl;
			design(pw,sf2); if(!pw.check()) return NULL;
			if(pw.debug) pw.dump_pdb("fa2A_"+tmp+".pdb");
			minimize(pw,sf2,bb-1,false,chi1>0);
			if(pw.debug) sf2->show(std::cerr,pose);
			if(pw.debug) pw.dump_pdb("fa2B_"+tmp+".pdb");
			// minimize(pw,sf2,bb-1,true);
			// if(pw.debug) sf2->show(std::cerr,pose);
			// if(pw.debug) pw.dump_pdb("fa2_"+tmp+".pdb");

			// std::cerr << "SCORE 3" << std::endl;
			design(pw,sf3);  if(!pw.check()) return NULL;
			if(pw.debug) pw.dump_pdb("fa3A_"+tmp+".pdb");
			minimize(pw,sf3,bb-1,false,chi1>0);
			if(pw.debug) sf3->show(std::cerr,pose);
			if(pw.debug) pw.dump_pdb("fa3B_"+tmp+".pdb");
			// minimize(pw,sf3,bb-1,true);
			// if(pw.debug) sf3->show(std::cerr,pose);
			// if(pw.debug) pw.dump_pdb("fa3_"+tmp+".pdb");

			// std::cerr << "SCORE 4" << std::endl;
			design(pw,sf4);  if(!pw.check()) return NULL;
			if(pw.debug) pw.dump_pdb("fa4A_"+tmp+".pdb");
			minimize(pw,sf4,bb-1,false,chi1>0);
			if(pw.debug) sf4->show(std::cerr,pose);
			if(pw.debug) pw.dump_pdb("fa4B_"+tmp+".pdb");
			// minimize(pw,sf4,bb-1,true);
			// if(pw.debug) sf4->show(std::cerr,pose);
			// if(pw.debug) pw.dump_pdb("fa4_"+tmp+".pdb");
			// TR << "deisgn/min1234   " << I(2,k) << " " << F(10,3,(*sf4)(pose)) << " " << pose.sequence().substr(0,pw.nres) << std::endl;
			// if( (*sf4)(best) >= (*sf4)(pose) ) best = pose;
		// }
		// pose = best;

		// sf4->set_weight(overlap_chainbreak,3.00);
		// sf4->set_weight(chainbreak,3.0);
			// std::cerr << "BB MIN" << std::endl;
			// sf4->set_weight(atom_pair_constraint,1.00);
			design(pw,sf4,true);  if(!pw.check()) return NULL;
			if(pw.debug) pw.dump_pdb("fa5A_"+tmp+".pdb");
			minimize(pw,sf4,bb,false,chi1>0);
			if(pw.debug) pw.dump_pdb("fa5_"+tmp+".pdb");
		// if(pw.debug) sf4->show(TR,pose);
		// TR << "deisgn/min cb3   " << F(10,3,(*sf4)(pose)) << " " << pose.sequence().substr(0,pw.nres) << std::endl;

		// sf4->set_weight(overlap_chainbreak,2.00);
		// sf4->set_weight(chainbreak,2.0);
		// sf4->set_weight(atom_pair_constraint,0.50);
		// design(pw,sf4); minimize(pw,sf4,true);
		// if(pw.debug) sf4->show(TR,pose);
		// TR << "deisgn/min cb3   " << F(10,3,(*sf4)(pose)) << " " << pose.sequence().substr(0,pw.nres) << std::endl;
		//
		// sf4->set_weight(overlap_chainbreak,1.00);
		// sf4->set_weight(chainbreak,1.0);
		// sf4->set_weight(atom_pair_constraint,0.20);
		// design(pw,sf4); minimize(pw,sf4,true);
		// if(pw.debug) sf4->show(TR,pose);
		// TR << "deisgn/min cb3   " << F(10,3,(*sf4)(pose)) << " " << pose.sequence().substr(0,pw.nres) << std::endl;

		// sf3->set_weight(coordinate_constraint,1.0);
		// sf4->set_weight(coordinate_constraint,1.0);
		// for(Size jj=1;jj<=100;++jj) TR << " " << std::endl;
		// design(pw,sf3); minimize(pw,sf3);
		// design(pw,sf4); minimize(pw,sf4);

		sf4->show(TR,pose);
		TR << std::endl;

		bb++;
		chi1++;
	}
	return sf4;
}




ScoreFunctionOP flxbb_nobu(PoseWrap & pw) {
	using namespace core::scoring;
	ScoreFunctionOP sf3,sf4;
	sf3 = get_score_function();
	sf4 = get_score_function();
	sf3->set_weight(fa_rep,0.100);
	sf4->set_weight(fa_rep,0.400);
	sf3 = new symmetry::SymmetricScoreFunction(*sf3);
	sf4 = new symmetry::SymmetricScoreFunction(*sf4);

	if( basic::options::option[basic::options::OptionKeys::smhybrid::fa_refine]() ) {
		protocols::flxbb::FlxbbDesign flxbb(sf3,sf4);
		flxbb.apply(pw.pose);
	} // else {
	 // 		protocols::relax::SimpleMultiRelax smr(sf3);
	 // 		smr.apply(pw.pose);
	 // 	}

	return sf4;
}

vector1<Real>
interface_energy_ratio( core::pose::Pose const & pose, ScoreFunctionOP sf ) {
	using namespace core;
	using namespace conformation::symmetry;
	using namespace scoring;
	SymmetricConformation const & SymmConf (dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	Energies const & energies( pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );

	vector1<Real> ife(symm_info->num_interfaces()+1,0.0);
	Real tot = 0.0;
	for ( Size i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
		// conformation::Residue const & resl( pose.residue( i ) );
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
			Size const j( edge.get_second_node_ind() );
			// conformation::Residue const & resu( pose.residue( j ) );
			Size const interface_number  = symm_info->interface_number(i,j);
			Size const score_mult = symm_info->score_multiply(i,j);
			Real const score = edge.dot(sf->weights());
			ife[interface_number] += score*score_mult;
			tot += score*score_mult;
		}
	}

	for(size_t i = 1; i <= ife.size(); ++i) {
		ife[i] /= tot;
	}

	return ife;
}

Real
symm_self_rep( PoseWrap & pw, ScoreFunctionOP sf, Size rsd ) {
	if(!sf)return 0;
	using namespace core;
	using namespace conformation::symmetry;
	using namespace scoring;
	SymmetricConformation const & SymmConf (dynamic_cast<SymmetricConformation const &> ( pw.pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	Energies const & energies( pw.pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );

	Real tot = 0.0;
	for ( Size i=1, i_end = pw.pose.total_residue(); i<= i_end; ++i ) {
		if( (i-1)%pw.nres+1 != rsd ) continue;
		// conformation::Residue const & resl( pose.residue( i ) );
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
			Size const j( edge.get_second_node_ind() );
			if( (j-1)%pw.nres+1 != rsd ) continue;
			// conformation::Residue const & resu( pose.residue( j ) );
		  	tot += symm_info->score_multiply(i,j) * edge[fa_rep];
		}
	}

	return tot;
}

std::map< std::pair<Size,Size>, Real >
find_clashes( PoseWrap & pw ) {
	using namespace core;
	using namespace conformation::symmetry;
	using namespace scoring;
	core::pose::Pose & pose(pw.pose);

	SymmetricConformation const & SymmConf (dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	Energies const & energies( pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );

	std::map< std::pair<Size,Size>, Real > clashes;
	for ( Size i=1, i_end = pw.nres; i<= i_end; ++i ) {
		// conformation::Residue const & resl( pose.residue( i ) );
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
			Size const j( edge.get_second_node_ind() );
			if(pose.residue(i).name3()=="BPY"&&pose.residue(j).name3()=="BPY") continue;
			// conformation::Residue const & resu( pose.residue( j ) );
			// if( edge[core::scoring::fa_rep] > 0 ) TR << "E EDGE rep " << i << " " << j << " " << edge[core::scoring::fa_rep] << std::endl;

			if(  edge[core::scoring::fa_rep] > 40.0 ) {
				clashes[std::pair<Size,Size>(i,j)] = edge[core::scoring::fa_rep];
			}
		}
	}

	return clashes;
}



void report( PoseWrap & pw, ScoreFunctionOP sf_fa, std::ostringstream & oss, Real censcore ) {
	using namespace core;
	using namespace ObjexxFCL::format;
	core::pose::Pose & pose(pw.pose);
	Size nres_mono = pw.nres;//(pose.n_residue()-4)/4;
	string tag = string_of(uniform());
	vector1<core::Real> ife;
	std::map< std::pair<Size,Size>, Real > clashes;
	if(sf_fa) {
	// core::pose::Pose posebefore = pose;
	// sf_fa->set_weight(scoring::holes_min,1.0);
		sf_fa->set_weight(scoring::rg,0.00001);
	// TR << "BEFORE HOLES MIN: " << (*sf4)(pose) << std::endl;
		(*sf_fa)(pose);
		sf_fa->show(TR,pose);
	// Real rholes = compute_rosettaholes_score(pw.pose).score();
	// core::Real before = pw.pose.energies().total_energies()[core::scoring::holes_min];
	// minimize(pw,sf_fa);
	// (*sf_fa)(pose);
	// TR << "holes_min " << before << " " << pw.pose.energies().total_energies()[core::scoring::holes_min] << std::endl;
	// sf4->show(pose);

	// for( Size i = 0; i < 4; ++i ) {
	// 	Size nres_mono = (pose.n_residue()-4)/4;
	// 	ref_rep += pose.energies().residue_total_energies(i*nres_mono+1)[scoring::fa_rep];
	//    	ref_rep += pose.energies().residue_total_energies(i*nres_mono+2)[scoring::fa_rep];
	// }


		ife = interface_energy_ratio(pose,sf_fa);
		clashes = find_clashes(pw);
	}
	using namespace core::scoring::constraints;
	ConstraintCOPs cs = pose.constraint_set()->get_all_constraints();
	// TR << "num cst: " << cs.size() << std::endl;
	for(Size i = 1; i <= cs.size(); i++) {
		// if( cs[i]->type() != "LocalCoordinateConstraint" ) continue;
		// LocalCoordinateConstraint const * cc = (LocalCoordinateConstraint const *)(cs[i]());
		// TR << "cst " << i << " " << cc->atom(1) << " " << cc->atom(2) << " " << cc->dist(pose) << std::endl;
	}

	core::conformation::symmetry::SymmetricConformation const & SymmConf (
		dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
	core::conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	Size nheavy = 0;
	for(Size i = 1; i <= pose.n_residue(); ++i) nheavy += pose.residue(i).nheavyatoms();

	Real sym_lig      = pose.energies().total_energies()[ scoring::sym_lig      ];
	Real floating_sc  = pose.energies().total_energies()[ scoring::atom_pair_constraint  ] / symm_info->score_multiply(1,1) / 6.0 / pw.floating_scs_.size();
	     floating_sc  += pose.energies().total_energies()[ scoring::angle_constraint   ] / symm_info->score_multiply(1,1) / 6.0 / pw.floating_scs_.size();
	Real lin_cb       = pose.energies().total_energies()[ scoring::linear_chainbreak  ] / 3.0 + pose.energies().total_energies()[ scoring::chainbreak  ] / 3.0;
	Real fa_atr       = pose.energies().total_energies()[ scoring::fa_atr       ];
	Real fa_rep       = pose.energies().total_energies()[ scoring::fa_rep       ];
	Real fa_sol       = pose.energies().total_energies()[ scoring::fa_sol       ];
	Real fa_intra_rep = pose.energies().total_energies()[ scoring::fa_intra_rep ];
	Real pro_close    = pose.energies().total_energies()[ scoring::pro_close    ];
	Real fa_pair      = pose.energies().total_energies()[ scoring::fa_pair      ];
	Real hbond_sr_bb  = pose.energies().total_energies()[ scoring::hbond_sr_bb  ];
	Real hbond_lr_bb  = pose.energies().total_energies()[ scoring::hbond_lr_bb  ];
	Real hbond_bb_sc  = pose.energies().total_energies()[ scoring::hbond_bb_sc  ];
	Real hbond_sc     = pose.energies().total_energies()[ scoring::hbond_sc     ];
	Real fa_dun       = pose.energies().total_energies()[ scoring::fa_dun       ];
	Real p_aa_pp      = pose.energies().total_energies()[ scoring::p_aa_pp      ];
	Real ref          = pose.energies().total_energies()[ scoring::ref          ];
	Real rg           = pose.energies().total_energies()[ scoring::rg           ];
	Real surface = 0, cavvol=0;
	if(sf_fa) surface = scoring::calc_total_sasa(pose,1.6);
	if(sf_fa) cavvol  = scoring::packstat::get_cavity_volume(pose);
	Real rholes       = 0;
	if(sf_fa) rholes = scoring::packstat::compute_packing_score( pw.pose , 2 );
	Real score = 0;
	if(sf_fa) score = (*sf_fa)(pose);
	if( sf_fa && pw.attach_rsd_[1]==1 && pw.attach_atom_[1]!="N" ) {
		Real symrep = symm_self_rep(pw,sf_fa,1);
		score -= sf_fa->get_weight(scoring::fa_rep) * symrep;
		fa_rep -= symrep;
		// std::cerr << "SYMREP: " << symrep << " " << score << " " << fa_rep << std::endl;
	}
	Real per_res_score = score/nres_mono/symm_info->score_multiply(1,1);
	oss << LJ(20,tag)          << " "
		<< pw.tag() << " "
	    << F(8,3,per_res_score)<< " "
	    << F(8,3,lin_cb    )<< " "
	    << F(8,3,floating_sc    )<< " "
		<< F(8,3,pw.rms_to_orig_subsubs() ) << " "
		<< F(8,3,pw.absrms_to_orig_subsubs() ) << " "
	    // << I(4, pw.get_closest_res_for_sc_attach(pw.floating_scs_[1]) )<< " "
        << F(8,3,ife[2]      ) << " "
        << F(8,3,ife[3]      ) << " "
        << F(8,3,ife[4]      ) << " "
		<< F(8,3, pose.energies().total_energies()[scoring::elec_dens_whole_structure_ca]/nres_mono/symm_info->score_multiply(1,1)) << " "
	    << F(8,3,score)        << " "
        << F(8,3,fa_atr      ) << " "
        << F(8,3,fa_rep      ) << " "
        << F(8,3,fa_sol      ) << " "
        << F(8,3,rholes      ) << " "
        << F(8,3,surface     ) << " "
        << F(8,3,cavvol      ) << " "
        << F(8,3,rg          ) << " "
		<< I(3,nres_mono)      << " "
		<< I(3,nheavy   )      << " "
	    << F(8,3,censcore)     << " "
		<< F(8,3,sym_lig     ) << " "
        << F(8,3,fa_intra_rep) << " "
        << F(8,3,pro_close   ) << " "
        << F(8,3,fa_pair     ) << " "
        << F(8,3,hbond_sr_bb ) << " "
        << F(8,3,hbond_lr_bb ) << " "
        << F(8,3,hbond_bb_sc ) << " "
        << F(8,3,hbond_sc    ) << " "
        << F(8,3,fa_dun      ) << " "
        << F(8,3,p_aa_pp     ) << " "
        << F(8,3,ref         ) << " "
		<< pose.sequence().substr(0,pw.nres);
	for(std::map< std::pair<Size,Size>, Real >::iterator i = clashes.begin(); i != clashes.end(); ++i) {
		oss << " CLASH " << i->first.first << "," << i->first.second << " " << i->second;
	}
	oss << std::endl;


	using namespace basic::options;
	string outdir = "./";
	if(option[OptionKeys::out::file::o].user()) outdir = option[OptionKeys::out::file::o]();
	if( per_res_score <= options::option[options::OptionKeys::in::file::silent_energy_cut]() ) {
		// cenpose.dump_pdb(string()+"/smhybrid"+tag+"_cen.pdb.gz");
		// posebefore .dump_pdb(outdir+"/smhybrid"+tag+"_fa.pdb.gz");
		{
			option[OptionKeys::out::file::output_virtual].value(true);
			string fn = outdir+"/smhybrid"+tag+"_virt.pdb";
			if(option[OptionKeys::out::pdb_gz]()) fn += ".gz";
			pw.dump_pdb(fn);
		}
		std::cout << oss.str();
		std::cout.flush();
		oss.clear();
		oss.str("");
	}

}




string make_rand_ss_hs(core::Size /*len*/) {
	string ss = "";
	Size nsheet = ((Size)(uniform()*3))*2+3;
	Size nhelix = uniform()*5+3;
	Size nloop = uniform()*4+0;
	for(Size i = 1; i <= nhelix; ++i) ss += "H";
	for(Size i = 1; i <= nloop+1; ++i) ss += "L";
	for(Size i = 1; i <= nhelix; ++i) ss += "H";
	for(Size i = 1; i <= nloop ; ++i) ss += "L";
	for(Size i = 1; i <= nsheet; ++i) ss += "E";
	// while(ss.size()<len) {
	// 	for(int i = 1; i <= uniform()*6   ; ++i) ss += "L";
	// 	for(int i = 1; i <= uniform()*15+5; ++i) ss += "H";
	// }
	return ss;
}


string make_rand_ss_h(core::Size len) {
	string ss = "";
	for(int i = 1; i <= uniform()*15+5; ++i) ss += "H";
	while(ss.size()<len) {
		for(int i = 1; i <= uniform()*6   ; ++i) ss += "L";
		for(int i = 1; i <= uniform()*15+5; ++i) ss += "H";
	}
	return ss;
}


// class SmhybridMover : public protocols::moves::Mover
// {
// public:
//
// public:
//
// 	// default constructor
// 	SmhybridMover() {}
//
// 	virtual ~SmhybridMover() {}
//
// 	// constructor with arguments
//
// 	virtual protocols::moves::MoverOP clone() const { return new SmhybridMover(*this); }
// 	virtual protocols::moves::MoverOP fresh_instance() const { return new SmhybridMover; }
//
// 	virtual void apply( core::pose::Pose & pose ){
void* doit(void* /*x = NULL*/) {
	using namespace core;
	using namespace pose;
	using namespace protocols;
	using namespace moves;
	using namespace ObjexxFCL::format;
	using numeric::random::uniform;

	PoseWrap pw = posewrap_from_command_line();


	std::map<string, vector1<core::fragment::FragDataOP> > fds = get_frags_map( true );
	std::ostringstream oss;
	Size Ncycle    = options::option[options::OptionKeys::smhybrid::abinitio_cycles]();

	for(Size iter = 1; iter <= (Size)options::option[options::OptionKeys::out::nstruct](); ++iter) {
		TR << "ITER " << iter << std::endl;


		PoseWrap pw = posewrap_from_command_line();

		// fold up
		// core::fragment::FragSetOP frags0 = make_frag_set_9mers(pw.nres);
		core::fragment::FragSetOP frags0 = NULL;//make_frag_set(pw,fds);
		core::fragment::FragSetOP frags3 = NULL;//make_frag_set(pw,fds);

		scoring::ScoreFunctionOP sf_cen = cen_fold(pw,Ncycle,frags0,frags3);

		if( options::option[basic::options::OptionKeys::smhybrid::filter].user() ) {
			if( pw.pose.energies().total_energies()[scoring::linear_chainbreak]    > 5.0 ||
		 		pw.pose.energies().total_energies()[scoring::atom_pair_constraint] > options::option[basic::options::OptionKeys::smhybrid::filter]()[5] )
			{
				TR << "chainbreak or floating_sc FAIL! redoing centroid" << std::endl;
				sf_cen->show(TR,pw.pose);
				--iter;
				continue;
			}
		}

		Real censcore = (*sf_cen)(pw.pose);
		// pw.pose.dump_pdb("cen.pdb");


		if(option[basic::options::OptionKeys::smhybrid::fa_refine]()) {
			ScoreFunctionOP sf_fa;
			// pw.dump_pdb("cen.pdb");
			if(!pw.switch_to_fa()) continue;
			sf_fa = fa_refine_and_design(pw,3);
			if(!sf_fa) {
				TR << "FA PW FAIL" << std::endl;
				continue;
			}
			report(pw,sf_fa,oss,censcore);
		} else {
			string tag = string_of(uniform());
			TR << "cen" << tag << " " << censcore << std::endl;
			if(censcore < option[basic::options::OptionKeys::in::file::silent_energy_cut]) {
				pw.pose.dump_pdb("cen"+tag+".pdb");
			}
		}

	}

	return NULL;
}


void* doit_refine(void* /*x = NULL*/) {
	using namespace core;
	using namespace pose;
	using namespace protocols;
	using namespace moves;
	using namespace ObjexxFCL::format;
	using numeric::random::uniform;

	PoseWrap pw = posewrap_from_command_line();

	std::ostringstream oss;

	vector1<utility::file::FileName> f = options::option[options::OptionKeys::smhybrid::refine]();

	for(Size ifile = 1; ifile <= f.size(); ++ifile) {
		utility::file::FileName const fname(f[ifile]);
		for(Size iter = 1; iter <= (Size)options::option[options::OptionKeys::out::nstruct](); ++iter) {
			TR << "FILE " << fname << "ITER " << iter << std::endl;

			PoseWrap pw = posewrap_from_command_line();
			pw.check_scattach_res();

			Pose ref;
			pose_from_pdb(ref,fname);
			for(Size ir = 1; ir <= pw.pose.n_residue(); ++ir) {
				pw.pose.replace_residue(ir,ref.residue(ir),false);
				for(Size ia = 1; ia <= ref.residue(ir).natoms(); ++ia) {
					pw.pose.set_xyz(AtomID(ia,ir),ref.xyz(AtomID(ia,ir)));
				}
			}
			for(Size ir = 1; ir <= pw.floating_scs_.size(); ++ir) {
				Size const rsd(pw.floating_scs_[ir]);
				pose::add_lower_terminus_type_to_pose_residue(pw.pose,rsd);
				pose::add_upper_terminus_type_to_pose_residue(pw.pose,rsd);
				core::pose::add_variant_type_to_pose_residue(pw.pose,"VIRTUAL_BBCB",rsd);
				for(Size ia = 1; ia <= ref.residue(rsd).natoms(); ++ia) {
					if(pw.pose.residue(rsd).has(ref.residue(rsd).atom_name(ia))) {
						pw.pose.set_xyz(AtomID(pw.pose.residue(rsd).atom_index(ref.residue(rsd).atom_name(ia)),rsd),ref.xyz(AtomID(ia,rsd)));
					}
				}
			}
			pw.check_scattach_res();
			pw.switch_to_fa();
			pw.check_scattach_res();
			pw.dump_pdb("test.pdb");

			ScoreFunctionOP sf_fa;

			sf_fa = fa_refine_and_design(pw,1);

			report(pw,sf_fa,oss,0.0);

		}
	}
	return NULL;
}


int
main( int argc, char * argv [] )
{

	try {

	devel::init(argc,argv);


	void* (*func)(void*) = &doit;
	if( basic::options::option[basic::options::OptionKeys::smhybrid::refine].user() ) func = &doit_refine;

	if (option[ basic::options::OptionKeys::parser::view ]()) {
		protocols::viewer::viewer_main( func );
	} else {
		func(NULL);
	}



	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
