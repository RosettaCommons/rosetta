// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// headers
 #define MAX_UINT16 65535
 #define MAX_UINT8    255

 #include <protocols/sic_dock/read_biounit.hh>

 #include <basic/Tracer.hh>
 #include <core/chemical/AtomType.hh>
 #include <core/chemical/ChemicalManager.hh>
 #include <core/import_pose/import_pose.hh>
 #include <core/pose/PDBInfo.hh>
 #include <core/pose/Pose.hh>
 #include <core/pose/annotated_sequence.hh>
 #include <core/pose/util.hh>
 #include <core/scoring/Energies.hh>
 #include <core/scoring/EnergyGraph.hh>
 #include <core/scoring/ScoreFunctionFactory.hh>
 #include <core/scoring/ScoreTypeManager.hh>
 #include <core/scoring/dssp/Dssp.hh>
 #include <core/scoring/hbonds/HBondOptions.hh>
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

// static data
 static basic::Tracer TR("protocols.sic_dock.read_biounit");

// types
 using std::string;
 typedef numeric::xyzVector<core::Real> Vec;
 typedef numeric::xyzMatrix<core::Real> Mat;

 using core::id::AtomID;
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
 using std::string;
 using utility::file_basename;
 using utility::vector1;
 using std::endl;
 using core::import_pose::pose_from_pdb;

 namespace protocols {
 namespace sic_dock {

using std::cerr;
using std::cout;

bool
read_biounit(
	std::string const & fname,
	core::pose::Pose & pose,
	int  const max_res,
	bool const debug
){
	vector1<core::Real> tmpvec1,tmpvec2;
	return read_biounit(fname,pose,tmpvec1,tmpvec2,max_res, debug);
}


bool
read_biounit(
	std::string const & fname,
	core::pose::Pose & pose,
	utility::vector1<core::Real> & bfactors,
	utility::vector1<core::Real> & occupancy,
	int  const max_res,
	bool const debug
){
	std::map<int,char> tmpmap;
	vector1<int> tmpvec;
	int tmpint;
	return read_biounit(fname,pose,bfactors,occupancy,tmpvec,tmpmap,tmpint,max_res, debug);
}


bool
read_biounit(
	string const & fname,
	Pose & pose,
	vector1<Real> & bfactors,
	vector1<Real> & occupancy,
	vector1<int>  & pdbres,
	std::map<int,char> & pdbchain,
	int        & nresmodel1,
	int  const max_res,
	bool const debug
){
	if(debug) cerr<<"my_read_biounit: reading pdb "<<fname<<endl;
	pdbchain.clear();
	pdbres.clear();
	try {
		if(pose.n_residue()!=0) utility_exit_with_message("must fill empty pose");
		core::chemical::ResidueTypeSetCAP fars = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		utility::vector1<Pose> poses;
		pose_from_pdb( poses, *fars, fname, false );
		bool only_one_model = -12345==nresmodel1; nresmodel1 = 9999999;
		for(vector1<Pose>::const_iterator ip = poses.begin(); ip != poses.end(); ++ip){
			if((int)ip->n_residue() > max_res) { cout<<"SKIP max_res "<<max_res<<" "<<fname<<endl; return false; }
			for(Size i = 1; i <= ip->n_residue(); ++i){
				if( i==1 || ip->residue(i).is_lower_terminus() || ip->residue(i).is_ligand() ){
					if(!ip->residue(i).is_protein()) continue;
					bool is_new_chain = ( 1 == i || ip->chain(i) != ip->chain(i-1) || ip->pdb_info()->chain(i) != ip->pdb_info()->chain(i-1) );
					pose.append_residue_by_jump(ip->residue(i),1,"","",is_new_chain);
					// if(is_new_chain) cout << "new chain " <<pose.n_residue() << " " << pose.chain(pose.n_residue()) << endl;
				} else {
					pose.append_residue_by_bond(ip->residue(i));
				}
				if((int)pose.n_residue() > max_res) { cout<<"SKIP max_res "<<max_res<<" "<<fname<<endl; return false; }
				pdbres.push_back(ip->pdb_info()->number(i));
				pdbchain[pose.chain(pose.n_residue())] = ip->pdb_info()->chain(i);
				Real bfac=0.0, occ=0.0;
				for(Size ia = 1; ia <= ip->residue(i).nheavyatoms(); ++ia){
					if(ip->pdb_info()->temperature(i,ia)<0.1) bfac = 9e9;
					bfac += ip->pdb_info()->temperature(i,ia);
					occ  += ip->pdb_info()->occupancy  (i,ia);
					// TR<<i <<" "<< ia <<" "<< ip->pdb_info()->temperature(i,ia)<<endl;
				}
				bfac /= (Real)ip->residue(i).nheavyatoms();
				occ  /= (Real)ip->residue(i).nheavyatoms();
				bfactors .push_back( sqrt(bfac/78.95684) );
				occupancy.push_back(      occ            );
			}
			nresmodel1 = min((int)pose.n_residue(),nresmodel1);
			if(only_one_model) break;
		}
		if(debug) pose.dump_pdb(utility::file_basename(fname)+"_ha.pdb");
	} catch(...) {
		return false;
	}
	// utility_exit_with_message("test biounit input");
	return true;
}

 } // sic_dock
 } // protocols
