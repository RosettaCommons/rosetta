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
	// #include <basic/database/open.hh>
	#include <basic/options/option_macros.hh>
	#include <basic/options/keys/mh.OptionKeys.gen.hh>

	#include <core/chemical/AtomType.hh>
	#include <core/chemical/orbitals/OrbitalType.hh>
	#include <core/chemical/ChemicalManager.hh>
	#include <core/chemical/ResidueTypeSet.hh>
	// #include <core/conformation/symmetry/util.hh>
	#include <core/conformation/ResidueFactory.hh>
	// #include <core/conformation/util.hh>
	// #include <core/import_pose/import_pose.hh>
	// #include <core/io/silent/SilentFileData.hh>
	// #include <core/pose/PDBInfo.hh>
	#include <core/pose/Pose.hh>
	// #include <core/pose/annotated_sequence.hh>
	#include <core/pose/motif/reference_frames.hh>
	// #include <core/pose/util.hh>
	// #include <core/pose/symmetry/util.hh>
	#include <core/io/pdb/pdb_writer.hh>
	// #include <core/kinematics/MoveMap.hh>
	// #include <core/scoring/Energies.hh>
	// #include <core/scoring/EnergyGraph.hh>
	#include <core/scoring/ScoreFunction.hh>
	// #include <core/scoring/ScoreFunctionFactory.hh>
	// #include <core/scoring/ScoreTypeManager.hh>
	// #include <core/scoring/dssp/Dssp.hh>
	#include <core/scoring/etable/Etable.hh>
	// #include <core/scoring/hbonds/HBondOptions.hh>
	// #include <core/scoring/hbonds/HBondSet.hh>
	// #include <core/scoring/hbonds/hbonds.hh>
	#include <core/pack/dunbrack/RotamerLibrary.hh>
	#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
	#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
	#include <core/pack/rotamer_set/RotamerSetFactory.hh>
	#include <core/pack/rotamer_set/RotamerSet.hh>
	#include <core/pack/task/TaskFactory.hh>
	#include <utility/graph/Graph.hh>
	#include <core/pack/packer_neighbors.hh>

	// #include <core/scoring/methods/EnergyMethodOptions.hh>
	// #include <core/scoring/packing/compute_holes_score.hh>
	// #include <core/scoring/rms_util.hh>
	// #include <core/scoring/sasa.hh>
	// #include <core/scoring/symmetry/SymmetricScoreFunction.hh>
	#include <devel/init.hh>
	#include <numeric/conversions.hh>
	#include <numeric/model_quality/rms.hh>
	#include <numeric/random/random.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyz.io.hh>
	#include <numeric/xyzVector.hh>
	// #include <protocols/idealize/IdealizeMover.hh>
	// #include <protocols/sicdock/Assay.hh>
	// #include <protocols/sic_dock/SICFast.hh>
	// #include <protocols/sic_dock/util.hh>
	// #include <protocols/sic_dock/read_biounit.hh>
	// #include <protocols/simple_moves/MinMover.hh>
	// #include <protocols/simple_moves/AddConstraintsToCurrentConformationMover.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	// #include <utility/fixedsizearray1.hh>
	// #include <utility/file/file_sys_util.hh>
	// #include <numeric/geometry/hashing/SixDHasher.hh>
	// #include <numeric/HomogeneousTransform.hh>

	// #include <apps/pilot/will/will_util.ihh>

	#include <boost/foreach.hpp>
	#include <boost/iterator/iterator_facade.hpp>

	#include <Eigen/Geometry>

	// #include <boost/archive/text_oarchive.hpp>
	// #include <boost/archive/text_iarchive.hpp>
	// #include <boost/archive/binary_oarchive.hpp>
	// #include <boost/archive/binary_iarchive.hpp>

	// // #include <boost/serialization/detail/stack_constructor.hpp>
	// #include <boost/serialization/hash_map.hpp>
	// #include <boost/serialization/hash_collections_save_imp.hpp>
	// #include <boost/serialization/hash_collections_load_imp.hpp>

static basic::Tracer TR("motif_hash_util");

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
	// using core::import_pose::pose_from_file;
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
	// using protocols::sic_dock::KMGT;


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

typedef utility::vector1<std::string> Strings;

OPT_1GRP_KEY( StringVector, mdhb, donres )
	OPT_1GRP_KEY( StringVector, mdhb, accres )
	OPT_1GRP_KEY( Real        , mdhb, tip_tol_deg )
	OPT_1GRP_KEY( Real        , mdhb, rot_samp_resl )


	void register_options() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		NEW_OPT(  mdhb::donres, "" , Strings() );
		NEW_OPT(  mdhb::accres, "" , Strings() );
		NEW_OPT(  mdhb::tip_tol_deg  , "" , 30.0 );
		NEW_OPT(  mdhb::rot_samp_resl, "" ,  7.0 );
	}

template<class T>
std::string str(T const & t, Size N=0){
	std::ostringstream oss;
	oss << t;
	std::string s = oss.str();
	while(s.size()<N) s = "0"+s;
	return s;
}

core::pack::rotamer_set::RotamerSetOP
get_rot_set(
	core::pose::Pose & pose,
	Size ir
){
	core::scoring::ScoreFunction dummy_sfxn;
	dummy_sfxn.score( pose );
	core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
	dummy_task->initialize_from_command_line();
	dummy_task->nonconst_residue_task( ir ).restrict_to_repacking();
	dummy_task->nonconst_residue_task( ir ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
	dummy_task->nonconst_residue_task( ir ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
	utility::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
	core::pack::rotamer_set::RotamerSetFactory rsf;
	core::pack::rotamer_set::RotamerSetOP rotset( rsf.create_rotamer_set( pose.residue( ir ) ) );
	rotset->set_resid( ir );
	rotset->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
	return rotset;
}

Real
get_rot_score(
	core::pose::Pose & pose,
	Size ir,
	core::pack::dunbrack::RotamerLibraryScratchSpace & scratch
){
	core::pack::rotamers::SingleResidueRotamerLibraryCAP rotlib =
		core::pack::dunbrack::RotamerLibrary::get_instance().get_rsd_library( pose.residue(ir).type() );
	return rotlib->rotamer_energy( pose.residue(ir), scratch );

}

inline void dump_pdb_atom(
	std::ostream & out,
	double x, double y, double z
){
	// std::string atomname,resname,elem;
	// int atomnum,resnum;
	// char chain;
	// bool ishet;
	// float occ,bfac;
	BOOST_VERIFY( x<10000 && x > -1000 );
	BOOST_VERIFY( y<10000 && y > -1000 );
	BOOST_VERIFY( z<10000 && z > -1000 );
	// cout << "ATOM   1604  C   GLU A 220       5.010  12.933   1.553  1.00 41.10           C" << endl;
	char buf[128];
	// std::string aname = a.atomname;
	// if( aname.size() == 1 ) aname = aname+"  ";
	// if( aname.size() == 2 ) aname = aname+" ";
	snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
		true ? "HETATM":"ATOM  ",
		1,
		"DOTS",
		"DOT",
		'A',
		1,
		x,y,z,
		1.0,
		1.0,
		"D"
	);
	out << buf;
}


/// @brief generate BCC lattice "hyper-ring" to sample rotationr around
/// the hbond axis as well as tilting the axis a up to tolang
utility::vector1<numeric::xyzMatrix<double> >
get_hbond_rotation_samples( double tolang, double reslang ){
	using namespace Eigen;
	typedef Map<Matrix3d> MMap;

	tolang = numeric::conversions::radians(tolang);
	reslang = numeric::conversions::radians(reslang);

	double spacing = sqrt(3.0)*reslang;
	int const n1 = std::ceil(M_PI*2.0/spacing);
	spacing = M_PI*2.0/n1;
	int const n2 = std::ceil(tolang/spacing);
	cout << "gen_hbond_rotations: spacing " << spacing << " n1 " << n1 << " n2 " << n2 << endl;

	// utility::io::ozstream out("test.pdb");
	// generate bcc "square rod" wrapped around in 3sphere wx plane
	utility::vector1<numeric::xyzMatrix<double> > rots;
	for(int i =  0 ; i <  n1; ++i){
	for(int j = -n2; j <= n2; ++j){
	for(int k = -n2; k <= n2; ++k){
	for(int o =  0;  o <=  1; ++o){
		// if(o && ( j==n2 || k==n2 ) ) continue;
		double const wx = i*spacing + (o?spacing/2.0:0.0) - M_PI;
		double const  w = sin( wx/2.0 );
		double const  x = cos( wx/2.0 );
		double const  y = ( j*spacing + (o?spacing/2.0:0.0) ) / 2.0;
		double const  z = ( k*spacing + (o?spacing/2.0:0.0) ) / 2.0;
		if( y*y + z*z > tolang*tolang/4.0 ) continue;
		// double const r = 1.0+y;
		// dump_pdb_atom(out,r*w*50.0,r*x*50.0,r*z*50.0);
		Quaterniond q(w,x,y,z);
		// printf("%4d %4d %4d %8.3f %8.3f %8.3f %8.3f\n",i,j,k,wx,y,z,q.norm());
		q.normalize();
		// dump_pdb_atom(out,wx*50.0,y*50.0,z*50.0);
		numeric::xyzMatrix<double> rot;
		MMap mm((double*)&rot);
		mm = q.matrix();
		// cout << rot << endl;
		// cout << q.matrix() << endl;
		// cout << endl;
		rots.push_back(rot);
	}}}}

	// out.close();

// std::exit(0);

	return rots;
}

class HBondedPairIterator
    : public boost::iterator_facade< HBondedPairIterator,
                                     Pose,
                                     boost::forward_traversal_tag,
                                     Pose const &
             						>
{
public:
	core::pose::Pose pose_,pose0_;
	core::pack::dunbrack::RotamerLibraryScratchSpace scratch_;
	core::pack::rotamer_set::RotamerSetOP rotset1_,rotset2_;
	Size irot1_,irot2_,idon_,iacc_,iorb_,ihbr_;
	Size nrots1_,nrots2_;
	utility::vector1<Size> donor_atoms_;
	utility::vector1<Size> donor_bases_;
	utility::vector1<Size> acceptor_atoms_;
	utility::vector1<utility::vector1<Size> > acceptor_orbitals_;
	utility::vector1<numeric::xyzMatrix<double> > rot_samples_;
	HBondedPairIterator(){ irot1_=0; irot2_=0; idon_=0; iacc_=0; iorb_=0; ihbr_=0; }
	HBondedPairIterator(std::string resn1, std::string resn2, double tol_ang=20.0, double tol_resl=5.0){
		rot_samples_ = get_hbond_rotation_samples(tol_ang,tol_resl);

		core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::chemical::ResidueType const & rtype1 = rts->name_map(resn1);
		core::chemical::ResidueType const & rtype2 = rts->name_map(resn2);
		core::conformation::ResidueOP res1op = core::conformation::ResidueFactory::create_residue(rtype1);
		core::conformation::ResidueOP res2op = core::conformation::ResidueFactory::create_residue(rtype2);
		pose_.append_residue_by_jump(*res1op,1);
		pose_.append_residue_by_jump(*res2op,1);
		rotset1_ = get_rot_set(pose_,1);
		rotset2_ = get_rot_set(pose_,2);
		nrots1_ = rotset1_->num_rotamers();
		nrots2_ = rotset2_->num_rotamers();
		donor_atoms_ = rtype1.Hpos_polar_sc();
		BOOST_FOREACH(Size i,donor_atoms_) donor_bases_.push_back( rtype1.atom_base(i) );
		acceptor_atoms_ = rtype2.accpt_pos_sc();
		acceptor_orbitals_.resize(acceptor_atoms_.size());
		for(Size iacc = 1; iacc <= acceptor_atoms_.size(); ++iacc){
			BOOST_FOREACH(Size j,rtype2.bonded_orbitals(acceptor_atoms_[iacc])){
				if(
					// rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::C_pi_sp2    ||
				    // rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::N_pi_sp2    ||
				    rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::N_p_sp2     ||
				    // rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_pi_sp2    ||
				    rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_p_sp2     ||
				    // rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::S_p_sp3     ||
				    // rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_pi_sp2_bb ||
				    // rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_p_sp2_bb  ||
				    rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_p_sp3
				){
					acceptor_orbitals_[iacc].push_back(j);
				}
			}
			if( acceptor_orbitals_[iacc].size()==0 ) utility_exit_with_message("no acceptor orbitals "+resn2);
		}
		if( donor_atoms_.size()==0 ) utility_exit_with_message("no donor atoms");
		if( acceptor_atoms_.size()==0 ) utility_exit_with_message("no acceptor atoms");
#ifdef DEBUG
		for(Size iacc=1; iacc<=acceptor_atoms_.size(); ++iacc)
			cout << "acceptor orbitals: " << iacc << " " << acceptor_atoms_[iacc] << " " << acceptor_orbitals_[iacc] << endl;
		cout << "donor atoms: " << donor_atoms_ << endl;
		cout << "donor nchi: " << pose_.residue(1).nchi() << endl;
		for(Size i = 1; i <= rotset1_->num_rotamers(); ++i){
			if( pose_.residue(1).nchi()==1 ) printf("donor chi %7.2f \n",rotset1_->rotamer(i)->chi(1));
			if( pose_.residue(1).nchi()==2 ) printf("donor chi %7.2f %7.2f \n",rotset1_->rotamer(i)->chi(1),rotset1_->rotamer(i)->chi(2));
			if( pose_.residue(1).nchi()==3 ) printf("donor chi %7.2f %7.2f %7.2f \n",
				rotset1_->rotamer(i)->chi(1),rotset1_->rotamer(i)->chi(2),rotset1_->rotamer(i)->chi(3));
			if( pose_.residue(1).nchi()==4 ) printf("donor chi %7.2f %7.2f %7.2f %7.2f \n",
				rotset1_->rotamer(i)->chi(1),rotset1_->rotamer(i)->chi(2),rotset1_->rotamer(i)->chi(3),rotset1_->rotamer(i)->chi(4));
		}
		cout << "acceptor nchi: " << pose_.residue(2).nchi() << endl;
		for(Size i = 1; i <= rotset2_->num_rotamers(); ++i){
			if( pose_.residue(2).nchi()==1 ) printf("acceptor chi %7.2f \n",rotset2_->rotamer(i)->chi(1));
			if( pose_.residue(2).nchi()==2 ) printf("acceptor chi %7.2f %7.2f \n",rotset2_->rotamer(i)->chi(1),rotset2_->rotamer(i)->chi(2));
			if( pose_.residue(2).nchi()==3 ) printf("acceptor chi %7.2f %7.2f %7.2f \n",
				rotset2_->rotamer(i)->chi(1),rotset2_->rotamer(i)->chi(2),rotset2_->rotamer(i)->chi(3));
			if( pose_.residue(2).nchi()==4 ) printf("acceptor chi %7.2f %7.2f %7.2f %7.2f \n",
				rotset2_->rotamer(i)->chi(1),rotset2_->rotamer(i)->chi(2),rotset2_->rotamer(i)->chi(3),rotset2_->rotamer(i)->chi(4));
		}
#endif
		irot1_=1; irot2_=1; idon_=1; iacc_=1; iorb_=1; ihbr_=1;
		update_chi();
		if( !update_hb() ) increment();

	}
private:
    friend class boost::iterator_core_access;
    Pose const & dereference() const {
    	return pose_;
    }
    void update_chi(){
		for(Size ichi=1; ichi <= pose_.residue(1).nchi(); ++ichi)
			if(pose_.residue(1).nchi()>0) pose_.set_chi(ichi,1,rotset1_->rotamer(irot1_)->chi(ichi));
		for(Size ichi=1; ichi <= pose_.residue(2).nchi(); ++ichi)
			if(pose_.residue(2).nchi()>0) pose_.set_chi(ichi,2,rotset2_->rotamer(irot2_)->chi(ichi));
    }
    bool update_hb(bool realign=true){
		// if(realign) cout << "REALIGN_HB " << irot1_ << " " << idon_ << " " << irot2_ << " " << iacc_ << " " << iorb_ << " " << ihbr_ << endl;
		// else        cout << "UPDATE_HB  " << irot1_ << " " << idon_ << " " << irot2_ << " " << iacc_ << " " << iorb_ << " " << ihbr_ << endl;
		Vec don  = Vec( pose_.residue(1).xyz(donor_atoms_[idon_]) );
		Vec donb = Vec( pose_.residue(1).xyz(donor_bases_[idon_]) );
		Vec acc  = Vec( pose_.residue(2).xyz(acceptor_atoms_[iacc_]));
		Vec acco = Vec( pose_.residue(2).orbital_xyz(acceptor_orbitals_[iacc_][iorb_]));

		Vec d = (don-donb).normalized();
		Vec a = (acc-acco).normalized();
		// check if already "aligned"
		if( realign ){
			Xform xd = Xform::align( Vec(1,0,0), don-donb );
			Xform xa = Xform::align( Vec(1,0,0), acc-acco );
			core::scoring::motif::xform_pose( pose_, xd, 1, 1 );
			core::scoring::motif::xform_pose( pose_, xa, 2, 2 );
			core::scoring::motif::xform_pose( pose_, Vec(-1.05,0,0)-pose_.residue(1).xyz(donor_atoms_   [idon_]), 1, 1 );
			core::scoring::motif::xform_pose( pose_, Vec( 0.65,0,0)-pose_.residue(2).xyz(acceptor_atoms_[iacc_]), 2, 2 );
			pose0_ = pose_;
		} else {
			for(Size ia = 1; ia <= pose_.residue(2).natoms(); ++ia){
				pose_.set_xyz( core::id::AtomID(ia,2) , rot_samples_[ihbr_] * pose0_.residue(2).xyz(ia) );
			}
		}

		for(Size ia = 1; ia <= pose_.residue(1).nheavyatoms(); ++ia){
			if(ia==4) continue;
			Real const lj1 = pose_.residue(1).atom_type(ia).lj_radius();
			for(Size ja = 1; ja <= pose_.residue(2).nheavyatoms(); ++ja){
				if( ja==4 || ( ia== donor_bases_[idon_] && ja==acceptor_atoms_[iacc_] ) ) continue;
				Real const lj2 = pose_.residue(2).atom_type(ja).lj_radius();
				Real const clash_dis = (lj1+lj2-0.5);
				Real const clash_dis2 = clash_dis*clash_dis;
				if( clash_dis2 > pose_.residue(1).xyz(ia).distance_squared( pose_.residue(2).xyz(ja) ) ){
					// cout << sqrt(clash_dis2) << " " << pose_.residue(1).xyz(ia).distance( pose_.residue(2).xyz(ja)) << endl;
					return false;
				}
			}
		}
		return true;
    }
    void increment(){
    	bool clash = true;
    	while(clash){
    		++ihbr_;
	    	if( ihbr_  > rot_samples_.size()              ){ ihbr_  = 1; ++iorb_ ; }
    		if( iorb_  > acceptor_orbitals_[iacc_].size() ){ iorb_  = 1; ++iacc_ ; }
    		if( iacc_  > acceptor_atoms_.size()           ){ iacc_  = 1; ++irot2_; }
			if( irot2_ > nrots2_                          ){ irot2_ = 1; ++idon_ ; }
			if( idon_  > donor_atoms_.size()              ){ idon_  = 1; ++irot1_; }
			if( irot1_ > nrots1_                          ){ irot1_ = irot2_ = idon_ = iacc_ = iorb_ = ihbr_ = 0;  return; }
			this->update_chi();
			clash = ! this->update_hb(ihbr_==1);
		}
    }
    bool equal(HBondedPairIterator const& o) const {
    	return o.irot1_==irot1_ && o.irot2_==irot2_ && o.iacc_==iacc_ && o.idon_==idon_ && iorb_==o.iorb_ && ihbr_==o.ihbr_;
    }

};

int main(int argc, char *argv[]) {
	register_options();
	devel::init(argc,argv);

	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::motif;

						// ballow2(i2,i1,ir) = plib->rotamer_energy( pose3.residue(2), scratch );
	ResPairMotifs motifs;
	utility::io::ozstream allout(option[mh::motif_out_file]()+".rpm.bin.gz");

	Size count = 0, totcount=0;
	BOOST_FOREACH( std::string resn1, option[mdhb::donres]() ){

		BOOST_FOREACH( std::string resn2, option[mdhb::accres]() ){
			cout << "motif_denovo_hb " << resn1 << " " << resn2 << endl;

			HBondedPairIterator end;
			HBondedPairIterator iter(
				resn1,
				resn2,
				option[mdhb::tip_tol_deg](),
				option[mdhb::rot_samp_resl]()
			);
			for(; iter != end; ++iter){
				Pose const & pose(*iter);
				++totcount;

				// Size irot1(iter.irot1_), irot2(iter.irot2_), iorb(iter.iorb_), iacc(iter.iacc_), idon(iter.idon_), ihbr(iter.ihbr_);
				// if( irot1==1 && irot2==1 && iacc==1 && idon==1 && count < 300 ){
				// if( numeric::random::uniform() < 0.001 && count < 100 ){
				// 	++count;
				// 	std::string tag = "test_"+resn1+"_"+resn2+"_"+str(irot1)+"_"+str(idon)+"_"+str(irot2)+"_"+str(irot2)+"_"+str(iacc)+"_"+str(iorb)+"_"+str(ihbr,2);
				// 	cout << tag << endl;
				// 	pose.dump_pdb(tag+".pdb");
				// }

				if(totcount%10000==9999) cout << "Samples: " << totcount+1 << endl;
				// if(totcount > 10000) break;
				Xform frm1 = core::pose::motif::get_backbone_reference_frame(pose,1);
				Xform frm2 = core::pose::motif::get_backbone_reference_frame(pose,2);
				Real6 rt = (~frm1*frm2).rt6();
				ResPairMotif motif(
					"DNOHB", // string tag,
					pose,    // Pose const & _pose,
					rt,      // Real6 const & _xform,
					1,       // Size const & _resi1,
					2,       // Size const & _resi2,
					1,       // Real const & _nbrs1,
					1,       // Real const & _nbrs2,
					0.0,     // Real const & _fa_atr,
					0.0,     // Real const & _fa_atr_sc_bb,
					0,       // Real const & _fa_atr_bb,
				   -1.5,     // Real const & _hb_sc,
					0.0,     // Real const & _hb_bb_sc,
					0,       // Real const & _hb_bb,
					0.0,     // Real const & _bfac1,
					0.0,     // Real const & _bfac2,
					false,   // bool is_sspaired,
					BB_BB    // RPM_Type _type
				);
				motifs.push_back(motif);

			}
		}
	}
	cout << motifs.size() << endl;
	if(!write_motifs_binary(allout,motifs)) utility_exit_with_message("error writing to file "+option[mh::motif_out_file]()+".rpm.bin.gz");


 }


