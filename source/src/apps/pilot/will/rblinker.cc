// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/rblinker.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/database/open.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <apps/pilot/will/mynamespaces.ihh>

using core::kinematics::Stub;
using core::conformation::Residue;

static THREAD_LOCAL basic::Tracer TR( "rblinker" );


inline void xform_pose( core::pose::Pose & pose, Stub const & s ) {
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, s.local2global(pose.xyz(aid)) );
		}
	}
}


void read_fragdata( vector1< core::fragment::FragDataOP > & fds, utility::io::izstream & in, bool /*design = false*/ ) {
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
	in.close();
}

std::map<string, vector1<core::fragment::FragDataOP> >
get_frags_map( bool design = false ) {
	using namespace core::fragment;
	TR << "reading frags" << std::endl;
	utility::io::izstream in;
	std::map<string,vector1<FragDataOP> > fds;
	basic::database::open(in,"sampling/ss_fragfiles/EEE.fragfile"); read_fragdata(fds["EEE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/EEH.fragfile"); read_fragdata(fds["EEH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/EEL.fragfile"); read_fragdata(fds["EEL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/EHH.fragfile"); read_fragdata(fds["EHH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/ELE.fragfile"); read_fragdata(fds["ELE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/ELH.fragfile"); read_fragdata(fds["ELH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/ELL.fragfile"); read_fragdata(fds["ELL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HEE.fragfile"); read_fragdata(fds["HEE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HEH.fragfile"); read_fragdata(fds["HEH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HEL.fragfile"); read_fragdata(fds["HEL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HHE.fragfile"); read_fragdata(fds["HHE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HHH.fragfile"); read_fragdata(fds["HHH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HHL.fragfile"); read_fragdata(fds["HHL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HLE.fragfile"); read_fragdata(fds["HLE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HLH.fragfile"); read_fragdata(fds["HLH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HLL.fragfile"); read_fragdata(fds["HLL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LEE.fragfile"); read_fragdata(fds["LEE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LEH.fragfile"); read_fragdata(fds["LEH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LEL.fragfile"); read_fragdata(fds["LEL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LHH.fragfile"); read_fragdata(fds["LHH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LLE.fragfile"); read_fragdata(fds["LLE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LLH.fragfile"); read_fragdata(fds["LLH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LLL.fragfile"); read_fragdata(fds["LLL"],in,design);
	return fds;
}

core::fragment::FragSetOP
make_frag_set(Pose const & pose, string ss, std::map<string, vector1<core::fragment::FragDataOP> > fds) {
	using namespace core::fragment;
	FragSetOP frags = new ConstantLengthFragSet();
	int const stop = pose.n_residue();
	if((int)1 >= stop) return NULL;
	for( Size i = 1; i <= (Size)stop; ++i ) {
		string ss3 = ss.substr(i-1,3);
		bool mkframe = true;
		for(Size j = 0; j < ss3.size(); ++j) if(ss3[j]!='H'&&ss3[j]!='E'&&ss3[j]!='L'&&ss3[j]!='*') mkframe = false;
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
		if(frame->nr_frags()) frags->add(frame);
		TR << "make frag " << i << ": " << ss3 << std::endl;
	}
	if(frags->size() == 0) return NULL;
	return frags;
}

struct ClashCheck {
	Pose const pose_;
	utility::vector1<Vec> points_;
	ObjexxFCL::FArray3D< vector1<Vec> > cubes_;
	Vec bbl_,bbu_;
	numeric::xyzTriple< platform::Size > cube_dim_;
	Real side_inv_, neighbor_cutoff_, neighbor_cutoff_sq_;
	ClashCheck(Pose const & pose_in, Real clash_dis_in = -1.0) : pose_(pose_in) {
		if(clash_dis_in > 0) init_clash_check( clash_dis_in );
		else init_clash_check( basic::options::option[basic::options::OptionKeys::rblinker::clash_dis]() );
	}
	void init_clash_check(Real neighbor_cutoff) {
		using numeric::min;
		using numeric::max;
		using numeric::square;
		typedef  numeric::xyzTriple< platform::Size >  CubeDim; // Cube dimensions
		typedef  numeric::xyzTriple< platform::Size >  CubeKey; // Cube index-triple key

		neighbor_cutoff_ = neighbor_cutoff;
		neighbor_cutoff_sq_ = ( neighbor_cutoff*neighbor_cutoff);
		points_.resize(pose_.n_residue()*5);
		for(Size i = 0; i < pose_.n_residue(); ++i) {
			for(Size j = 1; j <= 5; ++j) points_[5*i+j] = pose_.xyz(AtomID(j,i+1));
		}

		bbl_ = bbu_ = points_[1]; // Lower and upper corners of bounding box
		for ( Size ii = 2; ii <= points_.size(); ++ii ) { bbl_.min( points_[ ii ] ); bbu_.max( points_[ ii ] ); }
		bbl_ -= 10 * std::numeric_limits< Real >::epsilon();
		bbu_ += 10 * std::numeric_limits< Real >::epsilon();
		Real const side( neighbor_cutoff );
		assert( side > platform::Real( 0 ) );
		side_inv_ = ( platform::Real( 1 ) / side );
		cube_dim_ = CubeDim( // Cube dimensions
			platform::Size( std::ceil( ( bbu_.x() - bbl_.x() ) * side_inv_ ) ),             // Test that ceil values == platform::Size values
			platform::Size( std::ceil( ( bbu_.y() - bbl_.y() ) * side_inv_ ) ),
			platform::Size( std::ceil( ( bbu_.z() - bbl_.z() ) * side_inv_ ) )
		);

		cubes_.dimension( cube_dim_.x(), cube_dim_.y(), cube_dim_.z() );

		for ( Size i = 1; i <= points_.size(); ++i ) {
			Vec const pp( points_[ i ]);
			CubeKey const cube_key(
				platform::Size( ( pp.x() - bbl_.x() ) * side_inv_ ) + 1,
				platform::Size( ( pp.y() - bbl_.y() ) * side_inv_ ) + 1,
				platform::Size( ( pp.z() - bbl_.z() ) * side_inv_ ) + 1
			);
			assert( cube_key.x() <= cube_dim_.x() );
			assert( cube_key.y() <= cube_dim_.y() );
			assert( cube_key.z() <= cube_dim_.z() );
			Size i_index = cubes_.index( cube_key.x(), cube_key.y(), cube_key.z() );
			cubes_[ i_index ].push_back( pp );
		}
	}
	inline bool clash_check(Vec const & pp) {
		// doesn't help!
		// if( pp.x() + neighbor_cutoff_ < bbl_.x() ) return true;
		// if( pp.x() - neighbor_cutoff_ > bbu_.x() ) return true;
		// if( pp.y() + neighbor_cutoff_ < bbl_.y() ) return true;
		// if( pp.y() - neighbor_cutoff_ > bbu_.y() ) return true;
		// if( pp.z() + neighbor_cutoff_ < bbl_.z() ) return true;
		// if( pp.z() - neighbor_cutoff_ > bbu_.z() ) return true;
		Size const icx( Size( ( pp.x() - bbl_.x() ) * side_inv_ ) + 1 );
		Size const icy( Size( ( pp.y() - bbl_.y() ) * side_inv_ ) + 1 );
		Size const icz( Size( ( pp.z() - bbl_.z() ) * side_inv_ ) + 1 );
		for ( Size ix = max( icx, Size( 2 ) ) - 1,  ixe = min( icx + 1, cube_dim_.x() ); ix <= ixe; ++ix ) {
			for ( Size iy = max( icy, Size( 2 ) ) - 1, iye = min( icy + 1, cube_dim_.y() ); iy <= iye; ++iy ) {
				for ( Size iz = max( icz, Size( 2 ) ) - 1, ize = min( icz + 1, cube_dim_.z() ); iz <= ize; ++iz ) {
					Size cube_index = cubes_.index( ix, iy, iz );
					if ( cubes_[ cube_index ].size() != 0 ) { // Cube exists
						for ( vector1<Vec>::iterator ia = cubes_[ cube_index ].begin(), iae = cubes_[ cube_index ].end(); ia != iae; ++ia ) {
							Vec const j( *ia );
							Real const d_sq( distance_squared( pp, j ) );
							if ( d_sq <= neighbor_cutoff_sq_ ) {
								return false;
							}
						}
					}
				}
			}
		}
		return true;
	}
	inline bool clash_check_trimer(Pose & pose, Size refrsd) {
		Stub stubl(pose_.xyz(AtomID(2,1)),pose_.xyz(AtomID(2,2)),pose_.xyz(AtomID(2,3)));
		Stub stub1(pose .xyz(AtomID(2,1)),pose .xyz(AtomID(2,2)),pose .xyz(AtomID(2,3)));
		for(Size i = 1; i <= pose.residue(refrsd).nheavyatoms(); ++i) {
			if(i > 9) if( ! clash_check( stubl.local2global(stub1.global2local(pose.xyz(AtomID(i,refrsd+0*pose_.n_residue())))) ) ) return false;
			if( ! clash_check( stubl.local2global(stub1.global2local(pose.xyz(AtomID(i,refrsd+1*pose_.n_residue())))) ) ) return false;
			if( ! clash_check( stubl.local2global(stub1.global2local(pose.xyz(AtomID(i,refrsd+2*pose_.n_residue())))) ) ) return false;
		}
		Vec cen = stubl.local2global(stub1.global2local(Vec(0,0,0)));
		Vec axs = stubl.local2global(stub1.global2local(Vec(0,0,1)));
		axs = axs - cen;
		Mat rot = rotation_matrix_degrees(axs,120.0);
		for(vector1<Vec>::iterator i = points_.begin(); i != points_.end(); ++i) {
			if( ! clash_check( rot*(*i-cen)+cen ) ) return false;
		}
		return true;
	}
	inline bool clash_check(Pose & pose, Size refrsd) {
		Stub stubl(pose_.xyz(AtomID(2,1)),pose_.xyz(AtomID(2,2)),pose_.xyz(AtomID(2,3)));
		Stub stub1(pose .xyz(AtomID(2,1)),pose .xyz(AtomID(2,2)),pose .xyz(AtomID(2,3)));
		for(Size i = 9; i <= pose.residue(refrsd).nheavyatoms(); ++i) {
			if( ! clash_check( stubl.local2global(stub1.global2local(pose.xyz(AtomID(i,refrsd+0*pose_.n_residue())))) ) ) return false;
		}
		return true;
	}
	inline bool clash_check(Stub & stub, Vec pos) {
		Stub stubl(pose_.xyz(AtomID(2,1)),pose_.xyz(AtomID(2,2)),pose_.xyz(AtomID(2,3)));
		return clash_check( stubl.local2global(stub.global2local(pos)) );
	}
	inline bool clash_check_naive(Pose & pose) {
		for(Size i = 1; i <= pose.n_residue(); ++i) {
			for(Size j = 1; j <= 5; ++j) {
				Vec const xyz1( pose.xyz(AtomID(j,i)) );
				for(Size i2 = 1; i2 <= pose.n_residue(); ++i2) {
					for(Size j2 = 1; j2 <= 5; ++j2) {
						Vec const xyz2( pose.xyz(AtomID(j2,i2)) );
						Real const d_sq( distance_squared( xyz1, xyz2 ) );
						if ( d_sq <= neighbor_cutoff_sq_ ) {
							return false;
						}
					}
				}
			}
		}
		return true;
	}
};

// petf   37   42   45   75
// hyd   101  156  334  338
// psI  1492 1536 1489 1495
inline Real hyd_petf_sf4_dis(Pose const & hyd, Stub const & shyd, Pose const & petf, Stub const & spetf ) {
	Vec sf4A = ( hyd.xyz(AtomID(5, 101))+ hyd.xyz(AtomID(5, 156))+ hyd.xyz(AtomID(5, 334))+ hyd.xyz(AtomID(5, 338)))/4.0;
	Vec sf4B = (petf.xyz(AtomID(5,  37))+petf.xyz(AtomID(5,  42))+petf.xyz(AtomID(5,  45))+petf.xyz(AtomID(5,  75)))/4.0;
	return shyd.local2global(sf4A).distance(spetf.local2global(sf4B));
}
inline Real psI_petf_sf4_dis(Pose const & psI, Pose const & petf, Stub const & spetf ) {
	Vec sf4A = ( psI.xyz(AtomID(5,1492))+ psI.xyz(AtomID(5,1536))+ psI.xyz(AtomID(5,1489))+ psI.xyz(AtomID(5,1495)))/4.0;
	Vec sf4B = (petf.xyz(AtomID(5,  37))+petf.xyz(AtomID(5,  42))+petf.xyz(AtomID(5,  45))+petf.xyz(AtomID(5,  75)))/4.0;
	return sf4A.distance(spetf.local2global(sf4B));
}

// get stup that aligns r1 to r2
Stub getxform(Residue const & r1, Residue const & r2) {
	Stub s;
	s.M = alignVectorSets(r1.xyz(1)-r1.xyz(2),r1.xyz(3)-r1.xyz(2),r2.xyz(1)-r2.xyz(2),r2.xyz(3)-r2.xyz(2));
	s.v = r2.xyz(2)-s.M*r1.xyz(2);
	return s;
}

class BBMover : public protocols::moves::Mover {
	Size start_,stop_;
	Real mag_;
public:
	BBMover(Size start, Size stop, Real mag) : start_(start),stop_(stop),mag_(mag) {}
	void apply(core::pose::Pose & pose) {
		Size i = start_-1 + std::ceil(uniform()*(stop_-start_+1));
		if(uniform()<0.5) pose.set_phi(i,pose.phi(i)+gaussian()*mag_);
		else              pose.set_psi(i,pose.psi(i)+gaussian()*mag_);
	}
	std::string get_name() const { return "BBMover"; }
};

void* doit(void*) {
	using namespace core;
	using namespace chemical;
	using namespace conformation;
	using namespace pose;
	using namespace protocols;
	using namespace moves;
	using namespace ObjexxFCL::format;
	using numeric::random::uniform;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ResidueTypeSetCAP cenresset = ChemicalManager::get_instance()->residue_type_set( CENTROID );

	Pose lnk1,lnk2,hydinit,psIinit,petfinit;
	pose_from_pdb(psIinit ,*cenresset,"input/psI_0001_strip_0001.pdb");
	pose_from_pdb(hydinit ,*cenresset,"input/hyda1_0001_strip_0001.pdb");
	pose_from_pdb(petfinit,*cenresset,"input/PetF_0001.pdb");
	Pose const psI(psIinit),hyd(hydinit),petf(petfinit);
	ClashCheck clashpsI (psI ,3.5);
	ClashCheck clashhyd (hyd ,3.5);
	ClashCheck clashpetf(petf,3.5);
	Size linklen1 = option[rblinker::linker1_size]();
	Size linklen2 = option[rblinker::linker2_size]();
	Size psIroot1 = option[rblinker::linker1_root]();
	Size psIroot2 = option[rblinker::linker2_root]();
	Size linkstart1=9e9,linkstart2=9e9,hydroot=9e9,petfroot=9e9,linkend1=9e9,linkend2=9e9;
	bool nolinker = option[rblinker::nolinker]();
	bool clashcheck_hyd  = linklen1 > 0;
	bool clashcheck_petf = linklen2 > 0;
	string seq1, seq2;
	string ggs = "GGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGS";
	string lll = "***************************************************************************************";
	seq1 = "A" + ggs.substr(0,linklen1) + "A";
	seq2 = "A" + ggs.substr(0,linklen2) + "A";
	string ss1 = "_" + lll.substr(0,linklen1) + "_";
	string ss2 = "_" + lll.substr(0,linklen2) + "_";

	Size noutput = 0;

	if( psIroot1 > psI.n_residue() ) utility_exit_with_message("linker1_root outside of psI!");
	if( psIroot2 > psI.n_residue() ) utility_exit_with_message("linker2_root outside of psI!");

	TR << "linker1 size " << linklen1 << " psIroot1: " << psIroot1 << " seq1: " << seq1 << std::endl;
	TR << "linker2 size " << linklen2 << " psIroot2: " << psIroot2 << " seq2: " << seq2 << std::endl;

	// make lnk1
	if( linklen1 != 0 ) {
		if( psI.residue(psIroot1).is_upper_terminus() ) {
			make_pose_from_sequence(lnk1,"A",*cenresset,false);
			lnk1.replace_residue(1,psI.residue(psIroot1),false);
			remove_upper_terminus_type_from_pose_residue(lnk1,1);
			for(Size i = 2; i <= seq1.size(); ++i) {
				string name3 = name_from_aa(aa_from_oneletter_code(seq1[i-1]));
				lnk1.append_polymer_residue_after_seqpos(*ResidueFactory::create_residue(cenresset->name_map(name3)),lnk1.n_residue(),true);
			}
			lnk1.set_omega(1,180);
			for(Size i = 2; i <= lnk1.n_residue(); ++i) {	lnk1.set_phi(i,-60); lnk1.set_psi(i,-45); lnk1.set_omega(i,180); }
			//lnk1.set_phi  (lnk1.n_residue(),hyd.phi  (1));
			//lnk1.set_psi  (lnk1.n_residue(),hyd.psi  (1));
			//lnk1.set_omega(lnk1.n_residue(),hyd.omega(1));
			hydroot = 1;
			linkstart1 = 1;
			linkend1 = lnk1.n_residue();
		} else if ( psI.residue(psIroot1).is_lower_terminus() ) {
			make_pose_from_sequence(lnk1,"A",*cenresset,false);
			lnk1.replace_residue(lnk1.n_residue(),psI.residue(psIroot1),false);
			remove_lower_terminus_type_from_pose_residue(lnk1,1);
			for(Size i = seq1.size()-1; i > 0; --i) {
				string name3 = name_from_aa(aa_from_oneletter_code(seq1[i-1]));
				lnk1.prepend_polymer_residue_before_seqpos(*ResidueFactory::create_residue(cenresset->name_map(name3)),1,true);
			}

			//lnk1.set_omega(lnk1.n_residue(),180); // should already be ok...
			for(Size i = 1; i <= lnk1.n_residue()-1; ++i) {	lnk1.set_phi(i,-60); lnk1.set_psi(i,-45); lnk1.set_omega(i,180); }
			//lnk1.set_phi  (lnk1.n_residue(),hyd.phi  (hyd.n_residue()));
			//lnk1.set_psi  (lnk1.n_residue(),hyd.psi  (hyd.n_residue()));
			//lnk1.set_omega(lnk1.n_residue(),hyd.omega(hyd.n_residue()));
			hydroot = hyd.n_residue();
			linkstart1 = lnk1.n_residue();
			linkend1 = 1;
		} else {
			utility_exit_with_message("psIroot1 must be a terminus!");
		}
	}

	// make lnk2
	if( linklen2 != 0 ) {
		if( psI.residue(psIroot2).is_upper_terminus() ) {
			make_pose_from_sequence(lnk2,"A",*cenresset,false);
			lnk2.replace_residue(1,psI.residue(psIroot2),false);
			remove_upper_terminus_type_from_pose_residue(lnk2,1);
			for(Size i = 2; i <= seq2.size(); ++i) {
				string name3 = name_from_aa(aa_from_oneletter_code(seq2[i-1]));
				lnk2.append_polymer_residue_after_seqpos(*ResidueFactory::create_residue(cenresset->name_map(name3)),lnk2.n_residue(),true);
			}
			lnk2.set_omega(1,180);
			for(Size i = 2; i <= lnk2.n_residue(); ++i) {	lnk2.set_phi(i,-60); lnk2.set_psi(i,-45); lnk2.set_omega(i,180); }
			//lnk2.set_phi  (lnk2.n_residue(),petf.phi  (1));
			//lnk2.set_psi  (lnk2.n_residue(),petf.psi  (1));
			//lnk2.set_omega(lnk2.n_residue(),petf.omega(1));
			petfroot = 1;
			linkstart2 = 1;
			linkend2 = lnk2.n_residue();
		} else if ( psI.residue(psIroot2).is_lower_terminus() ) {
			make_pose_from_sequence(lnk2,"A",*cenresset,false);
			lnk2.replace_residue(lnk2.n_residue(),psI.residue(psIroot2),false);
			remove_lower_terminus_type_from_pose_residue(lnk2,1);
			for(Size i = seq2.size()-1; i > 0; --i) {
				string name3 = name_from_aa(aa_from_oneletter_code(seq2[i-1]));
				lnk2.prepend_polymer_residue_before_seqpos(*ResidueFactory::create_residue(cenresset->name_map(name3)),1,true);
			}
			//lnk2.set_omega(lnk2.n_residue(),180); // should already be ok...
			for(Size i = 1; i <= lnk2.n_residue()-1; ++i) {	lnk2.set_phi(i,-60); lnk2.set_psi(i,-45); lnk2.set_omega(i,180); }
			//lnk2.set_phi  (lnk2.n_residue(),petf.phi  (petf.n_residue())); // should be ok because copied residue
			//lnk2.set_psi  (lnk2.n_residue(),petf.psi  (petf.n_residue()));
			//lnk2.set_omega(lnk2.n_residue(),petf.omega(petf.n_residue()));
			petfroot = petf.n_residue();
			linkstart2 = lnk2.n_residue();
			linkend2 = 1;
		} else {
			utility_exit_with_message("psIroot2 must be a terminus!");
		}
	}

	// for(Size i = 1; i <= psI  .n_residue(); ++i ) psI  .residue(i).chain(1);
	// for(Size i = 1; i <= lnk1 .n_residue(); ++i ) lnk1 .residue(i).chain(2);
	// for(Size i = 1; i <= lnk2 .n_residue(); ++i ) lnk2 .residue(i).chain(3);
	// for(Size i = 1; i <= hyda1.n_residue(); ++i ) hyda1.residue(i).chain(4);
	// for(Size i = 1; i <= petf .n_residue(); ++i ) petf .residue(i).chain(5);

	if (option[ basic::options::OptionKeys::parser::view ]()) {
		protocols::viewer::add_conformation_viewer(lnk2.conformation(),"rblinker",1000,1000);
	}

	//std::map<string, vector1<core::fragment::FragDataOP> > fds = get_frags_map();
	//fragment::FragSetOP frags1 = make_frag_set(lnk1,ss1,fds);
	//protocols::moves::MoverOP fragins1 = new protocols::abinitio::ClassicFragmentMover(frags1);
	//fragment::FragSetOP frags2 = make_frag_set(lnk2,ss2,fds);
	//protocols::moves::MoverOP fragins2 = new protocols::abinitio::ClassicFragmentMover(frags2);
  protocols::moves::MoverOP fragins1 = new BBMover(2,lnk1.n_residue()-1,60.0);
  protocols::moves::MoverOP fragins2 = new BBMover(2,lnk2.n_residue()-1,60.0);

	scoring::ScoreFunctionOP sf = scoring::ScoreFunctionFactory::create_score_function( "score3" );
	sf->set_weight(core::scoring::rama,1.0);

	// psI.dump_pdb("test_psI.pdb");

	vector1<Size> histpsI(100);
	vector1<Size> histhyd(100);

	Pose last1 = lnk1;
	Pose last2 = lnk2;
	Real lastscore = 9e9;
	Size naccept = 0, ntries1 = 0, nclash1 = 0, ntries2 = 0, nclash2 = 0;
	Stub prevstubhyd (Mat::identity(),Vec(0,0,0));
	Stub prevstubpetf(Mat::identity(),Vec(0,0,0));
	for(int ITER = 1; ITER <= basic::options::option[basic::options::OptionKeys::out::nstruct](); ++ITER) {
		//TR << "fold linker1 " << ITER << std::endl;
		Stub shyd,spetf;
		while( nolinker || linklen1 > 0 ) { // always true if sampling lnk1
			//TR << "fold linker1" << std::endl;
			ntries1++;
			// move
			if( nolinker ) {
				shyd = prevstubhyd;
				shyd.v  += Vec( gaussian(), gaussian(), gaussian() );
				Vec axis = Vec( gaussian(), gaussian(), gaussian() ).normalized();
				shyd.M *= rotation_matrix_degrees(axis, 2*gaussian() );
			} else {
				if(ITER > 10) lnk1 = last1; // having this kills diversity??
				MonteCarloOP mc1 = new MonteCarlo( lnk1, *sf, 2.0 );
				mc1->set_autotemp( true, 2.0 ); mc1->set_temperature( 2.0 );
				RepeatMover( new TrialMover(fragins1,mc1), option[rblinker::moves_per_iter]() ).apply( lnk1 );
				shyd = getxform(hyd.residue(hydroot),lnk1.residue(linkend1));
			}

			// clash check hyd1a against linker
			bool clash = false;
			if( !nolinker && linklen1 > 0 ) {
				for(Size i = 1; i <= lnk1.n_residue(); ++i) {
					bool check_psI = true;
					bool check_hyd = true;
					if( linkstart1==1 ) {
						if(i <=                  2) check_psI = false;
						if(i >= lnk1.n_residue()-1) check_hyd = false;
					} else {
						if(i <=                  2) check_hyd = false;
						if(i >= lnk1.n_residue()-1) check_psI = false;
					}
					for(Size j = 1; j <= 4; ++j) {
						if( check_psI ) {
							if(!clashpsI.clash_check( lnk1.xyz(AtomID(j,i)) )) {
								clash = true;
								//TR << "     linker1 clashes with psI " << i << std::endl;
								//lnk1.dump_pdb("test.pdb");
								//std::exit(-1);
								break;
							}
						}
						if( check_hyd ) {
							if( !clashhyd.clash_check(shyd.global2local(lnk1.xyz(AtomID(j,i)))) ) {
								clash = true;
								//TR << "     linker1 clashes with hyd " << i << std::endl;
								break;
							}
						}
					}
					if(clash) break;
				}
				if(clash) continue;
			}
			if( clashcheck_hyd ) {
				for(Size i = 1; i <= hyd.n_residue(); ++i) {
					for(Size j = 1; j <= 4; ++j) {
						if( !clashpsI.clash_check(shyd.local2global(hyd.xyz(AtomID(j,i)))) ) clash = true;
						if(clash) break;
					}
					if(clash) break;
				}
				if(clash) {
					//TR << "     hyd/psI clash" << std::endl;
					continue;
				}
			}
			//TR << "linker1 no clash!" << std::endl;
			break;
		}

		while( option[rblinker::nolinker]() || linklen2 > 0 ) { // only true if sampling petf
			//TR << "fold linker2" << std::endl;
			ntries2++;
			// move
			if( nolinker ) {
				spetf = prevstubpetf;
				spetf.v += Vec( gaussian(), gaussian(), gaussian() );
				Vec axis = Vec( gaussian(), gaussian(), gaussian() ).normalized();
				spetf.M *= rotation_matrix_degrees(axis, 2*gaussian() );
			} else {
				if(ITER > 10) lnk2 = last2; // having this kills diversity??
				MonteCarloOP mc2 = new MonteCarlo( lnk2, *sf, 2.0 );
				mc2->set_autotemp( true, 2.0 ); mc2->set_temperature( 2.0 );
				RepeatMover( new TrialMover(fragins2,mc2), option[rblinker::moves_per_iter]() ).apply( lnk2 );
				spetf = getxform(petf.residue(petfroot),lnk2.residue(linkend2));
				if(spetf.v==prevstubpetf.v) continue;
			}

			bool clash = false;

			// clash check petf against linker
			if( !nolinker && linklen2 > 0 ) {
				for(Size i = 1; i <= lnk2.n_residue(); ++i) {
					bool check_psI  = true;
					bool check_petf = true;
					if( linkstart2==1 ) {
						if(i <=                  2) check_psI  = false;
						if(i >= lnk2.n_residue()-1) check_petf = false;
					} else {
						if(i <=                  2) check_petf = false;
						if(i >= lnk2.n_residue()-1) check_psI  = false;
					}
					for(Size j = 1; j <= 4; ++j) {
						if( check_psI ) {
							if(!clashpsI.clash_check( lnk2.xyz(AtomID(j,i)) )) {
								clash = true;
								//TR << "     linker2 clashes with psI " << i << std::endl;
								//lnk1.dump_pdb("test.pdb");
								//std::exit(-1);
								break;
							}
						}
						if( check_petf ) {
							if( !clashpetf.clash_check(spetf.global2local(lnk2.xyz(AtomID(j,i)))) ) {
								clash = true;
								//TR << "     linker2 clashes with petf " << i << std::endl;
								break;
							}
						}
					}
					if(clash) break;
				}
				if(clash) continue;
				//TR << "linker2 no clash" << std::endl;
			}

			// // clash check petf against linker
			// if( !nolinker && linklen2 > 0 ) {
			// 	for(Size i = 1; i <= lnk2.n_residue()-2; ++i) {
			// 		for(Size j = 1; j <= 4; ++j) {
			// 			if(i+1 < lnk2.n_residue()) clash = clash || !clashpsI .clash_check(                   lnk2.xyz(AtomID(j,i)) );
			// 			if(i   > 2)                clash = clash || !clashpetf.clash_check(spetf.global2local(lnk2.xyz(AtomID(j,i))));
			// 			if(clash) break;
			// 		}
			// 		if(clash) break;
			// 	}
			// 	if(clash) {
			// 		nclash++;
			// 		continue;
			// 	}
			// }

			// clash check petf against PsI
			if( clashcheck_petf ) {
				for(Size i = 1; i <= petf.n_residue(); ++i) {
					for(Size j = 1; j <= 4; ++j) {
						if( !clashpsI.clash_check(spetf.local2global(petf.xyz(AtomID(j,i)))) ) clash = true;
						if(clash) break;
					}
					if(clash) break;
				}
				if(clash) {
					//	nclash++;
					continue;
				}
			}
			// clash check petf against hyd1a
			if( clashcheck_petf && clashcheck_hyd ) {
				for(Size i = 1; i <= petf.n_residue(); ++i) {
					for(Size j = 1; j <= 4; ++j) {
						if( !clashhyd.clash_check(shyd.global2local(spetf.local2global(petf.xyz(AtomID(j,i))))) ) clash = true;
						if(clash) break;
					}
					if(clash) break;
				}
				if(clash) continue;
			}

			// clash check linkers against each other
			if( clashcheck_petf && clashcheck_hyd ) {
				for( Size ir = 1; ir <= lnk1.n_residue(); ++ir ) {
				// std::cerr << lnk1.residue(ir).name() << " " << lnk1.residue(ir).nheavyatoms() << std::endl;
				for( Size jr = 1; jr <= lnk2.n_residue(); ++jr ) {
				for( Size ia = 1; ia <= lnk1.residue(ir).nheavyatoms()-1; ++ia ) {
				for( Size ja = 1; ja <= lnk2.residue(jr).nheavyatoms()-1; ++ja ) {
				if( lnk1.residue(ir).atom(ia).xyz().distance_squared( lnk2.residue(jr).atom(ja).xyz() ) < 10.0 ) {
					clash = true;
					break;
				}
				}
				if(clash) break;
				}
				if(clash) break;
				}
				if(clash) break;
				}
				if(clash) break;
			}
			if(clash) continue;


			break;
		}


		Real score=0.0, psIdis=9e9, hyddis=9e9;
		if(clashcheck_hyd && clashcheck_petf) {
			hyddis = hyd_petf_sf4_dis(hyd,shyd,petf,spetf);
		 	score += -(min(0.0,hyddis-option[rblinker::dockscore_dis]()) *
		             min(0.0,hyddis-option[rblinker::dockscore_dis]()) ) * option[rblinker::dockscore_wt]()
					     + hyddis * option[rblinker::linearscore_wt]();
		} else if(clashcheck_petf) {
			psIdis = psI_petf_sf4_dis(psI,petf,spetf);
		 	score += -(min(0.0,psIdis-option[rblinker::dockscore_dis]()) *
		             min(0.0,psIdis-option[rblinker::dockscore_dis]()) ) * option[rblinker::dockscore_wt]()
					     + psIdis * option[rblinker::linearscore_wt]();
		}

		Real const boltz_factor = ( lastscore - score ) / option[rblinker::temp]();
		Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) );
		if ( numeric::random::rg().uniform() < probability ) {
			naccept++;
			lastscore = score;
			last1 = lnk1;
			last2 = lnk2;
			prevstubpetf = spetf;
			prevstubhyd  = shyd;
			//TR << "done accept " << ITER << " " << linklen << " " << lastscore << " " << (Real)naccept/(Real)ITER << " " << score /*<< " " << hyddis*/ << " " << psIdis << " ENDL" << std::endl;
			if(clashcheck_petf)                    histpsI[min((Size)100,(Size)std::ceil(psIdis))]++;
			if(clashcheck_petf && clashcheck_hyd ) histhyd[min((Size)100,(Size)std::ceil(hyddis))]++;
		} else {
			spetf = prevstubpetf; // redundant with above... but this is where it should go
			shyd  = prevstubhyd; // redundant with above... but this is where it should go
			//TR << "done reject " << ITER << " " << linklen << " " << lastscore << " " << (Real)naccept/(Real)ITER << " " << score /*<< " " << hyddis*/ << " " << psIdis << " ENDL" << std::endl;
		}

		// report
		if(ITER%1000==0) {
			TR << "ITER " << ITER << " ACCEPT " << naccept << " " << ntries1 << " " << (Real)naccept/(Real)ntries1 << " " << ntries2 << " " << (Real)naccept/(Real)ntries2 << " current dis: " << psIdis << " " << hyddis << " current score: " << score << " " << lastscore << std::endl;
			if( clashcheck_hyd && clashcheck_petf ) {
				TR << "HOUTHYD"; for(Size i = 1; i <= histhyd.size(); ++i) TR << " " << I(10,histhyd[i]); TR << " HISTHYD" << std::endl;
			} else if( clashcheck_petf ) {
				TR << "HOUTPSI"; for(Size i = 1; i <= histpsI.size(); ++i) TR << " " << I(10,histpsI[i]); TR << " HISTPSI" << std::endl;
			}
		}

		// output structure
		//std::cerr << "OUTTEST " << option[rblinker::output_pdb_bound].user() << " " << psIdis << " " << hyddis << std::endl;
		bool bound =  (psIdis <= 15.0 || hyddis <= 15.0);
		if( option[rblinker::output_pdb].user() || bound ) {
			if( ( option[rblinker::output_pdb]() && (uniform() < 0.01 || bound) ) ||
				 ( option[rblinker::output_pdb_bound]() && bound && uniform() < 1.0/(Real(noutput)+1.0) ) )
			{
				Pose tmphyd=hyd,tmppetf=petf;
				xform_pose(tmphyd ,shyd);
				xform_pose(tmppetf,spetf);
				string fn = option[out::file::o]()+"/test_"+string_of(linklen1)+"_"+string_of(linklen2)+"_"+ObjexxFCL::lead_zero_string_of(ITER,9);
				if(bound) fn += "_BOUND";
				if(bound) TR << "dumping PBD BOUND " << psIdis << " " << fn << std::endl;
				else      TR << "dumping PBD       " << psIdis << " " << fn << std::endl;
				// ozstream out(fn);
				if(linklen1 > 0) {
					lnk1   .dump_pdb(fn+"_lnk1.pdb");
					tmphyd .dump_pdb(fn+"_hyda.pdb");
				}
				if(linklen2 > 0) {
					lnk2   .dump_pdb(fn+"_lnk2.pdb");
					tmppetf.dump_pdb(fn+"_petf.pdb");
				}
				// out.close();
				noutput++;
			}
		}
	} // end ITER

	return NULL;
}


int main( int argc, char * argv [] ) {

	try {

	devel::init(argc,argv);

	void* (*func)(void*) = &doit;
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

/*

lnk1
1560 runs OK
1691 runs OK
1692 runs OK
1756 runs OK

lnk2
1560 runs OK
1691 runs OK
1692 runs OK
1756 runs OK


*/
