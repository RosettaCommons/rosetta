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

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/AtomID_Map.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/packing/PoseBalls.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/types.hh>
#include <numeric/model_quality/rms.hh>
#include <ObjexxFCL/FArray2D.hh>
// Auto-header: duplicate removed #include <core/scoring/dssp/Dssp.hh>
#include <pstream.h>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <numeric/model_quality/RmsData.hh>


std::map<std::string,utility::io::ozstream*> outs;

// this function assumes the poses are already aligned if localalign is false
core::id::AtomID_Map<core::Real>
local_rms( core::pose::Pose p1, core::pose::Pose p2, core::Real dist, bool localalign ) {
	using namespace core;
	using namespace core::id;
	AtomID_Map<core::Real> ret;
	core::pose::initialize_atomid_map_heavy_only(ret,p1,0.0);

	for( Size ir = 1; ir <= p1.n_residue(); ++ir ) {
		for( Size ia = 1; ia <= p1.residue(ir).nheavyatoms(); ++ia ) {
			ret[AtomID(ia,ir)] = 0.0;
			numeric::xyzVector<Real> const & xyz(p1.residue(ir).xyz(ia));
			int natoms = 0;
			ObjexxFCL::FArray2D< numeric::Real > p1a(3,999), p2a(3,999);
			for( Size jr = 1; jr <= p1.n_residue(); ++jr ) {
				for( Size ja = 1; ja <= p1.residue(jr).nheavyatoms(); ++ja ) {
					numeric::xyzVector<Real> const & xyz2(p1.residue(jr).xyz(ja));
					if( xyz.distance_squared(xyz2) <= dist*dist ) {
						numeric::xyzVector<Real> const & xyz2b(p2.residue(jr).xyz(ja));
						natoms++;
						p1a(1,natoms) = xyz2.x(); p2a(1,natoms) = xyz2b.x();
						p1a(2,natoms) = xyz2.y(); p2a(2,natoms) = xyz2b.y();
						p1a(3,natoms) = xyz2.z(); p2a(3,natoms) = xyz2b.z();
					}
				}
			}
			if( localalign ) {
				ret[AtomID(ia,ir)] = numeric::model_quality::rms_wrapper(natoms,p1a,p2a);
			} else {
				for( int i = 1; i<=natoms; i++ )
					for( int j = 1; j <= 3; j++ )
						ret[AtomID(ia,ir)] += (p1a(j,i)-p2a(j,i))*(p1a(j,i)-p2a(j,i));
				ret[AtomID(ia,ir)] /= natoms;
				ret[AtomID(ia,ir)] = sqrt(ret[AtomID(ia,ir)]);
			}// end if align
		}
	}
	return ret;
}


core::Real
local_rms_window(
	core::pose::Pose p1,
	core::pose::Pose p2,
	core::Size seq_st,
	core::Size seq_stop,
	bool localalign
) {
	using namespace core;
	using namespace core::id;
	core::Real ret = -9e9;

	Size natoms = 0;
	ObjexxFCL::FArray2D< numeric::Real > p1a(3,999), p2a(3,999);
	for( Size jr = seq_st; jr <= seq_stop; ++jr ) {
		for( Size ja = 1; ja <= p1.residue(jr).nheavyatoms(); ++ja ) {
			numeric::xyzVector<Real> const & xyz1(p1.residue(jr).xyz(ja));
			numeric::xyzVector<Real> const & xyz2(p2.residue(jr).xyz(ja));
			natoms++;
			p1a(1,natoms) = xyz1.x(); p2a(1,natoms) = xyz2.x();
			p1a(2,natoms) = xyz1.y(); p2a(2,natoms) = xyz2.y();
			p1a(3,natoms) = xyz1.z(); p2a(3,natoms) = xyz2.z();
		}
	}

	if( localalign ) {
		ret = numeric::model_quality::rms_wrapper(natoms,p1a,p2a);
	} else {
		ret = 0.0;
		for( Size i = 1; i<=natoms; i++ )
			for( Size j = 1; j <= 3; j++ )
				ret += (p1a(j,i)-p2a(j,i))*(p1a(j,i)-p2a(j,i));
		ret /= natoms;
		ret = sqrt(ret);
	}// end if align
	return ret;
}


// assigns SS by DSSP
utility::vector1<core::Size>
parse_ss_regions(
	core::pose::Pose & pose,
	utility::vector1<core::Size> & reg_begin_out,
	utility::vector1<core::Size> & reg_end_out
 )
{
	// using core::Size;
	// core::scoring::dssp::Dssp dssp(pose);
	// dssp.insert_ss_into_pose( pose );
	// utility::vector1<core::Size> ssreg(pose.n_residue());
	// Size regcount = 1;
	// ssreg[1] = regcount;
	// for( Size i = 2; i <= pose.n_residue()-1; ++i ) {
	// 	if( pose.secstruct(i) != pose.secstruct(i-1) ) {
	// 		if( pose.secstruct(i-1) != pose.secstruct(i+1) ) {
	// 			regcount++;
	// 		}
	// 	}
	// 	ssreg[i] = regcount;
	// }
	// ssreg[pose.n_residue()] = ssreg[pose.n_residue()-1];
	//
	// for( Size i = 1; i <= pose.n_residue(); i++ ) {
	// 	std::cout << pose.secstruct(i);
	// }
	// std::cout << std::endl;
	using core::Size;

	utility::vector1<Size> ssreg;
	Size i,ext = pose.n_residue() % 5 / 2;
	for( i = 1; i <= ext; i++ ) ssreg.push_back(1);
	Size count = 0;
	while( i+5 <= pose.n_residue() ) {
		count++;
		for( int j=0;j<5;j++ ) {
			ssreg.push_back(count);
			i++;
		}
	}
	for(i=i; i<= pose.n_residue(); i++) ssreg.push_back(count);

	reg_begin_out.push_back(1);
	for( Size i=2; i<=ssreg.size()-1;i++) {
		if( ssreg[i]!=ssreg[i+1] ) {
			reg_begin_out.push_back(i+1);
			reg_end_out.push_back(i);
		}
	}
	reg_end_out.push_back(ssreg.size());


	// for( Size i = 1; i <= ssreg.size(); ++i ) {
	// 	std::cout << pose.secstruct(i);
	// }
	// std::cout << std::endl;
	// std::string alpha = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	// alpha = alpha+alpha+alpha+alpha+alpha+alpha+alpha+alpha;
	// for( Size i = 1; i <= ssreg.size(); ++i ) {
	// 	std::cout << alpha[ssreg[i]];
	// }
	// std::cout << std::endl;

	return ssreg;
}


// must be scored!
utility::vector1<core::Real>
fill_hbond_e(
	core::pose::Pose & pose,
	core::scoring::packing::PoseBalls pb
) {
	utility::vector1<core::Real> hbond_e;
	hbond_e.resize(pb.nballs());
	core::scoring::hbonds::HBondSet hbset;
	core::scoring::hbonds::fill_hbond_set( pose, false, hbset, false );
	for( core::Size i = 1; i <= hbset.nhbonds(); ++i ) {
		core::scoring::hbonds::HBond const & hbond( hbset.hbond(i) );
		hbond_e[ pb.id_to_index(core::id::AtomID(hbond.acc_atm(),hbond.acc_res())) ] += hbond.energy();
		core::Size iabase = pose.residue(hbond.don_res()).atom_base(hbond.don_hatm());
		core::id::AtomID donparent(iabase,hbond.don_res());
		hbond_e[ pb.id_to_index(donparent) ] += hbond.energy();
		// std::cerr << "hbond " << donparent << " " << hbond.acc_atm() << "," << hbond.acc_res() << " " << hbond.energy() << std::endl;
	}
	return hbond_e;
}


std::string tag_from_pose( core::pose::Pose & pose ) {
	std::string tag( "empty_tag" );
	if ( pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		tag =	static_cast< basic::datacache::CacheableString const & >
			( pose.data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
	}
	return tag;
}


class HolesTrainingDataMover : public protocols::moves::Mover {

public:

void
apply(
	core::pose::Pose & pose
) {
	using namespace std;
	using namespace core;
	using namespace io::pdb;
	using namespace pose;
	using namespace scoring;
	using namespace packing;
	using namespace basic::options;
	using namespace utility;
	using core::Size;

	// Pose pose;
	// core::import_pose::pose_from_file(pose,fname, core::import_pose::PDB_file);

	string fname = tag_from_pose(pose);

	Size MAX_RES = 5000;
	if( pose.total_residue() > MAX_RES ) {
		std::cout << "nres > " << MAX_RES << ", skipping file " << fname << std::endl;
		return;
	}

	// core::scoring::dssp::Dssp dssp(pose); // done by PoseBalls now
	// dssp.insert_ss_into_pose( pose );
	int  hmode = basic::options::option[ OptionKeys::holes::h_mode ](); // default 0
	bool ignore_water = !basic::options::option[ OptionKeys::holes::water ]();
	PoseBalls pb( pose, hmode, ignore_water );

	vector1<Size> reg_begin,reg_end;
	vector1<Size> ssreg = parse_ss_regions(pose, reg_begin, reg_end);
	assert(ssreg.size()==pose.n_residue());

	bool HAVE_NATIVE = false;
	vector1<Real> reg_rms_abs, reg_rms_rel;
	Real gdt = -12345.6789;
	core::id::AtomID_Map<Real> absrms5,absrms7,absrms10,relrms5,relrms7,relrms10;
	if( basic::options::option[OptionKeys::in::file::native].user() ) {
		Pose native_pose;
		core::import_pose::pose_from_file(native_pose,core::options::option[OptionKeys::in::file::native](), core::import_pose::PDB_file);
		assert(native_pose.n_residue()==pose.n_residue());
		for( Size i = 1; i<= native_pose.n_residue(); ++i ) {
			assert( native_pose.residue(i).name3()==pose.residue(i).name3() );
			;
		}
		HAVE_NATIVE = true;
		// std::cout << pose.n_residue() << " " << native_pose.n_residue() << std::endl;
		core::scoring::calpha_superimpose_pose( pose, native_pose );
		gdt = core::scoring::CA_gdtmm( pose, native_pose );
		// absrms5  = local_rms( pose, native_pose, 5.0, false );
		absrms7  = local_rms( pose, native_pose, 7.0, false );
		// absrms10 = local_rms( pose, native_pose, 10., false );
		// relrms5  = local_rms( pose, native_pose, 5.0, true  );
		relrms7  = local_rms( pose, native_pose, 7.0, true  );
		// relrms10 = local_rms( pose, native_pose, 10., true  );

		for( Size i = 1; i<= reg_begin.size(); ++i ) {
			// std::cout << reg_begin[i] << "\t" << reg_end[i] << std::endl;
			reg_rms_abs.push_back( local_rms_window(pose,native_pose,reg_begin[i],reg_end[i], false ) );
			reg_rms_rel.push_back( local_rms_window(pose,native_pose,reg_begin[i],reg_end[i], true  ) );
		}

	}


	// for( core::Size ir = 1; ir <= pose.total_residue(); ir++ ) {
	// 	for( core::Size ia = 1; ia <= pose.residue(ir).nheavyatoms(); ia++ ) {
	// 		core::Size i = pb.id_to_index(id::AtomID(ia,ir));
	// 		if( "VIRT" != pose.residue(ir).atom_type(ia).name() ) {
	// 			printf("%5d %5d %3s %4s %3.3f\n",
	// 				(int)ir,(int)ia,
	// 				pb.res_name(i).c_str(),
	// 				pb.atom_name(i).c_str(),
	// 				pb.bfac(i)
	// 			);
	// 		} else {
	// 			printf("no atom for %d %d %d\n",(int)ir,(int)ia,(int)i);
	// 		}
	// 	}
	// }

	// for( int i = 1; i <= (int)pb.nballs(); i++ ) {
	// 	std::cerr << i <<" "<< pb.res_num(i) <<" "<< pb.res_name(i) <<" "<< pb.atom_name(i) <<" "<< pb.bfac(i) << std::endl;
	// }

	// score the pose
	core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
	sf->set_weight( core::scoring::fa_elec , 1.0 );
	(*sf)(pose);

	utility::vector1<Real> hbond_e = fill_hbond_e(pose,pb);

	time_t t = clock();
	std::cout << "starting surf calculation: " << std::endl;

	std::string cmd = basic::options::option[ OptionKeys::holes::dalphaball ]();
	redi::pstream proc( cmd + " alpha20_surf" );
	// ofstream proc("test.daball");
	proc << "NPOINTS" << endl << pb.nballs() << endl << "COORDS" << endl;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		Ball & b(pb.ball(i));
	// 	// cout << "PoseBalls " << i << " " << pb.index_to_id(i) << " " << b.x() << " " << b.y() << " " << b.z() << " " << b.r() <<  endl;
		proc << b.x() << " " << b.y() << " " << b.z() << " " << b.r() << " " << endl;
	}
	// proc << "WEIGHTS" << endl;
	// for( Size i = 1; i <= pb.nballs(); i++ ) {
	// 	proc << "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1" << endl;
	// }
	proc << "END" << endl << redi::peof;

	Size index,ialpha;
	Real val;
	for( Size a = 1; a <= 20; a++ ) {
		for( Size i = 1; i <= pb.nballs(); i++ ) {
			proc >> ialpha >> index >> val;
			if( i != index || a != ialpha ) {
				std::cout << "DAlphaBall output index mismatch " << ialpha << " " << a << "   " << i << " " << index << std::endl;
			}
			pb.set_surf(i,a,val);
		}
	}

	// output
	core::chemical::AtomTypeSetCAP atom_set(
		core::chemical::ChemicalManager::get_instance()->atom_type_set( core::chemical::FA_STANDARD ) );
	std::cout << "starting output" << std::endl;
	map<string,string> lines;
	std::string prefix = basic::options::option[ basic::options::OptionKeys::out::prefix ]();
	system(("mkdir -p "+prefix).c_str());

	string tag;
	string ERROR("");
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		// Size resnum = pb.res_num(i);
		core::Size atype = pb.atom_type(i);
		tag = string_of(atype);
		if( 18 <= atype && atype <= 21 ) tag += std::string("__") + pb.secstruct(i);
		if( 0  <  atype && atype < 100 ) {
			tag = tag + "__" + (*atom_set)[atype].name();
		} else if( atype > 100 ) { // no longer used?
			 tag = tag + "__" + (*atom_set)[atype-100].name() + "_Hpol";
		}
		if( pb.atom_type(i) < 10  ) tag = "0"+tag;
		if( pb.atom_type(i) < 100 ) tag = "0"+tag;
		lines[tag] += fname + " " + string_of(ssreg[pb.res_num(i)]) + " " + string_of(pb.res_num(i))
		              +" "+ pb.res_name(i) +" "+ pb.atom_name(i) +" ";
		if( HAVE_NATIVE ) {
			lines[tag] += (std::string)(basic::options::option[OptionKeys::in::file::native]()) + " ";
			lines[tag] += string_of(gdt) + " ";
			lines[tag] += string_of( reg_rms_abs[ ssreg[pb.res_num(i)] ] ) + " ";
			lines[tag] += string_of( reg_rms_rel[ ssreg[pb.res_num(i)] ] ) + " ";
			// lines[tag] += string_of(absrms10[pb.index_to_id(i)]) + " ";
			// lines[tag] += string_of(relrms10[pb.index_to_id(i)]) + " ";
			lines[tag] += string_of(absrms7 [pb.index_to_id(i)]) + " ";
			lines[tag] += string_of(relrms7 [pb.index_to_id(i)]) + " ";
			// lines[tag] += string_of(absrms5 [pb.index_to_id(i)]) + " ";
			// lines[tag] += string_of(relrms5 [pb.index_to_id(i)]) + " ";
		}
		lines[tag] += string_of(pb.bfac(i)) + " ";
		lines[tag] += string_of(hbond_e[i]) + " ";
		lines[tag] += string_of(pb.smooth_nb(i)) + " ";
		if( true ) {
			EnergyMap e;
			core::Real s = 0;
			if( pb.res_num(i) <= pose.n_residue() ) {
				e = pose.energies().residue_total_energies( pb.res_num(i) );
				s = pose.energies().residue_total_energy  ( pb.res_num(i) );
			}
			lines[tag] += string_of( s ) + " ";
			lines[tag] += string_of( e[ fa_atr    ] ) + " ";
			lines[tag] += string_of( e[ fa_rep    ] ) + " ";
			lines[tag] += string_of( e[ fa_sol    ] ) + " ";
			lines[tag] += string_of( e[ fa_dun    ] ) + " ";
			lines[tag] += string_of( e[ fa_pair   ] ) + " ";
			lines[tag] += string_of( e[ envsmooth ] ) + " ";
			lines[tag] += string_of( e[ p_aa_pp   ] ) + " ";
			lines[tag] += string_of( e[ fa_elec ] ) + " ";
		}
		for( Size a = 1; a <=20; a++ ) lines[tag] += " " + string_of(pb.surf(i,a));
		lines[tag] += "\n";
		if( pb.surf(i,1) < 1e-6 ) {
			ERROR += "BAD SA " + string_of(pb.surf(i,1)) + " for VDW sphere of atom ";
			ERROR += fname + " " + string_of(pb.res_num(i)) +" "+ pb.res_name(i) +" "+ pb.atom_name(i) +"\n";
		}
	}
	std::cout << "time: " << Real(clock()-t)/Real(CLOCKS_PER_SEC) << std::endl;
	// if( "" != ERROR ) {
	// 	std::cerr << ERROR << std::endl;
	// 	std::exit(-1);
	// }

	for( map<string,string>::iterator i = lines.begin(); i != lines.end(); i++ ) {
		cout << "writing: " << i->first << endl;
		utility::io::ozstream out( prefix+'/'+i->first+".lr", std::ios_base::app );
		out << i->second;
		out.close();
	}

	std::cout << "DONEDATAOUTPUT " << fname << std::endl;

}

};


int
main (int argc, char *argv[])
{

	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility;

	devel::init( argc, argv );

	HolesTrainingDataMover holes;

   protocols::jobdist::not_universal_main( holes );

	for( std::map<std::string,utility::io::ozstream*>::iterator i = outs.begin(); i != outs.end(); ++i ) {
		i->second->close();
		delete i->second;
	}

	return 0;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

// int
// main (int argc, char *argv[])
// {
//
//
// 	devel::init( argc, argv );
//
//   using namespace basic::options;
//   using namespace basic::options::OptionKeys;
//   using namespace utility;
//
// 	if( option[ in::file::s ].user() ) {
//   	vector1<file::FileName> files( option[ in::file::s ]() );
//   	for( size_t i = 1; i <= files.size(); ++i ) {
//     	test( files[i] );
//   	}
// 	} else if( option[ in::file::l ].user() ) {
//   		vector1<file::FileName> files( option[ in::file::l ]() );
//   		for( size_t i = 1; i <= files.size(); ++i ) {
// 			utility::io::izstream list( files[i] );
// 			std::string fname;
// 			while( list >> fname ) {
// 				// std::cerr << "'" << fname << "'" << std::endl;
//     		test( fname );
// 			}
//   		}
// 	}
//
// 	for( std::map<std::string,utility::io::ozstream*>::iterator i = outs.begin(); i != outs.end(); ++i ) {
// 		i->second->close();
// 		delete i->second;
// 	}
//
// 	return 0;
//
//
// }
