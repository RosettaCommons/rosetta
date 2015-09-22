/// @file
/// @brief

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <devel/init.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>

//options
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

/// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <fstream>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>

//hydrogen bond
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>

//kinematics
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

//others
#include <ObjexxFCL/FArray2D.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>

//boost
#include <boost/lexical_cast.hpp>

//Tracer
#include <basic/basic.hh>
#include <basic/Tracer.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot_apps.shilei.generate_hbond_geometry" );
        
using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;
using namespace ObjexxFCL;
using namespace conformation;
using namespace scoring;
using namespace protocols;
using utility::vector1;
using std::string;
using utility::to_string;
using numeric::constants::f::pi;
using numeric::constants::f::pi_2;


OPT_1GRP_KEY(Integer,generate_hbond_geometry,resid)
OPT_1GRP_KEY(String,generate_hbond_geometry,atomname)
OPT_1GRP_KEY(Real,generate_hbond_geometry,totalE)
OPT_1GRP_KEY(Real,generate_hbond_geometry,hbondE)

//check for clashes based on atom distances
//only check atoms in AtomID vector
//should be way faster than calculating entire score
bool
fast_clash_check(
  Pose const & pose,
  vector1< id::AtomID > const check_atids,
  Real const clash_dist_cut
)
{
  Real const clash_dist2_cut( clash_dist_cut * clash_dist_cut );
  for( Size iatid = 1; iatid <= check_atids.size(); ++iatid ){
    Vector const at1_xyz( pose.xyz( check_atids[ iatid ] ) );
    for( Size res2 = 1; res2 <= pose.total_residue(); ++res2 ){
      for( Size at2 = 1; at2 <= pose.residue( res2 ).natoms(); ++at2 ){
        //skip virtual atoms!
        if( pose.residue( res2 ).atom_type( at2 ).lj_wdepth() == 0.0 ) continue;
        id::AtomID atid2( at2, res2 );
        //skip if atid2 is in check_atids
        bool skip_at2( false );
        for( Size jatid = 1; jatid <= check_atids.size(); ++jatid ){
          if( atid2 == check_atids[ jatid ] ){ skip_at2 = true; break; }
        }
        if( skip_at2 ) continue;
        Real const dist2( at1_xyz.distance_squared( pose.xyz( atid2 ) ) );
        if( dist2 < clash_dist2_cut ){
          //TR_unsat << "CLASH!: " << check_atids[ iatid ] << " - " << atid2 <<
          //   " = " << dist2 << std::endl;
          return true;
        }
      }
    }
  }
  return false;
}

///////////////////////////////////////////////////////////////////////////////
class generate_hbond_geometry : public protocols::moves::Mover {
public:
	generate_hbond_geometry(){}
	void apply( pose::Pose & pose) {

	//input should be residue number and atom name
	core::Size seqpos=basic::options::option[ basic::options::OptionKeys::generate_hbond_geometry::resid];
	std::string atomname=basic::options::option[ basic::options::OptionKeys::generate_hbond_geometry::atomname];
	core::Real totalE=basic::options::option[ basic::options::OptionKeys::generate_hbond_geometry::totalE];
	core::Real hbondE=basic::options::option[ basic::options::OptionKeys::generate_hbond_geometry::hbondE];
	TR << "Find hydrogen bond to " << " residue " << seqpos << " atom name " << atomname << std::endl;


	//initialize
	scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function() ;
	core::scoring::hbonds::HBondDatabaseCOP hb_database_ = core::scoring::hbonds::HBondDatabase::get_database( "sp2_elec_params" );

	//copy from SemiExplicitWaterUnsatisfiedPolarsCalculator
        using namespace core;
        using namespace id;
        using namespace conformation;
        using namespace scoring;
        using namespace scoring::hbonds;
        using namespace kinematics;

	Residue rsd( pose.residue( seqpos ) );
	Pose ref_pose( pose );
	//TR << pose.total_residue() << std::endl;

	//Use ALA because of Cb
	Residue new_rsd( *ResidueFactory::create_residue( rsd.residue_type_set().name_map( "ALA" ) ) );
	Size atomno=rsd.atom_index(atomname);
	Size new_atomno=new_rsd.atom_index("H");

	pose.append_residue_by_jump( new_rsd, seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), true );

	Size new_seqpos( pose.total_residue() );
	Size jump_number( pose.fold_tree().num_jump() );
	Jump jump( pose.jump( jump_number ) );
	//TR << new_seqpos << "," << jump_number << std::endl;

	//store atom ids of ALA for clash check
	vector1< id::AtomID > clash_check_atids;
	for( Size iat = 1; iat <= new_rsd.natoms(); ++iat ){
	  	clash_check_atids.push_back( id::AtomID( iat, new_seqpos ) );
	}

       //which is acceptor, donor?
       Size aatm( 0 ), acc_pos( 0 ), hatm( 0 ), don_pos( 0 );
       bool wat_is_acc( false );
       //gly is acceptor
       if( rsd.atom_type( atomno ).is_polar_hydrogen() &&
         new_rsd.atom_type( new_atomno ).is_acceptor() ){
         don_pos = seqpos;
         hatm = atomno;
         acc_pos = new_seqpos;
         aatm = new_atomno;
         wat_is_acc = true;
       }
       //or gly is donor
       else if( new_rsd.atom_type( new_atomno ).is_polar_hydrogen() &&
         rsd.atom_type( atomno ).is_acceptor() ){
         don_pos = new_seqpos;
         hatm = new_atomno;
         acc_pos = seqpos;
         aatm = atomno;
	 //TR << "donor_position: " << don_pos << std::endl;
	 //TR << "hydro atom: " << hatm << std::endl;
	 //TR << "acceptor position: " << seqpos << std::endl;
	 //TR << "accep atom: " << aatm << std::endl;
       }
       else{ utility_exit_with_message( "ERROR: res " + to_string( seqpos ) + " atom " + to_string( atomno ) +
         " res " + to_string( new_seqpos ) + " atom " + to_string( new_atomno ) + " is not HB don/acc pair!!\n" );
       }


      //now get their base atoms to get datm and batm
      Size datm( pose.residue( don_pos ).atom_base( hatm ) );
      Size batm( pose.residue( acc_pos ).atom_base( aatm ) );
      Size b2atm( pose.residue( acc_pos ).abase2( aatm ) );
      Size dbatm( pose.residue( don_pos ).atom_base( datm ) ); 
	//TR << "donor_position: " << datm << std::endl;
	//TR << "hydro atom: " << batm << std::endl;
	//TR << "acceptor position: " << b2atm << std::endl;
	//TR << "accep atom: " << dbatm << std::endl;

      //add vrt res so final torsion exists
      //chemical::ResidueTypeSet const & rsd_set( rsd.residue_type_set() );
      //conformation::ResidueOP vrt_rsd( conformation::ResidueFactory::create_residue( rsd_set.name_map( "VRT" ) ) );
      //pose.append_residue_by_jump( *vrt_rsd, pose.total_residue() );
      FoldTree f_jump( pose.fold_tree() );
      //new naive fold tree
      FoldTree f_rot( pose.total_residue() );
      //switch to chem bond so can use bond angle defs directly
      f_rot.new_chemical_bond( seqpos, new_seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), pose.total_residue() - 2 );
      pose.fold_tree( f_rot );

      Size water_hb_states_tot( 0 );
      Size water_hb_states_good( 0 );
      hbonds::HBEvalTuple hbe_type( datm, pose.residue( don_pos ), aatm, pose.residue( acc_pos ) );
      Real AHdist_min(MIN_R), AHdist_max(MAX_R); Size AHdist_steps( 10 );
      Real BAHang_min( pi - std::acos( MIN_xH ) ), BAHang_max( pi - std::acos( MAX_xH ) ); Size BAHang_steps( 10 );
      Real AHDang_min( pi - std::acos( MIN_xD ) ), AHDang_max( pi - std::acos( MAX_xD ) ); Size AHDang_steps( 10 );
      Real B2BAHchi_min( 0 ), B2BAHchi_max( pi_2 ); Size B2BAHchi_steps( 20 );
      Real BAHDchi_min( 0 ), BAHDchi_max( pi_2 ); Size BAHDchi_steps( 6 );


      for( Real AHdist = AHdist_min; AHdist <= AHdist_max; AHdist += (AHdist_max - AHdist_min)/static_cast<Real>(AHdist_steps-1)){
        for( Real BAHang = BAHang_min - 0.0001; BAHang <= BAHang_max; BAHang += (BAHang_max - BAHang_min)/static_cast<Real>(BAHang_steps-1)){
          for( Real AHDang = AHDang_min - .0001; AHDang <= AHDang_max; AHDang += (AHDang_max - AHDang_min)/static_cast<Real>(AHDang_steps-1)){
            //this loop should be longer than the water internal orientation loops
            for( Real B2BAHchi = B2BAHchi_min; B2BAHchi <= B2BAHchi_max; B2BAHchi += (B2BAHchi_max - B2BAHchi_min)/static_cast<Real>(B2BAHchi_steps-1)){
              for( Real BAHDchi = BAHDchi_min; BAHDchi <= BAHDchi_max; BAHDchi += (BAHDchi_max - BAHDchi_min)/static_cast<Real>(BAHDchi_steps-1)){
    
     //            Real BAHang( pi - std::acos( cosBAH ) );
     //            Real AHDang( pi - std::acos( cosAHD ) );
                Real cosBAH( std::cos( pi - BAHang ) );
                Real cosAHD( std::cos( pi - AHDang ) );
                //call hbonds::hbond_compute_energy( ... ) with these angles
                // and skip if is zero hb energy
                // this allows computing exactly what frac of actual HB phase space
                // can be filled w/ a favorable water molecule
                Real hb_energy( 0.0 );
                scoring::hbonds::hbond_compute_energy( *( hb_database_ ), scorefxn->energy_method_options().hbond_options(), hbe_type, AHdist, cosAHD, cosBAH, B2BAHchi, hb_energy );
                //TR << "AHdist: " << AHdist << " cosAHD: " << numeric::conversions::degrees(cosAHD) << " cosBAH: " << numeric::conversions::degrees(cosBAH) << " B2BAHchi: " << numeric::conversions::degrees(B2BAHchi) << " hb_energy: " << hb_energy << std::endl;
                if( hb_energy >= hbondE ) continue; //was not actually an hbond
    
                ++water_hb_states_tot;
    
                //reset chem  bond ftree
                pose.fold_tree( f_rot );
    
                //Real water_ang( 0.0 );
                //store water internal bond angle
                //WARNING: hard-coded for TP3 water (H1-O-H2)
     //            if( wat_is_acc ){
     //              water_ang = pose.conformation().bond_angle( AtomID( 2, acc_pos ),
     //                  AtomID( 1, acc_pos  ),
     //                  AtomID( 3, acc_pos  ) );
     //            }
                pose.conformation().set_bond_angle( AtomID( batm, acc_pos ), AtomID( aatm, acc_pos  ), AtomID( hatm, don_pos  ), BAHang );
    
                pose.conformation().set_bond_angle( AtomID( aatm, acc_pos ), AtomID( hatm, don_pos  ), AtomID( datm, don_pos  ), AHDang );
    
                pose.conformation().set_torsion_angle( AtomID( batm, acc_pos ), AtomID( aatm, acc_pos  ), AtomID( hatm, don_pos  ), AtomID( datm, don_pos  ), BAHDchi );
    
                if( wat_is_acc ){
                  pose.conformation().set_torsion_angle( AtomID( dbatm, don_pos ), AtomID( datm, don_pos  ), AtomID( hatm, don_pos  ), AtomID( aatm, acc_pos  ), B2BAHchi );
                }
                else{
                  pose.conformation().set_torsion_angle( AtomID( b2atm, acc_pos ), AtomID( batm, acc_pos  ), AtomID( aatm, acc_pos  ), AtomID( hatm, don_pos  ), B2BAHchi );
                }
    
                pose.conformation().set_bond_length( AtomID( aatm, acc_pos ), AtomID( hatm, don_pos  ), AHdist );
    
                //do fast clash check, OH hbonds are only 0.8A!
                if( fast_clash_check( pose, clash_check_atids, 0.8 ) ) continue;
    
                pose.fold_tree( f_jump );
                scorefxn->score( pose );
    
                //get score
                Real wat_score( pose.energies().residue_total_energies( new_seqpos ).dot( scorefxn->weights() ) );
    
		//TR << "gly_score: " << wat_score << std::endl;
                if( wat_score <= totalE ){
                  ++water_hb_states_good;
                  //TR_unsat << "AHdist: " << AHdist << " ";
                  //TR_unsat << "\tBAHang: " << BAHang << " ";
                  //TR_unsat << "\tBAHDchi: " << BAHDchi << " ";
                  //TR_unsat << "\twat_score: " << wat_score << std::endl;
                  //TR_unsat << pose.energies().total_energies().weighted_to_string( scorefxn->weights() )
                  //  + " total_score: " + to_string( pose.energies().total_energies()[ total_score ] );
                  //pose.dump_pdb( "hbond." + to_string( seqpos ) + "." + to_string( atomno ) + "." + to_string(hb_energy) +"."+to_string(wat_score)+ "."+to_string(AHdist)+"." + to_string( Size( numeric::conversions::degrees( AHDang ) ) ) + "." + to_string( Size( numeric::conversions::degrees( BAHang ) ) ) + "." + to_string( Size( numeric::conversions::degrees( BAHDchi ) ) ) + "." + to_string( Size( numeric::conversions::degrees( B2BAHchi ) ) ) + ".pdb" );
		  ref_pose=pose;
		  ref_pose.conformation().delete_residue_range_slow(1,ref_pose.total_residue()-1);
                  ref_pose.dump_pdb( atomname + "_" + to_string( seqpos ) + "_" + to_string( atomno ) + "_" + to_string(hb_energy) +"_"+to_string(wat_score)+ ".pdb");
                }
              }
            }
          }
        }
      }
	TR << "hb_states good: " << water_hb_states_good << std::endl;
	TR << "hb_states tot: " << water_hb_states_tot << std::endl;

	TR << "Finished generate_hbond_geometry" << std::endl;
	}//end of apply

	virtual std::string get_name() const {
		return "generate_hbond_geometry";
	}
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;

	SequenceMoverOP seq( new SequenceMover() );
	seq->add_mover( new generate_hbond_geometry() );

	try{
		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
	NEW_OPT(generate_hbond_geometry::resid, "Residue name",0);
	NEW_OPT(generate_hbond_geometry::atomname, "Atom name","O");
	NEW_OPT(generate_hbond_geometry::totalE, "total energy for glycine",10.0);
	NEW_OPT(generate_hbond_geometry::hbondE, "hydrogen bond cutoff",-1.0);

	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
}
