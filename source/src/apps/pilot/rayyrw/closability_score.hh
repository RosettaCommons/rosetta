// 130425
// using loophash to determine whether two poses are closable or not given a residue separation
// WARNING: loophash doesn't know how to deal with terminus
//
//
#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

#include <core/id/NamedAtomID.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/conformation/Residue.hh>

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

// Real9
#include <numeric/HomogeneousTransform.hh>

#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>

#include <ObjexxFCL/format.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/owning_ptr.hh>

// data structure
#include <boost/unordered_map.hpp>

#include <basic/Tracer.hh>

#include <apps/pilot/rayyrw/rms_util.hh> // cal_distance

#ifndef apps_pilot_rayyrw_closability_score_hh
#define apps_pilot_rayyrw_closability_score_hh



//class ClosabilityScore;
//typedef utility::pointer::owning_ptr< ClosabilityScore > ClosabilityScoreOP;
//typedef utility::pointer::owning_ptr< ClosabilityScore const > ClosabilityScoreCOP;

class ClosabilityScore {
public:

    // constructor
    ClosabilityScore(){
        db_to_use_threshold_ = 25;
    }


    // using loop_hash
    // return counts and gap_size
    core::Real 
    closability_score( core::pose::Pose const &pose1, 
                        core::Size const pos1, 
                        core::pose::Pose const &pose2, 
                        core::Size const pos2,
                        core::Size &counts, 
                        int &gap_size, 
                        core::Size radius=2 );


    // using simply distance check
    // return gap_size
    core::Real 
    closability_score( core::pose::Pose const &pose1, 
                        core::Size const pos1, 
                        core::pose::Pose const &pose2, 
                        core::Size const pos2,
                        int &gap_size,
                        bool stringent=false );


    // if we want to control via command line
    void
    set_db_to_use_threshold( core::Size threshold );


    // get rigid body transform from a given pose, and i j residues
    bool 
    get_rt_over_leap( const core::pose::Pose & orig_pose, 
                    core::Size const ir, 
                    core::Size const jr, 
                    numeric::geometry::hashing::Real6 & rt_6 );


    void 
    miniPose_creator( core::pose::Pose &miniPose,
                    core::pose::Pose pose1,
                    core::pose::Pose pose2,
                    core::Size pos2_start=1 );


    core::Real
    gap_bound_cst_score( core::pose::Pose const &pose, 
                core::Size const gap_size,
                core::Real stdev=1,
                core::Real cst_weight=1 );



    core::Real 
    gap_upper_boundary( core::Size seq_gap );


    core::Size 
    loophash_query( core::pose::Pose const & miniPose, 
                    core::Size const db_to_use, 
                    core::Size const radius );


private:
  boost::unordered_map<core::Size, protocols::loophash::LoopHashLibraryOP> LH_Lib_Map_;
  core::Size db_to_use_threshold_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static basic::Tracer tr("closability_measurement_class");


void 
ClosabilityScore::
set_db_to_use_threshold( 
    core::Size threshold
){
    db_to_use_threshold_ = threshold;
}


// this function was taken directly from LoopHashMap
bool 
ClosabilityScore::
get_rt_over_leap( 
    const core::pose::Pose& orig_pose, 
    core::Size const ir, 
    core::Size const jr, 
    numeric::geometry::hashing::Real6 &rt_6 
){
	using namespace core;
	using namespace core::pose;
	using namespace conformation;
	using namespace kinematics;
	using namespace numeric::geometry::hashing;

	core::pose::Pose pose = orig_pose;

	//fpd vrt/ligand trim
	core::Size nres = pose.total_residue();
	while (!pose.residue_type(nres).is_polymer()) nres--;

	// get current cutpoints; don't try to connect these
	utility::vector1< core::Size > cuts_in = pose.fold_tree().cutpoints();
	std::sort( cuts_in.begin(), cuts_in.end() );

  // ryw commented this out because my Minipose only have i and j
  /*
	// bail if (ir,jr) crosses a cut
	for (Size i=1; i<=cuts_in.size(); ++i) {
		if (cuts_in[i]<=jr && cuts_in[i]>=ir) {
			TR.Error << "ERROR -- residue range crosses cut    IR: " << ir << "  JR: " << jr << "  CUT: " << cuts_in[i] << std::endl;
			return false;
		}
		//fpd insertions one position after the cut seem not to work ...
		//fpd perhaps if the foldtree for the local segment were reversed this might be ok
		if (cuts_in[i]==ir-1) {
			TR.Error << "ERROR -- startres immediately follows cut    IR: " << ir << "  CUT: " << cuts_in[i] << std::endl;
			return false;
		}
	}
  */

	// Create a fake foldtree with a jump from ir to jr, and the cutpoint just before jr. From that extract the
	// rigid body transfer. This is fairly hacky, but actually works v reliably and is easy to understand, without
	// having to mess deeply with the fold tree.
	// FoldTree f;
	// Size cutpoint= jr-1;
	// f.add_edge( 1, ir, Edge::PEPTIDE );
	// f.add_edge( ir, cutpoint, Edge::PEPTIDE );
	// f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
	// f.add_edge( jr, nres , Edge::PEPTIDE );
	// f.add_edge( ir, jr, 1 );  // this is the jump !!

	//fpd handle multiple chains/chainbreaks
	FoldTree f;
	core::Size last_cut=0, jump_num=1;
	Size cutpoint= jr-1;
	for (Size i=1; i<=cuts_in.size(); ++i) {
		if (cuts_in[i] >= nres) break;
		if (cutpoint > last_cut && cutpoint < cuts_in[i]) {
			f.add_edge( last_cut+1, ir, Edge::PEPTIDE );
			f.add_edge( ir, cutpoint, Edge::PEPTIDE );
			f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
			f.add_edge( jr, cuts_in[i] , Edge::PEPTIDE );
			f.add_edge( ir, jr, 1 );  // this is the jump !!
			if (last_cut!=0) f.add_edge( 1, last_cut+1, jump_num++);
		} else {
			if (last_cut+1 != cuts_in[i]) f.add_edge( last_cut+1, cuts_in[i], Edge::PEPTIDE );
			if (last_cut!=0) f.add_edge( 1, last_cut+1, jump_num++);
		}
		last_cut = cuts_in[i];
	}
	if (last_cut+1 <= nres) {
		if (cutpoint > last_cut && cutpoint < nres) {
			f.add_edge( last_cut+1, ir, Edge::PEPTIDE );
			f.add_edge( ir, cutpoint, Edge::PEPTIDE );
			f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
			f.add_edge( jr, nres , Edge::PEPTIDE );
			f.add_edge( ir, jr, 1 );  // this is the jump !!
			if (last_cut!=0) f.add_edge( 1, last_cut+1, jump_num++);
		} else {
			f.add_edge( last_cut+1, nres, Edge::PEPTIDE );
			if (last_cut!=0) f.add_edge( 1, last_cut+1, jump_num++);
		}
	}
	for (core::Size i=nres+1; i<=pose.total_residue(); ++i)
		f.add_edge( 1, i, jump_num++ );  // additional jumps

	core::Size theroot = 1;
	if( ir == 1 ) theroot = pose.total_residue();
	if( orig_pose.residue_type( orig_pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) theroot = orig_pose.fold_tree().root();  //fpd
	if( f.reorder(theroot) == false ){
		tr.Error << "ERROR During reordering of fold tree - am ignoring this LOOP ! bailing: The root: " << theroot << " NRES " << pose.total_residue() << "   IR: " << ir << "  JR: " << jr << std::endl;
		return false; // continuing leads to a segfault - instead ignore this loop !
	}

	// Apply this new foldtree to the pose.
	pose.fold_tree( f );

	// Now extract the rigid body transform!
	Jump myjump;
	myjump = pose.jump( 1 );

	// Aha, now you have the RT (RigidbodyTransform)
	RT rt = myjump.rt();

	// Create a 6 value representation: (just change the data format)
	numeric::HomogeneousTransform< Real > ht( rt.get_rotation() , rt.get_translation() );
	numeric::xyzVector < Real > euler_angles =  ht.euler_angles_rad();

	rt_6[1] = rt.get_translation().x();
	rt_6[2] = rt.get_translation().y();
	rt_6[3] = rt.get_translation().z();
	rt_6[4] = euler_angles.x()*180.0/numeric::constants::d::pi;
	rt_6[5] = euler_angles.y()*180.0/numeric::constants::d::pi;
	rt_6[6] = euler_angles.z()*180.0/numeric::constants::d::pi;

	// indicate success
	return true;
}

// this function is copied from /work/wangyr/rosetta/src/protocols/hybridization/FoldTreeHybridize.cc
// the number is provided as upper bound for bound constraints
// gap definition: 5 | 6 -> gap=1
// supposedly, we only allow
core::Real
ClosabilityScore::gap_upper_boundary( 
  core::Size seq_gap
){
  // delta = distribution I observed 
  // add 2 to the original upper bound
  core::Real gap_torr_1(6.0); //gap_torr_1( 5.0);  delta=~2.5A // used to be 4.0, but give it 5
  core::Real gap_torr_2(9.6); //gap_torr_2( 8.5);  delta=~4A
  core::Real gap_torr_3(13.0); //gap_torr_3(12.0);  delta=~7A
  core::Real gap_torr_4(16.0); //gap_torr_4(14.5);
  core::Real gap_torr_5(19.5); //gap_torr_5(18.0);
  core::Real gap_torr_6(22.0); //gap_torr_6(21.0);
  core::Real gap_torr_7(26.0); //gap_torr_7(24.5);
  core::Real gap_torr_8(29.0); //gap_torr_8(27.5);
  core::Real gap_torr_9(33.0); //gap_torr_9(31.0);
  core::Real gap_torr_10(35.0); //gap_torr_9(33.0);
  core::Real gap_torr_11(38.0); //gap_torr_9(35.9);
  core::Real gap_torr_12(40.0); //gap_torr_9(35.9);
  core::Real gap_torr_13(43.0); //gap_torr_9(35.9);
  core::Real gap_torr_14(45.0); //gap_torr_9(35.9);
  core::Real gap_torr_15(47.0); //gap_torr_9(35.9);

  switch (seq_gap) {
    case 1:
      return gap_torr_1; break;
    case 2:
      return gap_torr_2; break;
    case 3:
      return gap_torr_3; break;
    case 4:
      return gap_torr_4; break;
    case 5:
      return gap_torr_5; break;
    case 6:
      return gap_torr_6; break;
    case 7:
      return gap_torr_7; break;
    case 8:
      return gap_torr_8; break;
    case 9:
      return gap_torr_9; break;
    case 10:
      return gap_torr_10; break;
    case 11:
      return gap_torr_11; break;
    case 12:
      return gap_torr_12; break;
    case 13:
      return gap_torr_13; break;
    case 14:
      return gap_torr_14; break;
    case 15:
      return gap_torr_15; break;
    default:
      return 1.5;
  }
  return 1.5;
}


// this is probably not necessary anymore
core::Real
ClosabilityScore::
gap_bound_cst_score(
    core::pose::Pose const &pose, 
    core::Size const gap_size,
    core::Real stdev,
    core::Real cst_weight
){
  core::pose::Pose temp_pose;
  temp_pose = pose;
  core::scoring::func::FuncOP fx( new core::scoring::constraints::BoundFunc( 1.0, gap_upper_boundary( gap_size ), stdev, "gap" ) );
  temp_pose.add_constraint(
    core::scoring::constraints::ConstraintCOP( 
      core::scoring::constraints::ConstraintOP( 
        new core::scoring::constraints::AtomPairConstraint( 
          core::id::AtomID( temp_pose.residue(1).atom_index(" CA "), 1 ),
          core::id::AtomID( temp_pose.residue(2).atom_index(" CA "), 2 ),
          fx 
        )
      )
    )
  );
  tr << "energy function created" << std::endl;
  core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction() );
  scorefxn->set_weight( core::scoring::atom_pair_constraint, cst_weight );
  core::Real cst_score = (*scorefxn)(temp_pose);

  return cst_score;
}


// given two frag poses, the program takes ends residues and makes a two-residue pose out of that.
void
ClosabilityScore::
miniPose_creator(
    core::pose::Pose &miniPose,
    core::pose::Pose pose1,
    core::pose::Pose pose2,
    core::Size pos2_start
){
  // make sure a clean input here
  miniPose.clear(); 

  // create A mini pose out of two poses - poseA numbering must be further than poseB 
  // pose1 --------upper  lower--------pose2
  // pose1
  core::pose::remove_upper_terminus_type_from_pose_residue( pose1, pose1.n_residue() );
  core::pose::add_lower_terminus_type_to_pose_residue( pose1, pose1.n_residue() );
  core::conformation::ResidueCOP rsd1 = pose1.residue( pose1.n_residue() ).clone(); // copy residue from poseA

  // pose2
  core::pose::remove_lower_terminus_type_from_pose_residue( pose2, pos2_start );
  core::pose::add_upper_terminus_type_to_pose_residue( pose2, pos2_start );
  core::conformation::ResidueCOP rsd2 = pose2.residue( pos2_start ).clone();

  miniPose.append_residue_by_bond( *rsd1 );
  miniPose.append_residue_by_jump( *rsd2, 1, "", "", false ); // 1=jump_anchor residue -> should be rsd1
  //miniPose.append_residue_by_bond( *rsd2 );
  //miniPose.dump_pdb( "miniPose.pdb" );
  tr << "miniPose has" << miniPose.n_residue() << " residues" << std::endl;
}
  
 

core::Size
ClosabilityScore::
loophash_query(
    core::pose::Pose const & miniPose,
    core::Size const db_to_use,
    core::Size const radius
){
    if( ! basic::options::option[ basic::options::OptionKeys::lh::db_path ].user() ){
        utility_exit_with_message("user must specify -lh:db_path for loop hash database");
    }

    // get library query set up
    utility::vector1 < core::Size > loop_sizes;  // vector
    loop_sizes.push_back( db_to_use );

    // check whether the library has existed in LH_Lib_Map_[ db_to_use ] = library;
    boost::unordered_map< core::Size, protocols::loophash::LoopHashLibraryOP >::const_iterator exists = LH_Lib_Map_.find( db_to_use );
    if ( exists == LH_Lib_Map_.end() ){
        tr << "loop_hash library: " << db_to_use << " has not existed yet!" << std::endl;
        protocols::loophash::LoopHashLibraryOP library( new protocols::loophash::LoopHashLibrary( loop_sizes ) );
        library->load_db();
        LH_Lib_Map_[ db_to_use ] = library;
    } else {
        tr << "loop_hash library: " << db_to_use << " existed!" << std::endl;
    }

    // current_library is designed to have only one hash_size at a time, 
    // therefore the loop below should have only 1 iteration
    protocols::loophash::LoopHashLibraryOP current_library = LH_Lib_Map_[ db_to_use ];
    std::vector < core::Size > leap_index_bucket;
    for( std::vector< core::Size >::const_iterator jt = current_library->hash_sizes().begin(); 
            jt != current_library->hash_sizes().end(); ++jt ){

        core::Size loop_size = *jt;
        protocols::loophash::LoopHashMap &hashmap = current_library->gethash( loop_size );

        tr << "getting a loop_transform from miniPose" << std::endl;
        numeric::geometry::hashing::Real6 loop_transform;
        get_rt_over_leap( miniPose, 1, 2, loop_transform ); // 1, 2 are the numbering in miniPose

        hashmap.radial_lookup( radius, loop_transform, leap_index_bucket ); // note! should probably make the radius 2 as a variable
    }
    return leap_index_bucket.size();
}



// gap_size <= 0; return 0
// loop_hash version of closability score
core::Real
ClosabilityScore::
closability_score(
    core::pose::Pose const &pose1, 
    core::Size const pos1, 
    core::pose::Pose const &pose2, 
    core::Size const pos2, 
    core::Size & counts, 
    int & gap_size,
    core::Size radius
){
    if ( pos1 == pos2 ){ 
        tr << pos1 << " equals to " << pos2 << std::endl;
        return 0;
    }
    core::pose::Pose const &pose_lower = ( pos1 > pos2 ) ? pose2 : pose1;
    //core::pose::Pose const &pose_upper = ( pos1 > pos2 ) ? pose1 : pose2;
    core::Size pos_lower = std::min( pos1, pos2 );
    core::Size pos_upper = std::max( pos1, pos2 );

    // gap_size
    // therefore pos2 is always larger than pos1
    gap_size = pos_upper - ( pos_lower + pose_lower.n_residue() - 1 );
    tr << "gap_size: " << gap_size << " = " << pos2 << " -  (" << pos1 << " + " << pose1.n_residue() << " - 1 ) " << std::endl;

    // db_to_use starts with 3, which means gap_size starts with 2
    core::Size db_to_use = gap_size+1;
    //core::Real score( 0.05 );

    if ( gap_size <= 0 ){ // should return a proper value later on
        tr << "There is no gap for pos1: " << pos1 << " and pos2: " << pos2 << std::endl;
        return 0;  
    } else if ( db_to_use > db_to_use_threshold_ ){  
        // loophash library only has 3-25, db_to_use = gap_size+1
        // from 25 to 15, it shows certain 
        tr << "There is no corresponding lib for this size of gap: " << pos1 << " and pos2: " << pos2 << std::endl;
        return 0;  
    }

    core::pose::Pose miniPose;
    core::Real bound_cst_score;

    if ( db_to_use == 2  ){ // gap=1
        //miniPose_creator( miniPose, pose1, pose2 ); 
        miniPose_creator( miniPose, pose1, pose2 ); // 2 means for pose2 starts with res2

        tr << "distane: " << cal_distance( miniPose, 1, 2 ) << std::endl;
        bound_cst_score = gap_bound_cst_score( miniPose, gap_size, 0.7, 2.0 ); 
        tr << "bound_cst_score: " << bound_cst_score << std::endl;

        tr << "create a miniPose that consists of the end stub from pose1, and 2nd stub from pose2" << std::endl;
        miniPose_creator( miniPose, pose1, pose2, 2 ); // 2 means for pose2 starts with res2

        counts = loophash_query( miniPose, 3, radius );// db_to_use change into 3;
        tr << "loophash_query found: " << counts << std::endl;

    } else {
        tr << "create a miniPose that consists of two stubs from two poses" << std::endl;
        miniPose_creator( miniPose, pose1, pose2 );

        tr << "distane: " << cal_distance( miniPose, 1, 2 ) << std::endl;
        bound_cst_score = gap_bound_cst_score( miniPose, gap_size, 0.7, 2.0 ); 
        tr << "bound_cst_score: " << bound_cst_score << std::endl;

        counts = loophash_query( miniPose, db_to_use, radius );
        tr << "loophash_query found: " << counts << std::endl;
    } 


    // loophash query performs poorly on short gaps (1 or 2)
    // the optimum situation where loophash and bound_cst_score agree with each other
    // if gap_size > 2; bound_cst_score will return a absurd number, which is from upper and lower boundary = 1.5
    // for gap_size <= 2; I would trust bound_cst_score
    //if ( bound_cst_score == 0 && counts == 0 ){
    if ( bound_cst_score == 0 && counts == 0 ){ // closable
        return -1.0;
    } else if ( counts > 0 ){ // closable
        return -1.0;
    } else if ( counts == 0 ){ // not closable
        return 1.0;
    } else {
      return 1.0;  // no closable
    }
    //score = -1*weight*std::log( counts+score );
    //return score;
}


// simple distance count closability score
core::Real
ClosabilityScore::
closability_score(
    core::pose::Pose const &pose1, 
    core::Size const pos1, 
    core::pose::Pose const &pose2, 
    core::Size const pos2, 
    int &gap_size,
    bool stringent
){
    if ( pos1 == pos2 ) // although this logic is covered later on, this will be good if we return it right away
        return 0; 

    core::pose::Pose const &pose_lower = ( pos1 > pos2 ) ? pose2 : pose1;
    core::pose::Pose const &pose_upper = ( pos1 > pos2 ) ? pose1 : pose2;
    core::Size pos_lower = std::min( pos1, pos2 );
    core::Size pos_upper = std::max( pos1, pos2 );

    // gap_size
    gap_size = pos_upper - ( pos_lower + pose_lower.n_residue() - 1 );
    //tr << "gap_size: " << gap_size << " = " << pos_upper << " -  (" << pos_lower << " + " << pose_lower.n_residue() << " - 1 ) " << std::endl;

    if ( gap_size <= 0 ){ // should return a proper value later on
        //tr << "There is no gap for pos1: " << pos1 << " and pos2: " << pos2 << std::endl;
        return 0;  
    } else if ( gap_size > 15 ){
        //tr << "There is no corresponding data for this size of gap: "<< gap_size << ", " << pos1 << " and pos2: " << pos2 << std::endl;
        return 0;  
    } else {
        //core::pose::Pose miniPose;
        //core::Real bound_cst_score;

        //tr << "create a miniPose that consists of the end stub from pose1, and 1st stub from pose2" << std::endl;
        //miniPose_creator( miniPose, pose1, pose2 ); // 2 means for pose2 starts with res2
        //tr << "distane: " << cal_distance( miniPose, 1, 2 ) << std::endl;

        //bound_cst_score = gap_bound_cst_score( miniPose, gap_size, 1, 1.0 ); 
        //tr << "bound_cst_score: " << bound_cst_score << std::endl;

        //if ( bound_cst_score == 0 ){ // closable 
        //    return -1.0;
        //} else { // non-closable; at most, return 1
        //    return std::min( 1.0, bound_cst_score );
        //} 

        numeric::xyzVector< core::Real > atm_i = pose_lower.residue( pose_lower.total_residue() ).atom(2).xyz();
        numeric::xyzVector< core::Real > atm_j = pose_upper.residue( 1 ).atom(2).xyz();
        core::Real dist_ij = ( atm_i - atm_j ).length();
        core::Real tgt_length = gap_upper_boundary( gap_size );

        if ( dist_ij < tgt_length ){ // closable
            return -1.0;
        } else{
            // (dist-upper_bound)**2, 
            // at most return 1 140220, I reweight this later on via scripting
            if ( stringent == true )
                return 1.0;
            else
                return std::min( 1.0, ( ( dist_ij - tgt_length )*( dist_ij - tgt_length ) )); 
        }

    } // 0 < gap_size <= 15
}
#endif
