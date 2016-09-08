#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh> // have to add it! otherwise pose.residue_type(i) won't function

#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <apps/pilot/rayyrw/rms_util.hh>


#ifndef apps_pilot_rayyrw_util_hh
#define apps_pilot_rayyrw_util_hh

static basic::Tracer tr_util("rayyrw_util");

// read frag_filename and parse the string out of it
core::Size
get_pos_from_fragfn(
	utility::file::FileName const &fragfn
){
	//static basic::Tracer tr_get_pos_from_fragfn("get_pos_from_fragfn");
	if ( ! utility::file::file_exists( fragfn ) ) {
		utility_exit_with_message( "Unable to open file: " + fragfn.name() + '\n' );
	}

	// format after_rotation.mer.pos.picker_rank.shd_rank.pdb
	utility::vector1< std::string > fragfn_tag = utility::string_split( fragfn.base(), '.' );
	//tr_get_pos_from_fragfn << "fragfn" << fragfn << " read in." << std::endl;
	if ( fragfn_tag.size() < 3 )  return 0;
	core::Size pos = utility::string2int( fragfn_tag[3] );
	//tr_get_pos_from_fragfn << "fragfn" << pos << " read in." << std::endl;
	return pos;
}


// for cal_overlap_scores and cal_nonoverlap_scores apps to prevent intense IO
// output vectors: frag_poses,
//                 frag_pos
// the index of these vectors corresponds to that of pdb_files
void
read_pdbs(
	utility::vector1<utility::file::FileName> const &frag_files, //input
	utility::vector1<core::pose::Pose> &frag_poses,
	utility::vector1<core::Size> &frag_positions
){
	core::chemical::ResidueTypeSetCOP cen_rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );

	int counter=0;
	for ( core::Size i=1; i<=frag_files.size(); ++i ) {
		//  frag_poses
		core::pose::Pose frag_pose;
		core::import_pose::pose_from_file( frag_pose, *cen_rsd_set, frag_files[i] , core::import_pose::PDB_file);
		frag_poses[i] = frag_pose;

		core::Size pos = get_pos_from_fragfn( frag_files[i] );
		assert ( pos != 0 );
		frag_positions[i] = pos;

		counter++;
	}
	std::cout << "read in " << counter << " frags " << std::endl;
}

// for cal_overlap_scores and cal_nonoverlap_scores apps to constructing poses to often
// from read_pdbs()
// The index of fragpdb_files should match that in fragpdb_poses
void
get_frag_rmsd(
	utility::vector1<core::Real> &frag_rmsds,
	utility::vector1<core::Real> const &frag_positions,
	utility::vector1<core::pose::Pose> const &frag_poses,
	core::pose::Pose const &native_pose
){
	int counter=0;
	for ( core::Size i=1; i<=frag_positions.size(); ++i ) {
		frag_rmsds[i] = native_frag_CA_RMSD( native_pose, frag_poses[i], frag_positions[i] ); // return rmsd and assign it to frag_rmsd_vector
		counter++;
	}
	std::cout << "done calculate rmsd for " << counter << " frags " << std::endl;
}


std::string
int2str(
	int n
){
	std::stringstream Num;
	std::string str;
	Num << n;
	str = Num.str();
	return str;
}


// you get negative value when there two frags are overlapping
int
get_gap_size(
	core::pose::Pose const &pose1,
	core::Size const pos1,
	core::pose::Pose const &pose2,
	core::Size const pos2
){
	core::pose::Pose const &pose_lower = ( pos1 > pos2 ) ? pose2 : pose1;
	//core::pose::Pose const &pose_upper = ( pos1 > pos2 ) ? pose1 : pose2;
	core::Size pos_lower = std::min( pos1, pos2 );
	core::Size pos_upper = std::max( pos1, pos2 );

	return pos_upper - ( pos_lower + pose_lower.size() - 1 ); // gap_size
}


//  moved remove virtual residues to somewhere

// This function can only be called when you know pos1 > pos2;
// The rule is - pos2 is always further than pos1
// pass by reference to swap poses directly
void
swap_poses(
	core::pose::Pose & pose1,
	core::Size pos1,
	core::pose::Pose & pose2,
	core::Size pos2
){
	static basic::Tracer tr_swap_poses("swap_poses");
	runtime_assert ( pos1 > pos2 ); // "=" mean 1 overlap

	tr_swap_poses << "poses got swapped by reference" << std::endl;
	core::pose::Pose temp_pose = pose1;
	core::Size temp_pos = pos1;
	pose1 = pose2;
	pos1  = pos2;
	pose2 = temp_pose;
	pos2  = temp_pos;
}


#endif
