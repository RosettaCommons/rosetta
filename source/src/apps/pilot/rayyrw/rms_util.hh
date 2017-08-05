#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>  //make_pose_from_sequence
#include <core/pose/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/rms_util.tmpl.hh>

#include <iostream>
#include <string>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>


#ifndef apps_pilot_rayyrw_rms_util_hh
#define apps_pilot_rayyrw_rms_util_hh

static THREAD_LOCAL basic::Tracer tr_rms_util("rms_util");

// calculate distance between two residues from two poses
core::Real
cal_distance(
	core::pose::Pose const &pose1,
	core::Size const rsd1,
	core::pose::Pose const &pose2,
	core::Size const rsd2
){
	core::Vector diff = pose1.residue( rsd1 ).xyz("CA") - pose2.residue( rsd2 ).xyz("CA");
	core::Real dist = std::sqrt( diff.length_squared() );

	return dist;
}


core::Real
cal_distance(
	core::pose::Pose const & pose,
	core::Size const resi,
	core::Size const resj
){
	//std::string const atom_name( "CA" );
	core::Vector diff = pose.residue( resi ).xyz("CA") - pose.residue( resj ).xyz("CA"); // 2 == "CA"

	return std::sqrt( diff.length_squared() );
}


void
remove_all_virtual_residues(
	core::pose::Pose & pose
){
	//assert(pose);
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue_type(i).name() == "VRT" ) {
			pose.conformation().delete_residue_slow(i);
		}
	}
}


// rmsd_no_super
core::Real
cal_rmsd(
	core::pose::Pose & frag_pose,
	core::pose::Pose & native_frag_pose
)
{
	remove_all_virtual_residues( frag_pose );
	return core::scoring::rmsd_no_super( frag_pose, native_frag_pose, core::scoring::is_protein_CA );
}


// how to calculate RMSD when two frags are different
core::Real
nosuper_CA_rmsd(
	core::pose::Pose pose1,
	core::Size pos_1,
	core::pose::Pose pose2,
	core::Size pos_2
)
{
	tr_rms_util << pos_1 << " " << pos_2 << std::endl;
	if ( pos_1 > pos_2 ) {
		core::pose::Pose temp_pose = pose1;
		core::Size temp_pos = pos_1;
		tr_rms_util << pos_1 << " " << pos_2 << std::endl;
		pose1 = pose2;
		pos_1  = pos_2;
		pose2 = temp_pose;
		pos_2  = temp_pos;
	} // therefore pos2 is always larger than pos1

	tr_rms_util << pos_1 << " " << pos_2 << std::endl;
	assert( pos_2 < pos_1 );

	core::Real sum( 0.0 );
	core::Size natoms( 0 );
	std::string const atom_name( "CA" );

	core::Size offset = pos_2 - pos_1;
	tr_rms_util << "offset: " << offset << std::endl;

	for ( core::Size irsd=1; irsd <= pose1.size(); ++irsd ) {
		tr_rms_util << pose1.residue( irsd ) << std::endl;
		if ( irsd+offset > pose1.size() ) { break; }
		if ( pose1.residue( irsd ).is_virtual_residue() ) { continue; }

		core::Vector diff = pose1.residue( irsd + offset ).xyz( atom_name ) - pose2.residue( irsd ).xyz( atom_name );
		tr_rms_util << "i=" << irsd << " pose1 resn: " << pos_1 + irsd + offset - 1 << ", pose2 resn: " << pos_2 + irsd - 1 << " dist: " << std::sqrt( diff.length_squared() ) << std::endl;

		sum += diff.length_squared();
		natoms += 1;
	}
	tr_rms_util << "natoms: " << natoms << std::endl;
	//return std::sqrt(sum / natoms);
	return sum / natoms;
}



core::Real
native_frag_CA_RMSD(
	core::pose::Pose const & native_pose,
	core::pose::Pose const & frag_pose,
	core::Size seq_pos
){
	// get native sequence
	std::string sequence;
	sequence = native_pose.sequence();

	// get frag length
	core::Size Nmers_size;
	Nmers_size = frag_pose.size();

	//tr_rms_util << "native_sequence:" << sequence << std::endl;
	//tr_rms_util << "frag_pose Nmers_size:" << Nmers_size << std::endl;
	//tr_rms_util << "seq_pos: " << seq_pos << std::endl;
	std::string frag_seq = sequence.substr( seq_pos-1, Nmers_size ); // get the subsequce from native sequence

	utility::vector1< core::Size > positions;
	core::kinematics::FoldTree fold_tree( Nmers_size );
	core::pose::Pose native_frag_pose;
	core::pose::make_pose_from_sequence( native_frag_pose, frag_seq, *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )) );

	for ( core::Size irsd=seq_pos; irsd<seq_pos+Nmers_size; irsd++ ) {
		positions.push_back( irsd );
	}
	core::pose::create_subpose( native_pose, positions, fold_tree, native_frag_pose );

	// calculate rmsd_no_super
	core::Real rmsd_to_native;
	rmsd_to_native = core::scoring::rmsd_no_super( frag_pose, native_frag_pose, core::scoring::is_protein_CA );

	return rmsd_to_native;
}

#endif
