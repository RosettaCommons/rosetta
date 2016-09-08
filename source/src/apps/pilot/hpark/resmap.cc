#include <core/pose/PDBInfo.hh>

void
get_resmap( pose::Pose const &pose,
	    pose::Pose const &ref_pose,
	    std::map< Size, Size > &resmap,
	    std::map< Size, Size > &pose_resmap
	    )
{
  for ( Size ii = 1; ii <= pose.size(); ++ii ) {
    Size ii_pdb( pose.pdb_info()->number( ii ) );
    pose_resmap[ii_pdb] = ii;
    
    for ( Size jj = 1; jj <= ref_pose.size(); ++jj ) {
      Size jj_pdb( ref_pose.pdb_info()->number( jj ) );
      
      if( ii_pdb == jj_pdb ){
	id::AtomID id1( pose.residue(ii).atom_index( "CA" ), ii );
	id::AtomID id2( ref_pose.residue(jj).atom_index( "CA" ), jj );
	
	resmap[ii] = jj;
	std::cout << "Map: " << ii << " " << ii_pdb << " mapped to ";
	std::cout << jj << " " << jj_pdb << std::endl;
	break;
      }
    }
    
  }
}

