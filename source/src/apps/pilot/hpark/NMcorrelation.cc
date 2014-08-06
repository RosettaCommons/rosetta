// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/normalmode/NormalMode.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/scoring/rms_util.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/types.hh>
#include <devel/init.hh>
#include <sys/time.h>

OPT_1GRP_KEY(IntegerVector, fpd, modes)

namespace myspace{

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

void
CArmsd_and_superimpose_pose(
														pose::Pose & pose,
														pose::Pose & ref_pose,
														std::map< Size, Size > const resmap,
														Real &gdtmm,
														Real &rmsd
														)
{
  //scoring::superimpose_pose( pose, ref_pose, atom_map );
	rmsd = scoring::CA_rmsd( pose, ref_pose, resmap );
	gdtmm = scoring::CA_gdtmm( pose, ref_pose, resmap );
}

void
get_resmap( pose::Pose const &pose,
						pose::Pose const &ref_pose,
						std::map< Size, Size > &resmap,
						std::map< Size, Size > &pose_resmap
						)
{
  for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
 		Size ii_pdb( pose.pdb_info()->number( ii ) );
		pose_resmap[ii_pdb] = ii;

  	for ( Size jj = 1; jj <= ref_pose.total_residue(); ++jj ) {
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

void
find_best_projection( pose::Pose const &pose, pose::Pose const &pose_ref,
											std::map< Size, Size > const resmap,
											std::map< Size, Size > const pose_resmap,
											utility::vector1< utility::vector1< Vector > > const vecs,
											utility::vector1< Real > const importance,
											utility::vector1< Real > &proj_scale,
											utility::vector1< Vector > &proj_vec
											)
{
	// input "vec" should be normalized!
	// Both pose should be aligned to each other!
	//utility::vector1< Vector > diffvec;
	std::map< Size, Vector > diffvec;
	proj_vec.resize( pose.total_residue() );

	proj_scale.resize( vecs.size() );

	for( Size ires = 1; ires <= pose.total_residue(); ++ires ){
		proj_vec[ires] = Vector( 0.0, 0.0, 0.0 );
	}

	Real diffsum( 0.0 );
	//for( Size ires = 1; ires <= pose.total_residue(); ++ires ){

	std::map< Size, Size >::const_iterator it;
	std::map< Size, Vector >::iterator it2;

	for( it = resmap.begin(); it != resmap.end(); ++it ){
		Size ires( it->first );
		Size jres( it->second );

		assert( pose.residue(ires).type() == pose_ref.residue(jres) );

		Size i_ca = pose.residue( ires ).atom_index( " CA " );
		Vector dxyz( pose.residue( ires ).xyz( i_ca ) - pose_ref.residue( jres ).xyz( i_ca ) );
		//std::cout << "dxyz : " << ires << " " << std::sqrt(dxyz.dot_product( dxyz )) << std::endl;
		diffvec[ires] = dxyz;
		diffsum += dxyz.dot_product( dxyz );
	}
	diffsum = std::sqrt(diffsum);

	// Dot_product
	for( Size i_mode = 1; i_mode <= vecs.size(); ++i_mode ){
		Real dotsum( 0.0 );
		Real dotsum_n( 0.0 );
		utility::vector1< Vector > const &vec( vecs[i_mode] );
		
		for( it2 = diffvec.begin(); it2 != diffvec.end(); ++it2 ){
			Size resno( it2->first );
			Vector diffvec_i( it2->second );
			dotsum += diffvec_i.dot_product( vec[resno] );
		}

		proj_scale[i_mode] = dotsum;
		std::cout << "Best scale/Correlation for " << std::setw(3) << i_mode;
		std::cout << ", importance: " <<  std::setw(10) << importance[i_mode];
		std::cout << ": " << std::setw(10) << dotsum;
		std::cout << " " << std::setw(10) << dotsum/diffsum;
		std::cout << std::endl;

		for( Size i_res = 1; i_res <= pose.total_residue(); ++i_res ){
			proj_vec[i_res] += dotsum*vec[i_res];
			/*
			std::cout << "vecdiff: " << i_ca;
			printf(" %8.5f", diffvec[i_ca][0]);
			printf(" %8.5f", diffvec[i_ca][1]);
			printf(" %8.5f", diffvec[i_ca][2]);
			printf(" %8.5f", proj_vec[i_ca][0]);
			printf(" %8.5f", proj_vec[i_ca][1]);
			printf(" %8.5f", proj_vec[i_ca][2]);
			std::cout << std::endl;
			*/
		}
	}

	// Finally report for merged
	Real projsum( 0.0 );
	Real corr(0.0);

	//for( Size i_ca = 1; i_ca <= proj_vec.size(); ++i_ca ){
	for( it2 = diffvec.begin(); it2 != diffvec.end(); ++it2 ){
		Size resno( it2->first );
		projsum += proj_vec[resno].dot_product( proj_vec[resno] );
		corr += proj_vec[resno].dot_product( diffvec[resno] );
	}
	projsum = std::sqrt(projsum);

	corr /= (projsum*diffsum);
	std::cout << "Best scale/Correlation for merged mode: " << corr << std::endl;
}

utility::vector1< Vector >
project( Real const scale, 
				 utility::vector1< Real > const proj_scale,
				 utility::vector1< utility::vector1< Vector > > const vecs
				 )
{
	utility::vector1< Vector > proj_vec;
	proj_vec.resize( vecs[1].size() );

	for( Size ires = 1; ires <= vecs[1].size(); ++ires ){
		proj_vec[ires] = Vector( 0.0, 0.0, 0.0 );
	}

	for( Size i_mode = 1; i_mode <= vecs.size(); ++i_mode ){
		utility::vector1< Vector > const &vec( vecs[i_mode] );
		for( Size i_ca = 1; i_ca <= vec.size(); ++i_ca ){
			proj_vec[i_ca] += scale*proj_scale[i_mode]*vec[i_ca];
		}
	}

	return proj_vec;
}

pose::Pose
generate_proj_pose( pose::Pose const &pose_init,
										utility::vector1< Vector > proj_vec
										)
{
	pose::Pose pose_proj( pose_init );
	// This generates unreasonable structure but just do it simple
	for( Size ires = 1; ires <= pose_init.total_residue(); ++ires ){
		Vector newvec = proj_vec[ires]+pose_init.residue(ires).xyz( " CA " );

		id::AtomID CaID( pose_init.residue(ires).atom_index( " CA " ), ires );
		pose_proj.set_xyz( CaID, newvec );
	}
	return pose_proj;
}
}

int main( int argc, char *argv [] ){
	using namespace myspace;
  using namespace protocols::normalmode;

	NEW_OPT( fpd::modes, "modes", 0 );

	devel::init(argc, argv);

  core::chemical::ResidueTypeSetCAP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

  core::pose::Pose pose, pose_ref;

  core::import_pose::pose_from_pdb( pose_ref, *rsd_set, option[ in::file::s ](1) ); 
  core::import_pose::pose_from_pdb( pose, *rsd_set, option[ in::file::s ](2) ); 

	//timeval time0, time1;
	//gettimeofday(&time0, NULL );

	std::map< Size, Size > resmap, pose_resmap;
	get_resmap( pose, pose_ref, resmap, pose_resmap );

	NormalMode NM( "CA", 15.0 );
	NM.solve( pose );
	//NM.set_harmonic_constants( 10.0, 1.0, 1.0 );

	// First calculate RMSD and align each other
	Real gdtmm0, rmsd0;
	Real gdtmm, rmsd;
	utility::vector1< Vector > proj_vec;
	utility::vector1< Real > proj_scale;
	utility::vector1< utility::vector1< Vector > > vecs;
	utility::vector1< Real > importance;

	// Use modes defined in input arg
	if( option[ fpd::modes ].user() ){
		Size nmodes = option[ fpd::modes ]().size();
		vecs.resize( nmodes );
		importance.resize( nmodes );

		for( Size i_mode = 1; i_mode <= nmodes; ++i_mode ){
			Size modeno( option[ fpd::modes ]()[i_mode] );
			vecs[i_mode] =  NM.get_eigvec_cart( modeno );
			importance[i_mode] = NM.get_importance( modeno );
		}
		
		CArmsd_and_superimpose_pose( pose, pose_ref, resmap, gdtmm0, rmsd0 );

		// Get projection scale and vector
		find_best_projection( pose, pose_ref, resmap, pose_resmap, 
													vecs, importance,
													proj_scale, proj_vec );

	} else {
		Size nmodes( 20 );
		vecs.resize( nmodes );
		importance.resize( nmodes );

		// For align
		CArmsd_and_superimpose_pose( pose, pose_ref, resmap, gdtmm0, rmsd0 );

		std::cout << "Merging top " << nmodes << " modes..." << std::endl;
		vecs.resize( nmodes );
		for( Size i_mode = 1; i_mode <= nmodes; ++i_mode ){
			vecs[i_mode] =  NM.get_eigvec_cart( i_mode );
			importance[i_mode] = NM.get_importance( i_mode );
		}
		find_best_projection( pose, pose_ref, resmap, pose_resmap, 
													vecs, importance,
													proj_scale, proj_vec );

	}

	/*
	std::cout << "Before(rmsd/gdtmm): " << rmsd0 << " " << gdtmm0 << std::endl;

	// Scanning over proj_scale, report RMSD, GDT expected by best projection
	for( Size k = 1; k <= 4; ++k ){
		Real scale = 0.25*(Real)(k);

		utility::vector1< Vector > proj_vec_tmp;
		proj_vec_tmp = project( scale, proj_scale, vecs );

		std::cout << "At scale " << scale << std::endl;
		pose::Pose pose_proj = generate_proj_pose( pose, proj_vec_tmp );
		CArmsd_and_superimpose_pose( pose_proj, pose_ref, gdtmm, rmsd );

		std::stringstream pdbname;
		pdbname << option[ in::file::s ](2) << ".proj" << scale << ".pdb";
		pose_proj.dump_pdb( pdbname.str() );
		std::cout << "After (rmsd/gdtmm): " << rmsd << " " << gdtmm << std::endl;
	}
	*/

	return 0;
}

