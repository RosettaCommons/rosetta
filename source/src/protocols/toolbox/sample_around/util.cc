// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/sample_around/util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu

////////////////////////////////////////////////////
//
// Developed and used by DasLab -- see
//
//  apps/fcchou/adenosine_sample_around.cc
//  apps/rhiju/phosphate_sample_around.cc
//
////////////////////////////////////////////////////

#include <protocols/toolbox/sample_around/util.hh>
#include <protocols/toolbox/rigid_body/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/rna/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <utility/io/ozstream.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/init/init.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/viewer/viewers.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/sample_around.OptionKeys.gen.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <basic/Tracer.hh>

using namespace core;
using namespace protocols;
using namespace numeric::conversions;
using utility::vector1;
typedef  numeric::xyzMatrix< Real > Matrix;
using namespace basic::options::OptionKeys;
using namespace basic::options;

static basic::Tracer TR( "protocols.toolbox.sample_around.util" );

/////////////////////////////////////////////////////////////////////////////
//
// Contains grab-bag of functions for sampling probe atoms around nucleobases
//  ('nucleobase_sample_around' app) for water and phosphate probes, uses euler
//   rotation utils.
//
//                           -- rhiju, 2014
//
/////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace toolbox {
namespace sample_around {

/////////////////////////////////////////////////////////////////////////////
//FCC: Adding Virtual res
void
add_virtual_res ( core::pose::Pose & pose, bool set_res_as_root /*= true */ ) {
	int nres = pose.total_residue();

	// if already rooted on virtual residue , return
	if ( pose.residue ( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		std::cout << "addVirtualResAsRoot() called but pose is already rooted on a VRT residue ... continuing." << std::endl;
		return;
	}

	// attach virt res there
	//	bool fullatom = pose.is_fullatom();
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type ( 1 ).residue_type_set();
	core::chemical::ResidueTypeCOP rsd_type( residue_set.get_representative_type_name3( "VRT" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type ) );
	pose.append_residue_by_jump ( *new_res , 1 );

	// make the virt atom the root
	if ( set_res_as_root ) {
		kinematics::FoldTree newF ( pose.fold_tree() );
		newF.reorder ( nres + 1 );
		pose.fold_tree ( newF );
	}
}

void
add_another_virtual_res ( core::pose::Pose & pose ) {
	//	int nres = pose.total_residue();
	// attach virt res there
	//	bool fullatom = pose.is_fullatom();
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type ( 1 ).residue_type_set();
	core::chemical::ResidueTypeCOP rsd_type( residue_set.get_representative_type_name3( "VRT" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type ) );
	pose.append_residue_by_jump ( *new_res , pose.total_residue() );
}


/////////////////////////////////////////////////////////////////////////////////////////////
// Rhiju -- rotate to my favorite frame. Base centroid is now at origin.
//         X points to N1 atom. Z points normal to base. Y is orthonormal and points towards Hoogsteen edge, I think.
void
rotate_into_nucleobase_frame( core::pose::Pose & pose ){

	using namespace core::conformation;
	using namespace core::chemical::rna;
	using namespace core::id;

	// assuming pose has an RNA at residue 1 -- will rotate just that residue.
	Size const base_pos( 1 );
	Residue const & rsd = pose.residue( base_pos );

	Vector centroid = get_rna_base_centroid( rsd, false /*verbose*/ );
	Matrix M = get_rna_base_coordinate_system( rsd, centroid );
	kinematics::Stub stub( M, centroid );

	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		for (Size i = 1; i <= pose.residue(n).natoms(); i++ ){
			Vector xyz_new = stub.global2local( pose.residue(n).xyz( i ) ); // it is either this or M-inverse.
			pose.set_xyz( AtomID( i, n ), xyz_new );
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////
// rhiju's crazy frame
void
rotate_into_phosphate_frame( core::pose::Pose & pose, Size const n, bool const center_on_OP2 ) {
	using namespace core::conformation;
	using namespace core::id;

	// define my own kind of coordinate system with phosphate plopped into center.
	Vector x,y,z; // point x from P to OP2, y towards OP1 (but perpendicular to OP2), z perpendicular.
	Residue const & rsd = pose.residue( n );
	kinematics::Stub stub;

	if ( center_on_OP2 ) {
		x = rsd.xyz( " OP2") - rsd.xyz( " P  " ); x.normalize();
		y = rsd.xyz( " OP1") - rsd.xyz( " P  " );
		z = cross( x, y ); z.normalize();
		y = cross( z, x ); y.normalize();
		stub = kinematics::Stub( Matrix::cols(x,y,z), rsd.xyz( " OP2" ) );
	} else {
		x = 0.5 * ( rsd.xyz( " OP2") + rsd.xyz( " OP1" ) ) - rsd.xyz( " P  " ); x.normalize();
		y = rsd.xyz( " OP1") - rsd.xyz( " P  " );
		z = cross( x, y ); z.normalize();
		y = cross( z, x ); y.normalize();
		stub = kinematics::Stub( Matrix::cols(x,y,z), rsd.xyz( " P  " ) );
	}

	for (Size i = 1; i <= rsd.natoms(); i++ ){
		Vector xyz_new = stub.global2local( pose.residue(n).xyz( i ) );
		pose.set_xyz( AtomID( i, n ), xyz_new );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//Measure the centroid distance between two bases
Real
centroid_dist( core::pose::Pose & pose,
							 bool const sample_another_adenosine_ ) {

	using namespace core::conformation;
	using namespace core::chemical::rna;
	using namespace core::id;

	Residue const & rsd1 = pose.residue( 1 );
	Residue const & rsd2 = sample_another_adenosine_ ? pose.residue( 4 ) : pose.residue( 3 );

	Vector centroid1, centroid2;
	if ( rsd1.is_RNA() ) {
		centroid1 = get_rna_base_centroid( rsd1, false /*verbose*/ );
	} else {
		centroid1 = rsd1.nbr_atom_xyz(); //Just use the nbr atom if not RNA
	}

	if ( rsd2.is_RNA() ) {
		centroid2 = get_rna_base_centroid( rsd2, false /*verbose*/ );
	} else {
		centroid2 = rsd2.nbr_atom_xyz(); //Just use the nbr atom if not RNA
	}

	return (centroid1 - centroid2).length();
}


///////////////////////////////////////////////////////////////////////////////////////////
// This is imported from protocols/stepwise/RigidBodySampler.cco
Real
sample_all_rotations_at_jump( pose::Pose & pose, Size const num_jump, scoring::ScoreFunctionOP scorefxn /* = 0 */ ){

	Real alpha_, alpha_min_( 0 ), alpha_max_( 180.0 ), alpha_increment_( option[ OptionKeys::sample_around::alpha_increment ]() );
	Real beta_, cosbeta_min_( -1.0 ), cosbeta_max_( 1.0 ), cosbeta_increment_( option[ OptionKeys::sample_around::cosbeta_increment ]()  );
	Real gamma_, gamma_min_( 0 ), gamma_max_( 180.0 ), gamma_increment_( option[ OptionKeys::sample_around::gamma_increment ]() );

	Matrix M;
	Vector axis1( 1.0, 0.0, 0.0 ), axis2( 0.0, 1.0, 0.0 ), axis3( 0.0, 0.0, 1.0 );

	Real const kT = 0.5;
	Real partition_function = 0;
	Size  count( 0 );
	Real  score_min( 0.0 );
	kinematics::Jump  best_jump;

	for ( alpha_ = alpha_min_; alpha_ <= alpha_max_;  alpha_ += alpha_increment_ ){

		//std::cout << i++ << " out of " << N_SAMPLE_ALPHA << ". Current count: " << count_total_ <<
		//			". num poses that pass cuts: " << count_good_ << std::endl;

		for ( Real cosbeta = cosbeta_min_; cosbeta <= cosbeta_max_;  cosbeta += cosbeta_increment_ ){
			if ( cosbeta < -1.0 ){
				beta_ = -1.0 * degrees( std::acos( -2.0 - cosbeta ) );
			} else if ( cosbeta > 1.0 ){
				beta_ = -1.0 * degrees( std::acos( 2.0 - cosbeta ) );
			} else {
				beta_ = degrees( std::acos( cosbeta ) );
			}

			//std::cout << "BETA: " << beta_ << std::endl;

			// Try to avoid singularity at pole.
			Real gamma_min_local = gamma_min_;
			Real gamma_max_local = gamma_max_;
			Real gamma_increment_local = gamma_increment_;
			if ( (beta_<-179.999 || beta_>179.999) ){
				gamma_min_local = 0.0;
				gamma_max_local = 0.0;
				gamma_increment_local = 1.0;
			}

			for ( gamma_ = gamma_min_local; gamma_ <= gamma_max_local;  gamma_ += gamma_increment_local ){

				protocols::toolbox::rigid_body::create_euler_rotation( M, alpha_, beta_, gamma_, axis1, axis2, axis3 );

				kinematics::Jump jump = pose.jump( num_jump );
				jump.set_rotation( M );
				pose.set_jump( num_jump, jump );

				if ( scorefxn ) {
					Real const score = (*scorefxn)( pose );
					partition_function += exp( - score / kT );
					if ( score < score_min || count == 0 ) {
						score_min = score;
						best_jump = jump;
					}
				} else {
					// this is a test
					pose.dump_pdb( "S_" + ObjexxFCL::string_of( count ) + ".pdb" );
				}

				count++;

			} // gamma
		} // beta
	}// alpha

	pose.set_jump( num_jump, best_jump );

	Real const free_E = - log( partition_function / count );

//	std::cout << "Energies: " << free_E << ' ' << score_min << std::endl;
//	return score_min;
	return free_E;

}


/////////////////////////////////////////////////////////////////////////////////
Real
do_scoring( pose::Pose & pose,
						scoring::ScoreFunctionOP scorefxn,
						bool const & sample_rotations,
						Size const probe_jump_num ){

	if ( sample_rotations ){
		return sample_all_rotations_at_jump( pose, probe_jump_num, scorefxn );
	}

	return (*scorefxn)( pose );

}

//////////////////////////////////////////////////////////
void
do_xy_scan( pose::Pose & pose,
						scoring::ScoreFunctionOP scorefxn,
						std::string const & outfile,
						Real const z,
						Size const probe_jump_num,
						Real const box_bins,
						Real const translation_increment,
						bool const sample_rotations ){

	kinematics::Jump jump = pose.jump( probe_jump_num );

	utility::io::ozstream out;
	out.open( outfile );

	Size count( 0 );
	Real best_score( 0.0 );
	Vector best_translation( 0.0, 0.0, 0.0 );

	for (int i = -box_bins; i <= box_bins; ++i) {
		for (int j = -box_bins; j <= box_bins; ++j) {
			Real const x = j * translation_increment;
			Real const y = i * translation_increment;
			jump.set_translation( Vector( x, y, z ) ) ;
			pose.set_jump( probe_jump_num, jump );
			Real score = do_scoring( pose, scorefxn, sample_rotations, probe_jump_num );
			out << score << ' ' ;

			if (score < best_score || count++ == 0 ){
				best_translation = Vector( x, y, z );
				best_score = score;
			}

		}
		out << std::endl;
	}
	out.close();

	// return pose in lowest energy configuration that was found in scan.
	jump.set_translation( best_translation );
	pose.set_jump( probe_jump_num, jump );
	do_scoring( pose, scorefxn, sample_rotations, probe_jump_num );
}

} //sample_around
} //toolbox
} //protocols
