// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/CoordConstraintClaimer.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/LegacyRootClaim.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/RigidChunkClaimer.hh>

// Project Headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <protocols/toolbox/superimpose.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// #include <core/kinematics/MoveMap.hh>
// #include <core/fragment/FragSet.hh>
// #include <protocols/simple_moves/FragmentMover.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/LoopsFileIO.hh>


// ObjexxFCL Headers

// Utility headers
#include <protocols/topology_broker/Exceptions.hh>

//#include <utility/io/izstream.hh>
#include <utility>
#include <utility/io/ozstream.hh> //fur dump cst
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>
//#include <basic/options/option.hh>
#include <numeric/random/random.hh>
//// C++ headers
#include <fstream>

#include <core/id/NamedStubID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <ObjexxFCL/FArray2D.hh>

// option key includes

static THREAD_LOCAL basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;
using namespace scoring::constraints;
using namespace ObjexxFCL;

CoordConstraintClaimer::CoordConstraintClaimer() :
	filename_( "NO_FILE"),
	constraints_( /* NULL */ ),
	root_( 1 ),
	bRegenerateFromInputPose_ ( false ),
	cst_pose_ ( /* NULL */ ),
	bCentroid_( true ),
	bFullatom_( false ),
	bLocal_( false )
{}

CoordConstraintClaimer::CoordConstraintClaimer( std::string filename ) :
	filename_(std::move( filename )),
	constraints_( /* NULL */ ),
	root_( 1 ),
	bRegenerateFromInputPose_ ( false ),
	cst_pose_ ( /* NULL */ ),
	bCentroid_( true ),
	bFullatom_( false ),
	bLocal_( false )
{}

CoordConstraintClaimer::~CoordConstraintClaimer() = default;

void CoordConstraintClaimer::new_decoy() {
	if ( !constraints_ && !bUseXYZ_in_cstfile_ ) {
		if ( !cst_pose_ ) throw EXCN_Input( "CoordConstraintClaimer::new_decoy(): in broker setup provide PDB_FILE or say USE_XYZ_FROM_CSTFILE");
		generate_constraints( *cst_pose_ );
	}
}

void CoordConstraintClaimer::new_decoy( core::pose::Pose const& pose ) {
	if ( bRegenerateFromInputPose_ ) {
		sequence_ = ""; //force new stealing
		cst_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );
	}
	if ( bSuperimpose_ ) {
		sequence_ = ""; //force new superposition
	}
	new_decoy();
}

void CoordConstraintClaimer::generate_claims( claims::DofClaims& new_claims ) {
	if ( bLocal_ ) {
		Size new_root( root_ );
		RigidChunkClaimer::CM_SuggestFixResidue msg( root_from_label_ );
		broker().relay_message( msg );
		tr.Info << label() << " got a returned message: " << msg << std::endl;

		if ( msg.received() ) new_root=msg.good_fix_pos_;
		else throw EXCN_Input( "no fixed region (looked for "+root_from_label_+") found to use as reference for CoordinateConstraints");

		if ( !constraints_ && bUseXYZ_in_cstfile_ ) {
			//in this case we haven't been able to read the constraints file yet, since this is the first time a valid pose exists...
			cst_pose_ = core::pose::PoseOP( new pose::Pose( broker().current_pose() ) );
			read_constraints_from_file( *cst_pose_ );
		}

		if ( new_root != root_ ) {
			root_ = new_root;
			tr.Info << "set new root for CoordinateConstraints to " << root_ << std::endl;
			set_cst_root();
		}
		runtime_assert( root_ != 0 );
	}
	if ( !bLocal_ ) new_claims.push_back( claims::DofClaimOP( new claims::LegacyRootClaim( get_self_weak_ptr(), root_, claims::DofClaim::NEED_TO_KNOW ) ) );
}

void CoordConstraintClaimer::read_cst_pose() {
	cst_pose_ = core::pose::PoseOP( new core::pose::Pose );
	core::import_pose::pose_from_file( *cst_pose_,
		*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ),
		filename_, core::import_pose::PDB_file);
	sequence_ = "";
}

void CoordConstraintClaimer::add_constraints( core::pose::Pose& pose ) const {
	// Size new_root;

	std::string const new_sequence ( pose.annotated_sequence( true ) );
	bool fullatom( pose.is_fullatom() );

	//filter --- constraints can be active/inactive for full/centroid poses
	if ( fullatom && !bFullatom_ ) return;
	if ( !fullatom && !bCentroid_ ) return;

	/// new sequence ? --- re-configure constraints ( atomnumbers! )
	if ( sequence_ != new_sequence ) {
		tr.Debug << "new sequence -- remap constraints " << std::endl;

		/// we should have a cst_pose_ now
		if ( !cst_pose_ ) {
			throw EXCN_Input( "CoordConstraintClaimer::add_constraints(): "
				"in broker setup either provide PDB_FILE or set CST_FROM_INPUT_POSE");
		}

		//constraints_ should be defined, since new_decoy()
		runtime_assert( constraints_ != nullptr );

		//get constraints
		if ( !bUseXYZ_in_cstfile_ ) {
			// get new constrains from cst_pose   --- this takes care of atom-renumbering
			constraints_ = constraints_->steal_def_clone( *cst_pose_, pose );
		} else { //need to remap constraints or read from file
			constraints_ = constraints_->remapped_clone( *cst_pose_, pose );
		}

		if ( bSuperimpose_ ) {
			superimpose( pose );
		}
		// debug output
		if ( tr.Trace.visible() ) {
			utility::io::ozstream dump_cst( label()+"_coord_claimer.cst" );
			constraints_->show_definition( dump_cst, pose );
			cst_pose_->dump_pdb("cst_pose.pdb");
		}
	}

	//finally -- add constraints to pose
	tr.Debug << "add constraints "<< label() << std::endl;
	pose.add_constraints( constraints_->get_all_constraints() );

	//keep sequence -- to watch for changes
	sequence_ = new_sequence;
}

void CoordConstraintClaimer::set_cst_root() {
	runtime_assert( cst_pose_ != nullptr );
	id::StubID cst_fix_stub_ID( core::pose::named_stub_id_to_stub_id( id::NamedStubID( "N","CA","C", root_ ), *cst_pose_ ));
	runtime_assert( constraints_ != nullptr );
	ConstraintCOPs all_cst = constraints_->get_all_constraints();
	ConstraintSetOP new_set( new ConstraintSet );
	for ( ConstraintCOPs::const_iterator it = all_cst.begin(), eit = all_cst.end(); it!=eit; ++it ) {
		ConstraintOP new_cst = (*it)->clone();
		LocalCoordinateConstraintOP ll_cst = utility::pointer::dynamic_pointer_cast< core::scoring::constraints::LocalCoordinateConstraint > ( new_cst );
		runtime_assert( ll_cst != nullptr ); //only these should be in the constraint set!
		ll_cst->set_fixed_stub( cst_fix_stub_ID );
		new_set->add_constraint( ll_cst );
	}
	constraints_=new_set;
}

void CoordConstraintClaimer::generate_constraints( pose::Pose const& cst_pose ) const {

	//empty constraint set
	constraints_ = core::scoring::constraints::ConstraintSetOP( new ConstraintSet );

	//make local copy of region-definition
	loops::Loops rigid( rigid_ );

	//if no regions defined ---> replace it with 1..nres
	if ( !rigid.size() ) rigid.add_loop( 1, cst_pose.size(), 0 );

	//go thru regions and generate constraints
	 for ( auto const & it : rigid ) {
		for ( Size pos = it.start(); pos <= it.stop(); ++pos ) {

			//generate constraint for ( CA, pos )
			conformation::Residue const & rsd( cst_pose.residue( pos ) );
			id::AtomID cst_atomID( cst_pose.residue_type(pos).atom_index("CA"), pos );
			id::StubID cst_fix_stub_ID( pose::named_stub_id_to_stub_id( id::NamedStubID( "N","CA","C", root_ ), cst_pose ) );
			Vector ai(
				numeric::random::rg().uniform()*perturb_,
				numeric::random::rg().uniform()*perturb_,
				numeric::random::rg().uniform()*perturb_
			);
			Vector xyz( rsd.xyz( cst_atomID.atomno() ) + ai );
			if ( bLocal_ ) {
				constraints_->add_constraint( ConstraintCOP( ConstraintOP( new scoring::constraints::LocalCoordinateConstraint(
					cst_atomID,
					cst_fix_stub_ID,
					xyz,
					cst_func_) ) ) );

			} else {
				constraints_->add_constraint( ConstraintCOP( ConstraintOP( new scoring::constraints::CoordinateConstraint(
					cst_atomID,
					id::AtomID( 1, root_ ) /*this is completely ignored! */,
					xyz,
					cst_func_) ) ) );
			}
		}

	}
	if ( tr.Debug.visible() ) {
		tr.Debug << "CoordConstraintClaimer generate constraints" << std::endl;
		constraints_->show_definition( std::cout, cst_pose );
	}
}

void CoordConstraintClaimer::read_constraints_from_file( pose::Pose const& cst_pose ) const {
	constraints_ = ConstraintIO::get_instance()->read_constraints(
		cst_filename_, ConstraintSetOP( new ConstraintSet ), cst_pose );

	if ( tr.Debug.visible() ) {
		tr.Debug << "CoordConstraintClaimer: have read constraints from file:" << std::endl;
		constraints_->show_definition( std::cout, cst_pose );
	}
}

void CoordConstraintClaimer::superimpose( pose::Pose const& pose ) const {
	ConstraintCOPs all_cst = constraints_->get_all_constraints();
	ConstraintSetOP new_set( new ConstraintSet );

	core::Size const natoms( all_cst.size() );
	ObjexxFCL::FArray1D_double weights( natoms, 1.0 );
	ObjexxFCL::FArray2D_double ref_coords( 3, natoms, 0.0 );
	ObjexxFCL::FArray2D_double coords( 3, natoms, 0.0 );

	if ( superimpose_regions_.size() ) {
		for ( Size i = 1; i<=natoms; ++i ) {
			if ( superimpose_regions_.is_loop_residue( i ) ) {
				weights( i ) = 1.0;
			} else {
				weights( i ) = 0.0;
			}
		}
	}

	Size n = 1;
	if ( tr.Trace.visible() ) {
		utility::io::ozstream dump_cst( label()+"_before_fit.cst" );
		constraints_->show_definition( dump_cst, pose );
	}

	for ( ConstraintCOPs::const_iterator it = all_cst.begin(), eit = all_cst.end(); it!=eit; ++it, ++n ) {
		//  ConstraintOP new_cst = (*it)->clone();
		LocalCoordinateConstraintCOP ll_cst = utility::pointer::dynamic_pointer_cast< LocalCoordinateConstraint const > ( *it );
		runtime_assert( ll_cst != nullptr ); //only these should be in the constraint set!
		Vector xyz_ref( pose.xyz( ll_cst->atom( 1 ) ) );
		//n = ll_cst->atom( 1 ).rsd(); //uncomment for debugging
		tr.Trace << "constraint for atom " << ll_cst->atom( 1 ) << std::endl;
		Vector xyz_cst( ll_cst->xyz_target( pose ) );
		for ( Size d=1; d<=3; ++d ) {
			ref_coords( d, n ) = xyz_ref( d );
			coords( d, n ) = xyz_cst( d );
		}
	}

	//fitting
	FArray1D_double transvec( 3 );
	FArray1D_double ref_transvec( 3 );
	//   if ( tr.Trace.visible() ) { //if uncommented use test_n instead of n for filling of coordinates -- only works if all residues have a xzy!
	//    toolbox::dump_as_pdb( "cst_xyz_before_fit.pdb", natoms, coords );
	//    toolbox::dump_as_pdb( "cst_ref_fit.pdb", natoms, ref_coords );
	//   };

	toolbox::reset_x( natoms, coords, weights, transvec );
	toolbox::reset_x( natoms, ref_coords, weights, ref_transvec );
	toolbox::Matrix R;
	toolbox::fit_centered_coords( natoms, weights, ref_coords, coords, R );
	//   if ( tr.Trace.visible() ) {
	//    toolbox::dump_as_pdb( "cst_xyz_after_fit.pdb", natoms, coords, -ref_transvec );
	//   };

	//filling xyz_constraints with new coords
	n = 1;
	for ( ConstraintCOPs::const_iterator it = all_cst.begin(), eit = all_cst.end(); it!=eit; ++it, ++n ) {
		ConstraintOP new_cst = (*it)->clone();
		LocalCoordinateConstraintOP ll_cst = utility::pointer::dynamic_pointer_cast< core::scoring::constraints::LocalCoordinateConstraint > ( new_cst );
		Vector xyz_cst;
		//  n = ll_cst->atom( 1 ).rsd();
		for ( Size d=1; d<=3; ++d ) {
			xyz_cst( d ) = coords( d, n ) - ref_transvec( d );
		}
		ll_cst->set_xyz_target( xyz_cst, pose );
		new_set->add_constraint( ll_cst );
	}
	constraints_ = new_set;
}

/// @detail read setup file
/*
PDB_FILE   cst_pose
CST_FILE   definition of CoordinateConstraints: which atoms, which potential
LABEL      some name --- not important
ROOT        where to put the reference atom
CST_FROM_INPUT_POSE take cst-xyz from threading model, or silent-in model
POTENTIAL   a cst-func, e.g., BOUNDED 0 0 4
*/
void CoordConstraintClaimer::set_defaults() {
	Parent::set_defaults();
	cst_filename_ = "NoFile";
	filename_ = "NoFile";
	root_ = 1;
	root_from_label_ = "ALL";
	perturb_ = 0;
	bUseXYZ_in_cstfile_ = false;
	bSuperimpose_ = false;
	superimpose_regions_.clear();
}

bool CoordConstraintClaimer::read_tag( std::string tag, std::istream& is ) {
	loops::PoseNumberedLoopFileReader loop_file_reader;
	loop_file_reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
	bool const strict = false; /*prohibit_single_residue_loops ? no*/

	if ( tag == "pdb_file" || tag == "PDB_FILE" ) {
		is >> filename_;
		read_cst_pose();
	} else if ( tag == "CST_FILE" ) {
		is >> cst_filename_;
	} else if ( tag == "REGION" ) {
		loops::SerializedLoopList loops = loop_file_reader.read_pose_numbered_loops_file( is, type(), strict );
		rigid_ = loops::Loops( loops );
	} else if ( tag == "region_file" ) {
		std::string file;
		is >> file;
		std::ifstream infile( file.c_str() );

		if ( !infile.good() ) {
			utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + file + "'" );
		}
		loops::SerializedLoopList loops = loop_file_reader.read_pose_numbered_loops_file( is, file, strict );
		rigid_ = loops::Loops( loops ); // <==
	} else if ( tag == "ROOT" ) {
		is >> root_;
	} else if ( tag == "ASK_FOR_ROOT" ) {
		bLocal_ = true;
		//std::string label;
		is >> root_from_label_;
	} else if ( tag == "CST_FROM_INPUT_POSE" ) {
		bRegenerateFromInputPose_ = true;
	} else if ( tag == "USE_XYZ_FROM_CSTFILE" ) { // I made this an explicit option, since it is probably often a bad choice!
		bUseXYZ_in_cstfile_ = true;
	} else if ( tag == "NO_CENTROID" ) {
		bCentroid_ = false;
	} else if ( tag == "FULLATOM" ) {
		bFullatom_ = true;
	} else if ( tag == "PERTURB" ) {
		is >> perturb_;
	} else if ( tag == "SUPERIMPOSE" ) {
		bSuperimpose_ = true;
	} else if ( tag == "SUPERIMPOSE_REGION" ) {
		loops::SerializedLoopList loops = loop_file_reader.read_pose_numbered_loops_file( is, type(), strict );
		superimpose_regions_ = loops::Loops( loops );
	} else if ( tag == "SUPERIMPOSE_REGION_FILE" ) {
		std::string file;
		is >> file;
		std::ifstream infile( file.c_str() );

		if ( !infile.good() ) {
			utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + file + "'" );
		}
		loops::SerializedLoopList loops = loop_file_reader.read_pose_numbered_loops_file( infile, file, strict );
		superimpose_regions_ = loops::Loops( loops ); // <==
	} else if ( tag == "POTENTIAL" ) {
		std::string func_type;
		is >> func_type;
		cst_func_ = ConstraintIO::get_func_factory().new_func( func_type );
		cst_func_->read_data( is );
	} else return Parent::read_tag( tag, is );
	return true;
}

void CoordConstraintClaimer::init_after_reading() {
	if ( cst_pose_ ) {
		if ( cst_filename_ != "NoFile" ) {
			read_constraints_from_file( *cst_pose_ );
		} else {
			if ( !cst_func_ ) {
				throw EXCN_Input( "POTENTIAL not specified for " + type() );
			}
			generate_constraints( *cst_pose_ );
		}
	}
	sequence_ = ""; //so we refresh the constraint_set at first call to add_constraints

}


} //topology_broker
} //protocols
