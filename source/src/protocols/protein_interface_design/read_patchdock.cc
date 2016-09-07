// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/protein_interface_design/design_utils.cc
/// @brief various utilities for interface design.
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Project Headers
#include <protocols/protein_interface_design/read_patchdock.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/izstream.hh>

// Unit Headers

// C++ headers
#include <map>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <numeric/xyz.functions.hh>


using namespace core;
using namespace core::scoring;
using basic::options::option;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.read_patchdock" );

struct Transformation
{
	typedef core::Real Real;

	Real alpha, beta, gamma; // Euler angles
	numeric::xyzVector< Real > translation; // translation
};

namespace protocols {
namespace protein_interface_design {

PatchdockReader::PatchdockReader(){
	clear_internals();
	random_entry( option[ parser::patchdock_random_entry ].user() );
	if ( random_entry() ) {
		utility::vector1< core::Size > entry_num_extrema = option[ parser::patchdock_random_entry ]();
		runtime_assert( entry_num_extrema.size() == 2 );
		from_entry( entry_num_extrema[ 1 ] );
		to_entry( entry_num_extrema[ 2 ] );
	}//fi patchdock_random_entry
	if ( option[ parser::patchdock ].user() ) {
		patchdock_fname_ = option[ parser::patchdock ]();
	}
}

PatchdockReader::~PatchdockReader() = default;

void
PatchdockReader::clear_internals()
{
	patchdock_fname_ = "";
	patchdock_entry_num_ = 0;
	saved_input_tag_ = "";
	saved_native_tag_ = "";
	saved_input_pose_ = saved_native_pose_ = nullptr;
	saved_transformations_.clear();
}

/// @details how many entries are there in the patchdock file?
core::Size
PatchdockReader::number_of_patchdock_entries()
{
	core::Size const saved_patchdock_entry_num( patchdock_entry_num_ );
	patchdock_entry_num_ = 1; // to ensure validity of patchdock_entry_num_
	read_patchdock_entry(); // to ensure that saved_transformations_ is set
	patchdock_entry_num_ = saved_patchdock_entry_num;
	return( saved_transformations_.size() );
}

/// @details read a rigid-body transformation from a patchdock file.
/// caches transformations in saved_transformations_ to avoid multiple disk access
Transformation
PatchdockReader::read_patchdock_entry()
{
	runtime_assert( patchdock_entry_num_ );
	if ( saved_transformations_.size() ) {
		return( saved_transformations_[ patchdock_entry_num_ ] );
	}

	utility::io::izstream data( patchdock_fname_ );
	if ( !data ) {
		utility_exit_with_message( "Cannot open patchdock file: " + patchdock_fname_ );
	}

	std::string line;
	bool entries_found( false );
	while ( getline( data, line ) ) {
		using namespace std;

		Transformation t;
		t.alpha = t.beta = t.gamma = 0;
		t.translation.zero() ;

		istringstream line_stream( line );
		string first_field;
		line_stream >> first_field;

		if ( first_field == "#" ) { entries_found = true; continue; }
		if ( !entries_found ) continue;
		core::Size const wheres_pipe( line.find_first_of( "|" ) );
		if ( wheres_pipe == string::npos ) break; // no longer reading entries
		core::Size const transformation_begin( line.find_last_of( "||" ) + 2 );
		std::istringstream transData( line.substr( transformation_begin, 10000) );
		core::Real x,y,z;
		transData >> t.alpha >> t.beta >> t.gamma >> x >> y >> z;
		if ( transData.fail() ) {
			TR<<"Error parsing transformation data in line\n"<<line<<std::endl;
			runtime_assert( !transData.fail() );
		}
		t.translation.assign( x, y, z );
		saved_transformations_.push_back( t );
	}
	return( saved_transformations_[ patchdock_entry_num_ ] );
}

//@details transform a chain within the pose according to t. The transformation computed here is
//based on patchdock's transOutput.pl and pdb_trans
void
PatchdockReader::transform_pose( core::pose::Pose & pose, core::Size const chain, Transformation const & t )
{
	core::Size const chain_begin( pose.conformation().chain_begin( chain ) );
	core::Size const chain_end( pose.conformation().chain_end( chain ) );

	numeric::xyzMatrix< core::Real > rotation;
	{ //compute rotation matrix (taken from RBSMover), but here expecting radian rather than degrees
		core::Real const sa ( std::sin( t.alpha ));
		core::Real const ca ( std::cos( t.alpha ));
		core::Real const sb ( std::sin( t.beta  ));
		core::Real const cb ( std::cos( t.beta  ));
		core::Real const sg ( std::sin( t.gamma ));
		core::Real const cg ( std::cos( t.gamma ));
		// Adapted from code sent by Dina Schneidman of the Wolfson lab (Tel-Aviv U)
		rotation.xx( cg * cb ); rotation.xy( -sb*sa*cg - sg * ca ); rotation.xz(  -sb*ca*cg + sg * sa );
		rotation.yx( sg * cb ); rotation.yy(  -sb*sa*sg + ca*cg ); rotation.yz( -sb*ca*sg - sa*cg );
		rotation.zx( sb );            rotation.zy( cb*sa );            rotation.zz(  cb*ca );
	}//compute rotation

	//rotate each atom around the geometric centre of the chain
	for ( core::Size residue=chain_begin; residue<=chain_end; ++residue ) {
		core::Size const atom_begin( 1 );
		core::Size const atom_end( pose.residue( residue ).natoms() );

		numeric::xyzVector< core::Real > localX, localRX;
		for ( core::Size atom=atom_begin; atom<=atom_end; ++atom ) {
			id::AtomID const id( atom, residue );

			localX = pose.xyz( id );
			localRX = rotation * localX;
			pose.set_xyz( id, localRX );
		}
	}

	//translate
	for ( core::Size residue=chain_begin; residue<=chain_end; ++residue ) {
		core::Size const atom_begin( 1 );
		core::Size const atom_end( pose.residue( residue ).natoms() );

		for ( core::Size atom=atom_begin; atom<=atom_end; ++atom ) {
			id::AtomID const id( atom, residue );

			numeric::xyzVector< core::Real > const new_pos( pose.xyz( id ) + t.translation );
			pose.set_xyz( id, new_pos );
		}
	}
	// detect disulfides
	if ( option[ in::detect_disulf ].user() ?
			option[ in::detect_disulf ]() : // detect_disulf true
			pose.is_fullatom() // detect_disulf default but fa pose
			) {
		pose.conformation().detect_disulfides();
	}
}

void
PatchdockReader::read_poses( core::pose::Pose & input_pose, std::string & input_tag )
{
	core::pose::Pose dummy_pose;
	std::string dummy_tag( input_tag );

	read_poses( input_pose, dummy_pose, input_tag, dummy_tag );
}

void
PatchdockReader::read_patchdock( std::string & input_tag, std::string & native_tag )
{
	patchdock_entry_num_ = 0;
	if ( patchdock_fname_ == "" ) { // use default patchdock fname ( 1jjj_2xxx.pdb.gz -> 1jjj_2xxx.patchdock )
		core::Size const filename_end( input_tag.find_first_of( "." ) );

		patchdock_fname_ = input_tag.substr( 0, filename_end ) + ".patchdock";
	}
	TR<<"Reading from patchdock file name: "<<patchdock_fname_<<std::endl;
	core::Size const number_of_entries = number_of_patchdock_entries();
	if ( number_of_entries == 0 ) {
		utility_exit_with_message_status( "No patchdock entries found. Aborting", 0 );
	}
	if ( random_entry() ) {
		runtime_assert( to_entry() >= from_entry() );
		runtime_assert( from_entry() > 0 );

		TR<<number_of_entries<<" entries in patchdock file "<<patchdock_fname_<<std::endl;
		core::Size const actual_last_entry( std::min( to_entry(), number_of_entries ) );
		TR<<"sampling a number between "<<from_entry()<<" and "<<actual_last_entry<<std::endl;

		patchdock_entry_num_ = ( core::Size ) floor( numeric::random::rg().uniform() * ( actual_last_entry - from_entry() + 1 ) ) + from_entry();
		runtime_assert( patchdock_entry_num_ <= actual_last_entry );
		runtime_assert( patchdock_entry_num_ >= from_entry() );
		std::stringstream ss;
		ss << "." << patchdock_entry_num_;
		option[ out::user_tag ].value( ss.str() ); // to set the output tag
	} else {
		core::Size const entrynum_begin( input_tag.find_first_of( "." ) );
		core::Size const entrynum_end( input_tag.find_last_of( "." ) );
		std::stringstream ss( input_tag.substr( entrynum_begin+1, entrynum_end - entrynum_begin ) );
		ss >> patchdock_entry_num_;
		if ( patchdock_entry_num_ > number_of_entries ) {
			TR<<"number of patchdock entries exceeded. You've asked for entry "<< patchdock_entry_num_<<" but only "<<number_of_entries<<" entries were found"<<std::endl;
			utility_exit_with_message_status("aborting.", 0 );
		}
		if ( input_tag == native_tag ) {
			input_tag.replace( entrynum_begin, entrynum_end - entrynum_begin + 1, "." );
			native_tag.replace( entrynum_begin, entrynum_end - entrynum_begin + 1, "." );
		}
	}
}

//@details transform an input pdb according to patchdock parameters
//if patchdock option is not set, return as is.
/// If the native and input tags match those saved in the object, then the
/// pose will not be read again from disk. Only the patchdock entry will be read.
void
PatchdockReader::read_poses( core::pose::Pose & input_pose, core::pose::Pose & native_pose, std::string & input_tag, std::string & native_tag )
{
	if ( saved_input_tag_ == input_tag && saved_native_tag_ == native_tag ) {
		//we've already read this pose, do not go to disk again
		input_pose = *saved_input_pose_;
		native_pose = *saved_native_pose_;
		TR<<"Skipped reading pose from disk"<<std::endl;
	} else {
		clear_internals();
		saved_input_tag_ = input_tag;
		saved_native_tag_ = native_tag;
		TR<<"Reading pose from disk"<<std::endl;
		if ( option[ in::file::centroid_input ].user() ) {
			core::import_pose::centroid_pose_from_pdb( input_pose, input_tag );
			core::import_pose::centroid_pose_from_pdb( native_pose,  native_tag);
		} else {
			core::chemical::ResidueTypeSetCOP rsd_set(
				core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" )
			);

			core::import_pose::pose_from_file( input_pose, *rsd_set, input_tag , core::import_pose::PDB_file);
			core::import_pose::pose_from_file( native_pose, *rsd_set, native_tag , core::import_pose::PDB_file);
		}//else
		if ( option[ in::file::fold_tree ].user() ) {
			std::string const fold_tree_fname( option[ in::file::fold_tree ]() );
			utility::io::izstream data( fold_tree_fname );
			if ( !data ) {
				TR << "Cannot open  file: " << fold_tree_fname << std::endl;
				runtime_assert( data );
			}
			std::string line;
			bool ft_found( false );
			while ( getline( data, line ) ) {
				if ( line.substr(0,10) == "FOLD_TREE " ) {
					std::istringstream line_stream( line );
					kinematics::FoldTree f;
					line_stream >> f;
					input_pose.fold_tree( f );
					native_pose.fold_tree( f );
					ft_found = true;
					TR<<"Using user-defined fold-tree:\n"<<input_pose.fold_tree()<<std::endl;
					break;
				}//IF FOLD_TREE
			}//getline
			runtime_assert( ft_found );
		}//option foldtree

		saved_input_pose_ = core::pose::PoseOP( new core::pose::Pose( input_pose ) );
		saved_native_pose_ = core::pose::PoseOP( new core::pose::Pose( native_pose ) );
	}//else

	read_patchdock( input_tag, native_tag );

	if ( !patchdock_entry_num_ ) return; // no need for transformations
	TR<<"Reading patchdock entry "<<patchdock_entry_num_<<" from file: "<<patchdock_fname_<<std::endl;
	Transformation t( read_patchdock_entry() );

	transform_pose( input_pose, 2/*chain*/, t );
	transform_pose( native_pose, 2, t );
	TR.flush();
}

}//protein_interface_design
}//protocols

