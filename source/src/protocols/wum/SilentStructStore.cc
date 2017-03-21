// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/wum/SilentStructStore.cc
/// @brief
/// @author Mike Tyka

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/wum/SilentStructStore.hh>
//#include <protocols/frag_picker/VallChunk.hh>
//#include <protocols/frag_picker/VallProvider.hh>

#include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <basic/options/option.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/io/izstream.hh>


// C++ headers
//#include <cstdlib>

#include <iostream>
#include <fstream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

//Auto Headers
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <boost/unordered/unordered_map.hpp>


using namespace core;
using namespace kinematics;
//using namespace protocols::frag_picker;
using namespace core::io::silent;
using namespace core::pose;
using namespace core::scoring;
using namespace conformation;
using namespace protocols::moves;


namespace protocols {
namespace wum {

/// @details Auto-generated virtual destructor
SilentStructStore::~SilentStructStore() = default;


class sort_SilentStructOPs
{
public:
	sort_SilentStructOPs(std::string field = "score" ): field_(std::move(field)) {}

	bool operator () (const SilentStructOP& left, const SilentStructOP& right)
	{
		runtime_assert( left != nullptr );
		runtime_assert( right != nullptr );
		return left->get_energy( field_ ) < right->get_energy( field_ );
	}
private:
	std::string field_;
};


bool find_SilentStructOPs::operator () (const core::io::silent::SilentStructOP& check)
{
	if ( check->get_energy(field_) == value_ ) return true;
	return false;
}


static THREAD_LOCAL basic::Tracer TR( "SilentStructStore" );


void
SilentStructStore::clear()
{
	store_.clear();
}

void
SilentStructStore::add( const core::pose::Pose &pose ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentStructOP ss = option[ lh::bss]() ?
		core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary", opts) :
		core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts );
	ss->fill_struct( pose );
	add( ss );
}


void
SilentStructStore::add( SilentStructOP new_struct ){
	store_.push_back( new_struct );
}

void
SilentStructStore::add( const SilentStruct &new_struct ){
	SilentStructOP pss = new_struct.clone();
	store_.push_back( pss );
}

void
SilentStructStore::add( core::io::silent::SilentFileData const& sfd ) {
	using namespace core::io::silent;
	using namespace core::chemical;
	for ( SilentFileData::const_iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
		add( *it );
	}
}

void
SilentStructStore::add( SilentStructStore &mergestore ) {
	for ( std::vector < SilentStructOP >::const_iterator it = mergestore.store_.begin();
			it != mergestore.store_.end();
			++it ) {
		runtime_assert( *it != nullptr );
		store_.push_back( *it );
	}
}


// @brief This uses the pose stream to read in everything from -l, -s and -in:file:silent into this store.
void
SilentStructStore::read_from_cmd_line( ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::chemical::ResidueTypeSetCOP rsd_set;
	if ( option[ in::file::fullatom ]() ) {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	} else {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}

	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
	core::Size count = 0;
	while ( input.has_another_pose() && (count < 400 ) ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );
		add( pose );
		count ++;
	}
	TR.Info << "Read " << count << " structures from command line" << std::endl;
}


// @brief read from silent file
void
SilentStructStore::read_from_string( const std::string & input )  {
	std::istringstream iss(input);
	read_from_stream( iss );
}

// @brief read from silent file
void
SilentStructStore::read_from_stream( std::istream & input )  {
	SilentFileOptions opts;
	SilentFileData sfd( opts );
	utility::vector1< std::string > tags_wanted; // empty vector signals "read all" according to author of silent io
	sfd.read_stream( input, tags_wanted, false  );
	utility::vector1< std::string> comments = sfd.comment_lines();
	// utility::vector1< std::string >::iterator citer = comments.begin(); // Unused variable was causing a warning
	// Now loop over each structure in that silent file
	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
		add( *iter);
	}
}

void
SilentStructStore::read_from_file( const std::string &filename ){
	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message(
			"ERROR: Unable to open silent strcuture store file: '" + filename + "'"
		);
	}
	read_from_stream( data );
	data.close();
}


void
SilentStructStore::get_pose( core::Size index,  core::pose::Pose &pose ) const {
	runtime_assert( index < store_.size() );
	SilentStructCOP temp_struct = get_struct( index );
	temp_struct->fill_pose( pose );
}

// @brief GEt a random structure
SilentStructCOP SilentStructStore::get_struct_random() const{
	runtime_assert( store_.size() > 0 );
	core::Size choice=core::Size( numeric::random::rg().random_range(0,(store_.size()-1)));
	runtime_assert( choice < store_.size() );
	return store_[ choice ];
}

void SilentStructStore::serialize( std::ostream & out ) const {
	if ( store_.size() == 0 ) {
		TR.Warning << "Empty silent struct store serialized." << std::endl;
	} else {
		(*store_.begin())->print_header( out );
		SilentFileOptions opts;
		SilentFileData sfd( opts );
		for ( auto const & it : store_ ) {
			runtime_assert( it != nullptr );
			sfd.write_silent_struct( (*it), out );
		}
	}
}

void SilentStructStore::serialize( std::string & out ) const {
	std::ostringstream ss;
	serialize( ss );
	out = ss.str();
}


void SilentStructStore::serialize_to_file( const std::string &filename ) const  {
	std::ofstream out( filename.c_str() );
	if ( !out.good() ) {
		utility_exit_with_message( "ERROR: Unable to open output file : '" + filename + "'" );
	}
	serialize( out );
	out.close();
}


void SilentStructStore::print( std::ostream & out ) const {
	SilentFileOptions opts;
	SilentFileData sfd( opts );
	core::Size count=0;
	out << "----------------------------------------------" << std::endl;
	for ( auto const & it : store_ ) {
		out << count << " ";
		runtime_assert( it != nullptr );
		it->print_scores( out );
	}
	out << "----------------------------------------------" << std::endl;
}


core::Size
SilentStructStore::mem_footprint() const {
	core::Size total = 0;
	for ( auto const & it : store_ ) {
		total += it->mem_footprint();
	}
	return total;
}


void
SilentStructStore::sort_by( std::string field )
{
	if ( store_.size() == 0 ) return;
	sort_SilentStructOPs sort_by_field = field;
	std::sort( store_.begin(), store_.end(), sort_by_field );
}


void
SilentStructStore::all_add_energy( std::string scorename, core::Real value, core::Real weight )
{
	for ( auto & it : store_ ) {
		it->add_energy( scorename, value, weight );
	}
}

void
SilentStructStore::all_sort_silent_scores()
{
	for ( auto & it : store_ ) {
		it->sort_silent_scores( );
	}
}


// TOOLS


std::string
encode_alphanum(unsigned long number, int pad_width=0, char pad_char = '0')
{
	std::string code = "";
	const static char *codes = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	const static int codes_len = 52;

	int count_len = 0;
	while ( number > 0 ) {
		unsigned int digit = number % codes_len;
		code = codes[digit] + code;
		count_len++;
		number /= codes_len;
	}

	while ( count_len < pad_width ) {
		code = pad_char + code;
		count_len++;
	}

	return code;
}


// without locking this is not really threadsafe. But how to otherwise proved a *globally*
// unique number ? A singleton wont help either, just code overhead.

std::string
generate_unique_structure_id(){
	static long unique_count=0;
	int mpi_rank = 0;
	int mpi_npes = 0;
#ifdef USEMPI
		MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &mpi_rank ) );
		MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &mpi_npes ) );
#endif
	int width = int(floor(log(float(mpi_npes))/log(62.0)) + 1.0);  // 62 is the base of the coded number below.
	unique_count++;
	return encode_alphanum( unique_count ) + encode_alphanum( mpi_rank, width, '0' );
}


}
}

