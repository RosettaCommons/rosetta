// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/BackboneDB.cc
/// @brief
/// @author Mike Tyka
/// @author Ken Jung
#include <protocols/loophash/BackboneDB.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
#include <protocols/loophash/Exceptions.hh>
#include <boost/lexical_cast.hpp>
#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <iostream>
#include <fstream>
#include <sstream>

#include <protocols/frag_picker/VallChunk.hh>
#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace kinematics;
using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace loophash {

static THREAD_LOCAL basic::Tracer TR( "BackboneDB" );

short RealAngleToShort( core::Real angle ){
	while ( angle > 180.0 ) angle -= 360.0;
	while ( angle <-180.0 ) angle += 360.0;
	// range for short: -32768 to 32767
	short result = short( angle * 182.0 );
	return result;
}

core::Real ShortToRealAngle( short angle ){
	core::Real result = core::Real( angle ) / 182.0;
	return result;
}

void
BackboneSegment::apply_to_pose( core::pose::Pose &pose, core::Size ir, bool cut ) const
{
	core::Size length = phi_.size();
	core::Size jr = ir + length;
	if ( cut ) {
		//fpd vrt/ligand trim
		core::Size newroot=0;
		if ( pose.residue_type( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) newroot = pose.fold_tree().root();

		core::Size nres = pose.size();
		while ( !pose.residue_type(nres).is_polymer() ) nres--;

		// get current cutpoints; don't try to connect these
		utility::vector1< Size > cuts_in = pose.fold_tree().cutpoints();
		std::sort( cuts_in.begin(), cuts_in.end() );

		// bail if (ir,jr) crosses a cut
		for ( Size i=1; i<=cuts_in.size(); ++i ) {
			if ( cuts_in[i]<=jr && cuts_in[i]>=ir ) {
				TR.Error << "ERROR -- residue range crosses cut    IR: " << ir << "  JR: " << jr << "  CUT: " << cuts_in[i] << std::endl;
				return;
			}
			//fpd insertions one position after the cut seem not to work ...
			//fpd perhaps if the foldtree for the local segment were reversed this might be ok
			if ( cuts_in[i]==ir-1 ) {
				TR.Error << "ERROR -- startres immediately follows cut    IR: " << ir << "  CUT: " << cuts_in[i] << std::endl;
				return;
			}
		}

		//fpd handle multiple chains/chainbreaks
		FoldTree f;
		core::Size last_cut=0, jump_num=2;
		Size cutpoint= jr-1;
		for ( Size i=1; i<=cuts_in.size(); ++i ) {
			if ( cuts_in[i] >= nres ) break;
			if ( cutpoint > last_cut && cutpoint < cuts_in[i] ) {
				f.add_edge( last_cut+1, ir, Edge::PEPTIDE );
				f.add_edge( ir, cutpoint, Edge::PEPTIDE );
				f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
				f.add_edge( jr, cuts_in[i] , Edge::PEPTIDE );
				f.add_edge( ir, jr, 1 );  // this is the jump !!
				if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
			} else {
				f.add_edge( last_cut+1, cuts_in[i], Edge::PEPTIDE );
				if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
			}
			last_cut = cuts_in[i];
		}
		if ( last_cut+1 <= nres ) {
			if ( cutpoint > last_cut && cutpoint < nres ) {
				f.add_edge( last_cut+1, ir, Edge::PEPTIDE );
				f.add_edge( ir, cutpoint, Edge::PEPTIDE );
				f.add_edge( cutpoint + 1, jr, Edge::PEPTIDE );
				f.add_edge( jr, nres , Edge::PEPTIDE );
				f.add_edge( ir, jr, 1 );  // this is the jump !!
				if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
			} else {
				f.add_edge( last_cut+1, nres, Edge::PEPTIDE );
				if ( last_cut!=0 ) f.add_edge( 1, last_cut+1, jump_num++);
			}
		}
		for ( core::Size i=nres+1; i<=pose.size(); ++i ) {
			f.add_edge( 1, i, jump_num++ );  // additional jumps
		}

		core::Size theroot = 1;
		if ( ir == 1 ) theroot = pose.size();
		if ( newroot>0 ) theroot = newroot;  //fpd
		if ( f.reorder(theroot) == false ) {
			TR.Error << "ERROR During reordering of fold tree - am ignoring this LOOP ! bailing: The root: " << theroot << " NRES " << pose.size() << "   IR: " << ir << "  JR: " << jr << std::endl;
			return; // continuing leads to a segfault - instead ignore this loop !
		}
		pose.fold_tree(f);
	}

	for ( core::Size i = 0; i < length; i++ ) {
		core::Size ires = ir + i;
		if ( ires > pose.size() ) return;
		pose.set_phi( ires, phi_[i] );
		pose.set_psi( ires, psi_[i] );
		pose.set_omega( ires, omega_[i] );
	}
}

void
BackboneSegment::read_from_pose( core::pose::Pose const &pose, core::Size ir, core::Size length )
{
	phi_.clear();
	psi_.clear();
	omega_.clear();

	for ( core::Size i = 0; i < length; i++ ) {
		core::Size ires = ir + i;
		if ( ires > pose.size() ) return;

		phi_.  push_back( pose.phi(   ires ));
		psi_.  push_back( pose.psi(   ires ));
		omega_.push_back( pose.omega( ires ));
	}
}

void BackboneSegment::print() const {
	for ( double it : phi_ ) TR << it << "  " ;
	for ( double it : psi_ ) TR << it << "  " ;
	for ( double it : omega_ ) TR << it << "  " ;
	TR << std::endl;
}

bool BackboneSegment::compare(const BackboneSegment &bs1, core::Real tolerance) const {
	const BackboneSegment &bs2 = (*this);

	if ( bs1.phi().size() != bs2.phi().size() ) return false;
	if ( bs1.psi().size() != bs2.psi().size() ) return false;
	if ( bs1.omega().size() != bs2.omega().size() ) return false;

	for ( core::Size i = 0; i < bs1.phi().size(); i++ )    if ( (bs1.phi()[i]    - bs2.phi()[i]) > tolerance ) return false;
	for ( core::Size i = 0; i < bs1.psi().size(); i++ )    if ( (bs1.psi()[i]    - bs2.psi()[i])  > tolerance ) return false;
	for ( core::Size i = 0; i < bs1.omega().size(); i++ )  if ( (bs1.omega()[i]  - bs2.omega()[i])  > tolerance ) return false;

	return true;
}
bool BackboneSegment::operator==( const BackboneSegment &bs1 ) const {
	const BackboneSegment &bs2 = (*this);

	if ( bs1.phi().size() != bs2.phi().size() ) return false;
	if ( bs1.psi().size() != bs2.psi().size() ) return false;
	if ( bs1.omega().size() != bs2.omega().size() ) return false;

	for ( core::Size i = 0; i < bs1.phi().size(); i++ ) if ( bs1.phi()[i] == bs2.phi()[i] ) return false;
	for ( core::Size i = 0; i < bs1.psi().size(); i++ ) if ( bs1.psi()[i] == bs2.psi()[i] ) return false;
	for ( core::Size i = 0; i < bs1.omega().size(); i++ ) if ( bs1.omega()[i] == bs2.omega()[i] ) return false;

	return true;
}

core::Real get_rmsd( const BackboneSegment &bs1, const BackboneSegment &bs2 ){
	core::Real sumsqr = 0;
	core::Size count = 0;
	if ( bs1.phi().size() != bs2.phi().size() ) return -1;
	for ( core::Size i = 0; i < bs1.phi().size(); i++ ) {
		if ( bs1.phi()[i] == 0 || bs2.phi()[i] == 0 ) continue;
		core::Real diff = bs1.phi()[i] - bs2.phi()[i];
		while ( diff > 180  ) diff -= 360;
		while ( diff < -180  ) diff += 360;
		sumsqr += diff*diff;
		count++;
	}
	if ( bs1.psi().size() != bs2.psi().size() ) return -1;
	for ( core::Size i = 0; i < bs1.psi().size(); i++ ) {
		if ( bs1.psi()[i] == 0 || bs2.psi()[i] == 0 ) continue;
		core::Real diff = bs1.psi()[i] - bs2.psi()[i];
		while ( diff > 180  ) diff -= 360;
		while ( diff < -180  ) diff += 360;
		sumsqr += diff*diff;
		count++;
	}
	if ( bs1.omega().size() != bs2.omega().size() ) return -1;
	for ( core::Size i = 0; i < bs1.omega().size(); i++ ) {
		if ( bs1.omega()[i] == 0 || bs2.omega()[i] == 0 ) continue;
		core::Real diff = bs1.omega()[i] - bs2.omega()[i];
		while ( diff > 180  ) diff -= 360;
		while ( diff < -180  ) diff += 360;
		sumsqr += diff*diff;
		count++;
	}

	return sqrt( sumsqr/core::Real(count) );
}

core::Real
BackboneDB::angle( core::Size index, core::Size offset )
{
	if ( index >= data_.size() ) utility_exit_with_message( "Out of bounds error" );
	if ( offset >= data_[index].angles.size() ) utility_exit_with_message( "Out of bounds error" );

	return ShortToRealAngle( data_[ index ].angles[ offset ] );
}

void
BackboneDB::add_pose( const core::pose::Pose &pose, core::Size nres, core::Size &index, protocols::frag_picker::VallChunkOP chunk )
{
	if ( ! extra_ ) extra_ = true;
	index = data_.size(); // Index of protein
	BBData new_protein;
	for ( core::Size i = 0; i < nres; i++ ) {
		new_protein.angles.push_back( RealAngleToShort(pose.phi( 1 + i )));
		new_protein.angles.push_back( RealAngleToShort(pose.psi( 1 + i )));
		new_protein.angles.push_back( RealAngleToShort(pose.omega( 1 + i )));
	}
	BBExtraData extra_data;
	if ( chunk ) {
		//could modify this to move the chunk processing to LoopHashLibrary
		extra_data.sequence = chunk->get_sequence();
		extra_data.pdb_id = chunk->get_pdb_id() + chunk->get_chain_id();
	} else {
		extra_data.sequence = pose.sequence();
		std::string pose_id="";
		get_score_line_string( pose, "usid", pose_id );
		extra_data.pdb_id = pose_id;
	}
	new_protein.extra_key = extra_data_.size();
	extra_data_.push_back( extra_data );
	data_.push_back( new_protein );
}

// Maybe I should just overload the copy operator in the struct..
void BackboneDB::get_protein( core::Size index, BBData & protein ) const {
	protein.extra_key = data_[index].extra_key;
	protein.angles = data_[index].angles;
}

void BackboneDB::get_extra_data( core::Size index, BBExtraData & extra ) const {
	extra = extra_data_[index];
}

void BackboneDB::add_protein( BBData new_protein ) {
	data_.push_back( new_protein );
}

void BackboneDB::add_extra_data( BBExtraData extra ) {
	if ( ! extra_ ) extra_ = true;
	extra_data_.push_back( extra );
}

void
BackboneDB::get_backbone_segment(
	core::Size index,
	core::Size offset,
	core::Size len,
	BackboneSegment &bs
) const
{
	std::vector<core::Real> phi;
	std::vector<core::Real> psi;
	std::vector<core::Real> omega;
	core::Size pos = offset;
	for ( core::Size i = 0; i < len; i++ ) {
		phi.push_back( ShortToRealAngle(data_[index].angles[pos]) ); pos ++ ;
		psi.push_back( ShortToRealAngle(data_[index].angles[pos]) ); pos ++ ;
		omega.push_back( ShortToRealAngle(data_[index].angles[pos]) ); pos ++ ;
	}
	bs = BackboneSegment( phi, psi, omega );
}

void BackboneDB::write_db( std::string filename )
{
	std::ofstream file( filename.c_str() );
	if ( !file ) throw EXCN_DB_IO_Failed( filename, "write" );
	if ( data_.size() == 0 ) {
		file.close();
		return;
	}
	if ( ! extra_ ) throw EXCN_No_Extra_Data_To_Write();
	for ( auto & i : data_ ) {
		file << "pdb " << extra_data_[ i.extra_key ].pdb_id << std::endl;
		file << "seq " << extra_data_[ i.extra_key ].sequence << std::endl;
		file << "rot ";
		for (int j : extra_data_[ i.extra_key ].rotamer_id) {
			file << j << " ";
		}
		file << std::endl;
		file << "ang ";
		for (short angle : i.angles) {
			file << angle << " ";
		}
		file << std::endl;
	}
	file.close();
}

void
BackboneDB::read_legacydb( std::string filename )
{
	// use basic C input - C++ are too memory hungry to deal with these potentially v large files
	FILE *file = fopen( filename.c_str(), "r" );
	if ( file == nullptr ) throw EXCN_DB_IO_Failed( filename, "read" );

	data_.clear();
	BBData new_protein;
	data_.push_back( new_protein );
	unsigned count = 0;
	while ( !feof( file ) ) {
		count++;
		TR.Debug << "C: " << count << std::endl;
		const unsigned int bufsize = 16384;
		short bufferdata[16384];
		size_t readshorts = fread(&bufferdata[0],sizeof(short),bufsize,file);
		for ( unsigned i = 0; i< readshorts; i ++ ) {
			data_[0].angles.push_back( bufferdata[i] );
		}
	}
	fclose( file );
	TR.Debug << "End of read_db_from_binary" << std::endl;
}

void
BackboneDB::read_db( std::string filename, bool load_extra,
	core::Size num_partitions, core::Size assigned_num,
	std::pair< core::Size, core::Size > & loopdb_range,
	std::map< core::Size, bool > & homolog_index )
{
	std::ifstream file( filename.c_str() );
	if ( !file ) throw EXCN_DB_IO_Failed( filename, "read" );

	if ( option[ lh::exclude_homo ]() ) {
		TR << "Reading in homolog file" << std::endl;
		read_homologs();
	}

	extra_ = true;
	std::string line;

	core::Size num_lines = 0;
	// get number of lines in db
	while ( getline(file, line) ) {
		num_lines++;
	}

	// truncating to integer is good
	core::Size begin = assigned_num * num_lines / 4 / num_partitions;
	core::Size end = ( assigned_num + 1 ) * num_lines / 4 / num_partitions;
	//if( assigned_num == num_partitions - 1 ) end = 0;
	loopdb_range.first = begin;
	loopdb_range.second = end;

	TR.Info << "Reading in proteins " << begin << " to " << end << " out of " << num_lines / 4 << " , partition: "  << assigned_num+1 << "/"<< num_partitions << std::endl;
	// clear eof bit
	file.clear();
	file.seekg( 0, std::ios_base::beg );

	BBData new_protein;
	BBExtraData extra_data;
	std::string command;
	int line_counter = -1;
	bool is_homolog = false;

	unsigned int stat_count_protein = 0;
	while ( getline( file, line ) ) {
		line_counter++;
		if ( line_counter / 4 < int(begin) ) continue;
		if ( line_counter / 4 >= int(end) && int(end) != 0 ) continue;

		command = line.substr(0,3);
		if ( command == "" ) throw EXCN_Wrong_DB_Format( filename );
		// Even if we're not loading extra, still process pdb line
		// So we can use the pdb to filter homologs
		if ( command == "pdb" ) {
			new_protein.extra_key = extra_data_.size();
			extra_data.pdb_id = line.substr( 4 );
			stat_count_protein++; // count the number of proteins read in
			if ( homologs_.find( extra_data.pdb_id ) != homologs_.end() ) is_homolog = true;
		}
		if ( load_extra ) {
			if ( command == "seq" ) {
				extra_data.sequence = line.substr( 4 );
			}
			if ( command == "rot" ) {
				std::string buf;
				std::stringstream ss( line.substr( 4 ) );
				while ( ss >> buf )
						extra_data.rotamer_id.push_back( boost::lexical_cast< int >( buf ) );
			}
		}
		if ( command == "ang" ) {
			std::string buf;
			std::stringstream ss( line.substr( 4 ) );
			while ( ss >> buf ) new_protein.angles.push_back( boost::lexical_cast< short >( buf ) );
			if ( is_homolog ) {
				// Still leave a holder protein in data_ so indices in leapindex aren't messed up
				// but with no data so it doesn't take up space
				new_protein.angles.clear();
				// then add the index to the homolog map
				homolog_index[ data_.size() ] = true;
				TR << "Homolog " << extra_data.pdb_id << " rejected." << std::endl;
			}
			data_.push_back(new_protein);
			is_homolog = false;
			new_protein.angles.clear();
			if ( load_extra ) {
				// add extra data for holder proteins, since its not that much data
				// if this becomes too large, we can change later
				extra_data_.push_back( extra_data );
				extra_data.rotamer_id.clear();
			}
		}
	}
	TR.Info << "Read in " << stat_count_protein << " proteins" << std::endl;
	TR.Info << "Data_ size " << data_.size() << std::endl;
	file.close();
}

void BackboneDB::read_homologs()
{
	std::ifstream file( option[ lh::homo_file ]().c_str() );
	if ( !file ) throw EXCN_DB_IO_Failed( option[ lh::homo_file ](), "read" );
	std::string line;
	while ( getline( file, line) ) {
		utility::vector1< std::string > tokens ( utility::split( line ) );
		for ( utility::vector1< std::string >::const_iterator token = tokens.begin(); token != tokens.end(); ++token ) {
			std::string homolog_pdb_code_and_chain = (*token);
			TR << "Adding homolog: " << homolog_pdb_code_and_chain << std::endl;
			homologs_[homolog_pdb_code_and_chain] = true;
			if ( homolog_pdb_code_and_chain.size() == 5 ) {
				if ( homolog_pdb_code_and_chain[4] == 'A' || homolog_pdb_code_and_chain[4] == 'a' ) homologs_[ homolog_pdb_code_and_chain.replace( 4, 1, 1, '_' ) ] = true;
				if ( homolog_pdb_code_and_chain[4] == '_' ) homologs_[ homolog_pdb_code_and_chain.replace( 4, 1, 1, 'A' ) ] = true;
			}
		}
	}

	TR << "Homolog exclusion: ";
	for ( std::map< std::string, bool >::const_iterator hom = homologs_.begin(); hom != homologs_.end(); ++hom ) {
		TR << hom->first << " ";
	}
	TR << std::endl;

}

} // namespace loops
} // namespace protocols


