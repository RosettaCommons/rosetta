/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/raw_data/DisulfideFile.hh
/// @brief  A simple wrapper for a Disulfide File suitable for the -fix_disulf option.
/// @author Spencer Bliven <blivens@u.washington.edu>

#include <core/io/raw_data/DisulfideFile.hh>

#include <core/types.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/conformation/Conformation.hh>

// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <utility/vector1.hh>



namespace core {
namespace io {
namespace raw_data {

using namespace std;
using core::Size;
using utility::vector1;
using core::pose::Pose;

static basic::Tracer TR("core.io.raw_data.DisulfideFile");

/// @brief Get a list of disulfide bonds declared in the file
/// @details
/// The first time this is called it reads
/// the file and caches the resulting pairs of residues.
/// Subsequent calls to disulfides() are fast since they don't reparse the file
/// but merely reinterpret the results in terms of the specified.
///
/// This version ignores the possibility of PDB numbering and just reports on the file itself
///
/// @param[out] disulfides Appends pairs of residues specified by the file to this list
/// @param[in] none
///
/// @postcondition For each pair of residues (l,u) in this list,  l<u
void DisulfideFile::disulfides(vector1< pair<Size,Size> > & disulfides ) const {
	disulfides.resize(disulfides.size()+disulfides_.size());

	if(! up_to_date_ ) {
		parse_disulf_file();
		//mark cache as up to date
		up_to_date_ = true;
	}
	for(vector1< pair<ResNum,ResNum> >::const_iterator disulf = disulfides_.begin(),
			end_disulf = disulfides_.end();
			disulf != end_disulf; ++disulf)
	{
		Size l = disulf->first.n;
		Size u = disulf->second.n;
		disulfides.push_back(std::make_pair(l,u) );
	}
}


/// @brief Get a list of disulfide bonds declared in the file
/// @details
/// The first time this is called it reads
/// the file and caches the resulting pairs of residues.
/// Subsequent calls to disulfides() are fast since they don't reparse the file
/// but merely reinterpret the results in terms of the specified.
///
/// @param[out] disulfides Appends pairs of residues specified by the file to this list
/// @param[in] pose Only used if PDB numbering is used in the disulfide file
///
/// @postcondition For each pair of residues (l,u) in this list,  l<u
void DisulfideFile::disulfides(vector1< pair<Size,Size> > & disulfides, Pose const& pose) const {
	disulfides.resize(disulfides.size()+disulfides_.size());

	if(! up_to_date_ ) {
		parse_disulf_file();
		//mark cache as up to date
		up_to_date_ = true;
	}
	for(vector1< pair<ResNum,ResNum> >::const_iterator disulf = disulfides_.begin(),
			end_disulf = disulfides_.end();
			disulf != end_disulf; ++disulf)
	{
		Size l = resnum_to_rosetta_num(pose, disulf->first);
		Size u = resnum_to_rosetta_num(pose, disulf->second);
		disulfides.push_back(std::make_pair(l,u) );
	}
}

/// @brief Get a list of disulfide bonds declared in the file and manually set them in
/// the conformation
/// @details
/// The first time this is called it reads
/// the file and caches the resulting pairs of residues.
/// Subsequent calls to disulfides() are fast since they don't reparse the file
/// but merely reinterpret the results in terms of the specified.
/// (this is a necessary workaround for dealing with multiple disulfide specification files
/// in PyRosetta)  Because there is no equivalent of std::pair in python this method
/// avoids needing to provide a vector1 of std::pairs as arguments
///
/// @param[out] disulfides Appends pairs of residues specified by the file to this list
/// @param[in] pose is used if PDB numbering is used in the disulfide file and is also
/// needed to manually set the disulfides in the conformation
///
/// @postcondition For each pair of residues (l,u) in this list,  l<u
void DisulfideFile::read_in_and_set_disulfides( Pose &pose) {
    utility::vector1< std::pair<Size, Size> > disulfides;
        
        
    if(! up_to_date_ ) {
        parse_disulf_file();
        //mark cache as up to date
        up_to_date_ = true;
    }
    for(vector1< pair<ResNum,ResNum> >::const_iterator disulf = disulfides_.begin(),
        end_disulf = disulfides_.end();
        disulf != end_disulf; ++disulf)
    {
        Size l = resnum_to_rosetta_num(pose, disulf->first);
        Size u = resnum_to_rosetta_num(pose, disulf->second);
        disulfides.push_back(std::make_pair(l,u) );
    }
    pose.conformation().fix_disulfides( disulfides );
}
 
    
/// @brief Convert a ResNum object into the rosetta residue index
/// @details For ResNums with the rosetta_num type, this just return the n field.
///  For pdb_num ResNums the pose's PDB_Info is used to translate to rosetta numbering.
///
///  This function exits with an error message if it is unable to do the conversion.
Size DisulfideFile::resnum_to_rosetta_num(Pose const& pose, ResNum const& resnum) const
{
	using namespace core::pose;

	// Rosetta number
	if( resnum.type == rosetta_num ||
		(resnum.type == unknown_num && resnum.chain == 0) )
	{
		return resnum.n;
	}

	//PDB number
	PDB_InfoCOP info( pose.pdb_info() );
	if( info == 0 ) {
		TR.Error << "[ERROR] PDB Number expected from format, but no PDB Info present."
			<< std::endl;
		utility_exit();
	}
	Size n( info->pdb2pose(resnum.chain, resnum.n) );
	if( n == 0 ) {
		TR.Error << "[ERROR] PDB Number " << resnum.n << resnum.chain
			<< " does not correspond to a valid residue." << std::endl;
		utility_exit();
	}
	return n;

}

/// @brief Parses residue numbers out of a disulfide file
/// See \link core::io::raw_data::DisulfideFile DisulfideFile \endlink for more
/// details on the file format.
void DisulfideFile::parse_disulf_file() const {
	//Open disulfide file
	utility::io::izstream disulf_stm;
	disulf_stm.open( filename_ );
	if( disulf_stm.fail() ) {
		TR.Error << "[ERROR] Unable to open disulfide file " << filename_ << "."<< std::endl;
		utility_exit();
	}

	disulfides_.clear();

	//skip whitespace
	disulf_stm >> skipws;
	while(disulf_stm.good() ) {
		//Decide what to do based on the next character
		int next_char = disulf_stm.peek();
		if(disulf_stm.eof())
			break;
		if(disulf_stm.fail() || disulf_stm.bad() ) {
			TR.Error << "[ERROR] Error reading disulfide file " << filename_ << "." << std::endl;
			utility_exit();
		}
		switch( next_char ) {
			case '#': {
				//Ignore comments
				std::string line;
				disulf_stm.getline(line);
				continue;
			}
			default: {
				//Read pair of residue indices.
				Size l,u;
				disulf_stm >> l >> u >> std::ws;
				if(disulf_stm.fail() || disulf_stm.bad() ) {
					TR.Error << "[ERROR] Error reading disulfide file " << filename_ << "." << std::endl;
					utility_exit();
				}

				//Guarantee that l < u
				if( u < l) {
					Size tmp = u;
					u = l;
					l = tmp;
				}
				TR.Info << "Fixing a disulfide between "
					<< l << " and " << u << std::endl;
				ResNum lres = { l,0,rosetta_num };
				ResNum ures = { u,0,rosetta_num };
				disulfides_.push_back(make_pair(lres,ures));

				continue;
			}
		}
	}
	//clean up
	disulf_stm.close();

	up_to_date_ = true;
}

} //raw_data
} //io
} //core
