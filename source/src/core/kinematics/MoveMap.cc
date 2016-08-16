// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/MoveMap.cc
/// @brief  Method definitions for the MoveMap
/// @author Phil Bradley
/// @author Christopher Miles (cmiles@uw.edu)
/// @author Roland A. Pache

// Unit header
#include <core/kinematics/MoveMap.hh>

// Project headers
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>

// Utility headers
#include <utility/py/PyAssert.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// External headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <vector>


using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


namespace core {
namespace kinematics {

// Auto-generated virtual destructor
MoveMap::~MoveMap() {}

void MoveMap::set_ranges_unmodifiable(const std::vector<std::pair<Size, Size> >& ranges) {
	using std::pair;
	using std::vector;

	set_bb(true);
	vector<pair<Size, Size> >::const_iterator i;
	for ( i = ranges.begin(); i != ranges.end(); ++i ) {
		Size begin = i->first;
		Size end = i->second;

		debug_assert(begin > 0);
		debug_assert(end > 0);
		debug_assert(begin <= end);

		for ( Size i = begin; i <= end; ++i ) set_bb(i, false);
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @details all internal map data are cleared.
void
MoveMap::clear()
{
	torsion_type_map_.clear();
	move_map_torsion_id_map_.clear();
	torsion_id_map_.clear();

	dof_type_map_.clear();
	dof_id_map_.clear();

	jump_id_map_.clear();
}

/// set/get for JumpIDs --- fold-tree independent definition of jumps
void
MoveMap::set_jump( id::JumpID const & jump, bool const setting ) {
	jump_id_map_[ jump ] = setting;
}

///////////////////////////////////////////////////////////////////////////////
/// @details Set a specific TorsionType movable or not, e.g., "CHI".
/// Setting this TorsionType will clear data for individual MoveMapTorsionID and
/// TorsionID with this TorsionType to keep these three maps in sync.  Then query
/// for a specific TorsionID or MoveMapTorsionID will turn to setting for TorsionType.
void
MoveMap::set( TorsionType const & t, bool const setting )
{
	torsion_type_map_[t] = setting;

	{ // map by movemaptorsionid
		std::vector< MoveMapTorsionID > l;

		for ( MoveMapTorsionID_Map::const_iterator
				it=move_map_torsion_id_map_.begin(), it_end = move_map_torsion_id_map_.end(); it != it_end; ++it ) {
			// bad -- assumes knowledge of MoveMapTorsionID implementation as std::pair
			// but a whole object is probably overkill
			// compiler will catch type mismatch here
			if ( it->first.second == t ) l.push_back( it->first );
		}

		for ( std::vector< MoveMapTorsionID >::const_iterator it=l.begin(), it_end=l.end(); it != it_end; ++it ) {
			move_map_torsion_id_map_.erase( move_map_torsion_id_map_.find( *it ) );
		}
	}

	{ // map by torsionid
		std::vector< TorsionID > l;
		for ( TorsionID_Map::const_iterator
				it=torsion_id_map_.begin(), it_end = torsion_id_map_.end(); it != it_end; ++it ) {
			if ( it->first.type() == t ) l.push_back( it->first );
		}

		for ( std::vector< TorsionID >::const_iterator it=l.begin(), it_end=l.end(); it != it_end; ++it ) {
			torsion_id_map_.erase( torsion_id_map_.find( *it ) );
		}
	}

	// clear the map of JumpIDs
	if ( t == id::JUMP ) {
		jump_id_map_.clear();
	}
}

///////////////////////////////////////////////////////////////////////////////
// Set TorsionType flexible or fixed for one residue, e.g., BB torsions for residue 10
/// @details Setting this MoveMapTorsionID will clear data for individual TorsionID for this residue
/// with this TorsionType.
void
MoveMap::set( MoveMapTorsionID const & id, bool const setting )
{
	move_map_torsion_id_map_[ id ] = setting;

	// map by TorsionID
	Size const seqpos( id.first );
	TorsionType const & torsion_type( id.second );
	std::vector< TorsionID > l;
	for ( TorsionID_Map::const_iterator
			it=torsion_id_map_.begin(), it_end = torsion_id_map_.end(); it != it_end; ++it ) {
		if ( it->first.type() == torsion_type && it->first.rsd() == seqpos ) l.push_back( it->first );
	}

	for ( std::vector< TorsionID >::const_iterator it=l.begin(), it_end=l.end(); it != it_end; ++it ) {
		torsion_id_map_.erase( torsion_id_map_.find( *it ) );
	}
}

///////////////////////////////////////////////////////////////////////////////
// Set an individual Torsion movable for now, e.g., "BB torsion 2 of residue 4"
void
MoveMap::set( TorsionID const & id, bool const setting )
{
	torsion_id_map_[ id ] = setting;
}

///////////////////////////////////////////////////////////////////////////////
// the next two routines are identical between MoveMap and DOF_ID_Mask
//
// should think about how these two classes are related
//
// the thing is, the interface to this common data (dof-maps) is slightly
// different on the access side...
//
// set for this type of DOF, eg "PHI"
/// @details Setting this DOF type will also clear setting for individual DOF_ID of
/// this type in order to keep these two maps in sync, (query DOF_ID setting
/// will turn to DOF_type setting now)
void
MoveMap::set(
	DOF_Type const & t,
	bool const setting
)
{
	dof_type_map_[ t ] = setting;

	// now have to erase all the individual settings of this type, since
	// they are obliterated by this call
	std::vector< DOF_ID > l;
	for ( DOF_ID_Map::const_iterator it=dof_id_map_.begin(),
			it_end=dof_id_map_.end(); it != it_end; ++it ) {
		if ( it->first.type() == t ) {
			l.push_back( it->first );
		}
	}

	for ( std::vector< DOF_ID >::const_iterator it=l.begin(), it_end=l.end();
			it != it_end; ++it ) {
		dof_id_map_.erase( dof_id_map_.find( *it ) );
	}

}

///////////////////////////////////////////////////////////////////////////////
// set for an individual DoF, e.g., "PHI of Atom 3 in Residue 5"
void
MoveMap::set(
	DOF_ID const & id,
	bool const setting
)
{
	if ( get(id) == setting ) return;

	dof_id_map_[ id ] = setting;
}

bool
MoveMap::get_jump( id::JumpID const & jump ) const {
	JumpID_Map::const_iterator jid = jump_id_map_.find( jump );
	if ( jid != jump_id_map_.end() ) {
		return jid->second;
	} else {  //return global setting
		return get( id::JUMP );
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @details get setting for a specific TorsionType, such as "BB"
/// return false if no setting has been set for this TorsionType.
bool
MoveMap::get( TorsionType const & t ) const
{
	TorsionTypeMap::const_iterator i = torsion_type_map_.find( t );
	if ( i == torsion_type_map_.end() ) {
		return false;
	}
	return i->second;
}

///////////////////////////////////////////////////////////////////////////////
/// @details get TorsionType flexible or fixed for one residue, eg BB torsions for residue 10
/// if no setting for this MoveMapTorsionID, use the setting for the TorsionType
/// to which this MoveMapTorsionID belongs to.
bool
MoveMap::get( MoveMapTorsionID const & id ) const
{
	MoveMapTorsionID_Map::const_iterator i = move_map_torsion_id_map_.find( id );
	if ( i == move_map_torsion_id_map_.end() ) {
		TorsionType const & t( id.second );
		return get( t );
	}
	return i->second;
}

///////////////////////////////////////////////////////////////////////////////
/// @details get an individual torsion movable or not, eg BB torsion 2 of residue 4
/// if no setting for this specific TorsionID, use setting for the MoveMapTorsionID.
/// eg, no setting for BB torsion 2 of residue 4, check setting for BB torsions of
/// residue 4; if still not set, use setting for BB torsion type.
bool
MoveMap::get( TorsionID const & id ) const
{
	if ( !id.valid() ) return false;

	TorsionID_Map::const_iterator i = torsion_id_map_.find( id );
	if ( i == torsion_id_map_.end() ) {
		MoveMapTorsionID const move_map_torsion_id( id.rsd(), id.type() );
		return get( move_map_torsion_id );
	}
	return i->second;
}

/// @details get setting for this type of DOF, eg "PHI"
/// return false if no setting has been set to this DOF_type
bool
MoveMap::get(
	DOF_Type const & type
) const
{
	DOF_TypeMap::const_iterator iter( dof_type_map_.find( type ) );
	if ( iter == dof_type_map_.end() ) {
		return false;
	}
	return iter->second;
}

///////////////////////////////////////////////////////////////////////////////
/// @details get the setting for an individual dof, eg, PHI of Atom 3 in Residue 5
/// if no setting for this specific DOF_ID, get setting for the DOF type to
/// which this DOF_ID belongs to.
bool
MoveMap::get( DOF_ID const & id ) const
{
	if ( !id.valid() ) return false;

	DOF_ID_Map::const_iterator iter( dof_id_map_.find( id ) );
	if ( iter == dof_id_map_.end() ) {
		return get( id.type() );
	}
	return iter->second;
}

/// @brief find the explicit setting for the given TorsionType
/// @return iterator pointing to the TorsionType-bool pair, otherwise
///  torsion_type_end()
/// @warning Do not use this for general lookup, as it does not take
///  into account the stringency levels.  Only use this when you need
///  to check if a setting explicitly exists.
MoveMap::TorsionTypeMap::const_iterator
MoveMap::find( TorsionType const & t ) const {
	return torsion_type_map_.find( t );
}

/// @brief find the explicit setting for the given MoveMapTorsionID
/// @return iterator pointing to the MoveMapTorsionID-bool pair, otherwise
///  movemap_torsion_id_end()
/// @warning Do not use this for general lookup, as it does not take
///  into account the stringency levels.  Only use this when you need
///  to check if a setting explicitly exists.
MoveMap::MoveMapTorsionID_Map::const_iterator
MoveMap::find( MoveMapTorsionID const & id ) const {
	return move_map_torsion_id_map_.find( id );
}


/// @brief find the explicit setting for the given TorsionID
/// @return iterator pointing to the TorsionID-bool pair, otherwise torsion_id_end()
/// @warning Do not use this for general lookup, as it does not take
///  into account the stringency levels.  Only use this when you need
///  to check if a setting explicitly exists.
MoveMap::TorsionID_Map::const_iterator
MoveMap::find( TorsionID const & id ) const {
	return torsion_id_map_.find( id );
}

/// @brief find the explicit setting for the given JumpID
/// @return iterator pointing to the JumpID-bool pair, otherwise jump_id_end()
/// @warning Do not use this for general lookup, as it does not take
///  into account the stringency levels.  Only use this when you need
///  to check if a setting explicitly exists.
MoveMap::JumpID_Map::const_iterator
MoveMap::find( id::JumpID const & jump ) const {
	return jump_id_map_.find( jump );
}

/// @brief find the explicit setting for the given DOF_Type
/// @return iterator pointing to the DOF_Type-bool pair, otherwise dof_type_end()
/// @warning Do not use this for general lookup, as it does not take
///  into account the stringency levels.  Only use this when you need
///  to check if a setting explicitly exists.
MoveMap::DOF_TypeMap::const_iterator
MoveMap::find( DOF_Type const & t ) const {
	return dof_type_map_.find( t );
}

/// @brief find the explicit setting for the given DOF_ID
/// @return iterator pointing to the DOF_ID-bool pair, otherwise dof_id_end()
/// @warning Do not use this for general lookup, as it does not take
///  into account the stringency levels.  Only use this when you need
///  to check if a setting explicitly exists.
MoveMap::DOF_ID_Map::const_iterator
MoveMap::find( DOF_ID const & id ) const {
	return dof_id_map_.find( id );
}


/// @brief reads lines of format and set movemap
///    RESIDUE * CHI        # set all residues chi movable
///    RESIDUE 36 48 BBCHI  # set res 36-48 bb & chi movable
///    RESIDUE 89 NO        # set res 89 unmovable
///    JUMP * NO            # set all jumps unmovable
///    JUMP 1 YES           # set jump 1 movable
/// If a residue/default is not specified, mm defaults to current value.
/// If a value for a jum is not given (e.g. "JUMP 4\n"), it defaults to movable (YES)
/// Setting 'CHI' implies BB not movable, thus don't do:
///    RESIDUE * CHI
///    RESIDUE * BB
/// Instead:
///    RESIDUE * BBCHI
void
MoveMap::init_from_file( std::string const & filename ) {
	utility::io::izstream data( filename.c_str() );
	if ( !data ) {
		utility_exit_with_message("ERROR: could not open file " + filename );
	}

	std::string line;
	bool res_default_set=false, jump_default_set=false;
	while ( getline(data,line) ) {
		if ( line.substr(0,1) == "#" ) continue;
		if ( line == "" ) continue;

		utility::vector1< std::string > tokens ( utility::split( line ) );
		if ( tokens.size() > 0 ) {
			if ( tokens[1] == "JUMP" || tokens[1] == "jump" ) {
				if ( tokens.size() < 2 ) {
					utility_exit_with_message( "Error reading movemap at line: " + line );
				}

				if ( tokens[2] == "*" ) {
					if ( jump_default_set ) {
						utility_exit_with_message( "Error reading movemap: default jump set multiple times!" );
					} else {
						jump_default_set=true;
						if ( tokens.size() < 3 || tokens[3] == "YES" || tokens[3] == "yes" ) {
							set_jump( true );
						} else if ( tokens[3] == "NO" || tokens[3] == "no" ) {
							set_jump( false );
						} else {
							utility_exit_with_message( "Error reading movemap at line: " + line );
						}
					}
				} else {
					core::Size jump_num = (core::Size) atoi(tokens[2].c_str());
					if ( tokens.size() < 3 || tokens[3] == "YES" || tokens[3] == "yes" ) {
						set_jump( jump_num, true );
					} else if ( tokens[3] == "NO" || tokens[3] == "no" ) {
						set_jump( jump_num, false );
					} else {
						utility_exit_with_message( "Error reading movemap at line: " + line );
					}
				}
			} else if ( tokens[1] == "RESIDUE" || tokens[1] == "residue" ) {
				if ( tokens.size() < 3 ) {
					utility_exit_with_message( "Error reading movemap at line: " + line );
				}
				if ( tokens[2] == "*" ) {
					if ( res_default_set ) {
						utility_exit_with_message( "Error reading movemap: default residue set multiple times!" );
					} else {
						res_default_set=true;
						if ( tokens[3] == "BB" || tokens[3] == "bb" ) {
							set_bb( true ); set_chi( false );
						} else if ( tokens[3] == "CHI" || tokens[3] == "chi" ) {
							set_bb( false ); set_chi( true );
						} else if ( tokens[3] == "BBCHI" || tokens[3] == "bbchi" ) {
							set_bb( true ); set_chi( true );
						} else if ( tokens[3] == "NO" || tokens[3] == "no" ) {
							set_bb( false ); set_chi( false );
						} else {
							utility_exit_with_message( "Error reading movemap at line: " + line );
						}
					}
				} else if ( tokens.size() == 3 ) {
					core::Size start_res   = (core::Size) atoi(tokens[2].c_str());
					if ( tokens[3] == "BB" || tokens[3] == "bb" ) {
						set_bb( start_res, true ); set_chi( start_res, false );
					} else if ( tokens[3] == "CHI" || tokens[3] == "chi" ) {
						set_bb( start_res, false ); set_chi( start_res, true );
					} else if ( tokens[3] == "BBCHI" || tokens[3] == "bbchi" ) {
						set_bb( start_res, true ); set_chi( start_res, true );
					} else if ( tokens[3] == "NO" || tokens[3] == "no" ) {
						set_bb( start_res, false ); set_chi( start_res, false );
					} else {
						utility_exit_with_message( "Error reading movemap at line: " + line );
					}
				} else {
					core::Size start_res = (core::Size) atoi(tokens[2].c_str());
					core::Size end_res = (core::Size) atoi(tokens[3].c_str());
					if ( start_res > end_res ) {
						utility_exit_with_message( "Error reading movemap at line: " + line );
					}
					if ( tokens[4] == "BB" || tokens[4] == "bb" ) {
						for ( core::Size i = start_res; i<= end_res; ++i ) {
							set_bb( i, true ); set_chi( i, false );
						}
					} else if ( tokens[4] == "CHI" || tokens[4] == "chi" ) {
						for ( core::Size i = start_res; i<= end_res; ++i ) {
							set_bb( i, false ); set_chi( i, true );
						}
					} else if ( tokens[4] == "BBCHI" || tokens[4] == "bbchi" ) {
						for ( core::Size i = start_res; i<= end_res; ++i ) {
							set_bb( i, true ); set_chi( i, true );
						}
					} else if ( tokens[4] == "NO" || tokens[4] == "no" ) {
						for ( core::Size i = start_res; i<= end_res; ++i ) {
							set_bb( i, false ); set_chi( i, false );
						}
					} else {
						utility_exit_with_message( "Error reading movemap at line: " + line );
					}
				}
			} else {
				utility_exit_with_message( "Error reading movemap at line: " + line );
			}
		}
	}

} // init_from_file

void
MoveMap::show( std::ostream & out, Size n_residues_to_show ) const
{
	PyAssert( (n_residues_to_show>0), "MoveMap::show( std::ostream & out , Size n_residues_to_show ): "
		"input variable total_residue has a meaningless value");
	out << A(8, "resnum") << ' ';
	out << A(8, "BB") << ' ' << A(8, "CHI") << ' ' << A(8, "NU") << A(8, "BRANCH") << std::endl;
	for ( Size i = 1; i <= n_residues_to_show; ++i ) {
		std::string bb = "FALSE";
		std::string chi = "FALSE";
		std::string nu = "FALSE";
		std::string branches = "FALSE";
		if ( get_bb( i ) ) { bb = "TRUE "; }
		if ( get_chi( i ) ) { chi = "TRUE "; }
		if ( get_nu( i ) ) { nu = "TRUE "; }
		if ( get_branches( i ) ) { branches = "TRUE "; }
		out << I(8,3,i) << ' ' << A(8, bb) << ' ' << A(8, chi) << A(8, nu) << A(8, branches) << std::endl;
	}
}

void
MoveMap::show( std::ostream & out ) const
{
	out << "\n";
	out << "-------------------------------\n";
	out << A(8, "resnum") << ' ' << A(8, "Type") << ' ' << A(12, "TRUE/FALSE ") << "\n";
	out << "-------------------------------\n";
	// The general settings:
	out << A(8,"DEFAULT") <<' '<< A(7, id::to_string(id::BB) ) << "  " << A(8,( get(id::BB) ? "TRUE":"FALSE")) << "\n";
	out << A(8,"DEFAULT") <<' '<< A(7, id::to_string(id::CHI)) << "  " << A(8,( get(id::CHI) ? "TRUE":"FALSE")) << "\n";
	out << A(8,"DEFAULT") <<' '<< A(7, id::to_string(id::NU) ) << "  " << A(8,( get(id::NU) ? "TRUE":"FALSE")) << "\n";
	out << A(8,"DEFAULT") <<' '<< A(7, id::to_string(id::BRANCH) ) << "  "
		<< A(8,( get(id::BRANCH) ? "TRUE":"FALSE")) << "\n";
	// The overrides:
	Size prev_resnum = 0;
	utility::vector1< bool > jumpbool;
	utility::vector1< Size > jumpnum;
	for ( MoveMapTorsionID_Map::const_iterator it = movemap_torsion_id_begin(), it_end = movemap_torsion_id_end();
			it != it_end; ++it ) {
		MoveMapTorsionID mmtorsionID = it->first;
		bool boolean = it->second;
		Size res = mmtorsionID.first;
		TorsionType torsiontype = mmtorsionID.second;
		std::string type( id::to_string( torsiontype ) );

		// Jumps are handled under a separate heading.
		if ( torsiontype == id::JUMP ) {
			jumpbool.push_back(boolean);
			jumpnum.push_back(res);
			continue;
		}

		// Only show each residue once (and only if torsion type is BB, SC, or NU)
		if ( prev_resnum != mmtorsionID.first ) {
			out << I(8,3,res) <<' '<< A(7,type) << "  " << A(8, (boolean ? "TRUE":"FALSE")) << "\n";
		} else {
			out << A(8,' ') << ' ' << A(7,type) << "  " << A(8, (boolean ? "TRUE":"FALSE")) << "\n";
		}

		// Remember the previous residue/jump number
		prev_resnum = res;
	}
	out << "-------------------------------\n";
	out << A(8, "jumpnum") << ' ' << A(8, "Type") << ' ' << A(12, "TRUE/FALSE ") << "\n";
	out << "-------------------------------\n";
	// The general setting
	out << A(8,"DEFAULT")<<' '<< A(8,"JUMP") <<' '<< A(8,( get(id::JUMP) ? "TRUE":"FALSE"))<< "\n";
	// Jump overrides
	for ( Size i = 1; i <= jumpnum.size(); ++i ) {
		out << I(8,3,jumpnum[i])<<' '<< A(8,"JUMP") <<' '<< A(8,(jumpbool[i] ? "TRUE":"FALSE"))<< "\n";
	}

	out << "-------------------------------\n";
	out << A(8, "resnum") << ' ' << A(8, "atomnum") << ' ' << A(8, "Type") << ' ' << A(12, "TRUE/FALSE ") << "\n";
	out << "-------------------------------\n";
	// The defaults
	out << A(8,"DEFAULT") << ' ' << A(8, ' ') << ' ' << A(8,id::to_string(id::PHI)) <<' '<< A(8,( get(id::PHI) ? "TRUE":"FALSE"))<< "\n";
	out << A(8,"DEFAULT") << ' ' << A(8, ' ') << ' ' << A(8,id::to_string(id::THETA)) <<' '<< A(8,( get(id::THETA) ? "TRUE":"FALSE"))<< "\n";
	out << A(8,"DEFAULT") << ' ' << A(8, ' ') << ' ' << A(8,id::to_string(id::D)) <<' '<< A(8,( get(id::D) ? "TRUE":"FALSE"))<< "\n";
	out << A(8,"DEFAULT") << ' ' << A(8, ' ') << ' ' << A(8,id::to_string(id::RB1)) <<' '<< A(8,( get(id::RB1) ? "TRUE":"FALSE"))<< "\n";
	out << A(8,"DEFAULT") << ' ' << A(8, ' ') << ' ' << A(8,id::to_string(id::RB2)) <<' '<< A(8,( get(id::RB2) ? "TRUE":"FALSE"))<< "\n";
	out << A(8,"DEFAULT") << ' ' << A(8, ' ') << ' ' << A(8,id::to_string(id::RB3)) <<' '<< A(8,( get(id::RB3) ? "TRUE":"FALSE"))<< "\n";
	out << A(8,"DEFAULT") << ' ' << A(8, ' ') << ' ' << A(8,id::to_string(id::RB4)) <<' '<< A(8,( get(id::RB4) ? "TRUE":"FALSE"))<< "\n";
	out << A(8,"DEFAULT") << ' ' << A(8, ' ') << ' ' << A(8,id::to_string(id::RB5)) <<' '<< A(8,( get(id::RB5) ? "TRUE":"FALSE"))<< "\n";
	out << A(8,"DEFAULT") << ' ' << A(8, ' ') << ' ' << A(8,id::to_string(id::RB6)) <<' '<< A(8,( get(id::RB6) ? "TRUE":"FALSE"))<< "\n";
	prev_resnum = 0;
	for ( DOF_ID_Map::const_iterator it = dof_id_begin(), it_end = dof_id_end();
			it != it_end; ++it ) {
		DOF_ID const & dofID = it->first;
		bool boolean = it->second;
		Size res = dofID.rsd();
		Size atomno = dofID.atomno();
		DOF_Type doftype = dofID.type();
		std::string type( id::to_string( doftype ) );
		out << I(8,3,res) << ' ' << I(8,3,atomno) << ' ' << A(8,type) <<' '<< A(8,( boolean ? "TRUE":"FALSE"))<< "\n";
	}
	out << std::endl;
}

/// @brief import settings from another MoveMap
/// @param[in] rval The MoveMap to import settings from.
/// @param[in] import_true_settings Import True settings?
/// @param[in] import_false_settings Import False settings?
/// @return The total number of settings imported.
/// @remarks This function calls set() for each setting that exists in the
///  'rval' MoveMap in order from lowest to highest stringency.
Size MoveMap::import(
	MoveMap const & rval,
	bool const import_true_settings,
	bool const import_false_settings
)
{
	// Import settings from the other movemap.  Always add in order from
	// lowest stringency to highest.
	Size n = 0;

	// Step 1: torsions
	// TorsionType
	for ( TorsionTypeMap::const_iterator i = rval.torsion_type_begin(), ie = rval.torsion_type_end(); i != ie; ++i ) {
		if ( ( import_true_settings && i->second ) || ( import_false_settings && !i->second ) ) {
			set( i->first, i->second );
			++n;
		}
	}

	// MoveMapTorsionID
	for ( MoveMapTorsionID_Map::const_iterator i = rval.movemap_torsion_id_begin(), ie = rval.movemap_torsion_id_end(); i != ie; ++i ) {
		if ( ( import_true_settings && i->second ) || ( import_false_settings && !i->second ) ) {
			set( i->first, i->second );
			++n;
		}
	}

	// TorsionID
	for ( TorsionID_Map::const_iterator i = rval.torsion_id_begin(), ie = rval.torsion_id_end(); i != ie; ++i ) {
		if ( ( import_true_settings && i->second ) || ( import_false_settings && !i->second ) ) {
			set( i->first, i->second );
			++n;
		}
	}

	// Step 2: DOFs
	// DOF_Type
	for ( DOF_TypeMap::const_iterator i = rval.dof_type_begin(), ie = rval.dof_type_end(); i != ie; ++i ) {
		if ( ( import_true_settings && i->second ) || ( import_false_settings && !i->second ) ) {
			set( i->first, i->second );
			++n;
		}
	}

	// DOF_ID
	for ( DOF_ID_Map::const_iterator i = rval.dof_id_begin(), ie = rval.dof_id_end(); i != ie; ++i ) {
		if ( ( import_true_settings && i->second ) || ( import_false_settings && !i->second ) ) {
			set( i->first, i->second );
			++n;
		}
	}

	// Step 3: jumps
	// JumpID
	for ( JumpID_Map::const_iterator i = rval.jump_id_begin(), ie = rval.jump_id_end(); i != ie; ++i ) {
		if ( ( import_true_settings && i->second ) || ( import_false_settings && !i->second ) ) {
			set_jump( i->first, i->second );
			++n;
		}
	}

	return n;
}

} // namespace kinematics
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::MoveMap::save( Archive & arc ) const {
	arc( CEREAL_NVP( torsion_type_map_ ) ); // TorsionTypeMap
	arc( CEREAL_NVP( move_map_torsion_id_map_ ) ); // MoveMapTorsionID_Map
	arc( CEREAL_NVP( torsion_id_map_ ) ); // TorsionID_Map
	arc( CEREAL_NVP( dof_type_map_ ) ); // DOF_TypeMap
	arc( CEREAL_NVP( dof_id_map_ ) ); // DOF_ID_Map
	arc( CEREAL_NVP( jump_id_map_ ) ); // JumpID_Map
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::MoveMap::load( Archive & arc ) {
	arc( torsion_type_map_ ); // TorsionTypeMap
	arc( move_map_torsion_id_map_ ); // MoveMapTorsionID_Map
	arc( torsion_id_map_ ); // TorsionID_Map
	arc( dof_type_map_ ); // DOF_TypeMap
	arc( dof_id_map_ ); // DOF_ID_Map
	arc( jump_id_map_ ); // JumpID_Map
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::MoveMap );
CEREAL_REGISTER_TYPE( core::kinematics::MoveMap )

CEREAL_REGISTER_DYNAMIC_INIT( core_kinematics_MoveMap )
#endif // SERIALIZATION
