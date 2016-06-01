// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  Class to store symmetry data specified from input files
/// @file   core/conformation/symmetry/SymmData.cc
/// @author Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/VirtualCoordinate.hh>
#include <core/conformation/symmetry/SymDof.hh>

// Project headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/xyzTransform.hh>

// C++ headers
#include <iostream>

// Utility header
#include <utility/string_util.hh>

#include <utility/vector1.hh>
#include <fstream>
#include <iomanip>

// @details: This data structure contains all data necessary to generate a
// symmetrical system, score it and move it. It is used to initialize the
// SymmetryInfo object in SymmetricalConformation that carries the symmetry
// information pose need.

namespace core {
namespace conformation {
namespace symmetry {

using std::map;
using std::endl;
using std::string;
using std::pair;
using std::make_pair;
using std::vector;
using utility::vector1;
typedef numeric::xyzTransform<Real> Xform;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

//typedef vector1< Size > Clones;

typedef vector1< pair<Size,Real> > WtedClones;

static THREAD_LOCAL basic::Tracer TR( "core.conformation.symmetry.SymmData" );

static std::string const NOPARENT = "NOPARENT";

// @define: constructor
SymmData::SymmData() :
	utility::pointer::ReferenceCount(),
	subunits_( 0 ),
	num_components_( 1 ),
	interfaces_( 1 ),
	score_subunit_( 1 ),
	anchor_residue_( "1" ),
	recenter_( false ),
	root_( 1 )
{}

// @define: constructor
SymmData::SymmData( Size, Size ) :
	utility::pointer::ReferenceCount(),
	subunits_( 0 ),
	num_components_( 1 ),
	interfaces_( 1 ),
	score_subunit_( 1 ),
	anchor_residue_( "1" ),
	recenter_( false ),
	root_( 1 )
{}


// @define: constructor
SymmData::SymmData( SymmData const & src ) :
	utility::pointer::ReferenceCount(),
	symmetry_name_( src.symmetry_name_ ),
	symmetry_type_( src.symmetry_type_ ),
	subunits_( src.subunits_ ),
	num_components_( src.num_components_ ),
	interfaces_( src.interfaces_ ),
	score_subunit_( src.score_subunit_ ),
	anchor_residue_( src.anchor_residue_ ),
	recenter_( src.recenter_ ),
	root_(src.root_),
	slide_info_(src.slide_info_),
	slide_order_string_(src.slide_order_string_),
	symm_transforms_( src.symm_transforms_ ),
	rotation_matrices_( src.rotation_matrices_ ),
	translation_matrices_( src.translation_matrices_ ),
	virtual_coordinates_( src.virtual_coordinates_ ),
	jump_string_to_virtual_pair_( src.jump_string_to_virtual_pair_ ),
	jump_string_to_jump_num_( src.jump_string_to_jump_num_ ),
	virt_id_to_virt_num_( src.virt_id_to_virt_num_ ),
	virt_id_to_subunit_num_( src.virt_id_to_subunit_num_ ),
	virt_id_to_subunit_chain_( src.virt_id_to_subunit_chain_ ),
	virt_id_to_subunit_residue_( src.virt_id_to_subunit_residue_ ),
	virt_num_to_virt_id_( src.virt_num_to_virt_id_ ),
	subunit_num_to_virt_id_( src.subunit_num_to_virt_id_ ),
	jump_clones_( src.jump_clones_ ),
	dofs_( src.dofs_ ),
	allow_virtual_( src.allow_virtual_ ),
	score_multiply_subunit_( src.score_multiply_subunit_ ),
	include_subunit_( src.include_subunit_ ),
	output_subunit_( src.output_subunit_ ),
	components_(src.components_),
	name2component_(src.name2component_),
	jname2components_(src.jname2components_),
	jname2subunits_(src.jname2subunits_)
{}

// @define: deep copy of SymmData
SymmDataOP
SymmData::clone() const
{
	return SymmDataOP( new SymmData( *this ) );
}

// Accessor functions

//Size
//SymmData::get_nres_subunit() const
//{
//return nres_subunit_;
//}

//Size
//SymmData::get_njump_subunit() const
//{
//return njump_subunit_;
//}

string const &
SymmData::get_symmetry_name() const
{
	return symmetry_name_;
}

string const &
SymmData::get_symmetry_type() const
{
	return symmetry_type_;
}

core::Size
SymmData::get_subunits() const
{
	return subunits_;
}

core::Size
SymmData::get_num_components() const
{
	return num_components_;
}

core::Size
SymmData::get_interfaces() const
{
	return interfaces_;
}

core::Size
SymmData::get_score_subunit() const
{
	return score_subunit_;
}

string const &
SymmData::get_anchor_residue() const
{
	return anchor_residue_;
}

bool
SymmData::get_recenter() const
{
	return recenter_;
}

core::Size
SymmData::get_root() const
{
	return root_;
}

vector1< Size > const &
SymmData::get_score_multiply_subunit() const
{
	return score_multiply_subunit_;
}

vector1< Size > const &
SymmData::get_include_subunit() const
{
	return include_subunit_;
}

vector1< Size > const &
SymmData::get_output_subunit() const
{
	return output_subunit_;
}

vector< numeric::xyzMatrix< core::Real > > const &
SymmData::get_rotation_matrix() const
{
	return rotation_matrices_;
}

vector< numeric::xyzMatrix< core::Real > > const &
SymmData::get_translation_matrix() const
{
	return translation_matrices_;
}

map< string, VirtualCoordinate > const &
SymmData::get_virtual_coordinates() const
{
	return virtual_coordinates_;
}

core::Size
SymmData::get_num_virtual() const
{
	return virtual_coordinates_.size();
}

map< Size, SymDof > const &
SymmData::get_dofs() const
{
	return dofs_;
}

map< Size, WtedClones > const &
SymmData::get_jump_clones() const
{
	return jump_clones_;
}

map< string, Size > const &
SymmData::get_jump_string_to_jump_num() const
{
	return jump_string_to_jump_num_;
}

map< string, Size > const &
SymmData::get_virtual_id_to_num() const
{
	return virt_id_to_virt_num_;
}

map< string, Size > const &
SymmData::get_virt_id_to_subunit_num() const
{
	return virt_id_to_subunit_num_;
}

map< string, char > const &
SymmData::get_virt_id_to_subunit_chain() const
{
	return virt_id_to_subunit_chain_;
}

map< string, string > const &
SymmData::get_virt_id_to_subunit_residue() const
{
	return virt_id_to_subunit_residue_;
}

map< Size, string > const &
SymmData::get_subunit_num_to_virt_id() const
{
	return subunit_num_to_virt_id_;
}

map< Size, string > const &
SymmData::get_virtual_num_to_id() const
{
	return virt_num_to_virt_id_;
}

map< string, pair< string, string > > const &
SymmData::get_virtual_connects() const
{
	return jump_string_to_virtual_pair_;
}

SymSlideInfo const &
SymmData::get_slide_info() const
{
	return slide_info_;
}

void
SymmData::set_symmetry_name(
	string symm_name )
{
	symmetry_name_ = symm_name;
}

void
SymmData::set_symmetry_type(
	string symm_type )
{
	symmetry_type_ = symm_type;
}

void
SymmData::set_subunits(
	core::Size num_subunits )
{
	subunits_ = num_subunits;
}

void
SymmData::set_interfaces(
	core::Size interfaces )
{
	interfaces_ = interfaces;
}

void
SymmData::set_anchor_residue( string anchor )
{
	anchor_residue_ = anchor;
}

void
SymmData::set_score_multiply_subunit( vector1< Real > & score_multiply_subunit_vector )
{
	score_multiply_subunit_ = score_multiply_subunit_vector;
}

void
SymmData::set_slide_info( SymSlideInfo slide_info )
{
	slide_info_ = slide_info;
}

void
SymmData::set_rotation_matrix(
	vector< numeric::xyzMatrix< core::Real > > rotation_matrices )
{
	rotation_matrices_ = rotation_matrices;
}

void
SymmData::set_translation_matrix(
	vector< numeric::xyzMatrix< core::Real > > translation_matrices )
{
	translation_matrices_ = translation_matrices;
}

void
SymmData::set_symm_transforms(
	vector< vector< string> > symm_transforms )
{
	symm_transforms_ = symm_transforms;
}

// @define: destructor
SymmData::~SymmData(){}


// @details: Parse symmetry information from a textfile. This function fills all data necessary to generate a
// symmetrical system, score it and move it.
void
SymmData::read_symmetry_data_from_file(
	string filename
)
{
	std::ifstream infile( filename.c_str() );

	if ( !infile.good() ) {
		utility_exit_with_message( "[ERROR] Error opening symmetry file '" + filename + "'" );
	}

	read_symmetry_data_from_stream(infile);
}

void
SymmData::read_symmetry_data_from_stream(
	std::istream & infile
)
{
	string line;
	int linecount( 0 );
	bool read_virtual_coords (false);         // Find out when to start reading virtual coordinates
	bool read_transforms (false);             // Find out when to read symmetry transforms
	bool start_coordinates_found (false);     // Find out when start coordinates are found. Must be used
	// with read_transform
	bool have_virtual_coordinates( false );
	bool read_consecutive_transforms(false);
	bool read_subunits(false);
	bool set_jump_numbers_from_tags(true);
	bool connect_virtuals_specified(false);

	vector< Matrix > rot_matrix;                  // store rotation matrices during sym transform reading
	vector< Vector > trans_vector;                // store translation matrices during sym transform reading
	vector< pair< Size, string > > transform_type; // store type of symmetry transform. Rotation or translation
	Size num_transformations( 0 );                     // Number of sym transforms read
	core::Size N( 1 );                                 // Number of subunits in the system
	vector1< string > score_multiply_subunit_string;
	std::set<char> subchains;

	// Read the file
	while ( getline(infile,line) ) {
		linecount++;
		line = line.substr(0,line.find('#'));
		vector1< string > tokens( utility::split( line ) );
		if ( tokens.size()==0 ) continue;

		// if we're reading a coordinate block ...
		if ( read_virtual_coords ) {
			if ( tokens[1] == "xyz" || tokens[1] == "xyzM" ) {
				// Lines of virtual coordinate systems start with "xyz" and occur in
				//    virtual_coordinates_start ... virtual_coordinates_stop blocks
				// parse the coordinate system into a VirtualCoordinate object. Need to
				// have at least x and y axis. If no origin is specified we assume 0,0,0
				if ( tokens.size() < 5 ) {
					utility_exit_with_message("[ERROR] xyz lines in symm def file const contain an identifier and 3 coordinate vectors (X,Y,ORIG)" );
				}
				string identifier( tokens[2] );
				VirtualCoordinate coord;
				coord.set_mirror_z( tokens[1] == "xyzM" );
				coord.add_coordinate_from_string( tokens, 3 );

				if ( virtual_coordinates_.find(identifier) != virtual_coordinates_.end() ) {
					utility_exit_with_message("[ERROR] VRT identifier "+identifier+" defined twice in symmetry definition file!" );
				}

				virtual_coordinates_[identifier] = coord;
				Size jump_num (virtual_coordinates_.size());

				// vrt id to idx maps
				virt_id_to_virt_num_[identifier] = jump_num;
				virt_num_to_virt_id_[jump_num] = identifier;
			} else if ( tokens[1] == "virtual_coordinates_stop" ) {
				read_virtual_coords = false;
				have_virtual_coordinates = true;
			} else {
				utility_exit_with_message("[ERROR] Error reading symm def file while at '"+line+"'");
			}
		} else if ( read_transforms ) {
			// if we're reading a list of transforms
			// Read the start coordinate system upon which the first transformation operation work on.
			// Store the value in a VirtualCoordinate system
			if ( tokens[1] == "start" && tokens.size() >= 3 ) {
				string identifier = "VRT0001";
				VirtualCoordinate start_coord;
				start_coord.add_coordinate_from_string( tokens );
				virtual_coordinates_.insert( make_pair( identifier, start_coord ) );

				// The first one is always the master that controlls the clones
				string tag = "BASEJUMP";
				jump_string_to_jump_num_.insert( make_pair( tag, 1 ) );
				jump_string_to_virtual_pair_.insert( make_pair( tag, make_pair( identifier, "SUBUNIT" ) ) );
				virt_id_to_virt_num_.insert( make_pair( identifier, 1 ) );
				virt_num_to_virt_id_.insert( make_pair( 1, identifier ) );
				virt_id_to_subunit_num_.insert( make_pair( identifier, 1 ) );
				subunit_num_to_virt_id_.insert( make_pair( 1, identifier ) );
				start_coordinates_found = true;
			} else if ( tokens[1] == "rot" || tokens[1] == "trans" ) {
				if ( !start_coordinates_found ) {
					utility_exit_with_message( "[ERROR] When specifying Rotation operations a starting coordinate must be give" );
				}
				// Rotation around the x-axis
				if ( tokens[2] == "Rx" ) {
					// Check that we have both rot, Rx and a value specifying the number of rotations
					if ( tokens.size() != 3 ) {
						utility_exit_with_message( "[ERROR] Need to give two arguments..." );
					}
					// Read the N-fold rotation
					N = utility::string2int( tokens[2] );
					Vector axis (1,0,0);
					// Save the rotation matrix in the rot_matrix vector
					rot_matrix.push_back( numeric::rotation_matrix_degrees( axis, 360.0 / N ) );
					// We need to know what type of transformation this is
					transform_type.push_back( make_pair( ++num_transformations, "rot") );
				}
				if ( tokens[2] == "Ry" ) {
					if ( tokens.size() != 3 ) {
						utility_exit_with_message( "[ERROR] Need to give two arguments..." );
					}
					N = utility::string2int( tokens[2] );
					Vector axis (0,1,0);
					rot_matrix.push_back( numeric::rotation_matrix_degrees( axis, 360.0 / N ) );
					transform_type.push_back( make_pair( ++num_transformations, "rot") );
				}
				if ( tokens[2] == "Rz" ) {
					if ( tokens.size() != 3 ) {
						utility_exit_with_message( "[ERROR] Need to give two arguments..." );
					}
					N = utility::string2int( tokens[3] );
					Vector axis (0,0,1);
					rot_matrix.push_back( numeric::rotation_matrix_degrees( axis, 360.0 / N ) );
					transform_type.push_back( make_pair( ++num_transformations, "rot") );
				}
				// Rotate an angle "angle" around the x-axis
				if ( tokens[2] == "Rx_angle" ) {
					if ( tokens.size() != 3 ) {
						utility_exit_with_message( "[ERROR] Need to give two arguments..." );
					}
					Vector axis (1,0,0);
					Real angle = static_cast<core::Real>( std::atof( tokens[3].c_str() ) );
					rot_matrix.push_back( numeric::rotation_matrix_degrees( axis, angle ) );
					transform_type.push_back( make_pair( ++num_transformations, "rot") );
				}
				if ( tokens[2] == "Ry_angle" ) {
					if ( tokens.size() != 3 ) {
						utility_exit_with_message( "[ERROR] Need to give two arguments..." );
					}
					Vector axis (0,1,0);
					Real angle = static_cast<core::Real>( std::atof( tokens[3].c_str() ) );
					rot_matrix.push_back( numeric::rotation_matrix_degrees( axis, angle ) );
					transform_type.push_back( make_pair( ++num_transformations, "rot") );
				}
				if ( tokens[2] == "Rz_angle" ) {
					if ( tokens.size() != 3 ) {
						utility_exit_with_message( "[ERROR] Need to give two arguments..." );
					}
					Vector axis (0,0,1);
					Real angle = static_cast<core::Real>( std::atof( tokens[3].c_str() ) );
					rot_matrix.push_back( numeric::rotation_matrix_degrees( axis, angle ) );
					transform_type.push_back( make_pair( ++num_transformations, "rot") );
				}
				// read translations. Need to give a vector
				if ( tokens[1] == "trans" ) {
					vector1< string> split ( utility::string_split( tokens[2], ',' ) );
					runtime_assert( split.size() == 3 );
					Vector trans( ( static_cast<core::Real>( std::atof( split[1].c_str() ) ) ),
						( static_cast<core::Real>( std::atof( split[2].c_str() ) ) ),
						( static_cast<core::Real>( std::atof( split[3].c_str() ) ) ) );
					trans_vector.push_back( trans );
					transform_type.push_back( make_pair( ++num_transformations, "trans") );
				}
			} else if ( tokens[1] == "virtual_transforms_stop" ) {
				read_transforms = false;
				// For all the subunits generate a coordinate system. Go through all the tranformation
				// operations and operate on the previous coordinate system sequentially. The first
				// coordinate system is the start coordinates. All transformations are acted on for
				// virtual coordinate system. This may not make sense always, need to look at this more...
				Size index(0);
				Size num_rots(0), num_trans(0);
				// First find the start coordinates. They are labelled with VRT1
				if ( virtual_coordinates_.find( "VRT0001" ) == virtual_coordinates_.end()  ) {
					utility_exit_with_message( "[ERROR] start coordinates VRT0001 not found..." );
				}
				VirtualCoordinate virt_coord( virtual_coordinates_.find( "VRT0001" )->second );
				// read transforms in a consecutive manner
				if ( read_consecutive_transforms ) {
					for ( Size i=1; i < subunits_; ++i ) {
						Vector x_new( virt_coord.get_x() );
						Vector y_new( virt_coord.get_y() );
						Vector origin_new ( virt_coord.get_origin() );
						Size num_rots(0), num_trans(0);
						for ( vector< pair< Size, string > >::const_iterator it = transform_type.begin();
								it != transform_type.end(); ++it ) {
							if ( it->second == "rot" ) {
								Size matrix_num ( ++num_rots );
								x_new = rot_matrix[matrix_num -1]*x_new;
								y_new = rot_matrix[matrix_num -1]*y_new;
								origin_new = rot_matrix[matrix_num -1]*origin_new;
							} else if ( it->second == "trans" ) {
								Size vector_num ( ++num_trans );
								origin_new = trans_vector[vector_num -1] + origin_new;
							}
						}
						VirtualCoordinate coord( x_new, y_new, origin_new );
						Size jump_num ( virtual_coordinates_.size() + 1 );
						string tag;
						std::ostringstream stream;
						stream << std::setfill('0') << std::setw(4) << jump_num;
						tag = "CLONE" + stream.str();
						string identifier = "VRT" + stream.str();
						virtual_coordinates_.insert( make_pair( identifier, coord ) );
						jump_string_to_jump_num_.insert( make_pair( tag, jump_num ) );
						jump_string_to_virtual_pair_.insert( make_pair( tag, make_pair( identifier, "SUBUNIT" ) ) );
						virt_id_to_virt_num_.insert( make_pair( identifier, jump_num ) );
						virt_num_to_virt_id_.insert( make_pair( jump_num, identifier ) );
						virt_id_to_subunit_num_.insert( make_pair( identifier, virt_id_to_subunit_num_.size()+1 ) );
						subunit_num_to_virt_id_.insert( make_pair( subunit_num_to_virt_id_.size()+1, identifier ) );
						virt_coord = coord;
					}
				} else {
					for ( Size i=1; i < subunits_; ++i ) {
						if ( i > transform_type.size() ) {
							index = 0;
							num_rots = 0;
							num_trans = 0;
						}
						++index;
						Vector x_new( virt_coord.get_x() );
						Vector y_new( virt_coord.get_y() );
						Vector origin_new ( virt_coord.get_origin() );
						if ( transform_type[index-1].second == "rot" ) {
							Size matrix_num ( ++num_rots );
							x_new = rot_matrix[matrix_num -1]*x_new;
							y_new = rot_matrix[matrix_num -1]*y_new;
							origin_new = rot_matrix[matrix_num -1]*origin_new;
						} else if ( transform_type[i-1].second == "trans" ) {
							Size vector_num ( ++num_trans );
							origin_new = trans_vector[vector_num -1] + origin_new;
						}
						VirtualCoordinate coord( x_new, y_new, origin_new );
						Size jump_num ( virtual_coordinates_.size() + 1 );
						string tag;
						std::ostringstream stream;
						stream << std::setfill('0') << std::setw(4) << jump_num;
						tag = "CLONE" + stream.str();
						string identifier = "VRT" + stream.str();
						virtual_coordinates_.insert( make_pair( identifier, coord ) );
						jump_string_to_virtual_pair_.insert( make_pair( tag, make_pair( identifier, "SUBUNIT" ) ) );
						virt_id_to_virt_num_.insert( make_pair( identifier, jump_num ) );
						virt_num_to_virt_id_.insert( make_pair( jump_num, identifier ) );
						virt_id_to_subunit_num_.insert( make_pair( identifier, virt_id_to_subunit_num_.size()+1 ) );
						subunit_num_to_virt_id_.insert( make_pair( subunit_num_to_virt_id_.size()+1, identifier ) );
						virt_coord = coord;
					}
				}
			} else {
				utility_exit_with_message("[ERROR] Error reading symm def file while at '"+line+"'");
			}
		} else {
			// otherwise just reading std info
			if ( tokens[1] == "virtual_coordinates_start" ) {
				// Start virtual coordinate system block
				read_virtual_coords = true;
			} else if ( tokens[1] == "virtual_transforms_start" ) {
				// Start transformation block
				read_transforms = true;
				if ( tokens.size() >= 2 && tokens[2] == "consecutive" ) {
					read_consecutive_transforms = true;
				}
			} else if ( tokens[1] == "symmetry_name" ) {
				// parse the name of the symmetry type. Not required
				symmetry_name_ = tokens[2];
			} else if ( tokens[1] == "subunits" ) {
				// parse the number of subunits. Required.
				subunits_ = utility::string2int( tokens[2] );
				read_subunits = true;
				for ( Size i =1; i<=subunits_; i++ ) {
					include_subunit_.push_back(i);
					output_subunit_.push_back(i);
				}
			} else if ( tokens[1] == "number_of_interfaces" ) {
				// parse the number of different type of interfaces. Required
				interfaces_ = utility::string2int( tokens[2] );
			} else if ( tokens[1] == "anchor_residue" ) {
				// anchor residue id
				anchor_residue_ = tokens[2];
			} else if ( tokens[1] == "recenter" ) {
				recenter_ =  true;
			} else if ( tokens[1] == "E" && tokens[2] == "=" ) {
				// Define how the total energy should be calculated. The format is
				// E = X*E2 + Y*E3 + ...
				// where X and Y are integers. If not specified we assume that E = subunits_*E2.
				// Spaces between "E", "=" and the multipliers are necessary
				//for ( Size i = 2; i <= tokens.size(); i=i+2 ) {
				// score_multiply_subunit_string.push_back( tokens[i] );
				//}
				score_multiply_subunit_string.push_back( line );
			} else if ( tokens[1] == "connect_virtual" ) {
				// Read information on how virtual coordinate systems are connected by jumps. Not required
				connect_virtuals_specified = true;
				if ( tokens.size() < 4 ) {
					utility_exit_with_message( "[ERROR] You have to give jump identifier together with two virtual residue ids..." );
				}
				if ( jump_string_to_virtual_pair_.find( tokens[2] ) != jump_string_to_virtual_pair_.end() ) {
					utility_exit_with_message( "[ERROR] jump identifiers have to be unique..." + tokens[2] + " found twice " );
				}
				if ( virtual_coordinates_.find( tokens[3] ) == virtual_coordinates_.end() ) {
					utility_exit_with_message( "[ERROR] Virtual residue " + tokens[3] + " in jump " + tokens[2] + " not found" );
				}
				if ( tokens[4] != "SUBUNIT" && virtual_coordinates_.find( tokens[4] ) == virtual_coordinates_.end() ) {
					utility_exit_with_message( "[ERROR] Virtual residue " + tokens[4] + " in jump " + tokens[2] + " not found" );
				}

				jump_string_to_virtual_pair_[ tokens[2] ] = make_pair( tokens[3], tokens[4] );
				if ( tokens[4] == "SUBUNIT" ) {
					if ( tokens.size()==4 ) {
						if ( !read_subunits ) ++subunits_;
						Size newId = virt_id_to_subunit_num_.size()+1;
						virt_id_to_subunit_num_[ tokens[3] ] = newId;
						subunit_num_to_virt_id_[ newId ] = tokens[3];
					} else {
						// more than one component
						Size newId = virt_id_to_subunit_num_.size()+1;
						virt_id_to_subunit_num_[ tokens[3] ] = newId;
						subunit_num_to_virt_id_[ newId ] = tokens[3];
						if ( !have_virtual_coordinates ) utility_exit_with_message("multicomponent symmetry not supported with virtual_transforms_start/stop!");
						string subchain = tokens[5];
						if ( subchain.size() != 1 || subchain[0]==' ' ) utility_exit_with_message("[ERRIR] bad chain: "+subchain+" specified on SUBUNIT line:\n"+line);
						virt_id_to_subunit_chain_[ tokens[3] ] = subchain[0];
						subchains.insert(subchain[0]);
						if ( tokens.size() > 5 && tokens[5][0] != '#' ) {
							virt_id_to_subunit_residue_[ tokens[3] ] = tokens[6];
						} else {
							virt_id_to_subunit_residue_[ tokens[3] ] = get_anchor_residue(); // default subanchor
							continue;
						}
						if ( tokens.size() > 6 ) {
							if ( tokens[7][0] != '#' ) {
								utility_exit_with_message("[ERROR] extra input on SUBUNIT line:\n"+line);
							}
						}
					}
				}
			} else if ( tokens[1] == "set_dof" ) {
				// Store the degrees of freedom associated with a jump to a virtual residue. If not given
				// we assume that all dofs are allowed. Format:
				// set_dof x y z angle_x angle_x angle_y angle_z
				// If one of these tokens are missing it is assumed that the dof is disallowed
				// if we have no connect virtuals add connections by virtuals by default
				if ( !connect_virtuals_specified ) {
					for ( Size i = 2; i<= virtual_coordinates_.size(); ++i ) {
						std::ostringstream str, str2;
						str << std::setfill('0') << std::setw(4) << i;
						str2 << std::setfill('0') << std::setw(4) << i-1;
						string tag = "JUMP" + str2.str();
						string vrt_start = "VRT" + str2.str();
						string vrt_end = "VRT" + str.str();
						jump_string_to_virtual_pair_.insert( make_pair( tag, make_pair( vrt_start, vrt_end ) ) );
					}
				}
				// First convert jump tags to jump numbers
				// ASSUMES ALL JUMPS ALREADY DECLARED
				if ( set_jump_numbers_from_tags ) {
					Size insert_pos(0);
					map< string, pair< string, string > >::const_iterator itv_start = jump_string_to_virtual_pair_.begin();
					map< string, pair< string, string > >::const_iterator itv_end = jump_string_to_virtual_pair_.end();

					// jumps to the subunit first
					for ( map< string, pair< string, string > >::const_iterator itv = itv_start; itv != itv_end; ++itv ) {
						pair< string, string > connect( itv->second );
						string pos_id1( connect.first );
						string pos_id2( connect.second );
						if ( pos_id2 == "SUBUNIT" ) {
							++insert_pos;
							jump_string_to_jump_num_.insert( make_pair( itv->first, insert_pos ) );
						}
					}
					// then intra-VRT jumps
					for ( map< string, pair< string, string > >::const_iterator itv = itv_start; itv != itv_end; ++itv ) {
						pair< string, string > connect( itv->second );
						string pos_id2( connect.second );
						// We have already added the jumps from virtual residues to their corresponding subunits
						if ( pos_id2 == "SUBUNIT" ) continue;
						++insert_pos;
						jump_string_to_jump_num_.insert( make_pair( itv->first, insert_pos ) );
					}
					set_jump_numbers_from_tags = false;
				}

				if ( tokens.size() < 3 ) {
					utility_exit_with_message("[ERROR] Error reading set_dof line '"+line+"'");
				}

				string jump_id ( tokens[2] );
				if ( jump_string_to_jump_num_.find(jump_id) == jump_string_to_jump_num_.end() ) {
					string error( "[ERROR] Jump id is not found..." + jump_id );
					utility_exit_with_message( error );
				}
				Size jump_nbr( jump_string_to_jump_num_.find(jump_id)->second );
				SymDof dof;
				dof.add_dof_from_string(tokens);
				dofs_[ jump_nbr ] = dof;
			} else if ( tokens[1] == "include_subunit" ) {
				include_subunit_.clear();
				for ( Size i=2; i<=tokens.size(); i++ ) {
					core::Size j ( utility:: string2int( tokens[i] ) );
					include_subunit_.push_back(j);
				}
			} else if ( tokens[1] == "output_subunit" ) {
				output_subunit_.clear();
				for ( Size i =2; i<=tokens.size(); i++ ) {
					core::Size j ( utility:: string2int( tokens[i] ) );
					output_subunit_.push_back(j);
				}
			} else if ( tokens[1] == "set_jump_group" ) {
				// set groups of jumps that move together
				int master(-1);

				// First convert jump tags to jump numbers
				// ASSUMES ALL JUMPS ALREADY DECLARED
				if ( set_jump_numbers_from_tags ) {
					Size insert_pos(0);
					map< string, pair< string, string > >::const_iterator itv_start = jump_string_to_virtual_pair_.begin();
					map< string, pair< string, string > >::const_iterator itv_end = jump_string_to_virtual_pair_.end();

					// jumps to the subunit first
					for ( map< string, pair< string, string > >::const_iterator itv = itv_start; itv != itv_end; ++itv ) {
						pair< string, string > connect( itv->second );
						string pos_id1( connect.first );
						string pos_id2( connect.second );
						if ( pos_id2 == "SUBUNIT" ) {
							++insert_pos;
							jump_string_to_jump_num_.insert( make_pair( itv->first, insert_pos ) );
						}
					}
					// then intra-VRT jumps
					for ( map< string, pair< string, string > >::const_iterator itv = itv_start; itv != itv_end; ++itv ) {
						pair< string, string > connect( itv->second );
						string pos_id2( connect.second );
						// We have already added the jumps from virtual residues to their corresponding subunits
						if ( pos_id2 == "SUBUNIT" ) continue;
						++insert_pos;
						jump_string_to_jump_num_.insert( make_pair( itv->first, insert_pos ) );
					}
					set_jump_numbers_from_tags = false;
				}

				// The master is the jump that has a dof line
				vector1< Size > jump_nums;
				vector1< Real > jump_wts;
				for ( Size i = 3; i <= tokens.size(); ++i ) {
					// read tags of the form 'JUMP_0_0:0.333333'
					// which means this jump has a weight 0.333333
					vector1< string > jump_split_i ( utility::string_split( tokens[i], ':' ) );
					core::Real wt_i = 0.0;

					if ( jump_split_i.size() > 1 ) {
						wt_i = std::atof( jump_split_i[2].c_str() );
					}

					if ( jump_string_to_jump_num_.find( jump_split_i[1] ) == jump_string_to_jump_num_.end() ) {
						utility_exit_with_message("[ERROR] Undefined jump "+jump_split_i[1]+" in jump group '"+tokens[2]+"'");
					}

					Size jump_nbr( jump_string_to_jump_num_.find( jump_split_i[1] )->second );
					if ( dofs_.find(jump_nbr) != dofs_.end() ) {
						if ( master >= 0 ) {
							utility_exit_with_message("[ERROR] Multiple movable jumps specified in jump group '"+tokens[2]+"'");
						}
						master = jump_nbr;
					}

					jump_nums.push_back( jump_nbr );
					jump_wts.push_back( wt_i );
				}

				// If no dof is specified for any jump in this group just use the first group member
				if ( master < 0 ) {
					master = jump_nums[1];
				}
				runtime_assert (master >= 0);

				// Now set the jump group
				//vector1<Size> thisCloneList(0);
				vector1< pair<Size,Real> > thisCloneList(0);
				for ( Size i = 1; i <= jump_nums.size(); ++i ) {
					if ( (int)jump_nums[i] == master ) {
						if ( jump_wts[i] != 1.0 ) {
							TR.Warning << "Setting weight of master jump ( jump-id=" << master << " ) to 1.0 ";
							if ( jump_wts[i] == 0.0 ) {
								TR.Warning << "(was undefined)" << endl;
							} else {
								TR.Warning << "(was " << jump_wts[i] << ")" << endl;
							}
						}
					} else {
						thisCloneList.push_back( make_pair(jump_nums[i],jump_wts[i]) );
					}
				}
				jump_clones_[master] = thisCloneList;

				if ( jump_nums.size() > 1 ) {
					TR.Warning << "Setting jump_group " << tokens[2] << ": [master " << master << "] ";
					for ( Size i = 1; i <= jump_nums.size(); ++i ) {
						if ( (int)jump_nums[i] != master ) TR.Warning << " " << jump_nums[i] << ":" << jump_wts[i] << " ";
					}
					TR.Warning << endl;
				}
			} else if ( tokens[1] == "slide_type" ) {
				if ( tokens.size() != 2 ) utility_exit_with_message("[ERROR] Error reading slide_type '"+line+"'");
				if ( tokens[2] == "SEQUENTIAL" ) {
					slide_info_.set_slide_type( SEQUENTIAL );
				} else if ( tokens[2] == "ORDERED_SEQUENTIAL" ) {
					slide_info_.set_slide_type( ORDERED_SEQUENTIAL );
				} else if ( tokens[2] == "RANDOM" ) {
					slide_info_.set_slide_type( RANDOM );
				} else {
					utility_exit_with_message("[ERROR] Unknown slide_type '"+tokens[2]+"'");
				}
			} else if ( tokens[1] == "slide_criteria_type" ) {
				if ( tokens.size() != 2 ) utility_exit_with_message("[ERROR] Error reading slide_criteria_type '"+line+"'");
				if ( tokens[2] == "CEN_DOCK_SCORE" ) {
					slide_info_.set_SlideCriteriaType( CEN_DOCK_SCORE );
				} else if ( tokens[2] == "FA_REP_SCORE" ) {
					slide_info_.set_SlideCriteriaType( FA_REP_SCORE );
				} else if ( tokens[2] == "CONTACTS" ) {
					slide_info_.set_SlideCriteriaType( CONTACTS );
				} else {
					utility_exit_with_message("[ERROR] Unknown slide_type '"+tokens[2]+"'");
				}
			} else if ( tokens[1] == "slide_criteria_val" ) {
				if ( tokens.size() != 2 ) utility_exit_with_message("[ERROR] Error reading slide_criteria_val '"+line+"'");
				slide_info_.set_SlideCriteriaVal(tokens[2]);
			} else if ( tokens[1] == "slide_order" ) {
				for  ( Size record = 2; record <= tokens.size(); ++record ) {
					slide_order_string_.push_back( tokens[record] );
				}
			} else {
				utility_exit_with_message("[ERROR] Error reading symm def file while at '"+line+"'");
			}
		}
	}

	// postprocess multi-component related stuff // sheffler
	// count num subs per component, correct raw virt_id_to_subunit_num_
	if ( virt_id_to_subunit_chain_.size() > 0 ) {
		vector1<char> chains(subchains.begin(),subchains.end());
		if ( chains.size() <= 1 ) utility_exit_with_message("processing multicomponent symmetry, but only one chain!");
		std::sort(chains.begin(),chains.end());
		vector1<char>::const_iterator i = chains.begin();
		TR << "MULTICOMPONENT " << "You have specified the following chains for multi-component:";
		TR << " primary: " << *i;
		TR << ", secondary:";
		for ( ++i; i != chains.end(); ++i ) TR << " " << *i;
		TR << std::endl;
		char firstchain = chains[1];

		if ( virt_id_to_subunit_chain_.size() != virt_id_to_subunit_num_.size() ) utility_exit_with_message("missing component chains!");
		map<char,Size> chaincount;
		for ( map<string,char>::const_iterator i = virt_id_to_subunit_chain_.begin(); i != virt_id_to_subunit_chain_.end(); ++i ) {
			char const & chain = i->second;
			if ( chaincount.count(chain)==0 ) chaincount[chain] = 0;
			chaincount[chain]++;
		}
		if ( chaincount.size() == 1 ) utility_exit_with_message("For compatibility, don't use multicomponent format with only one component!!!");
		subunits_ = chaincount.begin()->second;
		num_components_ = chaincount.size();
		for ( map<char,Size>::const_iterator i = chaincount.begin(); i != chaincount.end(); ++i ) {
			if ( i->second != subunits_ ) {
				for ( map<char,Size>::const_iterator j = chaincount.begin(); j != chaincount.end(); ++j ) {
					TR << "SUBUNIT " << j->first << " num subs: " << j->second << endl;
				}
				utility_exit_with_message("[ERROR] bad number of subunits");
			}
		}

		// compute reference xforms
		map<pair<char,Size>,Xform> frames;
		for ( map<string,Size>::const_iterator i = virt_id_to_subunit_num_.begin(); i != virt_id_to_subunit_num_.end(); ++i ) {
			string const & virt_id = i->first;
			Size const & subnum = i->second;
			char const & chain = virt_id_to_subunit_chain_[virt_id];
			VirtualCoordinate const & vc( virtual_coordinates_[i->first] );

			Vec zaxis( vc.get_x().cross(vc.get_y()).normalized() );
			if ( vc.get_mirror_z() ) zaxis *= -1.0; // Mirror the z-axis if it says to do so in the VirtualCoordinate.
			Mat const R( Mat::cols( vc.get_x().normalized(), vc.get_y().normalized(), zaxis ) );
			Vec origin( vc.get_origin() );
			Xform toframe(R,origin);

			if ( vc.get_origin().distance(toframe  *Vec(0,0,0)             ) > 1e-10 ||
					vc.get_x().normalized().distance(toframe.R*Vec(1,0,0).normalized()) > 1e-10 ||
					vc.get_y().normalized().distance(toframe.R*Vec(0,1,0).normalized()) > 1e-10 ||
					zaxis.normalized().distance(toframe.R*Vec(0,0,1).normalized()) > 1e-10 ) {
				std::cerr << virt_id << " " << chain << " " << subnum << " " << "orig  " << vc.get_origin() << endl;
				std::cerr << virt_id << " " << chain << " " << subnum << " " << "xform " << toframe  *Vec(0,0,0) << endl;
				std::cerr << virt_id << " " << chain << " " << subnum << " " << "orig  " << vc.get_x() << endl;
				std::cerr << virt_id << " " << chain << " " << subnum << " " << "xform " << toframe.R*Vec(1,0,0) << endl;
				std::cerr << virt_id << " " << chain << " " << subnum << " " << "orig  " << vc.get_y() << endl;
				std::cerr << virt_id << " " << chain << " " << subnum << " " << "xform " << toframe.R*Vec(0,1,0) << endl;
				std::cerr << virt_id << " " << chain << " " << subnum << " " << "orig  " << zaxis << endl;
				std::cerr << virt_id << " " << chain << " " << subnum << " " << "xform " << toframe.R*Vec(0,0,1) << endl;
				utility_exit_with_message("computed frame xforms not correct!");
			}
			frames[make_pair(chain,(subnum-1)%subunits_+1)] = toframe;
		}

		// compute relative xforms
		map<pair<char,Size>,Xform> relxforms;
		for ( map<string,Size>::const_iterator i = virt_id_to_subunit_num_.begin(); i != virt_id_to_subunit_num_.end(); ++i ) {
			string const & virt_id = i->first;
			Size const & subnum = i->second;
			char const & chain = virt_id_to_subunit_chain_[virt_id];
			relxforms[make_pair(chain,subnum)] = frames[make_pair(chain,(subnum-1)%subunits_+1)] * ~frames[make_pair(chain,1)];
		}

		map<pair<char,Size>,Size> subperm;
		for ( map<string,Size>::const_iterator i = virt_id_to_subunit_num_.begin(); i != virt_id_to_subunit_num_.end(); ++i ) {
			string const & virt_id = i->first;
			Size const & subnum = i->second;
			char const & chain = virt_id_to_subunit_chain_[virt_id];
			if ( chain == firstchain ) continue;
			Xform const & xform1( relxforms[make_pair(chain,subnum)] );
			for ( map<pair<char,Size>,Xform>::const_iterator j = relxforms.begin(); j != relxforms.end(); ++j ) {
				Size const & subnum2(j->first.second);
				char const & chain2 (j->first.first);
				Xform const & xform2(j->second);
				if ( chain2!=firstchain ) continue;
				if ( xform1.distance(xform2) < 0.000001 ) {
					if ( subperm.count(make_pair(chain,subnum)) > 0 ) {
						utility_exit_with_message(string("multiple matching transforms for chain ")+
							chain+" in primary chain "+firstchain+"! probably a malformed or imprecise symmetry");
					}
					subperm[make_pair(chain,subnum)] = (subnum2-1)%subunits_+1;
				}
			}
			if ( subperm.find(make_pair(chain,subnum)) == subperm.end() ) {
				utility_exit_with_message(string("can't find matching transform for chain ")+
					chain+" in primary chain "+firstchain+"! probably a malformed or imprecise symmetry");
			}
		}

		// now map the extra component subunits
		for ( map<string,Size>::const_iterator i = virt_id_to_subunit_num_.begin(); i != virt_id_to_subunit_num_.end(); ++i ) {
			string const & virt_id = i->first;
			Size const & subnum = i->second;
			char const & chain = virt_id_to_subunit_chain_[virt_id];
			Size newsubnum = 0;
			if ( chain==firstchain ) {
				newsubnum = (subnum-1)%subunits_+1;
				TR << "MULTICOMPONENT " << "replace raw subunit "<< virt_id <<" num in primary chain " << chain << " with normalized subunit: "
					<< virt_id_to_subunit_num_[virt_id] << " -> " << newsubnum << ", " << virt_id << " -> " << subunit_num_to_virt_id_[newsubnum] << endl;
			} else {
				if ( subperm.find(make_pair(chain,subnum)) == subperm.end() ) {
					utility_exit_with_message("error computing correct subnum for component!");
				}
				newsubnum = subperm[make_pair(chain,subnum)];
				TR << "MULTICOMPONENT " << "replace raw subunit "<< virt_id <<" num in secondary chain " << chain << " with correctly permuted subunit: "
					<< virt_id_to_subunit_num_[virt_id] << " -> " << newsubnum << ", " << virt_id << " -> " << subunit_num_to_virt_id_[newsubnum] << endl;
			}
			if ( newsubnum==0 ) utility_exit_with_message("error computing correct subnum for component!");
			virt_id_to_subunit_num_[virt_id] = newsubnum;
		}

		// another sanity check
		map<char,string> chainres;
		for ( map<string,char>::const_iterator i = virt_id_to_subunit_chain_.begin(); i != virt_id_to_subunit_chain_.end(); ++i ) {
			string const & virt_id = i->first;
			char const & chain = i->second;
			// make sure [<res req>] are all same
			if ( chainres.count(chain)==0 ) chainres[chain] = virt_id_to_subunit_residue_[virt_id];
			if ( chainres[chain] != virt_id_to_subunit_residue_[virt_id] ) {
				std::cerr << "ERROR in chain " << chain << " " << virt_id << endl;
				utility_exit_with_message("SUBUNIT <chain> [<res>] must match!");
			}
		}

		// sanity check num subs pre component are all == num subunits (subunits_)
		components_ = chains;
		name2component_ = virt_id_to_subunit_chain_;
		map<Size,string> jnum2dofname_;
		for ( std::map<std::string,Size>::const_iterator i = jump_string_to_jump_num_.begin(); i != jump_string_to_jump_num_.end(); ++i ) {
			jnum2dofname_[i->second] = i->first;
		}

		for ( map< Size, SymDof >::iterator i = dofs_.begin(); i != dofs_.end(); ++i ) {
			string dofname = jnum2dofname_[i->first];
			jname2components_[dofname] = components_moved_by_jump(dofname);
			jname2subunits_[dofname] = subunits_moved_by_jump(dofname);
		}

	} else { // standard symmetry, one component
		num_components_ = 1;
	}
	// end sheffler mods

	// Done parsing the symm file, now do some post-processing
	// 1) find root of the system.
	//      All jumps are defined upstream to downstream.
	//      The only virtual residue that is not a downstream jump partner is the root
	int root(-1);
	map< string, pair< string, string > >::const_iterator it_start = jump_string_to_virtual_pair_.begin();
	map< string, pair< string, string > >::const_iterator it_end = jump_string_to_virtual_pair_.end();
	map< string, pair< string, string > >::const_iterator it, it2;

	Size nvrt( virt_num_to_virt_id_.size() );
	vector1<bool> downstream_targets( nvrt, false );
	for ( it = it_start; it != it_end; ++it ) {
		pair< string, string > connect( it->second );
		if ( connect.second == "SUBUNIT" ) continue;
		if ( downstream_targets[ virt_id_to_virt_num_[ connect.second ] ] ) {
			utility_exit_with_message( "[ERROR] Cycle found in connect_virtual" );
		}
		downstream_targets[ virt_id_to_virt_num_[ connect.second ] ] = true;
	}

	for ( int i=1; i<=(int)nvrt; ++i ) {
		if ( !downstream_targets[i] ) {
			if ( root >= 0 ) {
				utility_exit_with_message( "[ERROR] Foldtree not closed ("+virt_num_to_virt_id_[root]+","+ virt_num_to_virt_id_[i]+"...");
			}
			root = i;
		}
	}
	if ( root < 0 ) {
		utility_exit_with_message( "[ERROR] Cycle found in connect_virtual" );
	}
	root_ = root;

	vector1< Real > score_multiply_subunit_vector;

	// Initialize the first subunit with a factor corresponding to the number of
	// subunits. The rest of the values are set to 0
	for ( Size i = 1; i<= subunits_; ++i ) {
		score_multiply_subunit_vector.push_back(0.0);
	}
	// VRTs have a multiplier of 1
	for ( Size i = subunits_ + 1; i<= subunits_ + virtual_coordinates_.size(); ++i ) {
		score_multiply_subunit_vector.push_back(1.0);
	}

	if ( score_multiply_subunit_string.size() == 0 ) {
		utility_exit_with_message( "[ERROR] No total energy line specified!" );
	}
	// std::cout << "score_multiply_subunit_string " << line << std::endl;
	vector1< string> split_1 = utility::string_split( score_multiply_subunit_string[1], '=' );
	if ( split_1.size() < 2 ) {
		utility_exit_with_message( "[ERROR] Error parsing line '"+score_multiply_subunit_string[1]+"'" );
	}
	vector1< string> split_2 = utility::string_split( split_1[2], '+' );

	for ( Size i = 1; i<=split_2.size(); ++i ) {
		// trim whitespace
		utility::trim(split_2[i], " ");

		vector1< string> split_3 = utility::string_split( split_2[i], '*' );
		Real factor=1;
		string virtual_residues = split_3[1];
		if ( split_3.size() > 1 ) {
			factor = utility::string2Real( split_3[1] );
			virtual_residues = split_3[2];
		}
		utility::trim( virtual_residues, "()" );
		vector1< string> virtual_residues_split ( utility::string_split( virtual_residues, ':' ) );

		Size subunit(0);
		if ( virtual_residues_split.size() == 1 ) {
			utility::trim(virtual_residues_split[1], " ");
			if ( virt_id_to_subunit_num_.find( virtual_residues_split[1] ) == virt_id_to_subunit_num_.end() ) {
				utility_exit_with_message( "[ERROR] VRT " + virtual_residues_split[1] + " not attached to a subunit");
			}
			score_subunit_ = virt_id_to_subunit_num_.find( virtual_residues_split[1] )->second;
			subunit = score_subunit_;
		} else if ( virtual_residues_split.size() == 2 ) {
			utility::trim(virtual_residues_split[2], " ");
			if ( virt_id_to_subunit_num_.find( virtual_residues_split[2] ) == virt_id_to_subunit_num_.end() ) {
				utility_exit_with_message( "[ERROR] VRT " + virtual_residues_split[2] + " not attached to a subunit");
			}
			subunit = virt_id_to_subunit_num_.find( virtual_residues_split[2] )->second;
		} else {
			utility_exit_with_message( "[ERROR] Error parsing 'E =' line while at " + split_2[i] );
		}
		score_multiply_subunit_vector[subunit] = factor;
	}
	set_score_multiply_subunit( score_multiply_subunit_vector );

	// process slide order information
	vector< Size > slide_order;
	for ( vector< string >::iterator it = slide_order_string_.begin(); it != slide_order_string_.end(); ++it ) {
		string jump_id ( *it );
		if ( jump_string_to_jump_num_.find(jump_id) == jump_string_to_jump_num_.end() ) {
			string error( "[ERROR] Jump id is not found..." + jump_id );
			utility_exit_with_message( error );
		}
		Size jump_nbr( jump_string_to_jump_num_.find(jump_id)->second );
		slide_order.push_back(jump_nbr);
	}
	slide_info_.set_slide_order(slide_order);
	// Check that the information given makes sense and is enough to
	// specifiy the symmetry
	sanity_check();
	// Print data for symmetry
	show();
}

// @details function to check if input data makes sense. We need to add many more checks
// so that rosetta never exits without a error message at this stage
void
SymmData::sanity_check()
{
	if ( virtual_coordinates_.size() < 1 ) {
		utility_exit_with_message( "[ERROR] No virtual atoms specified..." );
	}
	if ( subunits_ < 1 ) {
		utility_exit_with_message( "[ERROR] Need to give number of subunits..." );
	}
	if ( interfaces_ < 1 ) {
		utility_exit_with_message( "[ERROR] Need to give number of interfaces..." );
	}
	if ( score_multiply_subunit_.size() < 1 ) {
		utility_exit_with_message( "[ERROR] Need to specify how to calculate symmetrical energy..." );
	}
	Size subunits(0);
	map< string, pair< string, string > >::const_iterator itv_start = jump_string_to_virtual_pair_.begin();
	map< string, pair< string, string > >::const_iterator itv_end = jump_string_to_virtual_pair_.end();
	for ( map< string, pair< string, string > >::const_iterator itv = itv_start; itv != itv_end; ++itv ) {
		pair< string, string > connect( itv->second );
		if ( connect.second == "SUBUNIT" ) ++subunits;
	}
	if ( subunits != num_components_*subunits_ ) {
		std::cerr << "subunits from connect_virtual: " << subunits << endl;
		std::cerr << "subunits declared: " << subunits_ << endl;
		std::cerr << "num symmetric components: " << num_components_ << endl;
		utility_exit_with_message( "[ERROR] The number of subunits is not equal to the number of jumps from virtual residues to subunits..." );
	}
	bool nullchain = false;
	for ( map<string,char>::const_iterator i = virt_id_to_subunit_chain_.begin(); i != virt_id_to_subunit_chain_.end(); ++i ) {
		if ( i->second == (char)0 ) nullchain = true;
		if ( nullchain && i->second != (char)0 ) utility_exit_with_message("[ERROR] all SUBUNITs must have chain if any have chain!");
	}
}

// @details print the symmetry data that we have initialized from file
void
SymmData::show()
{
	TR << "symmetry name: " << symmetry_name_ << endl;
	TR << "number of subunits: " << subunits_ << endl;
	TR << "number of interfaces: " << interfaces_ << endl;
	TR << "score subunit number: " << virt_num_to_virt_id_.find( score_subunit_ )->second << endl;
	TR << "anchor the subunits at residue: " << anchor_residue_ << endl;

	map< string, VirtualCoordinate >::iterator vit;
	map< string, VirtualCoordinate >::iterator vit_begin = virtual_coordinates_.begin();
	map<  string, VirtualCoordinate >::iterator vit_end = virtual_coordinates_.end();

	for ( vit = vit_begin; vit != vit_end; ++vit ) {
		string identifier( (*vit).first );
		VirtualCoordinate coord( (*vit).second );
		TR << " Virtual coordinate system " << identifier << endl;
		TR << "x: " << coord.get_x()(1) << " " << coord.get_x()(2) << " " << coord.get_x()(3)
			<< endl;
		TR << "y: " << coord.get_y()(1) << " " << coord.get_y()(2) << " " << coord.get_y()(3)
			<< endl;
		TR << "origin: " << coord.get_origin()(1) << " " << coord.get_origin()(2) << " " << coord.get_origin()(3)
			<< endl;
	}

	map< Size, SymDof >::iterator it;
	map< Size, SymDof >::iterator it_begin = dofs_.begin();
	map< Size, SymDof >::iterator it_end = dofs_.end();
	for ( it = it_begin; it != it_end; ++it ) {
		int jump_nbr ( (*it).first );
		SymDof dof ( (*it).second );
		TR << "Dof for jump: " << jump_nbr << endl;
		for ( Size i=1; i<=6; ++i ) {
			string dir ( "n2c" );
			if ( dof.jump_direction(i) == -1 ) dir = "c2n";
			if ( i == 1 ) TR << "x ";
			if ( i == 2 ) TR << "y ";
			if ( i == 3 ) TR << "z ";
			if ( i == 4 ) TR << "x_angle ";
			if ( i == 5 ) TR << "y_angle ";
			if ( i == 6 ) TR << "z_angle ";
			TR << dof.allow_dof(i) << ":" << dof.range1_lower(i) << "," << dof.range1_upper(i) << ":" << dof.range2_lower(i) << "," << dof.range2_upper(i)
				<< " " << dir << endl;
		}
	}

	map< string, pair< string, string > >::const_iterator itv_start = jump_string_to_virtual_pair_.begin();
	map< string, pair< string, string > >::const_iterator itv_end = jump_string_to_virtual_pair_.end();
	for ( map< string, pair< string, string > >::const_iterator itv = itv_start; itv != itv_end; ++itv ) {
		pair< string, string > connect( itv->second );
		string pos_id1( connect.first );
		string pos_id2( connect.second );
		TR << "Jump " << itv->first << " " << pos_id1 << " " << pos_id2 << endl;
	}
	TR << "Include subunit:";
	for ( vector1< Size >::iterator it = include_subunit_.begin(); it != include_subunit_.end(); ++it ) {
		TR << ' ' << (*it) ;
	}
	TR << endl;
	TR << "Output subunit:";
	for ( vector1< Size >::iterator it = output_subunit_.begin(); it != output_subunit_.end(); ++it ) {
		TR << ' ' << (*it) ;
	}
	TR << endl;
	TR << "SlideType: ";
	if ( slide_info_.get_slide_type() == SEQUENTIAL ) TR << "SEQUENTIAL" << endl;
	if ( slide_info_.get_slide_type() == ORDERED_SEQUENTIAL ) TR << "ORDERED_SEQUENTIAL" << endl;
	if ( slide_info_.get_slide_type() == RANDOM ) TR << "RANDOM" << endl;
	TR << "SlideCriteriaType: ";
	if ( slide_info_.get_SlideCriteriaType() == CEN_DOCK_SCORE ) TR << "CEN_DOCK_SCORE" << endl;
	if ( slide_info_.get_SlideCriteriaType() == FA_REP_SCORE ) TR << "FA_REP_SCORE" << endl;
	if ( slide_info_.get_SlideCriteriaType() == CONTACTS ) TR << "CONTACTS" << endl;
	TR << "SlideCriteriaVal: ";
	TR << slide_info_.get_SlideCriteriaVal() << endl;
	TR << "SlideOrder: ";
	if ( slide_order_string_.size() == 0 ) TR << "none";
	for ( vector< string >::iterator it = slide_order_string_.begin();
			it != slide_order_string_.end(); ++it ) {
		TR << ' ' << (*it) ;
	}

	TR << endl;
}


utility::vector1<char> const &
SymmData::get_components() const {
	return components_;
}

std::map<std::string,char> const &
SymmData::get_subunit_name_to_component() const {
	return name2component_;
}
std::map<std::string,utility::vector1<char> > const &
SymmData::get_jump_name_to_components() const {
	return jname2components_;
}

std::map<std::string,utility::vector1<Size> > const &
SymmData::get_jump_name_to_subunits() const {
	return jname2subunits_;
}

std::string const &
SymmData::get_parent_jump(std::string const & jname) const {
	if ( jump_string_to_virtual_pair_.find(jname) == jump_string_to_virtual_pair_.end() ) {
		utility_exit_with_message("can't find parent jump of jump"+jname);
	}
	std::string upvirt = jump_string_to_virtual_pair_.find(jname)->second.first;
	for ( map<string,pair<string,string> >::const_iterator i = jump_string_to_virtual_pair_.begin(); i != jump_string_to_virtual_pair_.end(); ++i ) {
		if ( i->second.second == upvirt ) return i->first;
	}
	// utility_exit_with_message("can't find parent jump of jump: "+jname);
	return NOPARENT;
}

std::string const &
SymmData::get_parent_virtual(std::string const & vname) const {
	for ( map<string,pair<string,string> >::const_iterator i = jump_string_to_virtual_pair_.begin(); i != jump_string_to_virtual_pair_.end(); ++i ) {
		if ( i->second.second == vname ) return i->second.first;
	}
	// utility_exit_with_message("can't find parent jump of virtual: "+vname);
	return NOPARENT;
}

bool
SymmData::is_ancestor_virtual(
	std::string const & ancestor,
	std::string const & child
) const {
	std::string const * p(&child);
	while ( *p != NOPARENT ) {
		if ( *p == ancestor ) return true;
		p = &get_parent_virtual(*p);
	}
	return false;
}

std::string const &
SymmData::get_least_common_ancertor_jump(std::string const & jname1, std::string const & jname2) const {
	std::set<string> seenit;
	std::string const * p(&jname1);
	while ( *p != NOPARENT ) {
		seenit.insert(*p);
		p = &get_parent_jump(*p);
	}
	p = &jname2;
	while ( *p != NOPARENT ) {
		if ( seenit.find(*p) != seenit.end() ) return *p;
		p = &get_parent_jump(*p);
	}
	return NOPARENT;
}

std::string const &
SymmData::get_least_common_ancertor_virtual(std::string const & vname1, std::string const & vname2) const {
	std::set<string> seenit;
	std::string const * p(&vname1);
	while ( *p != NOPARENT ) {
		seenit.insert(*p);
		p = &get_parent_virtual(*p);
	}
	p = &vname2;
	while ( *p != NOPARENT ) {
		if ( seenit.find(*p) != seenit.end() ) return *p;
		p = &get_parent_virtual(*p);
	}
	return NOPARENT;
}

vector1<string>
SymmData::leaves_of_jump(std::string const & jname) const {
	if ( jump_string_to_virtual_pair_.find(jname) == jump_string_to_virtual_pair_.end() ) {
		utility_exit_with_message("unknown jump name in leaves of jump");
	}
	string const & ancestor( jump_string_to_virtual_pair_.find(jname)->second.second );
	vector1<string> leaves;
	for ( map<std::string,Size>::const_iterator i = virt_id_to_subunit_num_.begin(); i != virt_id_to_subunit_num_.end(); ++i ) {
		if ( is_ancestor_virtual(ancestor,i->first) ) leaves.push_back(i->first);
	}
	return leaves;
}


vector1<char>
SymmData::components_moved_by_jump(std::string const & jname) const {
	vector1<string> leaves = leaves_of_jump(jname);
	std::set<char> components;
	for ( vector1<string>::const_iterator i = leaves.begin(); i != leaves.end(); ++i ) {
		if ( virt_id_to_subunit_chain_.find(*i) == virt_id_to_subunit_chain_.end() ) break;
		components.insert( virt_id_to_subunit_chain_.find(*i)->second );
	}
	return vector1<char>(components.begin(),components.end());
}

vector1<Size>
SymmData::subunits_moved_by_jump(std::string const & jname) const {
	vector1<string> leaves = leaves_of_jump(jname);
	std::set<Size> subunits;
	for ( vector1<string>::const_iterator i = leaves.begin(); i != leaves.end(); ++i ) {
		if ( virt_id_to_subunit_num_.find(*i) == virt_id_to_subunit_num_.end() ) break;
		subunits.insert( virt_id_to_subunit_num_.find(*i)->second );
	}
	return vector1<char>(subunits.begin(),subunits.end());
}

bool
operator==(
	SymmData const & a,
	SymmData const & b
) {
	if (
			(a.symmetry_name_ != b.symmetry_name_) ||
			(a.symmetry_type_ != b.symmetry_type_) ||
			(a.subunits_ != b.subunits_) ||
			(a.num_components_ != b.num_components_) ||
			(a.interfaces_ != b.interfaces_) ||
			(a.score_subunit_ != b.score_subunit_) ||
			(a.anchor_residue_ != b.anchor_residue_) ||
			(a.recenter_ != b.recenter_) ||
			(a.root_ != b.root_) ||
			(a.slide_info_ != a.slide_info_) ||
			!std::equal(
			a.slide_order_string_.begin(),
			a.slide_order_string_.end(),
			b.slide_order_string_.begin()) ||
			!std::equal(
			a.rotation_matrices_.begin(),
			a.rotation_matrices_.end(),
			b.rotation_matrices_.begin()) ||
			!std::equal(
			a.translation_matrices_.begin(),
			a.translation_matrices_.end(),
			b.translation_matrices_.begin()) ||
			!std::equal(
			a.virtual_coordinates_.begin(),
			a.virtual_coordinates_.end(),
			b.virtual_coordinates_.begin()) ||
			!std::equal(
			a.jump_string_to_virtual_pair_.begin(),
			a.jump_string_to_virtual_pair_.end(),
			b.jump_string_to_virtual_pair_.begin()) ||
			!std::equal(
			a.jump_string_to_jump_num_.begin(),
			a.jump_string_to_jump_num_.end(),
			b.jump_string_to_jump_num_.begin()) ||
			!std::equal(
			a.virt_id_to_virt_num_.begin(),
			a.virt_id_to_virt_num_.end(),
			b.virt_id_to_virt_num_.begin()) ||
			!std::equal(
			a.virt_id_to_subunit_num_.begin(),
			a.virt_id_to_subunit_num_.end(),
			b.virt_id_to_subunit_num_.begin()) ||
			!std::equal(
			a.virt_id_to_subunit_chain_.begin(),
			a.virt_id_to_subunit_chain_.end(),
			b.virt_id_to_subunit_chain_.begin()) ||
			!std::equal(
			a.virt_id_to_subunit_residue_.begin(),
			a.virt_id_to_subunit_residue_.end(),
			b.virt_id_to_subunit_residue_.begin()) ||
			!std::equal(
			a.virt_num_to_virt_id_.begin(),
			a.virt_num_to_virt_id_.end(),
			b.virt_num_to_virt_id_.begin()) ||
			!std::equal(
			a.subunit_num_to_virt_id_.begin(),
			a.subunit_num_to_virt_id_.end(),
			b.subunit_num_to_virt_id_.begin()) ||
			!std::equal(
			a.dofs_.begin(),
			a.dofs_.end(),
			b.dofs_.begin()) ||
			!std::equal(
			a.allow_virtual_.begin(),
			a.allow_virtual_.end(),
			b.allow_virtual_.begin()) ||
			!std::equal(
			a.score_multiply_subunit_.begin(),
			a.score_multiply_subunit_.end(),
			b.score_multiply_subunit_.begin()) ||
			!std::equal(
			a.include_subunit_.begin(),
			a.include_subunit_.end(),
			b.include_subunit_.begin()) ||
			!std::equal(
			a.output_subunit_.begin(),
			a.output_subunit_.end(),
			b.output_subunit_.begin()) ||
			!std::equal(
			a.components_.begin(),
			a.components_.end(),
			b.components_.begin()) ||
			!std::equal(
			a.name2component_.begin(),
			a.name2component_.end(),
			b.name2component_.begin()) ) {
		return false;
	}

	for (
			std::vector< std::vector< std::string > >::const_iterator
			ai = a.symm_transforms_.begin(), aie = a.symm_transforms_.end(),
			bi = b.symm_transforms_.begin();
			ai != aie; ++ai, ++bi ) {
		if ( bi == b.symm_transforms_.end() ||
				!std::equal(ai->begin(), ai->end(), bi->begin()) ) {
			return false;
		}
	}

	for (
			std::map< Size, WtedClones >::const_iterator
			ai = a.jump_clones_.begin(), aie = a.jump_clones_.end(),
			bi = b.jump_clones_.begin();
			ai != aie; ++ai, ++bi ) {
		if ( bi == b.jump_clones_.end() ||
				!std::equal(ai->second.begin(), ai->second.end(), bi->second.begin()) ) {
			return false;
		}
	}

	for (
			std::map< std::string, utility::vector1<char> >::const_iterator
			ai = a.jname2components_.begin(), aie = a.jname2components_.end(),
			bi = b.jname2components_.begin();
			ai != aie; ++ai, ++bi ) {
		if ( bi == b.jname2components_.end() ||
				!std::equal(ai->second.begin(), ai->second.end(), bi->second.begin()) ) {
			return false;
		}
	}

	for (
			std::map< std::string, utility::vector1<Size> >::const_iterator
			ai = a.jname2subunits_.begin(), aie = a.jname2subunits_.end(),
			bi = b.jname2subunits_.begin();
			ai != aie; ++ai, ++bi ) {
		if ( bi == b.jname2subunits_.end() ||
				!std::equal(ai->second.begin(), ai->second.end(), bi->second.begin()) ) {
			return false;
		}
	}

	return true;
}

} // symmetry
} // conformation
} // core
