// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/DownstreamAlgorithm.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/WriteUpstreamCoordinateKineamge.hh>

// Package headers
#include <protocols/match/downstream/DownstreamBuilder.hh>
#include <protocols/match/downstream/ClassicMatchAlgorithm.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/string_util.hh>

// C++ headers
#include <list>

#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {

/// TWO FUNCTIONS STOLEN FROM IAN: and slightly modified.
void print_node(
	std::ostream & out,
	int residue_num,
	int atom_num,
	core::Vector const & atom_xyz,
	core::chemical::ResidueType const & res,
	std::string const & extras = "" //< P for points, color, width/radius, etc.
)
{
	// atom_num is often 0 in fold tree, means no specific atom.
	// might as well use the first() one:
	if ( atom_num == 0 ) atom_num = 1;
	core::chemical::AtomType const & atom_type = res.atom_type(atom_num);
	// This info appears when you click on the point
	out << "{" << res.name3() << " " << residue_num
		<< " " << res.atom_name(atom_num) << " (" << atom_type.name() << ")"
		<< "}";
	// Color, width, etc. followed by coordinates
	//out << "col" << residue_num << " ";
	out << extras;
	out << " " << atom_xyz.x() << " " << atom_xyz.y() << " " << atom_xyz.z() << "\n";
}

void print_node(
	std::ostream & out,
	int residue_num,
	int atom_num,
	core::conformation::Residue const & res,
	std::string const & extras = "" //< P for points, color, width/radius, etc.
)
{
	print_node( out, residue_num, atom_num, res.xyz( atom_num ), res.type(), extras );
}


void print_node(
	std::ostream & out,
	int residue_num,
	std::string atom_name,
	core::conformation::Residue const & res,
	std::string const & extras = "" //< P for points, color, width/radius, etc.
)
{
	// atom_num is often 0 in fold tree, means no specific atom.
	// might as well use the first one:

	int atom_num;
	if ( atom_name == "" ) {
		atom_num = 1;
	} else {
		atom_num = res.atom_index( atom_name );
	}
	print_node( out, residue_num, atom_num, res, extras );
}

ResidueKinemageWriter::ResidueKinemageWriter() :
	master_( "" ),
	dominant_( true ),
	animate_( true ),
	group_( true ),
	write_virtual_atoms_( false )
{}

/// @details Don't write out multiple copies of the upstream coordinates if you're
/// drawing multiple instances of the downstream partner from a single rotamer.  Instead,
/// use the kinemage "instance" flag to point to the coordinates already written out
/// in this file.  This creates a smaller output file.
void
ResidueKinemageWriter::write_rsd_coords(
	std::ostream & ostr,
	Size const scaffold_build_point_id,
	Size const upstream_conf_id,
	core::conformation::Residue const & rsd,
	bool is_instance /* = false */
) const
{

	// intra-residue connections
	// do residues in different (~random) colors to help distinguish them
	//int const num_colors = 6;
	//std::string colors[num_colors] = {"pinktint", "peachtint", "yellowtint", "greentint", "bluetint", "lilactint"};
	//std::string color = colors[ rsd.seqpos() % num_colors ];
	std::string color = "bluetint";

	std::string tag = "";

	ostr << "@" << ( group_ ? "group" : "subgroup" ) << " { rot" << upstream_conf_id << " }";
	if ( animate_ ) ostr << " animate";
	if ( dominant_ ) ostr << " dominant";
	ostr << "\n";
	ostr << "@vectorlist {";
	if ( is_instance ) {
		ostr << rsd.name() << " " << scaffold_build_point_id << " " << upstream_conf_id;
	} else {
		ostr << rsd.name() << " " << scaffold_build_point_id;
	}
	ostr << "} color= " << color;
	if ( master_ != "" ) ostr << " master= {" << master_ << "}";
	if ( is_instance ) {
		ostr << " instance= {" << rsd.name() << " " << scaffold_build_point_id << "}";
	}
	ostr << "\n";

	if ( ! is_instance ) {
		for ( core::Size atom_i = 1; atom_i <= rsd.natoms(); ++atom_i ) {
			core::conformation::Residue::AtomIndices const & nbrs = rsd.nbrs(atom_i);
			for ( unsigned long atom_j : nbrs ) {
				if ( atom_j <= atom_i ) continue; // so we draw each bond just once, not twice
				bool const is_H = rsd.atom_is_hydrogen(atom_j) || rsd.atom_is_hydrogen(atom_i);

				if ( ! write_virtual_atoms_ && ( rsd.atom_type( atom_i ).element() == "X" || rsd.atom_type( atom_j ).element() == "X" ) ) continue;

				/// no backbone hydrogens...
				if ( rsd.atom_is_backbone( atom_i ) && rsd.atom_is_hydrogen( atom_i ) ) continue;
				if ( rsd.atom_is_backbone( atom_j ) && rsd.atom_is_hydrogen( atom_j ) ) continue;
				std::string const ptmaster = ( is_H ? " 'h'" : "" );
				print_node( ostr, rsd.seqpos(), atom_i, rsd, tag+"P"+ptmaster);
				print_node( ostr, rsd.seqpos(), atom_j, rsd, tag+ptmaster);
			}
		}
	}
}


void ResidueKinemageWriter::dominant( bool setting ) { dominant_ = setting; }
void ResidueKinemageWriter::animate(  bool setting ) { animate_  = setting; }
void ResidueKinemageWriter::group(    bool setting ) { group_    = setting; }
void ResidueKinemageWriter::master( std::string const & setting ) { master_ = setting; }
void ResidueKinemageWriter::write_virtual_atoms( bool setting ) { write_virtual_atoms_ = setting; }

WriteUpstreamCoordinateKinemage::WriteUpstreamCoordinateKinemage(
) :
	DownstreamAlgorithm( 1 ), // our geom_cst_id is irrelevant
	kinemage_file_name_( "NO_FILE" ),
	file_out_(),
	ostr_( file_out_ ),
	last_scaffold_build_point_( 0 ),
	nkins_( 0 ),
	n_downstream_to_output_( 0 ),
	return_pseudo_hits_( false )
{
}


WriteUpstreamCoordinateKinemage::WriteUpstreamCoordinateKinemage(
	std::string const & fname
) :
	DownstreamAlgorithm( 1 ), // our geom_cst_id is irrelevant
	kinemage_file_name_( fname ),
	file_out_(),
	ostr_( file_out_ ),
	last_scaffold_build_point_( 0 ),
	nkins_( 0 ),
	n_downstream_to_output_( 0 ),
	return_pseudo_hits_( false )
{
	file_out_.open( kinemage_file_name_.c_str()  );
}


WriteUpstreamCoordinateKinemage::WriteUpstreamCoordinateKinemage(
	std::ostream & ostr
) :
	DownstreamAlgorithm( 1 ), // our geom_cst_id is irrelevant
	kinemage_file_name_( "OSTREAM_MODE" ),
	file_out_(),
	ostr_( ostr ),
	last_scaffold_build_point_( 0 ),
	nkins_( 0 ),
	n_downstream_to_output_( 0 ),
	return_pseudo_hits_( false )
{
}

WriteUpstreamCoordinateKinemage::~WriteUpstreamCoordinateKinemage()
{
	if ( file_out_.is_open() ) file_out_.close();
}

downstream::DownstreamAlgorithmOP
WriteUpstreamCoordinateKinemage::clone() const
{
	return downstream::DownstreamAlgorithmOP( new WriteUpstreamCoordinateKinemage );
}

/// @brief To be invoked by derived classes, this function completes the building of
/// the downstream conformation once the coordinates of the upstream conformation
/// are known (and deemed non-colliding or, generally, pass any filter the upstream
/// builder would use).
std::list< Hit >
WriteUpstreamCoordinateKinemage::build(
	Size const scaffold_build_point_id,
	Size const upstream_conf_id,
	core::conformation::Residue const & rsd
) const
{
	std::list< Hit > empty_list;
	ResidueKinemageWriter writer;
	writer.master( "rotamers" );

	if ( last_scaffold_build_point_ != scaffold_build_point_id ) {
		ostr_ << "@kinemage {" << ++nkins_ << "}\n";
		ostr_ << "@title { matcher }\n";
		core::Vector const & ctr( rsd.has( "CB" ) ? rsd.xyz( "CB" ) : rsd.xyz("CA") );
		ostr_ << "@1center " << ctr.x() << " " << ctr.y() << " " << ctr.z() << "\n";
		ostr_ << "@1span 25\n";
		n_output_so_far_ = 0;
	}

	last_scaffold_build_point_ = scaffold_build_point_id;

	if ( match_algorithm_ ) {
		std::list< Hit > downstream_hits = match_algorithm_->build( scaffold_build_point_id, upstream_conf_id, rsd );
		if ( n_downstream_to_output_ != 0 && n_output_so_far_ >= n_downstream_to_output_ ) return downstream_hits;

		if ( dswriter_ ) {

			/// Write out the upstream residue coord once for each downstream hit so when we're
			/// animating we'll see the various downstream hits individually with their upstream residue.
			/// Take advantage of the kinemage "instance" flag so that the coordinates for the
			/// upstream residue are only written to the file once.
			for ( std::list< Hit >::const_iterator
					iter = downstream_hits.begin(), iter_begin = downstream_hits.begin(),
					iter_end = downstream_hits.end();
					iter != iter_end; ++iter ) {
				writer.write_rsd_coords( ostr_, scaffold_build_point_id, upstream_conf_id, rsd, (iter != iter_begin) );
				dswriter_->write_downstream_coordinates( *iter, ostr_ );

				++n_output_so_far_;
				if ( n_downstream_to_output_ != 0 && n_output_so_far_ >= n_downstream_to_output_ ) return downstream_hits;

			}
		} else {
			if ( n_downstream_to_output_ != 0 && n_output_so_far_ >= n_downstream_to_output_ ) return downstream_hits;

			/// write out the conformation if it generates at least one hit, but exclude it otherwise.
			if ( ! downstream_hits.empty() ) {
				writer.write_rsd_coords(  ostr_, scaffold_build_point_id, upstream_conf_id, rsd );
				++n_output_so_far_;
				if ( n_downstream_to_output_ != 0 && n_output_so_far_ >= n_downstream_to_output_ ) return downstream_hits;
			}
		}
		return downstream_hits; // return the set of found hits -- maybe someone listening wants them.
	} else {
		if ( n_downstream_to_output_ != 0 && n_output_so_far_ >= n_downstream_to_output_ ) return empty_list;

		writer.write_rsd_coords( ostr_, scaffold_build_point_id, upstream_conf_id, rsd );
		++n_output_so_far_;
	}

	if (  return_pseudo_hits_ ) {
		Hit hit;
		hit.first()[ 1 ] = scaffold_build_point_id;
		hit.first()[ 2 ] = upstream_conf_id;
		empty_list.push_back( hit );
	}

	return empty_list;
}

/// @details This is an acceptalbe return value for this debugging-purposed class.
bool
WriteUpstreamCoordinateKinemage::upstream_only() const
{
	return false;
}

bool
WriteUpstreamCoordinateKinemage::generates_primary_hits() const
{
	return true;
}


/// @details If this function is causing an exit, then there is a bug within the Matcher's
/// match-enumeration logic.  There is no meaningful way forward after this function is invoked.
/// It should not be invoked.  Truely, this class should not be used in match enumeration.
HitPtrListCOP
WriteUpstreamCoordinateKinemage::hits_to_include_with_partial_match( match_dspos1 const & ) const
{
	HitPtrListCOP empty;
	utility_exit_with_message( "Cannot invoke WriteUpstreamCoordinateKinemage::hits_to_include_with_partial_match()" );
	return empty;
}

WriteUpstreamCoordinateKinemage::Size
WriteUpstreamCoordinateKinemage::n_possible_hits_per_upstream_conformation() const
{
	if ( match_algorithm_ ) {
		return match_algorithm_->n_possible_hits_per_upstream_conformation();
	} else {
		return 0;
	}
}


void
WriteUpstreamCoordinateKinemage::set_kinemage_file_name( std::string const & filename )
{
	if ( kinemage_file_name_ == "OSTREAM_MODE" ) {
		utility_exit_with_message( "ERROR: Cannot set kinemage file name when in ostream mode" );
	}
	kinemage_file_name_ = filename;
	if ( file_out_.is_open() ) file_out_.close();
	file_out_.open( kinemage_file_name_.c_str()  );

}

void
WriteUpstreamCoordinateKinemage::set_match_algorithm( downstream::ClassicMatchAlgorithmCOP algorithm )
{
	match_algorithm_ = algorithm;
}

void
WriteUpstreamCoordinateKinemage::set_downstream_writer( DownstreamCoordinateKinemageWriterCOP dswriter )
{
	dswriter_ = dswriter;
}

void
WriteUpstreamCoordinateKinemage::set_n_downstream_to_output( Size n_downstream_to_output )
{
	n_downstream_to_output_ = n_downstream_to_output;
}

WriteUpstreamHitKinemage::WriteUpstreamHitKinemage()
:
	matches_output_count_( 0 ),
	use_default_master_( true ),
	animate_( false ),
	dominant_( false ),
	group_( false ),
	write_virtual_atoms_( false ),
	geom_id_( 0 ),
	kinemage_file_name_( "NO_FILE" ),
	file_out_(),
	ostr_( file_out_ )
{}

WriteUpstreamHitKinemage::WriteUpstreamHitKinemage(
	std::string const & fname
) :
	matches_output_count_( 0 ),
	use_default_master_( true ),
	animate_( false ),
	dominant_( false ),
	group_( false ),
	write_virtual_atoms_( false ),
	geom_id_( 0 ),
	kinemage_file_name_( fname ),
	file_out_(),
	ostr_( file_out_ )
{
	file_out_.open( kinemage_file_name_.c_str()  );
}


WriteUpstreamHitKinemage::WriteUpstreamHitKinemage( std::ostream & ostr ) :
	matches_output_count_( 0 ),
	use_default_master_( true ),
	animate_( false ),
	dominant_( false ),
	group_( false ),
	write_virtual_atoms_( false ),
	geom_id_( 0 ),
	kinemage_file_name_( "OSTREAM_MODE" ),
	file_out_(),
	ostr_( ostr )
{}

WriteUpstreamHitKinemage::~WriteUpstreamHitKinemage()
{
	if ( file_out_.is_open() ) file_out_.close();
}

void
WriteUpstreamHitKinemage::process_hit(
	Hit const & hit,
	core::conformation::Residue const & upstream_conformation
)
{
	output_hit( hit, upstream_conformation );
}

void
WriteUpstreamHitKinemage::output_hit(
	Hit const & hit,
	core::conformation::Residue const & upstream_conformation
)
{
	output_upstream_coordinates( hit, upstream_conformation );
	if ( dswriter_ ) {
		dswriter_->write_downstream_coordinates( hit, ostr_ );
	}
}

void
WriteUpstreamHitKinemage::output_upstream_coordinates(
	upstream_hit const & hit,
	core::conformation::Residue const & upstream_conformation
)
{
	ResidueKinemageWriter writer;
	writer.dominant( dominant_ ); writer.animate( animate_ ); writer.group( group_ );

	/// set the master
	if ( use_default_master_ ) {
		writer.master( "geom" + utility::to_string( geom_id_ ) );
	} else if ( master_ != "" ) {
		writer.master( master_ );
	}

	writer.write_rsd_coords( ostr_, hit.scaffold_build_id(), hit.upstream_conf_id(), upstream_conformation, false );

}


void
WriteUpstreamHitKinemage::start_new_match()
{
	if ( matches_output_count_ == 0 ) {
		ostr_ << "@kinemage 1\n";
	}
	++matches_output_count_;
	ostr_ << "@group { match " << matches_output_count_ << " } dominant animate\n";
}


void
WriteUpstreamHitKinemage::set_kinemage_file( std::string const & fname )
{
	if ( kinemage_file_name_ == "OSTREAM_MODE" ) {
		utility_exit_with_message( "ERROR: Cannot set kinemage file name when in ostream mode" );
	}
	kinemage_file_name_ = fname;
	if ( file_out_.is_open() ) file_out_.close();
	file_out_.open( kinemage_file_name_.c_str()  );
}

void
WriteUpstreamHitKinemage::set_dswriter( DownstreamCoordinateKinemageWriterOP dswriter )
{
	dswriter_ = dswriter;
}

void WriteUpstreamHitKinemage::geom_id( Size setting )
{
	geom_id_ = setting;
	if ( dswriter_ ) {
		std::string master( "geom" + utility::to_string( geom_id_ ));
		dswriter_->set_downstream_master( master );
	}
}

void
WriteUpstreamHitKinemage::set_master( std::string const & master )
{
	master_ = master;
	use_default_master_ = false;
}

void
WriteUpstreamHitKinemage::default_master( bool setting )
{
	use_default_master_ = setting;
}

/// @brief Returns whether or not the default master is being used.
bool
WriteUpstreamHitKinemage::default_master() const
{
	return use_default_master_;
}

void WriteUpstreamHitKinemage::animate( bool setting )
{
	animate_ = setting;
}

void WriteUpstreamHitKinemage::dominant( bool setting )
{
	dominant_ = setting;
}

void WriteUpstreamHitKinemage::group(    bool setting )
{
	group_ = setting;
}

void WriteUpstreamHitKinemage::write_virtual_atoms( bool setting )
{
	write_virtual_atoms_ = setting;
}


DownstreamCoordinateKinemageWriter::DownstreamCoordinateKinemageWriter() = default;
DownstreamCoordinateKinemageWriter::~DownstreamCoordinateKinemageWriter() = default;


SingleDownstreamResidueWriter::SingleDownstreamResidueWriter() = default;

SingleDownstreamResidueWriter::~SingleDownstreamResidueWriter() = default;


void
SingleDownstreamResidueWriter::write_downstream_coordinates(
	Hit const & hit,
	std::ostream & ostr
) const
{
	//std::cout << "outputting hit: ";
	//for ( Size ii = 1; ii <= 6; ++ii ) {
	// std::cout << hit.second[ ii ] << " ";
	//}
	//std::cout << std::endl;

	if ( ! restype_ ) {
		utility_exit_with_message( "ERROR: SingleDownstreamResidueWriter must have its residue type initialized!" );
	}

	utility::vector1< core::Vector > coords( restype_->natoms() );
	dsbuilder_->coordinates_from_hit( hit, all_atom_inds_, coords );

	// intra-residue connections
	// do residues in different (~random) colors to help distinguish them
	//int const num_colors = 6;
	//std::string colors[num_colors] = {"pinktint", "peachtint", "yellowtint", "greentint", "bluetint", "lilactint"};
	//std::string color = colors[ rsd.seqpos() % num_colors ];
	std::string color = "greentint";

	std::string tag = "";

	ostr << "@vectorlist {" << restype_->name() << " " << hit.first()[ 1 ] << " " << hit.first()[ 2 ] << " " << hit.first()[ 3 ];
	ostr << "} color= " << color << " master= {rotamers}";
	if ( master_ != "" ) {
		ostr << " master= {" << master_ << "}";
	}
	ostr << "\n";

	for ( core::Size atom_i = 1; atom_i <= restype_->natoms(); ++atom_i ) {
		core::chemical::AtomIndices const & nbrs = restype_->nbrs(atom_i);
		for ( unsigned long atom_j : nbrs ) {
			if ( atom_j <= atom_i ) continue; // so we draw each bond just once, not twice
			bool const is_H = restype_->atom_is_hydrogen(atom_j) || restype_->atom_is_hydrogen(atom_i);
			/// no backbone hydrogens...

			/// Don't write virtual atoms.
			if ( restype_->atom_type( atom_i ).element() == "X" || restype_->atom_type( atom_j ).element() == "X" ) continue;

			///if ( restype_->atom_is_backbone( atom_i ) && restype_->atom_is_hydrogen( atom_i ) ) continue;
			///if ( restype_->atom_is_backbone( atom_j ) && restype_->atom_is_hydrogen( atom_j ) ) continue;
			std::string const ptmaster = ( is_H ? " 'h'" : "" );
			print_node( ostr, hit.first()[ 1 ], atom_i, coords[ atom_i ], *restype_, tag+"P"+ptmaster);
			print_node( ostr, hit.first()[ 1 ], atom_j, coords[ atom_j ], *restype_, tag+ptmaster);
		}
	}

}

void
SingleDownstreamResidueWriter::set_restype( core::chemical::ResidueTypeCOP restype )
{
	restype_ = restype;
	all_atom_inds_.resize( restype_->natoms() );
	for ( Size ii = 1; ii <= all_atom_inds_.size(); ++ii ) all_atom_inds_[ ii ] = core::id::AtomID( ii, 1 );
}

void
SingleDownstreamResidueWriter::set_downstream_builder( downstream::DownstreamBuilderCOP dsbuilder )
{
	dsbuilder_ = dsbuilder;
}

void
SingleDownstreamResidueWriter::set_downstream_master( std::string const & master )
{
	master_ = master;
}

}
}
}
