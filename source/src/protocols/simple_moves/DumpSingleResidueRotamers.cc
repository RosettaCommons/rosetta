// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/DumpSingleResidueRotamers.cc
/// @brief Given a residue index, dump all of the rotamers to individual PDB files within 0-1 sd of the mean
/// @author Rebecca Alford (rfalford12@gmail.com)

// Unit headers
#include <protocols/simple_moves/DumpSingleResidueRotamers.hh>
#include <protocols/simple_moves/DumpSingleResidueRotamersCreator.hh>

// Protocol Headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/types.hh>

#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

// Basic/Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/ralford.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
#include <utility/graph/Graph.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <fstream>

static basic::Tracer TR( "protocols.simple_moves.DumpSingleResidueRotamers" );

namespace protocols {
namespace simple_moves {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
DumpSingleResidueRotamers::DumpSingleResidueRotamers():
	protocols::moves::Mover( DumpSingleResidueRotamers::mover_name() ),
	rsd_index_( 1 /* magic number */ ),
	prefix_( "output" ),
	write_rotamers_to_pdbs_( false ),
	all_positions_( false )
{
	init_from_options();
}

/// @brief Custom constructor
DumpSingleResidueRotamers::DumpSingleResidueRotamers( core::Size rsd_index ) :
	protocols::moves::Mover( DumpSingleResidueRotamers::mover_name() ),
	rsd_index_( rsd_index ),
	prefix_( "output" ),
	write_rotamers_to_pdbs_( false ),
	all_positions_( false )
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
DumpSingleResidueRotamers::DumpSingleResidueRotamers( DumpSingleResidueRotamers const & src ):
	protocols::moves::Mover( src ),
	rsd_index_( src.rsd_index_ ),
	prefix_( src.prefix_ ),
	write_rotamers_to_pdbs_( src.write_rotamers_to_pdbs_ ),
	all_positions_( src.all_positions_ )
{
	init_from_options();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
DumpSingleResidueRotamers::~DumpSingleResidueRotamers() = default;

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
DumpSingleResidueRotamers::apply( core::pose::Pose & pose ){

	using namespace core::pack::dunbrack;
	using namespace core::pack::rotamers;
	using namespace core::chemical;
	using namespace core::pose;

	runtime_assert( rsd_index_ >= 1 && rsd_index_ <= pose.total_residue() );

	// Set up the single residue dunbrack library
	ResidueTypeCOP concrete_residue( ResidueTypeOP( new ResidueType( pose.residue( rsd_index_ ).type() ) ) );
	SingleResidueRotamerLibraryCOP rsd_rl_temp1( SingleResidueRotamerLibraryFactory::get_instance()->get( *concrete_residue ) );
	SingleResidueDunbrackLibraryCOP rsd_rl = utility::pointer::dynamic_pointer_cast< const SingleResidueDunbrackLibrary >( rsd_rl_temp1 );

	if ( all_positions_ ) {

		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			utility::vector1< PoseOP > rotamer_poses = enumerate_aa_rotamer( rsd_rl, pose, ii );
			if ( write_rotamers_to_pdbs_ ) {
				write_poses_to_pdbs( rotamer_poses, ii );
			}
		}

	} else {

		utility::vector1< PoseOP > rotamer_poses = enumerate_aa_rotamer( rsd_rl, pose, rsd_index_ );
		if ( write_rotamers_to_pdbs_ ) {
			write_poses_to_pdbs( rotamer_poses, rsd_index_ );
		}

	}

} // apply

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
DumpSingleResidueRotamers::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
DumpSingleResidueRotamers::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	if ( tag->hasOption( "rsd_index" ) ) {
		rsd_index_ = tag->getOption< core::Size >( "rsd_index" );
	}

	if ( tag->hasOption( "all_positions" ) ) {
		all_positions_ = tag->getOption< bool >( "all_positions" );
	}

	if ( tag->hasOption( "write_rotamers_to_pdbs" ) ) {
		write_rotamers_to_pdbs_ = tag->getOption< bool >( "wrtie_rotamers_to_pdbs" );
	}

	if ( tag->hasOption( "prefix" ) ) {
		prefix_ = tag->getOption< std::string >( "prefix" );
	}

}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
DumpSingleResidueRotamers::fresh_instance() const
{
	return protocols::moves::MoverOP( new DumpSingleResidueRotamers );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
DumpSingleResidueRotamers::clone() const
{
	return protocols::moves::MoverOP( new DumpSingleResidueRotamers( *this ) );
}

std::string DumpSingleResidueRotamers::get_name() const {
	return mover_name();
}

std::string DumpSingleResidueRotamers::mover_name() {
	return "DumpSingleResidueRotamers";
}

void DumpSingleResidueRotamers::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "rsd_index", xsct_non_negative_integer, "Position at which to enumerate rotamers") +
		XMLSchemaAttribute( "prefix", xs_string, "Prefix for output rotamer info file and pdbs when specified" ) +
		XMLSchemaAttribute( "write_rotamers_to_pdbs", xsct_rosetta_bool, "Should we write the rotamers out to individual PDB files?" ) +
		XMLSchemaAttribute( "all_positions", xsct_rosetta_bool, "Dump rotamers at all pdb positions" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Dump all of the rotamers for a given amino acid type at the rsd_index position into individual PDB files", attlist );

}

/////////////// Creator ///////////////

protocols::moves::MoverOP
DumpSingleResidueRotamersCreator::create_mover() const
{
	return protocols::moves::MoverOP( new DumpSingleResidueRotamers );
}

std::string
DumpSingleResidueRotamersCreator::keyname() const
{
	return DumpSingleResidueRotamers::mover_name();
}

void DumpSingleResidueRotamersCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DumpSingleResidueRotamers::provide_xml_schema( xsd );
}

utility::vector1< core::pose::PoseOP >
DumpSingleResidueRotamers::enumerate_aa_rotamer(
	core::pack::dunbrack::SingleResidueDunbrackLibraryCOP rsd_rl,
	core::pose::Pose & pose,
	core::Size rsd_index ) const {

	using namespace core::pose;

	// Make a vector
	utility::vector1< core::pose::PoseOP > rotamer_poses;

	// Grab the current rotamer set given phi and psi
	utility::fixedsizearray1< core::Real, 5 > bbs;
	bbs[1] = pose.phi( rsd_index );
	bbs[2] = pose.psi( rsd_index );
	utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > drsd( rsd_rl->get_all_rotamer_samples( bbs) );

	// Write rotamer information to an output filename
	std::string rsd_aa( pose.residue( rsd_index ).name3() );
	std::string rsdpos( std::to_string( rsd_index ) );
	std::string rotamer_info_fn( prefix_ + "_" + rsd_aa + rsdpos + "_rotamer_info.txt" );
	std::ofstream rotamer_info( rotamer_info_fn );
	rotamer_info << "NRotamer Probability Stdev Chi1 Chi2 Chi3 Chi4" << std::endl;

	// Loop through the rotamers for a given residue
	for ( core::Size nrot = 1; nrot <= drsd.size(); ++nrot ) {

		// For each rotamer, store within one standard deviation of the mean
		for ( core::Real isd = -1; isd <= 1; ++isd ) {
			PoseOP pose_copy = PoseOP( new Pose( pose ) );
			core::Size nchi_in_aa( pose.residue( rsd_index ).type().nchi() -  pose.residue( rsd_index ).type().n_proton_chi() );
			core::Size nchi_tor( 4 ); // Max number of chi angles in canonical AAs
			rotamer_info << nrot << " " << drsd[nrot].probability() << " " << isd;
			for ( core::Size nchi = 1; nchi <= nchi_tor; ++nchi ) {

				core::Real stdev( drsd[ nrot ].chi_sd()[ nchi ] );
				core::Real chi( drsd[ nrot ].chi_mean()[ nchi ] + stdev*isd );

				if ( nchi <= nchi_in_aa ) {
					pose_copy->set_chi( nchi, rsd_index, chi );
				}

				rotamer_info << " " << chi;
				if ( nchi == nchi_tor ) {
					rotamer_info << std::endl;
				}
			}
			rotamer_poses.push_back( pose_copy );
		}
	}
	return rotamer_poses;
}

void
DumpSingleResidueRotamers::write_poses_to_pdbs(
	utility::vector1< core::pose::PoseOP > poses, core::Size position ) {

	// Stdev types
	utility::vector1< std::string > stdev_types;
	stdev_types.push_back( "minus1" );
	stdev_types.push_back( "zero" );
	stdev_types.push_back( "plus1" );

	for ( core::Size ii = 1; ii <= poses.size()/3; ++ii ) {
		std::string stdev( std::to_string( ii % 3 ) );
		std::string rsd_pos_name( std::to_string( position ) );
		for ( core::Size jj = 0; jj <= 2; ++jj ) {
			std::string aa_name( poses[ ii + jj  ]->residue( position ).name3() );
			std::string fn( prefix_ + "_" + aa_name + rsd_pos_name + "_" + std::to_string( ii ) + "_" + stdev_types[ jj + 1 ] + ".pdb");
			poses[ ii + jj ]->dump_pdb( fn );
		}
	}

}

void DumpSingleResidueRotamers::init_from_options() {

	using namespace basic::options;
	if ( option[ OptionKeys::ralford::dump_rotamers::rsd_index ].user() ) {
		rsd_index_ = option[ OptionKeys::ralford::dump_rotamers::rsd_index ];
	} else {
		utility_exit_with_message( "Residue index at which to sample rotamers not specified! Exiting..." );
	}

	if ( option[ OptionKeys::ralford::dump_rotamers::rotamer_info_prefix ].user() ) {
		prefix_ = option[ OptionKeys::ralford::dump_rotamers::rotamer_info_prefix ];
	}

	if ( option[ OptionKeys::ralford::dump_rotamers::write_rotamers_to_pdbs ].user() ) {
		write_rotamers_to_pdbs_ = option[ OptionKeys::ralford::dump_rotamers::write_rotamers_to_pdbs  ];
	}

	if ( option[ OptionKeys::ralford::dump_rotamers::all_positions ].user() ) {
		all_positions_ = option[ OptionKeys::ralford::dump_rotamers::all_positions ];
	}
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

std::ostream &
operator<<( std::ostream & os, DumpSingleResidueRotamers const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //simple_moves
