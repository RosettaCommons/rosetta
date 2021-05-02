// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/movers/RNA_Decoarsify.cc
/// @brief Make a coarse RNA pose into a fullatom representation.
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Unit headers
#include <protocols/rna/movers/RNA_Decoarsify.hh>
#include <protocols/rna/movers/RNA_DecoarsifyCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/rna/util.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/rna/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>

#include <protocols/rna/movers/RNA_Coarsify.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.rna.movers.RNA_Decoarsify" );

namespace protocols {
namespace rna {
namespace movers {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
RNA_Decoarsify::RNA_Decoarsify():
	protocols::moves::Mover( RNA_Decoarsify::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
RNA_Decoarsify::RNA_Decoarsify( RNA_Decoarsify const & src ):
	protocols::moves::Mover( src )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
RNA_Decoarsify::~RNA_Decoarsify(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
RNA_Decoarsify::apply( core::pose::Pose& pose_coarse ){

	// To decoarsify pose_coarse...
	// 1. make a fullatom pose with the same sequence: pose_scratch
	// 2. coarsify it: pose_scratch_coarsened
	// 3. transform pose_scratch_coarsened to look like pose_coarse
	// 4. apply the same transformations to pose_scratch to look like
	//    a fa version of pose_coarse

	using namespace core;
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace core::pose;

	typedef numeric::xyzVector< Real > Vector;

	auto rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );
	auto rsd_set_full = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// initialize helix pose
	// Pose pose_coarse;
	// std::string infile  = option[ in::file::s ][1];
	// import_pose::pose_from_file( pose_coarse, *rsd_set_coarse, infile , core::import_pose::PDB_file);
	core::pose::PoseOP pose = utility::pointer::make_shared< Pose >();

	core::pose::PoseOP pose_scratch = utility::pointer::make_shared< Pose >();
	core::pose::PoseOP pose_scratch_coarsened = utility::pointer::make_shared< Pose >();

	make_pose_from_sequence( *pose_scratch, pose_coarse.annotated_sequence(), *rsd_set_full );

	make_pose_from_sequence( *pose_scratch_coarsened, pose_coarse.annotated_sequence(), *rsd_set_full );
	auto coarsify = utility::pointer::make_shared< protocols::rna::movers::RNA_Coarsify >();
	coarsify->apply( *pose_scratch_coarsened );

	*pose = *pose_scratch;

	for ( Size ii = 1; ii <= pose->size(); ++ii ) {
		// A couple coordinate systems to allow easy superposition.
		Vector const origin1 =  pose_coarse.xyz( NamedAtomID( " P  ", ii ) );
		Vector const x1 =  pose_coarse.xyz( NamedAtomID( " P  ", ii ) );
		Vector const y1 =  pose_coarse.xyz( NamedAtomID( " S  ", ii ) );
		Vector const z1 =  pose_coarse.xyz( NamedAtomID( " CEN", ii ) );
		Stub stub1( origin1, x1, y1, z1 );

		// A couple coordinate systems to allow easy superposition.
		Vector const origin2 =  pose_scratch_coarsened->xyz( NamedAtomID( " P  ", ii ) );
		Vector const x2 =  pose_scratch_coarsened->xyz( NamedAtomID( " P  ", ii ) );
		Vector const y2 =  pose_scratch_coarsened->xyz( NamedAtomID( " S  ", ii ) );
		Vector const z2 =  pose_scratch_coarsened->xyz( NamedAtomID( " CEN", ii ) );
		Stub stub2( origin2, x2, y2, z2 );

		for ( Size jj = 1; jj <= pose->residue_type( ii ).natoms(); ++jj ) {
			pose->set_xyz( AtomID( jj, ii ), stub1.local2global( stub2.global2local( pose_scratch->xyz( AtomID( jj, ii) ) ) ) );
		}

	}

	//std::string outfile = "full.pdb";
	//if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();
	//std::cout<< "dumping to " << outfile << std::endl;
	pose_coarse = *pose;

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
RNA_Decoarsify::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
RNA_Decoarsify::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void RNA_Decoarsify::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
RNA_Decoarsify::fresh_instance() const
{
	return utility::pointer::make_shared< RNA_Decoarsify >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
RNA_Decoarsify::clone() const
{
	return utility::pointer::make_shared< RNA_Decoarsify >( *this );
}

std::string RNA_Decoarsify::get_name() const {
	return mover_name();
}

std::string RNA_Decoarsify::mover_name() {
	return "RNA_Decoarsify";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
RNA_DecoarsifyCreator::create_mover() const
{
	return utility::pointer::make_shared< RNA_Decoarsify >();
}

std::string
RNA_DecoarsifyCreator::keyname() const
{
	return RNA_Decoarsify::mover_name();
}

void RNA_DecoarsifyCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RNA_Decoarsify::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, RNA_Decoarsify const & mover )
{
	mover.show(os);
	return os;
}

} //movers
} //rna
} //protocols
