// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/rotamer_gen/RDKitRotamers.hh>
#include <protocols/rotamer_gen/RDKitRotamersCreator.hh>

#include <core/chemical/rdkit/util.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <protocols/chemistries/util.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <basic/Tracer.hh>

#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/DistGeomHelpers/Embedder.h>
#include <rdkit/GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <rdkit/GraphMol/Conformer.h>

namespace protocols {
namespace rotamer_gen {

static basic::Tracer TR("protocol.rotamer_gen.RDKitRotamers");

/////////////////////////////////////////////

std::string
RDKitRotamersCreator::keyname() const {
	return RDKitRotamers::class_name();
}

protocols::chemistries::ChemistryOP
RDKitRotamersCreator::create_chemistry() const {
	return protocols::chemistries::ChemistryOP( new RDKitRotamers );
}

void
RDKitRotamersCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	RDKitRotamers::provide_xml_schema( xsd );
}

////////////////////////////////////////////

void RDKitRotamers::apply(core::chemical::MutableResidueType & rsdtype) {

	//because fragment generation takes awhile, we only really want to do this once. If we have created
	//rotamers for this residuetype already, we should probably not do it again.
	//But if you're calling this Chemistry explicitly, you probably mean to
	if ( rsdtype.rotamer_library_specification() ) {
		TR << "Resetting rotamers for " << rsdtype.name() << " using RDKit." << std::endl;
	}

	// nconf logic adapted from Jean-Paul Ebejer's presentation
	// at the London RDKit User General Meeting
	// http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
	core::Size nconf( 50 );
	if ( rsdtype.nchi() > 12 ) {
		nconf = 300;
	} else if ( rsdtype.nchi() >= 8 ) {
		nconf = 200;
	}

	core::chemical::rotamers::StoredRotamerLibrarySpecificationOP rotamers_spec( new core::chemical::rotamers::StoredRotamerLibrarySpecification() );

	rotamers_spec->set_rotamers( generate_conformers( rsdtype, nconf ) );

	rsdtype.rotamer_library_specification( rotamers_spec );

	// Unless we crash, we always succeed, so there isn't a reason to change the success code.
}

utility::vector1< std::map< std::string, core::Vector > >
RDKitRotamers::generate_conformers( core::chemical::MutableResidueType const & rsdtype, core::Size nconf, bool minimize ) {

	utility::vector1< std::map< std::string, core::Vector > > retval;

	std::clock_t begin = std::clock();

	core::chemical::rdkit::RestypeToRDMol to_rdmol(rsdtype, /*neutralize = */ false, /*keep_hydro=*/ true);
	::RDKit::RWMOL_SPTR rdmol( to_rdmol.Mol() );

	if ( !rdmol ) {
		TR.Warning << "Can't convert restype " << rsdtype.name() << " to rdmol - skipping RDKit conformer generation." << std::endl;
		return retval;
	}

	core::chemical::VDIndexMapping const & rosetta_to_rdkit( to_rdmol.vd_to_index() );


	//::RDKit::MolOps::addHs( *rdmol, false, true ); // Not needed, because we keep them from above.

	// Generate the rotamers
	::RDKit::INT_VECT conf_ids; // List of generated conformations
	::RDKit::DGeomHelpers::EmbedMultipleConfs(*rdmol, conf_ids, nconf,
		/*numThreads*/ 1, /*maxIterations*/ 30, /*seed*/ 11111, /*clearConfs*/ true, /*useRandomCoords*/ false, /*boxSizeMult*/ 2.0,   // We use a fixed seed such that conformer generation is repeatable.
		/*randNegEig*/ true, /*numZeroFail*/ 1, /*pruneRmsThresh*/ -1.0, /*coordMap*/ nullptr, /*optimizerForceTol*/ 1e-3, /*ignoreSmoothingFailures*/ false,
		/*enforceChirality*/ true, /*useExpTorsionAnglePrefs*/ true, /*useBasicKnowledge*/ true  // These are the ones we're changing
	);

	TR << "Embedding took " << double(std::clock()-begin)/ CLOCKS_PER_SEC << "s " << std::endl;

	if ( minimize ) {
		std::vector< std::pair<int, double> > status; // Output status for each conf
		::RDKit::UFF::UFFOptimizeMoleculeConfs(*rdmol, status);
		//(We ignore the outputs status, because we're lazy.)
		TR << "Embedding && minimization took " << double(std::clock()-begin)/ CLOCKS_PER_SEC << "s " << std::endl;
	}

	// Now convert the coordinates to rotamers for Rosetta
	for ( core::Size ii(0); ii < conf_ids.size(); ++ii ) {

		::RDKit::Conformer const & conf( rdmol->getConformer(conf_ids[ii]) );
		std::map< std::string, core::Vector > single_rotamer_spec;

		for ( core::chemical::VD atm: rsdtype.all_atoms() ) {
			if ( rosetta_to_rdkit[ atm ] != rosetta_to_rdkit.invalid_entry() ) {
				::RDGeom::Point3D const & pos( conf.getAtomPos( rosetta_to_rdkit[ atm ] ) );
				//set the xyz coordinates in Rosetta
				single_rotamer_spec[ rsdtype.atom_name(atm) ] = core::Vector( pos.x, pos.y, pos.z );
			}
		}

		retval.push_back( single_rotamer_spec );
	}

	TR << "Generated " << retval.size() << " rotamers for " << rsdtype.name() << " in " << double(std::clock()-begin)/ CLOCKS_PER_SEC << "s total." << std::endl;

	return retval;
}

std::string
RDKitRotamers::class_name() {
	return "RDKitRotamers";
}

void
RDKitRotamers::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	protocols::chemistries::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"The RDKitRotamers chemistry will use the conformer generation protocol of "
		"Ebejer et al. (doi:10.1021/ci2004658) to generate a rotameric library "
		"and attach it to the given ResidueType.",
		attributes );
}

void
RDKitRotamers::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &)
{
	// Previous settings no longer used - ignore all.
}


}
}
