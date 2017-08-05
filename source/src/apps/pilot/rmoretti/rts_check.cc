// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/rmoretti/rts_check.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/residue_io.hh>

#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ResidueTypeKinWriter.hh>

#include <core/types.hh>

#include <devel/init.hh>

#include <utility/options/FileVectorOption.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <fstream>
#include <string>
#include <set>
using namespace basic::options;

static THREAD_LOCAL basic::Tracer TR("apps.rosetta_retype_check");

int
main( int argc, char * argv [] )
{

	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;

		devel::init(argc, argv);

		std::string typeset_name("fa_standard");
		if ( basic::options::option[ basic::options::OptionKeys::in::file::centroid ] ) {
			typeset_name = "centroid";
		}
		if ( basic::options::option[ basic::options::OptionKeys::in::file::extra_res_database ].user() ) {
			typeset_name = basic::options::option[ basic::options::OptionKeys::in::file::extra_res_database ];
		}

		TR << "Loading residue type set " << typeset_name << std::endl;
		core::chemical::ResidueTypeSetCAP rts( core::chemical::ChemicalManager::get_instance()->residue_type_set( typeset_name ) );

		//ResidueTypeSelector dna_selector;
		//dna_selector.set_property("DNA").exclude_variants();
		ResidueTypeCOPs dna_types = core::chemical::ResidueTypeFinder( *rts ).base_property( core::chemical::DNA ).get_possible_base_residue_types();

		TR << "Number of DNA types: " << dna_types.size();
		for ( core::Size j=1; j<=dna_types.size(); ++j ) {
			TR << "        " << dna_types[ j ]->name();
		}
		TR << std::endl;

		ResidueTypeCOPs dna_types2;
		dna_types2.push_back( rts->get_representative_type_aa( core::chemical::na_ade ) );
		dna_types2.push_back( rts->get_representative_type_aa( core::chemical::na_thy ) );
		dna_types2.push_back( rts->get_representative_type_aa( core::chemical::na_cyt ) );
		dna_types2.push_back( rts->get_representative_type_aa( core::chemical::na_gua ) );

		// TR << "Number of DNA types2: " << dna_types.size();
		// for (core::Size j=1; j<=dna_types2.size(); ++j ) {
		//  if( dna_types2[ j ] ) {
		//   TR << "        " << dna_types2[ j ]->name();
		//  } else {
		//   TR << "        " << "xxNULLxx";
		//  }
		// }
		// TR << std::endl;


		core::chemical::ResidueTypeCOPs const & all_rsd_types( core::chemical::ResidueTypeFinder( *rts ).get_all_possible_residue_types() );

		std::set< std::string > name3set;
		std::set< core::chemical::AA > aaset;
		std::set< char > name1set;

		for ( core::Size j=1; j<=all_rsd_types.size(); ++j ) {
			name3set.insert( all_rsd_types[j]->name3() );
			name1set.insert( all_rsd_types[j]->name1() );
			aaset.insert( all_rsd_types[j]->aa() );
			//  core::chemical::VariantTypeList vts( all_rsd_types[j]->variant_types() );
			//  while( vts.size() ) {
			//   core::chemical::VariantType vt( vts.back() );
			//   vts.pop_back();
			//   for( core::Size vv(1); vv <= vts.size(); ++vv ) {
			//    if( vt == vts[vv] ) {
			//  TR << all_rsd_types[j]->name() << " :: ";
			//     TR << vt << "  ";
			//  TR << std::endl;
			//    }
			//   }
			//  }
		}

		core::chemical::VariantTypeList vtl;

		TR << "Availible name3: ";
		for ( std::set< std::string >::const_iterator iter( name3set.begin() ); iter != name3set.end(); ++iter ) {
			TR << std::endl;
			TR << *iter << "::   ";
			core::chemical::ResidueTypeCOPs const & restypes( rts->get_all_types_with_variants_name3( *iter, vtl ) );

			for ( core::Size ii(1); ii <= restypes.size(); ++ii ) {
				TR << restypes[ii]->name() << "  ";
				//if( restypes[ii]->name().substr(0,5)=="HIS_D" ) { // restypes[ii]->name3() == *iter ) {
				// TR << restypes[ii]->name() << " ( " << restypes[ii]->name3() << " ) ";
				//}
			}
			TR << std::endl << std::endl;

		}
		TR << std::endl;


		TR << "Availible name1: ";
		for ( std::set< char >::const_iterator iter( name1set.begin() ); iter != name1set.end(); ++iter ) {
			TR << std::endl;
			TR << *iter << "::   ";
			core::chemical::ResidueTypeCOPs const & restypes( rts->get_all_types_with_variants_name1( *iter, vtl ) );

			for ( core::Size ii(1); ii <= restypes.size(); ++ii ) {
				TR << restypes[ii]->name() << "  ";
				//if( restypes[ii]->name().substr(0,5)=="HIS_D" ) { // restypes[ii]->name3() == *iter ) {
				// TR << restypes[ii]->name() << " ( " << restypes[ii]->name3() << " ) ";
				//}
			}
			TR << std::endl << std::endl;

		}
		TR << std::endl;


		TR << "Availible aa: ";
		for ( std::set< core::chemical::AA >::const_iterator iter( aaset.begin() ); iter != aaset.end(); ++iter ) {
			TR << std::endl;
			TR << *iter << "::   ";
			core::chemical::ResidueTypeCOPs const & restypes( rts->get_all_types_with_variants_aa( *iter, vtl ) );

			for ( core::Size ii(1); ii <= restypes.size(); ++ii ) {
				TR << restypes[ii]->name() << "  ";
				//if( restypes[ii]->name().substr(0,5)=="HIS_D" ) { // restypes[ii]->name3() == *iter ) {
				// TR << restypes[ii]->name() << " ( " << restypes[ii]->name3() << " ) ";
				//}
			}
			TR << std::endl << std::endl;

		}
		TR << std::endl;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

