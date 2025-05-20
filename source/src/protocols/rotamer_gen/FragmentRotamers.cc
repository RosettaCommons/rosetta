// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/rotamer_gen/FragmentRotamers.hh>
#include <protocols/rotamer_gen/FragmentRotamersCreator.hh>

#include <core/chemical/bcl/util.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/bcl/BCLFragmentHandler.hh>
#include <protocols/chemistries/util.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.hh>

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


#ifdef USEBCL
#include <io/bcl_io_file.h>
#include <chemistry/bcl_chemistry_fragment_ensemble.h>
#include <chemistry/bcl_chemistry_sample_conformations.h>
#include <linal/bcl_linal_vector_3d_operations.h>
#include <util/bcl_util_si_ptr_vector.h>
//#include <chemistry/bcl_chemistry_rotamer_library_interface.h>
#include <chemistry/bcl_chemistry_rotamer_library_file.h>

#endif

namespace protocols {
namespace rotamer_gen {

static basic::Tracer TR("protocol.rotamer_gen.FragmentRotamers");

/////////////////////////////////////////////

std::string
FragmentRotamersCreator::keyname() const {
	return FragmentRotamers::class_name();
}

protocols::chemistries::ChemistryOP
FragmentRotamersCreator::create_chemistry() const {
	return protocols::chemistries::ChemistryOP( new FragmentRotamers );
}

void
FragmentRotamersCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	FragmentRotamers::provide_xml_schema( xsd );
}

////////////////////////////////////////////

void FragmentRotamers::apply(core::chemical::MutableResidueType & rsdtype) {

#ifndef USEBCL
	TR.Error << "Cannot generate BCL Fragment rotamers for residue type " << rsdtype.name() << std::endl;
	core::chemical::bcl::require_bcl();
#else
	//because fragment generation takes awhile, we only really want to do this once. If we have created
	//rotamers for this residuetype already, we should probably not do it again.
	//But if you're calling this Chemistry explicitly, you probably mean to
	if( rsdtype.rotamer_library_specification() ) {
		TR << "Resetting rotamers for " << rsdtype.name() << " using the BCL." << std::endl;
	}
	std::clock_t begin = std::clock();

	core::chemical::bcl::BCLFragmentHandler frag_handler;
	::bcl::chemistry::FragmentComplete fragment = frag_handler.restype_to_fragment( rsdtype );

	core::chemical::IndexVDMapping bcl_to_rosetta( frag_handler.get_vd_to_index().reverse() );

	if(!conformation_sampler_.IsDefined()){
		::bcl::chemistry::RotamerLibraryFile rotamer_lib_file;
		rotamer_lib_file.AssertRead( ::bcl::util::ObjectDataLabel( "(number_of_files=100)"));
		// The default "SymmetryRMSD" comparer is *phenomenally* slow for some residues. Use an alternative comparer to get reasonable speeds.
		conformation_sampler_ = ::bcl::util::ShPtr< ::bcl::chemistry::SampleConformations>( new ::bcl::chemistry::SampleConformations( rotamer_lib_file, "DihedralBins", /*tolerance*/ 1, /*number_conf*/ 100, /*max_iter*/ 200, /*change_chirality*/ false ) );
	}

	::bcl::storage::Pair< ::bcl::chemistry::FragmentEnsemble, ::bcl::chemistry::FragmentEnsemble> output( conformation_sampler_->operator()(fragment));

	//use the conformation sampler to create all the fragments!
	::bcl::chemistry::FragmentEnsemble &all_fragments(  output.First());

	core::chemical::rotamers::StoredRotamerLibrarySpecificationOP rotamers_spec( new core::chemical::rotamers::StoredRotamerLibrarySpecification() );

	//This is where we convert the ensemble of rotamers created in the BCL to residues (rotamers) for Rosetta
	for
	(
			::bcl::storage::List< ::bcl::chemistry::FragmentComplete>::iterator itr_mols( all_fragments.Begin()), itr_mols_end( all_fragments.End());
			itr_mols != itr_mols_end;
			++itr_mols
	)
	{
		std::map< std::string, core::Vector > single_rotamer_spec;

		::bcl::util::SiPtrVector< const ::bcl::linal::Vector3D> atom_coords( itr_mols->GetAtomCoordinates());

		//keep track of what the atom number we are on
		core::Size i=0;

		for
		(
				::bcl::util::SiPtrVector< const ::bcl::linal::Vector3D>::const_iterator itr_coords( atom_coords.Begin()), itr_coords_end( atom_coords.End());
				itr_coords != itr_coords_end;
				++itr_coords, ++i
		)
		{
			//convert BCL fragment index base to Rosetta's VD based atoms
			//This relies on the fact that the BCL manipulation doesn't change indicies
			core::chemical::VD vd = bcl_to_rosetta[i];

			//get the xyz coordinates
			core::Real x, y, z;
			x= ( *itr_coords)->X();
			y = ( *itr_coords)->Y();
			z = ( *itr_coords)->Z();

			//set the xyz coordinates in Rosetta
			single_rotamer_spec[ rsdtype.atom(vd).name() ] = core::Vector( x, y, z );
		}
		rotamers_spec->add_rotamer( single_rotamer_spec );
	}

	rsdtype.rotamer_library_specification( rotamers_spec );

	TR << " number of rotamers created for "<< rsdtype.name() << ": " << rotamers_spec->coordinates().size() << std::endl;

	TR << "Rotamer generation took: " << double(std::clock()-begin)/ CLOCKS_PER_SEC << "s total" << std::endl;

#endif

	// Unless we crash, we always succeed, so there isn't a reason to change the success code.
}

std::string
FragmentRotamers::class_name() {
	return "FragmentRotamers";
}

void
FragmentRotamers::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	protocols::chemistries::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"The FragmentRotamers chemistry will call the BCL to generate a rotameric library "
		"for the ResidueType using the protocol of Kothiwale et al. (doi:10.1186/s13321-015-0095-1) "
		"and attach it to the ResidueType. "
		"Attempting to use this will result in a runtime error unless Rosetta was successfully compiled with `extras=bcl`.",
		attributes );
}

void
FragmentRotamers::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &)
{
	// Previous settings no longer used - ignore all.
}


}
}
