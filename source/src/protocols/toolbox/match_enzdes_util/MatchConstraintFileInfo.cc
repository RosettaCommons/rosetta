// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IO-functionality for enzyme and matching Constraints
/// @brief
/// @author Florian Richter, floric@u.washington.edu, may 2009

// Unit headers
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>
#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.hh>
#include <protocols/toolbox/match_enzdes_util/LigandConformer.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <core/pack/rotamer_set/bb_independent_rotamers.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.srlz.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/variant_util.hh>

// Basic headers
#include <basic/basic.hh>

// numeric headers
#include <numeric/HomogeneousTransform.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <basic/Tracer.hh>


// option key includes


#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.toolbox.match_enzdes_util.MatchConstraintFileIfo" );

#ifdef    SERIALIZATION
// Project serialization headers
#include <core/chemical/ResidueType.srlz.hh>
#include <core/chemical/ResidueTypeSet.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


/// @brief function to go through a list of restypes and
/// reduce them to chemically identical ones based on the same base_name
/// i.e. this function gets rid of the variant redundancy
void
add_relevant_restypes_to_subset(
	std::set< core::chemical::ResidueTypeCOP > & restype_subset,
	std::string const & name3,
	core::chemical::ResidueTypeSetCOP restype_set )
{
	using namespace core::chemical;
	core::chemical::ResidueTypeCOPs restypes( restype_set->get_base_types_name3( name3 ) );
	for ( core::Size i = 1; i <= restypes.size(); ++i ) {
		restype_subset.insert( restypes[i] );
	}
}


GeomSampleInfo::GeomSampleInfo(
	std::string tag ) :
	tag_(tag), function_tag_( "default" ), ideal_val_(0.0), tolerance_(0.0), periodicity_(360.0),
	force_const_(0.0), num_steps_(0), step_size_(0.0)
{}

GeomSampleInfo::GeomSampleInfo(
	core::Real ideal_val,
	core::Real tolerance,
	core::Real force_k,
	core::Real periodicity,
	core::Size num_steps
):
	tag_(""),
	function_tag_("default"),
	ideal_val_(ideal_val),
	tolerance_(tolerance),
	periodicity_(periodicity),
	force_const_(force_k),
	num_steps_(num_steps),
	step_size_(0.0)
{ if ( num_steps_ != 0 ) step_size_ = tolerance_ / num_steps_; }

// copy ctor
GeomSampleInfo::GeomSampleInfo( GeomSampleInfo const & gsi ) :
	tag_(gsi.tag_),
	function_tag_(gsi.function_tag_),
	ideal_val_(gsi.ideal_val_),
	tolerance_(gsi.tolerance_),
	periodicity_(gsi.periodicity_),
	force_const_(gsi.force_const_),
	num_steps_(gsi.num_steps_),
	step_size_(gsi.step_size_)
{}

GeomSampleInfo::~GeomSampleInfo() {}

bool
GeomSampleInfo::read_data( std::istringstream & line_stream )
{

	utility::vector1< std::string > fields;

	std::string buffer("");

	while ( true ) {
		line_stream >> buffer;
		if ( line_stream.fail() ) break;
		fields.push_back( buffer );
	}

	if ( fields.size() < 4 ) {

		tr << "Not enough fields detected for constraint " << tag_ << "." << std::endl;
		return false;
	};


	ideal_val_ = (core::Real ) atof( fields[1].c_str() );
	tolerance_ = (core::Real ) atof( fields[2].c_str() );
	force_const_ = (core::Real ) atof( fields[3].c_str() );
	periodicity_ = (core::Real ) atof( fields[4].c_str() );

	if ( (tag_ != "distanceAB:") && ((periodicity_ < 1.0) || (periodicity_ > 360.0)) ) {  //safeguard against stupid input
		std::cerr << "Error: illegal periodicity value of " << periodicity_ <<
			" requested for degree of freedom " << tag_ << "." << std::endl;
		utility_exit_with_message("Illegal periodicity value given. Value must be between 1.0 and 360.0 degrees.");
	}

	if ( fields.size() > 4 ) {

		//we have to test if fields[5] is an integer
		std::istringstream f5;
		f5.clear();
		f5.str( fields[5] );

		f5 >> num_steps_;
		if ( f5.bad() ) num_steps_ = 0;

		//num_steps_ = (core::Size ) atoi( fields[5].c_str() );
		//tr << "for tag " << tag_ << "fields[5] is |" << fields[5] << "| and num_steps is " << num_steps_ << std::endl;
	}

	if ( fields[ fields.size() ] == "PERIODIC" ) function_tag_ = "PERIODIC";

	if ( num_steps_ != 0 ) step_size_ = tolerance_ / num_steps_;

	tr.Debug  << "data read for GeomSampleInfo with tag " << tag_ << ": ideal_value(" << ideal_val_ << "), tolerance(" << tolerance_ << "), force_k(" << force_const_ << "), periodicity(" << periodicity_ << "), num_steps(" << num_steps_ << "), step size(" << step_size_ << "), function_tag(" << function_tag_ << ")." << std::endl;


	//some safeguards against retared user input
	if ( (tolerance_ == 0.0 ) && (num_steps_ != 0 ) ) {

		num_steps_ = 0;

		tr.Warning << "tolerance for constraint " << tag_ << " specified to be 0, yet num_steps specified to be non-0. Ignoring input and setting num_steps to 0." << std::endl;
	}

	return true;

}


utility::vector1< core::Real >
GeomSampleInfo::create_sample_vector() const
{

	core::Size num_ideal_val(1);

	bool distance( tag_ == "distanceAB:");

	//1. figure out the number of ideal values
	if ( (!distance ) && ( periodicity_ != 360.0 ) ) {
		num_ideal_val = static_cast< Size > ( 360.0 / periodicity_ );
	}


	//2. explicity create all the ideal values
	std::list< core::Real > ideal_values;

	//to make sure that there are no duplications
	//through user input of unforeseen weirdness
	std::set< core::Real > seen_values;

	for ( int i =  (int) -( num_ideal_val/2) ; i <= (int) ( num_ideal_val/2); ++i ) {

		core::Real val = ideal_val_ + ( i * periodicity_ );

		if ( !distance ) val = basic::unsigned_periodic_range( val, 360.0 );

		if ( seen_values.find( val ) == seen_values.end() ) {
			ideal_values.push_back( val );
			seen_values.insert( val );
		}
	}

	//sort from lowest to highest because apl sez so
	ideal_values.sort();

	//clear this bc it will be used again
	seen_values.clear();

	//3. build up the diversification samples around each dihedral value
	utility::vector1< core::Real > samples;
	samples.clear();

	tr.Debug << "ideal values for gsi with tag " << tag_ << ", ideal_val " << ideal_val_ << ", and periodicity " << periodicity_ << "are :";

	for ( std::list< core::Real >::const_iterator val_it = ideal_values.begin();
			val_it != ideal_values.end(); ++val_it ) {

		tr.Debug << *val_it << ", ";

		for ( int i = (int) -num_steps_; i <= (int) num_steps_; ++i ) {

			core::Real val =  *val_it + ( i * step_size_ );

			if ( !distance ) val = basic::unsigned_periodic_range( val, 360.0 );

			if ( seen_values.find( val ) == seen_values.end() ) {
				samples.push_back( val );
				seen_values.insert( val );
			}
		}
	} //over all ideal values

	tr.Debug << std::endl << " the generated samples are: ";

	for ( core::Size i = 1; i <= samples.size(); ++i ) {
		tr.Debug << samples[i] << ", ";
	}
	tr.Debug << std::endl;

	return samples;

}


MatchConstraintFileInfo::MatchConstraintFileInfo(
	core::Size index,
	core::chemical::ResidueTypeSetCOP restype_set )
:
	index_( index ),
	is_covalent_(false),
	restype_set_( restype_set ),
	native_ (false)
{
	allowed_seqpos_.clear();
	enz_template_res_.clear();
	constraints_.clear();
}

MatchConstraintFileInfo::~MatchConstraintFileInfo() {}


utility::vector1< core::chemical::ResidueTypeCOP > const
MatchConstraintFileInfo::allowed_restypes( core::Size which_cstres ) const
{
	EnzCstTemplateResCOP template_res = this->enz_cst_template_res( which_cstres );
	return sort_residue_type_pointers_by_name( template_res->allowed_res_types_pointers() );
}

utility::vector1< core::Size > const &
MatchConstraintFileInfo::template_atom_inds(
	core::Size which_cstres,
	core::Size which_template_atom,
	core::chemical::ResidueType const & restype
) const
{
	std::map< core::Size, EnzCstTemplateResOP >::const_iterator map_it =  enz_template_res_.find( which_cstres );
	if ( map_it == enz_template_res_.end() ) {
		utility_exit_with_message( "template res with code blabla not found in MatchConstraintFileInfo ");
	}
	return map_it->second->atom_inds_for_restype( which_template_atom, restype.get_self_ptr() );
}

EnzCstTemplateResCOP
MatchConstraintFileInfo::enz_cst_template_res( core::Size template_res ) const
{
	std::map< core::Size, EnzCstTemplateResOP >::const_iterator map_it =  enz_template_res_.find( template_res );
	if ( map_it == enz_template_res_.end() ) {
		utility_exit_with_message( "template res with code blabla not found in MatchConstraintFileInfo ");
	}
	return map_it->second;
}


//protocols::match::ExternalGeomSamplerCOP
//MatchConstraintFileInfo::exgs() const {
// return exgs_;
//}

std::string matcher_constraint_name_mangler( std::string const & foo ) {
	return "matcher_constraint_" + foo + "_type";
}

std::string matcher_constraint_combination_name_mangler( std::string const & foo ) {
	return "matcher_constraint_combination_" + foo + "_type";
}

void
setup_geometric_attribute_list( utility::tag::AttributeList & attlist ) {
	using namespace utility::tag;

	attlist
		+ XMLSchemaAttribute::required_attribute( "x0",  xsct_real, "x0 specifies the optimum distance x0 for the respective value." )
		+ XMLSchemaAttribute::required_attribute( "xtol",  xsct_real, "xtol specifies the allowed tolerance of the value" )
		+ XMLSchemaAttribute::required_attribute( "k",  xsct_real, "force constant k, or the strength of the parameter. If x is the value of the constrained parameter, the score penalty applied will be: 0 if |x - x0| is less than xtol and k * ( |x - x0| - xtol ) otherwise. This is only relevant for enzdes, and is not used by the matcher." )
		+ XMLSchemaAttribute::required_attribute( "periodicity",  xsct_real, "the periodicity of the constraint. For example, if x0 is 120 and per is 360, the constraint function will have a its minimum at 120 degrees. If x0 is 120 and per is 180, the constraint function will have two minima, one at 120 degrees and one at 300 degrees. If x0 is 120 and per is 120, the constraint function will have 3 minima, at 120, 240, and 360 degrees. In the case of distances, i.e. when specifying the DistanceAB parameter, this value has a special meaning: it indicates whether the constrained interaction is covalent or not with '1' meaning covalent and '0' meaning non-covalent. If the constraint is covalent, Rosetta will not evaluate the vdW term between DownstreamResidue:Atom1 and UpstreamResidue:Atom1 and their [1,3] neighbors." )
		+ XMLSchemaAttribute::required_attribute( "noSamples",  xsct_real, "specifies how many samples the matcher, if using the classic matching algorithm ( see the matcher documentation ), will place between the x0 and x0 +- tol value. If the value in this column is n, the matcher will sample 2n+1 points for the respective parameter. Note: the number of samples is also influenced by the periodicity, with the matcher sampling around every x0." );
}

void
add_subelement_for_constraint_combination(
	utility::tag::XMLSchemaSimpleSubelementList & ssl,
	utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaSimpleSubelementList combination_subelements;

	AttributeList
		combination_attributes, DistanceAB_attributes, AngleA_attributes, AngleB_attributes, TorsionA_attributes,
		TorsionB_attributes, TorsionAB_attributes;

	setup_geometric_attribute_list( DistanceAB_attributes );
	setup_geometric_attribute_list( AngleA_attributes );
	setup_geometric_attribute_list( AngleB_attributes );
	setup_geometric_attribute_list( TorsionA_attributes );
	setup_geometric_attribute_list( TorsionB_attributes );
	setup_geometric_attribute_list( TorsionAB_attributes );

	combination_subelements.add_simple_subelement( "DistanceAB", DistanceAB_attributes , "The distance between UpstreamResidue:Atom1 amd  UpstreamResidue:Atom1 in angstroms." );
	combination_subelements.add_simple_subelement( "AngleA", AngleA_attributes , "The angle formed by DownstreamResidue:Atom2, DownstreamResidue:Atom1, and UpstreamResidue:Atom1 in degrees." );
	combination_subelements.add_simple_subelement( "AngleB", AngleB_attributes , "The angle formed by DownstreamResidue:Atom1, UpstreamResidue:Atom1, UpstreamResidue:Atom2 in degrees." );
	combination_subelements.add_simple_subelement( "TorsionA", TorsionA_attributes , "The dihedral angle formed by DownstreamResidue:Atom3, DownstreamResidue:Atom2, DownstreamResidue:Atom1, and Upstream:Atom1 in degrees." );
	combination_subelements.add_simple_subelement( "TorsionB", TorsionA_attributes , "The dihedral angle formed by DownstreamResidue:Atom1, Upstream:Atom1, Upstream:Atom2, Upstream:Atom3" );
	combination_subelements.add_simple_subelement( "TorsionAB", TorsionA_attributes , "The dihedral angle formed by DownstreamResidue:Atom1, Upstream:Atom1, Upstream:Atom2, Upstream:Atom3" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & matcher_constraint_combination_name_mangler )
		.element_name( "Combination" )
		.description( "A fully-specified combinaiton of the six parameters needed to define an individual constraint. One or more Combinations can be specified in each MatcherConstraint tag." )
		.set_subelements_repeatable( combination_subelements ) //this is probably not the right call but it will work
		.add_attributes( combination_attributes ) //not that there are any
		.write_complex_type_to_schema( xsd ); //this makes it a permanent type

	ssl.add_already_defined_subelement( "Combination", & matcher_constraint_combination_name_mangler );
}

void MatchConstraintFileInfo::return_complex_type_for_MatcherConstraint(
	utility::tag::XMLSchemaSimpleSubelementList & ssl, utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	// attributes for both the Upstream and Downstream subelements
	AttributeList residue_subelement_attributes;
	residue_subelement_attributes
		+ XMLSchemaAttribute::required_attribute( "atom1", xs_string, "name of the first atom to use" )
		+ XMLSchemaAttribute( "atom2", xs_string, "name of the second atom to use" )
		+ XMLSchemaAttribute( "atom3", xs_string, "name of the third atom to use" )
		+ XMLSchemaAttribute( "name", xs_string, "name of the residue to consider for the constraint" );

	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & matcher_constraint_name_mangler );
	subelements.add_simple_subelement( "UpstreamResidue", residue_subelement_attributes,
		"Which atoms are constrained and what type of residue they are in" );
	subelements.add_simple_subelement( "DownstreamResidue", residue_subelement_attributes,
		"Which atoms are constrained and what type of residue they are in" );

	add_subelement_for_constraint_combination( subelements, xsd );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & matcher_constraint_name_mangler )
		.element_name( "MatcherConstraint" )
		.description( "A complete description of an interaction pair to sample, including residue types, atoms to use as reference points to build coordinate frames, and one or Combinations." )
		.set_subelements_repeatable( subelements ) //this is probably not the right call but it will work
		.add_attributes( attlist ) //not that there are any
		.write_complex_type_to_schema( xsd ); //this makes it a permanent type

	ssl.add_already_defined_subelement( "MatcherConstraint", & matcher_constraint_name_mangler );
}


void MatchConstraintFileInfo::add_enzyme_template( core::Size index, utility::tag::TagCOP const tag ) {
	std::string a1 = tag->getOption< std::string >( "atom1", "" );
	std::string a2 = tag->getOption< std::string >( "atom2", "" );
	std::string a3 = tag->getOption< std::string >( "atom3", "" );

	// TODO: edit the XSD to support a comma separated list of residue names
	utility::vector1< std::string > allowed_3res;
	allowed_3res.push_back( tag->getOption< std::string >( "name", "" ) );

	EnzCstTemplateResOP enz_cst_template( new EnzCstTemplateRes( restype_set_ ) );
	enz_cst_template->initialize_from_params( a1, a2, a3, allowed_3res);
	std::pair< core::Size, EnzCstTemplateResOP > to_insert( index, enz_cst_template );
	enz_template_res_.insert( to_insert );
}

void MatchConstraintFileInfo::initialize_from_tag( utility::tag::TagCOP const tag ) {

	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );

	for ( auto const & branch_tag : branch_tags ) {
		if ( branch_tag->getName() == "DownstreamResidue" ) { add_enzyme_template( 1, branch_tag ); }
		else if ( branch_tag->getName() == "UpstreamResidue" ) { add_enzyme_template( 2, branch_tag ); }
		else if ( branch_tag->getName() == "Combination" ) {
			// Fill out the constraint struct!
			SingleConstraint constraint;

			// iterate over the components of the combination to create GeomSampleInfo instancees with the ctor:
			utility::vector1< utility::tag::TagCOP > const combination_tags( branch_tag->getTags() );
			for ( auto const & combination_tag : combination_tags ) {
				core::Real x0 = combination_tag->getOption< core::Real >( "x0", 0. );
				core::Real xtol = combination_tag->getOption< core::Real >( "xtol", 0. );
				core::Real k = combination_tag->getOption< core::Real >( "k", 0. );
				core::Real periodicity = combination_tag->getOption< core::Real >( "periodicity", 0. );
				core::Size noSamples = combination_tag->getOption< core::Size >( "noSamples", 0 );
				GeomSampleInfoOP gsi( new GeomSampleInfo( x0, xtol, k, periodicity, noSamples ) );

				gsi->tag( combination_tag->getName() );  // this will let us keep track of which gsi is which

				if ( combination_tag->getName() == "DistanceAB" ) { gsi->tag( "distanceAB:" ); constraint.dis_U1D1 = gsi; }
				else if ( combination_tag->getName() == "AngleA" ) { constraint.ang_U1D2 = gsi; }
				else if ( combination_tag->getName() == "AngleB" ) { constraint.ang_U2D1 = gsi; }
				else if ( combination_tag->getName() == "TorsionA" ) { constraint.tor_U1D3 = gsi; }
				else if ( combination_tag->getName() == "TorsionB" ) { constraint.tor_U3D1 = gsi; }
				else if ( combination_tag->getName() == "TorsionAB" ) { constraint.tor_U2D2 = gsi; }
			}
			constraints_.push_back( constraint );
		}
	}
}

bool
MatchConstraintFileInfo::read_data( utility::io::izstream & data )
{
	std::istringstream line_stream;

	std::string line, key(""), tag(""),res3;
	core::Size map_id(0);

	//std::cerr << "calling read data for mcfi " << std::endl;

	SingleConstraint current_constraint;

	while ( !data.eof() ) {

		key = ""; tag = "";
		getline(data,line);

		utility::vector1< std::string > comment_split = utility::string_split( line, '#' );
		if ( comment_split[1] == "" ) continue;
		line_stream.clear();
		line_stream.str( comment_split[1] );
		line_stream >> key;

		//std::cerr << "reading things, line is " << line << ", key is " << key;
		//Kui Native 110809
		if ( key == "NATIVE" ) {
			native_ = true;
		} else if ( key == "TEMPLATE::" ) {
			line_stream >> tag;
			//tr.Info << "tag is: " << tag << " ";
			if ( tag == "ATOM_MAP:" ) {

				line_stream >> map_id;

				std::map< core::Size, EnzCstTemplateResOP >::iterator map_it = enz_template_res_.find( map_id );

				if ( map_it == enz_template_res_.end() ) {

					std::pair< core::Size, EnzCstTemplateResOP > to_insert( map_id, EnzCstTemplateResOP( new EnzCstTemplateRes( restype_set_ ) ) );

					enz_template_res_.insert( to_insert );

					map_it = enz_template_res_.find( map_id );
					map_it->second->set_param_index( map_id );
				}

				map_it->second->read_params( line_stream );
			}

			//std::cerr << "  end of file line, tag was " << tag << std::endl;
		} else if ( key == "CONSTRAINT::" ) {
			line_stream >> tag;

			GeomSampleInfoOP gs_info( new GeomSampleInfo( tag ) );

			if ( !gs_info->read_data( line_stream ) ) return false;

			if ( tag == "distanceAB:" ) {

				current_constraint.dis_U1D1 = gs_info;

				//old convention to declare covalency in file
				if ( current_constraint.dis_U1D1->periodicity() == 1.0 ) is_covalent_ = true;
				else is_covalent_ = false;
			} else if ( tag == "angle_A:" ) current_constraint.ang_U1D2 = gs_info;

			else if ( tag == "angle_B:" ) current_constraint.ang_U2D1 = gs_info;

			else if ( tag == "torsion_A:" ) current_constraint.tor_U1D3 = gs_info;

			else if ( tag == "torsion_AB:" ) current_constraint.tor_U2D2 = gs_info;

			else if ( tag == "torsion_B:" ) current_constraint.tor_U3D1 = gs_info;

			else {
				std::cerr << "The following line in the cst file with key " << key << " was not recognized and will be ignored: " << std::endl << line << std::endl;
			}

			//std::cerr << "  end of file line, tag was " << tag << std::endl;

		} else if ( key == "ALGORITHM_INFO::" ) { //if key==CONSTRAINT

			line_stream >> tag;

			if ( !this->process_algorithm_info( tag, data ) ) return false;

		} else if ( key == "CST::END" ) {
			constraints_.push_back( current_constraint );
			current_constraint = SingleConstraint();
			return true;
		} else if ( key != "" ) {
			std::cerr << "The following line in the cst file with key " << key << " was not recognized and will be ignored: " << std::endl << line << std::endl;
		}

	} //while ( !data.eof )

	//if we get to here, that means the cstfile is corrupted

	return false;

} //read_data


void
MatchConstraintFileInfo::process_data()
{

	core::chemical::ResidueTypeSetCOP restype_set( restype_set_ );

	for ( std::map< core::Size, EnzCstTemplateResOP >::iterator map_it = enz_template_res_.begin();
			map_it != enz_template_res_.end(); ++map_it ) {

		utility::vector1< std::string > const & res_name3s =
			map_it->second->allowed_res_types();
		std::set< core::chemical::ResidueTypeCOP > restypes_this_res;
		for ( core::Size j = 1; j <= res_name3s.size(); ++j ) {
			add_relevant_restypes_to_subset( restypes_this_res, res_name3s[j], restype_set_ );
		}

		for ( std::set< core::chemical::ResidueTypeCOP >::iterator
				set_it  = restypes_this_res.begin();
				set_it != restypes_this_res.end(); ++set_it ) {
			if ( map_it->second->compatible_restype( *set_it ) ) {
				map_it->second->determine_atom_inds_for_restype( *set_it );
			}
		}

		if ( map_it->second->atom_inds_for_restype_begin() == map_it->second->atom_inds_for_restype_end() ) {
			tr.Warning << "No appropriate atoms found for entry " << map_it->first << " in constraint " << index_ << std::endl;
		}
	} //loop over all Enz_cst_template-res

}


std::list< core::conformation::ResidueCOP >
MatchConstraintFileInfo::inverse_rotamers_against_residue(
	core::Size const target_template,
	core::conformation::ResidueCOP target_conf
) const
{
	runtime_assert( enz_template_res_.size() == 2 );
	core::Size const invrot_template( target_template == 1 ? 2 : 1 );

	//the exgs created based on the cst file info might be wrong if we have to create inverse rotamers
	//of a residue that is the upstream residue in the cstfile
	bool flip_exgs_upstream_downstream_samples( false );
	if ( invrot_template == this->upstream_res() ) flip_exgs_upstream_downstream_samples = true;

	std::list< core::conformation::ResidueCOP > to_return;
	core::Size rotcount_buffer(0);

	//1. loop over all allowed restypes of the inverse template
	utility::vector1< core::chemical::ResidueTypeCOP > invrot_restypes(this->allowed_restypes( invrot_template ));

	//if we're dealing with backbone interaction, only build glycine rotamers
	bool backbone_interaction(false);
	if ( this->is_backbone( invrot_template ) ) {
		core::chemical::ResidueTypeSetCOP restype_set( restype_set_ );
		backbone_interaction = true;
		invrot_restypes.clear();
		invrot_restypes.push_back( restype_set->name_mapOP("ALA") );
		tr << "Only Ala inverse rotamers will be built because it is a backbone interaction." << std::endl;
	}
	for ( core::Size ii =1; ii <= invrot_restypes.size(); ++ii ) {

		//2 get the relevant atoms through wich the orientation of the two residues
		//is defined in the constraint file
		utility::vector1< utility::vector1< core::Size > > target_template_atom_inds(3), invrot_template_atom_inds(3);
		for ( core::Size atct = 1; atct <= 3; ++atct ) {
			target_template_atom_inds[atct] = this->template_atom_inds( target_template, atct, target_conf->type() );
			invrot_template_atom_inds[atct] = this->template_atom_inds( invrot_template, atct, *(invrot_restypes[ii]) );
		}

		//3. loop over all possible combinations of atoms in the present residue and the inverse rotamer
		for ( core::Size jj = 1; jj <= target_template_atom_inds[1].size(); ++jj ) {
			utility::vector1< core::Size > targ_ats(3);
			targ_ats[1] = target_template_atom_inds[1][jj]; targ_ats[2] = target_template_atom_inds[2][jj]; targ_ats[3] = target_template_atom_inds[3][jj];
			for ( core::Size kk = 1; kk <= invrot_template_atom_inds[1].size(); ++kk ) {

				utility::vector1< core::Size > invrot_ats(3);
				invrot_ats[1] = invrot_template_atom_inds[1][kk]; invrot_ats[2] = invrot_template_atom_inds[2][kk]; invrot_ats[3] = invrot_template_atom_inds[3][kk];

				//4. hand off to other function so code stays readable
				std::list<core::conformation::ResidueCOP > inv_rots_this_combo = this->inverse_rotamers_against_residue( *target_conf, invrot_restypes[ii], targ_ats, invrot_ats, flip_exgs_upstream_downstream_samples, backbone_interaction );
				to_return.splice( to_return.end(), inv_rots_this_combo  );

			} // kk loop over all possible atoms in the inverse rotamer
		} //jj loop over all possible atoms in the present residue
		tr << to_return.size() - rotcount_buffer << " inverse rotamers were created for restype " << invrot_restypes[ii]->name() << "." << std::endl;
		rotcount_buffer = to_return.size();
	} //ii loop over all allowed restypes of the inverse template
	return to_return;
}

std::list< core::conformation::ResidueCOP >
MatchConstraintFileInfo::inverse_rotamers_against_residue(
	core::conformation::Residue const & target_conf,
	core::chemical::ResidueTypeCOP invrot_restype,
	utility::vector1< core::Size > const & target_ats,
	utility::vector1< core::Size > const & invrot_ats,
	bool const flip_exgs_upstream_downstream_samples,
	bool const backbone_interaction
) const
{
	//using namespace protocols::match::downstream;

	std::list< core::conformation::ResidueCOP > to_return;

	utility::vector1< core::conformation::ResidueCOP > rotamers( core::pack::rotamer_set::bb_independent_rotamers( invrot_restype, true ) ); //This is getting the residue specific inverse rotamers

	// bool no_theozyme_inverse_rotamers( basic::options::option[ basic::options::OptionKeys::enzdes::no_theozyme_inverse_rotamers ]() );
	// if ( no_theozyme_inverse_rotamers ) get rid of rotamers ;

	runtime_assert( rotamers.size() > 0 );
	tr << rotamers.size() << " bbindependent rotamers for Residue " << rotamers[1]->type().name() << "." << std::endl;

	//note: if we have a backbone interaction, this means we need to diversify
	//the phi value of the rotamer
	if ( backbone_interaction ) {
		runtime_assert( (rotamers.size() == 1) && (rotamers[1]->name3() == "ALA") );
		this->diversify_backbone_only_rotamers( rotamers );
	}
	core::Size inv_oat1(0), inv_oat2(0), inv_oat3(0);
	rotamers[1]->select_orient_atoms( inv_oat1, inv_oat2, inv_oat3 );

	utility::vector1< LigandConformer > invrot_conformers;
	for ( core::Size rotcount(1); rotcount <= rotamers.size(); ++rotcount ) {
		invrot_conformers.push_back( LigandConformer() );
		invrot_conformers[ rotcount ].initialize_from_residue( invrot_ats[1], invrot_ats[2], invrot_ats[3], inv_oat1, inv_oat2, inv_oat3, *(rotamers[rotcount]));
	}

	utility::vector1< ExternalGeomSamplerOP > exgs_list( this->create_exgs() );

	core::Real atom1_atom2_distance = invrot_conformers[1].atom1_atom2_distance();
	core::Real atom2_atom3_distance = invrot_conformers[1].atom2_atom3_distance();
	core::Real atom1_atom2_atom3_angle = invrot_conformers[1].atom1_atom2_atom3_angle();

	for ( auto & exgs : exgs_list ) {

		//apparently we have to do some stuff with the sampler
		exgs->set_dis_D1D2( atom1_atom2_distance );
		exgs->set_dis_D2D3( atom2_atom3_distance );
		exgs->set_ang_D1D2D3( atom1_atom2_atom3_angle );
		if ( flip_exgs_upstream_downstream_samples ) exgs->flip_upstream_downstream_samples();
		exgs->precompute_transforms();

		HTReal ht_start( target_conf.xyz(target_ats[3]), target_conf.xyz(target_ats[2]), target_conf.xyz(target_ats[1]) );

		for ( Size ii = 1; ii <= exgs->n_tor_U3D1_samples(); ++ii ) {
			HTReal ht_ii = ht_start * exgs->transform( HT_tor_U3D1, ii );
			for ( Size jj = 1; jj <= exgs->n_ang_U2D1_samples(); ++jj ) {
				HTReal ht_jj = ht_ii * exgs->transform( HT_ang_U2D1, jj );
				for ( Size kk = 1; kk <= exgs->n_dis_U1D1_samples(); ++kk ) {
					HTReal ht_kk = ht_jj;
					ht_kk.walk_along_z( exgs->dis_U1D1_samples()[ kk ] );
					for ( Size ll = 1; ll <= exgs->n_tor_U2D2_samples(); ++ll ) {
						HTReal ht_ll = ht_kk * exgs->transform( HT_tor_U2D2, ll );
						for ( Size mm = 1; mm <= exgs->n_ang_U1D2_samples(); ++mm ) {
							HTReal ht_mm = ht_ll * exgs->transform( HT_ang_U1D2, mm );
							for ( Size nn = 1; nn <= exgs->n_tor_U1D3_samples(); ++nn ) {
								HTReal ht_nn = ht_mm * exgs->transform( HT_tor_U1D3, nn );

								for ( core::Size rotcount(1); rotcount <= rotamers.size(); ++rotcount ) {
									core::conformation::ResidueOP rot( new core::conformation::Residue( *(rotamers[rotcount]) ) );
									for ( core::Size atm = 1; atm <= rot->natoms(); ++atm ) {
										rot->set_xyz( atm, invrot_conformers[rotcount].coordinate_in_D3_frame( atm, ht_nn ) );
									}
									to_return.push_back( rot );
								} //loop over all rotamers
							} //nn sampler
						} //mm sampler
					}// ll sampler
				} //kk sampler
			} // jj sampler
		} //ii sampler
	} // iteration over all ExternalGeomSamplers
	return to_return;
}

/// @details
/// helper function to keep code readable
/// for rotamers that make backbone interactions,
/// as opposed to sidechain interactions only, the default
/// phi (-150) that comes out of the bb-indep rotamers function
/// has an influence on what fragments in sampling can overlap
/// with this rotamer. thus we'll put in more samples to allow
/// for more diversity
/// the implementation is quite clumsy, make a one residue pose
/// add the chainbreak variant, set the chi, return the residue
/// but there's no easier way to simply rotate around a bond.
/// additional samples will be put at a phi of -60 and 70,
/// i.e. other regions observed in ramachandran plot
void
MatchConstraintFileInfo::diversify_backbone_only_rotamers( utility::vector1< core::conformation::ResidueCOP > & rotamers ) const
{
	//core::conformation::ResidueOP changeres( rotamers[1]->clone() );
	core::pose::Pose dummy_pose;
	dummy_pose.append_residue_by_jump( *(rotamers[1]), (core::Size) 0 );
	core::pose::add_variant_type_to_pose_residue( dummy_pose, core::chemical::CUTPOINT_UPPER, 1 );
	dummy_pose.set_phi( 1, -60.0 );
	rotamers.push_back( core::pose::remove_variant_type_from_residue( dummy_pose.residue(1), core::chemical::CUTPOINT_UPPER, dummy_pose ) );
	dummy_pose.set_phi( 1, 70.0 );
	rotamers.push_back( core::pose::remove_variant_type_from_residue( dummy_pose.residue(1), core::chemical::CUTPOINT_UPPER, dummy_pose ) );
}

bool
MatchConstraintFileInfo::process_algorithm_info(
	std::string tag,
	utility::io::izstream & data
) {

	// according to apl's request, only allow prespecified tags
	// if you want to read in data for additional tags, you have
	// to specify those here
	if ( ( tag != "match") && ( tag != "match_positions" ) && ( tag != "test") && ( tag != "invrot_tree" ) ) {
		utility_exit_with_message("Tag "+tag+" not a legal option for ALGORITHM_INFO block.");
	}

	if ( algorithm_inputs_.find( tag ) != algorithm_inputs_.end() ) {
		tr.Error << "tag " << tag << " was found twice in the same cstfile block." << std::endl;
		return false;
	}

	utility::vector1< std::string > alg_strings;

	while ( !data.eof() ) {

		std::string line("");
		getline(data,line);

		//if ( line == "ALGORITHM_INFO::END") {
		if ( utility::trimmed_compare( line, "ALGORITHM_INFO::END") ) {
			if ( alg_strings.size() != 0 ) {
				algorithm_inputs_.insert( std::pair< std::string, utility::vector1< std::string > >( tag, alg_strings ) );
			} else tr.Warning << "ALGORITHM_INFO block for " << tag << " seemed to contain no information." << std::endl;
			return true;
		}
		utility::vector1< std::string > comment_split = utility::string_split( line, '#' );
		if ( comment_split[1] != "" ) alg_strings.push_back( comment_split[1] );

	} //while ( !data.eof() ) {

	tr << "Error, when reading algorithm info block with tag " << tag << ", no ALGORITHM_INFO::END line was found." << std::endl;

	return false;

} //process algorithm_info


utility::vector1< ExternalGeomSamplerOP >
MatchConstraintFileInfo::create_exgs() const
{
	// So if we're doing this with a list of explicit samples, we should return a vector of EXGSs
	utility::vector1< ExternalGeomSamplerOP > exgs_list;
	for ( auto const & constraint : constraints_ ) {
		ExternalGeomSamplerOP exgs( new ExternalGeomSampler() );

		utility::vector1< std::string > tags_undefined_gsi;

		if ( constraint.dis_U1D1 ) exgs->set_dis_U1D1_samples( constraint.dis_U1D1->create_sample_vector() );
		else tags_undefined_gsi.push_back( "distanceAB:" );

		if ( constraint.ang_U1D2 ) exgs->set_ang_U1D2_samples( constraint.ang_U1D2->create_sample_vector() );
		else tags_undefined_gsi.push_back( "angle_A:" );

		if ( constraint.ang_U2D1 ) exgs->set_ang_U2D1_samples( constraint.ang_U2D1->create_sample_vector() );
		else tags_undefined_gsi.push_back( "angle_B:" );

		if ( constraint.tor_U1D3 ) exgs->set_tor_U1D3_samples( constraint.tor_U1D3->create_sample_vector() );
		else tags_undefined_gsi.push_back( "torsion_A:" );

		if ( constraint.tor_U2D2 ) exgs->set_tor_U2D2_samples( constraint.tor_U2D2->create_sample_vector() );
		else tags_undefined_gsi.push_back( "torsion_AB:" );

		if ( constraint.tor_U3D1 ) exgs->set_tor_U3D1_samples( constraint.tor_U3D1->create_sample_vector() );
		else tags_undefined_gsi.push_back( "torsion_B:" );

		if ( tags_undefined_gsi.size() != 0 ) {

			tr.Warning << "could not create external geom sampler from file input because not all 6 necessary degrees of freedom are specified.\n The following DOFs are missing specifications: ";

			for ( core::Size i = 1; i <= tags_undefined_gsi.size(); ++i ) {
				tr.Warning << tags_undefined_gsi[i] << ", ";
			}
			tr.Warning << "." << std::endl;

			//std::cerr << "setting external geom sampler to null pointer" << std::endl;
			exgs = NULL;
		}
		if ( exgs ) {
			exgs_list.push_back( exgs );
		}
	}
	return exgs_list;
}


MatchConstraintFileInfoList::MatchConstraintFileInfoList(
	core::chemical::ResidueTypeSetCOP restype_set
) :
	restype_set_( restype_set )
{
	mcfis_.clear();
}

MatchConstraintFileInfoList::~MatchConstraintFileInfoList() {}

/// @brief temporary implementation for now, only one MCFI supported
bool
MatchConstraintFileInfoList::read_data( utility::io::izstream & data )
{

	//std::cerr << "calling read data for mcfi list " << std::endl;
	//active_mcfi_ = 1;

	//mcfis_.clear();

	core::Size new_index = mcfis_.size() + 1;
	MatchConstraintFileInfoOP mcfi( new MatchConstraintFileInfo( new_index, restype_set_ ) );

	if ( mcfi->read_data( data ) ) {
		add_mcfi( mcfi );
		return true;
	}
	return false;
}

void MatchConstraintFileInfoList::add_mcfi( MatchConstraintFileInfoOP mcfi )
{
	mcfi->process_data();
	mcfis_.push_back( mcfi );
	determine_upstream_restypes();
}


utility::vector1< MatchConstraintFileInfoCOP > const &
MatchConstraintFileInfoList::mcfis_for_upstream_restype( core::chemical::ResidueTypeCOP restype ) const
{
	std::map< core::chemical::ResidueTypeCOP, utility::vector1< MatchConstraintFileInfoCOP > >::const_iterator
		mcfi_it = mcfis_for_restype_.find( restype );
	if ( mcfi_it == mcfis_for_restype_.end() ) {
		utility_exit_with_message( " could not find mcfi list for given restype" );
	}
	return mcfi_it->second;
}

std::list< core::conformation::ResidueCOP >
MatchConstraintFileInfoList::inverse_rotamers_against_residue(
	core::Size const target_template,
	core::conformation::ResidueCOP target_conf ) const
{
	std::list< core::conformation::ResidueCOP > to_return;
	for ( core::Size i = 1; i <= mcfis_.size(); ++i ) {
		if ( mcfis_[i]->num_enz_cst_template_res() != 2 ) {
			tr << "Can't create inverse rotamers for mcfi " << i << " because it has more or less than 2 template res." << std::endl;
			continue;
		}
		//core::Size const invrot_template( target_template == 1 ? 2 : 1 );
		if ( std::find(
				mcfis_[i]->allowed_res_name3s( target_template ).begin(),
				mcfis_[i]->allowed_res_name3s( target_template ).end(),
				target_conf->name3() )
				== mcfis_[i]->allowed_res_name3s( target_template ).end() ) {
			tr << "Can't create inverse rotamers for mcfi " << i << " because it doesn't contain target template for residue " << target_conf->name3() << "." << std::endl;
			continue;
		}
		std::list< core::conformation::ResidueCOP > mcfi_invrots(
			mcfis_[i]->inverse_rotamers_against_residue( target_template, target_conf ) );
		to_return.splice( to_return.end(), mcfi_invrots  );
	}
	return to_return;
}


void
MatchConstraintFileInfoList::determine_upstream_restypes()
{

	upstream_restypes_.clear();
	mcfis_for_restype_.clear();

	utility::vector1< std::string > gly_vec;
	gly_vec.push_back("GLY");

	utility::vector1< core::chemical::ResidueTypeCOP > restype_temp_set;
	core::chemical::ResidueTypeSetCOP restype_set( restype_set_ );

	for ( core::Size i =1 ; i <= mcfis_.size(); ++i ) {
		bool is_backbone( mcfis_[i]->is_backbone( mcfis_[i]->upstream_res() ) );

		//in case the mcfi is a backbone interaction, only allow glycine as the restype
		//for this mcfi
		utility::vector1< std::string > const & res_name3s( is_backbone ? gly_vec : mcfis_[i]->allowed_res_name3s( mcfis_[i]->upstream_res() ) );

		std::set< core::chemical::ResidueTypeCOP > restypes_this_mcfi;

		for ( core::Size j = 1; j <= res_name3s.size(); ++j ) {
			add_relevant_restypes_to_subset( restypes_this_mcfi, res_name3s[j], restype_set_ );
		}

		//add the restypes_this_mcfi to the total set of upstream residue
		//(if they haven't already been added)
		for ( std::set< core::chemical::ResidueTypeCOP >::iterator set_it = restypes_this_mcfi.begin();
				set_it != restypes_this_mcfi.end(); ++set_it ) {
			//build the restype_to_mcfi mapping
			std::map< core::chemical::ResidueTypeCOP, utility::vector1< MatchConstraintFileInfoCOP > >::iterator
				res_mcfi_it = mcfis_for_restype_.find( *set_it );
			if ( res_mcfi_it == mcfis_for_restype_.end() ) {
				std::pair<core::chemical::ResidueTypeCOP, utility::vector1< MatchConstraintFileInfoCOP > >
					to_insert ( *set_it, utility::vector1< MatchConstraintFileInfoCOP >() );
				mcfis_for_restype_.insert( to_insert );
				res_mcfi_it = mcfis_for_restype_.find( *set_it );
			}

			res_mcfi_it->second.push_back( mcfis_[i] );
			//restype_to_mcfi mapping updated

			if ( std::find( restype_temp_set.begin(), restype_temp_set.end(), *set_it ) == restype_temp_set.end() ) {
				restype_temp_set.push_back( *set_it );
			}

		} //loop over restypes this mcfi

	} //loop over all mcfis


	// This is a new sort which will hopefully reduce instability in match_* integration tests.
	upstream_restypes_ = sort_residue_type_pointers_by_name( restype_temp_set );
}


}
}//enzdes
}//protocols


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::toolbox::match_enzdes_util::MatchConstraintFileInfoList::MatchConstraintFileInfoList() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::MatchConstraintFileInfoList::save( Archive & arc ) const {
	arc( CEREAL_NVP( mcfis_ ) ); // utility::vector1<MatchConstraintFileInfoOP>
	arc( CEREAL_NVP( upstream_restypes_ ) );

	arc( mcfis_for_restype_.size() );
	for ( std::map<core::chemical::ResidueTypeCOP, utility::vector1<MatchConstraintFileInfoCOP> >::const_iterator
			iter = mcfis_for_restype_.begin(), iter_end = mcfis_for_restype_.end(); iter != iter_end; ++iter ) {
		arc( iter->first );
		arc( iter->second );
	}

	arc( CEREAL_NVP( restype_set_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::MatchConstraintFileInfoList::load( Archive & arc ) {
	arc( mcfis_ ); // utility::vector1<MatchConstraintFileInfoOP>
	arc( upstream_restypes_ );

	core::Size nrestypes; arc( nrestypes );
	for ( core::Size ii = 1; ii <= nrestypes; ++ii ) {
		core::chemical::ResidueTypeCOP restype;
		arc( restype );
		utility::vector1< MatchConstraintFileInfoOP > mcfis;
		arc( mcfis );
		mcfis_for_restype_[ restype ] = mcfis;
	}

	arc( restype_set_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::match_enzdes_util::MatchConstraintFileInfoList );
CEREAL_REGISTER_TYPE( protocols::toolbox::match_enzdes_util::MatchConstraintFileInfoList )


/// @brief Default constructor required by cereal to deserialize this class
protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo::MatchConstraintFileInfo() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( index_ ) ); // core::Size
	arc( CEREAL_NVP( allowed_seqpos_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( enz_template_res_ ) ); // std::map<core::Size, EnzCstTemplateResOP>
	arc( CEREAL_NVP( is_covalent_ ) ); // _Bool
	arc( CEREAL_NVP( constraints_ ) ); // utility::vector1< SingleConstraint >
	arc( CEREAL_NVP( algorithm_inputs_ ) ); // std::map<std::string, utility::vector1<std::string> >
	core::chemical::serialize_residue_type_set( arc, restype_set_ );
	arc( CEREAL_NVP( native_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo::load( Archive & arc ) {
	arc( index_ ); // core::Size
	arc( allowed_seqpos_ ); // utility::vector1<core::Size>
	arc( enz_template_res_ ); // std::map<core::Size, EnzCstTemplateResOP>
	arc( is_covalent_ ); // _Bool
	arc( constraints_ ); // utility::vector1< SingleConstraint >
	arc( algorithm_inputs_ ); // std::map<std::string, utility::vector1<std::string> >

	core::chemical::deserialize_residue_type_set( arc, restype_set_ );

	arc( native_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo );
CEREAL_REGISTER_TYPE( protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo )


/// @brief SingleConstraint serialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::SingleConstraint::save( Archive & arc ) const {
	arc( CEREAL_NVP( dis_U1D1 ) ); // GeomSampleInfoCOP
	arc( CEREAL_NVP( ang_U1D2 ) ); // GeomSampleInfoCOP
	arc( CEREAL_NVP( ang_U2D1 ) ); // GeomSampleInfoCOP
	arc( CEREAL_NVP( tor_U1D3 ) ); // GeomSampleInfoCOP
	arc( CEREAL_NVP( tor_U2D2 ) ); // GeomSampleInfoCOP
	arc( CEREAL_NVP( tor_U3D1 ) ); // GeomSampleInfoCOP
}

/// @brief SingleConstraint deserialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::SingleConstraint::load( Archive & arc ) {
	// make a local non-const pointer to a GeomSampleInfo
	GeomSampleInfoOP local_gis;

	arc( local_gis ); // GeomSampleInfoCOP masquerading as a GeomSampleInfoOP
	dis_U1D1 = local_gis; // copy the address of the non-const pointer into the const pointer

	// repeat the process for each data member
	arc( local_gis );
	ang_U1D2 = local_gis;

	arc( local_gis );
	ang_U2D1 = local_gis;

	arc( local_gis );
	tor_U1D3 = local_gis;

	arc( local_gis );
	tor_U2D2 = local_gis;

	arc( local_gis );
	tor_U3D1 = local_gis;
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::match_enzdes_util::SingleConstraint );

/// @brief Default constructor required by cereal to deserialize this class
protocols::toolbox::match_enzdes_util::GeomSampleInfo::GeomSampleInfo() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::GeomSampleInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( tag_ ) ); // std::string
	arc( CEREAL_NVP( function_tag_ ) ); // std::string
	arc( CEREAL_NVP( ideal_val_ ) ); // core::Real
	arc( CEREAL_NVP( tolerance_ ) ); // core::Real
	arc( CEREAL_NVP( periodicity_ ) ); // core::Real
	arc( CEREAL_NVP( force_const_ ) ); // core::Real
	arc( CEREAL_NVP( num_steps_ ) ); // core::Size
	arc( CEREAL_NVP( step_size_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::GeomSampleInfo::load( Archive & arc ) {
	arc( tag_ ); // std::string
	arc( function_tag_ ); // std::string
	arc( ideal_val_ ); // core::Real
	arc( tolerance_ ); // core::Real
	arc( periodicity_ ); // core::Real
	arc( force_const_ ); // core::Real
	arc( num_steps_ ); // core::Size
	arc( step_size_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::match_enzdes_util::GeomSampleInfo );
CEREAL_REGISTER_TYPE( protocols::toolbox::match_enzdes_util::GeomSampleInfo )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_match_enzdes_util_MatchConstraintFileInfo )
#endif // SERIALIZATION

