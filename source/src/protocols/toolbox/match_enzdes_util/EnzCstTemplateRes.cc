// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file IO-functionality for enzyme Constraints
/// @brief
/// @author Florian Richter, floric@u.washington.edu

// Unit headers
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh> //needed to convert one letter AA codes
#include <core/chemical/ResidueTypeSet.hh> //have to include complete file
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <basic/options/option.hh>
#include <core/id/SequenceMapping.hh>


// Utility Headers
#include <utility/string_util.hh>
#include <iostream>
#include <string>
#include <sstream>

#include <basic/Tracer.hh>


// option key includes

#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.toolbox.match_enzdes_util.EnzCstTemplateRes" );

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

/// @details Auto-generated virtual destructor
EnzCstTemplateRes::~EnzCstTemplateRes() {}

/// @details Auto-generated virtual destructor
EnzCstTemplateResAtoms::~EnzCstTemplateResAtoms() {}


void
EnzCstTemplateResAtoms::remap_resid( core::id::SequenceMapping const & smap )
{
	remap_atomid_vector( smap, atom1_ );
	remap_atomid_vector( smap, atom2_ );
	remap_atomid_vector( smap, atom3_ );
}

void
EnzCstTemplateResAtoms::remap_atomid_vector(
	core::id::SequenceMapping const & smap,
	std::vector< core::id::AtomID > & atomid_vec
)
{
	for ( std::vector< core::id::AtomID >::iterator vec_it = atomid_vec.begin();
			vec_it != atomid_vec.end(); ++vec_it ) {

		core::Size newpos = smap[ vec_it->rsd() ];
		if ( newpos == 0 ) utility_exit_with_message("A catalytic residue is apparently missing from the pose");

		*vec_it = core::id::AtomID( vec_it->atomno(), newpos );
	}
}


EnzCstTemplateRes::EnzCstTemplateRes(
	core::chemical::ResidueTypeSetCAP src_restype_set
) : rb_minimizable_(true), is_backbone_(false),
	identical_tag_found_(false), corresponding_res_block_(0),
	corresponding_res_num_in_block_(0),
	restype_set_(src_restype_set), enz_io_param_( /* NULL */ )
{
	clear_all();
}


EnzCstTemplateRes::EnzCstTemplateRes(
	core::chemical::ResidueTypeSetCAP src_restype_set,
	EnzConstraintParametersCAP src_enzio_param )
:
	rb_minimizable_(true),
	is_backbone_(false),
	identical_tag_found_(false),
	corresponding_res_block_(0),
	corresponding_res_num_in_block_(0),
	restype_set_(src_restype_set),
	enz_io_param_( src_enzio_param)
{
	clear_all();
}


/// @brief WARNING: currently doesn't copy the template atoms in the respos map
EnzCstTemplateRes::EnzCstTemplateRes(
	EnzCstTemplateResCOP other,
	EnzConstraintParametersCAP new_ref_param)
:
	//respos_map_(other->respos_map_), //old, remove
	atom1_(other->atom1_),
	atom2_(other->atom2_),
	atom3_(other->atom3_),
	at1_type_(other->at1_type_),
	at2_type_(other->at2_type_),
	at3_type_(other->at3_type_),
	allowed_res_types_(other->allowed_res_types_),
	allowed_res_types_pointers_(other->allowed_res_types_pointers_),
	atom_inds_for_restype_( other->atom_inds_for_restype_),
	rb_minimizable_( other->rb_minimizable_ ),
	is_backbone_( other->is_backbone_ ),
	respos_from_external_(other->respos_from_external_),
	identical_tag_found_(other->identical_tag_found_),
	corresponding_res_block_(other->corresponding_res_block_),
	corresponding_res_num_in_block_(other->corresponding_res_num_in_block_),
	restype_set_(other->restype_set_),
	enz_io_param_(new_ref_param),
	param_index_(other->param_index_)
{
}

/// @brief read and set up the general information for one residue of one pair block in a cstfile
void
EnzCstTemplateRes::read_params(std::istringstream & line_stream)
{
	using namespace core::chemical;

	rb_minimizable_ = true; //true by default

	//input from cst file
	std::string tag = "";
	std::string allowed_1res_raw;
	std::vector< std::string > allowed_3res_raw;
	std::string buffer = "";
	core::Size size_buffer(0);


	line_stream >> tag;
	//tr.Info << "template_res.read tag is: " << tag << " ";
	if ( tag == "atom_name:" ) {
		std::string a1, a2, a3;
		line_stream >> a1 >> a2 >> a3 ;
		atom1_.push_back(a1);
		atom2_.push_back(a2);
		atom3_.push_back(a3);
		// else if (tag == "atom_name") line_stream >> skip(1) >> bite(4, parameters->res1_at1name_) >> skip(1);
	} else if ( tag == "atom_type:" ) {
		line_stream >> at1_type_ >> at2_type_ >> at3_type_;
	} else if ( tag == "residue1:" ) {
		line_stream >> allowed_1res_raw;
	} else if ( tag == "residue3:" ) {
		while ( !line_stream.fail() ) { line_stream >> buffer; allowed_3res_raw.push_back(buffer);}
	} else if ( tag == "identical:" ) {
		line_stream >> corresponding_res_block_ >> corresponding_res_num_in_block_;
		identical_tag_found_ = true;
	} else if ( tag == "seqpos:" ) {
		while ( !line_stream.fail() ) {
			line_stream >> size_buffer;
			respos_from_external_.push_back( size_buffer );
		}
	} else if ( tag == "no_rb_min" ) rb_minimizable_ = false;
	else if ( tag == "is_backbone" ) is_backbone_ = true;

	else {
		std::cerr << "Line in cstfile specifying template residue with tag " << tag << " was not recognized and will be ignored. " << std::endl;
	}

	//done reading parameters, now the processing starts:
	//1. make sure the ResidueTypeSet knows about the residues specified in the input file
	//2. if atom types have been put in, make sure they only have 4 letters.
	//3. if this is an identical line, make sure that the two residues declared identical
	//   have the same allowed res types

	core::chemical::ResidueTypeSetCOP restype_set( restype_set_ );

	//first process 3-letter code input
	for ( std::vector< std::string >::iterator it = allowed_3res_raw.begin(); it != allowed_3res_raw.end(); ++it ) {
		if ( it->size() == 2 ) *it = " " + *it;
		if ( restype_set->has_name3( *it ) ) {
			allowed_res_types_.push_back( *it );
		} else {
			utility_exit_with_message("Error in cstfile: Residue with 3-letter code "+*it+" is unknown.");
		}
	}
	//and now do the same for 1-letter code input (only canonical aas will work)
	for ( Size ii = 0; ii != allowed_1res_raw.size(); ii++ ) {

		if ( oneletter_code_specifies_aa( allowed_1res_raw[ii] ) ) {

			std::string cur_res_name3 = name_from_aa( aa_from_oneletter_code( allowed_1res_raw[ii] ) );
			if ( restype_set->has_name3( cur_res_name3 ) ) allowed_res_types_.push_back( cur_res_name3 );
			else {
				utility_exit_with_message("Unexpected error in program: Residue "+cur_res_name3+" is unknown. ResidueTypeSet setup seems to not have worked properly.");
			}
		} else {
			//std::string fock( (char*) &allowed_1res_raw[ii] );
			std::string fock( 1,allowed_1res_raw[ii] );
			utility_exit_with_message("Error in cstfile: Residue with one letter code "+ fock +" is unknown.");
		}
	}

	//make sure atom names are capped to four characters
	if ( at1_type_ != "" ) { at1_type_ = at1_type_.substr(0,4);}
	if ( at2_type_ != "" ) { at2_type_ = at2_type_.substr(0,4);}
	if ( at3_type_ != "" ) { at3_type_ = at3_type_.substr(0,4);}

} //EnzCstTemplateRes::read_params


/// @brief show the contents of a particular instance
void EnzCstTemplateRes::show_params() const
{
	using namespace core::chemical;

	tr.Info << "Parameters read for template residue: atom ids are :";
	for ( core::Size i = 1; i <= atom1_.size(); ++i ) tr.Info << " (" << atom1_[i] << " " << atom2_[i] << " " << atom3_[i] << ")";
	tr.Info << ", first atom name is " << at1_type_  <<" and allowed residues are: ";
	for ( utility::vector1< std::string >::const_iterator res_it = allowed_res_types_.begin(); res_it != allowed_res_types_.end(); ++res_it ) {
		tr.Info << *res_it << ", ";
	}
	tr.Info <<  std::endl ;
}


/// @brief check whether the data gathered from the cst file and the pdbfile/pose is consistent in itself and with one another
void EnzCstTemplateRes::get_pose_data(core::pose::Pose & pose) const {

	using namespace core::chemical;
	EnzCstTemplateResCacheOP template_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( enz_io_param_.lock()->cst_block() )->template_res_cache( param_index_ ) );

	//first check whether we actually have a residue to process new
	if ( template_cache->seqpos_map_.size() == 0 ) {
		std::cerr << "Error: no residues for this template." << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}


	//then we check if this acutally needs to be done (cause if it's done twice we get fucked up behaviour)
	if ( template_cache->pose_data_uptodate_ == true ) return;

	//then check whether residues specified in the pdbfile were in the list specified in the cstfile,
	//and for each residue, check whether the atoms specified in the cst file actually match atoms in this residue type
	for ( std::map< Size, EnzCstTemplateResAtomsOP >::iterator respos_it = template_cache->seqpos_map_.begin(); respos_it != template_cache->seqpos_map_.end(); ++respos_it ) {

		ResidueTypeCOP cur_res = pose.residue_type( respos_it->first ).get_self_ptr();
		std::string cur_res_name3 = pose.residue( respos_it->first ).name3();
		//utility::trim( cur_res_name3 );

		//tr.Info << "resis in pose " << pose.residue( respos_it->first ).name3() << " " << respos_it->first <<";";

		if ( std::find( allowed_res_types_.begin(), allowed_res_types_.end(), cur_res_name3 ) == allowed_res_types_.end() ) {
			std::cerr << "Error: residue " << pose.residue( respos_it->first).name3() << respos_it->first << "found in pdb header is not allowed by data in cstfile." << std::endl;
			std::cerr << "Allowed restypes:";
			for ( utility::vector1< std::string >::const_iterator iter = allowed_res_types_.begin(); iter != allowed_res_types_.end(); ++iter ) {
				std::cerr << " " << *iter;
			}
			std::cerr << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}

		//now assign the atoms based on the cst file input, input of atom name takes priority  over input of atom type
		std::map< core::chemical::ResidueTypeCOP, utility::vector1< utility::vector1< core::Size > > >::iterator res_aid_it = atom_inds_for_restype_.find( cur_res );

		if ( res_aid_it == atom_inds_for_restype_.end() ) {
			this->determine_atom_inds_for_restype( cur_res );
			res_aid_it = atom_inds_for_restype_.find( cur_res );
		}

		for ( core::Size at1_ct = 1; at1_ct <= res_aid_it->second[1].size(); ++at1_ct ) {
			core::id::AtomID at1id( res_aid_it->second[1][at1_ct], respos_it->first);
			respos_it->second->atom1_.push_back( at1id );
			//tr << "for pose residue " << respos_it->first << " atom1, pushed back id " << res_aid_it->second[1][at1_ct] << std::endl;
		}

		for ( core::Size at2_ct = 1; at2_ct <= res_aid_it->second[2].size(); ++at2_ct ) {
			core::id::AtomID at2id( res_aid_it->second[2][at2_ct], respos_it->first);
			respos_it->second->atom2_.push_back( at2id );

			//tr << "for pose residue " << respos_it->first << " atom2, pushed back id " << res_aid_it->second[2][at2_ct] << std::endl;
		}

		for ( core::Size at3_ct = 1; at3_ct <= res_aid_it->second[3].size(); ++at3_ct ) {
			core::id::AtomID at3id( res_aid_it->second[3][at3_ct], respos_it->first);
			respos_it->second->atom3_.push_back( at3id );
			//tr << "for pose residue " << respos_it->first << " atom3, pushed back id " << res_aid_it->second[3][at3_ct] << std::endl;
		}

	} //loop over all residue positions in the pose

	template_cache->pose_data_uptodate_ = true;
} //get pose data function


EnzCstTemplateResAtomsCOP
EnzCstTemplateRes::get_template_atoms_at_pos( core::pose::Pose const & pose, core::Size seqpos ) const
{

	EnzCstTemplateResCacheCOP template_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( enz_io_param_.lock()->cst_block() )->template_res_cache( param_index_ ) );

	std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator at_it = template_cache->seqpos_map_.find( seqpos );

	if ( at_it == template_cache->seqpos_map_.end() ) {
		utility_exit_with_message("Error: could not find template atoms in EnzCstTemplateRes.\n");
	}

	return at_it->second;;
}

bool
EnzCstTemplateRes::find_in_pose_if_missing_from_header( core::pose::Pose & pose)
{
	EnzCstTemplateResCacheOP template_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( enz_io_param_.lock()->cst_block() )->template_res_cache( param_index_ ) );
	//there are three strategies to find a residue that wasn't mentioned in the PDB REMARKs
	//1. if it was given externally, i.e. through a 'seqpos' entry in the cst file or
	// through being speficially set through the accessor function
	//
	//2. if this residue is set to corresponding/identical to another residue in the cstfile
	//
	//3. if there is only a single residue in the pose that is in the allowed residue types
	// of this TemplateRes. This may sound somewhat arbitrary, but in most enzyme design cases
	// where there is only one ligand this will yield that ligand

	if ( template_cache->seqpos_map_.size() != 0 ) {
		utility_exit_with_message("Error: function find_in_pose_... was called even though there is stuff in the respos_map.\n");
	}

	//first option 1
	if ( respos_from_external_.size() != 0 ) {
		for ( utility::vector1< core::Size >::const_iterator ex_it = respos_from_external_.begin();
				ex_it != respos_from_external_.end(); ++ex_it ) {
			template_cache->add_position_in_pose( *ex_it );
		}
		template_cache->not_in_pose_ = false;
		//not_in_pose_ = false;
		return true;
	}

	//then option 2
	if ( (corresponding_res_block_ != 0 ) && ( corresponding_res_num_in_block_ != 0 ) ) {

		EnzCstTemplateResCOP corresponding_res;
		if ( corresponding_res_num_in_block_ == 1 ) {
			corresponding_res = enz_io_param_.lock()->enz_io().lock()->enz_cst_params( corresponding_res_block_ )->resA();
		} else {
			corresponding_res = enz_io_param_.lock()->enz_io().lock()->enz_cst_params( corresponding_res_block_ )->resB();
		}
		EnzCstTemplateResCacheCOP corresponding_res_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( corresponding_res->enz_io_param_.lock()->cst_block() )->template_res_cache( corresponding_res->param_index_ ) );

		//found the right template residue, now get the positions in the pose where it is
		for ( std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator pos_it = corresponding_res_cache->seqpos_map_begin();
				pos_it != corresponding_res_cache->seqpos_map_end(); ++pos_it ) {
			template_cache->add_position_in_pose( pos_it->first );
		}

		if ( template_cache->seqpos_map_.size() == 0 ) return false;
		else {
			template_cache->not_in_pose_ = false;
			return true;
		}
	}

	//then option 3
	utility::vector1< core::Size > found_positions;
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		std::string cur_res_name3 = pose.residue( i ).name3();
		utility::vector1< std::string >::iterator resfind = find( allowed_res_types_.begin(), allowed_res_types_.end(), cur_res_name3);

		if ( resfind != allowed_res_types_.end() ) found_positions.push_back( i );
	}
	if ( found_positions.size() == 1 ) {
		template_cache->add_position_in_pose( found_positions[1] );
		tr << "Found residue " << pose.residue( found_positions[1] ).name3() << " for CstBlock " << enz_io_param_.lock()->cst_block() << " without REMARK line in pose at position " << found_positions[1] << "." << std::endl;
		template_cache->not_in_pose_ = false;
		return true;
	}

	return false;

} //find_in_pose_if_missing_from_header


void
EnzCstTemplateRes::clear_all(){

	atom1_.clear(); atom2_.clear(); atom3_.clear();
	at1_type_ = ""; at2_type_ = ""; at3_type_ = "";
	allowed_res_types_.clear();
	allowed_res_types_pointers_.clear();
	corresponding_res_block_ = 0;
	corresponding_res_num_in_block_ = 0;
	respos_from_external_.clear();
	atom_inds_for_restype_.clear();
	param_index_ = 0;
}

void
EnzCstTemplateRes::set_external_position( core::Size resnum ){
	respos_from_external_.clear();
	respos_from_external_.push_back( resnum );
}


void
EnzCstTemplateRes::remap_resid( core::id::SequenceMapping const & smap )
{
	for ( core::Size i = 1; i <= respos_from_external_.size(); ++i ) {

		core::Size newpos = smap[ respos_from_external_[i] ];
		if ( newpos == 0 ) utility_exit_with_message("A catalytic residue is apparently missing from the pose");

		respos_from_external_[i] = newpos;
	}
}

bool
EnzCstTemplateRes::compatible_restype(
	core::chemical::ResidueTypeCOP restype
) const
{
	if ( !atom1_.empty() ) {
		for ( core::Size ii = 1; ii <= atom1_.size(); ++ii ) {
			if ( ! restype->has( atom1_[ii] ) ) {
				return false;
			}
		}
	}
	if ( !atom2_.empty() ) {
		for ( core::Size ii = 1; ii <= atom2_.size(); ++ii ) {
			if ( ! restype->has( atom2_[ii] ) ) {
				return false;
			}
		}
	}
	if ( !atom3_.empty() ) {
		for ( core::Size ii = 1; ii <= atom3_.size(); ++ii ) {
			if ( ! restype->has( atom3_[ii] ) ) {
				return false;
			}
		}
	}
	return true;
}


void
EnzCstTemplateRes::determine_atom_inds_for_restype(
	core::chemical::ResidueTypeCOP restype
) const
{

	Size natoms = restype->natoms();

	utility::vector1< core::Size > at1_ids, at2_ids, at3_ids;

	if ( !atom1_.empty() ) {
		for ( core::Size i = 1; i <= atom1_.size(); ++i ) {

			at1_ids.push_back( restype->atom_index(atom1_[i]) );

			if ( atom2_[i] != "" ) {
				at2_ids.push_back( restype->atom_index(atom2_[i]) );
			}

			if ( atom3_[i] != "" ) {
				at3_ids.push_back( restype->atom_index(atom3_[i]) );
			}
		}
	} else if ( at1_type_ != "" ) {

		// we need to know the atom indices based on the atom_type_name.
		// there is no function in the ResidueType class to do this yet, so we have to
		//do the following hacky implementation for now

		for ( Size ii = 1; ii <= natoms; ii++ ) {
			core::chemical::AtomType at1type = restype->atom_type( ii );

			if ( at1type.name() == at1_type_ ) {
				if ( basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] ) {
					tr.Info << "Adding atom " << ii << " with name " << at1type.name() << " for restype " << restype->name() << " to atom1 vector." << std::endl;
				}

				at1_ids.push_back( ii );
				//core::id::AtomID temp_at1id(ii,respos_it->first);
				//respos_it->second->atom1_.push_back( temp_at1id );
			}
		}

		if ( at1_ids.size() == 0 ) {
			std::cerr << "Error: ResidueType " << restype->name() <<  " does not have an atom of type " << at1_type_  << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}


		//now we have to find the base atoms for each possible atoms for atom1
		for ( utility::vector1< core::Size >::iterator at_it = at1_ids.begin(); at_it != at1_ids.end(); ++at_it ) {

			Size id_firstbase = restype->atom_base( *at_it );
			Size id_secondbase = restype->atom_base( id_firstbase );

			//ATTENTION: BACKBONE CONSTRAINT SPECIAL CASE: if an nbb or a CA atom are constrained
			//the secondbase will be identical to the first atom, i.e. this torsion will be undefined
			//to circumvent this, we set the id_secondbase to 3, the C atom
			if ( ( *at_it == 1) && ( id_firstbase == 2 ) && (id_secondbase == 1 ) ) id_secondbase = 3;
			else if ( ( *at_it == 2) && ( id_firstbase == 1) && (id_secondbase == 2 ) ) id_secondbase = 3;

			//ATTENTION: HISTIDINE SPECIAL CASE: if the first atom happens to be the ND1, we set the second base
			//to be the CD2. this ensures that torsion B is identical for ND1 and NE2, so only one cst file needed
			if ( restype->name3() == "HIS" && restype->atom_name( *at_it ) == " ND1" ) {
				id_secondbase = restype->atom_index( " CD2");
			}

			at2_ids.push_back( id_firstbase );
			at3_ids.push_back( id_secondbase );

			//core::id::AtomID temp_at2id( id_firstbase , respos_it->first );
			//core::id::AtomID temp_at3id( id_secondbase , respos_it->first  );
			//respos_it->second->atom2_.push_back( temp_at2id );
			//respos_it->second->atom3_.push_back( temp_at3id );
		}
	} else {
		std::cerr << "Error: cstfile did not specify any atoms for Residue " << restype->name3() << "." << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}


	utility::vector1< utility::vector1< core::Size > > at_ids;
	at_ids.push_back( at1_ids );
	at_ids.push_back( at2_ids );
	at_ids.push_back( at3_ids );

	atom_inds_for_restype_[ restype ] = at_ids;
	allowed_res_types_pointers_.push_back( restype );


} //determine_atom_ids_for_restype


utility::vector1< core::Size > const &
EnzCstTemplateRes::atom_inds_for_restype(
	core::Size template_atom,
	core::chemical::ResidueTypeCOP restype ) const
{

	//utility::vector1< std::string > const & atom_names;
	//std::string const & atom_type;


	// this might change in the future
	if ( template_atom == 1 ) {
		//atom_names = atom1_;
		//atom_type = at1_type_;
	} else if ( template_atom == 2 ) {
		//atom_names = atom2_;
		//atom_type = at2_type_;
	} else if ( template_atom == 3 ) {
		//atom_names = atom3_;
		//atom_type = at3_type_;
	} else utility_exit_with_message("When trying to find template_atom atom_ids for restype, illegal parameter was passed in for the template_atom. this has to be either 1,2, or 3" );

	RestypeToTemplateAtomsMap::const_iterator res_it = atom_inds_for_restype_.find( restype );
	if ( res_it == atom_inds_for_restype_.end() ) {
		this->determine_atom_inds_for_restype( restype );
		res_it = atom_inds_for_restype_.find( restype );
	}
	utility::vector1< utility::vector1< core::Size > > const & template_atoms = res_it->second;
	return template_atoms[ template_atom ];

}

bool
EnzCstTemplateRes::residue_conformations_redundant(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2
) const {
	core::Real const sqdist_cutoff(0.2*0.2); //hardcoded for now

	RestypeToTemplateAtomsMap::const_iterator res1_it(atom_inds_for_restype_.find( res1.type().get_self_ptr() ) ), res2_it(atom_inds_for_restype_.find( res2.type().get_self_ptr() ) );
	if ( res1_it == atom_inds_for_restype_.end() ) utility_exit_with_message("Residue of type "+res1.type().name()+" is not part of EnzCstTemplateRes.");
	if ( res2_it == atom_inds_for_restype_.end() ) utility_exit_with_message("Residue of type "+res2.type().name()+" is not part of EnzCstTemplateRes.");

	utility::vector1< utility::vector1< core::Size > > const & res1_at_ids( res1_it->second );
	utility::vector1< utility::vector1< core::Size > > const & res2_at_ids( res2_it->second );
	//tr << "beginning redundancy check, there are " << res1_at_ids[1].size() << " atom definitions for res1." << std::endl;
	//first redundancy check: if there are different numbers of
	//atom1s, the residues are non-redundant
	if ( res1_at_ids[1].size() != res2_at_ids[1].size() ) return false;

	for ( core::Size i = 1; i <= res1_at_ids[1].size(); ++i ) {
		//find the closest set from res2
		core::Real smallest_large_deviation_this_set( sqdist_cutoff + 1.0 );

		//tr << "ATOM1 res1 atom " << res1.atom_name( res1_at_ids[1][i] ) << " has sqdist of " << closest_atom1_sqdist << " from res2 atom " << res2.atom_name( res2_at_ids[1][1] ) << std::endl;

		for ( core::Size j = 1; j <= res2_at_ids[1].size(); ++j ) {
			core::Real large_deviation_this_pair(  res1.atom( res1_at_ids[1][i] ).xyz().distance_squared( res2.atom( res2_at_ids[1][j] ).xyz() ) );

			core::Real atom2_deviation(  res1.atom( res1_at_ids[2][i] ).xyz().distance_squared( res2.atom( res2_at_ids[2][j] ).xyz() ) );
			if ( atom2_deviation > large_deviation_this_pair ) large_deviation_this_pair = atom2_deviation;

			core::Real atom3_deviation(  res1.atom( res1_at_ids[3][i] ).xyz().distance_squared( res2.atom( res2_at_ids[3][j] ).xyz() ) );

			if ( atom3_deviation > large_deviation_this_pair ) large_deviation_this_pair = atom3_deviation;

			if ( large_deviation_this_pair < smallest_large_deviation_this_set ) smallest_large_deviation_this_set = large_deviation_this_pair;
			//tr << "res 1 set" <<i << " and res 2 set" << j << " have large_deviation of " << large_deviation_this_pair << std::endl;
		}
		if ( smallest_large_deviation_this_set > sqdist_cutoff ) return false;

	} //loop over the different template atom sets

	//if we make it to here, that means no atoms were
	//farther apart than the cutoff
	return true;
}

void
EnzCstTemplateRes::identical_info_consistency_check() const
{

	if ( !identical_tag_found_ ) return;

	//first some safety checks
	if ( (corresponding_res_block_ == 0) || ( ( corresponding_res_num_in_block_ != 1) && ( corresponding_res_num_in_block_ != 2) )
			|| ( enz_io_param_.lock()->cst_block() == corresponding_res_block_ ) ) {

		utility_exit_with_message("Cstfile has wrong format: 'identical' tag for Cstblock "+utility::to_string( enz_io_param_.lock()->cst_block() ) + " not formatted properly.\n");
	}

	if ( enz_io_param_.lock()->cst_block() < corresponding_res_block_ ) {
		utility_exit_with_message("Cstfile has wrong format: 'identical' tag for Cstblock "+utility::to_string( enz_io_param_.lock()->cst_block() ) + " is referring to a block of higher number. Please rewrite Cstfile such that identity tags only refer to blocks of higher numbers.\n");
	}

	EnzConstraintParametersCOP corresponding_pair = enz_io_param_.lock()->enz_io().lock()->enz_cst_params( corresponding_res_block_);

	utility::vector1< std::string > residue_names_to_match;

	if ( corresponding_res_num_in_block_ == 1 ) residue_names_to_match = corresponding_pair->resA()->allowed_res_types();
	else residue_names_to_match = corresponding_pair->resB()->allowed_res_types();

	for ( utility::vector1< std::string >::const_iterator restype_it = allowed_res_types_.begin();  restype_it != allowed_res_types_.end(); ++restype_it ) {
		utility::vector1< std::string >::const_iterator resfind = find(residue_names_to_match.begin(), residue_names_to_match.end(), *restype_it );

		if ( resfind == residue_names_to_match.end() ) {
			utility_exit_with_message("Error in cstfile 'identical' tag for Cstblock "+utility::to_string( enz_io_param_.lock()->cst_block() ) + " Allowed residue types for this cst block is not identical to the allowed residue types for the other cst_block.\n");
		}
	}
}

}
}//enzdes
}//protocols
