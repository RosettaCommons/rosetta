// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @author Eva-Maria Strauch (evas01@u.washington.edu)
// Unit headers
#include <protocols/seeded_abinitio/CAcstGenerator.hh>
#include <protocols/seeded_abinitio/CAcstGeneratorCreator.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/types.hh>
#include <core/pose/selection.hh>
#include <core/scoring/func/HarmonicFunc.hh>


#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <utility/vector0.hh>
#include <basic/options/keys/OptionKeys.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace seeded_abinitio {

using namespace core;
using namespace scoring::constraints;
using namespace protocols::moves;

static THREAD_LOCAL basic::Tracer TR( "protocols.seeded_abinitio.CAcstGenerator" );

// XRW TEMP std::string
// XRW TEMP CAcstGeneratorCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return CAcstGenerator::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP CAcstGeneratorCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new CAcstGenerator );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CAcstGenerator::mover_name()
// XRW TEMP {
// XRW TEMP  return "CAcstGenerator";
// XRW TEMP }

CAcstGenerator::~CAcstGenerator() = default;

CAcstGenerator::CAcstGenerator() :
	protocols::moves::Mover( CAcstGenerator::mover_name() )
{
	stddev_ = 3.0;
	seed_exceptions_.clear();
	distance_cutoff_ = 6;
}

protocols::moves::MoverOP
CAcstGenerator::clone() const {
	return( protocols::moves::MoverOP( new CAcstGenerator( *this ) ) );
}

protocols::moves::MoverOP
CAcstGenerator::fresh_instance() const {
	return protocols::moves::MoverOP( new CAcstGenerator );
}

bool
is_part( utility::vector1<core::Size> & cut_points,
	core::Size & residue){
	bool res_cut = false;
	for ( Size i = 1; i <= cut_points.size(); ++i ) {
		if ( cut_points[i] == residue ) {
			res_cut = true;
		}
	}
	return res_cut;
}

///this method will get de-convoluted soon....
void add_dist_constraints(
	pose::Pose & pose,
	pose::PoseOP & pose_of_int,
	core::Size start_relevant_chain,
	core::scoring::constraints::ConstraintSetOP & cst,
	protocols::loops::Loops & seeds,// referring to template if one is given! or just individual chain numbering
	protocols::loops::Loops & clear_area,
	utility::vector1< core::Size > cut_points,
	utility::vector1< core::Size > seed_exceptions,
	bool add_cst_seed,
	core::Real stddev,
	core::Size seq_separation,
	core::Real distance_cutoff
){

	using namespace scoring::constraints;
	using namespace id;
	//using namespace basic::options;
	//using namespace basic::options::OptionKeys;

	TR<<"stddev for harmonic constraints: " << stddev <<std::endl;

	if ( cut_points.empty() ) {
		TR<<"there are no cut point registered" <<std::endl;
	}

	TR.Debug << "start_relevant_chain: " << start_relevant_chain << std::endl;
	TR.Debug << "seeds: " << seeds << std::endl;
	TR.Debug << "clear_area " << clear_area << std::endl;

	for ( Size i=1; i<=cut_points.size(); i++ ) {
		TR <<"cutpoints: " << cut_points[i] << std::endl;
	}

	//adjust cutpoints first to relevant numbering
	for ( Size i = 1 ; i <= cut_points.size(); ++i ) {
		TR.Debug <<"rosetta numbering: cutpoint: "<< cut_points[i] <<std::endl;
		cut_points[i] = cut_points[i] - (start_relevant_chain -1);
		TR.Debug <<"adjusted cutpoint: "<< cut_points[i] <<std::endl;
	}

	for ( Size i=1; i<=cut_points.size(); i++ ) {
		TR <<"cutpoints: " << cut_points[i] << std::endl;
	}


	for ( Size pos = 1; pos <= pose_of_int->size(); pos++ ) {
		for ( Size pos_2 = 1; pos_2 <=pose_of_int->size(); pos_2++   ) {

			bool res_is_loop = false;
			bool res2_is_loop = false;

			if ( !is_part( seed_exceptions, pos) ) {
				//if residues are part of the loop, set to true
				res_is_loop = seeds.is_loop_residue( pos );
			}

			if ( !is_part( seed_exceptions, pos_2 ) ) {
				res2_is_loop = seeds.is_loop_residue( pos_2 );
			}

			bool cut_point_pos = false;
			bool cut_point_pos_2 = false;
			bool res_cst_free = false;

			//add the cst of the segment that will be replaced by the seed(s)
			if ( add_cst_seed ) {
				res_is_loop = false;
				res2_is_loop = false;
			}

			//user specified area that is allowed to float
			if ( clear_area.size() > 0 ) {
				if ( clear_area.is_loop_residue(pos) || clear_area.is_loop_residue(pos_2) ) {
					res_cst_free = true;
				}
			}

			//mark cut points
			if ( cut_points.size() > 1 ) {
				cut_point_pos = is_part(cut_points, pos );
				cut_point_pos_2 = is_part(cut_points, pos_2 );
			}

			//avoiding doubling of constraints
			Size seq_sep = 0;

			if ( pos > pos_2 ) {
				seq_sep = pos - pos_2;
			} else {
				seq_sep = 0; //to avoid repetition
			}
			//seq_sep = pos_2 - pos;

			if ( seq_sep >= seq_separation  ) {
				if ( !res_is_loop && !res2_is_loop ) {
					if ( !cut_point_pos && !cut_point_pos_2 ) {

						if ( !res_cst_free ) {
							core::conformation::Residue res_pos = pose_of_int->residue(pos);
							core::conformation::Residue res_pos_2 = pose_of_int->residue(pos_2);

							Real const distance_ca( res_pos.xyz( res_pos.atom_index("CA") ).distance( res_pos_2.xyz( res_pos_2.atom_index("CA") )));
							TR.Debug  <<"distance contraints for: " << pos_2 <<" "<< pos << " "<<distance_ca <<std::endl;
							TR.Debug <<"updated: "<<pos_2 + start_relevant_chain -1 <<" " << pos + start_relevant_chain - 1<< " " <<distance_ca << std::endl;

							if ( distance_ca > distance_cutoff ) {
								//adjust numbering to the current input pose!
								core::conformation::Residue res_in_pose = pose.residue( pos + start_relevant_chain - 1 );
								core::conformation::Residue res_in_pose2= pose.residue( pos_2 + start_relevant_chain - 1);

								cst->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint ( AtomID(res_in_pose.atom_index("CA"), pos + start_relevant_chain - 1), AtomID(res_in_pose2.atom_index("CA"),pos_2 + start_relevant_chain - 1), core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc( distance_ca, stddev ) ) ) ) ) );
							}
						}
					}
				}
			}//end seq_sep
		}
	}
}

void
CAcstGenerator::apply( pose::Pose & pose ){

	using namespace core;

	utility::vector1< Size > cutpoints = pose.fold_tree().cutpoints();

	for ( Size i=1; i<=cutpoints.size(); i++ ) {
		TR.Debug <<"cutpoints: " << cutpoints[i] << std::endl;
	}

	TR.Debug << "foldtree: " << pose.fold_tree();

	pose::PoseOP donor_poseOP;
	Size start_recipient_chain = 0;
	ca_cst_ = core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() );

	//this all has to be checked during run time since parse time input will be inaccurate if the pose has been grown before
	if ( from_chain_ == 0 ) {
		TR<<"user did not specify for which chain constraints should be derrived, defaulting to all chains"<<std::endl;
		if ( template_presence_ ) {
			donor_poseOP = template_pdb_;
			TR<<"derriving CA distance constraints from the template pdb, since a template was given"<<std::endl;
		} else {
			donor_poseOP = pose.get_self_ptr();
		}
	}

	//statements to prevent bogus from happening
	if ( to_chain_ == 0 && template_presence_ ) {
		if ( pose.size() != template_pdb_->size() ) {
			utility_exit_with_message("chain(s) to derrive constraints form and chain(s) to apply it to, do NOT have the same residue number. NOT supported." );
		}
	}

	if ( from_chain_ != 0 ) {
		if ( template_presence_ ) {
			donor_poseOP = template_pdb_->split_by_chain( from_chain_ );
			TR<<"derriving CA distance constraints from the template pdb"<<std::endl;
		} else {
			donor_poseOP = pose.split_by_chain( from_chain_ );
		}
	}

	if ( to_chain_ != 0 ) {
		start_recipient_chain = pose.conformation().chain_begin( to_chain_ );
		TR<<"adding constraints starting from index residue: " <<start_recipient_chain <<std::endl;
	}

	if ( from_chain_ != 0 && to_chain_ != 0 ) {
		TR<<"donor chain length: " << donor_poseOP->size() <<", recipient chain length: " << pose.conformation().chain_end( to_chain_ ) - pose.conformation().chain_begin( to_chain_) +1 <<std::endl;
		if ( donor_poseOP->size() != (pose.conformation().chain_end( to_chain_ ) - pose.conformation().chain_begin( to_chain_) + 1) ) {
			utility_exit_with_message("donor pose and recipient chain do not have the same residue numbers, check your template or input pdb");
		}
	}


	add_dist_constraints( pose, donor_poseOP , start_recipient_chain , ca_cst_, all_seeds_ , clear_seeds_ , cutpoints,
		seed_exceptions_, add_cst_seed_ , stddev_, seq_separation_, distance_cutoff_);

	if ( replace_ ) {
		TR<<"replacing all constraints with newly generated constraint set" <<std::endl;
		pose.constraint_set( ca_cst_ );
	} else {
		utility_exit_with_message("ADDing new constraints to pose, is currently not supported, try just replacing"); //<<std::endl;
	}
	TR.flush();
}

// XRW TEMP std::string
// XRW TEMP CAcstGenerator::get_name() const {
// XRW TEMP  return CAcstGenerator::mover_name();
// XRW TEMP }


void
CAcstGenerator::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	Movers_map const &,
	Pose const & ){

	TR<<"CAcstGenerator has been invoked"<<std::endl;

	/// the constraints can be either derrived from a chain of the input pose or from a template pose
	/// If a template pose is provided, the constraints will be derrived from that pdb, defaulting to chain 1 of it
	/// otherwise, the user can specify for which chain of the input pose constraints will be derrived
	/// currently it only allows one chain at a time or all chains....

	template_presence_ = false;

	stddev_ = tag->getOption<core::Real>( "stddev", 3.0);

	TR<<"setting constraint standard deviation to "<< stddev_<< std::endl;

	if ( tag->hasOption( "template_pdb" ) ) {
		std::string const template_pdb_fname( tag->getOption< std::string >( "template_pdb" ));
		template_pdb_ = core::pose::PoseOP( new core::pose::Pose ) ;
		core::import_pose::pose_from_file( *template_pdb_, template_pdb_fname , core::import_pose::PDB_file);
		TR<<"read in a template pdb with " <<template_pdb_->size() <<"residues"<<std::endl;
		template_presence_ = true;
	}

	add_cst_seed_ = tag->getOption< bool >("add_cst_seed", 0 ); ///header

	replace_ = tag->getOption< bool >("replace", 1 );

	seq_separation_ = tag->getOption< core::Size >( "seq_separation", 6 );

	distance_cutoff_ = tag->getOption< core::Real >("distance", 6.0 );

	if ( tag->hasOption( "add_seed_residues" ) ) {
		std::string residues_string = tag->getOption< std::string > ("add_seed_residues" );
		utility::vector1< std::string > const residue_keys( utility::string_split( residues_string, ',' ) );
		for ( std::string const & key : residue_keys ) {
			Size const res( utility::string2int( key ) );
			TR  << "add constraints to residues  within seed, residue: "<< key <<std::endl;
			seed_exceptions_.push_back( res );
		}
	}

	//parsing branch tags
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );

	for ( TagCOP const btag : branch_tags ) {
		//parse the pdb of interest, which is either the template or the input pdb depending on the users specificiation
		if ( template_presence_ ) {
			curr_pose_ = template_pdb_;
		}

		if ( btag->getName() == "Seeds" ) { //need an assertion for the presence of these or at least for the option file
			//needs some assertions to avoid bogus input
			std::string const beginS( btag->getOption<std::string>( "begin" ) );
			std::string const endS( btag->getOption<std::string>( "end" ) );
			core::Size const begin( core::pose::parse_resnum( beginS, *curr_pose_ ) );
			core::Size const end( core::pose::parse_resnum( endS, *curr_pose_ ) );

			TR.Debug <<"parsing seeds: \n"<< begin <<" and " << end <<std::endl;
			TR.Debug <<"seeds: "<< all_seeds_ <<std::endl;

			all_seeds_.add_loop( begin , end , 0, 0, false );

		}//end seed tags

		if ( btag->getName() == "Clear_cst_segment" ) { //need an assertion for the presence of these or at least for the option file

			std::string const begin_str( btag->getOption<std::string>( "begin" ) );
			std::string const end_str( btag->getOption<std::string>( "end" ) );
			core::Size const begin( core::pose::parse_resnum( begin_str, *curr_pose_ ) );
			core::Size const end( core::pose::parse_resnum( end_str, *curr_pose_ ) );
			clear_seeds_.add_loop( begin , end , 0, 0, false );

		}//end seed tags
	}//end branch tags

	///could be eventually a vector of chains if desired
	///but currently just the simple version...

	if ( tag->hasOption( "from_chain" ) ) {
		from_chain_ = tag->getOption< core::Size >( "from_chain", 1 );
		TR<<"chain to derrive constraints from: "<< from_chain_ <<std::endl;
	}

	if ( tag->hasOption( "to_chain" ) ) {
		to_chain_ = tag->getOption< core::Size >( "to_chain", 1 );
		TR<<"chain to apply constraints to: "<< to_chain_ <<std::endl;
	}

}//end parse my tag

std::string CAcstGenerator::get_name() const {
	return mover_name();
}

std::string CAcstGenerator::mover_name() {
	return "CAcstGenerator";
}

void CAcstGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	// Main attributes
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default("stddev", xsct_real, "Contraint strength as stddev.", "3.0")
		+ XMLSchemaAttribute("template_pdb", xs_string, "Template pdb to derive coordinate contraints from.")
		+ XMLSchemaAttribute::attribute_w_default("add_cst_seed", xsct_rosetta_bool, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default("replace", xsct_rosetta_bool, "XRW TO DO", "1")
		+ XMLSchemaAttribute::attribute_w_default("seq_separation", xsct_non_negative_integer, "XRW TO DO", "6")
		+ XMLSchemaAttribute::attribute_w_default("distance", xsct_real, "XRW TO DO", "6.0")
		+ XMLSchemaAttribute("add_seed_residues", xs_string, "Residues to which to apply contraints.")
		+ XMLSchemaAttribute::attribute_w_default("from_chain", xsct_non_negative_integer, "Chain from which to derive the constraints from.","1")
		+ XMLSchemaAttribute::attribute_w_default("to_chain", xsct_non_negative_integer, "Chain which to apply the contraints to.","1");

	// Subelements
	AttributeList subelement_attributes;
	subelement_attributes
		+ XMLSchemaAttribute::required_attribute("begin", xs_string, "First residue of a segment.")
		+ XMLSchemaAttribute::required_attribute("end", xs_string, "Last residue of a segment.");
	XMLSchemaSimpleSubelementList subelement_list;
	subelement_list.add_simple_subelement("Seeds", subelement_attributes, "Segment to which to apply the contraints?");
	subelement_list.add_simple_subelement("Clear_cst_segments", subelement_attributes, "Remove contraints from this segment?");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(),
		"Generate coordinate constraints from the input or a template pdb", attlist, subelement_list );
}

std::string CAcstGeneratorCreator::keyname() const {
	return CAcstGenerator::mover_name();
}

protocols::moves::MoverOP
CAcstGeneratorCreator::create_mover() const {
	return protocols::moves::MoverOP( new CAcstGenerator );
}

void CAcstGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CAcstGenerator::provide_xml_schema( xsd );
}


}//CAcstGenerator
}//protocol
