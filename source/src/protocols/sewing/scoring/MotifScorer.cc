// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MotifScorer.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/scoring/MotifScorer.hh>
#include <protocols/sewing/scoring/MotifScorerCreator.hh>
//Project headers
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <protocols/sewing/scoring/AssemblyScorerFactory.hh>
//Core headers
#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/motif/reference_frames.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/motif/util.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzTransform.hh>
#include <utility/io/ozstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
namespace protocols {
namespace sewing  {
namespace scoring {

static basic::Tracer TR("protocols.sewing.scoring.MotifScorer");

MotifScorer::MotifScorer():
	AssemblyScorer(),
	mman_(*core::scoring::motif::MotifHashManager::get_instance()),
	res_type_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) )
{
	weight_ = 1.0;
}

MotifScorer::MotifScorer( MotifScorer const & src ):
	AssemblyScorer( src ),
	mman_( src.mman_ ),
	res_type_set_( src.res_type_set_ ),
	weight_( src.weight_ )
{
}
///@details Use the negative normalized motif score
core::Real
MotifScorer::score(
	data_storage::SmartAssemblyCOP assembly
) {
	this->set_old_last_score( this->get_last_score() );
	core::Real nscore = this->norm_motif_score( assembly );
	this->set_last_score(  -1.0 * nscore );
	TR.Debug << "Final motif score: " << this->get_last_score() << std::endl;
	return this->get_last_score();
}


///@details Motif score of the entire Assembly divided by total residue
core::Real
MotifScorer::norm_motif_score(
	data_storage::SmartAssemblyCOP assembly
) {
	core::Real fscore = this->full_motif_score( assembly );
	return fscore / core::Real (assembly->get_length() );
}


///@details use Will's Motif score to calculate the motif score for the entire Assembly
core::Real
MotifScorer::full_motif_score(
	data_storage::SmartAssemblyCOP assembly
) {
	core::Real score = 0.0;
	data_storage::SmartSegmentOP current_segment = assembly->get_n_terminal_segment();
	runtime_assert( current_segment != nullptr );
	//Note that this segment is pulled from the assembly's set of segments, not the master set
	while ( current_segment != nullptr ) {
		//TR << "Scoring segment " << current_segment->get_segment_id() << " of length " << current_segment->get_length() << std::endl;
		//DSSP actually lives in segment and won't change, so just get it here
		char ss1 = current_segment->get_dssp_code();
		//Iterate through the residues in this segment
		//I need the indices to avoid double counting
		//  for( data_storage::SmartSewingResidueOP current_res : current_segment->get_residue_vector() ){
		for ( core::Size current_res_index = 1; current_res_index <=current_segment->get_length() ; ++current_res_index ) {
			//New version of get_stub just takes the SmartSewingResidueOP
			data_storage::SmartSewingResidueCOP current_res = current_segment->get_const_residue_vector().at( current_res_index );
			numeric::xyzTransform< core::Real > stub1 = get_stub( current_res );
			//We already have the DSSP, just get the aa type
			char aa1 = res_type_set_->name_map( current_res->get_amino_acid_type() ).name1();

			//Now loop through all of the segments starting at this segment
			data_storage::SmartSegmentOP second_segment = current_segment;
			while ( second_segment != nullptr ) {
				//TR << "Scoring against segment " << second_segment->get_segment_id() << " of length " << second_segment->get_length() << std::endl;
				char ss2 = second_segment->get_dssp_code();
				//for( data_storage::SmartSewingResidueOP second_res : second_segment->get_residue_vector() ){
				for ( core::Size second_res_index = 1; second_res_index <= second_segment->get_length(); ++second_res_index ) {
					//Double check that we're not comparing the residue to itself
					if ( current_segment->get_segment_id() == second_segment->get_segment_id() && second_res_index == current_res_index ) { //These two pointers actually point to the same residue
						continue;
					}
					data_storage::SmartSewingResidueCOP second_res = second_segment->get_const_residue_vector().at( second_res_index );
					numeric::xyzTransform< core::Real > stub2 = get_stub( second_res );
					char aa2 = res_type_set_->name_map( second_res->get_amino_acid_type() ).name1();
					//TR << "Scoring atom pair with ss1 " << ss1 << " ss2 " << ss2 << " aa1 " << aa1 << " aa2 " << aa2 << std::endl;
					score += get_score(stub1, ss1, aa1, stub2, ss2, aa2);
				} //end for second_res
				second_segment = second_segment->get_c_terminal_neighbor();
			} //end while second_segment
		} //end for current_res
		//if( !( current_segment->is_c_terminus_fixed() ) ){
		// break;
		//}
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	return score;
}


//    core::pose::motif::get_backbone_reference_frame( atom1.xyz(), atom2.xyz(), atom3.xyz() )
//This will change to just take the SmartSewingResidueOP (that's all it needs)
numeric::xyzTransform<core::Real>
MotifScorer::get_stub(
	data_storage::SmartSewingResidueCOP residue) const
{
	//It doesn't like trying to do this all on one line
	utility::vector1< core::conformation::Atom > const & atoms = residue->get_const_atom_vector();
	return core::pose::motif::get_backbone_reference_frame(
		atoms[ 1 ].xyz(), atoms[ 2 ].xyz(), atoms[ 3 ].xyz() );
}


//This really shouldn't need to change at all!
core::Real
MotifScorer::get_score(
	numeric::xyzTransform<core::Real> stub1,
	char ss1,
	char aa1,
	numeric::xyzTransform<core::Real> stub2,
	char ss2,
	char aa2
) const {
	core::Real score = 0;
	core::scoring::motif::XformScoreCOP xs_bb_fxn1 = mman_.get_xform_score_BB_BB(ss1,ss2,aa1,aa2);
	if ( !xs_bb_fxn1 ) {
		std::stringstream err;
		err << "Null XformScore pointer!" << std::endl;
		err << "ss1 " << ss1 << " aa1" << aa1 << std::endl;
		err << "ss2 " << ss2 << " aa2" << aa2 << std::endl;
		utility_exit_with_message(err.str());
	}
	numeric::xyzTransform<core::Real> const Xbb = stub1.inverse() * stub2;
	if ( Xbb.x() < 16 && Xbb.x() > -16
			&& Xbb.y() < 16 && Xbb.y() > -16
			&& Xbb.z() < 16 && Xbb.z() > -16
			) {
		score += xs_bb_fxn1->score_of_bin(Xbb);
	}

	core::scoring::motif::XformScoreCOP xs_bb_fxn2 = mman_.get_xform_score_BB_BB(ss2,ss1,aa2,aa1);
	if ( !xs_bb_fxn2 ) {
		std::stringstream err;
		err << "Null XformScore pointer!" << std::endl;
		err << "ss1 " << ss1 << " aa1" << aa1 << std::endl;
		err << "ss2 " << ss2 << " aa2" << aa2 << std::endl;
		utility_exit_with_message(err.str());
	}
	numeric::xyzTransform<core::Real> const Xbbi = Xbb.inverse();
	if ( Xbbi.x() < 16 && Xbbi.x() > -16
			&& Xbbi.y() < 16 && Xbbi.y() > -16
			&& Xbbi.z() < 16 && Xbbi.z() > -16
			) {
		score += xs_bb_fxn2->score_of_bin(Xbbi);
	}
	return score;
}


void
MotifScorer::set_options_from_tag(
	utility::tag::TagCOP scorer,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const &)
{
	//No data members to set!

	if ( !scorer->hasOption( "weight" ) ) {
		TR << "No weight provided for scorer " << scorer->getName() << " Using default value" << std::endl;
	}

	weight_ = scorer->getOption<core::Real>( "weight", 1.0 );

	TR.Debug << "Created MotifScorer from tag!" << std::endl;

}


std::string
MotifScorer::get_name() const {
	return "MotifScorer";
}

core::Real
MotifScorer::get_weight() const{
	return weight_;
}

void
MotifScorer::set_weight( core::Real weight ){
	weight_ = weight;
}

void
MotifScorer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	//only attribute is the weight

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "How heavily will this term be weighted during scoring?", "1.0"  );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyScorerFactory::assembly_scorer_ct_namer )
		.element_name( MotifScorer::type_name() )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.description( "Basic Motif score among all helices" )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations



}

std::string
MotifScorer::type_name(){
	return "MotifScorer";
}



//////Creator methods//////////
AssemblyScorerOP
MotifScorerCreator::create_assembly_scorer() const{
	return AssemblyScorerOP( new MotifScorer );
}



std::string
MotifScorerCreator::keyname() const{
	return MotifScorer::type_name();
}

void
MotifScorerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	MotifScorer::provide_xml_schema( xsd );


}

core::Real
MotifScorer::get_last_score() const{
	return last_score_;
}

void
MotifScorer::set_last_score( core::Real score ) {
	last_score_ = score;
}

core::Real
MotifScorer::get_old_last_score() const{
	return old_last_score_;
}

void
MotifScorer::set_old_last_score( core::Real score ){
	old_last_score_ = score;
}

/////End Creator methods///////



//Not really sure when we even use this. Why would we convert to a pose here?
/*
void
MotifScorer::dump_motif(
AssemblyCOP assembly
) const {
core::pose::Pose allpose(assembly->to_pose(core::chemical::FA_STANDARD, false));
core::pose::Pose assembledpose = *(allpose.split_by_chain(1));
core::scoring::dssp::Dssp dssp(assembledpose);
dssp.insert_ss_into_pose(assembledpose);
core::scoring::motif::ResPairMotifQuery opt(assembledpose);
opt.interface_only() = false;
core::scoring::motif::MotifHits hits;
core::scoring::motif::MotifHashManager::get_instance()->get_matching_motifs(opt,hits);
if ( hits.size() ) {
std::string outfile = "motif_test.pdb";
std::cout << "dump " << hits.size() << " (before pruning) hits to " << outfile << std::endl;
utility::io::ozstream pdbout( outfile );
pdbout << "MODEL MAIN" << std::endl;
assembledpose.dump_pdb(pdbout);
pdbout << "ENDMDL" << std::endl;
hits.dump_motifs_pdb(pdbout);
pdbout.close();
}
}
*/

} //scoring namespace
} //sewing namespace
} //protocols namespace
