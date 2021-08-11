// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/movers/OmnibusDisulfideAnalysisLabelerMover.cc
/// @brief reports information about disulfides
/// @author Frank Teets (frankdt@email.unc.edu)

// Unit headers
#include <protocols/pose_sewing/movers/OmnibusDisulfideAnalysisLabelerMover.hh>
#include <protocols/pose_sewing/movers/OmnibusDisulfideAnalysisLabelerMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/io/Remarks.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/methods/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/util/disulfide_util.hh>

#include <protocols/pose_sewing/util.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
//#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <numeric/conversions.hh>
#include <basic/svd/SVD_Solver.hh>
#include <numeric/PCA.hh>
#include <protocols/sewing/scoring/MotifScorer.hh>

#include <protocols/rosetta_scripts/util.hh>

static basic::Tracer TR( "protocols.pose_sewing.movers.OmnibusDisulfideAnalysisLabelerMover" );

namespace protocols {
namespace pose_sewing {
namespace movers {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
OmnibusDisulfideAnalysisLabelerMover::OmnibusDisulfideAnalysisLabelerMover():
	protocols::moves::Mover( OmnibusDisulfideAnalysisLabelerMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
OmnibusDisulfideAnalysisLabelerMover::~OmnibusDisulfideAnalysisLabelerMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
OmnibusDisulfideAnalysisLabelerMover::apply( core::pose::Pose& pose){

	core::scoring::dssp::Dssp dssp(pose);
	dssp.insert_ss_into_pose(pose);
	core::pose::PDBInfoOP new_info = core::pose::PDBInfoOP(new core::pose::PDBInfo(pose,true));
	pose.pdb_info(new_info);
	core::io::Remarks pdb_remarks = pose.pdb_info()->remarks();
	core::io::RemarkInfo new_remark;
	//scorefxn_->score(pose);

	utility::vector1<std::pair<core::Size, core::Size>> respairs;
	std::pair<core::Size, core::Size> current_pair;
	//utility::vector1<core::scoring::ScoreType> score_types_to_use_only;
	//utility::vector1<core::scoring::ScoreType> score_types_to_ignore;
	//std::map< core::Size, core::Real> pair_energies;
	//core::Real relevant_score;

	std::string label;

	//pose.pdb_info()->add_reslabel(1,label);
	for ( core::Size upstream_res = 1; upstream_res <pose.size(); ++upstream_res ) {
		if ( pose.residue(upstream_res).type().is_disulfide_bonded() ) {
			for ( core::Size downstream_res = upstream_res+2; downstream_res <= pose.size(); ++downstream_res ) {
				if ( pose.residue(downstream_res).type().is_disulfide_bonded() && pose.residue(upstream_res).is_bonded(pose.residue(downstream_res)) ) {
					respairs.clear();
					current_pair.first = upstream_res;
					current_pair.second = downstream_res;
					respairs.push_back(current_pair);
					//TR << respairs << std::endl;
					// pair_energies = core::scoring::calculate_srlr_interaction_energy( pose, respairs, score_types_to_use_only, score_types_to_ignore);
					label = "DSSP_DISULFIDE_" + std::to_string(upstream_res) + "_" + std::to_string(downstream_res) + "_" + std::string(1,pose.secstruct(upstream_res)) + "_" + std::string(1,pose.secstruct(downstream_res));
					pose.pdb_info()->add_reslabel(upstream_res,label);
					new_remark.num = pdb_remarks.size();
					new_remark.value = label;
					pdb_remarks.push_back(new_remark);
					//TR << label << std::endl;
					// relevant_score = pair_energies[upstream_res];

					// pose.pdb_info()->add_reslabel(upstream_res,std::to_string(relevant_score));
				}
			}
		}
	}
	//
	// DisulfidizablePairs
	//

	std::map< core::Size, core::Size > element_blocks;
	protocols::pose_sewing::calculate_helices(element_blocks, pose, min_helix_length_); // key is residue, value is block

	std::set<std::pair<core::Size, core::Size>> helix_pairs;

	std::set<std::pair<core::Size, core::Size>> disulfidizable_helix_pairs;

	std::set<std::pair<core::Size, core::Size>> disulfide_bound_helix_pairs;

	std::pair<core::Size, core::Size> helix_pair;

	for ( core::Size upstream_res = 1; upstream_res < pose.size(); ++upstream_res ) {
		if ( pose.secstruct(upstream_res) == 'H' ) {
			for ( core::Size downstream_res = upstream_res+2; downstream_res <= pose.size(); ++downstream_res ) {
				if ( pose.secstruct(downstream_res) == 'H' ) {
					if ( element_blocks[upstream_res] != element_blocks[downstream_res] && pose.residue(upstream_res).xyz(2).distance(pose.residue(upstream_res).xyz(2)) <= crit_dist_ ) {
						helix_pair.first = element_blocks[upstream_res];
						helix_pair.second = element_blocks[downstream_res];
						helix_pairs.insert(helix_pair);
						if ( pose.residue(upstream_res).type().is_disulfide_bonded() && pose.residue(downstream_res).type().is_disulfide_bonded() && pose.residue(upstream_res).is_bonded(pose.residue(downstream_res)) ) {
							disulfide_bound_helix_pairs.insert(helix_pair);
						}
					}
				}

			}
		}
	}
	std::map<core::Size, std::map< core::Size, core::Real >> pair_energies;
	core::scoring::disulfides::DisulfideMatchingPotential disulfPot;
	core::Size pair_count = 0;
	for ( auto current_pair : helix_pairs ) {
		bool has_disulf = false;
		utility::vector1< core::Size > selection1;
		utility::vector1< core::Size > selection2;

		for ( core::Size current_res = 1; current_res < pose.size(); ++current_res ) {
			if ( element_blocks[current_res] == current_pair.first ) {
				selection1.push_back(current_res);
			}
			if ( element_blocks[current_res] == current_pair.second ) {
				selection2.push_back(current_res);
			}
		}

		for ( core::Size res1 : selection1 ) {
			std::map< core::Size, core::Real > energies;
			for ( core::Size res2 : selection2 ) {
				if ( res1 == res2 || has_disulf ) continue;

				//TR <<"Selection1/2 " << res1<<" " << res2 << std::endl;
				core::Energy match_t = 0.0;
				core::Energy match_r = 0.0;
				core::Energy match_rt = 0.0;
				disulfPot.score_disulfide( pose.residue(res1), pose.residue(res2), match_t, match_r, match_rt, false /* mirror */);
				if ( match_rt <= cutoff_ ) {
					has_disulf = true;
					disulfidizable_helix_pairs.insert(current_pair);
				}
			}
		}
		if ( has_disulf ) {
			++pair_count;
		}
	}
	new_remark.num = pdb_remarks.size();
	new_remark.value = "DISULFIDIZABLE_PAIRS: " + utility::to_string(pair_count);
	pdb_remarks.push_back(new_remark);
	//
	//  HelixPairStats
	//
	protocols::sewing::scoring::MotifScorer motif_scorer = protocols::sewing::scoring::MotifScorer();

	//first, find all pairs
	//std::set<std::pair<core::Size, core::Size>> helix_pairs; // need to split into disulfidizable_helix_pairs and all_helix_pairs
	//helix_pairs.clear();

	//std::pair<core::Size, core::Size> helix_pair;

	//for(core::Size upstream_res = 1; upstream_res < pose.size();++upstream_res){
	//        if(pose.secstruct(upstream_res) == 'H' && (!disulf_ || pose.residue(upstream_res).type().is_disulfide_bonded())){
	//                for(core::Size downstream_res = upstream_res+2; downstream_res <= pose.size();++downstream_res){
	//                        if(pose.secstruct(downstream_res) == 'H' && (!disulf_ || (pose.residue(downstream_res).type().is_disulfide_bonded() && pose.residue(upstream_res).is_bonded(pose.residue(downstream_res))))){
	//                                if(element_blocks[upstream_res] != element_blocks[downstream_res] && pose.residue(upstream_res).xyz(2).distance(pose.residue(upstream_res).xyz(2)) <= crit_dist_ ){
	//                                        helix_pair.first = element_blocks[upstream_res];
	//                                        helix_pair.second = element_blocks[downstream_res];
	//                                        helix_pairs.insert(helix_pair);
	//                                }
	//                        }

	//                }
	//        }
	//}
	// end pair find
	//
	// now, prune helix pairs that are insufficiently large
	std::map<std::pair<core::Size,core::Size>,core::Real> outmap;
	// now, prune helix pairs that are insufficiently packable
	core::Real running_score;
	for ( auto current_pair : helix_pairs ) {
		core::Size pair_count = 0;
		for ( core::Size N_resnum = 1; N_resnum <= pose.size(); N_resnum++ ) {
			if ( !pose.residue(N_resnum).is_virtual_residue() && element_blocks[N_resnum] == current_pair.first ) {
				char ss1 = pose.secstruct(N_resnum);
				core::conformation::Residue const & N_res = pose.residue(N_resnum);

				for ( core::Size C_resnum = 1; C_resnum <= pose.size(); C_resnum++ ) {
					if ( !pose.residue(C_resnum).is_virtual_residue() && element_blocks[C_resnum] == current_pair.second ) {

						char ss2 = pose.secstruct(C_resnum);
						core::conformation::Residue const & C_res = pose.residue(C_resnum);

						if ( use_motifs_ ) {
							running_score = -1 * protocols::pose_sewing::calculate_motif_score_bt_residues(motif_scorer, N_res, C_res, ss1, ss2);
						} else {
							running_score = -1 * protocols::pose_sewing::calculate_distance_score_bt_residues(N_res,C_res,-2.0, 4.0, 0.1);
						}

						//running_score = -1 * motif_scorer->get_score(stub1, ss1, aa1, stub2, ss2, aa2);

						if ( running_score < crit_score_ ) {
							++ pair_count;
							if ( outmap.count(current_pair) == 0 ) {
								outmap[current_pair] =  running_score;
							} else {
								outmap[current_pair] = outmap[current_pair] + running_score;
							}
						}
					}
				}
			}
		}
		if ( pair_count > 0 ) {
			outmap[current_pair] = outmap[current_pair] / pair_count;
		} else {
			outmap[current_pair] = 0;
		}

	}

	// end prune
	//std::map<std::pair<core::Size,core::Size>,core::Real> outmap;
	for ( auto current_pair : helix_pairs ) {
		if ( outmap[current_pair] > qualifying_score_ ) {
			continue;
		}
		//below taken from HelixHelixAngleFilter
		utility::vector1< numeric::xyzVector <core::Real > > helix_vector_1;
		utility::vector1< numeric::xyzVector <core::Real > > helix_vector_2;

		utility::vector1< numeric::xyzVector< core::Real> > bb_coords_1;
		bb_coords_1.clear();
		utility::vector1< numeric::xyzVector< core::Real> > bb_coords_2;
		bb_coords_2.clear();

		for ( core::Size resnum = 1; resnum <= pose.size(); resnum++ ) {
			core::conformation::Residue rsd( pose.residue( resnum ) );
			if ( element_blocks[resnum] == current_pair.first ) {
				for ( core::Size j=1; j <= rsd.nheavyatoms(); ++j ) {
					if ( rsd.atom_is_backbone( j ) ) {
						core::Vector const & bbvec( rsd.xyz( j ) );
						bb_coords_1.push_back( bbvec );
					}
				}
			}
			if ( element_blocks[resnum] == current_pair.second ) {
				for ( core::Size j=1; j <= rsd.nheavyatoms(); ++j ) {
					if ( rsd.atom_is_backbone( j ) ) {
						core::Vector const & bbvec( rsd.xyz( j ) );
						bb_coords_2.push_back( bbvec );
					}
				}
			}


		}
		numeric::xyzVector<core::Real> com_1 = numeric::center_of_mass( bb_coords_1 );
		numeric::xyzVector<core::Real> com_2 = numeric::center_of_mass( bb_coords_2 );

		numeric::xyzVector<core::Real> first_principal_component_1 = numeric::first_principal_component( bb_coords_1 );
		numeric::xyzVector<core::Real> first_principal_component_2 = numeric::first_principal_component( bb_coords_2 );

		numeric::xyzVector<core::Real> com_principal_component_1 = first_principal_component_1 += com_1;
		numeric::xyzVector<core::Real> com_principal_component_2 = first_principal_component_2 += com_2;

		//p0 is point on principal component vector closest to N-term
		numeric::xyzVector<core::Real> p0_1 = numeric::closest_point_on_line( com_1, com_principal_component_1, bb_coords_1[1]);
		numeric::xyzVector<core::Real> p0_2 = numeric::closest_point_on_line( com_2, com_principal_component_2, bb_coords_2[1]);

		//p1 is point on principal component vector closest to C-term
		numeric::xyzVector<core::Real> p1_1 = numeric::closest_point_on_line( com_1, com_principal_component_1, bb_coords_1.back());
		numeric::xyzVector<core::Real> p1_2 = numeric::closest_point_on_line( com_2, com_principal_component_2, bb_coords_2.back());

		helix_vector_1.push_back( p0_1 );
		helix_vector_1.push_back( p1_1 );

		helix_vector_2.push_back( p0_2 );
		helix_vector_2.push_back( p1_2 );

		std::pair< numeric::xyzVector<core::Real>, numeric::xyzVector<core::Real>> closest_points;

		numeric::xyzVector< core::Real > u( helix_vector_1[2] - helix_vector_1[1] );
		numeric::xyzVector< core::Real > v( helix_vector_2[2] - helix_vector_2[1] );
		numeric::xyzVector< core::Real > w( helix_vector_1[1] - helix_vector_2[1] );

		core::Real a( u.dot( u ) );
		core::Real b( u.dot( v ) );
		core::Real c( v.dot( v ) );
		core::Real d( u.dot( w ) );
		core::Real e( v.dot( w ) );
		core::Real D( ( a * c ) - ( b * b ) );

		core::Real sc( 0.0 );
		core::Real tc( 0.0 );

		if ( D < 0.000001 ) {
			sc = 0.0;
			tc = ( b > c ? d / b : e / c );
		} else {
			sc = ( ( b * e ) - ( c * d ) ) / D;
			tc = ( ( a * e ) - ( b * d ) ) / D;
		}

		numeric::xyzVector< core::Real > v1( helix_vector_1[1] + sc * ( helix_vector_1[2] - helix_vector_1[1] ) );
		numeric::xyzVector< core::Real > v2( helix_vector_2[1] + tc * ( helix_vector_2[2] - helix_vector_2[1] ) );
		//TR << "v1 " << v1.x() << " " << v1.y() << " " << v1.z() << std::endl;
		//TR << "v2 " << v2.x() << " " << v2.y() << " " << v2.z() << std::endl;

		closest_points.first = v1;
		closest_points.second = v2;

		core::Real cross_angle = numeric::dihedral_degrees( helix_vector_1[1], closest_points.first, closest_points.second, helix_vector_2[1] );

		new_remark.num = pdb_remarks.size();
		new_remark.value = "ALL_HELIX-HELIX_ANGLE: " + utility::to_string(cross_angle);
		pdb_remarks.push_back(new_remark);

		if ( disulfidizable_helix_pairs.find(current_pair) != disulfidizable_helix_pairs.end() ) {
			new_remark.num = pdb_remarks.size();
			new_remark.value = "DISULFIDIZABLE_HELIX-HELIX_ANGLE: " + utility::to_string(cross_angle);
			pdb_remarks.push_back(new_remark);
		}
		if ( disulfide_bound_helix_pairs.find(current_pair) != disulfide_bound_helix_pairs.end() ) {
			new_remark.num = pdb_remarks.size();
			new_remark.value = "DISULFIDE_BOUND_HELIX-HELIX_ANGLE: " + utility::to_string(cross_angle);
			pdb_remarks.push_back(new_remark);
		}

		//out[std::to_string(current_pair.first)+":"+std::to_string(current_pair.second)] =  cross_angle;

	}

	//for (auto outpair : out){
	//        //TR << "OUTPAIR " << outpair.second << std::endl;
	//}
	//
	pose.pdb_info()->remarks(pdb_remarks);

}


////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
OmnibusDisulfideAnalysisLabelerMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
OmnibusDisulfideAnalysisLabelerMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap
) {

	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );
	}

}
void OmnibusDisulfideAnalysisLabelerMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
OmnibusDisulfideAnalysisLabelerMover::fresh_instance() const
{
	return utility::pointer::make_shared< OmnibusDisulfideAnalysisLabelerMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
OmnibusDisulfideAnalysisLabelerMover::clone() const
{
	return utility::pointer::make_shared< OmnibusDisulfideAnalysisLabelerMover >( *this );
}

std::string OmnibusDisulfideAnalysisLabelerMover::get_name() const {
	return mover_name();
}

std::string OmnibusDisulfideAnalysisLabelerMover::mover_name() {
	return "OmnibusDisulfideAnalysisLabelerMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
OmnibusDisulfideAnalysisLabelerMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< OmnibusDisulfideAnalysisLabelerMover >();
}

std::string
OmnibusDisulfideAnalysisLabelerMoverCreator::keyname() const
{
	return OmnibusDisulfideAnalysisLabelerMover::mover_name();
}

void OmnibusDisulfideAnalysisLabelerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	OmnibusDisulfideAnalysisLabelerMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns Frank Teets as its author.
void
OmnibusDisulfideAnalysisLabelerMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"OmnibusDisulfideAnalysisLabelerMover", basic::citation_manager::CitedModuleType::Mover,
		"Frank Teets",
		"TODO: institution",
		"frankdt@email.unc.edu",
		"Wrote the OmnibusDisulfideAnalysisLabelerMover."
		)
	);
}


////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, OmnibusDisulfideAnalysisLabelerMover const & mover )
{
	mover.show(os);
	return os;
}


} //movers
} //pose_sewing
} //protocols
