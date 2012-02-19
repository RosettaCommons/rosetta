// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file BundlePairRmsdWorkUnit.cc
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

//Unit
#include <devel/helixAssembly/BundlePairRmsdWorkUnit.hh>

//Basic
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>

//Utility
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

//C++
#include <string>
#include <map>

static basic::Tracer TR("BundlePairRmsdWorkUnit");

using namespace std;

BundlePairRmsdWorkUnit::BundlePairRmsdWorkUnit(utility::sql_database::sessionOP db_session):
DatabaseEntryWorkUnit(db_session)
{}

BundlePairRmsdWorkUnit::BundlePairRmsdWorkUnit( std::map<std::string,std::string> row_map ):
DatabaseEntryWorkUnit(row_map)
{
    TR << "Constructing bundle pair RMSD work unit!" << endl;
}

/// @brief Calculate the pair RMSD using data from the results_map_
void
BundlePairRmsdWorkUnit::run(){
    using namespace cppdb;
    using namespace utility;
    
    TR << "BundlePairRMSD work unit recieved following key/value pairs:" << endl;
        
    for( map<string,string>::const_iterator it = row_map_.begin();
        it != row_map_.end(); ++it){
        TR << it->first << ", " << it->second << endl;
    }
    
    core::Size struct_id, bundle_id;
    struct_id=string2int(row_map_["struct_id"]);
    bundle_id=string2int(row_map_["bundle_id"]);
    
    std::string select_string = 
    "SELECT * FROM bundle_helices\n"
    "WHERE bundle_id = ?";
        
    statement query_statement(basic::database::safely_prepare_statement(select_string,db_session_));
    
    query_statement.bind(1,bundle_id);
    
    result res(basic::database::safely_read_from_database(query_statement));
    
    int test_counter=0;
    
    while(res.next()){
        ++test_counter;
    }
    
    //create the result query - this will be sent back to master and executed on the DB
    result_query_string_ = "Test query returned " + utility::to_string(test_counter) + " rows";
}

//    THINK ABOUT CREATING  A TEMPORARY TABLE FROM BELOW STATEMENT AND DISTRIBUTING A SELECTION FROM THAT

//    std::string select_string =
//    "SELECT\n"
//    "   bundles.struct_id,\n"
//    "	bundles.bundle_id,\n"
//    "   helices1.helix_id,\n"
//    "	helices1.residue_begin AS helix_1_begin,\n"
//    "	helices1.residue_end AS helix_1_end,\n"
//    "   helices2.helix_id,\n"
//    "	helices2.residue_begin AS helix_2_begin,\n"
//    "	helices2.residue_end AS helix_2_end\n"
//    "FROM\n"
//    "	helix_bundles AS bundles\n"
//    "JOIN bundle_helices AS helices1 ON\n"
//    "   bundles.bundle_id = helices1.bundle_id\n"
//    "JOIN bundle_helices AS helices2 ON\n"
//    "   bundles.bundle_id = helices2.bundle_id\n"
//    "WHERE\n"
//    "   helices1.bundle_id = helices2.bundle_id AND\n"
//    "   helices1.helix_id <> helices2.helix_id";
//    
//	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session_));
//	result res1(basic::database::safely_read_from_database(select_statement));
//    result res2(basic::database::safely_read_from_database(select_statement));
//    
//    while(res1.next()){
//        Size pair1_struct_id, pair1_bundle_id, pair1_helix1_id, pair1_helix1_begin, pair1_helix1_end, 
//        pair1_helix2_id, pair1_helix2_begin, pair1_helix2_end;
//        
//        res1 >> pair1_struct_id >> pair1_bundle_id >> pair1_helix1_id >> pair1_helix1_begin >> pair1_helix1_end >>
//        pair1_helix2_id >> pair1_helix2_begin >> pair1_helix2_end;
//        
//        //Create a DB Work Unit out of res1 and res2?
//        
//        std::pair<HelicalFragment, HelicalFragment> pair1(HelicalFragment(pair1_helix1_begin, pair1_helix1_end), 
//                                                          HelicalFragment(pair1_helix2_begin, pair1_helix2_end));
//        
//        //This is needlessly expensive since we'll be recreating poses even when struct_id doesn't change
//        core::pose::Pose pose1;
//        pose_conformation_features_->load_into_pose(db_session_, pair1_struct_id, pose1);
//        protein_residue_conformation_features_->load_into_pose(db_session_, pair1_struct_id, pose1);
//        
//        utility::vector1< std::pair<HelicalFragment,HelicalFragment> > bundle_pairs;
//        while(res2.next()){
//            
//            Size pair2_struct_id, pair2_bundle_id, pair2_helix1_id, pair2_helix1_begin, pair2_helix1_end, 
//            pair2_helix2_id, pair2_helix2_begin, pair2_helix2_end;
//            
//            res2 >> pair2_struct_id >> pair2_bundle_id >> pair2_helix1_id >> pair2_helix1_begin >> pair2_helix1_end >>
//            pair2_helix2_id >> pair2_helix2_begin >> pair2_helix2_end;
//            
//            std::pair<HelicalFragment, HelicalFragment> pair2(HelicalFragment(pair2_helix1_begin, pair2_helix1_end), HelicalFragment(pair2_helix2_begin, pair2_helix2_end));
//            
//            
//            //if these pairs are part of the same bundle, don't add their RMSD to the table
//            if(pair1_struct_id != pair2_struct_id && pair1_bundle_id != pair2_bundle_id){
//                
//                //This is needlessly expensive since we'll be recreating poses even when struct_id doesn't change
//                core::pose::Pose pose2;
//                pose_conformation_features_->load_into_pose(db_session_, pair2_struct_id, pose2);
//                protein_residue_conformation_features_->load_into_pose(db_session_, pair2_struct_id, pose2);
//                
//                Real rmsd = BundlePairFeatures::calculate_pair_rmsd(pose1, pair1, pose2, pair2);
//                
//                std::string insert_string = "INSERT INTO bundle_pair_comparisons VALUES (?,?,?,?,?,?,?,?,?)";
//                statement insert_statement(basic::database::safely_prepare_statement(insert_string,db_session_));
//                insert_statement.bind(1,pair1_struct_id);
//                insert_statement.bind(2,pair1_bundle_id);
//                insert_statement.bind(3,pair1_helix1_id);
//                insert_statement.bind(4,pair1_helix2_id);
//                
//                insert_statement.bind(5,pair2_struct_id);
//                insert_statement.bind(6,pair2_bundle_id);
//                insert_statement.bind(7,pair2_helix1_id);
//                insert_statement.bind(8,pair2_helix2_id);
//                
//                insert_statement.bind(9,rmsd);
//                basic::database::safely_write_to_database(insert_statement);
//            }
//        }
//    }




//namespace protocols {
//    namespace features {
//        namespace helixAssembly {
//            
//            using namespace std;
//            using namespace core;
//            using core::pose::Pose;
//            using utility::vector1;
//            using utility::sql_database::sessionOP;
//            using cppdb::statement;
//            using cppdb::result;
//            
//            BundlePairFeatures::BundlePairFeatures():
//            pose_conformation_features_( new protocols::features::PoseConformationFeatures() ),
//            protein_residue_conformation_features_( new protocols::features::ProteinResidueConformationFeatures() )
//            {}
//            
//            string
//            BundlePairFeatures::schema() const {
//                return
//                
//                "CREATE TABLE IF NOT EXISTS bundle_pair_comparisons (\n"
//                "   struct_id_1 INTEGER,\n"
//                "   bundle_id_1 INTEGER,\n"
//                "   bundle_1_helix_1 INTEGER,\n"
//                "   bundle_1_helix_2 INTEGER,\n"
//                "   struct_id_2 INTEGER,\n"
//                "   bundle_id_2 INTEGER,\n"
//                "   bundle_2_helix_1 INTEGER,\n"
//                "   bundle_2_helix_2 INTEGER,\n"
//                "   bb_rmsd INTEGER,\n"
//                "	FOREIGN KEY(struct_id_1, bundle_id_1, struct_id_2, bundle_id_2)\n"
//                "		REFERENCES helix_bundles(struct_id, bundle_id, struct_id, bundle_id)\n"
//                "		DEFERRABLE INITIALLY DEFERRED,\n"
//                "	FOREIGN KEY(bundle_1_helix_1, bundle_1_helix_2, bundle_2_helix_1, bundle_2_helix_2)\n"
//                "		REFERENCES helix_bundles(helix_id,helix_id,helix_id,helix_id)\n"
//                "		DEFERRABLE INITIALLY DEFERRED,\n"
//                "	PRIMARY KEY(struct_id_1, bundle_id_1, struct_id_2, bundle_id_2));"
//                ;
//            }
//            
//            //@brief calculate the minimum RMSD for the two pairs of helices
//            Real BundlePairFeatures::calculate_pair_rmsd(core::pose::Pose pose1, std::pair<HelicalFragment,HelicalFragment> pair1, core::pose::Pose pose2,
//                                                         std::pair<HelicalFragment,HelicalFragment> pair2){
//                
//                //get the minimum helix size - this will be used as the RMSD window size
//                core::Size helix_1_min_size = min(pair1.first.get_size(), pair2.first.get_size());
//                core::Size helix_2_min_size = min(pair1.second.get_size(), pair2.second.get_size());
//                
//                std::map< core::id::AtomID, core::id::AtomID > atom_map;
//                atom_map.clear();
//                
//                Real min_rmsd=1000;
//                //iterate through all 4 helices to find minimal pair RMSD using window size
//                for(Size i=pair1.first.get_start(); i<=pair1.first.get_end()-helix_1_min_size; ++i){
//                    for(Size j=pair2.first.get_start(); j<=pair2.first.get_end()-helix_1_min_size; ++j){
//                        
//                        for(Size offset=0; offset<=helix_1_min_size; ++offset){
//                            core::id::AtomID const id1( pose1.residue(i+offset).atom_index("CA"), i+offset);
//                            core::id::AtomID const id2( pose2.residue(j+offset).atom_index("CA"), j+offset);
//                            atom_map.insert(std::pair<id::AtomID, id::AtomID>(id1, id2));
//                            
//                            core::id::AtomID const id3( pose1.residue(i+offset).atom_index("C"), i+offset);
//                            core::id::AtomID const id4( pose2.residue(j+offset).atom_index("C"), j+offset);
//                            atom_map.insert(std::pair<id::AtomID, id::AtomID>(id3, id4));
//                            
//                            core::id::AtomID const id5( pose1.residue(i+offset).atom_index("N"), i+offset);
//                            core::id::AtomID const id6( pose2.residue(j+offset).atom_index("N"), j+offset);
//                            atom_map.insert(std::pair<id::AtomID, id::AtomID>(id5, id6));
//                            
//                            core::id::AtomID const id7( pose1.residue(i+offset).atom_index("O"), i+offset);
//                            core::id::AtomID const id8( pose2.residue(j+offset).atom_index("O"), j+offset);
//                            atom_map.insert(std::pair<id::AtomID, id::AtomID>(id7, id8));
//                        }
//                        
//                        
//                        for(Size k=pair1.second.get_start(); k<=pair1.second.get_end()-helix_2_min_size; ++k){
//                            for(Size l=pair2.second.get_start(); l<=pair2.second.get_end()-helix_2_min_size; ++l){
//                                
//                                for(Size offset=0; offset<=helix_1_min_size; ++offset){
//                                    core::id::AtomID const id9( pose1.residue(k+offset).atom_index("CA"), k+offset);
//                                    core::id::AtomID const id10( pose2.residue(l+offset).atom_index("CA"), l+offset);
//                                    atom_map.insert(std::pair<id::AtomID, id::AtomID>(id9, id10));
//                                    
//                                    core::id::AtomID const id11( pose1.residue(k+offset).atom_index("C"), k+offset);
//                                    core::id::AtomID const id12( pose2.residue(l+offset).atom_index("C"), l+offset);
//                                    atom_map.insert(std::pair<id::AtomID, id::AtomID>(id11, id12));
//                                    
//                                    core::id::AtomID const id13( pose1.residue(k+offset).atom_index("N"), k+offset);
//                                    core::id::AtomID const id14( pose2.residue(l+offset).atom_index("N"), l+offset);
//                                    atom_map.insert(std::pair<id::AtomID, id::AtomID>(id13, id14));
//                                    
//                                    core::id::AtomID const id15( pose1.residue(k+offset).atom_index("O"), k+offset);
//                                    core::id::AtomID const id16( pose2.residue(l+offset).atom_index("O"), l+offset);
//                                    atom_map.insert(std::pair<id::AtomID, id::AtomID>(id15, id16));
//                                }
//                                
//                                Real rmsd = core::scoring::rms_at_all_corresponding_atoms(pose1, pose2, atom_map);
//                                if(rmsd < min_rmsd){
//                                    min_rmsd=rmsd;
//                                }
//                            }
//                        }
//                        
//                    }
//                }
//                return min_rmsd;
//            }
//            
//            ///@brief collect all the feature data for the pose
//            core::Size
//            BundlePairFeatures::report_features(
//                                                core::pose::Pose const & pose,
//                                                utility::vector1<bool> const & relevant_residues,
//                                                core::Size struct_id,
//                                                utility::sql_database::sessionOP db_session
//                                                ){
//                
//                std::string select_string =
//                "SELECT\n"
//                "   bundles.struct_id,\n"
//                "	bundles.bundle_id,\n"
//                "   helices1.helix_id,\n"
//                "	helices1.residue_begin AS helix_1_begin,\n"
//                "	helices1.residue_end AS helix_1_end,\n"
//                "   helices2.helix_id,\n"
//                "	helices2.residue_begin AS helix_2_begin,\n"
//                "	helices2.residue_end AS helix_2_end\n"
//                "FROM\n"
//                "	helix_bundles AS bundles\n"
//                "JOIN bundle_helices AS helices1 ON\n"
//                "   bundles.bundle_id = helices1.bundle_id\n"
//                "JOIN bundle_helices AS helices2 ON\n"
//                "   bundles.bundle_id = helices2.bundle_id\n"
//                "WHERE\n"
//                "   helices1.bundle_id = helices2.bundle_id AND\n"
//                "   helices1.helix_id <> helices2.helix_id";
//                
//                statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
//                result res1(basic::database::safely_read_from_database(select_statement));
//                result res2(basic::database::safely_read_from_database(select_statement));
//                
//                while(res1.next()){
//                    Size pair1_struct_id, pair1_bundle_id, pair1_helix1_id, pair1_helix1_begin, pair1_helix1_end, 
//                    pair1_helix2_id, pair1_helix2_begin, pair1_helix2_end;
//                    
//                    res1 >> pair1_struct_id >> pair1_bundle_id >> pair1_helix1_id >> pair1_helix1_begin >> pair1_helix1_end >>
//                    pair1_helix2_id >> pair1_helix2_begin >> pair1_helix2_end;
//                    
//                    std::pair<HelicalFragment, HelicalFragment> pair1(HelicalFragment(pair1_helix1_begin, pair1_helix1_end), 
//                                                                      HelicalFragment(pair1_helix2_begin, pair1_helix2_end));
//                    
//                    //This is needlessly expensive since we'll be recreating poses even when struct_id doesn't change
//                    core::pose::Pose pose1;
//                    pose_conformation_features_->load_into_pose(db_session, pair1_struct_id, pose1);
//                    protein_residue_conformation_features_->load_into_pose(db_session, pair1_struct_id, pose1);
//                    
//                    utility::vector1< std::pair<HelicalFragment,HelicalFragment> > bundle_pairs;
//                    while(res2.next()){
//                        
//                        Size pair2_struct_id, pair2_bundle_id, pair2_helix1_id, pair2_helix1_begin, pair2_helix1_end, 
//                        pair2_helix2_id, pair2_helix2_begin, pair2_helix2_end;
//                        
//                        res2 >> pair2_struct_id >> pair2_bundle_id >> pair2_helix1_id >> pair2_helix1_begin >> pair2_helix1_end >>
//                        pair2_helix2_id >> pair2_helix2_begin >> pair2_helix2_end;
//                        
//                        std::pair<HelicalFragment, HelicalFragment> pair2(HelicalFragment(pair2_helix1_begin, pair2_helix1_end), HelicalFragment(pair2_helix2_begin, pair2_helix2_end));
//                        
//                        
//                        //if these pairs are part of the same bundle, don't add their RMSD to the table
//                        if(pair1_struct_id != pair2_struct_id && pair1_bundle_id != pair2_bundle_id){
//                            
//                            //This is needlessly expensive since we'll be recreating poses even when struct_id doesn't change
//                            core::pose::Pose pose2;
//                            pose_conformation_features_->load_into_pose(db_session, pair2_struct_id, pose2);
//                            protein_residue_conformation_features_->load_into_pose(db_session, pair2_struct_id, pose2);
//                            
//                            Real rmsd = BundlePairFeatures::calculate_pair_rmsd(pose1, pair1, pose2, pair2);
//                            
//                            std::string insert_string = "INSERT INTO bundle_pair_comparisons VALUES (?,?,?,?,?,?,?,?,?)";
//                            statement insert_statement(basic::database::safely_prepare_statement(insert_string,db_session));
//                            insert_statement.bind(1,pair1_struct_id);
//                            insert_statement.bind(2,pair1_bundle_id);
//                            insert_statement.bind(3,pair1_helix1_id);
//                            insert_statement.bind(4,pair1_helix2_id);
//                            
//                            insert_statement.bind(5,pair2_struct_id);
//                            insert_statement.bind(6,pair2_bundle_id);
//                            insert_statement.bind(7,pair2_helix1_id);
//                            insert_statement.bind(8,pair2_helix2_id);
//                            
//                            insert_statement.bind(9,rmsd);
//                            basic::database::safely_write_to_database(insert_statement);
//                        }
//                    }
//                }
//                return 0;
//            }
//            
//        } //namespace helixAssembly
//    } //namespace features
//} //namespace protocols

