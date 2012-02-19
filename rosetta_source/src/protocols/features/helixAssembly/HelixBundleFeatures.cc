// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file protocols/features/helixAssembly/HelixBundleFeatures.cc
/// @brief Search through a pose for sets of 3 helices and generate RMSD calculations between all pairs of them
/// @author Tim Jacobs

//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

//Devel
#include <protocols/features/helixAssembly/HelixBundleFeatures.hh>
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility and basic
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//C++
#include <string>
#include <math.h>

//External Headers
#include <cppdb/frontend.h>

//Basic
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/helixAssembly.OptionKeys.gen.hh>


static basic::Tracer TR("protocols.features.helixAssembly.HelixBundleFeatures");

namespace protocols {
namespace features {
namespace helixAssembly {

using namespace std;
using namespace core;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

HelixBundleFeatures::HelixBundleFeatures() :
helix_cap_dist_cutoff_(12.0),
helix_contact_dist_cutoff_(10.0),
min_helix_size_(14)
{
    init_from_options();
}

void HelixBundleFeatures::init_from_options(){ 
    using namespace basic::options;

    if(option[OptionKeys::helixAssembly::min_helix_size].user()){
        min_helix_size_ = option[OptionKeys::helixAssembly::min_helix_size];
    }
    if(option[OptionKeys::helixAssembly::helix_contact_dist_cutoff].user()){
        helix_contact_dist_cutoff_ = option[OptionKeys::helixAssembly::helix_contact_dist_cutoff];
    }
    if(option[OptionKeys::helixAssembly::helix_cap_dist_cutoff].user()){
        helix_cap_dist_cutoff_ = option[OptionKeys::helixAssembly::helix_cap_dist_cutoff];
    }
}
    
string
HelixBundleFeatures::schema() const {
	return

    "CREATE TABLE IF NOT EXISTS helix_bundles (\n"
    "   bundle_id INTEGER PRIMARY KEY,\n"
    "   struct_id INTEGER,\n"
    "	FOREIGN KEY(struct_id)\n"
    "		REFERENCES structures(struct_id)\n"
    "		DEFERRABLE INITIALLY DEFERRED);"
    
    //does this need to be keyed on struct_id too?
    "CREATE TABLE IF NOT EXISTS bundle_helices (\n"
    "   helix_id INTEGER PRIMARY KEY,\n"
    "   bundle_id INTEGER,\n"
    "	residue_begin INTEGER,\n"
    "	residue_end INTEGER,\n"
    "	FOREIGN KEY(bundle_id)\n"
    "		REFERENCES helix_bundles(bundle_id)\n"
    "		DEFERRABLE INITIALLY DEFERRED,\n"
    "	FOREIGN KEY(residue_begin,residue_end)\n"
    "		REFERENCES residues(resNum,resNum)\n"
    "		DEFERRABLE INITIALLY DEFERRED);"
    ;
}

//Select all helical segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<HelicalFragment> HelixBundleFeatures::get_full_helices(Size struct_id, sessionOP db_session){
    std::string select_string =
    "SELECT\n"
    "	helices.helix_id,\n"
    "	helices.residue_begin,\n"
    "	helices.residue_end\n"
    "FROM\n"
    "	helix_segments as helices\n"
    "WHERE\n"
    "	helices.struct_id = ?;";
    
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));
    
    utility::vector1<HelicalFragment> all_helices;
	while(res.next()){
		Size helix_id, residue_begin, residue_end;
		res >> helix_id >> residue_begin >> residue_end;
        all_helices.push_back(HelicalFragment(residue_begin, residue_end));
    }
    
    return all_helices;
}

//Check three helical fragments to ensure they make contacts with eachother. The requirement here is that each helical residue must have at least
//half of its residues within 10 angstroms to any residue on each of the other two helices
bool HelixBundleFeatures::checkHelixContacts(Pose const & pose, HelicalFragment helix_1, HelicalFragment helix_2, HelicalFragment helix_3){
    
    Size dist_sq_cutoff = (Size)pow((double)helix_contact_dist_cutoff_, 2);
    
    core::Size helix_1_and_2_contacts = 0;
    core::Size helix_1_and_3_contacts = 0;
    
    for(Size i=helix_1.get_start(); i<=helix_1.get_end(); i++){
        
        //bool helix_2_pass=false;       
        for(Size j=helix_2.get_start(); j<=helix_2.get_end(); j++){
            if(pose.residue(i).atom("CA").xyz().distance_squared(pose.residue(j).atom("CA").xyz()) < dist_sq_cutoff){
                helix_1_and_2_contacts++;
                break;
            }
        }
        
        //bool helix_3_pass=false;
        for(Size j=helix_3.get_start(); j<helix_3.get_end(); j++){
            if(pose.residue(i).atom("CA").xyz().distance_squared(pose.residue(j).atom("CA").xyz()) < dist_sq_cutoff){
                helix_1_and_3_contacts++;
                break;
            }
        }
    }
    
    if(helix_1_and_2_contacts < (helix_1.get_size()/2) || helix_1_and_3_contacts < (helix_1.get_size()/2)){
        return false;
    }
    
    //Now check interactions between helices 2 and 3
    core::Size helix_2_and_3_contacts = 0;
    for(Size i=helix_2.get_start(); i<=helix_2.get_end(); i++){
        
        //bool helix_2_pass=false;       
        for(Size j=helix_3.get_start(); j<=helix_3.get_end(); j++){
            if(pose.residue(i).atom("CA").xyz().distance_squared(pose.residue(j).atom("CA").xyz()) < dist_sq_cutoff){
                helix_2_and_3_contacts++;
                break;
            }
        }
    }
    if(helix_2_and_3_contacts < (helix_2.get_size()/2)){
        return false;
    }
    return true;
}

///@brief collect all the feature data for the pose
core::Size
HelixBundleFeatures::report_features(
                                     core::pose::Pose const & pose,
                                     utility::vector1<bool> const & relevant_residues,
                                     core::Size struct_id,
                                     utility::sql_database::sessionOP db_session
                                     ){
    
    utility::vector1<HelicalFragment> all_helices = get_full_helices(struct_id, db_session);
        
    //This function gets really ugly right here and should probably be further functionalized
    core::Size bundle_counter=0;
    Size dist_sq_cutoff = (Size)pow((double)helix_cap_dist_cutoff_, 2);
    
    for(Size i=1; i<=all_helices.size(); ++i){
        for(Size j=i+1; j<=all_helices.size(); ++j){
            for(Size k=j+1; k<=all_helices.size(); ++k){
                               
                core::Size max_bundle_residues=0;
                
                //Helical fragment holders for final 3-helix bundle
                HelicalFragment helix_i(0,0);
                HelicalFragment helix_j(0,0);
                HelicalFragment helix_k(0,0);
                
                
                //If all three helices are parallel (Note: This code could potentially be much faster if I only loop over each helical residue once,
                //but for ease and because this will only be run once, I'm not doing it that way)
                bool found_start=false;
                bool found_end=false;
                bool broke=false;
                
                core::Size helix_i_start=0;
                core::Size helix_j_start=0;
                core::Size helix_k_start=0;
                
                core::Size helix_i_end=0;
                core::Size helix_j_end=0;
                core::Size helix_k_end=0;
                
                for(Size helix_i_offset=0; helix_i_offset<all_helices[i].get_size(); helix_i_offset++){
                    if(broke){break;}
                    for(Size helix_j_offset=0; helix_j_offset<all_helices[j].get_size(); helix_j_offset++){
                        if(broke){break;}
                        for(Size helix_k_offset=0; helix_k_offset<all_helices[k].get_size(); helix_k_offset++){
                            
                            //If we haven't found a suitable start position for the 3 helices, check distances between current helical residues
                            if(!found_start){
                                //distance between n-term of helix-i and n-term of helix-j
                                DistanceSquared start_distance_1 = pose.residue(all_helices[i].get_start()+helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[j].get_start()+helix_j_offset).atom("CA").xyz());
                                
                                //distance between n-term of helix-i and n-term of helix-k
                                DistanceSquared start_distance_2 = pose.residue(all_helices[i].get_start()+helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_start()+helix_k_offset).atom("CA").xyz());
                                
                                //distance between n-term of helix-j and n-term of helix-k
                                DistanceSquared start_distance_3 = pose.residue(all_helices[j].get_start()+helix_j_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_start()+helix_k_offset).atom("CA").xyz());
                                
                                if(start_distance_1 <= dist_sq_cutoff && start_distance_2 <= dist_sq_cutoff && start_distance_3 <= dist_sq_cutoff){
                                    found_start = true;
                                    helix_i_start=(all_helices[i].get_start()+helix_i_offset);
                                    helix_j_start=(all_helices[j].get_start()+helix_j_offset);
                                    helix_k_start=(all_helices[k].get_start()+helix_k_offset);
                                }
                            }
                            
                            if(!found_end){
                                //distance between c-term of helix-i and c-term of helix-j
                                DistanceSquared end_distance_1 = pose.residue(all_helices[i].get_end()-helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[j].get_end()-helix_j_offset).atom("CA").xyz());
                                
                                //distance between c-term of helix-i and c-term of helix-k
                                DistanceSquared end_distance_2 = pose.residue(all_helices[i].get_end()-helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_end()-helix_k_offset).atom("CA").xyz());
                                
                                //distance between c-term of helix-j and c-term of helix-k
                                DistanceSquared end_distance_3 = pose.residue(all_helices[j].get_end()-helix_j_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_end()-helix_k_offset).atom("CA").xyz());
                                
                                if(end_distance_1 <= dist_sq_cutoff && end_distance_2 <= dist_sq_cutoff && end_distance_3 <= dist_sq_cutoff){
                                    found_end = true;
                                    helix_i_end=(all_helices[i].get_end()-helix_i_offset);
                                    helix_j_end=(all_helices[j].get_end()-helix_j_offset);
                                    helix_k_end=(all_helices[k].get_end()-helix_k_offset);
                                }
                            }
                            
                            //If we found start and end residues that satisfy cap distance requirements, then check helix contacts
                            if(found_start && found_end){
                                HelicalFragment temp_helix_i(min(helix_i_start, helix_i_end), max(helix_i_start, helix_i_end));
                                HelicalFragment temp_helix_j(min(helix_j_start, helix_j_end), max(helix_j_start, helix_j_end));
                                HelicalFragment temp_helix_k(min(helix_k_start, helix_k_end), max(helix_k_start, helix_k_end));
                                
                                TR.Debug << "Unchecked parallel helix found: \ni:" << temp_helix_i.get_start() << "," << temp_helix_i.get_end() 
                                << "\n j:" << temp_helix_j.get_start() << "," << temp_helix_j.get_end()
                                << "\n k:" << temp_helix_k.get_start() << "," << temp_helix_k.get_end() << endl;
                                
                                if(temp_helix_i.get_size() >= min_helix_size_ &&
                                   temp_helix_j.get_size() >= min_helix_size_ &&
                                   temp_helix_k.get_size() >= min_helix_size_){
                                    
                                    TR.Debug << "Size check passed!" << endl;
                                    
                                    //If the found start and ends don't satisfy this closeness check then we want to keep looking in these
                                    //helices for different starts and ends
                                    if(!checkHelixContacts(pose, temp_helix_i, temp_helix_j, temp_helix_k)){
                                        found_start=false;
                                        found_end=false;
                                        
                                        TR.Debug << "Close helix contacts failed!" << endl;
                                    }
                                    
                                    //Only take these helices if they are bigger than helices from other 
                                    else if(temp_helix_i.get_size() + temp_helix_j.get_size() + temp_helix_k.get_size() > max_bundle_residues){
                                        helix_i=temp_helix_i;
                                        helix_j=temp_helix_j;
                                        helix_k=temp_helix_k;
                                        max_bundle_residues=temp_helix_i.get_size() + temp_helix_j.get_size() + temp_helix_k.get_size();
                                        TR.Debug << "Found new parallel bundle i-" << temp_helix_i.get_size() << " j-" << temp_helix_j.get_size() << " k-" << temp_helix_k.get_size() << endl;
                                        broke=true;
                                        break;
                                    }
                                }
                                else{
                                    broke=true;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                
                
                //If helix-j is antiparallel to helix-i and helix-k
                found_start=false;
                found_end=false;
                broke=false;
                
                helix_i_start=0;
                helix_j_start=0;
                helix_k_start=0;
                
                helix_i_end=0;
                helix_j_end=0;
                helix_k_end=0;
                
                for(Size helix_i_offset=0; helix_i_offset<all_helices[i].get_size(); helix_i_offset++){
                    if(broke){break;}
                    for(Size helix_j_offset=0; helix_j_offset<all_helices[j].get_size(); helix_j_offset++){
                        if(broke){break;}
                        for(Size helix_k_offset=0; helix_k_offset<all_helices[k].get_size(); helix_k_offset++){
                            
                            //If we haven't found a suitable start position for the 3 helices, check distances between current helical residues
                            if(!found_start){
                                //distance between n-term of helix-i and c-term of helix-j
                                DistanceSquared start_distance_1 = pose.residue(all_helices[i].get_start()+helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[j].get_end()-helix_j_offset).atom("CA").xyz());
                                
                                //distance between n-term of helix-i and n-term of helix-k
                                DistanceSquared start_distance_2 = pose.residue(all_helices[i].get_start()+helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_start()+helix_k_offset).atom("CA").xyz());
                                
                                //distance between c-term of helix-j and n-term of helix-k
                                DistanceSquared start_distance_3 = pose.residue(all_helices[j].get_end()-helix_j_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_start()+helix_k_offset).atom("CA").xyz());
                                
                                if(start_distance_1 <= dist_sq_cutoff && start_distance_2 <= dist_sq_cutoff && start_distance_3 <= dist_sq_cutoff){
                                    found_start = true;
                                    helix_i_start=(all_helices[i].get_start()+helix_i_offset);
                                    helix_j_start=(all_helices[j].get_end()-helix_j_offset);
                                    helix_k_start=(all_helices[k].get_start()+helix_k_offset);
                                }
                            }
                            
                            if(!found_end){
                                //distance between c-term of helix-i and n-term of helix-j
                                DistanceSquared end_distance_1 = pose.residue(all_helices[i].get_end()-helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[j].get_start()+helix_j_offset).atom("CA").xyz());
                                
                                //distance between c-term of helix-i and c-term of helix-k
                                DistanceSquared end_distance_2 = pose.residue(all_helices[i].get_end()-helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_end()-helix_k_offset).atom("CA").xyz());
                                
                                //distance between c-term of helix-j and c-term of helix-k
                                DistanceSquared end_distance_3 = pose.residue(all_helices[j].get_end()-helix_j_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_end()-helix_k_offset).atom("CA").xyz());
                                
                                if(end_distance_1 <= dist_sq_cutoff && end_distance_2 <= dist_sq_cutoff && end_distance_3 <= dist_sq_cutoff){
                                    found_end = true;
                                    helix_i_end=(all_helices[i].get_end()-helix_i_offset);
                                    helix_j_end=(all_helices[j].get_start()+helix_j_offset);
                                    helix_k_end=(all_helices[k].get_end()-helix_k_offset);
                                }
                            }
                            
                            //If we found start and end residues that satisfy cap distance requirements, then check helix contacts
                            if(found_start && found_end){
                                HelicalFragment temp_helix_i(min(helix_i_start, helix_i_end), max(helix_i_start, helix_i_end));
                                HelicalFragment temp_helix_j(min(helix_j_start, helix_j_end), max(helix_j_start, helix_j_end));
                                HelicalFragment temp_helix_k(min(helix_k_start, helix_k_end), max(helix_k_start, helix_k_end));
                                
                                TR.Debug << "Unchecked j-flipped helix found: \ni:" << temp_helix_i.get_start() << "," << temp_helix_i.get_end() 
                                << "\n j:" << temp_helix_j.get_start() << "," << temp_helix_j.get_end()
                                << "\n k:" << temp_helix_k.get_start() << "," << temp_helix_k.get_end() << endl;
                                
                                if(temp_helix_i.get_size() >= min_helix_size_ &&
                                   temp_helix_j.get_size() >= min_helix_size_ &&
                                   temp_helix_k.get_size() >= min_helix_size_){
                                    
                                    TR.Debug << "Size check passed!" << endl;
                                    
                                    //If the found start and ends don't satisfy this closeness check then we want to keep looking in these
                                    //helices for different starts and ends
                                    if(!checkHelixContacts(pose, temp_helix_i, temp_helix_j, temp_helix_k)){
                                        found_start=false;
                                        found_end=false;
                                        TR.Debug << "Close helix contacts failed!" << endl;
                                    }
                                    
                                    //Only take these helices if they are bigger than helices from other 
                                    else if(temp_helix_i.get_size() + temp_helix_j.get_size() + temp_helix_k.get_size() > max_bundle_residues){
                                        helix_i=temp_helix_i;
                                        helix_j=temp_helix_j;
                                        helix_k=temp_helix_k;
                                        max_bundle_residues=temp_helix_i.get_size() + temp_helix_j.get_size() + temp_helix_k.get_size();
                                        TR.Debug << "Found new j-flipped bundle i-" << temp_helix_i.get_size() << " j-" << temp_helix_j.get_size() << " k-" << temp_helix_k.get_size() << endl;
                                        broke=true;
                                        break;
                                    }
                                }
                                else{
                                    broke=true;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                
                //If helix-k is antiparallel to helix-i and helix-j
                found_start=false;
                found_end=false;
                broke=false;
                
                helix_i_start=0;
                helix_j_start=0;
                helix_k_start=0;
                
                helix_i_end=0;
                helix_j_end=0;
                helix_k_end=0;
                
                for(Size helix_i_offset=0; helix_i_offset<all_helices[i].get_size(); helix_i_offset++){
                    if(broke){break;}
                    for(Size helix_j_offset=0; helix_j_offset<all_helices[j].get_size(); helix_j_offset++){
                        if(broke){break;}
                        for(Size helix_k_offset=0; helix_k_offset<all_helices[k].get_size(); helix_k_offset++){
                            
                            //If we haven't found a suitable start position for the 3 helices, check distances between current helical residues
                            if(!found_start){
                                //distance between n-term of helix-i and n-term of helix-j
                                DistanceSquared start_distance_1 = pose.residue(all_helices[i].get_start()+helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[j].get_start()+helix_j_offset).atom("CA").xyz());
                                
                                //distance between n-term of helix-i and c-term of helix-k
                                DistanceSquared start_distance_2 = pose.residue(all_helices[i].get_start()+helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_end()-helix_k_offset).atom("CA").xyz());
                                
                                //distance between n-term of helix-j and c-term of helix-k
                                DistanceSquared start_distance_3 = pose.residue(all_helices[j].get_start()+helix_j_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_end()-helix_k_offset).atom("CA").xyz());
                                
                                if(start_distance_1 <= dist_sq_cutoff && start_distance_2 <= dist_sq_cutoff && start_distance_3 <= dist_sq_cutoff){
                                    found_start = true;
                                    helix_i_start=(all_helices[i].get_start()+helix_i_offset);
                                    helix_j_start=(all_helices[j].get_start()+helix_j_offset);
                                    helix_k_start=(all_helices[k].get_end()-helix_k_offset);
                                }
                            }
                            
                            if(!found_end){
                                //distance between c-term of helix-i and c-term of helix-j
                                DistanceSquared end_distance_1 = pose.residue(all_helices[i].get_end()-helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[j].get_end()-helix_j_offset).atom("CA").xyz());
                                
                                //distance between c-term of helix-i and n-term of helix-k
                                DistanceSquared end_distance_2 = pose.residue(all_helices[i].get_end()-helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_start()+helix_k_offset).atom("CA").xyz());
                                
                                //distance between c-term of helix-j and n-term of helix-k
                                DistanceSquared end_distance_3 = pose.residue(all_helices[j].get_end()-helix_j_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_start()+helix_k_offset).atom("CA").xyz());
                                
                                if(end_distance_1 <= dist_sq_cutoff && end_distance_2 <= dist_sq_cutoff && end_distance_3 <= dist_sq_cutoff){
                                    found_end = true;
                                    helix_i_end=(all_helices[i].get_end()-helix_i_offset);
                                    helix_j_end=(all_helices[j].get_end()-helix_j_offset);
                                    helix_k_end=(all_helices[k].get_start()+helix_k_offset);
                                }
                            }
                            
                            //If we found start and end residues that satisfy cap distance requirements, then check helix contacts
                            if(found_start && found_end){
                                HelicalFragment temp_helix_i(min(helix_i_start, helix_i_end), max(helix_i_start, helix_i_end));
                                HelicalFragment temp_helix_j(min(helix_j_start, helix_j_end), max(helix_j_start, helix_j_end));
                                HelicalFragment temp_helix_k(min(helix_k_start, helix_k_end), max(helix_k_start, helix_k_end));
                                
                                TR.Debug << "Unchecked k-flipped helix found: \ni:" << temp_helix_i.get_start() << "," << temp_helix_i.get_end() 
                                << "\n j:" << temp_helix_j.get_start() << "," << temp_helix_j.get_end()
                                << "\n k:" << temp_helix_k.get_start() << "," << temp_helix_k.get_end() << endl;
                                
                                if(temp_helix_i.get_size() >= min_helix_size_ &&
                                   temp_helix_j.get_size() >= min_helix_size_ &&
                                   temp_helix_k.get_size() >= min_helix_size_){
                                    
                                    TR.Debug << "Size check passed!" << endl;
                                    
                                    //If the found start and ends don't satisfy this closeness check then we want to keep looking in these
                                    //helices for different starts and ends
                                    if(!checkHelixContacts(pose, temp_helix_i, temp_helix_j, temp_helix_k)){
                                        found_start=false;
                                        found_end=false;
                                        TR.Debug << "Close helix contacts failed!" << endl;
                                    }
                                    
                                    //Only take these helices if they are bigger than helices from other 
                                    else if(temp_helix_i.get_size() + temp_helix_j.get_size() + temp_helix_k.get_size() > max_bundle_residues){
                                        helix_i=temp_helix_i;
                                        helix_j=temp_helix_j;
                                        helix_k=temp_helix_k;
                                        max_bundle_residues=temp_helix_i.get_size() + temp_helix_j.get_size() + temp_helix_k.get_size();
                                        TR.Debug << "Found new k-flipped bundle i-" << temp_helix_i.get_size() << " j-" << temp_helix_j.get_size() << " k-" << temp_helix_k.get_size() << endl;
                                        broke=true;
                                        break;
                                    }
                                }
                                else{
                                    broke=true;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                
                //If helix-j & helix-k are both antiparallel to helix-i
                found_start=false;
                found_end=false;
                broke=false;
                
                helix_i_start=0;
                helix_j_start=0;
                helix_k_start=0;
                
                helix_i_end=0;
                helix_j_end=0;
                helix_k_end=0;
                
                for(Size helix_i_offset=0; helix_i_offset<all_helices[i].get_size(); helix_i_offset++){
                    if(broke){break;}
                    for(Size helix_j_offset=0; helix_j_offset<all_helices[j].get_size(); helix_j_offset++){
                        if(broke){break;}
                        for(Size helix_k_offset=0; helix_k_offset<all_helices[k].get_size(); helix_k_offset++){
                            
                            //If we haven't found a suitable start position for the 3 helices, check distances between current helical residues
                            if(!found_start){
                                //distance between n-term of helix-i and c-term of helix-j
                                DistanceSquared start_distance_1 = pose.residue(all_helices[i].get_start()+helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[j].get_end()-helix_j_offset).atom("CA").xyz());
                                
                                //distance between n-term of helix-i and c-term of helix-k
                                DistanceSquared start_distance_2 = pose.residue(all_helices[i].get_start()+helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_end()-helix_k_offset).atom("CA").xyz());
                                
                                //distance between c-term of helix-j and c-term of helix-k
                                DistanceSquared start_distance_3 = pose.residue(all_helices[j].get_end()-helix_j_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_end()-helix_k_offset).atom("CA").xyz());
                                
                                if(start_distance_1 <= dist_sq_cutoff && start_distance_2 <= dist_sq_cutoff && start_distance_3 <= dist_sq_cutoff){
                                    found_start = true;
                                    helix_i_start=(all_helices[i].get_start()+helix_i_offset);
                                    helix_j_start=(all_helices[j].get_end()-helix_j_offset);
                                    helix_k_start=(all_helices[k].get_end()-helix_k_offset);
                                }
                            }
                            
                            if(!found_end){
                                //distance between c-term of helix-i and n-term of helix-j
                                DistanceSquared end_distance_1 = pose.residue(all_helices[i].get_end()-helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[j].get_start()+helix_j_offset).atom("CA").xyz());
                                
                                //distance between c-term of helix-i and n-term of helix-k
                                DistanceSquared end_distance_2 = pose.residue(all_helices[i].get_end()-helix_i_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_start()+helix_k_offset).atom("CA").xyz());
                                
                                //distance between c-term of helix-j and n-term of helix-k
                                DistanceSquared end_distance_3 = pose.residue(all_helices[j].get_start()+helix_j_offset).atom("CA").xyz().distance_squared(                                                                                                                               pose.residue(all_helices[k].get_start()+helix_k_offset).atom("CA").xyz());
                                
                                if(end_distance_1 <= dist_sq_cutoff && end_distance_2 <= dist_sq_cutoff && end_distance_3 <= dist_sq_cutoff){
                                    found_end = true;
                                    helix_i_end=(all_helices[i].get_end()-helix_i_offset);
                                    helix_j_end=(all_helices[j].get_start()+helix_j_offset);
                                    helix_k_end=(all_helices[k].get_start()+helix_k_offset);
                                }
                            }
                            
                            //If we found start and end residues that satisfy cap distance requirements, then check helix contacts
                            if(found_start && found_end){
                                HelicalFragment temp_helix_i(min(helix_i_start, helix_i_end), max(helix_i_start, helix_i_end));
                                HelicalFragment temp_helix_j(min(helix_j_start, helix_j_end), max(helix_j_start, helix_j_end));
                                HelicalFragment temp_helix_k(min(helix_k_start, helix_k_end), max(helix_k_start, helix_k_end));
                                
                                TR.Debug << "Unchecked j&k flipped helix found: \ni:" << temp_helix_i.get_start() << "," << temp_helix_i.get_end() 
                                << "\n j:" << temp_helix_j.get_start() << "," << temp_helix_j.get_end()
                                << "\n k:" << temp_helix_k.get_start() << "," << temp_helix_k.get_end() << endl;
                                
                                //Quit here if any of the helices are too small (they only get smaller)
                                if(temp_helix_i.get_size() >= min_helix_size_ &&
                                   temp_helix_j.get_size() >= min_helix_size_ &&
                                   temp_helix_k.get_size() >= min_helix_size_){
                                    
                                    TR.Debug << "Size check passed!" << endl;
                                    
                                    //If the found start and ends don't satisfy this closeness check then we want to keep looking in these
                                    //helices for different starts and ends
                                    if(!checkHelixContacts(pose, temp_helix_i, temp_helix_j, temp_helix_k)){
                                        found_start=false;
                                        found_end=false;
                                        TR.Debug << "Close helix contacts failed!" << endl;
                                    }
                                    
                                    //Only take these helices if they are bigger than helices from other 
                                    else if(temp_helix_i.get_size() + temp_helix_j.get_size() + temp_helix_k.get_size() > max_bundle_residues){
                                        helix_i=temp_helix_i;
                                        helix_j=temp_helix_j;
                                        helix_k=temp_helix_k;
                                        max_bundle_residues=temp_helix_i.get_size() + temp_helix_j.get_size() + temp_helix_k.get_size();
                                        TR.Debug << "Found new j&k flipped bundle i-" << temp_helix_i.get_size() << " j-" << temp_helix_j.get_size() << " k-" << temp_helix_k.get_size() <<endl;
                                        broke=true;
                                        break;
                                    }
                                }
                                else{
                                    broke=true;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                //If any of the different orientations produced a bundle, print it
                if(helix_i.get_size() > 0){
                    
                    bundle_counter++;
                    TR << "saving bundle: " << bundle_counter << endl;
                                     
                    string bundle_insert =  "INSERT INTO helix_bundles VALUES (?,?);";
                    statement bundle_insert_stmt(basic::database::safely_prepare_statement(bundle_insert,db_session));
                    bundle_insert_stmt.bind(1,NULL);//auto-increment
                    bundle_insert_stmt.bind(2,struct_id);
                    basic::database::safely_write_to_database(bundle_insert_stmt);
                    
                    //Get bundle primary key               
                    core::Size bundle_id(bundle_insert_stmt.last_insert_id());
                                        
                    string helix_insert =  "INSERT INTO bundle_helices VALUES (?,?,?,?);";
                    statement helix_1_insert_stmt(basic::database::safely_prepare_statement(helix_insert,db_session));
                    helix_1_insert_stmt.bind(1,NULL);
                    helix_1_insert_stmt.bind(2,bundle_id);
                    helix_1_insert_stmt.bind(3,helix_i.get_start());
                    helix_1_insert_stmt.bind(4,helix_i.get_end());
                    basic::database::safely_write_to_database(helix_1_insert_stmt); 
                    
                    statement helix_2_insert_stmt(basic::database::safely_prepare_statement(helix_insert,db_session));
                    helix_2_insert_stmt.bind(1,NULL);
                    helix_2_insert_stmt.bind(2,bundle_id);
                    helix_2_insert_stmt.bind(3,helix_j.get_start());
                    helix_2_insert_stmt.bind(4,helix_j.get_end());
                    basic::database::safely_write_to_database(helix_2_insert_stmt);

                    statement helix_3_insert_stmt(basic::database::safely_prepare_statement(helix_insert,db_session));
                    helix_3_insert_stmt.bind(1,NULL);
                    helix_3_insert_stmt.bind(2,bundle_id);
                    helix_3_insert_stmt.bind(3,helix_k.get_start());
                    helix_3_insert_stmt.bind(4,helix_k.get_end());
                    basic::database::safely_write_to_database(helix_3_insert_stmt);
                }
            }
        }
    }
    
    TR << "Done saving helices" << endl;
    
    return 0;
}

} //namespace helixAssembly
} //namespace features
} //namespace protocols
