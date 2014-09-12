Real ligand_rmsd (pose::Pose const & pose1, pose::Pose const & pose2) {

	core::Size ref_res_num = 0;
  for ( int j = 1, resnum = pose1.total_residue(); j <= resnum; ++j ) {
    if (!pose1.residue(j).is_protein()){
      ref_res_num = j;
      break;
    }
  }
  if (ref_res_num == 0){
		std::cout<<"Error, no reference ligand for RMSD calculation" << std::endl;
    exit(1);
  }
	core::conformation::Residue const & pose1_rsd = pose1.conformation().residue(ref_res_num);

	core::Size inp_res_num = 0;
  for ( int j = 1, resnum = pose2.total_residue(); j <= resnum; ++j ) {
    if (!pose2.residue(j).is_protein()){
      inp_res_num = j;
      break;
    }
  }
  if (inp_res_num == 0){
		std::cout<<"Error, no input ligand for RMSD calculation" << std::endl;
    exit(1);
  }
	core::conformation::Residue const & pose2_rsd = pose2.conformation().residue(inp_res_num);


	//CALCULATE RMSD
	//IMPORTANT:this rmsd calculation does not consider symmetry
	//assert input & reference ligand should have same num of atoms
	if(!pose2_rsd.nheavyatoms() == pose1_rsd.nheavyatoms()){
		std::cout<<"Error, input & reference ligand should have same no. of atoms in the same order" << std::endl;
    exit(1);
  }
	else{
		core::Real rmsd, dist_sum(0.);
		Size j = 1;
		for(Size i = 1; i <= pose1_rsd.nheavyatoms(); ++i ) {
			if ( j <= pose2_rsd.nheavyatoms() ) {
				assert ( pose1_rsd.atom_name(i) == pose2_rsd.atom_name(j) );
				core::Real x_dist =  ( (pose1_rsd.atom(i).xyz()(1) - pose2_rsd.atom(j).xyz()(1)) * (pose1_rsd.atom(i).xyz()(1) - pose2_rsd.atom(j).xyz()(1)) );
				core::Real y_dist =  ( (pose1_rsd.atom(i).xyz()(2) - pose2_rsd.atom(j).xyz()(2)) * (pose1_rsd.atom(i).xyz()(2) - pose2_rsd.atom(j).xyz()(2)) );
				core::Real z_dist =  ( (pose1_rsd.atom(i).xyz()(3) - pose2_rsd.atom(j).xyz()(3)) * (pose1_rsd.atom(i).xyz()(3) - pose2_rsd.atom(j).xyz()(3)) );
				dist_sum += x_dist + y_dist + z_dist;
			}
			++j;
		}
		rmsd = sqrt(dist_sum/pose1_rsd.nheavyatoms());
		std::cout<<"Ligand RMSD:	" << rmsd << std::endl;
	}
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
	return 0;
}
