//struct xyzmatrix;
//struct xyzvector;


float
calculate_helix_clash_score(
  __global struct xyzvector * node1_other_coords,
  int node1_n_other_coords,
  __global int * node1_center_coords,
  int node1_ncenters,
  __global struct xyzvector * node2_other_coords,
  int node2_n_other_coords,
  __global int * node2_center_coords,
  int node2_ncenters,
  __private struct xyzmatrix * UU//,
  //__global float * all_square_distances
)
{
  float closest_contact = 12345;

  //int count_square_distances = -1;
#define GROUPSIZE 2

  int n1_ncenter_groups = ( node1_ncenters + GROUPSIZE - 1 ) / GROUPSIZE;
  for ( int ii = 0; ii < n1_ncenter_groups; ++ii ) {
    struct xyzvector iicoords[ GROUPSIZE ];
    int iioffset = ii*GROUPSIZE;
    for ( int jj = 0; jj < GROUPSIZE; ++jj ) {
      if ( jj+iioffset < node1_ncenters ) { iicoords[jj] = node1_other_coords[ node1_center_coords[ jj+iioffset ]]; }
    }
    for ( int jj = 0; jj < GROUPSIZE; ++jj ) {
      if ( jj+iioffset >= node1_ncenters ) continue;
      for ( int kk = 0; kk < node2_n_other_coords; ++kk ) {
	struct xyzvector kkcoord = node2_other_coords[ kk ];
	struct xyzvector kkrotated = xyzmatrix_xyzvector_multiply_private( UU, &kkcoord );
	float dist2 = xyzvector_square_distance_private( &iicoords[jj+iioffset], & kkrotated );
	closest_contact = dist2 < closest_contact ? dist2 : closest_contact;
	//all_square_distances[ ++count_square_distances ] = dist2;
      }
    }
  }

  int n2_ncenter_groups = ( node2_ncenters + GROUPSIZE - 1 ) / GROUPSIZE;
  for ( int ii = 0; ii < n2_ncenter_groups; ++ii ) {
    struct xyzvector iirotated[ GROUPSIZE ];
    int iioffset = ii*GROUPSIZE;
    for ( int jj = 0; jj < GROUPSIZE; ++jj ) {
      if ( jj+iioffset < node2_ncenters ) {
	struct xyzvector jjcoord = node2_other_coords[ node2_center_coords[ jj+iioffset ]];
	iirotated[jj] = xyzmatrix_xyzvector_multiply_private( UU, &jjcoord );
      }
    }
    for ( int jj = 0; jj < GROUPSIZE; ++jj ) {
      if ( jj+iioffset >= node2_ncenters ) continue;
      for ( int kk = 0; kk < node1_n_other_coords; ++kk ) {
	struct xyzvector kkcoord = node1_other_coords[ kk ];
	float dist2 = xyzvector_square_distance_private( & iirotated[jj], & kkcoord );
	closest_contact = dist2 < closest_contact ? dist2 : closest_contact;
	//all_square_distances[ ++count_square_distances ] = dist2;
      }
    }
  }

#undef GROUPSIZE

  // TEMP DEBUG
  //for ( int ii = 0; ii < node2_n_other_coords; ++ii ) {
  //  struct xyzvector iicoord = node2_other_coords[ ii ];
  //  node2_other_coords[ ii ] = xyzmatrix_xyzvector_multiply_private( UU, &iicoord );
  //}
  return closest_contact;
}

void
calculate_rmsd_and_clash_score_for_bundles(
  __global struct xyzvector * node1_coords_at_origin,
  __global float * node1_coords_at_origin_moment_of_inertia,
  __global struct xyzvector * node2_coords_at_origin,
  __global float * node2_coords_at_origin_moment_of_inertia,
  int node_npoints,
  __global struct xyzvector * node1_other_coords,
  int node1_n_other_coords,
  __global int * node1_center_coords,
  int node1_ncenters,
  __global struct xyzvector * node2_other_coords,
  int node2_n_other_coords,
  __global int * node2_center_coords,
  int node2_ncenters,
  __private float * rms,
  __private float * clash//,
  //__global float * all_square_distances // <-- either debug collision calculation
  //__global struct xyzmatrix * UUout, /// <--- or use the next three to debug the findUU calculation
  //__global float * sigma3out ,
  //__global struct xyzmatrix * m_moment_out
)
{
  struct xyzmatrix UU;
  float sigma3;
  findUU_no_translation_to_origin(
    node1_coords_at_origin,
    node2_coords_at_origin,
    node_npoints, &UU, &sigma3//, m_moment_out
  );

  //*UUout = UU;
  //*sigma3out = sigma3;

  *rms = calculate_rmsd_fast(
    node1_coords_at_origin_moment_of_inertia,
    node2_coords_at_origin_moment_of_inertia,
    node_npoints, sigma3 );

  *clash = calculate_helix_clash_score(
    node1_other_coords, node1_n_other_coords, node1_center_coords, node1_ncenters,
    node2_other_coords, node2_n_other_coords, node2_center_coords, node2_ncenters,
    &UU//, all_square_distances
  );
}

__kernel
void
compute_rmsd_and_clash_scores(
  int block_size,
  int n_nodes_1,
  int n_nodes_2,
  int n_atoms_in_rms_calc,
  __global int * node1_bundle_inds,
  __global int * node2_bundle_inds,
  __global struct xyzvector * rmsd_coords_1,
  __global struct xyzvector * rmsd_coords_2,
  __global float * rmsd_coords_moment_of_inertia_1,
  __global float * rmsd_coords_moment_of_inertia_2,
  __global struct xyzvector * col_coords_1,
  __global struct xyzvector * col_coords_2,
  __global int * col_coords_offsets_1,
  __global int * col_coords_offsets_2,
  __global int * n_coords_for_col_calc_1,
  __global int * n_coords_for_col_calc_2,
  __global int * comparison_coord_ind_list_1,
  __global int * comparison_coord_ind_list_2,
  __global int * comparison_coord_ind_offset_1,
  __global int * comparison_coord_ind_offset_2,
  __global int * n_comparison_coords_1,
  __global int * n_comparison_coords_2,
  __global float * rmsd_table,
  __global float * collision_table,
  __global unsigned char * calc_rmsd_table//,
  //__global float * all_square_distances
  //__global struct xyzmatrix * UUout,
  //__global float * sigma3out,
  //__global struct xyzmatrix * m_moment_out
)
{
  int n1 = get_global_id(0);
  int n2 = get_global_id(1);

  if ( n1 < n_nodes_1 && n2 < n_nodes_2 ) {
    unsigned char calculate_this_pair = node1_bundle_inds[ n1 ] < node2_bundle_inds[ n2 ] ? 1 : 0;
    calc_rmsd_table[ n1*block_size + n2 ] = calculate_this_pair;

    if ( calculate_this_pair ) {
      float rms, collision;
      //printf( "%x %x :: %x %x\n", rmsd_coords_1, rmsd_coords_2, &rmsd_coords_1[ n1 * n_atoms_in_rms_calc ], &rmsd_coords_2[ n2 * n_atoms_in_rms_calc ] );
      calculate_rmsd_and_clash_score_for_bundles(
        &rmsd_coords_1[ n1 * n_atoms_in_rms_calc ], &rmsd_coords_moment_of_inertia_1[ n1 ],
        &rmsd_coords_2[ n2 * n_atoms_in_rms_calc ], &rmsd_coords_moment_of_inertia_2[ n2 ],
        n_atoms_in_rms_calc,
        &col_coords_1[ col_coords_offsets_1[ n1 ] ], n_coords_for_col_calc_1[ n1 ],
        &comparison_coord_ind_list_1[ comparison_coord_ind_offset_1[ n1 ] ], n_comparison_coords_1[ n1 ],
        &col_coords_2[ col_coords_offsets_2[ n2 ] ], n_coords_for_col_calc_2[ n2 ],
        &comparison_coord_ind_list_2[ comparison_coord_ind_offset_2[ n2 ] ], n_comparison_coords_2[ n2 ],
        &rms, &collision
        //, UUout, sigma3out, m_moment_out // <-- debugging variables
        //, all_square_distances
      );
  
      rmsd_table[ n1*block_size + n2 ] = rms;
      collision_table[ n1*block_size + n2 ] = collision;
    }
  } else if ( n1 < block_size && n2 < block_size ) {
    calc_rmsd_table[ n1*block_size + n2 ] = 0;
    rmsd_table[ n1*block_size + n2 ] = -1234; // debug
    collision_table[ n1*block_size + n2 ] = -1234; // debug
  }
}

/// find the center of mass of a list of coordinates used in a future RMS calculation
/// and translate them so that their center of mass is at the origin.  Next compute
/// the moment of inertia for those translated points.  Finally, translate a second
/// set of coordinates using the same translation vector as for the first set of
/// coordinates.
__kernel
void
translate_rmscoords_to_origin(
  int n_nodes,
  int n_rms_coords_per_node,
  __global struct xyzvector * rms_coords,
  __global float * node_moments_of_inertia,
  __global struct xyzvector * collision_coords,
  __global int * col_coords_offsets,
  __global int * n_coords_for_col_calc
)
{
  if ( get_global_id(0) < n_nodes ) {
    int node_id = get_global_id(0);
    struct xyzvector rms_coords_center_of_mass;
    float moment_of_inertia = transform_center_of_mass_to_origin(
      & rms_coords[ node_id * n_rms_coords_per_node ],
      & rms_coords_center_of_mass,
      n_rms_coords_per_node );
    int node_num_collision_coords = n_coords_for_col_calc[ node_id ];
    int offset = col_coords_offsets[ node_id ];
    for ( int ii = 0; ii < node_num_collision_coords; ++ii ) {
      struct xyzvector iicolcoord = collision_coords[ ii+offset ];
      iicolcoord.x_ -= rms_coords_center_of_mass.x_;
      iicolcoord.y_ -= rms_coords_center_of_mass.y_;
      iicolcoord.z_ -= rms_coords_center_of_mass.z_;
      collision_coords[ ii+offset ] = iicolcoord;
    }
    node_moments_of_inertia[ node_id ] = moment_of_inertia;
  }
}

__kernel
void
filter_good_rmsd_pairs(
  int block_size,
  int n_nodes_1,
  int n_nodes_2,
  float rmsd_filter,
  float clash_filter,
  __global float * rmsd_table,
  __global float * collision_table,
  __global unsigned char * calc_rmsd_table
)
{
  int n1 = get_global_id(0);
  int n2 = get_global_id(1);  
  if ( n1 < n_nodes_1 && n2 < n_nodes_2 && calc_rmsd_table[ n1*block_size + n2 ] == 1 ) {
    if ( rmsd_table[ n1*block_size + n2 ] > rmsd_filter || collision_table[ n1*block_size + n2 ] < clash_filter ) {
      calc_rmsd_table[ n1*block_size + n2 ] = 0;
    }
  }
}

/// Look at all the values computed and stored in the rmsd_table and clash_table and move only
/// the good rmsd/clash scores (and their indices) to the beginning of their respective tables so that
/// the CPU need only read back the good results from GPU.
__kernel
void
scan_good_rmsd_clash_pairs(
  int block_size,
  int n_nodes_1,
  int n_nodes_2,
  float max_acceptible_rms,
  float minimum_acceptible_clash_d2,
  __global unsigned char * calc_rmsd_table,
  __global float * rmsd_table,
  __global float * clash_table,
  __global int * good_indices_table,
  __global int * ngood_indices//,
  //__global int * scan_result_table, // <-- debug results
  //__global int * start_scan_results_table
)
{
  __local int values[ 256 ];

  int lid = get_local_id(0);
  int local_size  = get_local_size(0);
  int niterations = ( block_size * block_size + local_size - 1 ) / local_size;
  int last_offset = 0;

  for ( int ii = 0; ii < niterations; ++ii ) {   
    int iiid = ii*local_size + lid;
    int n1   = iiid / block_size;
    int n2   = iiid % block_size;
    if ( lid > 0 ) last_offset = 0; // reset every thread but thread 0

    bool iigood = false;
    float iirms = 0.0;
    float iiclash = 0.0;
    if ( n1 < n_nodes_1 && n2 < n_nodes_2 && calc_rmsd_table[ iiid ] != 0 ) {
      iirms = rmsd_table[ iiid ];
      iiclash = clash_table[ iiid ];
      iigood = iirms <= max_acceptible_rms && iiclash >= minimum_acceptible_clash_d2;
    }

    values[ lid ] = last_offset + ( iigood ? 1 : 0 );

    // debug: save the starting scan results
    //if ( iiid < block_size * block_size ) start_scan_results_table[ iiid ] = values[ lid ];
   
    // run a quck inclusive scan to get the compacted indices for the good rmsd/clash pairs
    for ( int offset = 1; offset < local_size; offset *= 2 ) {
      barrier(CLK_LOCAL_MEM_FENCE);
      int val;
      if ( lid-offset >= 0 ) val = values[lid-offset];

      barrier(CLK_LOCAL_MEM_FENCE);
      if ( lid-offset >= 0 ) values[ lid ] += val;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    // to get the exclusive scan result, take the value from position lid-1
    int my_val = last_offset; // thread 0 takes its value from the last value in the previous iteration
    if ( lid > 0 ) {
      my_val = values[ lid - 1 ];
    }


    // write the good results to the output rmsd, clash, and pair-index tables
    if ( iigood ) {
      int good_index = my_val;
      rmsd_table[ good_index ] = iirms;
      clash_table[ good_index ] = iiclash;
      good_indices_table[ good_index ] = iiid;
    }
    
    // debug: also save the scan result for this index
    //if ( iiid < block_size * block_size ) scan_result_table[ iiid ] = my_val;

    // node 0 keeps the last of the offsets
    if ( lid == 0 ) last_offset = values[ local_size - 1 ];
  }


  if ( lid == 0 ) *ngood_indices = last_offset;
}
