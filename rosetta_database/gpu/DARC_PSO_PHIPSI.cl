__kernel void Check_for_intersection(                                     
		                     __global float4 *rays,
                                     __global float4 *atoms,
                                     __global float *ray_scores,
                                     const unsigned int NUM_ATOMS,
                                     __global float *weights,
				     const unsigned int NUM_RAYS,
                                     __global float4 *atoms_maxmin_phipsi
                                     )
{
  int rayID = get_global_id(0);
  if(rayID >= NUM_RAYS) return;
  unsigned int atomID;
  float4 ray = rays[rayID];
  //scoring variables
  float const missing_point_weight = weights[0];
  float const steric_weight = weights[1];
  float const extra_point_weight = weights[2];
  int const NUM_PARTICLES = (int)weights[3];

  const int large_dist = 999.;
  float best_rho_sq, min_intersect_SQ, dist_deviation;
  //quadratic setup starts for ray
  float dirX = large_dist*sin(ray.x)*cos(ray.y);
  float dirY = large_dist*sin(ray.x)*sin(ray.y);
  float dirZ = large_dist*cos(ray.x);
  // setup our quadratic equation
  float a = (dirX*dirX) + (dirY*dirY) + (dirZ*dirZ);

  for(int particleID = 0; particleID < NUM_PARTICLES; particleID++){

    int ray_score_ID = (particleID*NUM_RAYS+rayID);
    best_rho_sq = 9999.;
    int atom_ID_start = particleID*NUM_ATOMS;
    int atom_ID_stop = (particleID+1)*NUM_ATOMS;
    for(atomID = atom_ID_start; atomID < atom_ID_stop; atomID ++) {

      float4 atom = atoms[atomID];
      float4 maxmin_phipsi = atoms_maxmin_phipsi[atomID];
      float curr_phi = ray.x;
      float curr_psi = ray.y;
      while ( curr_phi < maxmin_phipsi.z ) {
        curr_phi += 2*M_PI;
      }
      while ( curr_phi > maxmin_phipsi.x ) {
        curr_phi -= 2*M_PI;
      }
      while ( curr_psi < maxmin_phipsi.w ) {
        curr_psi += 2*M_PI;
      }
      while ( curr_psi > maxmin_phipsi.y ) {
        curr_psi -= 2*M_PI;
      }
      if ( curr_phi < maxmin_phipsi.z ) continue;
      if ( curr_psi < maxmin_phipsi.w ) continue;
      if ( curr_phi > maxmin_phipsi.x ) continue;
      if ( curr_psi > maxmin_phipsi.y ) continue;

      min_intersect_SQ = 9999.;
      //quadratic setup starts for atom
      float b = 2.0 * ( (dirX*(-atom.x)) + (dirY*(-atom.y)) + (dirZ*(-atom.z)) );
      float c = atom.x*atom.x + atom.y*atom.y + atom.z*atom.z - (atom.w * atom.w);
      // test for intersection
      float inside_sq = ( b * b ) - ( 4.0 * a * c );
      if(inside_sq > 0.0) {
        float inside = sqrt(inside_sq);
        float mu1 = -(b-inside) / ( 2.0 * a);
        float x1 =  mu1 * dirX;
        float y1 =  mu1 * dirY;
        float z1 =  mu1 * dirZ;
        float dist1_sq = x1*x1 + y1*y1 + z1*z1;
        float mu2 = -(b+inside) / ( 2.0 * a);
        float x2 = mu2 * dirX;
        float y2 = mu2 * dirY;
        float z2 = mu2 * dirZ;
        float dist2_sq = x2*x2 + y2*y2 + z2*z2;
	if(dist2_sq < dist1_sq){
	  min_intersect_SQ = dist2_sq;
	    }else{
	        min_intersect_SQ = dist1_sq;
		    }
      }
      if(min_intersect_SQ < best_rho_sq ){
        best_rho_sq  = min_intersect_SQ;
      }
    }
    float plaid_rho = 9999.;
    if(best_rho_sq < 9998.0){
      plaid_rho = sqrt(best_rho_sq );
    }
    if( (plaid_rho > 9998.0) && (ray.z > 0.001) ) {
      ray_scores[ray_score_ID] = missing_point_weight;
    }else if( (plaid_rho < 9998.0) && (ray.z < 0.001) ) {
      ray_scores[ray_score_ID] = extra_point_weight;
    }else if( (plaid_rho < 9998.0) && (ray.z > 0.001) ) {
      dist_deviation = plaid_rho - ray.z;
      if(dist_deviation < 0.0){
        dist_deviation = (ray.z - plaid_rho)*steric_weight;
      }
      ray_scores[ray_score_ID] = dist_deviation;
    }else{
      ray_scores[ray_score_ID] = 0.0;
    }
  }
}

__kernel void Get_scores(
                         __global float *ray_scores,
                         __global float *particle_scores,
                         const unsigned int NUM_RAYS,
							 const unsigned int num_particles
                         )
{
  int particleID = get_global_id(0);
  if(get_global_id(0) >= num_particles) return;
  int start = particleID*NUM_RAYS;
  int stop = start + NUM_RAYS;
  int ray_score_ID;
  float score = 0;
  int num = 0;
  for(ray_score_ID=start; ray_score_ID<stop; ray_score_ID++){
    if(ray_scores[ray_score_ID]>0.0){
      score += ray_scores[ray_score_ID];
      num++;
    }
  }
  particle_scores[particleID] = score/num;
}
