__kernel void Check_for_intersection(
                                     __global float4 *rays,
                                     __global float4 *atoms,
                                     __global float *ray_scores,
                                     const unsigned int NUM_ATOMS,
                                     __global float *weights,
                                     const unsigned int NUM_RAYS
                                     )
{
  const int large_dist = 999.0;
  int rayID = get_global_id(0);
  if(rayID>=NUM_RAYS) return;
  unsigned int atomID;
  float dirX, dirY, dirZ;
  float min_intersect_SQ, best_rho_sq;
  float4 ray = rays[rayID];
  //scoring variables
  float missing_point_weight = weights[0];
  float steric_weight = weights[1];
  float extra_point_weight = weights[2];
  int NUM_PARTICLES = (int)weights[3];
  //quadratic setup starts for ray
  dirX = large_dist*sin(ray.x)*cos(ray.y);
  dirY = large_dist*sin(ray.x)*sin(ray.y);
  dirZ = large_dist*cos(ray.x);
  // setup our quadratic equation
  float a = (dirX*dirX) + (dirY*dirY) + (dirZ*dirZ);
  for(int particleID = 0; particleID < NUM_PARTICLES; particleID++){
    //ray_score array: 0-99,999 particle 1, 100,000-199,999 particle 2, etc
    int ray_score_ID = (particleID*7596+rayID);
    best_rho_sq = 9999.;
    for(atomID = particleID*NUM_ATOMS; atomID < (particleID+1)*NUM_ATOMS; atomID ++) {
      min_intersect_SQ = 9999.;
      //quadratic setup starts for atom
      float4 atom = atoms[atomID];
      float b = 2.0 * ( (dirX*(-atom.x)) + (dirY*(-atom.y)) + (dirZ*(-atom.z)) );
      float c = atom.x*atom.x + atom.y*atom.y + atom.z*atom.z - (atom.w * atom.w);
      // test for intersection
      float inside_sq = ( b * b ) - ( 4.0 * a * c );
      if (inside_sq > 0) {
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
        if(dist2_sq < dist1_sq) min_intersect_SQ = dist2_sq;
        if(dist1_sq < dist2_sq) min_intersect_SQ = dist1_sq;
      }
      if(min_intersect_SQ<best_rho_sq) best_rho_sq = min_intersect_SQ;
    }
    
    float plaid_rho = 9999.;
    if(best_rho_sq<9998.){
      plaid_rho = sqrt(best_rho_sq);
    }
    ray_scores[ray_score_ID] = 0.0;
    if ( (plaid_rho > 9998.0) && (ray.z > 0.001) ) {
      ray_scores[ray_score_ID] = missing_point_weight;
    }
    if ( (plaid_rho < 9999.0) && (ray.z < 0.001) ) {
      ray_scores[ray_score_ID] = extra_point_weight;
    }
    if ( (plaid_rho < 9999.0) && (ray.z > 0.001) ) {
      float dist_deviation = plaid_rho - ray.z;
      if (dist_deviation < 0.0){
        dist_deviation = (ray.z - plaid_rho)*steric_weight;
      }
      ray_scores[ray_score_ID] = dist_deviation;
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
  if(particleID>=num_particles) return;
  int ray_score_ID;
  float score = 0;
  int num = 0;
  for(ray_score_ID=(particleID*NUM_RAYS); ray_score_ID<((particleID*NUM_RAYS)+NUM_RAYS); ray_score_ID++){
    if(ray_scores[ray_score_ID]>0){
        score = score + ray_scores[ray_score_ID];
      num = num + 1;
    }
  }
  particle_scores[particleID] = score/num;
}
