__kernel void Check_for_intersection(
                                     __global const float4 *rays,
                                     __global const float4 *atoms,
                                     __global float *ray_scores,
                                     const unsigned int NUM_ATOMS,
                                     __global float *weights
                                     )
{
  const int large_dist = 999.0;
  int rayID = get_global_id(0);
  unsigned int atomID;
  float dirX, dirY, dirZ;
  float distance, smallest;
  float4 ray;
  //scoring variables
  float missing_point_weight = weights[0];
  float steric_weight = weights[1];
  float extra_point_weight = weights[2];
  int NUM_PARTICLES = (int)weights[3];
  //quadratic setup starts for ray
  ray = rays[rayID];
  dirX = large_dist*sin(ray.x)*cos(ray.y);
  dirY = large_dist*sin(ray.x)*sin(ray.y);
  dirZ = large_dist*cos(ray.x);
  // setup our quadratic equation
  float a = (dirX*dirX) + (dirY*dirY) + (dirZ*dirZ);

  for(int particleID = 0; particleID < NUM_PARTICLES; particleID++){
    //ray_score array: 0-99,999 particle 1, 100,000-199,999 particle 2, etc
    int ray_score_ID = (particleID*100000+rayID);
    distance = 9999.;
    smallest = 9999.;
    float dis_deviation = 0;
    int atom_ID_start = particleID*NUM_ATOMS;
    int atom_ID_stop = (particleID+1)*NUM_ATOMS;
    for(atomID = atom_ID_start; atomID < atom_ID_stop; atomID ++) {
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
        if(dist1_sq < dist2_sq)
          distance = sqrt(dist1_sq);
        if(dist2_sq < dist1_sq)
          distance = sqrt(dist2_sq);
      }
      if(distance<smallest) smallest = distance;
    }
    ray_scores[ray_score_ID] = 0.0;
    if ( (ray.z > 0.001) && (smallest > 9998.0) ) {
      ray_scores[ray_score_ID] = missing_point_weight;
    }
    if ( (ray.z < 0.001) && (smallest < 9999.0) ) {
      ray_scores[ray_score_ID] = extra_point_weight;
    }
    if ( (ray.z > 0.001) && (smallest < 9999.0) ) {
      dist_deviation = smallest - ray.z;
      if (dist_deviation < 0){
        dist_deviation = (ray.z - smallest)*steric_weight;
      }
      ray_scores[ray_score_ID] = dist_deviation;
    }
  }
}
__kernel void Get_scores(
                         __global float *ray_scores,
                         __global float4 *particle_scores,
                         const unsigned int NUM_RAYS
                         )
{
  int particleID = get_global_id(0);
  particle_scores[particleID].x = 0;
  particle_scores[particleID].y = 0;
  int start = particleID*100000;
  int stop = start + NUM_RAYS;
  int ray_score_ID;
  float score = 0;
  int num = 0;
  for(ray_score_ID=start; ray_score_ID<stop; ray_score_ID++){
    if(ray_scores[ray_score_ID]>0){
      score = score + ray_scores[ray_score_ID];
      num = num + 1;
    }
  }
  particle_scores[particleID].z = score/num;
}
