__kernel void Check_for_intersection(__global const float *phiC,
                                     __global const float *psiC,
                                     __global const float *phiS, 
                                     __global const float *psiS, 
                                     __global const float *atomX, 
                                     __global const float *atomY, 
                                     __global const float *atomZ, 
                                     __global const float *radii, 
                                     __global float *distances,
                                     __global const unsigned int *NUM_ATOMS){

    const int large_dist = 999.0;
        
    register int rayID = get_global_id(0);
    register int atomID, i, j, l;
    register float dirX, dirY, dirZ;
    register float distance = 9999.0;
    register float smallest = 9999.0;
    int Number_of_atoms =  NUM_ATOMS[0];
        
    dirX = large_dist*phiS[rayID]*psiC[rayID];
    dirY = large_dist*phiS[rayID]*psiS[rayID];
    dirZ = large_dist*phiC[rayID];
        
    // setup our quadratic equation
    float a = (dirX*dirX) + (dirY*dirY) + (dirZ*dirZ);
        
    for(atomID = 0; atomID < Number_of_atoms; ++atomID) {
        float atom_X = atomX[atomID];
        float atom_Y = atomY[atomID];
        float atom_Z = atomZ[atomID];
        float atomR = radii[atomID];
        float b = 2.0 * ( (dirX*(-atom_X)) + (dirY*(-atom_Y)) + (dirZ*(-atom_Z)) );
        float c = atom_X*atom_X + atom_Y*atom_Y + atom_Z*atom_Z - (atomR * atomR);
                
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
                    
            if(dist1_sq < distance)
                distance = dist1_sq;
            if(dist2_sq < distance)
                distance = dist2_sq;
        }	
        if(distance<smallest){
            smallest = distance;
        }
    }
    distances[rayID] = smallest;
}
