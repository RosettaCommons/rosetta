__kernel void test_fold(
                         __global half8  const * input,
                         __global float4 const * output
                         ){
  float4 ma = input[get_global_id(0)];

}
