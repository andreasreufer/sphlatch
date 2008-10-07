#version 3.6;
global_settings {  assumed_gamma 1.0 }
#default{ finish{ ambient 0.1 diffuse 0.9}}
#include "colors.inc"
camera {location <0.7 , 0.7 ,-1.2>
        look_at  <0.0 , 0.0 , 0.0>}
light_source{<1500,2500,-2500> color White}
