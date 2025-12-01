#ifndef STRUCTURES
#define STRUCTURES

struct vec3
{
    float x, y, z;
};

struct particle
{
    struct vec3 pos;
    struct vec3 vel;
    float charge;
    float mass;
};

#endif