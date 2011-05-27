/* physics.h
   Physics engine for spring embedding.

   Copyright (c) 2008,2009,2010,2011 Matthias Kramm <kramm@quiss.org> 
 
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#ifndef __physics_h__
#define __physics_h__

#define bool unsigned char

typedef struct _vec3
{
    float x,y,z;
} vec3_t;

vec3_t vec3_diff(vec3_t*v1, vec3_t*v2);
void vec3_zero(vec3_t*v);
vec3_t vec3_sum(vec3_t*v1, vec3_t*v2, float c);
void vec3_add(vec3_t*v1, vec3_t*v2, float c);
float vec3_dot(vec3_t*v1, vec3_t*v2);
float vec3_length(vec3_t*l);

typedef struct _particle
{
    vec3_t position;
    vec3_t velocity;
    vec3_t force;
    float mass;
    bool free;
} particle_t;

void particle_init(particle_t*p, float x, float y, float z);
void particle_move(particle_t*p, float x, float y, float z);
void particle_accelerate(particle_t*p, float x, float y, float z);

typedef struct _force
{
    void (*apply)(struct _force*f);
    bool active;
    particle_t*p1;
    particle_t*p2;
} force_t;

typedef struct _attraction // extends _force
{
    void (*apply)(struct _force*f);
    bool active;
    particle_t*p1;
    particle_t*p2;
    
    float c;
    float len;
} attraction_t;

void attraction_apply(force_t*a);
void attraction_init(attraction_t*a, particle_t*p1, particle_t*p2, float c, float len);
void attraction_apply(force_t*_a);

typedef struct _spring // extends _force
{
    void (*apply)(struct _force*f);
    bool active;
    particle_t*p1;
    particle_t*p2;
    
    float c;
    float len;
    float damp;
} spring_t;

void spring_apply(force_t*s);
void spring_init(spring_t*s, particle_t*p1, particle_t*p2, float c, float len, float damp);
void spring_apply(force_t*_s);

typedef struct _forcelist
{
    force_t*force;
    struct _forcelist*next;
} forcelist_t;

typedef struct _particlelist
{
    particle_t*particle;
    struct _particlelist*next;
} particlelist_t;

typedef struct _physics
{
    forcelist_t*forces;
    particlelist_t*particles;
    float friction;
    float gravity;
    float stepsize;
} physics_t;

void physics_init(physics_t*physics);
void physics_add_force(physics_t*physics,force_t*force);
void physics_add_particle(physics_t*physics,particle_t*particle);
void physics_step(physics_t*physics);
vec3_t*physics_get_forces(physics_t*physics, particle_t**particles, int num_particles);

#endif //__physics_h__
