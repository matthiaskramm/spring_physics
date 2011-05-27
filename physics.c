/* physics.c
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

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "physics.h"

vec3_t vec3_diff(vec3_t*v1, vec3_t*v2)
{
    vec3_t result;
    result.x = v1->x - v2->x;
    result.y = v1->y - v2->y;
    result.z = v1->z - v2->z;
    return result;
}

void vec3_zero(vec3_t*v)
{
    v->x = 0;
    v->y = 0;
    v->z = 0;
}

vec3_t vec3_sum(vec3_t*v1, vec3_t*v2, float c)
{
    vec3_t result;
    result.x = v1->x + v2->x * c;
    result.y = v1->y + v2->y * c;
    result.z = v1->z + v2->z * c;
    return result;
}

void vec3_add(vec3_t*v1, vec3_t*v2, float c)
{
    v1->x += v2->x * c;
    v1->y += v2->y * c;
    v1->z += v2->z * c;
}

float vec3_dot(vec3_t*v1, vec3_t*v2)
{
    return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z;
}

float vec3_length(vec3_t*l)
{
    return sqrt(l->x*l->x + l->y*l->y + l->z*l->z);
}

void physics_init(physics_t*physics)
{
    memset(physics, 0, sizeof(physics_t));
    physics->forces = 0;
    physics->stepsize = 0.1;
}
void physics_set_friction(physics_t*physics, float friction)
{
    physics->friction = friction;
}
void physics_set_gravity(physics_t*physics, float gravity)
{
    physics->gravity = gravity;
}
void physics_add_force(physics_t*physics,force_t*force)
{
    forcelist_t*fl = (forcelist_t*)malloc(sizeof(forcelist_t));
    fl->force = force;
    fl->next = physics->forces;
    physics->forces = fl;
}
void physics_remove_force(physics_t*physics,force_t*force)
{
    forcelist_t*it = physics->forces;
    forcelist_t*prev = 0;
    while(it) {
        if(it->force->p1 == force->p1 && it->force->p2 == force->p2 ||
           it->force->p2 == force->p1 && it->force->p1 == force->p2) 
        {
            if(prev) {
                prev->next = it->next;
            } else {
                physics->forces = it->next;
            }
            free(it->force);it->force = 0;
            it->next = 0;
            free(it);
            return;
        }
        prev = it;
        it = it->next;
    }
}
void physics_add_particle(physics_t*physics,particle_t*particle)
{
    particlelist_t*fl = (particlelist_t*)malloc(sizeof(particlelist_t));
    fl->particle = particle;
    fl->next = physics->particles;
    physics->particles = fl;
}
vec3_t*physics_get_forces(physics_t*physics, particle_t**particles, int num_particles)
{
    forcelist_t*fl = physics->forces;
    vec3_t*v = (vec3_t*)malloc(sizeof(vec3_t)*num_particles);
    int t;
    for(t=0;t<num_particles;t++) {
        vec3_zero(&particles[t]->force);
    }
    while(fl) {
        if(fl->force->active) {
            fl->force->apply(fl->force);
        }
        fl = fl->next;
    }
    /* apply friction*/
    if(physics->friction) {
        for(t=0;t<num_particles;t++) {
            vec3_add(&particles[t]->force, &particles[t]->velocity, -physics->friction);
        }
    }
    /* apply gravity */
    if(physics->gravity) {
        vec3_t down = {0,1,0};
        for(t=0;t<num_particles;t++) {
            vec3_add(&particles[t]->force, &down, physics->gravity);
        }
    }
    for(t=0;t<num_particles;t++) {
        v[t] = particles[t]->force;
    }
    return v;
}
void integrate(particle_t**particles, int num_particles, physics_t*physics)
{
    vec3_t*startpos= (vec3_t*)malloc(sizeof(vec3_t)*num_particles);
    vec3_t*startvelocity = (vec3_t*)malloc(sizeof(vec3_t)*num_particles);
    vec3_t*velocity1 = (vec3_t*)malloc(sizeof(vec3_t)*num_particles);
    vec3_t*velocity2 = (vec3_t*)malloc(sizeof(vec3_t)*num_particles);
    vec3_t*velocity3 = (vec3_t*)malloc(sizeof(vec3_t)*num_particles);
    vec3_t*velocity4 = (vec3_t*)malloc(sizeof(vec3_t)*num_particles);
    int t;
    for(t=0;t<num_particles;t++) {
        startpos[t] = particles[t]->position;
        startvelocity[t] = particles[t]->velocity;
        velocity1[t] = particles[t]->velocity;
    }
    vec3_t*force1 = physics_get_forces(physics, particles, num_particles);
    for(t=0;t<num_particles;t++) {
        particles[t]->position = vec3_sum(&startpos[t], &velocity1[t], 0.5*physics->stepsize);
        particles[t]->velocity = vec3_sum(&startvelocity[t], &force1[t], 0.5*physics->stepsize / particles[t]->mass);
        velocity2[t] = particles[t]->velocity;
    }

    vec3_t*force2 = physics_get_forces(physics, particles, num_particles);
    for(t=0;t<num_particles;t++) {
        particles[t]->position = vec3_sum(&startpos[t], &velocity2[t], 0.5*physics->stepsize);
        particles[t]->velocity = vec3_sum(&startvelocity[t], &force2[t], 0.5*physics->stepsize / particles[t]->mass);
        velocity3[t] = particles[t]->velocity;
    }
    vec3_t*force3 = physics_get_forces(physics, particles, num_particles);
    for(t=0;t<num_particles;t++) {
        particles[t]->position = vec3_sum(&startpos[t], &velocity3[t], physics->stepsize);
        particles[t]->velocity = vec3_sum(&startvelocity[t], &force3[t], physics->stepsize / particles[t]->mass);
        velocity4[t] = particles[t]->velocity;
    }
    vec3_t*force4 = physics_get_forces(physics, particles, num_particles);
    for(t=0;t<num_particles;t++) {
        vec3_t*v1 = &velocity1[t];
        vec3_t*v2 = &velocity2[t];
        vec3_t*v3 = &velocity3[t];
        vec3_t*v4 = &velocity4[t];
        particles[t]->position.x = startpos[t].x + (physics->stepsize/6) * (v1->x + 2*v2->x + 2*v3->x + v4->x);
        particles[t]->position.y = startpos[t].y + (physics->stepsize/6) * (v1->y + 2*v2->y + 2*v3->y + v4->y);
        particles[t]->position.z = startpos[t].z + (physics->stepsize/6) * (v1->z + 2*v2->z + 2*v3->z + v4->z);
        vec3_t*f1 = &force1[t];
        vec3_t*f2 = &force2[t];
        vec3_t*f3 = &force3[t];
        vec3_t*f4 = &force4[t];
        particles[t]->velocity.x = startvelocity[t].x + (physics->stepsize/6) * (f1->x + 2*f2->x + 2*f3->x + f4->x);
        particles[t]->velocity.y = startvelocity[t].y + (physics->stepsize/6) * (f1->y + 2*f2->y + 2*f3->y + f4->y);
        particles[t]->velocity.z = startvelocity[t].z + (physics->stepsize/6) * (f1->z + 2*f2->z + 2*f3->z + f4->z);
    }
    free(force1);
    free(force2);
    free(force3);
    free(force4);
    free(velocity1);
    free(velocity2);
    free(velocity3);
    free(velocity4);
    free(startpos);
    free(startvelocity);
}

void physics_step(physics_t*physics)
{
    particlelist_t*pl;
    int num=0;
    pl = physics->particles;
    while(pl) {
        num++;
        pl = pl->next;
    }
    particle_t**list = (particle_t**)malloc(sizeof(particle_t*)*num);

    num = 0;
    pl = physics->particles;
    while(pl) {
        list[num] = pl->particle;
        num++;
        pl = pl->next;
    }
    integrate(list,num,physics);
}

void particle_init(particle_t*p, float x, float y, float z)
{
    memset(p, 0, sizeof(particle_t));
    p->position.x = x;
    p->position.y = y;
    p->position.z = z;
    p->mass = 1.0;
    p->free = 1;
}
void particle_move(particle_t*p, float x, float y, float z)
{
    p->position.x += x;
    p->position.y += y;
    p->position.z += z;
}
void particle_accelerate(particle_t*p, float x, float y, float z)
{
    p->velocity.x += x;
    p->velocity.y += y;
    p->velocity.z += z;
}

void attraction_init(attraction_t*a, particle_t*p1, particle_t*p2, float c, float len)
{
    a->p1 = p1;
    a->p2 = p2;
    a->c = c;
    a->len = len;
    a->active = 1;
    a->apply = attraction_apply;
}
void attraction_apply(force_t*_a)
{
    attraction_t*a = (attraction_t*)_a;

    vec3_t diff = vec3_diff(&a->p1->position, &a->p2->position);
    float dist = vec3_length(&diff);
    if(dist == 0) {
        return;
    }
    /* normalize */
    float ddist = 1.0/dist;
    diff.x *= ddist;
    diff.y *= ddist;
    diff.z *= ddist;

    float force = a->c * a->p1->mass * a->p2->mass;
    if(dist < a->len) {
        ddist = 1/a->len;
    }
    force *= ddist*ddist;

    if(a->p1->free)
        vec3_add(&a->p1->force, &diff, -force);
    if(a->p2->free)
        vec3_add(&a->p2->force, &diff, force);
}

void spring_init(spring_t*s, particle_t*p1, particle_t*p2, float c, float len, float damp)
{
    s->p1 = p1;
    s->p2 = p2;
    s->c = c;
    s->len = len;
    s->active = 1;
    s->damp = damp;
    s->apply = spring_apply;
}

void spring_apply(force_t*_s)
{
    spring_t*s = (spring_t*)_s;

    vec3_t diff = vec3_diff(&s->p1->position, &s->p2->position);
    float dist = vec3_length(&diff);
    float ddist = 0;
    if(dist > 0) {
        /* normalize */
        float ddist = 1.0/dist;
        diff.x *= ddist;
        diff.y *= ddist;
        diff.z *= ddist;
    }

    float springForce = -(dist - s->len) * s->c;
    vec3_t veldiff = vec3_diff(&s->p1->velocity, &s->p2->velocity);
    float dampingForce = -s->damp * vec3_dot(&diff,&veldiff);
    if(dampingForce < -s->damp)
        dampingForce = -s->damp;
    if(dampingForce > s->damp)
        dampingForce = s->damp;
    float r = springForce + dampingForce;

    diff.x *= r;
    diff.y *= r;
    diff.z *= r;

    //printf("diff=(%f %f %f), springForce=%f dampingForce=%f, r=%f\n", diff.x,diff.y,diff.z, springForce, dampingForce, r);

    if(s->p1->free)
        vec3_add(&s->p1->force,&diff,-r);
    if(s->p2->free)
        vec3_add(&s->p2->force,&diff,r);
}

