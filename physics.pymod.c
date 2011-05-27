/* physics.pymod.c
   Physics engine for spring embedding- Python interface.

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

#include <stdio.h>
#undef _POSIX_C_SOURCE
#include <Python.h>
#include "physics.h"

static char* strf(char*format, ...)
{
    char buf[1024];
    int l;
    va_list arglist;
    va_start(arglist, format);
    vsprintf(buf, format, arglist);
    va_end(arglist);
    return strdup(buf);
}

extern PyTypeObject ParticleClass;
extern PyTypeObject ParticleSystemClass;

typedef struct {
    PyObject_HEAD
    particle_t particle;
} ParticleObject;

typedef struct {
    PyObject_HEAD
    physics_t system;
} ParticleSystemObject;

#define PY_ERROR(s,args...) (PyErr_SetString(PyExc_Exception, strf(s, ## args)),NULL)
#define PY_NONE Py_BuildValue("s", 0)

//---------------------------------------------------------------------

static PyMethodDef ParticleMethods[] = {
    {0,0,0,0}
};
static int particle_setattr(PyObject * _self, char* a, PyObject* o)
{
    ParticleObject* self = (ParticleObject*)_self;
    if(!strcmp(a, "mass")) {
	return !PyArg_Parse(o, "f", &self->particle.mass);
    } else if(!strcmp(a, "x")) { 
	return !PyArg_Parse(o, "f", &self->particle.position.x);
    } else if(!strcmp(a, "y")) { 
	return !PyArg_Parse(o, "f", &self->particle.position.y);
    } else if(!strcmp(a, "z")) { 
	return !PyArg_Parse(o, "f", &self->particle.position.z);
    }
    return 1;
}
static PyObject* particle_getattr(PyObject * _self, char* a)
{
    ParticleObject* self = (ParticleObject*)_self;
    if(!strcmp(a, "x")) return PyFloat_FromDouble(self->particle.position.x);
    else if(!strcmp(a, "y")) return PyFloat_FromDouble(self->particle.position.y);
    else if(!strcmp(a, "z")) return PyFloat_FromDouble(self->particle.position.z);
    else if(!strcmp(a, "pos")) return Py_BuildValue("(fff)",self->particle.position.x,self->particle.position.y,self->particle.position.z);
    else if(!strcmp(a, "force")) return Py_BuildValue("(fff)",self->particle.force.x,self->particle.force.y,self->particle.force.z);
    else if(!strcmp(a, "velocity")) return Py_BuildValue("(fff)",self->particle.velocity.x,self->particle.velocity.y,self->particle.velocity.z);
    return Py_FindMethod(ParticleMethods, _self, a);
}

static void particle_dealloc(PyObject* self)
{
    PyObject_Del(self);
}
//---------------------------------------------------------------------

PyObject* f_ParticleSystem(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    static char *kwlist[] = {"stepsize", NULL};
    float stepsize = 0;
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "|f", kwlist, &stepsize))
        return NULL;
    ParticleSystemObject*system= PyObject_New(ParticleSystemObject, &ParticleSystemClass);
    physics_init(&system->system);
    if(stepsize) {
        system->system.stepsize = stepsize;
    }
    return (PyObject*)system;
}
static int particlesystem_setattr(PyObject * _self, char* a, PyObject* o)
{
    ParticleSystemObject* self = (ParticleSystemObject*)_self;
    if(!strcmp(a, "friction")) {
	return !PyArg_Parse(o, "f", &self->system.friction);
    } else if(!strcmp(a, "gravity")) {
	return !PyArg_Parse(o, "f", &self->system.gravity);
    }
    return 1;
}
void particlesystem_dealloc(PyObject* self)
{
    PyObject_Del(self);
}
PyObject* particlesystem_createParticle(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ParticleSystemObject* self = (ParticleSystemObject*)_self;
    static char *kwlist[] = {"x","y","z","free",NULL};
    float x,y,z;
    char free=1;
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "fff|i", kwlist, &x, &y, &z, &free))
        return NULL;
    ParticleObject*particle= PyObject_New(ParticleObject, &ParticleClass);
    particle_init(&particle->particle, x,y,z);
    particle->particle.free = free;
    physics_add_particle(&self->system,&particle->particle);
    return (PyObject*)particle;
}

PyObject* particlesystem_step(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ParticleSystemObject* self = (ParticleSystemObject*)_self;
    physics_step(&self->system);
    return PY_NONE;
}

PyObject* particlesystem_addSpring(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ParticleSystemObject* self = (ParticleSystemObject*)_self;
    static char *kwlist[] = {"p1","p2","c", "len", "damp", NULL};
    float c=0,len=0,damp=0;
    ParticleObject*p1=0,*p2=0;
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!fff", kwlist, &ParticleClass, &p1, &ParticleClass, &p2, &c,&len,&damp))
        return NULL;

    spring_t* s = malloc(sizeof(spring_t));
    spring_init(s,&p1->particle,&p2->particle,c,len,damp);
    physics_add_force(&self->system,(force_t*)s);
    
    return PY_NONE;
}

PyObject* particlesystem_addAttraction(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ParticleSystemObject* self = (ParticleSystemObject*)_self;
    static char *kwlist[] = {"p1","p2","c", "len", NULL};
    float c=0,len=0;
    ParticleObject*p1=0,*p2=0;
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!ff", kwlist, &ParticleClass, &p1, &ParticleClass, &p2, &c,&len))
        return NULL;

    attraction_t* a = malloc(sizeof(attraction_t));
    attraction_init(a,&p1->particle,&p2->particle,c,len);
    physics_add_force(&self->system,(force_t*)a);
    return PY_NONE;
}

PyObject* particlesystem_remove(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ParticleSystemObject* self = (ParticleSystemObject*)_self;
    static char *kwlist[] = {"p1","p2",NULL};
    ParticleObject*p1=0,*p2=0;
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!", kwlist, &ParticleClass, &p1, &ParticleClass, &p2))
        return NULL;

    force_t f;
    f.p1 = &p1->particle;
    f.p2 = &p2->particle;
    physics_remove_force(&self->system,&f);
    return PY_NONE;
}

static PyMethodDef ParticleSystemMethods[] = 
{
    {"addSpring", (PyCFunction)particlesystem_addSpring, METH_KEYWORDS, "Adds a spring"},
    {"addAttraction", (PyCFunction)particlesystem_addAttraction, METH_KEYWORDS, "Adds an attraction"},
    {"createParticle", (PyCFunction)particlesystem_createParticle, METH_KEYWORDS, "Adds a new particle"},
    {"remove", (PyCFunction)particlesystem_remove, METH_KEYWORDS, "Removes a force"},
    {"step", (PyCFunction)particlesystem_step, METH_KEYWORDS, "Performs a time step"},
    {0,0,0,0}
};

static PyObject* particlesystem_getattr(PyObject * _self, char* a)
{
    ParticleSystemObject* self = (ParticleSystemObject*)_self;
    return Py_FindMethod(ParticleSystemMethods, _self, a);
}
//---------------------------------------------------------------------

PyTypeObject ParticleClass = 
{
    PyObject_HEAD_INIT(NULL)
    0,
    tp_name: "Particle",
    tp_basicsize: sizeof(ParticleObject),
    tp_itemsize: 0,
    tp_dealloc: particle_dealloc,
    tp_print: 0,
    tp_getattr: particle_getattr,
    tp_setattr: particle_setattr,
};

PyTypeObject ParticleSystemClass = 
{
    PyObject_HEAD_INIT(NULL)
    0,
    tp_name: "ParticleSystem",
    tp_basicsize: sizeof(ParticleSystemObject),
    tp_itemsize: 0,
    tp_dealloc: particlesystem_dealloc,
    tp_print: 0,
    tp_getattr: particlesystem_getattr,
    tp_setattr: particlesystem_setattr,
};

static PyMethodDef physics_methods[] =
{
    /* module functions */
    {"ParticleSystem", (PyCFunction)f_ParticleSystem, METH_KEYWORDS, ""},
    {0, 0, 0, 0}
};

void initphysics(void)
{
    PyObject*module = Py_InitModule("physics", physics_methods);
    ParticleSystemClass.ob_type = &PyType_Type;
    ParticleClass.ob_type = &PyType_Type;
}
