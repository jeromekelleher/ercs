/*
** Copyright (C) 2012 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
**  
** This file is part of ercs.
** 
** ercs is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** ercs is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with ercs.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <Python.h>

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

#include <float.h>
#include "lib/ercs.h"

#define DISC_EVENT_CLASS 0
#define GAUSSIAN_EVENT_CLASS 1

#define MODULE_DOC \
"Simple Python interface to the ``ercs`` C library. Supports Python 2 and 3."

static PyObject *ErcsInputError;
static PyObject *ErcsLibraryError;

static int
pyercs_check_parameters(ercs_t *sim)
{
    int v;
    if (sim->torus_diameter <= 0.0) {
        PyErr_SetString(ErcsInputError, "torus_diameter must be > 0.0"); 
        goto out;
    }
    if (sim->max_time < 0.0) {
        PyErr_SetString(ErcsInputError, "torus_diameter must be >= 0.0"); 
        goto out;
    } 
    if (sim->num_parents <= 0) {
        PyErr_SetString(ErcsInputError, "num_parents must be > 0"); 
        goto out;
    }
    v = sim->kdtree_bucket_size;
    if (v <= 0) {
        PyErr_SetString(ErcsInputError, "kdtree_bucket_size must be > 0"); 
        goto out;
    }
    if ((v & (v - 1)) != 0) {
        PyErr_SetString(ErcsInputError, 
                "kdtree_bucket_size must be a power of 2"); 
        goto out;
    }
    if (sim->max_kdtree_insertions < 0) {
        PyErr_SetString(ErcsInputError, "max_kdtree_insertions must be >= 0");
        goto out;
    }
    if (sim->max_lineages < 2) {
        PyErr_SetString(ErcsInputError, "max_lineages must be >= 2");
        goto out;
    }
out:   
    return PyErr_Occurred() == NULL;
}

static int
pyercs_parse_sample(PyObject *py_sample, ercs_t *sim)
{
    int size;
    int j, k;
    double v;
    PyObject *item, *value;    
    int n = PyList_Size(py_sample);
    if (n == 0) {
        PyErr_SetString(ErcsInputError, "Empty sample"); 
        goto out;
    }
    sim->sample_size = n; 
    sim->sample = xcalloc(2 * n, sizeof(double));
    for (j = 0; j < n; j++) {
        item = PyList_GetItem(py_sample, j);
        size = 0;
        if (!PyTuple_Check(item)) {
            PyErr_SetString(ErcsInputError, "Samples must be 2-tuples"); 
            goto out;
        } else {
            size = PyTuple_Size(item);
            if (size != 2) {
                PyErr_SetString(ErcsInputError, "Dimension != 2 not supported"); 
                goto out;
            }
            for (k = 0; k < 2; k++) {
                value = PyTuple_GetItem(item, k);
                if (!PyNumber_Check(value)) {
                    PyErr_SetString(ErcsInputError, 
                            "Locations must be numeric");
                    goto out;
                }
                v = PyFloat_AsDouble(value);
                sim->sample[j * 2 + k] = v;
                if (v < 0.0 || v >= sim->torus_diameter) {
                    PyErr_SetString(ErcsInputError, 
                            "sample location: must have 0 <= v < L"); 
                    goto out;
                }
            }
        }
    }
out:   
    return PyErr_Occurred() == NULL;
}

static int  
pyercs_parse_recombination(PyObject *py_recombination, ercs_t *sim) 
{
    int j;
    int size;
    double v;
    PyObject *item;
    size = PyList_Size(py_recombination);
    sim->num_loci = size + 1; 
    sim->recombination_probabilities = xcalloc(sim->num_loci, sizeof(double));
    for (j = 0; j < size; j++) {
        item = PyList_GetItem(py_recombination, j);
        if (!PyNumber_Check(item)) {
            PyErr_SetString(ErcsInputError, "Not a number"); 
            goto out;
        }
        v = PyFloat_AsDouble(item);
        sim->recombination_probabilities[j] = v;
        if (v < 0.0 || v > 1) {
            PyErr_SetString(ErcsInputError, "invalid rho: not a probability"); 
            goto out;
        }
    }
out: 
    return PyErr_Occurred() == NULL;
}

/*
 * Retrieves a number value with the specified key from the specified 
 * dictionary.
 */
static PyObject *
get_dict_number(PyObject *dict, const char *key_str)
{
    PyObject *ret = NULL;
    PyObject *value;
    PyObject *key = Py_BuildValue("s", key_str);
    if (!PyDict_Contains(dict, key)) {
        PyErr_Format(ErcsInputError, "'%s' not specified", key_str); 
        goto out;
    }
    value = PyDict_GetItem(dict, key);
    if (!PyNumber_Check(value)) {
        PyErr_Format(ErcsInputError, "'%s' is not number", key_str); 
        goto out;
    }
    ret = value;
out:
    Py_DECREF(key);
    return ret;
}


static int
pyercs_parse_events(PyObject *py_events, ercs_t *sim) 
{
    int j, size;
    long type;
    double rate, u, r, alpha, theta;
    PyObject *item, *value;
    size = PyList_Size(py_events);
    if (size == 0) {
        PyErr_SetString(ErcsInputError, "must have > 0 events"); 
        goto out;
    }
    sim->num_event_classes = size; 
    sim->event_classes = xmalloc(size * sizeof(event_class_t));
    for (j = 0; j < size; j++) {
        item = PyList_GetItem(py_events, j);
        if (!PyDict_Check(item)) {
            PyErr_SetString(ErcsInputError, "not a dictionary"); 
            goto out;
        }
        value = get_dict_number(item, "rate");
        if (value == NULL) {
            goto out;
        }
        rate = PyFloat_AsDouble(value);
        value = get_dict_number(item, "type");
        if (value == NULL) {
            goto out;
        }
        type = PyLong_AsLong(value);
        if (type == DISC_EVENT_CLASS) {
            value = get_dict_number(item, "r");
            if (value == NULL) {
                goto out;
            }
            r = PyFloat_AsDouble(value);
            value = get_dict_number(item, "u");
            if (value == NULL) {
                goto out;
            }
            u = PyFloat_AsDouble(value);
            alloc_disc_event_class(&sim->event_classes[j], rate, r, u);
        } else  if (type == GAUSSIAN_EVENT_CLASS) {
            value = get_dict_number(item, "theta");
            if (value == NULL) {
                goto out;
            }
            theta = PyFloat_AsDouble(value);
            value = get_dict_number(item, "alpha");
            if (value == NULL) {
                goto out;
            }
            alpha = PyFloat_AsDouble(value);
            value = get_dict_number(item, "u0");
            if (value == NULL) {
                goto out;
            }
            u = PyFloat_AsDouble(value);
            alloc_gaussian_event_class(&sim->event_classes[j], rate, theta,
                    alpha, u);

        } else {
            PyErr_SetString(ErcsInputError, "Unknown event type");
            goto out;
        }
    }
out:
    return PyErr_Occurred() == NULL;
}


#define SIMULATE_DOC \
"Allocates an ercs object from the C library, calls the simulate function "\
"and returns the resulting genealogy. All arguments must be specified and "\
"be in the correct order.\n\n"\
"   :param random_seed: The seed used for the GSL random generator\n"\
"   :type random_seed: long integer\n"\
"   :param torus_diameter: The diameter of the torus\n"\
"   :type torus_diameter: double\n"\
"   :param num_parents: The number of parents in each event\n"\
"   :type num_parents: integer\n"\
"   :param sample: The (zero indexed) sample of 2D locations\n"\
"   :type sample: list of numeric (x, y) tuples\n"\
"   :param event_classes: The list of event classes and their rates\n"\
"   :type event_classes: list of dictionaries\n"\
"   :param recombination_probabilities: probability of recombination between"\
" adjacent loci\n"\
"   :type recombination_probabilities: list of doubles\n"\
"   :param kdtree_bucket_size: The number of points in a kdtree bucket\n"\
"   :type kdtree_bucket_size: integer\n"\
"   :param max_kdtree_insertions: The maximum number of insertions before "\
" the kdtree is rebuilt; if 0, the kdtree is never rebuilt\n"\
"   :type max_kdtree_insertions: integer\n"\
"   :param max_lineages: The maximum number of extant lineages\n"\
"   :type max_lineages: integer\n"\
"   :param max_time: the maximum time we simulate back into the past;"\
" if 0.0, simulate until coalescence\n"\
"   :type max_time: double\n"\
"   :param ancestry_algorithm: currently unused\n"\
"   :type ancestry_algorithm: int\n"\
"   :return: the simulated history of the sample, (pi, tau).\n"\
"   :rtype: a tuple (pi, tau); pi is a list of lists of integers, and "\
"        tau is a list of lists of doubles\n"\
"   :raises: InputError when the input is not correctly formed.\n"\
"   :raises: LibraryError when the C library encounters an error\n"\

static PyObject *
pyercs_simulate(PyObject *self, PyObject *args)
{
    int err, ercs_ret, not_done;
    PyObject *ret = NULL;
    PyObject *pi = NULL;
    PyObject *tau = NULL;
    int j, l, n, ancestry_algorithm;
    PyObject *py_sample, *py_recombination, *py_events, *pi_locus, *tau_locus;
    ercs_t *sim = xcalloc(1, sizeof(ercs_t));
    sim->sample = NULL;
    sim->event_classes = NULL;
    sim->recombination_probabilities = NULL;
    if (!PyArg_ParseTuple(args, "ldiO!O!O!iiidi", 
                &sim->random_seed,
                &sim->torus_diameter,
                &sim->num_parents,
                &PyList_Type, &py_sample,
                &PyList_Type, &py_events,
                &PyList_Type, &py_recombination,
                &sim->kdtree_bucket_size,
                &sim->max_kdtree_insertions,
                &sim->max_lineages,
                &sim->max_time,
                &ancestry_algorithm /* unused */)) {
        goto out; 
    }
    if (!pyercs_check_parameters(sim)) {
        goto out;
    }
    if (!pyercs_parse_sample(py_sample, sim)) {
        goto out;
    }
    if (!pyercs_parse_events(py_events, sim)) {
        goto out;
    }
    if (!pyercs_parse_recombination(py_recombination, sim)) {
        goto out;
    }
    /* We now have valid input values */
    if (sim->max_time == 0.0) {
        sim->max_time = DBL_MAX;
    }
    if (sim->max_kdtree_insertions == 0) {
        sim->max_kdtree_insertions = INT_MAX;
    }
    ercs_ret = ercs_initialise(sim);
    ERCS_ERROR_CHECK(ercs_ret, cleanup);
    not_done = 1; 
    while (not_done) {
        ercs_ret = ercs_simulate(sim, 1u<<20);
        ERCS_ERROR_CHECK(ercs_ret, cleanup);
        not_done = ercs_ret == ERCS_SIM_NOT_DONE;
        if (PyErr_CheckSignals() < 0) {
            goto out;
        }
    }
    /* output */
    pi = PyList_New(sim->num_loci);
    if (pi == NULL) {
        goto cleanup;
    }
    tau = PyList_New(sim->num_loci);
    if (tau == NULL) {
        goto cleanup;
    }
    n = 2 * sim->sample_size;
    for (l = 0; l < sim->num_loci; l++) {
        pi_locus = PyList_New(n);
        if (pi_locus == NULL) {
            goto cleanup;
        }
        err = PyList_SetItem(pi, l, pi_locus);
        if (err < 0) {
            goto cleanup;
        }
        tau_locus = PyList_New(n);
        if (tau_locus == NULL) {
            goto cleanup;
        }
        err = PyList_SetItem(tau, l, tau_locus);
        if (err < 0) {
            goto cleanup;
        }
        for (j = 0; j < n; j++) {
            err = PyList_SetItem(pi_locus, j, PyLong_FromLong(sim->pi[l][j])); 
            if (err < 0) {
                goto cleanup;
            }
            err = PyList_SetItem(tau_locus, j, 
                    PyFloat_FromDouble(sim->tau[l][j])); 
            if (err < 0) {
                goto cleanup;
            }
        }
    }
    ret = Py_BuildValue("(O, O)", pi, tau);
cleanup:
    if (pi != NULL) {
        Py_DECREF(pi);
    }
    if (tau != NULL) {
        Py_DECREF(tau);
    }
    ercs_free(sim);
    if (ercs_ret < 0) {
        PyErr_SetString(ErcsLibraryError, ercs_error_str(ercs_ret)); 
    }
out:
    if (sim->sample != NULL) {
        free(sim->sample);
    }
    if (sim->event_classes != NULL) {
        free(sim->event_classes);
    }
    if (sim->recombination_probabilities != NULL) {
        free(sim->recombination_probabilities);
    }
    free(sim);
    return ret; 
}

static PyMethodDef ErcsMethods[] = {
    {"simulate",  pyercs_simulate, METH_VARARGS, SIMULATE_DOC},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


/* Initialisation code supports Python 2.x and 3.x. The framework uses the 
 * recommended structure from http://docs.python.org/howto/cporting.html. 
 * I've ignored the point about storing state in globals, as the examples 
 * from the Python documentation still use this idiom. 
 */

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef ercsmodule = {
    PyModuleDef_HEAD_INIT,
    "_ercs",   /* name of module */
    MODULE_DOC, /* module documentation, may be NULL */
    -1,    
    ErcsMethods
};

#define INITERROR return NULL

PyObject * 
PyInit__ercs(void)

#else
#define INITERROR return

void
init_ercs(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&ercsmodule);
#else
    PyObject *module = Py_InitModule3("_ercs", ErcsMethods, MODULE_DOC);
#endif
    if (module == NULL) {
        INITERROR;
    }
    ErcsInputError = PyErr_NewException("_ercs.InputError", NULL, NULL);
    Py_INCREF(ErcsInputError);
    PyModule_AddObject(module, "InputError", ErcsInputError);
    ErcsLibraryError = PyErr_NewException("_ercs.LibraryError", NULL, NULL);
    Py_INCREF(ErcsLibraryError);
    PyModule_AddObject(module, "LibraryError", ErcsLibraryError);
    
    PyModule_AddIntConstant(module, "DISC_EVENT_CLASS", DISC_EVENT_CLASS);
    PyModule_AddIntConstant(module, "GAUSSIAN_EVENT_CLASS", 
            GAUSSIAN_EVENT_CLASS);

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}






