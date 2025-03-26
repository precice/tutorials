#include "model.h"
#include <float.h> // for DBL_MIN
#include <math.h>  // for fabs()
#include "config.h"

#define V_MIN (0.1)
#define EVENT_EPSILON (1e-10)

void setStartValues(ModelInstance *comp)
{

  M(r)     = 0;
  M(e_ls)  = 0;
  M(kp)    = 0;
  M(ki)    = 0;
  M(kd)    = 0;
  M(P)     = 0;
  M(I)     = 0;
  M(D)     = 0;
  M(y_1)   = 0;
  M(y_2)   = 0;
  M(u_1)   = 0;
  M(u_2)   = 0;
  M(I_max) = 100;

  M(use_implicit_method) = false;
  M(compute_u_1)         = true;
}

Status calculateValues(ModelInstance *comp)
{

  UNUSED(comp);
  return OK;
}

Status getFloat64(ModelInstance *comp, ValueReference vr, double *value, size_t *index)
{
  switch (vr) {
  case vr_time:
    value[(*index)++] = comp->time;
    return OK;
  case vr_r:
    value[(*index)++] = M(r);
    return OK;
  case vr_e:
    value[(*index)++] = M(e);
    return OK;
  case vr_e_ls:
    value[(*index)++] = M(e_ls);
    return OK;
  case vr_kp:
    value[(*index)++] = M(kp);
    return OK;
  case vr_ki:
    value[(*index)++] = M(ki);
    return OK;
  case vr_kd:
    value[(*index)++] = M(kd);
    return OK;
  case vr_P:
    value[(*index)++] = M(P);
    return OK;
  case vr_I:
    value[(*index)++] = M(I);
    return OK;
  case vr_D:
    value[(*index)++] = M(D);
    return OK;
  case vr_y_1:
    value[(*index)++] = M(y_1);
    return OK;
  case vr_y_2:
    value[(*index)++] = M(y_2);
    return OK;
  case vr_u_1:
    value[(*index)++] = M(u_1);
    return OK;
  case vr_u_2:
    value[(*index)++] = M(u_2);
    return OK;
  case vr_I_max:
    value[(*index)++] = M(I_max);
    return OK;
  default:
    logError(comp, "Get Float64 is not allowed for value reference %u.", vr);
    return Error;
  }
}

Status getBoolean(ModelInstance *comp, ValueReference vr, bool *value, size_t *index)
{

  switch (vr) {
  case vr_use_implicit_method:
    value[(*index)++] = M(use_implicit_method);
    break;
  case vr_compute_u_1:
    value[(*index)++] = M(compute_u_1);
    break;
  default:
    logError(comp, "Get Boolean is not allowed for value reference %u.", vr);
    return Error;
  }

  return OK;
}

Status setFloat64(ModelInstance *comp, ValueReference vr, const double *value, size_t *index)
{
  switch (vr) {

  case vr_r:
    M(r) = value[(*index)++];
    return OK;
  case vr_kp:
    M(kp) = value[(*index)++];
    return OK;
  case vr_ki:
    M(ki) = value[(*index)++];
    return OK;
  case vr_kd:
    M(kd) = value[(*index)++];
    return OK;
  case vr_y_1:
    M(y_1) = value[(*index)++];
    return OK;
  case vr_y_2:
    M(y_2) = value[(*index)++];
    return OK;
  case vr_I_max:
    M(I_max) = value[(*index)++];
    return OK;
  default:
    logError(comp, "Unexpected value reference: %u.", vr);
    return Error;
  }
}

Status setBoolean(ModelInstance *comp, ValueReference vr, const bool *value, size_t *index)
{

  switch (vr) {
  case vr_use_implicit_method:
    M(use_implicit_method) = value[(*index)++];
    break;
  case vr_compute_u_1:
    M(compute_u_1) = value[(*index)++];
    break;
  default:
    logError(comp, "Set Boolean is not allowed for value reference %u.", vr);
    return Error;
  }

  comp->isDirtyValues = true;

  return OK;
}

Status getOutputDerivative(ModelInstance *comp, ValueReference valueReference, int order, double *value)
{

  if (order != 1) {
    logError(comp, "The output derivative order %d for value reference %u is not available.", order, valueReference);
    return Error;
  }

  switch (valueReference) {

  default:
    logError(comp, "The output derivative for value reference %u is not available.", valueReference);
    return Error;
  }

  UNUSED(value);
}

void eventUpdate(ModelInstance *comp)
{

  if (false) {

  } else {
    comp->valuesOfContinuousStatesChanged = false;
  }

  comp->nominalsOfContinuousStatesChanged = false;
  comp->terminateSimulation               = false;
  comp->nextEventTimeDefined              = false;
}

void getContinuousStates(ModelInstance *comp, double x[], size_t nx)
{
  UNUSED(comp);
  UNUSED(nx);
  UNUSED(x);
}

void setContinuousStates(ModelInstance *comp, const double x[], size_t nx)
{
  UNUSED(comp);
  UNUSED(nx);
  UNUSED(x);
}

void getDerivatives(ModelInstance *comp, double dx[], size_t nx)
{
  UNUSED(comp);
  UNUSED(nx);
  UNUSED(dx);
}

void getEventIndicators(ModelInstance *comp, double z[], size_t nz)
{

  UNUSED(nz);
  UNUSED(comp);
  UNUSED(z);
}
