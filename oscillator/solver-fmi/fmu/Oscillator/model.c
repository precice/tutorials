#include "model.h"
#include <float.h> // for DBL_MIN
#include <math.h>  // for fabs()
#include "config.h"

#define V_MIN (0.1)
#define EVENT_EPSILON (1e-10)

void setStartValues(ModelInstance *comp)
{
  M(mass_m)          = 1;
  M(mass_u)          = 0;
  M(mass_v)          = 0;
  M(spring_fixed_c)  = 4 * 3.1416;
  M(spring_middle_c) = 16 * 3.1416;
  M(force_in)        = 0;
  M(force_out)       = 0;
  M(alpha_f)         = 0;
  M(alpha_m)         = 0;
  M(mass_a)          = (M(force_in) - (M(spring_fixed_c) + M(spring_middle_c)) * M(mass_u)) / M(mass_m);
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
  case vr_mass_m:
    value[(*index)++] = M(mass_m);
    return OK;
  case vr_mass_u:
    value[(*index)++] = M(mass_u);
    return OK;
  case vr_mass_v:
    value[(*index)++] = M(mass_v);
    return OK;
  case vr_mass_a:
    value[(*index)++] = M(mass_a);
    return OK;
  case vr_spring_fixed_c:
    value[(*index)++] = M(spring_fixed_c);
    return OK;
  case vr_spring_middle_c:
    value[(*index)++] = M(spring_middle_c);
    return OK;
  case vr_force_in:
    value[(*index)++] = M(force_in);
    return OK;
  case vr_force_out:
    value[(*index)++] = M(force_out);
    return OK;
  case vr_alpha_f:
    value[(*index)++] = M(alpha_f);
    return OK;
  case vr_alpha_m:
    value[(*index)++] = M(alpha_m);
    return OK;
  default:
    logError(comp, "Get Float64 is not allowed for value reference %u.", vr);
    return Error;
  }
}

Status setFloat64(ModelInstance *comp, ValueReference vr, const double *value, size_t *index)
{
  switch (vr) {

  case vr_mass_m:
    M(mass_m) = value[(*index)++];
    return OK;
  case vr_mass_u:
    M(mass_u) = value[(*index)++];
    return OK;
  case vr_mass_v:
    M(mass_v) = value[(*index)++];
    return OK;
  case vr_mass_a:
    M(mass_a) = value[(*index)++];
    return OK;
  case vr_spring_fixed_c:
    M(spring_fixed_c) = value[(*index)++];
    return OK;
  case vr_spring_middle_c:
    M(spring_middle_c) = value[(*index)++];
    return OK;
  case vr_force_in:
    M(force_in) = value[(*index)++];
    return OK;
  case vr_force_out:
    M(force_out) = value[(*index)++];
    return OK;
  case vr_alpha_f:
    M(alpha_f) = value[(*index)++];
    return OK;
  case vr_alpha_m:
    M(alpha_m) = value[(*index)++];
    return OK;
  default:
    logError(comp, "Unexpected value reference: %u.", vr);
    return Error;
  }
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
