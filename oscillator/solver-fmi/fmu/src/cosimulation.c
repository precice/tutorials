#include "cosimulation.h"
#include <float.h> // for DBL_EPSILON
#include <math.h>  // for fabs()
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h> // for calloc(), free()
#include <string.h>
#include "config.h"

#if FMI_VERSION == 3
#include "fmi3Functions.h"
#endif

#ifdef _MSC_VER
#define strdup _strdup
#endif

ModelInstance *createModelInstance(
    loggerType             cbLogger,
    intermediateUpdateType intermediateUpdate,
    void                  *componentEnvironment,
    const char            *instanceName,
    const char            *instantiationToken,
    const char            *resourceLocation,
    bool                   loggingOn,
    InterfaceType          interfaceType)
{

  ModelInstance *comp = NULL;

  if (!instanceName || strlen(instanceName) == 0) {
    if (cbLogger) {
#if FMI_VERSION < 3
      cbLogger(componentEnvironment, "?", Error, "error", "Missing instance name.");
#else
      cbLogger(componentEnvironment, Error, "error", "Missing instance name.");
#endif
    }
    return NULL;
  }

  if (!instantiationToken || strlen(instantiationToken) == 0) {
    if (cbLogger) {
#if FMI_VERSION < 3
      cbLogger(componentEnvironment, instanceName, Error, "error", "Missing GUID.");
#else
      cbLogger(componentEnvironment, Error, "error", "Missing instantiationToken.");
#endif
    }
    return NULL;
  }

  if (strcmp(instantiationToken, INSTANTIATION_TOKEN)) {
    if (cbLogger) {
#if FMI_VERSION < 3
      cbLogger(componentEnvironment, instanceName, Error, "error", "Wrong GUID.");
#else
      cbLogger(componentEnvironment, Error, "error", "Wrong instantiationToken.");
#endif
    }
    return NULL;
  }

  comp = (ModelInstance *) calloc(1, sizeof(ModelInstance));

  if (comp) {
    comp->componentEnvironment = componentEnvironment;
    comp->logger               = cbLogger;
    comp->intermediateUpdate   = intermediateUpdate;
    comp->lockPreemtion        = NULL;
    comp->unlockPreemtion      = NULL;
    comp->instanceName         = strdup(instanceName);
    comp->resourceLocation     = resourceLocation ? strdup(resourceLocation) : NULL;
    comp->status               = OK;
    comp->logEvents            = loggingOn;
    comp->logErrors            = true; // always log errors
    comp->nSteps               = 0;
    comp->solverStepSize       = 1e-3;
    comp->earlyReturnAllowed   = false;
    comp->eventModeUsed        = false;
  }

  if (!comp || !comp->instanceName) {
    logError(comp, "Out of memory.");
    return NULL;
  }

  comp->time = 0.0; // overwrite in fmi*SetupExperiment, fmi*SetTime
  comp->type = interfaceType;

  comp->state = Instantiated;

  comp->newDiscreteStatesNeeded           = false;
  comp->terminateSimulation               = false;
  comp->nominalsOfContinuousStatesChanged = false;
  comp->valuesOfContinuousStatesChanged   = false;
  comp->nextEventTimeDefined              = false;
  comp->nextEventTime                     = 0;

  setStartValues(comp);

  comp->isDirtyValues = true;

  return comp;
}

void freeModelInstance(ModelInstance *comp)
{
  free((void *) comp->instanceName);
  free(comp);
}

void reset(ModelInstance *comp)
{
  comp->state          = Instantiated;
  comp->startTime      = 0.0;
  comp->time           = 0.0;
  comp->nSteps         = 0;
  comp->solverStepSize = 1e-3;
  comp->status         = OK;
  setStartValues(comp);
  comp->isDirtyValues = true;
}

bool invalidNumber(ModelInstance *comp, const char *f, const char *arg, size_t actual, size_t expected)
{

  if (actual != expected) {
    comp->state = modelError;
    logError(comp, "%s: Invalid argument %s = %d. Expected %d.", f, arg, actual, expected);
    return true;
  }

  return false;
}

bool invalidState(ModelInstance *comp, const char *f, int statesExpected)
{

  UNUSED(f);
  UNUSED(statesExpected);

  if (!comp) {
    return true;
  }

  // TODO: add missing states and check state
  return false;

  //    if (!(comp->state & statesExpected)) {
  //        comp->state = modelError;
  //        logError(comp, "%s: Illegal call sequence.", f);
  //        return true;
  //    }
  //
  //    return false;
}

bool nullPointer(ModelInstance *comp, const char *f, const char *arg, const void *p)
{

  if (!p) {
    comp->state = modelError;
    logError(comp, "%s: Invalid argument %s = NULL.", f, arg);
    return true;
  }

  return false;
}

Status setDebugLogging(ModelInstance *comp, bool loggingOn, size_t nCategories, const char *const categories[])
{

  if (loggingOn) {
    for (size_t i = 0; i < nCategories; i++) {
      if (categories[i] == NULL) {
        logError(comp, "Log category[%d] must not be NULL", i);
        return Error;
      } else if (strcmp(categories[i], "logEvents") == 0) {
        comp->logEvents = true;
      } else if (strcmp(categories[i], "logStatusError") == 0) {
        comp->logErrors = true;
      } else {
        logError(comp, "Log category[%d] must be one of logEvents or logStatusError but was %s", i, categories[i]);
        return Error;
      }
    }
  } else {
    // disable logging
    comp->logEvents = false;
    comp->logErrors = false;
  }

  return OK;
}

static void logMessage(ModelInstance *comp, int status, const char *category, const char *message, va_list args)
{

  if (!comp->logger) {
    return;
  }

  va_list args1;
  size_t  len = 0;
  char   *buf = "";

  va_copy(args1, args);
  len = vsnprintf(buf, len, message, args1);
  va_end(args1);

  va_copy(args1, args);
  buf = (char *) calloc(len + 1, sizeof(char));
  vsnprintf(buf, len + 1, message, args);
  va_end(args1);

  // no need to distinguish between FMI versions since we're not using variadic arguments
#if FMI_VERSION < 3
  comp->logger(comp->componentEnvironment, comp->instanceName, status, category, buf);
#else
  comp->logger(comp->componentEnvironment, status, category, buf);
#endif

  free(buf);
}

void logEvent(ModelInstance *comp, const char *message, ...)
{

  if (!comp || !comp->logEvents)
    return;

  va_list args;
  va_start(args, message);
  logMessage(comp, OK, "logEvents", message, args);
  va_end(args);
}

void logError(ModelInstance *comp, const char *message, ...)
{

  if (!comp || !comp->logErrors)
    return;

  va_list args;
  va_start(args, message);
  logMessage(comp, Error, "logStatusError", message, args);
  va_end(args);
}

// default implementations
#if NZ < 1
void getEventIndicators(ModelInstance *comp, double z[], size_t nz)
{
  UNUSED(comp);
  UNUSED(z);
  UNUSED(nz);
  // do nothing
}
#endif

#define GET_NOT_ALLOWED(t)                           \
  do {                                               \
    UNUSED(vr);                                      \
    UNUSED(value);                                   \
    UNUSED(index);                                   \
    logError(comp, "Getting " t " is not allowed."); \
    return Error;                                    \
  } while (false)

#ifndef GET_FLOAT32
Status getFloat32(ModelInstance *comp, ValueReference vr, float value[], size_t *index)
{
  GET_NOT_ALLOWED("Float32");
}
#endif

#ifndef GET_INT8
Status getInt8(ModelInstance *comp, ValueReference vr, int8_t value[], size_t *index)
{
  GET_NOT_ALLOWED("Int8");
}
#endif

#ifndef GET_UINT8
Status getUInt8(ModelInstance *comp, ValueReference vr, uint8_t value[], size_t *index)
{
  GET_NOT_ALLOWED("UInt8");
}
#endif

#ifndef GET_INT16
Status getInt16(ModelInstance *comp, ValueReference vr, int16_t value[], size_t *index)
{
  GET_NOT_ALLOWED("Int16");
}
#endif

#ifndef GET_UINT16
Status getUInt16(ModelInstance *comp, ValueReference vr, uint16_t value[], size_t *index)
{
  GET_NOT_ALLOWED("UInt16");
}
#endif

#ifndef GET_INT32
Status getInt32(ModelInstance *comp, ValueReference vr, int32_t value[], size_t *index)
{
  GET_NOT_ALLOWED("Int32");
}
#endif

#ifndef GET_UINT32
Status getUInt32(ModelInstance *comp, ValueReference vr, uint32_t value[], size_t *index)
{
  GET_NOT_ALLOWED("UInt32");
}
#endif

#ifndef GET_INT64
Status getInt64(ModelInstance *comp, ValueReference vr, int64_t value[], size_t *index)
{
  GET_NOT_ALLOWED("Int64");
}
#endif

#ifndef GET_UINT64
Status getUInt64(ModelInstance *comp, ValueReference vr, uint64_t value[], size_t *index)
{
  GET_NOT_ALLOWED("UInt64");
}
#endif

#ifndef GET_BOOLEAN
Status getBoolean(ModelInstance *comp, ValueReference vr, bool value[], size_t *index)
{
  GET_NOT_ALLOWED("Boolean");
}
#endif

#ifndef GET_STRING
Status getString(ModelInstance *comp, ValueReference vr, const char *value[], size_t *index)
{
  GET_NOT_ALLOWED("String");
}
#endif

#ifndef GET_BINARY
Status getBinary(ModelInstance *comp, ValueReference vr, size_t size[], const char *value[], size_t *index)
{
  UNUSED(size);
  GET_NOT_ALLOWED("Binary");
}
#endif

#define SET_NOT_ALLOWED(t)                           \
  do {                                               \
    UNUSED(vr);                                      \
    UNUSED(value);                                   \
    UNUSED(index);                                   \
    logError(comp, "Setting " t " is not allowed."); \
    return Error;                                    \
  } while (false)

#ifndef SET_FLOAT32
Status setFloat32(ModelInstance *comp, ValueReference vr, const float value[], size_t *index)
{
  SET_NOT_ALLOWED("Float32");
}
#endif

#ifndef SET_FLOAT64
Status setFloat64(ModelInstance *comp, ValueReference vr, const double value[], size_t *index)
{
  SET_NOT_ALLOWED("Float64");
}
#endif

#ifndef SET_INT8
Status setInt8(ModelInstance *comp, ValueReference vr, const int8_t value[], size_t *index)
{
  SET_NOT_ALLOWED("Int8");
}
#endif

#ifndef SET_UINT8
Status setUInt8(ModelInstance *comp, ValueReference vr, const uint8_t value[], size_t *index)
{
  SET_NOT_ALLOWED("UInt8");
}
#endif

#ifndef SET_INT16
Status setInt16(ModelInstance *comp, ValueReference vr, const int16_t value[], size_t *index)
{
  SET_NOT_ALLOWED("Int16");
}
#endif

#ifndef SET_UINT16
Status setUInt16(ModelInstance *comp, ValueReference vr, const uint16_t value[], size_t *index)
{
  SET_NOT_ALLOWED("UInt16");
}
#endif

#ifndef SET_INT32
Status setInt32(ModelInstance *comp, ValueReference vr, const int32_t value[], size_t *index)
{
  SET_NOT_ALLOWED("Int32");
}
#endif

#ifndef SET_UINT32
Status setUInt32(ModelInstance *comp, ValueReference vr, const uint32_t value[], size_t *index)
{
  SET_NOT_ALLOWED("UInt32");
}
#endif

#ifndef SET_INT64
Status setInt64(ModelInstance *comp, ValueReference vr, const int64_t value[], size_t *index)
{
  SET_NOT_ALLOWED("Int64");
}
#endif

#ifndef SET_UINT64
Status setUInt64(ModelInstance *comp, ValueReference vr, const uint64_t value[], size_t *index)
{
  SET_NOT_ALLOWED("UInt64");
}
#endif

#ifndef SET_BOOLEAN
Status setBoolean(ModelInstance *comp, ValueReference vr, const bool value[], size_t *index)
{
  SET_NOT_ALLOWED("Boolean");
}
#endif

#ifndef SET_STRING
Status setString(ModelInstance *comp, ValueReference vr, const char *const value[], size_t *index)
{
  SET_NOT_ALLOWED("String");
}
#endif

#ifndef SET_BINARY
Status setBinary(ModelInstance *comp, ValueReference vr, const size_t size[], const char *const value[], size_t *index)
{
  UNUSED(size);
  SET_NOT_ALLOWED("Binary");
}
#endif

#ifndef ACTIVATE_CLOCK
Status activateClock(ModelInstance *comp, ValueReference vr)
{
  UNUSED(comp);
  UNUSED(vr);
  return Error;
}
#endif

#ifndef GET_CLOCK
Status getClock(ModelInstance *comp, ValueReference vr, bool *value)
{
  UNUSED(comp);
  UNUSED(vr);
  UNUSED(value);
  return Error;
}
#endif

#ifndef GET_INTERVAL
Status getInterval(ModelInstance *comp, ValueReference vr, double *interval, int *qualifier)
{
  UNUSED(comp);
  UNUSED(vr);
  UNUSED(interval);
  UNUSED(qualifier);
  return Error;
}
#endif

#ifndef ACTIVATE_MODEL_PARTITION
Status activateModelPartition(ModelInstance *comp, ValueReference vr, double activationTime)
{
  UNUSED(comp);
  UNUSED(vr);
  UNUSED(activationTime);
  return Error;
}
#endif

#if NX < 1
void getContinuousStates(ModelInstance *comp, double x[], size_t nx)
{
  UNUSED(comp);
  UNUSED(x);
  UNUSED(nx);
}

void setContinuousStates(ModelInstance *comp, const double x[], size_t nx)
{
  UNUSED(comp);
  UNUSED(x);
  UNUSED(nx);
}

void getDerivatives(ModelInstance *comp, double dx[], size_t nx)
{
  UNUSED(comp);
  UNUSED(dx);
  UNUSED(nx);
}
#endif

#ifndef GET_PARTIAL_DERIVATIVE
Status getPartialDerivative(ModelInstance *comp, ValueReference unknown, ValueReference known, double *partialDerivative)
{
  UNUSED(comp);
  UNUSED(unknown);
  UNUSED(known);
  UNUSED(partialDerivative);
  return Error;
}
#endif

void *getFMUState(ModelInstance *comp)
{

  ModelInstance *fmuState = (ModelInstance *) calloc(1, sizeof(ModelInstance));

  memcpy(fmuState, comp, sizeof(ModelInstance));

  return fmuState;
}

void setFMUState(ModelInstance *comp, void *FMUState)
{

  ModelInstance *s = (ModelInstance *) FMUState;

  comp->startTime                         = s->startTime;
  comp->time                              = s->time;
  comp->solverStepSize                    = s->solverStepSize;
  comp->status                            = s->status;
  comp->state                             = s->state;
  comp->newDiscreteStatesNeeded           = s->newDiscreteStatesNeeded;
  comp->terminateSimulation               = s->terminateSimulation;
  comp->nominalsOfContinuousStatesChanged = s->nominalsOfContinuousStatesChanged;
  comp->valuesOfContinuousStatesChanged   = s->valuesOfContinuousStatesChanged;
  comp->nextEventTimeDefined              = s->nextEventTimeDefined;
  comp->nextEventTime                     = s->nextEventTime;
  comp->clocksTicked                      = s->clocksTicked;
  comp->isDirtyValues                     = s->isDirtyValues;
  comp->modelData                         = s->modelData;
#if NZ > 0
  memcpy(comp->z, s->z, NZ * sizeof(double));
#endif
  comp->nSteps = s->nSteps;
}

void doFixedStep(ModelInstance *comp, bool *stateEvent, bool *timeEvent)
{

  UNUSED(comp);
  UNUSED(stateEvent);
  UNUSED(timeEvent);
  logError(comp, "Fixed step not yet implemented. Please use doAlphaStep or implement method.");
}

void doAlphaStep(ModelInstance *comp, bool *stateEvent, bool *timeEvent)
{

  // get current time step size
  double dt = comp->solverStepSize;

  // create local variables
  double t         = comp->time;
  double m[3]      = {0};
  double force_in  = M(force_in);
  double force_out = 0;
  double gamma     = 0.5 - M(alpha_m) + M(alpha_f);
  double beta      = 0.25 * (gamma + 0.5);
  double u_new     = 0;
  double v_new     = 0;
  double a_new     = 0;
  double k_bar     = 0;
  double stiffness = M(spring_fixed_c) + M(spring_middle_c);

  // do generalized alpha step
  m[0] = (1 - M(alpha_m)) / (beta * dt * dt);
  m[1] = (1 - M(alpha_m)) / (beta * dt);
  m[2] = (1 - M(alpha_m) - 2 * beta) / (2 * beta);

  k_bar = stiffness * (1 - M(alpha_f)) + m[0] * M(mass_m);
  u_new = (force_in - M(alpha_f) * stiffness * M(mass_u) + M(mass_m) * (m[0] * M(mass_u) + m[1] * M(mass_v) + m[2] * M(mass_a))) / k_bar;
  a_new = 1.0 / (beta * dt * dt) * (u_new - M(mass_u) - dt * M(mass_v)) - (1 - 2 * beta) / (2 * beta) * M(mass_a);
  v_new = M(mass_v) + dt * ((1 - gamma) * M(mass_a) + gamma * a_new);

  // compute force
  force_out = M(spring_middle_c) * u_new;

  // store new values
  M(mass_u)    = u_new;
  M(mass_v)    = v_new;
  M(mass_a)    = a_new;
  M(force_out) = force_out;

  // advance time
  comp->nSteps++;
  comp->time = t + comp->solverStepSize;

  // state event
  *stateEvent = false;
  UNUSED(timeEvent);
}
