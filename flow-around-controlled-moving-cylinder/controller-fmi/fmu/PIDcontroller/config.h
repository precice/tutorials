#ifndef config_h
#define config_h

#include <stdbool.h> // for bool

// define class name and unique id
#define MODEL_IDENTIFIER PIDcontroller
#define INSTANTIATION_TOKEN "{1AE5E10D-9521-4DE3-80B9-D0EAAA7D5AF1}"

#define CO_SIMULATION
#define MODEL_EXCHANGE

// define model size
#define NX 2
#define NZ 1

#define SET_FLOAT64
#define SET_BOOLEAN
#define GET_OUTPUT_DERIVATIVE
#define GET_BOOLEAN
#define EVENT_UPDATE

#define FIXED_SOLVER_STEP 1e-3
#define DEFAULT_STOP_TIME 3

typedef enum {

  vr_time,
  vr_r,
  vr_e,
  vr_e_ls,
  vr_kp,
  vr_ki,
  vr_kd,
  vr_P,
  vr_I,
  vr_D,
  vr_y_1,
  vr_y_2,
  vr_u_1,
  vr_u_2,
  vr_I_max,
  vr_use_implicit_method,
  vr_compute_u_1

} ValueReference;

typedef struct {

  double r;
  double e;
  double e_ls;
  double kp;
  double ki;
  double kd;
  double P;
  double I;
  double D;
  double y_1;
  double y_2;
  double u_1;
  double u_2;
  double I_max;
  bool   use_implicit_method;
  bool   compute_u_1;

} ModelData;

#endif /* config_h */
