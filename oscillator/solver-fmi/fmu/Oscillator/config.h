#ifndef config_h
#define config_h

// define class name and unique id
#define MODEL_IDENTIFIER Oscillator
#define INSTANTIATION_TOKEN "{1AE5E10D-9521-4DE3-80B9-D0EAAA7D5AF1}"

#define CO_SIMULATION
#define MODEL_EXCHANGE

// define model size
#define NX 2
#define NZ 1

#define SET_FLOAT64
#define GET_OUTPUT_DERIVATIVE
#define EVENT_UPDATE

#define FIXED_SOLVER_STEP 1e-3
#define DEFAULT_STOP_TIME 3

typedef enum {
  vr_time,
  vr_mass_m,
  vr_mass_u,
  vr_mass_v,
  vr_mass_a,
  vr_spring_fixed_c,
  vr_spring_middle_c,
  vr_force_in,
  vr_force_out,
  vr_alpha_f,
  vr_alpha_m
} ValueReference;

typedef struct {

  double mass_m;
  double mass_u;
  double mass_v;
  double mass_a;
  double spring_fixed_c;
  double spring_middle_c;
  double force_in;
  double force_out;
  double alpha_f;
  double alpha_m;

} ModelData;

#endif /* config_h */
