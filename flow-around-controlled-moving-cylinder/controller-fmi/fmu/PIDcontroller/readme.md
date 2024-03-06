# PID Controller

The PID Controller implements the following system of equations:

```
explicit: Forward euler
e = r - y
u = kp*e + ki* (I + e*dt) + kd*(e-e_ls)*(1/dt)

implicit: Trapezoid
e = r - y
u = kp*e + ki * (I + (e + e_ls) * dt) + kd*(e-e_ls)*(1/dt)

```

with the variables

| Variable 		| Start | Causality | Variability | Description
|:---------------------|------:|-----------|-------------|:---------------
| u_1      		|     0 | output    | continuous  | Control output 1
| u_2      		|     0 | output    | continuous  | Control output 2
| y_1      		|     0 | parameter | tunable     | Control input 1 
| y_2      		|     0 | parameter | tunable     | Control input 2
| r        		|     0 | parameter | tunable     | Reference value
| e        		|     - | local     | continuous  | Error between input and reference
| e_ls     		|     - | local     | continuous  | Error between input and reference from last time step
| kp       		|     0 | parameter | tunable     | Proportional gain
| ki       		|     0 | parameter | tunable     | Integral gain
| kd       		|     0 | parameter | tunable     | Derivative gain
| dt       		|  1e-2 | parameter | tunable     | Solver step size
| P        		|     - | local     | continuous  | Proportional term
| I        		|     - | local     | continuous  | Integral term
| D        		|     - | local     | continuous  | Differential term
| I_max    		|   100 | local     | continuous  | Maximumim for I to avoid integrator windup
| compute_u_1   	|  true | local     | continuous  | Flag to compute either u_1 or u_2
| use_implicit_method  | false | local     | continuous  | Flag to use explicit or implicit method

