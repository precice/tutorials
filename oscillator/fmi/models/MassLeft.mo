model MassLeft
  Modelica.Mechanics.Translational.Components.Spring spring1(c = 4*3.1416*3.1416, s_rel(fixed = false))  annotation(
    Placement(visible = true, transformation(origin = {-58, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.Mass mass1(L = 0, m = 1, s(fixed = true, start = 0), v(fixed = true, start = 0))  annotation(
    Placement(visible = true, transformation(origin = {-22, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.Fixed fixed annotation(
    Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Mechanics.Translational.Components.Spring spring12(c = 16*3.1416*3.1416, s_rel0 = 0) annotation(
    Placement(visible = true, transformation(origin = {14, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -180)));
  Modelica.Blocks.Interfaces.RealOutput disp1 annotation(
    Placement(visible = true, transformation(origin = {52, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {22, -52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Sources.Position position(v(fixed = true, start = 0))  annotation(
    Placement(visible = true, transformation(origin = {54, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Sensors.PositionSensor positionSensor annotation(
    Placement(visible = true, transformation(origin = {-2, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput disp2 annotation(
    Placement(visible = true, transformation(origin = {14, 64}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {14, 64}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
equation
  connect(spring1.flange_a, fixed.flange) annotation(
    Line(points = {{-68, 0}, {-80, 0}}, color = {0, 127, 0}));
  connect(mass1.flange_a, spring1.flange_b) annotation(
    Line(points = {{-32, 0}, {-48, 0}}, color = {0, 127, 0}));
  connect(mass1.flange_b, positionSensor.flange) annotation(
    Line(points = {{-12, 0}, {-12, -54}}, color = {0, 127, 0}));
  connect(disp1, positionSensor.s) annotation(
    Line(points = {{52, -44}, {52, -55}, {9, -55}, {9, -54}}, color = {0, 0, 127}));
  connect(spring12.flange_a, mass1.flange_b) annotation(
    Line(points = {{4, 0}, {-12, 0}}, color = {0, 127, 0}));
  connect(spring12.flange_b, position.flange) annotation(
    Line(points = {{24, 0}, {44, 0}}, color = {0, 127, 0}));
  connect(disp2, position.s_ref) annotation(
    Line(points = {{14, 64}, {66, 64}, {66, 0}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "4.0.0")),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "trapezoid"));
end MassLeft;