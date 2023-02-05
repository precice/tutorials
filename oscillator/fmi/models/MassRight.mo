model MassRight
  Modelica.Mechanics.Translational.Components.Fixed fixed annotation(
    Placement(visible = true, transformation(origin = {82, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Mechanics.Translational.Components.Spring spring2(c = 4*3.1416*3.1416, s_rel(fixed = false), s_rel0 = 0) annotation(
    Placement(visible = true, transformation(origin = {44, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.Spring spring12(c = 16*3.1416*3.1416, s_rel(start = 0), s_rel0 = 0) annotation(
    Placement(visible = true, transformation(origin = {-30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.Mass mass2(L = 0, m = 1, s(fixed = true, start = 0), v(displayUnit = "Gm/s", fixed = true, start = 0)) annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Sensors.PositionSensor positionSensor annotation(
    Placement(visible = true, transformation(origin = {0, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Sources.Position position(v(fixed = true, start = 0))  annotation(
    Placement(visible = true, transformation(origin = {-68, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput disp1 annotation(
    Placement(visible = true, transformation(origin = {-46, 62}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-46, 62}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput disp2 annotation(
    Placement(visible = true, transformation(origin = {80, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(spring2.flange_b, fixed.flange) annotation(
    Line(points = {{54, 0}, {54, 23}, {82, 23}, {82, 0}}, color = {0, 127, 0}));
  connect(mass2.flange_a, spring12.flange_b) annotation(
    Line(points = {{-10, 0}, {-20, 0}}, color = {0, 127, 0}));
  connect(position.flange, spring12.flange_a) annotation(
    Line(points = {{-58, 0}, {-40, 0}}, color = {0, 127, 0}));
  connect(disp1, position.s_ref) annotation(
    Line(points = {{-46, 62}, {-80, 62}, {-80, 0}}, color = {0, 0, 127}));
  connect(positionSensor.s, disp2) annotation(
    Line(points = {{11, -58}, {47.5, -58}, {47.5, -60}, {80, -60}}, color = {0, 0, 127}));
  connect(mass2.flange_a, positionSensor.flange) annotation(
    Line(points = {{-10, 0}, {-10, -58}}, color = {0, 127, 0}));
  connect(mass2.flange_b, spring2.flange_a) annotation(
    Line(points = {{10, 0}, {34, 0}}, color = {0, 127, 0}));
  annotation(
    uses(Modelica(version = "4.0.0")),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "trapezoid"));
end MassRight;