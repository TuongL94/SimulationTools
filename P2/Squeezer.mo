model Squeezer
  // Ledernas positioner
  constant Real xa = -0.06934;
  constant Real ya = -0.00227;
  constant Real xb = -0.03635;
  constant Real yb = 0.03273;
  constant Real xc = 0.014;
  constant Real yc = 0.072;

  // Massor
  constant Real m1 = 0.04325;
  constant Real m2 = 0.00365;
  constant Real m3 = 0.02373;
  constant Real m4 = 0.00706;
  constant Real m5 = 0.07050;
  constant Real m6 = 0.00706;
  constant Real m7 = 0.05498;

  // Tröghetsmoment
  constant Real I1 = 2.194e-6;
  constant Real I2 = 4.410e-7;
  constant Real I3 = 5.255e-6;
  constant Real I4 = 5.667e-7;
  constant Real I5 = 1.169e-5;
  constant Real I6 = 5.667e-7;
  constant Real I7 = 1.912e-5;



  // Geometri
  constant Real d = 0.028;
  constant Real e = 0.02;
  constant Real rr = 0.007;
  constant Real ta = 0.02308;
  constant Real tb = 0.00916;
  constant Real u = 0.04;
  constant Real ua = 0.01228;
  constant Real ub = 0.00449;
  constant Real zf = 0.02;
  constant Real zt = 0.04;


  function getCM
    input Real L1;
    input Real L2;
    input Real L3;
    output Real CM_loc;
  protected 
    Real w1;
    Real w2;
  algorithm 
    w1 := sqrt(L1^2 + L3^2);
    w2 := sqrt(L2^2 + L3^2);
    CM_loc := (L1*w1 + L2*w2)/(2*w1 + 2*w2);
  end getCM;

  inner Modelica.Mechanics.MultiBody.World world(axisLength=0.01, nominalLength=
       1/10)
    annotation (Placement(transformation(extent={{6,6},{-6,-6}},
        rotation=180,
        origin={84,-26})));
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(r={xa,ya,
        0}, animation=false)
            annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-30,-6})));
  Modelica.Mechanics.MultiBody.Joints.Revolute A(w(fixed=false), phi(fixed=false,
        start=-0.5235987755983))
    annotation (Placement(transformation(extent={{5,-5},{-5,5}},
        rotation=90,
        origin={-87,-39})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K7(
    m=m7,
    I_33=I7,
    shapeType="beam",
    v_0(fixed=true),
    r={u,0,0},
    animateSphere=true,
    r_CM={getCM(
        ua,
        u - ua,
        ub),-ub/2,0})
    annotation (Placement(transformation(extent={{-52,-54},{-32,-34}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K6(
    shapeType="beam",
    v_0(fixed=true),
    animateSphere=true,
    m=m6,
    I_33=I6,
    r={zf,0,0},
    r_CM={zf/2,0,0})
             annotation (Placement(transformation(extent={{6,-46},{26,-26}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute(phi(fixed=false, start=1.0471975511966))
    annotation (Placement(transformation(extent={{-26,-50},{-16,-40}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K5(
    shapeType="beam",
    v_0(fixed=true),
    animateSphere=true,
    r={zt,0,0},
    r_CM={getCM(
        ta,
        zt - ta,
        tb),tb/2,0},
    m=m5,
    I_33=I5) annotation (Placement(transformation(extent={{-82,4},{-62,24}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute A1(phi(fixed=false), w(fixed=false))
    annotation (Placement(transformation(extent={{5,-5},{-5,5}},
        rotation=90,
        origin={-93,-13})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute1
    annotation (Placement(transformation(extent={{-54,10},{-44,20}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K4(
    shapeType="beam",
    v_0(fixed=true),
    animateSphere=true,
    r={e,0,0},
    r_CM={e/2,0,0},
    m=m4,
    I_33=I4) annotation (Placement(transformation(extent={{-18,2},{2,22}})));
  Modelica.Mechanics.MultiBody.Joints.RevolutePlanarLoopConstraint revolute2
    annotation (Placement(transformation(extent={{28,8},{38,18}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K2(
    shapeType="beam",
    animateSphere=true,
    r={d,0,0},
    r_CM={d/2,0,0},
    m=m2,
    I_33=I2,
    v_0(fixed=false))
             annotation (Placement(transformation(extent={{68,6},{88,26}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(phi(fixed=false, start=
         0))
    annotation (Placement(transformation(extent={{48,12},{58,22}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute4(phi(fixed=false, start=
         0))
    annotation (Placement(transformation(extent={{98,10},{108,20}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K1(
    shapeType="beam",
    v_0(fixed=true),
    animateSphere=true,
    r={rr,0,0},
    r_CM={rr/2,0,0},
    m=m1,
    I_33=I1) annotation (Placement(transformation(extent={{124,-24},{144,-4}})));
  Modelica.Mechanics.MultiBody.Joints.RevolutePlanarLoopConstraint revolute5
    annotation (Placement(transformation(extent={{102,-24},{112,-14}})));
equation 
  connect(world.frame_b, fixedTranslation.frame_a) annotation (Line(
      points={{90,-26},{90,-6},{34,-6},{-20,-6}},
      color={95,95,95},
      thickness=0.5));
  connect(fixedTranslation.frame_b, A.frame_a) annotation (Line(
      points={{-40,-6},{-48,-6},{-48,-24},{-87,-24},{-87,-34}},
      color={95,95,95},
      thickness=0.5));
  connect(K7.frame_a, A.frame_b) annotation (Line(
      points={{-52,-44},{-87,-44}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute.frame_a, K7.frame_b) annotation (Line(
      points={{-26,-45},{-30,-45},{-30,-44},{-32,-44}},
      color={95,95,95},
      thickness=0.5));
  connect(K6.frame_a, revolute.frame_b) annotation (Line(
      points={{6,-36},{-12,-36},{-12,-45},{-16,-45}},
      color={95,95,95},
      thickness=0.5));
  connect(A1.frame_b, A.frame_a) annotation (Line(
      points={{-93,-18},{-90,-18},{-90,-34},{-87,-34}},
      color={95,95,95},
      thickness=0.5));
  connect(K5.frame_a, A1.frame_a) annotation (Line(
      points={{-82,14},{-88,14},{-88,-8},{-93,-8}},
      color={95,95,95},
      thickness=0.5));
  connect(K5.frame_b, revolute1.frame_a) annotation (Line(
      points={{-62,14},{-64,14},{-64,15},{-54,15}},
      color={95,95,95},
      thickness=0.5));
  connect(K4.frame_a, revolute1.frame_b) annotation (Line(
      points={{-18,12},{-34,12},{-34,15},{-44,15}},
      color={95,95,95},
      thickness=0.5));
  connect(K4.frame_b, revolute2.frame_a) annotation (Line(
      points={{2,12},{16,12},{16,13},{28,13}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute2.frame_b, K6.frame_b) annotation (Line(
      points={{38,13},{38,14},{38,-10},{38,-36},{26,-36}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute3.frame_a, revolute2.frame_b) annotation (Line(
      points={{48,17},{44,17},{44,13},{38,13}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute3.frame_b, K2.frame_a) annotation (Line(
      points={{58,17},{62,17},{62,16},{68,16}},
      color={95,95,95},
      thickness=0.5));
  connect(K2.frame_b, revolute4.frame_a) annotation (Line(
      points={{88,16},{98,16},{98,15}},
      color={95,95,95},
      thickness=0.5));
  connect(K1.frame_b, revolute4.frame_b) annotation (Line(
      points={{144,-14},{150,-14},{150,15},{108,15}},
      color={95,95,95},
      thickness=0.5));
  connect(K1.frame_a, revolute5.frame_b) annotation (Line(
      points={{124,-14},{120,-14},{120,-19},{112,-19}},
      color={95,95,95},
      thickness=0.5));
  connect(world.frame_b, revolute5.frame_a) annotation (Line(
      points={{90,-26},{96,-26},{96,-19},{102,-19}},
      color={95,95,95},
      thickness=0.5));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Modelica(version="3.2.2")));


end Squeezer;
