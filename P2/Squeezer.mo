model Squeezer

  // parameter Real Theta = 0.0; // Theta = -60 fungerar


  constant Real PI = Modelica.Constants.pi;
  Real beta(start = 0.0);
  //   Real Theta(start = 0.0, fixed=true);
  Real Theta(start = 0.0);
  Real gamma(start = 0.0);
  Real Phi(start = 0.0);
  Real delta(start = 0.0);
  Real Omega(start = 0.0);
  Real epsilon(start = 0.0);


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
  constant Real c_0 = 4530;
  constant Real d = 0.028;
  constant Real da = 0.0115;
  constant Real e = 0.02;
  constant Real ea = 0.01421;
  constant Real fa = ea;
  constant Real l_0 = 0.07785;
  constant Real ra = 0.00092;
  constant Real rr = 0.007;
  constant Real sa = 0.01874;
  constant Real sb = 0.01043;
  constant Real sc = 0.018;
  constant Real sd = 0.02;
  constant Real ss = 0.035;
  constant Real ta = 0.02308;
  constant Real tb = 0.00916;
  constant Real u = 0.04;
  constant Real ua = 0.01228;
  constant Real ub = 0.00449;
  constant Real zf = 0.02;
  constant Real zt = 0.04;


  // Moment
  constant Real mom = 0.033;

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
       1/10) annotation (Placement(transformation(
        extent={{6,6},{-6,-6}},
        rotation=180,
        origin={84,-26})));
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation OA(r={xa,ya,0}, animation=
       false) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-30,-6})));
  Modelica.Mechanics.MultiBody.Joints.Revolute A(w(fixed=true), phi(fixed=false,
        start=-0.26179938779915))
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
    r_CM={ua,-ub,0})
    annotation (Placement(transformation(extent={{-52,-54},{-32,-34}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K6(
    shapeType="beam",
    animateSphere=true,
    m=m6,
    I_33=I6,
    r={zf,0,0},
    r_CM={zf - fa,0,0},
    v_0(fixed=false))
             annotation (Placement(transformation(extent={{8,-46},{28,-26}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute K6_K7(phi(fixed=false, start=
          1.0471975511966))
    annotation (Placement(transformation(extent={{-26,-50},{-16,-40}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K5(
    shapeType="beam",
    v_0(fixed=true),
    animateSphere=true,
    r={zt,0,0},
    m=m5,
    I_33=I5,
    r_CM={ta,tb,0})
             annotation (Placement(transformation(extent={{-82,4},{-62,24}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute A1(                  w(fixed=
          false), phi(fixed=false, start=0.69813170079773))
    annotation (Placement(transformation(extent={{-5,-5},{5,5}},
        rotation=90,
        origin={-93,-13})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K4(
    shapeType="beam",
    animateSphere=true,
    r={e,0,0},
    m=m4,
    I_33=I4,
    r_CM={e - ea,0,0},
    v_0(fixed=false))
             annotation (Placement(transformation(extent={{-18,4},{2,24}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K2(
    shapeType="beam",
    animateSphere=true,
    r={d,0,0},
    m=m2,
    I_33=I2,
    v_0(fixed=false),
    r_CM={d - da,0,0})
             annotation (Placement(transformation(extent={{116,14},{136,34}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K1(
    shapeType="beam",
    v_0(fixed=true),
    animateSphere=true,
    m=m1,
    I_33=I1,
    r={rr,0,0},
    r_CM={ra,0,0})
             annotation (Placement(transformation(extent={{124,-24},{144,-4}})));
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation OB(animation=false, r={xb,
        yb,0}) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={28,42})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K3(
    shapeType="beam",
    animateSphere=true,
    v_0(fixed=false),
    r={ss,0,0},
    r_CM={sa,sb,0},
    m=m3,
    I_33=I3) annotation (Placement(transformation(extent={{62,54},{82,74}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute B(phi(fixed=false, start=-1.0471975511966),
      w(fixed=false))
    annotation (Placement(transformation(extent={{28,58},{38,68}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(useAxisFlange=true,
      phi(fixed=false, start=1.5707963267949))
    annotation (Placement(transformation(extent={{106,-22},{116,-12}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(phi(fixed=true, start=0),
      w(fixed=false))
    annotation (Placement(transformation(extent={{-5,-5},{5,5}},
        rotation=90,
        origin={151,7})));
  Modelica.Mechanics.MultiBody.Joints.RevolutePlanarLoopConstraint revolute4
    annotation (Placement(transformation(extent={{100,18},{110,28}})));
  Modelica.Mechanics.MultiBody.Forces.Spring spring(c=c_0, s_unstretched=l_0)
    annotation (Placement(transformation(extent={{124,80},{144,100}})));
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation OC(animation=false, r={xc,
        yc,0}) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={170,24})));
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation BD(r={sc,sd,0}, animation=
       false) annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=90,
        origin={63,83})));
  Modelica.Mechanics.Rotational.Sources.ConstantTorque constantTorque(
      tau_constant=mom, useSupport=true)
    annotation (Placement(transformation(extent={{76,-74},{96,-54}})));
  Modelica.Mechanics.Rotational.Components.Fixed fixed
    annotation (Placement(transformation(extent={{76,-98},{96,-78}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute6(phi(fixed=false, start
        =4.7123889803847))
    annotation (Placement(transformation(extent={{-46,8},{-36,18}})));
  Modelica.Mechanics.MultiBody.Joints.RevolutePlanarLoopConstraint K4_K2
    annotation (Placement(transformation(extent={{30,6},{40,16}})));
  Modelica.Mechanics.MultiBody.Joints.RevolutePlanarLoopConstraint revolute2
    annotation (Placement(transformation(extent={{48,-28},{58,-18}})));
equation 
  connect(world.frame_b, OA.frame_a) annotation (Line(
      points={{90,-26},{90,-6},{34,-6},{-20,-6}},
      color={95,95,95},
      thickness=0.5));
  connect(OA.frame_b, A.frame_a) annotation (Line(
      points={{-40,-6},{-48,-6},{-48,-24},{-87,-24},{-87,-34}},
      color={95,95,95},
      thickness=0.5));
  connect(K7.frame_a, A.frame_b) annotation (Line(
      points={{-52,-44},{-87,-44}},
      color={95,95,95},
      thickness=0.5));
  connect(K6_K7.frame_a, K7.frame_b) annotation (Line(
      points={{-26,-45},{-30,-45},{-30,-44},{-32,-44}},
      color={95,95,95},
      thickness=0.5));
  connect(K6.frame_a, K6_K7.frame_b) annotation (Line(
      points={{8,-36},{-12,-36},{-12,-45},{-16,-45}},
      color={95,95,95},
      thickness=0.5));
  connect(OB.frame_a, world.frame_b) annotation (Line(
      points={{28,32},{90,32},{90,-26}},
      color={95,95,95},
      thickness=0.5));
  connect(OB.frame_b, B.frame_a) annotation (Line(
      points={{28,52},{28,63}},
      color={95,95,95},
      thickness=0.5));
  connect(B.frame_b, K3.frame_a) annotation (Line(
      points={{38,63},{46,63},{46,64},{62,64}},
      color={95,95,95},
      thickness=0.5));
  connect(world.frame_b, revolute1.frame_a) annotation (Line(
      points={{90,-26},{98,-26},{98,-17},{106,-17}},
      color={95,95,95},
      thickness=0.5));
  connect(K1.frame_a, revolute1.frame_b) annotation (Line(
      points={{124,-14},{118,-14},{118,-17},{116,-17}},
      color={95,95,95},
      thickness=0.5));
  connect(K1.frame_b, revolute3.frame_a) annotation (Line(
      points={{144,-14},{148,-14},{148,2},{151,2}},
      color={95,95,95},
      thickness=0.5));
  connect(K2.frame_b, revolute3.frame_b) annotation (Line(
      points={{136,24},{140,24},{140,12},{151,12}},
      color={95,95,95},
      thickness=0.5));
  connect(K3.frame_b, revolute4.frame_a) annotation (Line(
      points={{82,64},{92,64},{92,23},{100,23}},
      color={95,95,95},
      thickness=0.5));
  connect(K2.frame_a, revolute4.frame_b) annotation (Line(
      points={{116,24},{110,23}},
      color={95,95,95},
      thickness=0.5));
  connect(OC.frame_a, world.frame_b) annotation (Line(
      points={{170,14},{170,14},{170,-26},{90,-26}},
      color={95,95,95},
      thickness=0.5));
  connect(OC.frame_b, spring.frame_b) annotation (Line(
      points={{170,34},{170,34},{170,90},{144,90}},
      color={95,95,95},
      thickness=0.5));
  connect(K3.frame_a, BD.frame_a) annotation (Line(
      points={{62,64},{62,76},{63,76}},
      color={95,95,95},
      thickness=0.5));
  connect(spring.frame_a, BD.frame_b) annotation (Line(
      points={{124,90},{63,90}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute1.axis, constantTorque.flange) annotation (Line(points={{111,-12},
          {104,-12},{104,-64},{96,-64}}, color={0,0,0}));
  connect(constantTorque.support, fixed.flange)
    annotation (Line(points={{86,-74},{86,-88}}, color={0,0,0}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Modelica(version="3.2.2")));

algorithm 

  // Alla vinklar uttrycks i radianer

  beta := revolute1.phi;
  Theta := revolute3.phi;
  gamma := -B.phi;
  Phi := 7*PI/2 - revolute6.phi;
  Omega :=   K6_K7.phi - PI/2;
  epsilon := A.phi + PI/2;
  delta := A1.phi;




equation 
  connect(K6.frame_b, revolute2.frame_a) annotation (Line(
      points={{28,-36},{38,-36},{38,-23},{48,-23}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute2.frame_b, revolute4.frame_a) annotation (Line(
      points={{58,-23},{66,-23},{66,23},{100,23}},
      color={95,95,95},
      thickness=0.5));
  connect(A1.frame_b, K5.frame_a) annotation (Line(
      points={{-93,-8},{-88,-8},{-88,14},{-82,14}},
      color={95,95,95},
      thickness=0.5));
  connect(A1.frame_a, A.frame_a) annotation (Line(
      points={{-93,-18},{-90,-18},{-90,-24},{-87,-24},{-87,-34}},
      color={95,95,95},
      thickness=0.5));
  connect(K4.frame_b, K4_K2.frame_a) annotation (Line(
      points={{2,14},{16,14},{16,11},{30,11}},
      color={95,95,95},
      thickness=0.5));
  connect(K4_K2.frame_b, revolute4.frame_a) annotation (Line(
      points={{40,11},{54,11},{54,22},{66,22},{66,23},{100,23}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute6.frame_b, K4.frame_a) annotation (Line(
      points={{-36,13},{-28,13},{-28,14},{-18,14}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute6.frame_a, K5.frame_b) annotation (Line(
      points={{-46,13},{-54,13},{-54,14},{-62,14}},
      color={95,95,95},
      thickness=0.5));
end Squeezer;
