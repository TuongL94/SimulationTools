within ;
model WoodPecker

  // Constants
  constant Real m_S(unit="kg") = 3.0e-4;
  constant Real m_B(unit="kg") = 4.5e-3;
  constant Real J_S(unit="kg.m2") = 5.0e-9;
  constant Real J_B(unit="kg.m2") = 7.0e-7;
  constant Real r_0(unit="m") = 2.5e-3;
  constant Real r_S(unit="m") = 3.1e-3;
  constant Real l_S(unit="m") = 1.0e-2;
  constant Real l_G(unit="m") = 1.5e-2;
  constant Real l_B(unit="m") = 2.01e-2;
  constant Real h_S(unit="m") = 1.16e-2;
  constant Real h_B(unit="m") = 2.0e-2;
  constant Real c_p(unit="N.m/rad") = 5.6e-3;
  constant Real g(unit="m/s2") = 9.81;
  constant Real epsilon = 1e-8;

  Integer state(start = 1,fixed = true);
  Integer nbrOfBeakHits(start = 0, fixed = true);

  Real phi_S(start = 0.0, unit="rad");
  Real phi_S_dot(start = 0.5, unit="rad/s", fixed = true);
  Real phi_B(start = 0.0, unit="rad", fixed = true);
  Real phi_B_dot(start = 0.5, unit="rad/s");
  Real z(start = 0.0, unit="m");
  Real z_dot(start = 0.0, unit="m/s", fixed = true);
  Real lambda_1(start = 0);
  Real lambda_2(start = 0);

   function momentum
     input Real z_dot_minus;
     input Real phi_S_dot_minus;
     input Real phi_B_dot_minus;
     output Real phi_B_dot_plus;
   algorithm
        phi_B_dot_plus := (m_B*l_G*z_dot_minus + m_B*l_S*l_G*phi_S_dot_minus +
        (J_B+m_B*l_G^2)*phi_B_dot_minus)/(J_B+m_B*l_G^2);
   end momentum;

equation
  der(phi_S) = phi_S_dot;
  der(phi_B) = phi_B_dot;
  der(z) = z_dot;

  if state == 1 then
    (m_S + m_B)*der(z_dot) + m_B*l_S*der(phi_S_dot) + m_B*l_G*der(phi_B_dot) = -
      (m_S + m_B)*g;
    (m_B*l_S)*der(z_dot) + (J_S + m_B*l_S^2)*der(phi_S_dot) + (m_B*l_S*l_G)*der(
       phi_B_dot) = c_p*(phi_B - phi_S) - m_B*l_S*g - lambda_1;
    m_B*l_G*der(z_dot) + (m_B*l_S*l_G)*der(phi_S_dot) + (J_B + m_B*l_G^2)*der(
      phi_B_dot) = c_p*(phi_S - phi_B) - m_B*l_G*g - lambda_2;
    lambda_1 = 0;
    lambda_2 = 0;
  elseif state == 2 then
    (m_S + m_B)*der(z_dot) + m_B*l_S*der(phi_S_dot) + m_B*l_G*der(phi_B_dot) = -
      (m_S + m_B)*g - lambda_2;
    (m_B*l_S)*der(z_dot) + (J_S + m_B*l_S^2)*der(phi_S_dot) + (m_B*l_S*l_G)*der(
       phi_B_dot) = c_p*(phi_B - phi_S) - m_B*l_S*g - h_S*lambda_1 - lambda_2*
      r_S;
    m_B*l_G*der(z_dot) + (m_B*l_S*l_G)*der(phi_S_dot) + (J_B + m_B*l_G^2)*der(
      phi_B_dot) = c_p*(phi_S - phi_B) - m_B*l_G*g;
    der(phi_S_dot) = 0;
    der(z_dot) + r_S*der(phi_S_dot) = 0;
  else
    (m_S + m_B)*der(z_dot) + m_B*l_S*der(phi_S_dot) + m_B*l_G*der(phi_B_dot) = -
      (m_S + m_B)*g - lambda_2;
    (m_B*l_S)*der(z_dot) + (J_S + m_B*l_S^2)*der(phi_S_dot) + (m_B*l_S*l_G)*der(
       phi_B_dot) = c_p*(phi_B - phi_S) - m_B*l_S*g + h_S*lambda_1 - lambda_2*
      r_S;
    m_B*l_G*der(z_dot) + (m_B*l_S*l_G)*der(phi_S_dot) + (J_B + m_B*l_G^2)*der(
      phi_B_dot) = c_p*(phi_S - phi_B) - m_B*l_G*g;
    der(phi_S_dot) = 0;
    der(z_dot) + r_S*der(phi_S_dot) = 0;
  end if;

algorithm
  // State 1 to 2
  when {h_S*phi_S + (r_S - r_0) < 0,h_S*phi_S + (r_S - r_0) > 0} then
    if state == 1 and phi_B_dot < 0 then
      state := 2;
    end if;
  // State 1 to 3
  elsewhen {h_S*phi_S - (r_S - r_0) < 0,h_S*phi_S - (r_S - r_0) > 0} then
    if state == 1 and phi_B_dot > 0 then
      state := 3;
    end if;
  elsewhen {lambda_1 < -epsilon,lambda_1 > epsilon} then
    // State 2 to 1
    if state == 2 then
      state := 1;
      // State 3 to 1
    elseif state == 3 and phi_B_dot < 0 then
      state := 1;
    end if;
  elsewhen {h_B*phi_B - (l_S + l_G - l_B - r_0) > 0 and phi_B_dot > 0} then
    nbrOfBeakHits := nbrOfBeakHits + 1;
  end when;

equation
   when {h_B*phi_B - (l_S + l_G - l_B - r_0) > 0 and phi_B_dot > 0} then
    reinit(phi_B_dot, -pre(phi_B_dot));
   elsewhen state == 2 then
    reinit(phi_B_dot,momentum(pre(z_dot),pre(phi_S_dot),pre(phi_B_dot)));
    reinit(z_dot, 0);
    reinit(phi_S_dot, 0);
   elsewhen state == 3 then
    reinit(phi_B_dot, momentum(
      pre(z_dot),
      pre(phi_S_dot),
      pre(phi_B_dot)));
    reinit(z_dot, 0);
    reinit(phi_S_dot, 0);
  end when;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end WoodPecker;
