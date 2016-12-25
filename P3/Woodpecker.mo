model Woodpecker

  // Konstanter
  constant Real m_S(unit="kg") = 3.0e-4;
  constant Real m_B(unit="kg") = 4.5e-3;
  constant Real J_S(unit="kg.m") = 5.0e-9;
  constant Real J_B(unit="kg.m") = 7.0e-7;
  constant Real r_0(unit="m") = 2.5e-3;
  constant Real r_S(unit="m") = 3.1e-3;
  constant Real l_S(unit="m") = 1.0e-2;
  constant Real l_G(unit="m") = 1.5e-2;
  constant Real l_B(unit="m") = 2.01e-2;
  constant Real h_B(unit="m") = 2.0e-2;
  constant Real c_p(unit="N/rad") = 5.6e-3;
  constant Real g(unit="m/s2") = 9.81;

  Integer state(start = 1, fixed = true);

  Real phi_S(start = 0, unit="rad");
  Real phi_S_dot(start = 0, unit="rad/s");
  Real phi_B(start = 0, unit="rad");
  Real phi_B_dot(start = 0, unit="rad/s");
  Real z(start = 0, unit="m");
  Real z_dot(start = 0, unit="m/s");
  
  Real lambda_2(start = 0);

  function KD
    input Integer a;
    input Integer b;
    output Real res;
  algorithm 
    if a == b then
      res = 1;
    else
      res = 0;
    end if;
  end KD;

equation 
  der(phi_S) = phi_S_dot;
  der(phi_B) = phi_B_dot;
  der(z) = z_dot;
  (m_S + m_B)*der(z_dot) + m_B*l_S*der(phi_S_dot) + m_B*l_G*der(phi_B_dot) = -(m_S+m_B)*g + (KD(state,1) + KD(state,2))*(-lambda_2); // 1-3a
  
algorithm
  when h_S*phi_S == -(r_S - r_0) then
    if state == 1 and phi_B_dot < 0 then
      state := 2;
    end if;
  end when;
    
  
end Woodpecker;
