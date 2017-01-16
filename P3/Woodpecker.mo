model Woodpecker_splitEq

  // test: lambda1 och 2 till "protected"

  // Konstanter
  constant Real m_S(unit="kg") = 3.0e-4;
  constant Real m_B(unit="kg") = 4.5e-3;
  constant Real J_S(unit="kg.m2") = 5.0e-9;
  constant Real J_B(unit="kg.m2") = 7.0e-7;
  constant Real r_0(unit="m") = 2.5e-3;
  constant Real r_S(unit="m") = 3.1e-3;
  constant Real h_S(unit="m") = 2.0e-2;
  constant Real l_S(unit="m") = 1.0e-2;
  constant Real l_G(unit="m") = 1.5e-2;
  constant Real l_B(unit="m") = 2.01e-2;
  constant Real h_B(unit="m") = 2.0e-2;
  constant Real c_p(unit="N.m/rad") = 5.6e-3;
  constant Real g(unit="m/s2") = 9.81;

  discrete Integer state(start = 1, fixed = true);

  Real phi_S(start = 0.00, unit="rad");
  Real phi_S_dot(start = .5, unit="rad/s", fixed=true);
  Real phi_B(start = 0, unit="rad");
  Real phi_B_dot(start = .0,unit="rad/s");
  Real z(start = 0, unit="m", fixed=true);
  Real z_dot(start = 0, unit="m/s");

  Real mom_pre(start=0);

  Boolean sleeveBlocking(start=false);

  Real testcond2(start = 0);
  Real testcond3(start = 0);
  Real testcond4(start = 0);
  
  

  Real lambda_1(start = 0.0,fixed=true);
  Real lambda_2(start = 0.0,fixed=true);
  
  discrete Integer nbrOfBeakHits(start = 0);
  
  function momentum
    input Real zd;
    input Real phiSd;
    input Real phiBd;
    output Real I;
  algorithm 
    I := m_B*l_G*zd + m_B*l_S*l_G*phiSd +( J_B+m_B*l_G^2)*phiBd;
  end momentum;

  function lambda_1Change
  algorithm 
    if state == 2 or (state == 3 and phi_B_dot < 0) then
      state := 1;
    end if;
  end lambda_1Change;

//  Real mom_pre;

protected 
  
equation 
  der(phi_S) = phi_S_dot;
  der(phi_B) = phi_B_dot;
  der(z) = z_dot;
  
  if state == 1 then // lambda_1 = lambda_2 = 0
  
  
  (m_S + m_B)*der(z_dot) + m_B*l_S*der(phi_S_dot) + m_B*l_G*der(phi_B_dot) = -(m_S+m_B)*g; // 1a
  m_B*l_S*der(z_dot) +( J_S + m_B*l_S^2)*der(phi_S_dot) + m_B*l_S*l_G*der(phi_B_dot) = c_p*(phi_B - phi_S) - m_B*l_S*g; //2a
  m_B*l_G*der(z_dot) + m_B*l_S*l_G*der(phi_S_dot) +( J_B + m_B*l_G^2)*der(phi_B_dot) = c_p*(phi_S - phi_B) - m_B*l_G*g; //3a
    0.0 = lambda_1;  //1d
    0.0 = lambda_2; //1e
  else // z_dot = phi_S_dot = 0
        
    m_B*l_G*der(phi_B_dot) = -(m_S+m_B)*g - lambda_2; // 2-3a
    if state == 2 then
      // m_B*l_S*l_G*der(phi_B_dot) = c_p*(phi_B - phi_S) - m_B*l_S*g + lambda_1*h_S - lambda_2*r_S;// 2b
      m_B*l_S*l_G*der(phi_B_dot) = c_p*(phi_B - (r_0-r_S)/h_S) - m_B*l_S*g + lambda_1*h_S - lambda_2*r_S;// 2b
      (J_B + m_B*l_G^2)*der(phi_B_dot) = c_p*((r_0-r_S)/h_S - phi_B) - m_B*l_G*g; //2c
      phi_S = (r_0-r_S)/h_S;// 0 = (r_S-r_0) + h_S*phi_S;
    else
      // m_B*l_S*l_G*der(phi_B_dot) = c_p*(phi_B - phi_S) - m_B*l_S*g - lambda_1*h_S - lambda_2*r_S;// 3b
      m_B*l_S*l_G*der(phi_B_dot) = c_p*(phi_B - phi_S) - m_B*l_S*g - lambda_1*h_S - lambda_2*r_S;// 3b
      (J_B + m_B*l_G^2)*der(phi_B_dot) = c_p*((r_S-r_0)/h_S - phi_B) - m_B*l_G*g; //3c
      phi_S = (r_S-r_0)/h_S; //0 = (r_S-r_0) - h_S*phi_S;
    end if;
    //0.0 = phi_S_dot;
    //(J_B + m_B*l_G^2)*der(phi_B_dot) = c_p*(phi_S - phi_B) - m_B*l_G*g;
    0.0 = z_dot;//2-3e
  end if;
  //1-3c //Här har lagts till en faktor r_S före lambda_2 då det i annat fall blir dimensionsfel
  // Spelar eventuellt inte ngn större roll eftersom lambda_2 ska vara noll

  // Här blir det problem. När tillståndet lämnar tillstånd 1 försöker Dymola kanske dela med noll?

  //when state == 4 then
    //reinit(phi_B_dot,-pre(phi_B_dot));
  //end when;

  when sleeveBlocking then
    reinit(z_dot,0);
    reinit(phi_S_dot,0);
    reinit(phi_B_dot, pre(mom_pre)/(J_B+m_B*l_G^2));
  end when;

algorithm 

  mom_pre := momentum(z_dot,phi_S_dot,phi_B_dot);
  testcond2 := h_S*phi_S + (r_S - r_0);
  testcond3 := h_S*phi_S - (r_S - r_0);
  testcond4 := h_B*phi_B - (l_S + l_G - l_B - r_0);

  when h_S*phi_S < -(r_S - r_0) then
    if state == 1 then
      state := 2;
    end if;
  end when;

  when h_S*phi_S > (r_S - r_0) then
    if state == 1 then
      state := 3;
    end if;
  end when;

  when lambda_1 > 0 then
    if state == 2 or (state == 3 and phi_B_dot < 0) then
     state := 1;
    end if;
  end when;

  when lambda_1 < 0 then
    if state == 2 or (state == 3 and phi_B_dot < 0) then
     state := 1;
    end if;
  end when;

  when h_B*phi_B > l_S + l_G - l_B - r_0 then
    if state == 3 then //if state == 3 and phi_B_dot > 0 then
      state := 4;
    end if;
  end when;

  when state == 4 then
    nbrOfBeakHits := nbrOfBeakHits + 1;
    reinit(phi_B_dot,-pre(phi_B_dot));
    state := 3;
  end when;

  when state == 2 then
      sleeveBlocking :=true;
  end when;

  when state == 3 then
    sleeveBlocking :=true;
  end when;

  when sleeveBlocking then
    sleeveBlocking :=false;
  end when;

end Woodpecker_splitEq;
