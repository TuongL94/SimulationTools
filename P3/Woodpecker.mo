model Woodpecker

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

  Integer state(start = 1, fixed = true);

  Real phi_S(start = 0.05, unit="rad", fixed=true);
  Real phi_S_dot(start = 0, unit="rad/s");
  Real phi_B(start = 0, unit="rad");
  Real phi_B_dot(start = 0.5,unit="rad/s", fixed=true);
  Real z(start = 0, unit="m");
  Real z_dot(start = 0, unit="m/s");

  Real lambda_1(start = 0);
  Real lambda_2(start = 0);
  
  function KD
    // Kroneckerdelta
    input Integer a;
    input Integer b;
    output Real res;
  algorithm 
    if a == b then
      res :=1;
    else
      res :=0;
    end if;
  end KD;

  function momentum
    input Real zd;
    input Real phiSd;
    input Real phiBd;
    output Real I;
  algorithm 
    I := m_B*l_G*zd + m_B*l_S*l_G*phiSd +( J_B+m_B*l_G^2)*phiBd;
  end momentum;

//  function sleeveBlocking
                          //sleeveBlocking(z_dot,phi_S_dot,phi_B_dot)

            // Lokala variabler
  //protected 
    //Real mom_pre;
//  algorithm 
    
                              //sleeve blocking

 // end sleeveBlocking;

  function lambda_1Change
  algorithm 
    if state == 2 or (state == 3 and phi_B_dot < 0) then
      state := 1;
    end if;
  end lambda_1Change;
  
protected
  Real mom_pre;
  
equation 
  der(phi_S) = phi_S_dot;
  der(phi_B) = phi_B_dot;
  der(z) = z_dot;

  (m_S + m_B)*der(z_dot) + m_B*l_S*der(phi_S_dot) + m_B*l_G*der(phi_B_dot) = -(m_S+m_B)*g + (KD(state,1) + KD(state,2))*(-lambda_2); // 1-3a
  m_B*l_S*der(z_dot) +( J_S + m_B*l_S^2)*der(phi_S_dot) + m_B*l_S*l_G*der(
    phi_B_dot) = c_p*(phi_B - phi_S) - m_B*l_S*g + lambda_1*(-KD(state, 1) +
    h_S*(-KD(state, 2) + KD(state, 3))) - lambda_2*r_S*(KD(state, 2) + KD(state, 3));
  // 1-3b
  m_B*l_G*der(z_dot) + m_B*l_S*l_G*der(phi_S_dot) +( J_B + m_B*l_G^2)*der(
    phi_B_dot) = c_p*(phi_S - phi_B) - m_B*l_G*g - lambda_2*KD(state, 1);
  //1-3c

  if state == 1 then
    0 = lambda_1;  //1d
  else
    0 = (r_S-r_0) + h_S*phi_S*(KD(state,2) - KD(state,3)); //2-3d
  end if;

  if state == 1 then
    0 = lambda_2; //1e
  else
    0 = z_dot + r_S*phi_S_dot; //2-3e
  end if;
  
  
  when lambda_1 > 0 then
    if state == 2 or (state == 3 and phi_B_dot < 0) then
     //state := 1;
     //sleeveBlocking(z_dot,phi_S_dot,phi_B_dot);
     // sleeve blocking
         
    reinit(z_dot,0);
    reinit(phi_S_dot,0);
    reinit(phi_B_dot, mom_pre/(J_B+m_B*l_G^2));
    end if;
  end when;
  
    when lambda_1 < 0 then
      if state == 2 or (state == 3 and phi_B_dot < 0) then
     //state := 1;
     //sleeveBlocking(z_dot,phi_S_dot,phi_B_dot);
     // sleeve blocking
        //mom_pre :=momentum(z_dot,phi_S_dot,phi_B_dot);
        reinit(z_dot,0);
        reinit(phi_S_dot,0);
        reinit(phi_B_dot, mom_pre/(J_B+m_B*l_G^2));
      else
      end if;
    end when;
    
    
      when h_B*phi_B > l_S + l_G - l_B - r_0 then
    if state == 3 and phi_B_dot > 0 then
      //state := 4;
     // sleeve blocking
      //mom_pre :=momentum(z_dot,phi_S_dot,phi_B_dot);
      reinit(z_dot,0);
      reinit(phi_S_dot,0);
      reinit(phi_B_dot, mom_pre/(J_B+m_B*l_G^2));
    else
    end if;
  end when;
    
algorithm 
  
  mom_pre :=momentum(z_dot,phi_S_dot,phi_B_dot);
  
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
     //sleeveBlocking(z_dot,phi_S_dot,phi_B_dot);
     // sleeve blocking
//         mom_pre :=momentum(z_dot,phi_S_dot,phi_B_dot);
//    reinit(z_dot,0);
 //   reinit(phi_S_dot,0);
//    reinit(phi_B_dot, mom_pre/(J_B+m_B*l_G^2));
    end if;
  end when;
  
    when lambda_1 < 0 then
    if state == 2 or (state == 3 and phi_B_dot < 0) then
     state := 1;
 //        mom_pre :=momentum(z_dot,phi_S_dot,phi_B_dot);
 //   reinit(z_dot,0);
  //  reinit(phi_S_dot,0);
  //  reinit(phi_B_dot, mom_pre/(J_B+m_B*l_G^2));
    end if;
  end when;


  when h_B*phi_B > l_S + l_G - l_B - r_0 then
    if state == 3 and phi_B_dot > 0 then
      state := 4;
      //sleeveBlocking(z_dot,phi_S_dot,phi_B_dot);
    end if;
  end when;

  when state == 4 then
    state := 3;
  end when;



end Woodpecker;
