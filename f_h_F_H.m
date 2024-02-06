    clc
    clearvars
    close all

    N = 12;                                                                 % custom number of submodules (SMs) in one arm.
    syms('v_pv',[N 1])                                                      % vector of SM voltages in one arm [V]
    syms('i_pv',[N 1])                                                      % vector of SM PV currents [A]
    syms('duty',[N 1])                                                      % vector of SM duty cycles [/]
    syms('gate',[N 1])                                                      % vector of SM gate signals [/]
    syms('G',[N 1])                                                         % vector of SM irradiances (in p.u. w.r.t. 1000 W/m^2)
    syms i_arm                                                              % MMC arm current [A]
    syms dt                                                                 % Integration time step [s]
    syms Csm                                                                % Submodule capacitances [F]
    

    % Expressions of functions f and g as they appear in (3) and (5),
    % without replacing the specific formula of i_pv (4).
    f = [v_pv + dt/Csm*(i_arm.*duty + i_pv); G]
    h = gate.'*v_pv

    %%

    % We now replace i_pv with its specific formula (4). 
    syms('Rsh',[N 1])
    syms('iph',[N 1])
    syms('Rs')
    syms('vT',[N 1])
    syms('i0',[N 1])

    A    = (Rs*Rsh.*i0)./((Rsh+Rs).*vT);
    B    = (v_pv.*Rsh + Rsh*Rs.*(i0+iph))./((Rsh+Rs).*vT);
    i_pv = (Rsh.*(iph+i0)-v_pv)./(Rsh + Rs) - ...
            vT/Rs.*lambertw(A.*exp(B));

    % Expressions of functions f and g as they appear in (3) and (5),
    % by replacing the specific formula of i_pv (4) but not those of Rsh,
    % Rs, iph, i0, iph and vT.
    f = [v_pv + dt/Csm*(i_arm.*duty + i_pv); G]
    h = gate.'*v_pv

    %%
    % We now replace Rsh, Rs, iph, i0, and vT with their specific
    % formulas in (1)

    syms Nser                                                               % number of PV cells in series in a module
    syms Npar                                                               % number of PV cells in parallel in a module
    syms('T',[N 1])                                                         % vector of SM temperatures [C]
    syms T_ref                                                              % reference temperature [K]
    syms Rs_ref 
    syms Rsh_ref
    syms i0_ref
    syms vT_ref
    syms Eg_ref
    syms alpha_isc                                                          % ?
    syms iL
    syms k1
    syms dEdT

    G_ref  = 1000;                                                          % refenrence irradiance [W/m^2]
    T_k    = T + 273.15;                                                    % vector of SM temperatures [K]
    Rs     = Rs_ref*Nser/Npar;
    Rsh    = Rsh_ref*Nser/Npar*G_ref./(G*G_ref);
    iph    =  G*G_ref.*(alpha_isc + iL*(T_k - T_ref));
    Eg     = Eg_ref*((T_k-T_ref)*dEdT + 1);
    i0     = i0_ref*Npar*(T_k/T_ref).^3.*exp(Eg_ref/(k1*T_ref) - Eg/(k1*T_ref));
    vT     = vT_ref*Nser*T_k/T_ref;
    
    A    = (Rs*Rsh.*i0)./((Rsh+Rs).*vT);
    B    = (v_pv.*Rsh + Rsh*Rs.*(i0+iph))./((Rsh+Rs).*vT);
    i_pv = (Rsh.*(iph+i0)-v_pv)./(Rsh + Rs) - ...
            vT/Rs.*lambertw(A.*exp(B));

    % Expressions of functions f and g as they appear in (3) and (5),
    % by replacing the specific formula of i_pv (4) but not those of Rsh,
    % Rs, iph, i0, iph and vT.
    f = simplify([v_pv + dt/Csm*(i_arm.*duty + i_pv); G])
    h = gate.'*v_pv


    % Jacobian matrices F and H. Not visible entries due to excessive
    % length can be inspected individually
    F = simplify(jacobian(f,[v_pv; G]))
    H = simplify(jacobian(h,[v_pv; G]))

    matlabFunction(f,F,'File','prediction_EKF',...
        'Vars',{[v_pv;G],T,gate,duty,dt,i_arm,[i_arm dt Csm Nser ...
        Npar T_ref Rs_ref Rsh_ref i0_ref vT_ref Eg_ref alpha_isc iL k1 dEdT]});

    matlabFunction(h,H,'File','correction_EKF',...
        'Vars',{[v_pv;G],T,gate,duty,dt,i_arm,[i_arm dt Csm Nser ...
        Npar T_ref Rs_ref Rsh_ref i0_ref vT_ref Eg_ref alpha_isc iL k1 dEdT]});



   
