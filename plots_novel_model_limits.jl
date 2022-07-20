import Pkg;
Pkg.add("Plots");
Pkg.add("Printf");
Pkg.add("LinearAlgebra");
Pkg.add("SpecialFunctions");
Pkg.add("Interpolations");
Pkg.add("ForwardDiff");
Pkg.add("LaTeXStrings");
Pkg.add("LinearSolve");
Pkg.add("Trapz");
using Trapz;
using Plots;
using Printf;
using LinearAlgebra;
using SpecialFunctions;
using Interpolations;
using ForwardDiff;
using Statistics;
using LinearSolve;

default(fontfamily="Computer Modern");

include("projection.jl");

function describe_domain(x_I, n_points, L)
    Ω = collect(range(0, L, n_points - 1));
    idx = n_points;
    for i in 1 : n_points - 1
        if x_I < Ω[i]
            idx = i;
            break;
        end
    end
    insert!(Ω, idx, x_I);
    Ω[idx + 1] = Ω[idx    ] + (Ω[idx + 2] - Ω[idx    ]) / 2.0;
    Ω[idx - 1] = Ω[idx - 2] + (Ω[idx    ] - Ω[idx - 2]) / 2.0;
    return Ω, idx;
end

function central_FD_second(h_1, h_2, u_m, u, u_p)
    return 2 / (h_1 + h_2) * ((u_p - u) / h_2 - (u - u_m) / h_1);
end

function upwind_first(h, u, u_p, u_pp)
    return (-u_pp + 4*u_p - 3*u) / (2*h);
end

function downwind_first(h, u_mm, u_m, u)
    return (3*u - 4*u_m + u_mm) / (2*h);
end

function p_sat_compr(T, p::Dict)
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=4&Type=ANTOINE&Plot=on
    # Antoine Equations for 155 <= T <= 195
    A = 6.81228;
    B = 1301.679;
    C = -3.494;
    return 10^(A-B/(C+T))*100000;
end

function residual_compr(U, p::Dict, s::Dict, Ω, index_I, Uᶜ, Uʰ, θ, dt, t)
    y  = zero(U);

    x_I_cur = Ω[index_I] + U[index_I+1] * dt;
    pʰ = ((s["L"] - s["x_I_0"]) * s["ρ_h_0"] - (x_I_cur - s["x_I_0"]) * p["rho"]) / (s["L"] - x_I_cur) * p["R"] * 1 / (trapz(Ω[index_I:end], U[index_I+2:end].^(-1))) * (Ω[end]-Ω[index_I]);
    rho_h = pʰ / (p["R"] * 1 / (trapz(Ω[index_I:end], U[index_I+2:end].^(-1))) * (Ω[end]-Ω[index_I]));

    #aʰ = p["a_h"];
    # Trivial IC
    #aʰ = p["k_h"] / (p["rho"] * p["c_h"]);
    # Advanced IC
    aʰ  = p["k_h"] / (s["ρ_h_0"] * p["c_h"]);

    # Cold BC
    if (true)
        # Robin type
        hₛ = Ω[2] - Ω[1];
        F_1 = 1 + θ * dt * p["a_c"] * 2 * (1 + s["C"] * hₛ) / hₛ^2;
        F_1_p = 1 - 2 * dt * (1-θ) * p["a_c"] * (1 / hₛ^2 + s["C"] / hₛ);
        F_2 = 2 * dt * θ * p["a_c"] / hₛ^2;
        F_2_p = 2 * dt * (1-θ) * p["a_c"] / hₛ^2;
        y[1] = U[1] * F_1 - U[2] * F_2 - Uᶜ[1] * F_1_p - Uᶜ[2] * F_2_p - 2 * dt * p["a_c"] * s["C"] * s["T_ext"] / hₛ;
    else
        # Dirichlet
        y[1] = U[1] - s["T_ext"];
    end

    # Intermediate cold
    for i in 2 : index_I - 1
        h₁ = Ω[i]     - Ω[i - 1];
        h₂ = Ω[i + 1] - Ω[i];
        FD_cur  = central_FD_second(h₁, h₂, U[i-1], U[i], U[i+1]);
        FD_prev = central_FD_second(h₁, h₂, Uᶜ[i-1], Uᶜ[i], Uᶜ[i+1]);
        y[i] = U[i] - dt * θ * p["a_c"] * FD_cur - Uᶜ[i] + dt * (1-θ) * p["a_c"] * FD_prev;
    end

    # Intermediate hot
    for i in index_I + 3 : s["n_points"] + 1
        h₁ = Ω[i - 2] - Ω[i - 3];
        h₂ = Ω[i - 1] - Ω[i - 2];
        FD_cur  = central_FD_second(h₁, h₂, U[i-1], U[i], U[i+1]);
        FD_prev = central_FD_second(h₁, h₂, Uʰ[i - index_I - 2], Uʰ[i - index_I - 1], Uʰ[i - index_I]);
        y[i] = U[i] - dt * θ * aʰ * FD_cur - Uʰ[i - index_I - 1] + dt * (1-θ) * aʰ * FD_prev;
    end

    #T_ext = (1/4*sin(2*π*t/3600) + 1)* s["T_ext_h"];
    if (true)
        T_ext = s["T_ext_h"];
    else
        T_ext = (abs(1/8*sin(2*π*t/3600)) + 1)* s["T_ext_h"];
    end 

    if (true)
        hₛ = Ω[2] - Ω[1];
        F_end = 1 + dt * aʰ * 2 * (1 + s["C_h"] * hₛ) / hₛ^2;
        F_end_p = 1;
        F_end_1 = 2 * dt * aʰ / hₛ^2;
        y[end] = U[end] * F_end - U[end-1] * F_end_1 - Uʰ[end] * F_end_p - 2 * dt * aʰ * s["C_h"] * T_ext / hₛ;
    else
        y[s["n_points"] + 2] = U[s["n_points"] + 2] - T_ext;
    end
    # Hot BC (Dirichlet)
    #y[s["n_points"] + 2] = U[s["n_points"] + 2] - s["T_1"];

    # Hot BC (Robin type)
    #F_end = 1 + dt * aʰ * 2 * (1 + s["C_h"] * hₛ) / hₛ^2;
    #F_end_p = 1;
    #F_end_1 = 2 * dt * aʰ / hₛ^2;
    #y[end] = U[end] * F_end - U[end-1] * F_end_1 - Uʰ[end] * F_end_p - 2 * dt * aʰ * s["C_h"] * s["T_ext_h"] / hₛ;


    # Helpful definitions
    h₁ = Ω[index_I]     - Ω[index_I - 1];
    h₂ = Ω[index_I + 1] - Ω[index_I];
    FD_c = downwind_first(h₁, U[index_I - 2], U[index_I - 1], U[index_I    ]);
    FD_h =   upwind_first(h₂, U[index_I + 2], U[index_I + 3], U[index_I + 4]);
    qʰ = -FD_h * p["k_h"];
    Tʰ = U[index_I+2];
    vᴵ = U[index_I+1];
    Tᶜ = U[index_I];
    R  = p["R"];
    ρᶜ = p["rho"];
    pᶜ = p_sat_compr(Tᶜ, p);
    γ = p["gamma"];
    ϑ = p["theta"];
    # Stefan Condition
    y[index_I + 1] = ρᶜ * p["H"] * vᴵ - p["k_c"] * FD_c + p["k_h"] * FD_h;


    # Trivial interface conditions
    #y[index_I]     = U[index_I] - p["T_f"];
    #y[index_I + 2] = U[index_I] - U[index_I + 2];

    # Interface conditions
    #y[index_I] = qʰ - 2*ϑ/(2-ϑ)*sqrt(2/π) * (pᶜ * √(p["R"] * U[index_I]) - pʰ * √(p["R"] * U[index_I + 2]));
    #y[index_I] = -q_h / (p["R"]*Tʰ) + 2*√(2/π)*(1/4*ϑ/(2-ϑ)*(4*Tᶜ-5*Tʰ)/Tʰ*(pᶜ/√(p["R"]*Tᶜ) - pʰ/√(p["R"]*Tʰ)) - pʰ/√(p["R"]Tʰ)*(Tʰ-Tᶜ)/Tʰ);
    #y[index_I+2] = p["rho"] * vᴵ + ϑ/(2-ϑ)*√(2/π) * (pᶜ / √(p["R"]*Tᶜ) - pʰ / √(p["R"]*Tʰ));
    # With gamma=1
    #y[index_I]   = qʰ/(R*Tʰ) - ϑ/(2-ϑ) * (4*Tᶜ-5*Tʰ)/Tʰ * (pᶜ/√(2*π*R*Tᶜ) - pʰ/√(2*π*R*Tʰ)) + 4*pʰ/√(2*π*R*Tʰ) * (Tʰ-Tᶜ)/Tʰ;
    # With speculiar reflection coefficient
    y[index_I]   = qʰ/(R*Tʰ) - ϑ/(2-ϑ) * (4*Tᶜ-5*Tʰ)/Tʰ * (pᶜ/√(2*π*R*Tᶜ) - pʰ/√(2*π*R*Tʰ)) + 4*pʰ/√(2*π*R*Tʰ) * (Tʰ-Tᶜ)/Tʰ * (2/(γ*ϑ-γ-ϑ+2)-1);
    y[index_I+2] = ρᶜ*vᴵ + 2*ϑ/(2-ϑ) * (pᶜ/√(2*π*R*Tᶜ) - pʰ/√(2*π*R*Tʰ));
    return y;
end

function Θ_init(Θ, p::Dict, s::Dict)
    # Θ = [xᴵ, ∂ₓTᶜ, ∂ₓTʰ, Tᶜ, Tʰ, T₁, T₀]
    y = zero(Θ)
    #γ = p["gamma"];
    #ϑ = p["theta"];
    # USE ONLY FOR KINETIC COMPARISON
    γ = 1;
    ϑ = 1;
    if (ϑ>0)
        y[1] = Θ[6] - s["T_1"];
        y[2] = Θ[7] - s["T_0_init"];
        y[3] = Θ[2] * p["k_c"] - Θ[3] * p["k_h"];
        y[4] = Θ[2] - (Θ[4] - Θ[7]) / Θ[1];
        y[5] = Θ[3] - (Θ[6] - Θ[5]) / (s["L"] - Θ[1]);
        y[6] = (s["p_0"] / √Θ[5] - p_sat_compr(Θ[4], p) / √Θ[4]);
        #y[7] = Θ[3] - 2 / p["k_h"] * √(2/π) * s["p_0"] * √p["R"] * (Θ[5] - Θ[4]) / √Θ[5];
        y[7] = Θ[3] - 2 / p["k_h"] * √(2/π) * s["p_0"] * (2/(γ*ϑ-γ-ϑ+2)-1) * √p["R"] * (Θ[5] - Θ[4]) / √Θ[5];
    else
        y[1] = Θ[6] - s["T_1"];
        y[2] = Θ[7] - s["T_0_init"];
        y[3] = Θ[2] * p["k_c"] - Θ[3] * p["k_h"];
        y[4] = Θ[2] - (Θ[4] - Θ[7]) / Θ[1];
        y[5] = Θ[3] - (Θ[6] - Θ[5]) / (s["L"] - Θ[1]);
        y[6] = Θ[1]-s["L"]/2;
        #y[7] = Θ[3] - 2 / p["k_h"] * √(2/π) * s["p_0"] * √p["R"] * (Θ[5] - Θ[4]) / √Θ[5];
        y[7] = Θ[3] - 2 / p["k_h"] * √(2/π) * s["p_0"] * (2/(γ*ϑ-γ-ϑ+2)-1) * √p["R"] * (Θ[5] - Θ[4]) / √Θ[5];
    end
    return y;
end

function Θ_init_robin_pressure(Θ, p::Dict, s::Dict, xᴵ_0, T_0)
    # Θ = [xᴵ, ∂ₓTᶜ, ∂ₓTʰ, Tᶜ, Tʰ, T₁, T₀]
    #      1 ,   2 ,  3  , 4 ,  5,  6,  7
    y = zero(Θ);
    T̃ʰ = (s["L"] - Θ[1]) / (1/Θ[3] * log(Θ[6]/(Θ[3]*(Θ[1]-s["L"]) + Θ[6])));
    pʰ = ((s["L"] - xᴵ_0) * s["p_0"] / (p["R"] * T_0) - (Θ[1] - xᴵ_0) * p["rho"]) / (s["L"] - Θ[1]) * p["R"] * T̃ʰ;
    γ = p["gamma"];
    ϑ = p["theta"];
    if (ϑ>0)
        y[1] = Θ[3] - s["C_h"] * (s["T_ext_h"] - Θ[6]);
        y[2] = Θ[2] - s["C"] * (Θ[7] - s["T_ext"]);
        y[3] = Θ[2] * p["k_c"] - Θ[3] * p["k_h"];
        y[4] = Θ[2] - (Θ[4] - Θ[7]) / Θ[1];
        y[5] = Θ[3] - (Θ[6] - Θ[5]) / (s["L"] - Θ[1]);
        y[6] = (pʰ / √Θ[5] - p_sat_compr(Θ[4], p) / √Θ[4]);
        #y[7] = Θ[3] - 2 / p["k_h"] * √(2/π) * pʰ * √p["R"] * (Θ[5] - Θ[4]) / √Θ[5];
        y[7] = Θ[3] - 2 / p["k_h"] * √(2/π) * pʰ * (2/(γ*ϑ-γ-ϑ+2)-1) * √p["R"] * (Θ[5] - Θ[4]) / √Θ[5];
    else
        y[1] = Θ[3] - s["C_h"] * (s["T_ext_h"] - Θ[6]);
        y[2] = Θ[2] - s["C"] * (Θ[7] - s["T_ext"]);
        y[3] = Θ[2] * p["k_c"] - Θ[3] * p["k_h"];
        y[4] = Θ[2] - (Θ[4] - Θ[7]) / Θ[1];
        y[5] = Θ[3] - (Θ[6] - Θ[5]) / (s["L"] - Θ[1]);
        y[6] = Θ[1]-xᴵ_0;
        y[7] = Θ[3] - 2 / p["k_h"] * √(2/π) * pʰ * (2/(γ*ϑ-γ-ϑ+2)-1) * √p["R"] * (Θ[5] - Θ[4]) / √Θ[5];
    end
    return y;
end

function non_uniform_scale_theta_adv_bc_newton(p::Dict, s::Dict)
    # Calculate initial condition

    # First, good starting values for Newtons method by assuming trivial IC
      xᴵ_init = s["L"] / (1 + p["k_h"] / p["k_c"] * (s["T_1"] - p["T_f"]) / (p["T_f"] - s["T_0_init"]));
    ∂ₓTᶜ_init = (p["T_f"] - s["T_0_init"]) / xᴵ_init;
    ∂ₓTʰ_init = (s["T_1"] - p["T_f"]) / (s["L"] - xᴵ_init);
      Tᶜ_init = p["T_f"];
      Tʰ_init = p["T_f"];
      T₁_init = s["T_1"];
      T₀_init = s["T_0_init"];
    Θ = [xᴵ_init, ∂ₓTᶜ_init, ∂ₓTʰ_init, Tᶜ_init, Tʰ_init, T₁_init, T₀_init];

    # Newton's Method for actual starting values
    res = norm(Θ_init(Θ, p, s), 2);
    n_newton = 0;
    while (res > 1.0e-8 && n_newton < 100)
        A = ForwardDiff.jacobian(Θ -> Θ_init(Θ, p, s), Θ);
        b = -Θ_init(Θ, p, s);
        dΘ = A\b;
        λ = 2;
        damped_res = typemax(Float64);
        for l in 1 : 10
            λ /= 2;
            damped_res = norm(Θ_init(Θ + λ*dΘ, p, s), 2);
            if damped_res < res
                break;
            end
        end
        Θ += λ*dΘ;
        res = damped_res;
        n_newton += 1;
        println("  Newton steps: ", n_newton, ", residual = ", res);
    end
    xᴵ_init, ∂ₓTᶜ_init, ∂ₓTʰ_init, Tᶜ_init, Tʰ_init, T₁_init, T₀_init = Θ;
    println("xᴵ_init, ∂ₓTᶜ_init, ∂ₓTʰ_init, Tᶜ_init, Tʰ_init, T₁_init, T₀_init: (starting values)");
    println(Θ);
    Ω, index_I = describe_domain(xᴵ_init, s["n_points"], s["L"]);
    


    U_c = ∂ₓTᶜ_init .* Ω[1:index_I]   .+ T₀_init;
    U_h = ∂ₓTʰ_init .* Ω[index_I:end] .+ (Tʰ_init*s["L"] - T₁_init*xᴵ_init) / (s["L"] - xᴵ_init);
    xᴵ_max = xᴵ_init + s["p_0"]/(1 / (trapz(Ω[index_I:end], U_h.^(-1))) * (Ω[end]-Ω[index_I])*p["R"]*p["rho"])*(s["L"]-xᴵ_init);
    println(xᴵ_max);

    U_c_init = U_c;
    U_h_init = U_h;
    plot_init_c = Ω[1:index_I];
    plot_init_h = Ω[index_I:end];

    # Stationary result for whole system
    xᴵ_0 = xᴵ_init;
    T_0 = (s["L"] - xᴵ_init) / (1/∂ₓTʰ_init * log(T₁_init/(∂ₓTʰ_init*(xᴵ_init-s["L"]) + T₁_init)));

    res = norm(Θ_init_robin_pressure(Θ, p, s, xᴵ_0, T_0), 2);
    n_newton = 0;
    while (res > 1.0e-7 && n_newton < 0)
        A = ForwardDiff.jacobian(Θ -> Θ_init_robin_pressure(Θ, p, s, xᴵ_0, T_0), Θ);
        b = -Θ_init_robin_pressure(Θ, p, s, xᴵ_0, T_0);
        dΘ = A\b;

        λ = 2;
        damped_res = typemax(Float64);
        for l in 1 : 5
            λ /= 2;
            damped_res = norm(Θ_init_robin_pressure(Θ + λ*dΘ, p, s, xᴵ_0, T_0), 2);
            if damped_res < res
                break;
            end
        end
        Θ += λ*dΘ;
        res = damped_res;
        n_newton += 1;
        println(Θ);
    end

    plot_stat_c = collect(range(0, Θ[1], 300));
    plot_stat_h = collect(range(Θ[1], s["L"], 300));
    U_c_stat = Θ[2] .* plot_stat_c   .+ Θ[7];
    U_h_stat = Θ[3] .* plot_stat_h .+ (Θ[5]*s["L"] - Θ[6]*Θ[1]) / (s["L"] - Θ[1]);

    # xᴵ_init, ∂ₓTᶜ_init, ∂ₓTʰ_init, Tᶜ_init, Tʰ_init, T₁_init, T₀_init = Θ;
    #println("xᴵ_init, ∂ₓTᶜ_init, ∂ₓTʰ_init, Tᶜ_init, Tʰ_init, T₁_init, T₀_init: (stationary result)");
    #println(Θ);
    p_h_tmp = ((s["L"] - xᴵ_0) * s["p_0"] / (p["R"] * T_0) - (Θ[1] - xᴵ_0) * p["rho"]) / (s["L"] - Θ[1]) * p["R"] * (s["L"] - Θ[1]) / (1/Θ[3] * log(Θ[6]/(Θ[3]*(Θ[1]-s["L"]) + Θ[6])));
    # s["p_0"] = p_h_tmp;
    println("Pressure at stationary limit: ", p_h_tmp)

    # Describe domain
    Ω, index_I = describe_domain(xᴵ_init, s["n_points"], s["L"]);

    U_c = ∂ₓTᶜ_init .* Ω[1:index_I]   .+ T₀_init;
    U_h = ∂ₓTʰ_init .* Ω[index_I:end] .+ (Tʰ_init*s["L"] - T₁_init*xᴵ_init) / (s["L"] - xᴵ_init);

    xᴵ_max = s["p_0"]/(1 / (trapz(Ω[index_I:end], U_h.^(-1))) * (Ω[end]-Ω[index_I])*p["R"]*p["rho"])*(s["L"]-xᴵ_init);
    println("Maximum ice growth: ", xᴵ_max);
    println("Maximum interface position: ", xᴵ_max + xᴵ_init);
    maximum_interface = xᴵ_max + xᴵ_init;
    println("Initial interface at x = ", xᴵ_init, " with index ", index_I, " from ", s["n_points"]);

    # Solving parameters
    dt = (s["T"][2] - s["T"][1]) / s["n_timesteps"];
    println("Solving from ", s["T"][1], " to ", s["T"][2],
            " in ", s["n_timesteps"], " steps (dt = ", dt, ")");

    x_I = zeros(s["n_timesteps"] + 1);
    Δp = zeros(s["n_timesteps"] + 1);
    ΔT = zeros(s["n_timesteps"] + 1);
    p_gas = zeros(s["n_timesteps"] + 1);
    Kn1 = zeros(s["n_timesteps"] + 1);
    Kn2 = zeros(s["n_timesteps"] + 1);
    x_I[1] = xᴵ_init;

    color_init = 1;
    θ = 1;
    U = vcat(U_c, xᴵ_max / dt, U_h);
    # Compressibility parameters
    s["x_I_0"] = x_I[1];
    s["ρ_h_0"] = s["p_0"] / (p["R"] * 1 / (trapz(Ω[index_I:end], U_h.^(-1))) * (Ω[end]-Ω[index_I]));
    p_H = ((s["L"] - s["x_I_0"]) * s["ρ_h_0"] - (Ω[index_I] - s["x_I_0"]) * p["rho"]) / (s["L"] - Ω[index_I]) * p["R"] * 1 / (trapz(Ω[index_I:end], U_h.^(-1))) * (Ω[end]-Ω[index_I]);
    Δp[1] = (p_H - p_sat_compr(U_c[end], p)) / p_H;
    ΔT[1] = (U_h[1] - U_c[end])/U_h[1];
    p_gas[1] = p_H;
    temp_plot = plot();
    Kn1[1] = (1.380649e-3 * U_h[1])/(sqrt(2)*π*3.34^2*abs(p_H)*(s["L"]-x_I[1]));
    #Kn2[1] = 3 * p["k_h"] / (p["c_h"] * sqrt(1/(1-p["R"]/p["c_h"]) * abs(p_H)/s["ρ_h_0"])) / (s["L"]-x_I[1]);
    #Kn2[1] = 3 * p["k_h"] * sqrt(p["R"] * U_h[1]) / (p["c_h"] * abs(p_H) * sqrt(1/(1-p["R"]/p["c_h"]))) / (s["L"]-x_I[1]);
    Kn2[1] = 1e-6/abs(p_H)*sqrt((p["R"]*U_h[1]*π)/2) / (s["L"]-x_I[1]);
    for k in 1 : s["n_timesteps"]
        t = s["T"][1] + k * dt;
        n_newton_max = 10;
        res = norm(residual_compr(U, p, s, Ω, index_I, U_c, U_h, θ, dt, t), 2);
        n_newton = 0;


        while (res > 1.0e-7 && n_newton < n_newton_max)
            A = ForwardDiff.jacobian(U -> residual_compr(U, p, s, Ω, index_I, U_c, U_h, θ, dt, t), U);
            #println("Condition number of matrix (with ", s["n_points"], " points and h = ", s["L"]/(s["n_points"]-1) ,"): ", cond(A));
            b = -residual_compr(U, p, s, Ω, index_I, U_c, U_h, θ, dt, t);
            dU = A\b;
            λ = 2;
            damped_res = typemax(Float64);
            for l in 1 : 5
                λ /= 2;
                damped_res = norm(residual_compr(U + λ*dU, p, s, Ω, index_I, U_c, U_h, θ, dt, t), 2);
                if damped_res < res
                    break;
                end
            end
            U += λ*dU;
            res = damped_res;
            n_newton += 1;
        end
        println("Iteration ", k, "/", s["n_timesteps"]);
        println("  Newton steps: ", n_newton, ", residual = ", res);
        Ω_prev = Ω;
        index_I_prev = index_I;

        #press_h = ((s["L"] - s["x_I_0"]) * s["ρ_h_0"] - (Ω[index_I] - s["x_I_0"]) * p["rho"]) / (s["L"] - Ω[index_I]) * p["R"] * mean(U[index_I+2:end]);
        press_h = ((s["L"] - s["x_I_0"]) * s["ρ_h_0"] - (Ω[index_I] - s["x_I_0"]) * p["rho"]) / (s["L"] - Ω[index_I]) * p["R"] * 1 / (trapz(Ω[index_I:end], U[index_I+2:end].^(-1))) * (Ω[end]-Ω[index_I]);

        x_I[k+1] = x_I[k] + dt * U[index_I + 1];
        Δp[k+1]  = (press_h - p_sat_compr(U[index_I], p))/press_h;
        ΔT[k+1]  = (U[index_I+2]-U[index_I])/U[index_I+2];
        p_gas[k+1] = press_h;
        Kn1[k+1] = (1.380649e-3 * U[index_I+2]) / (sqrt(2)*π*3.34^2*abs(press_h) * (s["L"]-x_I[k+1]));
        #Kn2[k+1] = 3 * p["k_h"] / (p["c_h"] * sqrt(1/(1-p["R"]/p["c_h"]) * abs(press_h)/s["ρ_h_0"])) / (s["L"]-x_I[1]);
        #Kn2[k+1] = 3 * p["k_h"] * sqrt(p["R"] * U[index_I+2]) / (p["c_h"] * abs(press_h) * sqrt(1/(1-p["R"]/p["c_h"]))) / (s["L"]-x_I[1]);
        Kn2[k+1] = (1e-6)/abs(press_h)*sqrt((p["R"]*U[index_I+2]*π)/2) / (s["L"]-x_I[k+1]);
        println("  -- Interface --")
        println("    (Tʰ-Tᶜ)/Tʰ = ", ΔT[k+1], ", Tʰ = ", U[index_I+2], " Tᶜ = ", U[index_I]);
        println("    (pʰ-pᶜ)/pʰ = ", Δp[k+1], ", pʰ = ", press_h, " pᶜ = ", p_sat_compr(U[index_I], p));
        println("    vᴵ = ", U[index_I+1]);
        println("  ---------------")

        Ω, index_I = describe_domain(x_I[k+1], s["n_points"], s["L"]);

        U_c = l2_projection(Ω_prev[1:index_I_prev], Ω[1:index_I], U[1:index_I_prev]);
        U_h = l2_projection(Ω_prev[index_I_prev:end], Ω[index_I:end], U[index_I_prev+2:end]);

        
        if (k==10 || k==80 || k==200 || k==250 || k==500)
            dist = 1.4*0.0036;
            dist_y = -2;
            if (k == 10)
                plot!(
                [Ω_prev[1:index_I_prev],
                Ω_prev[index_I_prev:end]],
                [U[1:index_I_prev],
                U[index_I_prev+2:end]],
                color = [color_init color_init],
                label = ["" ""],
                legend = :bottomright);
                annotate!(0.06+dist, 163.5+dist_y, L"1", annotationfontsize = 10);
                ylabel!(L"Temperature $T$ (K)");
                xlabel!(L"Position $x$ (m)");
                color_init = color_init + 1;
                #plot!([x_I[k+1], x_I[k+1]], [U[index_I_prev], U[index_I_prev+2]], label = "");
            else
                    if (k == 500)
                        plot!(
                            [Ω_prev[1:index_I_prev],
                            Ω_prev[index_I_prev:end]],
                            [U[1:index_I_prev],
                            U[index_I_prev+2:end]],
                            color = [:goldenrod2 :goldenrod2],
                            label = "",
                            legend = :bottomright,
                            linestyle = [:dash :dash]);

                    else
                        plot!(
                            [Ω_prev[1:index_I_prev],
                            Ω_prev[index_I_prev:end]],
                            [U[1:index_I_prev],
                            U[index_I_prev+2:end]],
                            color = [color_init color_init],
                            label = "",
                            legend = :bottomright);
                    end
                    color_init = color_init + 1;
            end
            if (k == 500)
                plot!([plot_init_c, plot_init_h], [U_c_init, U_h_init],
                    color = :darkgrey,
                    linestyle = [:solid :solid],
                    label = ["initial condition" ""])
                annotate!(0.06+dist, 106+dist_y, L"2", annotationfontsize = 10);
                annotate!(0.06+dist, 83+dist_y, L"3", annotationfontsize = 10);
                annotate!(0.06+dist, 57+dist_y, L"4", annotationfontsize = 10);
                annotate!(0.06+dist, 22.5+dist_y, L"5", annotationfontsize = 10);
                lens!([0.1, 0.11], [184.5, 185.1], inset = (1, bbox(0.05, 0.3, 0.3, 0.4, :top, :right)),
                    subplot=2,
                    framestyle=:box,
                    xticks = ([0.1025, 0.105, 0.1075], ["0.1025", "0.105", "0.1075"]),
                    yticks = [184.6, 184.7, 184.8, 184.9, 185.0],
                    ls=:dot,
                    lc=:black)
            end
        end

        U = vcat(U_c, U[index_I_prev+1], U_h);
        println();
        println(".........................................");
        println();
    end
    display(temp_plot); savefig(temp_plot,"temperature_plot_novel_co2_limit_case.pdf");
    if (false)
        interface_plot = plot(s["T"][1]:dt:s["T"][2], x_I,
                                ylabel = "Interface position in [m]",
                                xlabel = "Time in [s]",
                                title = "Interface propagation",
                                legend = :bottomright);
        #savefig(propagation_plot,"non_uniform.pdf");
        pressure_plot = plot(s["T"][1]:dt:s["T"][2], Δp,
                                ylabel = "Pressure jump (pʰ-pᶜ)/pʰ",
                                xlabel = "Time in [s]",
                                title = "Pressure difference",
                                legend = :bottomright);
        temperature_plot = plot(s["T"][1]:dt:s["T"][2], ΔT,
                                ylabel = "Temperature jump (Tʰ-Tᶜ)/Tʰ",
                                xlabel = "Time in [s]",
                                title = "Temperature difference",
                                legend = :bottomright);
        display(interface_plot);
        display(pressure_plot);
        display(temperature_plot);
        savefig(interface_plot,"interface_plot3.pdf");
        savefig(pressure_plot,"pressure_plot3.pdf");
        savefig(temperature_plot,"temperature_plot3.pdf");
        return x_I;
    end
    if (true)
        # Plot: Interface with marked time points
        ticks = [(i-1) * 3600 for i = 1 : 2:31]
        labels = [(i-1) for i = 1 :2: 31]
        interface_plot = plot(
            [s["T"][1]:dt:s["T"][2]], [x_I],
            ylabel = L"Interface Position $x_I$ (m)",
            xlabel = L"Time $t$ (h)",
            legend = :none,
            xticks = (ticks, labels),
            color = :blue,
            linestyle = :solid,
        );
        dist = 3*900;
        dist2 = 0.000001
        iI = LinearInterpolation(collect(s["T"][1]:dt:s["T"][2]), x_I);
        scatter!([2160], [iI(2160)], markershape = :rect, color = :black, label = "");
        annotate!(2160+dist, iI(2160)-dist2, L"1", annotationfontsize = 10);
        scatter!([17280], [iI(17280)], markershape = :rect, color = :black, label = "");
        annotate!(17280+dist, iI(17280)-dist2, L"2", annotationfontsize = 10);
        scatter!([43200], [iI(43200)], markershape = :rect, color = :black, label = "");
        annotate!(43200+dist, iI(43200)-dist2, L"3", annotationfontsize = 10);
        scatter!([54000], [iI(54000)], markershape = :rect, color = :black, label = "");
        annotate!(54000+dist, iI(54000)-dist2, L"4", annotationfontsize = 10);   
        scatter!([108000], [iI(108000)], markershape = :rect, color = :black, label = "");
        annotate!(108000+dist, iI(108000)-dist2, L"5", annotationfontsize = 10);
        savefig(interface_plot,"interface_plot_novel_co2_limit_case.pdf");
        display(interface_plot);
        println(x_I[end]);
    end
    if (true)
        ticks = [(i-1) * 3600 for i = 1 : 2 : 31]
        labels = [(i-1) for i = 1 : 2 : 31]
        # Temperature difference at interface
        temperature_plot = plot([s["T"][1]:dt:s["T"][2]], [ΔT],
        ylabel = "\n"*L"Temperature jump $\frac{T_h-T_c}{T_h}$ at $x_I$",
        xlabel = "Time t (s)",
        legend = :none,
        color = :blue,
        linestyle = :solid,
        xticks = (ticks, labels));
        display(temperature_plot);
        savefig(temperature_plot,"novel_co2_temperature_diff_limit_case.pdf");
    end
    if (true)
        ticks = [(i-1) * 3600 for i = 1 : 2:31]
        labels = [(i-1) for i = 1 :2: 31]
        # Relative pressure in gas phase
        T̃ʰ = (s["L"] - Θ[1]) / (1/Θ[3] * log(Θ[6]/(Θ[3]*(Θ[1]-s["L"]) + Θ[6])));
        pʰ = ((s["L"] - xᴵ_init) * s["p_0"] / (p["R"] * T_0) - (Θ[1] - xᴵ_init) * p["rho"]) / (s["L"] - Θ[1]) * p["R"] * T̃ʰ;
        press_plot = plot([s["T"][1]:dt:s["T"][2]], [p_gas ./ s["p_0"]],
        ylabel = "\n"*L"Relative pressure $\frac{p_h}{p_{\mathrm{init}}}$",
        xlabel = "Time t (s)",
        color = :blue,
        legend = :none,
        linestyle = :solid,
        #ylims = (-0.05,1.05),
        xticks = (ticks, labels));
        display(press_plot);
        savefig(press_plot,"novel_co2_press_prop_limit_case.pdf");
    end
    if (true)
        ticks = [(i-1) * 3600 for i = 1 : 2 : 31]
        labels = [(i-1) for i = 1 : 2 : 31]
        # Knudsen plot at interface
        temperature_plot = plot([s["T"][1]:dt:s["T"][2]], [Kn2],
        ylabel = "\n"*L"Knudsen number $\mathit{Kn}$ at $x_I$",
        xlabel = "Time t (s)",
        legend = :bottomright,
        label = "",
        color = :blue,
        #ylims = (0, 100),
        yaxis = :log,
        yticks = ([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6], [L"10^{-7}" L"10^{-6}" L"10^{-5}" L"10^{-4}" L"10^{-3}" L"10^{-2}" L"10^{-1}" L"1" L"10^{1}" L"10^{2}" L"10^{3}" L"10^{4}" L"10^{5}" L"10^{6}"]),
        linestyle = :solid,
        xticks = (ticks, labels));
        display(temperature_plot);
        savefig(temperature_plot,"novel_co2_Kn_limit_case.pdf");
    end
    println("Used gas: ", (x_I[end]-xᴵ_init)/(xᴵ_max)*100, "%");
    temp_frac=(Θ[5]-Θ[4])/Θ[5];
    return x_I, ΔT, Θ[1], temp_frac, maximum_interface;
end

p_co2 = Dict(
    # Material parameters, water liquid <--> solid
    "H"   => 571000.0,    # [J / kg]     latent heat of sublimation
    "T_f" => 180.0,       # [K]          phase-shift temperature (at 1000 Pa)
    "rho" => 1600,      # [kg / m^3]   density #1600

    "c_h" => 777.05,      # [J / (kg K)] specific heat capacity
    "c_c" => 1239.49,     # [J / (kg K)] specific heat capacity
    "k_h" => 1,   # [W / (m K)]  thermal conductivity # 1
    "k_c" => 0.5, # [W / (m K)]  thermal conductivity # 0.025
    "R"   => 188.98,       # [J / (kg K)]  specific gas constant

    "gamma" => 1,  # Accomodation coefficient
    "theta" => 1  # Condensation coeff
);
p_co2["a_h"] = p_co2["k_h"] / (p_co2["rho"] * p_co2["c_h"]); # [m^2 / s]    thermal diffusivity
p_co2["a_c"] = p_co2["k_c"] / (p_co2["rho"] * p_co2["c_c"]); # [m^2 / s]    thermal diffusivity

s_co2 = Dict(
    "T"           => [0, 30*3600],  # [[s], [s]] Solving time
    "n_timesteps" => 500,        # [-]        number of timesteps
    "n_points"    => 100,        # [-]        number of spacial discretization points
    "L"           => 0.2,        # [m]        canal length
    "T_1"         => 185,        # [K]        fixed temperature at hot surface (BC, IC)
    "T_0_init"    => 160,        # [K]        initial temperature at cold surface (IC)
    "C"           => 100,         # [1 / m]    constant of Newton's cooling (BC)
    "T_ext"       => 20,  #20    # [K]        ambient temperature for Newton's cooling (BC)
    "C_h"         => 100,         # [1 / m]    constant of Newton's cooling on hot side
    "T_ext_h"     => 185,  #220  # [K]        ambient temperature hot side
    "p_0"         => 20000,      # [Pa]       initial gas pressure in hot phase
);

I6, ΔT6, stat_I, stat_T = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);

if (false)
    p_co2["gamma"] = 0;
    p_co2["theta"] = 0.1;
    I1, ΔT1, I1_s, ΔT1_s, maximum_interface = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    p_co2["gamma"] = 0;
    p_co2["theta"] = 0.05;
    I2, ΔT2, I2_s, ΔT2_s, maximum_interface = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    p_co2["gamma"] = 0;
    p_co2["theta"] = 0;
    I3, ΔT3, I3_s, ΔT3_s, maximum_interface = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    dt=(s_co2["T"][2]-s_co2["T"][1])/s_co2["n_timesteps"];
    ticks = [(i-1) * 3600 for i = 1 : 2:31]
    labels = [(i-1) for i = 1 :2: 31]
    lens_area = [0, 30*3600];
    comparison_plot = plot([s_co2["T"][1]:dt:s_co2["T"][2], s_co2["T"][1]:dt:s_co2["T"][2], s_co2["T"][1]:dt:s_co2["T"][2], lens_area, lens_area, lens_area], [(I1 .- I1[1])./(maximum_interface-I1[1]), (I2 .- I1[1])./(maximum_interface-I1[1]), (I3 .- I1[1])./(maximum_interface-I1[1]), ones(2)*(I1_s .- I1[1])./(maximum_interface-I1[1]), ones(2)*(I2_s .- I1[1])./(maximum_interface-I1[1]), ones(2)*(I3_s .- I1[1])./(maximum_interface-I1[1])],
    ylabel = "\n"*L"Relative interface position $\frac{x_I-x_I^{(0)}}{x_{I,\mathrm{max}}-x_I^{(0)}}$",
    xlabel = L"Time $t$ (h)",
    legend = :right,
    label = [L"\gamma=0, \vartheta=0.1" L"\gamma=0, \vartheta=0.05" L"\gamma=0, \vartheta=0" "Interface positions at stationary limit" "" ""],
    color = [:blue :orange :red :darkgrey :darkgrey :darkgrey],
    linestyle = [:solid :solid :solid :dash :dash :dash],
    xticks = (ticks, labels));
    display(comparison_plot);
    savefig(comparison_plot, "novel_comp_gamma0theta0_interface.pdf");
    dt=(s_co2["T"][2]-s_co2["T"][1])/s_co2["n_timesteps"];
    ticks = [(i-1) * 3600 for i = 1 : 2:31]
    labels = [(i-1) for i = 1 :2: 31]
    comparison_plot = plot(s_co2["T"][1]:dt:s_co2["T"][2], [ΔT1, ΔT2, ΔT3, ones(s_co2["n_timesteps"]+1)*ΔT1_s, ones(s_co2["n_timesteps"]+1)*ΔT2_s, ones(s_co2["n_timesteps"]+1)*ΔT3_s],
    ylabel = "\n"*L"Temperature jump $\frac{T_h-T_c}{T_h}$ at $x_I$",
    xlabel = L"Time $t$ (h)",
    legend = :right,
    label = [L"\gamma=0, \vartheta=0.1" L"\gamma=0, \vartheta=0.05" L"\gamma=0, \vartheta=0" "Temperature jumps at stationary limit" "" ""],
    xticks = (ticks, labels),
    color = [:blue :orange :red :darkgrey :darkgrey :darkgrey],
    linestyle = [:solid :solid :solid :dash :dash :dash]);
    display(comparison_plot);
    savefig(comparison_plot, "novel_comp_gamma0theta0_temp.pdf");
end

if (false)
    p_co2["gamma"] = 0;
    p_co2["theta"] = 1;
    I1, ΔT1, I1_s, ΔT1_s, maximum_interface = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    p_co2["gamma"] = 0;
    p_co2["theta"] = 0.5;
    I2, ΔT2, I2_s, ΔT2_s, maximum_interface = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    p_co2["gamma"] = 0;
    p_co2["theta"] = 0.1;
    I3, ΔT3, I3_s, ΔT3_s, maximum_interface = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    dt=(s_co2["T"][2]-s_co2["T"][1])/s_co2["n_timesteps"];
    ticks = [(i-1) * 3600 for i = 1 : 2:31]
    labels = [(i-1) for i = 1 :2: 31]
    lens_area = [28*3600-900, 29*3600+900];
    comparison_plot = plot([s_co2["T"][1]:dt:s_co2["T"][2], s_co2["T"][1]:dt:s_co2["T"][2], s_co2["T"][1]:dt:s_co2["T"][2], lens_area, lens_area, lens_area], [(I1 .- I1[1])./(maximum_interface-I1[1]), (I2 .- I1[1])./(maximum_interface-I1[1]), (I3 .- I1[1])./(maximum_interface-I1[1]), ones(2)*(I1_s .- I1[1])./(maximum_interface-I1[1]), ones(2)*(I2_s .- I1[1])./(maximum_interface-I1[1]), ones(2)*(I3_s .- I1[1])./(maximum_interface-I1[1])],
    ylabel = "\n"*L"Relative interface position $\frac{x_I-x_I^{(0)}}{x_{I,\mathrm{max}}-x_I^{(0)}}$",
    xlabel = L"Time $t$ (h)",
    legend = :bottomright,
    label = [L"\gamma=0, \vartheta=1" L"\gamma=0, \vartheta=0.5" L"\gamma=0, \vartheta=0.1" "" "" ""],
    color = [:blue :orange :red :blue :orange :red],
    linestyle = [:solid :solid :solid :dash :dash :dash],
    xticks = (ticks, labels));
    lens!(lens_area, [0.94125, 0.95375], inset = (1, bbox(0.1, 0.2, 0.4, 0.4, :top, :right)),
    xticks = ([28*3600, 28.5*3600, 29*3600], ["29", "29.5", "30"]),
    yticks = ([0.9425, 0.945, 0.9475, 0.950, 0.9525]),
    subplot=2,
    framestyle=:box,
    ls=:dot,
    lc=:black)
    annotate!(29*3600-1800, 0.95, "Dashed lines are stationary limits", subplot=2, annotationfontsize=8);
    display(comparison_plot);
    savefig(comparison_plot, "novel_comp_gamma0_interface.pdf");
    dt=(s_co2["T"][2]-s_co2["T"][1])/s_co2["n_timesteps"];
    ticks = [(i-1) * 3600 for i = 1 : 2:31]
    labels = [(i-1) for i = 1 :2: 31]
    comparison_plot = plot(s_co2["T"][1]:dt:s_co2["T"][2], [ΔT1, ΔT2, ΔT3, ones(s_co2["n_timesteps"]+1)*ΔT1_s, ones(s_co2["n_timesteps"]+1)*ΔT2_s, ones(s_co2["n_timesteps"]+1)*ΔT3_s],
    ylabel = "\n"*L"Temperature jump $\frac{T_h-T_c}{T_h}$ at $x_I$",
    xlabel = L"Time $t$ (h)",
    legend = :right,
    label = [L"\gamma=0, \vartheta=1" L"\gamma=0, \vartheta=0.5" L"\gamma=0, \vartheta=0.1" "Temperature jumps at stationary limit" "" ""],
    xticks = (ticks, labels),
    color = [:blue :orange :red :darkgrey :darkgrey :darkgrey],
    linestyle = [:solid :solid :solid :dash :dash :dash]);
    display(comparison_plot);
    savefig(comparison_plot, "novel_comp_gamma0_temp.pdf");
end

if (false)
    p_co2["gamma"] = 0.5;
    p_co2["theta"] = 1;
    I1, ΔT1, I1_s, ΔT1_s, maximum_interface = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    p_co2["gamma"] = 0.5;
    p_co2["theta"] = 0.5;
    I2, ΔT2, I2_s, ΔT2_s = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    p_co2["gamma"] = 0.5;
    p_co2["theta"] = 0.1;
    I3, ΔT3, I3_s, ΔT3_s = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    dt=(s_co2["T"][2]-s_co2["T"][1])/s_co2["n_timesteps"];
    ticks = [(i-1) * 3600 for i = 1 : 2:31]
    labels = [(i-1) for i = 1 :2: 31]
    lens_area = [28*3600-900, 29*3600+900];
    comparison_plot = plot([s_co2["T"][1]:dt:s_co2["T"][2], s_co2["T"][1]:dt:s_co2["T"][2], s_co2["T"][1]:dt:s_co2["T"][2], lens_area, lens_area, lens_area], [(I1 .- I1[1])./(maximum_interface-I1[1]), (I2 .- I1[1])./(maximum_interface-I1[1]), (I3 .- I1[1])./(maximum_interface-I1[1]), ones(2)*(I1_s .- I1[1])./(maximum_interface-I1[1]), ones(2)*(I2_s .- I1[1])./(maximum_interface-I1[1]), ones(2)*(I3_s .- I1[1])./(maximum_interface-I1[1])],
    ylabel = "\n"*L"Relative interface position $\frac{x_I-x_I^{(0)}}{x_{I,\mathrm{max}}-x_I^{(0)}}$",
    xlabel = L"Time $t$ (h)",
    legend = :bottomright,
    label = [L"\gamma=0.5, \vartheta=1" L"\gamma=0.5, \vartheta=0.5" L"\gamma=0.5, \vartheta=0.1" "" "" ""],
    color = [:blue :orange :red :blue :orange :red],
    linestyle = [:solid :solid :solid :dash :dash :dash],
    xticks = (ticks, labels));
    lens!(lens_area, [0.94175, 0.94375], inset = (1, bbox(0.1, 0.2, 0.4, 0.4, :top, :right)),
    xticks = ([28*3600, 28.5*3600, 29*3600], ["29", "29.5", "30"]),
    yticks = ([0.942, 0.9425, 0.943, 0.9435]),
    subplot=2,
    framestyle=:box,
    ls=:dot,
    lc=:black)
    annotate!(29*3600-1800, 0.94325, "Dashed lines are stationary limits", subplot=2, annotationfontsize=8);
    display(comparison_plot);
    savefig(comparison_plot, "novel_comp_gamma05_interface.pdf");
    dt=(s_co2["T"][2]-s_co2["T"][1])/s_co2["n_timesteps"];
    ticks = [(i-1) * 3600 for i = 1 : 2:31]
    labels = [(i-1) for i = 1 :2: 31]
    comparison_plot = plot(s_co2["T"][1]:dt:s_co2["T"][2], [ΔT1, ΔT2, ΔT3, ones(s_co2["n_timesteps"]+1)*ΔT1_s, ones(s_co2["n_timesteps"]+1)*ΔT2_s, ones(s_co2["n_timesteps"]+1)*ΔT3_s],
    ylabel = "\n"*L"Temperature jump $\frac{T_h-T_c}{T_h}$ at $x_I$",
    xlabel = L"Time $t$ (h)",
    legend = :bottomright,
    label = [L"\gamma=0.5, \vartheta=1" L"\gamma=0.5, \vartheta=0.5" L"\gamma=0.5, \vartheta=0.1" "Temperature jumps at stationary limit" "" ""],
    xticks = (ticks, labels),
    color = [:blue :orange :red :darkgrey :darkgrey :darkgrey],
    linestyle = [:solid :solid :solid :dash :dash :dash]);
    display(comparison_plot);
    savefig(comparison_plot, "novel_comp_gamma05_temp.pdf");
end

if (false)
    # Gamma=1, mehrer theta
    p_co2["gamma"] = 1;
    p_co2["theta"] = 1;
    I1, ΔT1, I1_s, ΔT1_s, maximum_interface = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    p_co2["theta"] = 1.0e-5;
    I2, ΔT2, I2_s, ΔT2_s = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    p_co2["theta"] = 1.0e-6;
    I3, ΔT3, I3_s, ΔT3_s = non_uniform_scale_theta_adv_bc_newton(p_co2, s_co2);
    dt=(s_co2["T"][2]-s_co2["T"][1])/s_co2["n_timesteps"];
    ticks = [(i-1) * 3600 for i = 1 : 2:31]
    labels = [(i-1) for i = 1 :2: 31]
    comparison_plot = plot(s_co2["T"][1]:dt:s_co2["T"][2], [(I1 .- I1[1])./(maximum_interface-I1[1]), (I2 .- I1[1])./(maximum_interface-I1[1]), (I3 .- I1[1])./(maximum_interface-I1[1]), ones(s_co2["n_timesteps"]+1)*(I1_s .- I1[1])./(maximum_interface-I1[1])],
    ylabel = "\n"*L"Relative interface position $\frac{x_I-x_I^{(0)}}{x_{I,\mathrm{max}}-x_I^{(0)}}$",
    xlabel = L"Time $t$ (h)",
    legend = :bottomright,
    label = [L"\gamma=1, \vartheta=1" L"\gamma=1, \vartheta=10^{-5}" L"\gamma=1, \vartheta=10^{-6}" "Interface position at stationary limit"],
    color = [:blue :orange :red :darkgrey],
    linestyle = [:solid :solid :solid :dash],
    xticks = (ticks, labels));
    display(comparison_plot);
    savefig(comparison_plot, "novel_comp_gamma1_interface.pdf");
    dt=(s_co2["T"][2]-s_co2["T"][1])/s_co2["n_timesteps"];
    ticks = [(i-1) * 3600 for i = 1 : 2:31]
    labels = [(i-1) for i = 1 :2: 31]
    comparison_plot = plot(s_co2["T"][1]:dt:s_co2["T"][2], [ΔT1, ΔT2, ΔT3, ones(s_co2["n_timesteps"]+1)*ΔT1_s],
    ylabel = "\n"*L"Temperature jump $\frac{T_h-T_c}{T_h}$ at $x_I$",
    xlabel = L"Time $t$ (h)",
    legend = :bottomright,
    label = [L"\gamma=1, \vartheta=1" L"\gamma=1, \vartheta=10^{-5}" L"\gamma=1, \vartheta=10^{-6}" "Temperature jump at stationary limit"],
    xticks = (ticks, labels),
    color = [:blue :orange :red :darkgrey],
    linestyle = [:solid :solid :solid :dash]);
    display(comparison_plot);
    savefig(comparison_plot, "novel_comp_gamma1_temp.pdf");
end


if (false)
    ticks = [(i-1) * 3600 for i = 1 : 2:31]
    labels = [(i-1) for i = 1 :2: 31]
    # Temperature difference at interface
    temperature_plot = plot(s["T"][1]:dt:s["T"][2], (interface_1-interface_2)./interface_1,
    ylabel = "\n"*L"Temperature jump $\frac{T_h-T_c}{T_h}$ at $x_I$",
    xlabel = "Time t (s)",
    legend = :bottomright,
    label = ["Temperature jump over time" "Temperature jump at stationary limit"],
    color = [:blue :darkgrey],
    linestyle = [:solid :dash],
    xticks = (ticks, labels));
    display(temperature_plot);
    dt = (s_co2["T"][2] - s_co2["T"][1]) / s_co2["n_timesteps"]
end