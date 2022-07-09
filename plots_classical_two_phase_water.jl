#import Pkg;
Pkg.add("Plots");
Pkg.add("Printf");
Pkg.add("LinearAlgebra");
Pkg.add("SpecialFunctions");
Pkg.add("Interpolations");
Pkg.add("ForwardDiff");
Pkg.add("LaTeXStrings");
Pkg.add("LinearSolve");
Pkg.add("Trapz");
Pkg.add("ColorSchemes");
using Trapz;
using Plots;
using Printf;
using LinearAlgebra;
using SpecialFunctions;
using Interpolations;
using ForwardDiff;
using Statistics;
using LinearSolve;
using LaTeXStrings;
using ColorSchemes;

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

    aʰ = p["k_h"] / (p["rho"] * p["c_h"]);

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


    T_ext = s["T_ext_h"];

    if (true)
        hₛ = Ω[2] - Ω[1];
        F_end = 1 + dt * aʰ * 2 * (1 + s["C_h"] * hₛ) / hₛ^2;
        F_end_p = 1;
        F_end_1 = 2 * dt * aʰ / hₛ^2;
        y[end] = U[end] * F_end - U[end-1] * F_end_1 - Uʰ[end] * F_end_p - 2 * dt * aʰ * s["C_h"] * T_ext / hₛ;
    else
        y[s["n_points"] + 2] = U[s["n_points"] + 2] - T_ext;
    end

    # Helpful definitions
    h₁ = Ω[index_I]     - Ω[index_I - 1];
    h₂ = Ω[index_I + 1] - Ω[index_I];
    FD_c = downwind_first(h₁, U[index_I - 2], U[index_I - 1], U[index_I    ]);
    FD_h = upwind_first(h₂, U[index_I + 2], U[index_I + 3], U[index_I + 4]);
    vᴵ = U[index_I+1];
    ρᶜ = p["rho"];

    # Stefan Condition
    y[index_I + 1] = ρᶜ * p["H"] * vᴵ - p["k_c"] * FD_c + p["k_h"] * FD_h;

    # Trivial interface conditions
    y[index_I]     = U[index_I] - p["T_f"];
    y[index_I + 2] = U[index_I] - U[index_I + 2];
    return y;
end

function non_uniform_scale_theta_adv_bc_newton_analytical_comparison(p::Dict, s::Dict)

    # First, trivial IC
    xᴵ_init = s["L"] / (1 + p["k_h"] / p["k_c"] * (s["T_1"] - p["T_f"]) / (p["T_f"] - s["T_0_init"]));
    ∂ₓTᶜ_init = (p["T_f"] - s["T_0_init"]) / xᴵ_init;
    ∂ₓTʰ_init = (s["T_1"] - p["T_f"]) / (s["L"] - xᴵ_init);
      Tᶜ_init = p["T_f"];
      Tʰ_init = p["T_f"];
      T₁_init = s["T_1"];
      T₀_init = s["T_0_init"];

    # Describe domain
    Ω, index_I = describe_domain(xᴵ_init, s["n_points"], s["L"]);
    println("initial interface at index ", index_I, " from ", s["n_points"]);

    U_c = ∂ₓTᶜ_init .* Ω[1:index_I]   .+ T₀_init;
    U_h = ∂ₓTʰ_init .* Ω[index_I:end] .+ (Tʰ_init*s["L"] - T₁_init*xᴵ_init) / (s["L"] - xᴵ_init);

    # Solving parameters
    dt = (s["T"][2] - s["T"][1]) / s["n_timesteps"];
    println("Solving from ", s["T"][1], " to ", s["T"][2],
            " in ", s["n_timesteps"], " steps (dt = ", dt, ")");

    x_I = zeros(s["n_timesteps"] + 1);
    ΔT = zeros(s["n_timesteps"] + 1);

    x_I[1] = xᴵ_init;

    xᴵ_stat = (s["C"]*p["k_c"]*(1+s["C_h"]*s["L"])*(s["T_ext"]-p["T_f"])+s["C_h"]*p["k_h"]*(s["T_ext_h"]-p["T_f"])) / (s["C_h"]*s["C"]*(p["k_c"]*(s["T_ext"]-p["T_f"])-p["k_h"]*(s["T_ext_h"]-p["T_f"])));
    println(xᴵ_stat);
    ∂ₓTᶜ_stat = p["k_h"]/p["k_c"]*(s["C_h"]*(s["T_ext_h"]-p["T_f"]))/(1+s["C_h"]*(s["L"]-xᴵ_stat));
    ∂ₓTʰ_stat = (s["C_h"]*(s["T_ext_h"]-p["T_f"]))/(1+s["C_h"]*(s["L"]-xᴵ_stat));
      Tᶜ_stat = p["T_f"];
      Tʰ_stat = p["T_f"];
      T₁_stat = s["T_ext_h"] + (p["T_f"]-s["T_ext_h"])/(1+s["C_h"]*(s["L"]-xᴵ_stat));
      T₀_stat = s["T_ext"] + p["k_h"]/p["k_c"]*s["C_h"]/s["C"]*(s["T_ext_h"]-p["T_f"])/(1+s["C_h"]*(s["L"]-xᴵ_stat));

    # Plotting preparation
    plot_stat_c = collect(range(0,xᴵ_stat,300));
    plot_stat_h = collect(range(xᴵ_stat,s["L"],300));
    U_c_stat = ∂ₓTᶜ_stat .* plot_stat_c   .+ T₀_stat;
    U_h_stat = ∂ₓTʰ_stat .* plot_stat_h .+ (Tʰ_stat*s["L"] - T₁_stat*xᴵ_stat) / (s["L"] - xᴵ_stat);
    U_c_init = U_c;
    U_h_init = U_h;
    plot_init_c = Ω[1:index_I];
    plot_init_h = Ω[index_I:end];

    θ = 1;
    U = vcat(U_c, 1e-6, U_h);
    # Compressibility parameters
    s["x_I_0"] = x_I[1];

    ΔT[1] = (U_h[1] - U_c[end])/U_h[1];

    temp_plot = plot();
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
        Ω_prev = Ω;
        index_I_prev = index_I;

        x_I[k+1] = x_I[k] + dt * U[index_I + 1];
        ΔT[k+1]  = (U[index_I+2]-U[index_I])/U[index_I+2];

        Ω, index_I = describe_domain(x_I[k+1], s["n_points"], s["L"]);
        index_diff = index_I - index_I_prev;


        U_c = l2_projection(Ω_prev[1:index_I_prev], Ω[1:index_I], U[1:index_I_prev]);
        U_h = l2_projection(Ω_prev[index_I_prev:end], Ω[index_I:end], U[index_I_prev+2:end]);


        # Temperature Plot
        if (k == 10 || k == 50 || (mod(k, 100) == 0 && k != 700 && k != 600 && k != 500))
            if (k == 100)
                plot!(
                [Ω_prev[1:index_I_prev],
                Ω_prev[index_I_prev:end]],
                [U[1:index_I_prev],
                U[index_I_prev+2:end]],
                color = [:blue :orange],
                label = [L"T_{c}" L"T_{h}"],
                legend = :bottomright);
                scatter!([0.009], [273.15], markershape = :rect, color = :black, label = "");
                annotate!(0.009+0.0018, 273.15, L"1", annotationfontsize = 10);
                ylabel!(L"Temperature $T$ (K)");
                xlabel!(L"Position $x$ (m)");
            else
                plot!(
                    [Ω_prev[1:index_I_prev],
                    Ω_prev[index_I_prev:end]],
                    [U[1:index_I_prev],
                    U[index_I_prev+2:end]],
                    color = [:blue :orange],
                    label = "",
                    legend = :bottomright);
            end
            if (k == 800)
                plot!([plot_init_c, plot_init_h, plot_stat_c, plot_stat_h], [U_c_init, U_h_init, U_c_stat, U_h_stat],
                    color = :darkgrey,
                    linestyle = [:solid :solid :dash :dash],
                    label = ["initial condition" "" "stationary limit" ""])
                scatter!([0.026], [273.15], markershape = :rect, color = :black, label = "");
                annotate!(0.026+0.0018, 273.15, L"2", annotationfontsize = 10);
                scatter!([0.039], [273.15], markershape = :rect, color = :black, label = "");
                annotate!(0.039+0.0018, 273.15, L"3", annotationfontsize = 10);
                scatter!([0.0555], [273.15], markershape = :rect, color = :black, label = "");
                annotate!(0.0555+0.0018, 273.15, L"4", annotationfontsize = 10);
                scatter!([0.0655], [273.15], markershape = :rect, color = :black, label = "");
                annotate!(0.0655+0.0018, 273.15, L"5", annotationfontsize = 10);
                scatter!([0.0705], [273.15], markershape = :rect, color = :black, label = "");
                annotate!(0.0705+0.0018, 273.15, L"6", annotationfontsize = 10);
                scatter!([0.0745], [273.15], markershape = :rect, color = :black, label = "");
                annotate!(0.0745+0.0018, 273.15, L"7", annotationfontsize = 10);
            end
        end

        U = vcat(U_c, U[index_I_prev+1], U_h);
    end
    display(temp_plot); savefig(temp_plot,"temperature_plot_water.pdf");

    # Interface Plot
    ticks = [(i-1) * 3600 for i = 1 : 21]
    labels = [(i-1) for i = 1 : 21]
    interface_plot = plot(s["T"][1]:dt:s["T"][2], x_I,
                            ylabel = L"Interface Position $x_I$ (m)",
                            xlabel = L"Time $t$ (h)",
                            legend = :bottomright,
                            label = :none,
                            color = :blue,
                            xticks = (ticks, labels));
    dist = 2*900;
    iI = LinearInterpolation(collect(s["T"][1]:dt:s["T"][2]), x_I);
    scatter!([900], [iI(900)], markershape = :rect, color = :black, label = "");
    annotate!(900+dist, iI(900), L"1", annotationfontsize = 10);
    scatter!([4500], [iI(4500)], markershape = :rect, color = :black, label = "");
    annotate!(4500+dist, iI(4500), L"2", annotationfontsize = 10);
    scatter!([9000], [iI(9000)], markershape = :rect, color = :black, label = "");
    annotate!(9000+dist, iI(9000), L"3", annotationfontsize = 10);
    scatter!([18000], [iI(18000)], markershape = :rect, color = :black, label = "");
    annotate!(18000+dist, iI(18000), L"4", annotationfontsize = 10);   
    scatter!([27000], [iI(27000)], markershape = :rect, color = :black, label = "");
    annotate!(27000+dist, iI(27000), L"5", annotationfontsize = 10);
    scatter!([36000], [iI(36000)], markershape = :rect, color = :black, label = "");
    annotate!(36000+dist, iI(36000), L"6", annotationfontsize = 10);
    scatter!([72000], [iI(72000)], markershape = :rect, color = :black, label = "");
    annotate!(72000+dist, iI(72000), L"7", annotationfontsize = 10);
    savefig(interface_plot,"interface_plot_water.pdf");
    display(interface_plot);
    println("Interface position at t_end: ", x_I[end]);
end


p = Dict(
    # Material parameters, water liquid <--> solid
    "H"   => 333500,    # [J / kg]     latent heat of sublimation
    "T_f" => 273.15,       # [K]          phase-shift temperature (at 1000 Pa)
    "rho" => 1000,      # [kg / m^3]   density #1600

    "c_h" => 4200,      # [J / (kg K)] specific heat capacity
    "c_c" => 2060,     # [J / (kg K)] specific heat capacity
    "k_h" => 0.6,     # [W / (m K)]  thermal conductivity
    "k_c" => 2.33,    # [W / (m K)]  thermal conductivity
);
p["a_h"] = p["k_h"] / (p["rho"] * p["c_h"]); # [m^2 / s]    thermal diffusivity
p["a_c"] = p["k_c"] / (p["rho"] * p["c_c"]); # [m^2 / s]    thermal diffusivity

s = Dict(
    "T"           => [0, 20*3600],  # [[s], [s]] Solving time
    "n_timesteps" => 800,        # [-]        number of timesteps
    "n_points"    => 100,        # [-]        number of spacial discretization points
    "L"           => 0.08,        # [m]        canal length
    "T_1"         => 290,        # [K]        fixed temperature at hot surface (BC, IC)
    "T_0_init"    => 273.1,        # [K]        initial temperature at cold surface (IC)
    "C"           => 100,         # [1 / m]    constant of Newton's cooling (BC)
    "T_ext"       => 250,     # [K]        ambient temperature for Newton's cooling (BC)
    "C_h"         => 100,         # [1 / m]    constant of Newton's cooling on hot side
    "T_ext_h"     => 290,  # [K]        ambient temperature hot side
);

interface_1 = non_uniform_scale_theta_adv_bc_newton_analytical_comparison(p, s);
