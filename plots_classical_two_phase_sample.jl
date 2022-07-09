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
    if (false)
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
    if (false)
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
    FD_h =   upwind_first(h₂, U[index_I + 2], U[index_I + 3], U[index_I + 4]);
    vᴵ = U[index_I+1];
    ρᶜ = p["rho"];

    # Stefan Condition
    y[index_I + 1] = ρᶜ * p["H"] * vᴵ - p["k_c"] * FD_c + p["k_h"] * FD_h;

    # Trivial interface conditions
    y[index_I]     = U[index_I] - p["T_f"];
    y[index_I + 2] = U[index_I] - U[index_I + 2];

    return y;
end

# Analytical solutions
function xI(t, lambda, a_c)
    return 2*sqrt(a_c)*lambda*sqrt.(t);
end

function Tc(T0, Tsat, a_c, lambda, t, x)
    return T0 .+ (Tsat-T0) .* erf.(x ./ (2*sqrt(t*a_c))) ./ erf(xI(t, lambda, a_c)/(2*sqrt(t*a_c)));
end

function Th(T1, Tsat, a_c, lambda, t, x, a_h)
    return T1 .+ (Tsat-T1) .* erfc.(x ./ (2*sqrt(t*a_h))) ./ erfc(xI(t, lambda, a_c)/(2*sqrt(t*a_h)));
end

function non_uniform_scale_theta_adv_bc_newton_analytical_comparison(p::Dict, s::Dict)
    # Initial condition for analytical comparison
    t_init = s["T"][1];
    # Calculate lambda with some linear algebra program for the paramters in p and s
    lambda = 0.205428;
    xᴵ_init = xI(t_init, lambda, p["a_c"]);

    # Describe domain
    Ω, index_I = describe_domain(xᴵ_init, s["n_points"], s["L"]);
    println("Initial interface at index ", index_I, " from ", s["n_points"]);

    U_c = Tc(s["T_0_init"], p["T_f"], p["a_c"], lambda, t_init, Ω[1:index_I]);
    U_h = Th(s["T_1"], p["T_f"], p["a_c"], lambda, t_init, Ω[index_I:end], p["a_h"]);

    # Solving parameters
    dt = (s["T"][2] - s["T"][1]) / s["n_timesteps"];
    println("Solving from ", s["T"][1], " to ", s["T"][2],
            " in ", s["n_timesteps"], " steps (dt = ", dt, ")");

    x_I = zeros(s["n_timesteps"] + 1);
    ΔT = zeros(s["n_timesteps"] + 1);

    x_I[1] = xᴵ_init;

    θ = 1;
    U = vcat(U_c, 10e-6, U_h);

    # Compressibility parameters
    s["x_I_0"] = x_I[1];
    ΔT[1] = (U_h[1] - U_c[end])/U_h[1];

    # Plotting preparation
    tmp_plot = plot();
    Omega_plot = zeros(s["n_points"], s["n_timesteps"]+1);
    Omega_plot[:,1] = Ω;
    single_T_plot = zeros(s["n_timesteps"] + 1);
    single_T_plot_int = LinearInterpolation(Ω, vcat(U_c, U_h[2:end]));
    single_T_plot[1] = single_T_plot_int(0.2);

    # Time loop
    for k in 1 : s["n_timesteps"]
        t = s["T"][1] + k * dt;
        n_newton_max = 10;
        res = norm(residual_compr(U, p, s, Ω, index_I, U_c, U_h, θ, dt, t), 2);
        n_newton = 0;

        while (res > 1.0e-7 && n_newton < n_newton_max)
            A = ForwardDiff.jacobian(U -> residual_compr(U, p, s, Ω, index_I, U_c, U_h, θ, dt, t), U);
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
        Omega_plot[:,k+1] = Ω;

        index_diff = index_I - index_I_prev;

        # Time stepping: Interpolation
        #interp_linear_c = LinearInterpolation(Ω_prev[1:index_I_prev],   U[1:index_I_prev], extrapolation_bc=Line());
        #interp_linear_h = LinearInterpolation(Ω_prev[index_I_prev:end], U[index_I_prev+2:end], extrapolation_bc=Line());
        #U_c = interp_linear_c(Ω[1:index_I]);
        #U_h = interp_linear_h(Ω[index_I:end]);

        # Time stepping: Projection
        U_c = l2_projection(Ω_prev[1:index_I_prev], Ω[1:index_I], U[1:index_I_prev]);
        U_h = l2_projection(Ω_prev[index_I_prev:end], Ω[index_I:end], U[index_I_prev+2:end]);

        # Time stepping: Scaling transformation
        if (false)
            if (index_diff >= 0)
                U_c = zeros(index_I);
                U_c[1:index_I_prev-2] = U[1:index_I_prev-2];
                Ω_sc = [Ω_prev[index_I_prev - 2], (Ω[index_I] + Ω_prev[index_I_prev - 2]) / 2, Ω[index_I]];
                U_c[index_I_prev-1:index_I] = l2_projection(Ω_sc, Ω[index_I_prev-1:index_I], U[index_I_prev-2:index_I_prev]);
                U_c[index_I] = U[index_I_prev];

                U_h = zeros(s["n_points"] - index_I + 1);
                Ω_sc = zeros(3 + index_diff);
                Ω_sc[1] = Ω[index_I];
                h_2_tilde = (Ω[index_I+2] - Ω[index_I]) / (2 + index_diff * (Ω_prev[index_I_prev+3] - Ω_prev[index_I_prev+2]) / (Ω_prev[index_I_prev+1] - Ω_prev[index_I_prev]));
                h_2 = h_2_tilde * (Ω_prev[index_I_prev+3] - Ω_prev[index_I_prev+2]) / (Ω_prev[index_I_prev+1] - Ω_prev[index_I_prev]);
                Ω_sc[2] = Ω_sc[1] + h_2_tilde;
                Ω_sc[3] = Ω_sc[2] + h_2_tilde;
                for k in 4 : 2 + index_diff
                    Ω_sc[k] = Ω_sc[k-1] + h_2;
                end
                Ω_sc[end] = Ω[index_I+2];
                U_h[1:3] = l2_projection(Ω_sc, Ω[index_I:index_I+2], U[index_I_prev+2:index_I+4]);
                U_h[1] = U[index_I_prev+2];
                U_h[3:end] = U[index_I+4:end];
            else
                U_h = zeros(s["n_points"] - index_I + 1);
                U_h[3-index_diff:end] = U[index_I_prev+4:end];
                Ω_sc = [Ω[index_I], (Ω[index_I] + Ω_prev[index_I_prev + 2]) / 2, Ω_prev[index_I_prev+2]];
                U_h[1:2-index_diff] = l2_projection(Ω_sc, Ω[index_I:index_I_prev+1], U[index_I_prev+2:index_I_prev+4]);
                U_h[1] = U[index_I_prev+2];

                U_c = zeros(index_I);
                U_c[1:index_I-2] = U[1:index_I-2];
                Ω_sc = zeros(3-index_diff);
                Ω_sc[end] = Ω[index_I];
                h_2_tilde = (Ω[index_I] - Ω[index_I-2]) / (2 - index_diff * (Ω_prev[index_I_prev-2] - Ω_prev[index_I_prev-3]) / (Ω_prev[index_I_prev] - Ω_prev[index_I_prev-1]));
                h_2 = h_2_tilde * (Ω_prev[index_I_prev-2] - Ω_prev[index_I_prev-3]) / (Ω_prev[index_I_prev] - Ω_prev[index_I_prev-1]);
                Ω_sc[end-1] = Ω_sc[end] - h_2_tilde;
                Ω_sc[end-2] = Ω_sc[end-1] - h_2_tilde;
                for k in 4 : 2 - index_diff
                    Ω_sc[end-k+1] = Ω_sc[end-k+2] - h_2;
                end
                Ω_sc[1] = Ω[index_I-2];
                U_c[end-2:end] = l2_projection(Ω_sc, Ω[index_I-2:index_I], U[index_I-2:index_I_prev]);
                U_c[end] = U[index_I_prev];
            end
        end

        # Plotting
        if (t == 0.002 && (false))
            # When set to true, add 25 to all temperatures
            # (saturation, boundary), otherwise the error calculation is not
            # correct
            L2_err = abs(x_I[k+1] - xI(t, lambda, p["a_c"])) / abs(xI(t, lambda, p["a_c"]));
            println("L2 error of Interface Position: ", L2_err, " (this was with n_timesteps = ", s["n_timesteps"], "and with n_points = ", s["n_points"], ")!", "at time ", t);
            T_a = vcat(Tc(s["T_0_init"], p["T_f"], p["a_c"], lambda, t, Ω_prev[1:index_I_prev]),Th(s["T_1"], p["T_f"], p["a_c"], lambda, t, Ω_prev[index_I_prev:end], p["a_h"]));
            T_n = vcat(U[1:index_I_prev], U[index_I_prev+2:end]);
            L2_err = sqrt(1/s["n_points"] * transpose((T_n - T_a) ./ T_a) * ((T_n - T_a) ./ T_a));
            println("L2 error of Temperature: ", L2_err, " (this was with n_timesteps = ", s["n_timesteps"], "and with n_points = ", s["n_points"], ")! at time ", t);
            return;
        end

        if ((mod(k, 41) == 0 && mod(k,2)!=0) || k == 2)
            # Temperature Plot
            x_for_analytical1 = collect(range(0,xI(t, lambda, p["a_c"]), 400));
            x_for_analytical2 = collect(range(xI(t, lambda, p["a_c"]), s["L"], 400));

            if (k == 2)
                plot!(
                [Ω_prev[1:index_I_prev],
                Ω_prev[index_I_prev:end],
                x_for_analytical1,
                x_for_analytical2],
                [U[1:index_I_prev],
                U[index_I_prev+2:end],
                Tc(s["T_0_init"], p["T_f"], p["a_c"], lambda, t, x_for_analytical1),
                Th(s["T_1"], p["T_f"], p["a_c"], lambda, t, x_for_analytical2, p["a_h"])],
                color = [:blue :blue :orange :orange],
                linestyle = [:solid :solid :dash :dash],
                label = [L"T_{\textrm{numerical}}" "" L"T_{\textrm{analytical}}" ""],
                legend = :bottomright);
                scatter!([xI(t, lambda, p["a_c"])], [0], markershape = :rect, color = :black, label = "");
                annotate!(xI(t, lambda, p["a_c"])+0.018, 0, L"1", annotationfontsize = 10);
                ylabel!(L"T");
                xlabel!(L"x");
                xlims!((0,1));
            else
                plot!(
                    [Ω_prev[1:index_I_prev],
                    Ω_prev[index_I_prev:end],
                    x_for_analytical1,
                    x_for_analytical2],
                    [U[1:index_I_prev],
                    U[index_I_prev+2:end],
                    Tc(s["T_0_init"], p["T_f"], p["a_c"], lambda, t, x_for_analytical1),
                    Th(s["T_1"], p["T_f"], p["a_c"], lambda, t, x_for_analytical2, p["a_h"])],
                    color = [:blue :blue :orange :orange],
                    linestyle = [:solid :solid :dash :dash],
                    label = "",
                    legend = :bottomright);
                    scatter!([xI(0.042, lambda, p["a_c"])], [0], markershape = :rect, color = :black, label = "");
                    annotate!(xI(0.042, lambda, p["a_c"])+0.018, 0, L"2", annotationfontsize = 10);
                    scatter!([xI(0.124, lambda, p["a_c"])], [0], markershape = :rect, color = :black, label = "");
                    annotate!(xI(0.124, lambda, p["a_c"])+0.018, 0, L"3", annotationfontsize = 10);
                    scatter!([xI(0.206, lambda, p["a_c"])], [0], markershape = :rect, color = :black, label = "");
                    annotate!(xI(0.206, lambda, p["a_c"])+0.018, 0, L"4", annotationfontsize = 10);
                    scatter!([xI(0.288, lambda, p["a_c"])], [0], markershape = :rect, color = :black, label = "");
                    annotate!(xI(0.288, lambda, p["a_c"])+0.018, 0, L"5", annotationfontsize = 10);
            end
        end

        single_T_plot_int = LinearInterpolation(Ω_prev, vcat(U[1:index_I_prev], U[index_I_prev+3:end]));
        single_T_plot[k+1] = single_T_plot_int(0.2);

        U = vcat(U_c, U[index_I_prev+1], U_h);
    end
    display(tmp_plot);
    savefig(tmp_plot,"temperature_classical_stefan.pdf");
    # Interface Plot
    interface_plot = plot(s["T"][1]:dt:s["T"][2], [x_I, xI(s["T"][1]:dt:s["T"][2], lambda, p["a_c"])],
                            ylabel = L"x_I",
                            xlabel = L"t",
                            legend = :bottomright,
                            label = [L"x_{I,\textrm{numerical}}" L"x_{I,\textrm{analytical}}"],
                            linestyle = [:solid :dash],
                            color = [:blue :orange]);
    scatter!([s["T"][1]+dt], [xI(s["T"][1]+dt, lambda, p["a_c"])], markershape = :rect, color = :black, label = "");
    annotate!(s["T"][1]+dt+0.008, xI(s["T"][1]+dt, lambda, p["a_c"]), L"1", annotationfontsize = 10);
    scatter!([0.042], [xI(0.042, lambda, p["a_c"])], markershape = :rect, color = :black, label = "");
    annotate!(0.042+0.008, xI(0.042, lambda, p["a_c"]), L"2", annotationfontsize = 10);
    scatter!([0.124], [xI(0.124, lambda, p["a_c"])], markershape = :rect, color = :black, label = "");
    annotate!(0.124+0.008, xI(0.124, lambda, p["a_c"]), L"3", annotationfontsize = 10);
    scatter!([0.206], [xI(0.206, lambda, p["a_c"])], markershape = :rect, color = :black, label = "");
    annotate!(0.206+0.008, xI(0.206, lambda, p["a_c"]), L"4", annotationfontsize = 10);
    scatter!([0.288], [xI(0.288, lambda, p["a_c"])], markershape = :rect, color = :black, label = "");
    annotate!(0.288+0.008, xI(0.288, lambda, p["a_c"]), L"5", annotationfontsize = 10);
    savefig(interface_plot,"interface_plot.pdf");
    display(interface_plot);

    # Domain Plot
    domain_plot = plot(Omega_plot', s["T"][1]:dt:s["T"][2],
        legend = :none,
        color_palette = palette([:blue, :orange, :darkred]),
        xlabel = L"x",
        ylabel = L"t",
        xlim = (-0.01,0.3));
    display(domain_plot);
    savefig(domain_plot,"domain_plot.pdf");

    # Single Temp Plot
    single_T_anal = zeros(s["n_timesteps"]+1);
    for i in 1 : s["n_timesteps"]+1
        t = s["T"][1]+(i-1)*dt;
        x_I = xI(t, lambda, p["a_c"]);
        T_c_anal = Tc(s["T_0_init"], p["T_f"], p["a_c"], lambda, t, 0.2);
        T_h_anal = Th(s["T_1"], p["T_f"], p["a_c"], lambda, t, 0.2, p["a_h"]);
        if (x_I < 0.2)
            single_T_anal[i] = T_h_anal;
        else
            single_T_anal[i] = T_c_anal;
        end
    end
    single_T = plot(s["T"][1]:dt:s["T"][2], [single_T_plot single_T_anal],
    label = [L"T_{\textrm{numerical}}(0.2, t)" L"T_{\textrm{analytical}}(0.2, t)"],
                            linestyle = [:solid :dash],
                            color = [:blue :orange],
                            xlabel = L"t",
                            ylabel = L"T");
    display(single_T);
    savefig(single_T,"single_T.pdf");
end

p = Dict(
    # Material parameters, water liquid <--> solid
    "H"   => 338,    # [J / kg]     latent heat of sublimation
    "T_f" => 0,       # [K]          phase-shift temperature (at 1000 Pa)
    "rho" => 1,      # [kg / m^3]   density #1600

    "c_h" => 4.226,      # [J / (kg K)] specific heat capacity
    "c_c" => 1.762,     # [J / (kg K)] specific heat capacity
    "k_h" => 0.556,   # [W / (m K)]  thermal conductivity # 1
    "k_c" => 2.22, # [W / (m K)]  thermal conductivity # 0.025
);
p["a_h"] = p["k_h"] / (p["rho"] * p["c_h"]); # [m^2 / s]    thermal diffusivity
p["a_c"] = p["k_c"] / (p["rho"] * p["c_c"]); # [m^2 / s]    thermal diffusivity

s = Dict(
    "T"           => [0.001, 0.288],  # [[s], [s]] Solving time
    "n_timesteps" => 287,        # [-]        number of timesteps
    "n_points"    => 100,        # [-]        number of spacial discretization points
    "L"           => 1,        # [m]        canal length
    "T_1"         => 10,        # [K]        fixed temperature at hot surface (BC, IC)
    "T_0_init"    => -20,        # [K]        initial temperature at cold surface (IC)
    "C"           => 100,         # [1 / m]    constant of Newton's cooling (BC)
    "T_ext"       => -20,  #20    # [K]        ambient temperature for Newton's cooling (BC)
    "C_h"         => 100,         # [1 / m]    constant of Newton's cooling on hot side
    "T_ext_h"     => 10,  #220  # [K]        ambient temperature hot side
);

interface_1 = non_uniform_scale_theta_adv_bc_newton_analytical_comparison(p, s);

if (true)
    # Plot for error sample problem. To reproduce, set corresponding 'if' in
    # non_uniform_scale_theta_adv_bc_newton_analytical_comparison to 'true'
    # (line 275)
    cond_N = [100, 200, 500, 1000, 2000, 5000];
    #cond_k = [0.01394777792796869, 0.0027911733872067837, 0.0004254370334512404, 0.00010417764666112757, 2.5697755234483897e-5, 3.939775721202375e-6];
    cond_k = [0.04435602535767886, 0.005661603667922837, 0.0008080062274602715, 0.00019193008764232646, 5.038833348557257e-5, 7.88857602904146e-6]
    cond_t = [0.10274624564134395, 0.01829149516233783, 0.0027915593956100125, 0.0006903761313142608, 0.00016987681940149537, 2.4868803096266185e-5];
    slope_x = [100, 5000];
    slope_3 = [100^(-3)*0.01^(-3), 5000^(-3)*0.01^(-3)];
    slope_2 = [100^(-2)*350, 5000^(-2)*350];
    slope_1 = [100^(-2), 10000^(-2)];
    ticklabels = [L"100", L"200", L"500", L"1000", L"2000", L"5000"];

    yticks = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
    yticklabels = [L"10^{-1}", L"10^{-2}", L"10^{-3}", L"10^{-4}", L"10^{-5}", L"10^{-6}"];
    default(fontfamily="Computer Modern");
    condition_number = plot([cond_N, cond_N, slope_x], [cond_t, cond_k, slope_2],
        xaxis = :log,
        yaxis = :log,
        xticks = (cond_N, ticklabels),
        yticks = (yticks, yticklabels),
        ylims = (yticks[end]/2, 2*yticks[1]),
        label = ["error interface position" "error temperature" L"slope $2$"],
        xlabel = L"n_x",
        ylabel = "error",
        linestyle = [:solid :solid :dash],
        markershape = [:circle :circle :none],
        markersize = 4,
        color = [:blue :orange :darkgrey]
    )
    display(condition_number);
    savefig(condition_number,"convergence_two_phase.pdf");
end


if (true)
    # Plot for error sample problem with different time stepping. To reproduce,
    # set corresponding 'if' in non_uniform_scale_theta_adv_bc_newton_analytical_comparison
    # to 'true' (line 275) and use desired time stepping method (lines 217-273)
    cond_N = [100, 200, 500, 1000, 2000];
    #cond_k = [0.01394777792796869, 0.0027911733872067837, 0.0004254370334512404, 0.00010417764666112757, 2.5697755234483897e-5, 3.939775721202375e-6];
    proj_T = [0.04435602535767886, 0.005661603667922837, 0.0008080062274602715, 0.00019193008764232646, 5.038833348557257e-5]
    proj_t = [0.10274624564134395, 0.01829149516233783, 0.0027915593956100125, 0.0006903761313142608, 0.00016987681940149537];
    scale_T = [0.04435602535767886, 0.004642582344313742, 0.0004220822208193931, 1.803321879852518e-5, 4.409884700982907e-5];
    scale_t = [0.10274624564134395, 0.026232540359088837, 0.006216921763188181, 0.0024939205145511124, 0.001089293271963888];
    int_T = [0.04435602535767886, 0.005575941624863174, 0.0008002421958919414, 0.0001899599503041007, 4.6872473493938944e-5];
    int_t = [0.10274624564134395, 0.018988549420947062, 0.002860518977458971, 0.000709584894749687, 0.00017492228509802398];
    slope_x = [100, 2000];
    slope_3 = [100^(-3), 2000^(-3)];
    slope_2 = 400*[100^(-2), 2000^(-2)];
    slope_1 = 4*[100^(-1), 2000^(-1)];
    ticklabels = [L"100", L"200", L"500", L"1000", L"2000"];

    yticks = [0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5];
    yticklabels = ["", L"10^{-1}", L"10^{-2}", L"10^{-3}", L"10^{-4}", L"10^{-5}"];
    default(fontfamily="Computer Modern");
    condition_number = plot([cond_N, cond_N, cond_N, cond_N, cond_N, cond_N, slope_x, slope_x], [int_t, int_T, proj_t, proj_T, scale_t, scale_T, slope_1, slope_2],
        xaxis = :log,
        yaxis = :log,
        xticks = (cond_N, ticklabels),
        yticks = (yticks, yticklabels),
        ylims = (yticks[end]/2, 2*yticks[1]),
        label = ["interpolation: error interface" "interpolation: error temperature" "projection: error interface" "projection: error temperature" "scaling: error interface" "scaling: error temperature" L"slope $1$" L"slope $2$"],
        xlabel = L"n_x",
        ylabel = "error",
        linestyle = [:solid :dash :solid :dash :solid :dash :dot :solid],
        markershape = [:circle :circle :circle :circle :circle :circle :none :none],
        markersize = 4,
        color = [:blue :blue :orange :orange :red :red :darkgrey :darkgrey]
    )
    display(condition_number);
    savefig(condition_number,"convergence_two_phase_test_bla.pdf");
end
