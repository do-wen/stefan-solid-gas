import Pkg;
Pkg.add("Plots");
Pkg.add("SpecialFunctions");
Pkg.add("LaTeXStrings");
using SpecialFunctions;
using Plots;
using LaTeXStrings;

default(fontfamily="Computer Modern");

function plot_classical_stefan_problem(p::Dict)
    #Plots.scalefontsizes(1.25); # side-by-side plot
    # Position
    # Determine λ from Mathematica or so
    λ = 0.37813;
    L = 0.06;
    n = 10; # stunden
    t_end = 60*60*n;
    t = collect(range(0, t_end, 500));
    T0 = 298.15;
    ticks = [60*60*(i-1) for i = 1 : n + 1];
    ticklabels = [string(i-1) for i = 1 : n + 1];
    interface_plot = plot(t, vcat(0, 2*√p["a_h"]*λ*sqrt.(t[2:end])),
        xticks = (ticks, ticklabels),
        legend = :none,
        color = :blue);
    ylabel!(L"Interface Position $x_I$ (m)");
    xlabel!(L"Time $t$ (h)");
    distx = 900;
    disty = -0.0009;
    scatter!([60*60*2], [2*√p["a_h"]*λ*sqrt(60*60*2)], markershape = :rect, color = :black, label = "");
    annotate!(60*60*2+distx, 2*√p["a_h"]*λ*sqrt(60*60*2)+disty, L"1", annotationfontsize = 10);
    scatter!([60*60*4], [2*√p["a_h"]*λ*sqrt(60*60*4)], markershape = :rect, color = :black, label = "");
    annotate!(60*60*4+distx, 2*√p["a_h"]*λ*sqrt(60*60*4)+disty, L"2", annotationfontsize = 10);
    scatter!([60*60*6], [2*√p["a_h"]*λ*sqrt(60*60*6)], markershape = :rect, color = :black, label = "");
    annotate!(60*60*6+distx, 2*√p["a_h"]*λ*sqrt(60*60*6)+disty, L"3", annotationfontsize = 10);
    scatter!([60*60*8], [2*√p["a_h"]*λ*sqrt(60*60*8)], markershape = :rect, color = :black, label = "");
    annotate!(60*60*8+distx, 2*√p["a_h"]*λ*sqrt(60*60*8)+disty, L"4", annotationfontsize = 10);
    scatter!([60*60*10], [2*√p["a_h"]*λ*sqrt(60*60*10)], markershape = :rect, color = :black, label = "");
    annotate!(60*60*10+distx, 2*√p["a_h"]*λ*sqrt(60*60*10)+disty, L"5", annotationfontsize = 10);
    display(interface_plot);
    savefig(interface_plot,"interface_classical_stefan.pdf");

    x = collect(range(0, L, 499));
    t = 60*60*2
    x_I = 2*√p["a_h"]*λ*sqrt(t);
    push!(x, x_I);
    sort!(x);
    idx = findall(y->y<x_I+0.000001 &&y > x_I - 0.000001, x)[1];
    tmp_plot = plot([x[1:idx], x[idx:end]], [T0.+(p["T_f"]-T0)*(erf.(x[1:idx]./(2*sqrt.(p["a_h"]*t))))/(erf(2*√p["a_h"]*λ*sqrt.(t)/(2*sqrt.(p["a_h"]*t)))), ones(500-idx+1) * p["T_f"]],
    legend = :topright, color = [:orange :blue], label = ["Temperature in liquid phase" "Temperature in solid phase"]);

    x = collect(range(0, L, 499));
    t = 60*60*4
    x_I = 2*√p["a_h"]*λ*sqrt(t)
    push!(x, x_I);
    sort!(x);
    idx = findall(y->y<x_I+0.000001 &&y > x_I - 0.000001, x)[1];
    plot!([x[1:idx], x[idx:end]], [T0.+(p["T_f"]-T0)*(erf.(x[1:idx]./(2*sqrt.(p["a_h"]*t))))/(erf(2*√p["a_h"]*λ*sqrt.(t)/(2*sqrt.(p["a_h"]*t)))), ones(500-idx+1) * p["T_f"]], color = [:orange :blue], label = ["" ""]);

    x = collect(range(0, L, 499));
    t = 60*60*6
    x_I = 2*√p["a_h"]*λ*sqrt(t)
    push!(x, x_I);
    sort!(x);
    idx = findall(y->y<x_I+0.000001 &&y > x_I - 0.000001, x)[1];
    plot!([x[1:idx], x[idx:end]], [T0.+(p["T_f"]-T0)*(erf.(x[1:idx]./(2*sqrt.(p["a_h"]*t))))/(erf(2*√p["a_h"]*λ*sqrt.(t)/(2*sqrt.(p["a_h"]*t)))), ones(500-idx+1) * p["T_f"]], color = [:orange :blue], label = ["" ""]);

    x = collect(range(0, L, 499));
    t = 60*60*8
    x_I = 2*√p["a_h"]*λ*sqrt(t)
    push!(x, x_I);
    sort!(x);
    idx = findall(y->y<x_I+0.000001 &&y > x_I - 0.000001, x)[1];
    plot!([x[1:idx], x[idx:end]], [T0.+(p["T_f"]-T0)*(erf.(x[1:idx]./(2*sqrt.(p["a_h"]*t))))/(erf(2*√p["a_h"]*λ*sqrt.(t)/(2*sqrt.(p["a_h"]*t)))), ones(500-idx+1) * p["T_f"]], color = [:orange :blue], label = ["" ""]);

    x = collect(range(0, L, 499));
    t = 60*60*10
    x_I = 2*√p["a_h"]*λ*sqrt(t)
    push!(x, x_I);
    sort!(x);
    idx = findall(y->y<x_I+0.000001 &&y > x_I - 0.000001, x)[1];
    plot!([x[1:idx], x[idx:end]], [T0.+(p["T_f"]-T0)*(erf.(x[1:idx]./(2*sqrt.(p["a_h"]*t))))/(erf(2*√p["a_h"]*λ*sqrt.(t)/(2*sqrt.(p["a_h"]*t)))), ones(500-idx+1) * p["T_f"]], color = [:orange :blue], label = ["" ""]);

    distx = 0.0001;
    disty = -0.7;
    annotate!(0.0173+distx, 280+disty, L"1", annotationfontsize = 10);
    annotate!(0.0243+distx, 280+disty, L"2", annotationfontsize = 10);
    annotate!(0.0297+distx, 280+disty, L"3", annotationfontsize = 10);
    annotate!(0.0345+distx, 280+disty, L"4", annotationfontsize = 10);
    annotate!(0.0388+distx, 280+disty, L"5", annotationfontsize = 10);

    ticks = [298.15, 290, 280, 270, 273.15];
    yticks!(ticks);
    ylabel!(L"Temperature $T$ (K)");
    xlabel!(L"Position $x$ (m)");

    display(tmp_plot);
    savefig(tmp_plot,"temperature_classical_stefan.pdf");
    Plots.scalefontsizes();
end



p = Dict(
    # Material parameters, water liquid <--> solid
    "H"   => 333500.0,          # [J / kg]      latent heat of melting
    #"H"   => 2501000.0,         # [J / kg]      latent heat of vaporization
    "T_f" => 273.15,            # [K]           phase-shift temperature
    "rho" => 1000.0,            # [kg / m^3]    density

    "c_h" => 4200.0,            # [J / (kg K)]  specific heat capacity

    "k_h" => 0.6,               # [W / (m K)]   thermal conductivity

);
p["a_h"] = p["k_h"] / (p["rho"] * p["c_h"]); # [m^2 / s]    thermal diffusivity

plot_classical_stefan_problem(p);
