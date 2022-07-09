import Pkg;
Pkg.add("Plots");
Pkg.add("Printf");
Pkg.add("LinearAlgebra");
Pkg.add("SpecialFunctions");
Pkg.add("Interpolations");
Pkg.add("LaTeXStrings");
using Plots;
using SpecialFunctions;
using Interpolations;
using LinearAlgebra;
using Printf;
using LaTeXStrings;
default(fontfamily="Computer Modern");

include("projection.jl");

function Plotting_MA_INT()
    #INTERPOLATION
    Omega1   = [0, 1, 2, 3, 3.8, 4.6, 5.3, 6, 7, 8, 9, 10, 11];
    Omega2 = [0, 1, 2, 3, 4, 5, 5.6, 6.2, 7.1, 8, 9, 10, 11];
    U1 = collect(range(1,3,14));
    U1 .*= -U1.*U1;
    U1 .+= 283;
    reverse!(U1);
    U1_c = U1[1:6];
    U1_h = U1[7:14];

    interp_linear_c = LinearInterpolation(Omega1[4:6], U1_c[4:6], extrapolation_bc=Flat());
    interp_linear_h = LinearInterpolation(Omega1[6:10], U1_h[1:5], extrapolation_bc=Flat());
    U2_c_int = vcat(U1_c[1:3], interp_linear_c(Omega2[4:8]));
    U2_h_int = vcat(interp_linear_h(Omega2[8:10]), U1_h[6:8]);

    U2_c_proj = vcat(U1_c[1:3], l2_projection(Omega1[4:6], Omega2[4:8], U1_c[4:6]));
    U2_h_proj = vcat(l2_projection(Omega1[6:10], Omega2[8:10], U1_h[1:5]), U1_h[6:8]);

    # Scale
    x_I_1 = 6;
    x_I_2 = 8;

    x_I_diff = x_I_2 - x_I_1;
    U2_c_sc = zeros(x_I_2);
    U2_c_sc[1:x_I_1-2] = U1_c[1:x_I_1 - 2];
    Omega1_sc = [Omega1[x_I_1-2], (Omega2[x_I_2]+Omega1[x_I_1-2])/2, Omega2[x_I_2]];
    U2_c_sc[x_I_1-1:x_I_2] = l2_projection(Omega1_sc, Omega2[x_I_1-1:x_I_2], U1_c[x_I_1-2:x_I_1]);
    U2_c_sc[x_I_2] = U1_c[x_I_1];
    println(Omega1_sc);

    U2_h_sc = zeros(13-x_I_2+1);
    Omega1_sc = zeros(3+x_I_2-x_I_1);
    Omega1_sc[1] = Omega2[x_I_2];
    h_2_tilde = (Omega2[x_I_2+2]-Omega2[x_I_2]) / (2 + x_I_diff * (Omega1[x_I_1+3]-Omega1[x_I_1+2])/(Omega1[x_I_1+1]-Omega1[x_I_1]));
    h_2 = h_2_tilde * (Omega1[x_I_1+3]-Omega1[x_I_1+2])/(Omega1[x_I_1+1]-Omega1[x_I_1]);
    Omega1_sc[2] = Omega1_sc[1]+h_2_tilde;
    Omega1_sc[3] = Omega1_sc[2]+h_2_tilde;
    for k in 4:3+x_I_2-x_I_1-1
        Omega1_sc[k] = Omega1_sc[k-1] + h_2;
    end
    Omega1_sc[end] = Omega2[x_I_2 + 2];
    U2_h_sc[1:3] = l2_projection(Omega1_sc, Omega2[x_I_2:x_I_2+2], U1_h[1:3+x_I_2-x_I_1]);
    U2_h_sc[1] = U1_h[1];
    U2_h_sc[3:end] = U1_h[3+x_I_2-x_I_1:end];

    # Select either U2_c_sc/U2_h_sc OR U2_c_proj/U2_h_proj OR U2_c_int/U2_h_int for plotting
    U_c_plt = U2_c_sc;
    U_h_plt = U2_h_sc;

    ticks = vcat(Omega1, Omega2);
    Omega2_Str = String[];
    for elem in Omega2
        push!(Omega2_Str, "\n"*string(elem));
    end
    ticklabels = vcat(Omega1, Omega2_Str);
    plot_plt = plot();
    plot!([4.0, 4.0], [265, 280], color=:lightblue, label="", linestyle=:dash);
    plot!([5.0, 5.0], [265, 280], color=:lightblue, label="", linestyle=:dash);
    plot!([5.6, 5.6], [265, 280], color=:lightblue, label="", linestyle=:dash);
    plot!([6.2, 6.2], [265, 280], color=:green, label="", linestyle=:dash);
    plot!([7.1, 7.1], [265, 280], color=:orange, label="", linestyle=:dash);
    plot!([Omega1[1:6], Omega1[6:13], Omega2[1:8], Omega2[8:13]], [U1_c, U1_h, U_c_plt, U_h_plt],
        markeralpha=[0.5 0.5 1 1],
        linealpha=[0.5 0.5 1 1],
        markershape=[:circle :circle :diamond :diamond],
        label = [L"$T_c^{(n)}$" L"$T_h^{(n)}$" L"$T_{c, scale}^{(n\rightarrow n+1)}$" L"$T_{h, scale}^{(n\rightarrow n+1)}$"],
        linestyle = [:dash :dash :dash :dash],
        markercolor = [:blue :red],
        linecolor = [:blue :red],
        xlabel = L"Position $x$",
        ylabel = L"Temperature $T$",
        legend = :topleft,
        xticks = (ticks, ticklabels),
        annotations = [([4.35],[262.2], text(L"Interface at time step $n$", rotation=90.0, pointsize=10, color=:blue, "Computer Modern")), ([5.95],[263.5], text(L"Interface at time step $n+1$", rotation=90.0, pointsize=10, color=:green, "Computer Modern"))]);
    vline!([4.6], color=:blue, label="", linestyle=:dot);
    vline!([6.2], color=:green, label="", linestyle=:dot);

    display(plot_plt);
    savefig(plot_plt,"time_stepping.pdf");

end

# Select desired time stepping method in the function itself
Plotting_MA_INT();