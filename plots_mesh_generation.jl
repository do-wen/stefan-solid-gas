import Pkg;
Pkg.add("Plots");
Pkg.add("SpecialFunctions");
Pkg.add("LaTeXStrings");
using SpecialFunctions;
using Plots;
using LaTeXStrings;

default(fontfamily="Computer Modern");

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

function plot_mesh_generation()
    default(fontfamily="Computer Modern");
    L=10;
    N = 12;
    x_I_1 = 4.4;
    x_I_2 = 4.6;
    x_I_3 = 4.9;
    x_I_4 = 5.2;

    initial_mesh = plot(collect(range(0,L,N-1)), ones(N-1),
        markershape = :circle,
        markersize = 4,
        ylims = (0.5,1.5),
        grid = :none,
        legend = :none,
        yticks = :none,
        showaxis = :x,
        xticks = collect(range(0,L,N-1)),
        size = (600, 100),
        color = :blue
        )
    display(initial_mesh);
    savefig(initial_mesh, "initial_mesh.pdf");
    initial_mesh_w_interface = plot(collect(range(0,L,N-1)), ones(N-1), 
        markershape = :circle,
        markersize = 4,
        ylims = (0.99,1.01),
        grid = :none,
        legend = :none,
        yticks = :none,
        showaxis = :x,
        xticks = collect(range(0,L,N-1)),
        size = (600, 100),
        color = :blue
        )
    plot!([x_I_1], [1], markershape = :circle,
    markersize = 4, color = :orange);
    display(initial_mesh_w_interface);
    savefig(initial_mesh_w_interface, "initial_mesh_w_interface.pdf");
    
    Omega, idx = describe_domain(x_I_1, N, L);

    initial_mesh_done = plot(Omega, ones(N), 
        markershape = :circle,
        markersize = 4,
        ylims = (0.99,1.01),
        grid = :none,
        legend = :none,
        yticks = :none,
        showaxis = :x,
        xticks = collect(range(0,L,N-1)),
        size = (600, 100),
        color = :blue
        )
    plot!([x_I_1], [1], markershape = :circle,
    markersize = 4, color = :orange);
    display(initial_mesh_done);
    savefig(initial_mesh_done, "initial_mesh_done.pdf");
    Omega2, idx = describe_domain(x_I_2, N, L);
    Omega3, idx = describe_domain(x_I_3, N, L);
    Omega4, idx = describe_domain(x_I_4, N, L);

    mesh_prop = plot([Omega, Omega2, Omega3, Omega4], [ones(N), 2*ones(N), 3*ones(N), 4*ones(N)], 
        markershape = :circle,
        markersize = 4,
        ylims = (0.5,4.5),
        ylabel = "\nTime Step",
        grid = :none,
        legend = :none,
        xticks = collect(range(0,L,N-1)),
        yticks = ([1, 2, 3, 4], [L"n-2", L"n-1", L"n", L"n+1"]),
        size = (600, 200),
        color = :blue
        )
    plot!([x_I_1], [1], markershape = :circle,
    markersize = 4, color = :orange);
    plot!([x_I_2], [2], markershape = :circle,
    markersize = 4, color = :orange);
    plot!([x_I_3], [3], markershape = :circle,
    markersize = 4, color = :orange);
    plot!([x_I_4], [4], markershape = :circle,
    markersize = 4, color = :orange);
    display(mesh_prop);
    savefig(mesh_prop, "mesh_prop.pdf");
end


plot_mesh_generation();
