import Pkg;
Pkg.add("LinearAlgebra");
Pkg.add("SpecialFunctions");
Pkg.add("Interpolations");
using SpecialFunctions;
using Interpolations;
using LinearAlgebra;


function gauss_quadrature(n::Int)
    if n == 1
        return 1, [0], [1];
    elseif n == 2
        return 2, [-sqrt(1/3) sqrt(1/3)], [1 1];
    elseif n == 3
        return 3, [-sqrt(3/5) 0 sqrt(3/5)], [5/9 8/9 5/9];
    else
        error("n>3 not supported");
    end
end

function fun2(x, U1, Omega1)
    interp_linear = LinearInterpolation(Omega1, U1, extrapolation_bc=Flat());
    return interp_linear(x);
end

function fun(x, U, Omega)
    if x > Omega[end]
        return U[end];
    end
    if x < Omega[1]
        return U[1];
    end
    low = 1;
    high = length(Omega);
    idx = 0;
    while low < high && low != high-1
        idx = floor(Int, (low + high) / 2);
        if Omega[idx] > x
            high = idx;
        else
            low = idx;
        end
    end
    idx = floor(Int, (low + high) / 2);
    return U[idx:idx+1] ⋅ phi(x, Omega[idx], Omega[idx+1]);
end

function phi(x, h1, h2)
    return [1-(x-h1)/(h2-h1) (x-h1)/(h2-h1)];
end

function elem_mat_sys(h1, h2)
    A_elem = [h2/3 - h1/3 h2/6 - h1/6
              h2/6 - h1/6 h2/3 - h1/3];
    return A_elem;
end
 
function elem_vec_sys(h1, h2, U, Omega)
    quad_points, quad_coords, quad_weight = gauss_quadrature(2);
    V_elem = [0 0];
    for i in 1 : quad_points
        coord = (h2-h1)/2 * quad_coords[i] + (h1 + h2)/2; # Transformation
        V_elem += quad_weight[i] * fun(coord, U, Omega) * phi(coord, h1, h2);
    end
    V_elem *= (h2 - h1) / 2; # Transformation
    return V_elem;
end

function l2_projection(Omega1, Omega2, U1)
    N = length(Omega2);
    M = zeros(N, N);
    b = zeros(N);

    for i in 1 : N - 1
        A_elem = elem_mat_sys(Omega2[i], Omega2[i+1]);
        M[i,   i]   += A_elem[1, 1];
        M[i+1, i+1] += A_elem[2, 2];
        M[i+1, i]   += A_elem[2, 1];
        M[i,   i+1] += A_elem[1, 2];

        V_elem = elem_vec_sys(Omega2[i], Omega2[i+1], U1, Omega1);
        b[i]   += V_elem[1];
        b[i+1] += V_elem[2];
    end
    U2 = M\b;

    for i in 1:length(Omega2)
        if Omega2[i] > Omega1[end]
            U2[i] = U1[end];
        end
    end
    return U2;
end
