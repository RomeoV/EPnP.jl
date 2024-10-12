using LinearAlgebra
using StaticArrays
using StatsBase
using Rotations
using SimpleNonlinearSolve

# runway params
Δx = 2_000e0 # m
Δy = 50e0    # m

corners = [
    SVector(0, -Δy/2, 0),
    SVector(0, +Δy/2, 0),
    SVector(Δx, -Δy/2, 1),
    SVector(Δx, +Δy/2, 2),
    # let's for now get at least six points.
    SVector(10*rand(), 10*rand(), 10*rand()),
    SVector(10*rand(), 10*rand(), 10*rand()),
    SVector(10*rand(), 10*rand(), 10*rand()),
    SVector(10*rand(), 10*rand(), 10*rand()),
    SVector(10*rand(), 10*rand(), 10*rand()),
]


# some arbitrarily chosen reference points in world coordinate system
c_w = (
    corners[1],
    corners[2],
    mean([corners[3], corners[4]]),
    mean(corners) + SVector(0, 0, 10)
)

function compute_barycentric(p, refs)
  M = [stack(refs, dims=2); ones(length(refs))']
  M \ [p;1]
end

αs = compute_barycentric.(corners, [c_w])

cam_pos_true = SVector(-500e0, 0e0, 50e0);
cam_rot_true = RotXYZ(0.1, 0.2, 0.3)

function project(cam_pos, pt; focal_length=40e-3)
    pt′ = cam_rot_true * (pt - cam_pos)
    focal_length .* SVector(pt′[2] / pt′[1], pt′[3] / pt′[1])
end

projs = project.([cam_pos_true], corners)
us = getindex.(projs, 1)
vs = getindex.(projs, 2)

u_c = 0; v_c = 0;
f_u = f_v = focal_length = 40e-3;
M = vcat([
    [hcat([[(αs[i][j] * -(us[i] - u_c)) (αs[i][j] * f_u) 0] for j in 1:4]...);
     hcat([[(αs[i][j] * -(vs[i] - v_c)) 0 (αs[i][j] * f_v)] for j in 1:4]...)]
    for i in eachindex(αs)]...)

v_flat = nullspace(M)[:,1]
v = reshape(v_flat, 3, :) |> eachcol

β = (sum(norm(v[i] - v[j]) * norm(c_w[i] - c_w[j]) for i in 1:4, j in 1:4)/
     sum(norm(v[i] - v[j])^2                       for i in 1:4, j in 1:4))


c_c = β .* v
C_c = stack(c_c; dims=2)
C_w = stack(c_w; dims=2)

# R_t * [C_c ; 1] = C_w
# but we can't invert from the left. so we transpose:
# [C_c ; 1]' * R_t' = C_w'
# solve for R_t', and then transpose back.

R_t = ([C_c' ones(4)] \ C_w')'
R = R_t[1:3, 1:3]
t = R_t[1:3, end]

Rotations.params(RotXYZ(R')) ./ (2*pi)


# function f(x, ps)
#     R = RotXYZ(x[1:3]...)
#     t = x[4:6]
#     mean.([R] .* (c_w .- [t]) .- c_c)
# end

# prob = NonlinearLeastSquaresProblem{false}(f, [0;0;0;-500;0;50])
# res = solve(prob, SimpleNewtonRaphson(); maxiters=1e5)
