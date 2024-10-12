using EPnP
using StaticArrays
using Rotations
using Test

# runway params
Δx = 2_000e0 # m
Δy = 50e0    # m

corners = (
    SVector(0, -Δy/2, 0),
    SVector(0, +Δy/2, 0),
    SVector(Δx, -Δy/2, 1),
    SVector(Δx, +Δy/2, 2),
    # let's for now get at least six points.
    SVector(10*rand(), 10*rand(), 10*rand()),
    SVector(10*rand(), 10*rand(), 10*rand()),
)


# some arbitrarily chosen reference points in world coordinate system
c_w = (
    corners[1],
    corners[2],
    sum([corners[3], corners[4]])/2,
    sum(corners) + SVector(0, 0, 10)
)

αs = compute_barycentric.(corners, [c_w])

cam_pos_true = SVector(-500e0 + 100*randn(), 0e0, 50e0+10*randn());
cam_rot_true_ = SVector(0.1*2pi +randn(), 0.2*2pi+randn(), 0.3*2pi+randn()) .* 0
cam_rot_true = RotXYZ(cam_rot_true_...)

projs = project.([cam_pos_true], [cam_rot_true], corners)
us = getindex.(projs, 1)
vs = getindex.(projs, 2)

rot, pos = compute_pose(us, vs, c_w, αs)
@test rot ≈ cam_rot_true_
@test pos ≈ cam_pos_true
