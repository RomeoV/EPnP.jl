module EPnP
using LinearAlgebra
using StaticArrays
using Rotations

export compute_barycentric
export project
export compute_pose

sones(N) = @SVector ones(N)

function compute_barycentric(p, refs)
  M = [hcat(refs...); sones(4)']
  M \ [p;1]
end


function project(cam_pos, cam_rot, pt; focal_length=40e-3)
    pt′ = cam_rot * (pt - cam_pos)
    focal_length .* SVector(pt′[2] / pt′[1], pt′[3] / pt′[1])
end

function compute_pose(us, vs, c_w, αs)
    u_c = 0; v_c = 0;
    f_u = f_v = 40e-3;
    # TODO: There are several strange things with this matrix.
    # Basically, it only seems to work (sometimes) if we put the z component first.
    # Then, there's still some ambiguity with the sign...
    # Very strange. Would be good to figure out soon.
    M = vcat([
        [hcat((SVector((αs[i][j] * +(us[i] - u_c)), (αs[i][j] * f_u), 0)' for j in 1:4)...);
         hcat((SVector((αs[i][j] * +(vs[i] - v_c)), 0, (αs[i][j] * f_v))' for j in 1:4)...)]
        for i in eachindex(αs)]...)
    # M = vcat([
    #     [hcat((SVector((αs[i][j] * f_u), 0, (αs[i][j] * +(us[i] - u_c)) )' for j in 1:4)...);
    #      hcat((SVector(0, (αs[i][j] * f_v), (αs[i][j] * +(vs[i] - v_c)) )' for j in 1:4)...)]
    #     for i in eachindex(αs)]...)

    # @info size(nullspace(M))
    v_flat = SVector{12}(nullspace(M)[:,1])
    v = reshape(v_flat, Size(3, 4)) |> eachcol

    β = (sum(norm(v[i] - v[j]) * norm(c_w[i] - c_w[j]) for i in 1:4, j in 1:4)/
         sum(norm(v[i] - v[j])^2                       for i in 1:4, j in 1:4))


    c_c = β .* v
    C_c = hcat(c_c...)
    C_w = hcat(c_w...)

    R_t = ([C_c' sones(4)] \ C_w')'
    idx = SVector((1:3)...)
    R = R_t[idx, idx]'
    t = R_t[idx, 4]

    rots = Rotations.params(RotXYZ(R))
    return rots, t
end

end
