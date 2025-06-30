function transform_to_direct_forcing(ps::KitBase.PSpace2D, ctr::KitBase.AM, pd_out::Matrix, ib::IBM2D, Δt)
    # 与ib的参考变量    
    tem_ref = ib.rv.temperature
    rho_ref = ib.rv.density
    l_ref = ps.dx[1] * ib.rv.length
    t_ref = Δt
    u_ref = l_ref / t_ref
    force_ref = rho_ref * l_ref / t_ref^2

    xf = OffsetArray(vec(ps.x), 0:size(ctr,1) * size(ctr,2)-1) .* ib.rv.length ./ l_ref
    yf = OffsetArray(vec(ps.y), 0:size(ctr,1) * size(ctr,2)-1) .* ib.rv.length ./ l_ref
    rhof = OffsetArray(zeros(size(ctr,1) * size(ctr,2)), 0:size(ctr,1) * size(ctr,2)-1)
    ux = OffsetArray(zeros(size(ctr,1) * size(ctr,2)), 0:size(ctr,1) * size(ctr,2)-1)
    uy = OffsetArray(zeros(size(ctr,1) * size(ctr,2)), 0:size(ctr,1) * size(ctr,2)-1)
    Tf = OffsetArray(zeros(size(ctr,1) * size(ctr,2)), 0:size(ctr,1) * size(ctr,2)-1)
    for i in axes(ctr, 1)
        for j in axes(ctr, 2)
            n = i + j * size(ctr, 1)
            rhof[n] = ctr[i, j].prim[1] * ib.rv.density / rho_ref
            ux[n] = ctr[i, j].prim[2] * ib.rv.velocity / u_ref
            uy[n] = ctr[i, j].prim[3] * ib.rv.velocity / u_ref
            Tf[n] = ctr[i, j].prim[4] * ib.rv.temperature / tem_ref
        end
    end

    xpd = OffsetArray(pd_out[1, :] ./ l_ref, 0:size(pd_out, 2)-1)
    ypd = OffsetArray(pd_out[2, :] ./ l_ref, 0:size(pd_out, 2)-1)
    fxpd = OffsetArray(pd_out[9, :] ./ force_ref, 0:size(pd_out, 2)-1)
    fypd = OffsetArray(pd_out[10, :] ./ force_ref, 0:size(pd_out, 2)-1)
    tempd = OffsetArray((pd_out[7,:] .+ ib.rv.temperature) ./ tem_ref, 0:size(pd_out, 2)-1)

    return xf, yf, rhof, ux, uy, Tf, xpd, ypd, fxpd, fypd, tempd
end