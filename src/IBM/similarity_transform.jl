function compute_similarity_transform!(data::Dict{Int, Vector{Vector{Vector{Float64}}}}, trans_para::Vector{Float64})
    length(trans_para) == 4 || error("变换参数需要 [s, θ, tx, ty]")
    s, θ, tx, ty = trans_para
    
    # 构造 3D 旋转矩阵
    R = [cos(θ) -sin(θ) 0;
         sin(θ)  cos(θ) 0;
         0       0      1]
    
    # 对每个边界边应用变换
    for (_, edges) in data
        for edge in edges
            for point in edge
                # 转换为列向量
                v = [point[1], point[2], point[3]]
                
                # 应用变换
                v_transformed = s * R * v + [tx, ty, 0.0]
                
                # 更新原数据
                point[1], point[2], point[3] = v_transformed[1], v_transformed[2], v_transformed[3]
            end
        end
    end
    return data
end

function compute_similarity_transform(data::Dict{Int, Vector{Vector{Vector{Float64}}}}, trans_para::Vector{Float64})
    new_data = deepcopy(data)
    compute_similarity_transform!(new_data, trans_para)
    return new_data
end

function compute_similarity_transform(points::Matrix{Float64}, trans_para::Vector{Float64})
    length(trans_para) == 4 || error("变换参数需要 [s, θ, tx, ty]")
    s, θ, tx, ty = trans_para
    
    # 构造 3D 旋转矩阵（绕 Z 轴旋转）
    R = [cos(θ) -sin(θ) 0;
         sin(θ)  cos(θ) 0;
         0       0      1]
    
    # 应用变换：缩放 → 旋转 → 平移
    sR = s * R
    transformed = sR * points .+ [tx, ty, 0]

    return transformed
end
