struct ReferenceVariables
    length::Float64
    velocity::Float64
    time::Float64
    density::Float64
    temperature::Float64
    pressure::Float64
    R::Float64
    μ::Float64
    k::Float64
    geo::Vector{Float64}

    function ReferenceVariables(len, vel, time, density, temp, press, R, μ, k, geo)
        len > 0 || error("需给定正数的参考长度")
        density > 0 || error("需给定正数的参考密度")
        temp > 0 || error("需给定正数的参考温度")
        press > 0 || error("需给定正数的参考压强")
        new(len, vel, time, density, temp, press, R, μ, k, geo)
    end
end

# 外部构造函数：方便通过关键字参数创建
function ReferenceVariables(;
    length = -1.0,
    velocity = -1.0,
    time = -1.0,
    density = -1.0,
    temperature = -1.0,
    pressure = -1.0,
    R = -1.0,
    m = -1.0,
    μ = -1.0,
    k = -1.0,
    geo = nothing,
)
    if velocity < 0.0
        if R > 0
            velocity = (2 * R * temperature) ^ 0.5
        elseif m > 0
            velocity = (2 * 8.314 * temperature / m) ^ 0.5
        else
            error("需给定气体常数或分子质量以计算参考速度")
        end
    end

    if R < 0.0 && m > 0.0
        R = 8.31 / m
    elseif R < 0.0 && m < 0.0
        m = 29 * 0.001
        R = 8.31 / m
    end

    if μ < 0.0
        μ = density * velocity * length
    end

    if k < 0.0
        k = μ * 3.5 * R
    end

    if time < 0
        time = length / velocity
    end

    if isnothing(geo)
        geo = [1.0, 0.0, 0.0, 0.0]
    end

    ReferenceVariables(length, velocity, time, density, temperature, pressure, R, μ, k, geo)
end