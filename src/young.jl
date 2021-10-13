module young

import Pkg
using Plots
using DataFrames
using Unitful
import Unitful: cm, mm, kg, N
import PhysicalConstants.CODATA2018: g_n as g

struct Rod
    material
    l
    ul
    d
    ud
end

struct ExperimentData
   rod::Rod
   df::DataFrame
end

function computedata(df::DataFrame)
    # compute force

    insertcols!(df, 2,
        :F => map(m->uconvert(N,m*g), df[!, :m])
    )

    # compute absolute deformation (up and down)

    xdown = [df[1, "Δx↓"]]
    
    for Δx in df[2:end, "Δx↓"]
        push!(xdown, last(xdown) + Δx)
    end

    xup = [last(xdown)]

    for Δx in df[end:-1:begin+1, "Δx↑"]
        push!(xup, last(xup) - Δx)
    end

    reverse!(xup)

    insertcols!(df, "x↑" => xup, "x↓" => xdown)

    # compute average deformation

    transform!(df,
        ["x↑", "x↓"] => ((xu, xd)->(xu+xd)/4) => :Δl
    )

    df
end

function main()
    steelrod = Rod(
        "steel", # material
        106.5cm, # length
        1mm,     # u(l)
        0.81mm,  # diameter
        0.01mm   # u(d)
    )

    brassrod = Rod(
        "brass", # material
        107cm,   # length
        1mm,     # u(l)
        1.21mm,  # diameter
        0.01mm   # u(d)
    )

    exp1 = ExperimentData(
        steelrod,
        DataFrame(
            [
                1kg 1.00mm 1.00mm
                2kg 0.50mm 0.50mm
                3kg 0.38mm 0.37mm
                4kg 0.30mm 0.30mm
                5kg 0.28mm 0.30mm
                6kg 0.29mm 0.28mm
                7kg 0.24mm 0.28mm
                8kg 0.26mm 0.23mm
                9kg 0.24mm 0.28mm
               10kg 0.23mm 0.24mm
            ],
            ["m", "Δx↑", "Δx↓"]
        ) |> computedata
    )

    exp2 = ExperimentData(
        brassrod,
        DataFrame(
            [
                1kg 0.72mm 0.75mm
                2kg 0.35mm 0.41mm
                3kg 0.32mm 0.31mm
                4kg 0.25mm 0.25mm
                5kg 0.25mm 0.25mm
                6kg 0.24mm 0.23mm
            ],
            ["m", "Δx↑", "Δx↓"]
        ) |> computedata
    )

    plts = map([exp1, exp2]) do exp
        plot(
            exp.df.F .|> ustrip,
            exp.df.Δl .|> ustrip,
            xlabel="F [N]",
            ylabel="Δl [mm]",
            legend=false,
            title=exp.rod.material,
        )
    end

    println(typeof(plts))

    plot(plts...)
end

function prepareenv()
    Pkg.add([
        "DataFrames",
        "Measurements",
        "PhysicalConstants",    
        "Plots",
        "Unitful",
    ])
end

end # module
