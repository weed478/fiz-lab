begin
    println("Adding packages, please wait...")
    import Pkg
    Pkg.add([
        "DataFrames",
        "GLM",
        "Measurements",
        "PhysicalConstants",    
        "Plots",
        "Unitful",
    ])
end

module young

import Pkg
using Plots
using DataFrames
using GLM
using Unitful
using Measurements
import Unitful: mm, cm, m, kg, N, P, GP
import PhysicalConstants.CODATA2018: g_n as g

struct Rod
    material
    l
    d
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

    # compute average deformation (/2)

    transform!(df,
        ["x↑", "x↓"] => ((xu, xd)->(xu+xd)/4) => :Δl
    )

    df
end

function linearregression(df::DataFrame, skip::Integer)
    # convert to SI floats, remove uncertainty

    df = select(df,
        :F => (x->Measurements.value.(ustrip.(N, x))) => :F,
        :Δl => (x->Measurements.value.(ustrip.(m, x))) => :Δl
    )

    # linear regression

    olm = lm(@formula(Δl ~ F), df[1+skip:end, :])
    a = coef(olm)[2] ± stderror(olm)[2]
    b = coef(olm)[1] ± stderror(olm)[1]

    (a)u"m/N", (b)u"m"
end

function getyoung(exp::ExperimentData)
    l = exp.rod.l
    d = exp.rod.d
    a, _ = linearregression(exp.df, 2)

    # convert to SI floats
    l = ustrip(u"m", l)
    d = ustrip(u"m", d)
    a = ustrip(u"m/N", a)

    # result is in pascals, convert to GP
    uconvert(GP, (4l / (pi * d^2 * a))P)
end

function getframehookcoef(exp::ExperimentData, Eᵣ)
    l = exp.rod.l
    d = exp.rod.d
    S = pi * d^2 / 4

    l = ustrip(m, l)
    S = ustrip(m^2, S)
    Eᵣ = ustrip(P, Eᵣ)

    a, _ = linearregression(exp.df, 2)
    λ = a - (l/(Eᵣ*S))m/N
    println("λₖ($(exp.rod.material)) = $λ")
    λ
end

function getcorrectedyoung(expᵣ::ExperimentData, Eᵣ, exp::ExperimentData)
    λₖ = getframehookcoef(expᵣ, Eᵣ)

    mass, F = exp.df[end, [:m, :F]]
    Δl = uconvert(mm, λₖ * F)
    println("Δlₖ($(exp.rod.material), $(Measurements.value(mass))) = $Δl")

    l = exp.rod.l
    d = exp.rod.d
    S = pi * d^2 / 4
    a, _ = linearregression(exp.df, 2)

    E = l / (S * (a - λₖ))

    uconvert(GP, ustrip(N/m^2, E)P)
end

function generateoutputs(exps::Vector{ExperimentData})
    rm("output", recursive=true, force=true)
    mkpath("output")

    for exp in exps
        println("E($(exp.rod.material)) = $(getyoung(exp))")

        a, b = linearregression(exp.df, 2)
        a = ustrip(mm/N, a) |> Measurements.value
        b = ustrip(mm, b) |> Measurements.value

        X = Measurements.value.(ustrip.(N, exp.df.F))
        Y = Measurements.value.(ustrip.(mm, exp.df.Δl))

        scatter(
            X,
            Y,
            color=:blue,
            xlabel="F [N]",
            ylabel="Δl [mm]",
            legend=false,
            title=exp.rod.material,
            xlims=(0, 1.2maximum(X)),
            ylims=(0, 1.2maximum(Y)),
            size=(400, 500),
        )

        plot!(
            (F -> a*F + b),
            color=:red,
        )

        for ext in ["png", "svg"]
            savefig("output/$(exp.rod.material)-Fvdl.$ext")
        end
    end
end

function adduncert(M::Matrix, us::Vector)
    M′ = similar(M, Any)
    for col in 1:size(M)[2]
        M′[:, col] .= M[:, col] .± us[col]
    end
    M′
end

function main()
    steelrod = Rod(
        "steel",
        106.5cm ± 1mm,   # length
        0.81mm ± 0.01mm, # diameter
    )

    brassrod = Rod(
        "brass",
        107cm ± 1mm,     # length
        1.21mm ± 0.01mm, # diameter
    )

    exp1 = ExperimentData(
        steelrod,
        DataFrame(
            adduncert(
                [
                    1kg 1.06mm 1.00mm
                    2kg 1.56mm 1.50mm
                    3kg 1.94mm 1.87mm
                    4kg 2.24mm 2.17mm
                    5kg 2.52mm 2.47mm
                    6kg 2.81mm 2.75mm
                    7kg 3.05mm 3.03mm
                    8kg 3.31mm 3.26mm
                    9kg 3.55mm 3.54mm
                   10kg 3.78mm 3.78mm
                ],
                [0kg, 0.01mm, 0.01mm]
            ),
            ["m", "x↑", "x↓"]
        ) |> computedata
    )

    exp2 = ExperimentData(
        brassrod,
        DataFrame(
            adduncert(
                [
                    1kg 0.79mm 0.75mm
                    2kg 1.14mm 1.16mm
                    3kg 1.46mm 1.47mm
                    4kg 1.71mm 1.72mm
                    5kg 1.96mm 1.97mm
                    6kg 2.20mm 2.20mm
                ],
                [0kg, 0.01mm, 0.01mm]
            ),
            ["m", "x↑", "x↓"]
        ) |> computedata
    )

    generateoutputs([exp1, exp2])

    println("E_corrected(steel)) = $(getcorrectedyoung(exp2, 100GP, exp1))")
    println("E_corrected(brass)) = $(getcorrectedyoung(exp1, 215GP, exp2))")

    nothing
end

end # module

young.main()
