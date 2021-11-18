module kompas

using Plots
using Statistics: mean
using Measurements
using Measurements: value
using Unitful
using DataFrames
import CSV
import Unitful: mm, mA, °, μT
import PhysicalConstants.CODATA2018: μ_0 as μ₀

B₀(R, N, I, α) = uconvert(μT, μ₀ * (N * I) / (2 * R * tan(α)))

function plotstuff(df)
    N = unique(df.N)
    dfs = map(N) do n
        filter(df) do row
            row.N == n
        end
    end

    plot(
        xlabel="I [mA]",
        ylabel="α [°]",
        ylims=(-Inf, 100),
    )
    for df in dfs |> reverse
        plot!(
            ustrip.(mA, df.I .|> value),
            ustrip.(°, df.α .|> value),
            label="N=$(df.N[1])",
        )
    end
    savefig("output/I-vs-alpha")

    plot(
        xlabel="I [mA]",
        ylabel="B₀ [μT]",
        ylims=(0, 35),
    )
    for df in dfs
        scatter!(
            ustrip.(mA, df.I .|> value),
            ustrip.(μT, df.B₀ .|> value),
            label="N=$(df.N[1])",
        )
    end
    savefig("output/I-vs-B_0")
    nothing
end

function main()
    d = 260mm ± 3mm
    R = d / 2

    df = DataFrame(
        [
            12	100	16	15
            12	200	27	28
            12	300	40	39
            12	400	47	46
            12	500	54	51
                        
            16	100	19	18
            16	200	35	34
            16	300	46	44
            16	400	55	54
            16	500	62	60
                        
            24	100	28	28
            24	200	46	46
            24	300	58	58
            24	400	66	65
            24	500	70	69
                        
            36	100	40	37
            36	200	60	58
            36	300	69	66
            36	400	72	73
            36	500	76	77
                        
            40	100	42	44
            40	200	60	61
            40	300	69	71
            40	400	73	74
            40	500	77	78
        ] .* [1 mA ° °],
        [:N, :I, :left, :right]
    )

    transform!(df,
        [:I] => (I -> I .± 3mA) => :I,
        [:left] => (a -> a .± 2°) => :left,
        [:right] => (a -> a .± 2°) => :right,
    )

    transform!(df, [:left, :right] => ((a, b) -> (a .+ b) ./ 2) => :α)

    transform!(df, [:N, :I, :α] => ((N, I, α) -> B₀.(R, N, I, α)) => :B₀)

    display(df)

    plotstuff(df)

    println("B₀ = $(mean(df.B₀))")

    CSV.write("output/kompas.csv",
        select(df,
            :N => :N,
            :I => ByRow(x -> round(Int, ustrip(mA, value(x)))) => "I [mA]",
            :left => ByRow(x -> round(Int, ustrip(°, value(x)))) => "α lewo [°]",
            :right => ByRow(x -> round(Int, ustrip(°, value(x)))) => "α prawo [°]",
            :α => ByRow(x -> round(ustrip(°, value(x)), digits=2)) => "α [°]",
            :B₀ => ByRow(x -> round(ustrip(μT, value(x)), digits=1)) => "B₀ [μT]",
            :B₀ => ByRow(x -> round(ustrip(μT, Measurements.uncertainty(x)), digits=1)) => "u(B₀) [μT]",
        ),
        newline=" \\\\\n\\midrule\n",
        delim=" & ",
    )

    nothing
end

end # module
