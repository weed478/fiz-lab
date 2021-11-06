module azot

using Plots
using DataFrames
using CSV
using Unitful
using Measurements
using Polynomials
using GLM
import PhysicalConstants.CODATA2018: R

function calcQp(dTdp)
    p = 1u"atm"
    T = 77.3u"K"

    ρ = 0.808u"g/cm^3"
    m = 1u"kg"
    V₁ = m/ρ

    M = 28u"g"
    n = (m/M)u"mol"
    V₂ = n*R*T/p

    Q = T*(V₂ - V₁)/dTdp

    Qₚ = uconvert(u"J/g", Q/m)
    @show Qₚ
end

function makecsv(atmpress, p, tup, tdown, t, pvac, tvac)
    p = @. Measurements.value(p)
    p = @. round(ustrip(u"bar", p - atmpress), digits=2)

    pvac = @. Measurements.value(pvac)
    pvac = @. round(ustrip(u"bar", pvac - atmpress), digits=2)

    atmpress = ustrip(u"bar", atmpress)
    
    tup = @. round(ustrip(u"K", tup), digits=2)
    tdown = @. round(ustrip(u"K", tdown), digits=2)
    t = @. round(ustrip(u"K", t), digits=2)
    tvac = @. round(ustrip(u"K", tvac), digits=2)

    CSV.write(
        "output/wyniki-pos.csv",
        DataFrame(
            "\$\\Delta\$p [bar]" => p,
            "p [bar]" => round.(p .+ atmpress, digits=2),
            "T\\,\$\\uparrow\$ [K]" => tup,
            "T\\,\$\\downarrow\$ [K]" => tdown,
            "T [K]" => t,
        ),
        delim=" & ",
        newline=" \\\\\n",
    )
    CSV.write(
        "output/wyniki-vac.csv",
        DataFrame(
            "\$\\Delta\$p [bar]" => pvac,
            "p [bar]" => round.(pvac .+ atmpress, digits=2),
            "T [K]" => tvac,
        ),
        delim=" & ",
        newline=" \\\\\n",
    )
end

function main()
    atmpress = 980u"hPa"

    Mpress = @. atmpress + [
        0.0 
        0.1
        0.2
        0.3
        0.4
        0.5
        0.6
        0.7
        0.8
        0.9
        1.0
    ]u"bar" ± 0.05u"bar"

    Mtempup = [
        76.7 # 20.0
        78.4 # 20.7
        79.1 # 21.0
        79.8 # 21.3
        80.5 # 21.6
        80.9 # 21.8
        81.6 # 22.1
        82.1 # 22.3
        82.6 # 22.5
        83.1 # 22.7
        83.8 # 23.0
    ]u"K"

    Mtempdown = [
        77.2 # 20.2
        78.1 # 20.6
        78.8 # 20.9
        79.5 # 21.2
        80.5 # 21.6
        80.9 # 21.8
        81.6 # 22.1
        82.1 # 22.3
        82.8 # 22.6
        83.3 # 22.8
        83.8 # 23.0
    ]u"K"

    Mtemp = @. (Mtempup + Mtempdown) / 2

    Mpressvac = @. atmpress + [
         0.0
        -0.1
        -0.2
        -0.3
        -0.4
        -0.5
        -0.6
        -0.7
        -0.8
        -0.9
        -1.0
    ]u"bar" ± 0.05u"bar"

    Mtempvac = [
        77.2 # 20.2
        76.5 # 19.9
        75.5 # 19.5
        74.6 # 19.1
        73.6 # 18.7
        72.2 # 18.1
        70.8 # 17.5
        68.6 # 16.6
        66.2 # 15.6
        64.3 # 14.8
        58.1 # 12.3
    ]u"K"

    Tice = 63.6u"K" # 14.5

    makecsv(
        atmpress,
        Mpress,
        Mtempup,
        Mtempdown,
        Mtemp,
        Mpressvac,
        Mtempvac,
    )

    scatter(
        ustrip.(u"bar", Mpress),
        ustrip.(u"K", Mtemp),
        xlabel="p [bar]",
        ylabel="T [K]",
        legend=false,
    )
    savefig("output/T-vs-p-high-press.png")

    scatter(
        ustrip.(u"bar", Mpressvac),
        ustrip.(u"K", Mtempvac),
        xlabel="p [bar]",
        ylabel="T [K]",
        legend=false,
    )
    savefig("output/T-vs-p-vacuum.png")

    @assert Mpress[1] == Mpressvac[1]
    Mpresscombined = [reverse(Mpressvac); Mpress[2:end]]
    Mtempcombined = [reverse(Mtempvac[2:end]); (Mtemp[1]+Mtempvac[1])/2; Mtemp[2:end]]

    scatter(
        ustrip.(u"bar", Mpresscombined),
        ustrip.(u"K", Mtempcombined),
        xlabel="p [bar]",
        ylabel="T [K]",
        legend=false,
    )
    savefig("output/T-vs-p-combined.png")

    x0 = ustrip(u"bar", 1u"atm")

    df = DataFrame(
        :X => ustrip.(u"bar", Measurements.value.(Mpresscombined)) .- x0,
        :Y => ustrip.(u"K", Mtempcombined),
    )
    olm = lm(@formula(Y ~ X + X^2 + X^3 + X^4 + X^5), df)
    poly = Polynomial(coef(olm))
    f(x) = poly(x - x0)

    plot!(
        ustrip.(u"bar", Measurements.value.(Mpresscombined)),
        f.(ustrip.(u"bar", Measurements.value.(Mpresscombined))),
    )

    a = poly[1]
    b = poly[0] - a*x0
    line(x) = a*x + b

    xs = [-0.1x0, 1.9x0]
    ys = line.(xs)

    plot!(
        ustrip.(u"bar", (xs)u"bar"),
        ys,
    )
    
    savefig("output/T-vs-p-combined-polynomial.png")

    a = (a)u"K/bar" ± stderr(olm)[2]u"K/bar"
    println("a₁ = $a")

    calcQp(a)

    nothing
end

end
