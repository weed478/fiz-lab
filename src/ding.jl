module ding

using Statistics: mean
using Measurements
using Unitful
import Unitful: mm, cm, g, Hz

struct Experiment
    material
    area
    l
    m
    f
end

circle(d) = pi * (d/2)^2
square(a) = a^2

function getyoung(exp)
    f = exp.f
    ρ = exp.m / (exp.area * exp.l)

    λ = 2exp.l ./ (1:length(f))
    v = λ .* f

    E = ρ * mean(v)^2
    E = uconvert(u"GPa", E)

    println("E($(exp.material)) = $E")
end

function main()
    ud = 0.05mm
    ul = 2mm
    um = 1g
    uf = 6Hz

    exps = [
        Experiment(
            "copper",
            circle(12.05mm ± ud),
            202cm ± ul,
            2038g ± um,
            [
                943
                1893
                2836
                3779
                4723
                5672
            ]Hz .± uf,
        ),
        Experiment(
            "alu",
            square(10.00mm ± ud),
            100cm ± ul,
            273g ± um,
            [
                2484
                4963
                7441
                9926
            ]Hz .± uf,
        ),
        Experiment(
            "steel",
            circle(10.25mm ± ud),
            198cm ± ul,
            1139g ± um,
            [
                1441
                2883
                4318
                5766
            ]Hz .± uf,
        ),
        Experiment(
            "brass",
            circle(12.00mm ± ud),
            196cm ± ul,
            1869g ± um,
            [
                896
                1793
                2689
                3586
                4482
                5379
            ]Hz .± uf,
        ),
    ]

    getyoung.(exps)

    nothing
end

end
