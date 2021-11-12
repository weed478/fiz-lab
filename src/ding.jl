module ding

using Statistics: mean
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
    exps = [
        Experiment(
            "copper",
            circle(12.05mm),
            202cm,
            2038g,
            [
                943
                1893
                2836
                3779
                4723
                5672
            ]Hz,
        ),
        Experiment(
            "alu",
            square(10.00mm),
            100cm,
            273g,
            [
                2484
                4963
                7441
                9926
            ]Hz,
        ),
        Experiment(
            "steel",
            circle(10.25mm),
            198cm,
            1139g,
            [
                1441
                2883
                4318
                5766
            ]Hz,
        ),
        Experiment(
            "brass",
            circle(12.00mm),
            196cm,
            1869g,
            [
                896
                1793
                2689
                3586
                4482
                5379
            ]Hz,
        ),
    ]

    getyoung.(exps)

    nothing
end

end
