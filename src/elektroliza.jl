module elektroliza

using Unitful
using Measurements
import PhysicalConstants.CODATA2018: N_A

function main()
    μ = 63.58u"g"
    w = 2
    t = 30u"minute" ± 5u"s"
    I = 0.5u"A" ± 3u"mA"
    m₁ = 110.246u"g" ± 0.1u"g"
    m₂ = 110.532u"g" ± 0.1u"g"

    m = m₂ - m₁
    @show m

    k = uconvert(u"mg/C", m/I/t)
    @show k

    F = uconvert(u"C", (μ/w/k))
    @show F

    e = ustrip(u"C*mol", F/N_A)u"C"
    @show e

    nothing
end

end
