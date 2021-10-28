module pole

using PlotlyJS
using Statistics
using Interpolations

function fillmissing!(M)
    for xi=1:size(M)[1]
        val = mean(filter(x->!isnan(x), M[xi, :]))
        for yi=1:size(M)[2]
            if isnan(M[xi,yi])
                M[xi,yi] = val
            end
        end
    end
end

function main()
    data1 = [
        NaN  NaN  1.04 NaN  NaN  0.86 NaN  NaN  0.75 NaN  NaN  0.94 NaN  NaN  0.95 NaN 
        NaN  NaN  1.79 NaN  NaN  1.71 NaN  NaN  1.69 NaN  NaN  1.72 NaN  NaN  1.74 NaN
        2.84 2.73 2.58 2.49 2.36 2.44 2.56 2.42 2.50 2.56 2.52 2.47 2.52 2.53 2.61 2.78
        NaN  NaN  3.27 NaN  NaN  3.15 NaN  NaN  3.26 NaN  NaN  3.30 NaN  NaN  3.31 NaN
        NaN  NaN  3.92 NaN  NaN  3.92 NaN  NaN  4.09 NaN  NaN  4.06 NaN  NaN  4.08 NaN
        4.64 4.67 4.71 4.65 4.63 4.80 4.88 4.86 4.95 4.94 4.97 4.94 4.82 4.82 4.85 4.88
        NaN  NaN  5.35 NaN  NaN  5.49 NaN  NaN  5.80 NaN  NaN  5.78 NaN  NaN  5.64 NaN
        NaN  NaN  6.05 NaN  NaN  6.27 NaN  NaN  6.56 NaN  NaN  6.56 NaN  NaN  6.36 NaN
        6.50 6.58 6.68 7.03 7.07 7.11 7.17 7.19 7.33 7.36 7.43 7.31 7.24 7.24 7.14 7.03
        NaN  NaN  7.43 NaN  NaN  7.88 NaN  NaN  8.08 NaN  NaN  8.17 NaN  NaN  7.95 NaN
        NaN  NaN  8.11 NaN  NaN  8.62 NaN  NaN  8.90 NaN  NaN  8.93 NaN  NaN  8.65 NaN
    ]

    fillmissing!(data1)

    xs = 1:size(data1)[1]
    ys = 1:size(data1)[2]
    zs = data1

    itp = LinearInterpolation((xs, ys), zs)
    xd = range(minimum(xs), maximum(xs), length=100)
    yd = range(minimum(ys), maximum(ys), length=100)
    zd = [itp(x, y) for x=xd, y=yd]

    heatmap(
        x=xd,
        y=yd,
        z=zd,
    ) |> plot |> display

    nothing
end

end