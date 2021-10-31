module pole

using Plots
using Statistics
using Interpolations
using DataFrames
using CSV
using Unitful

function fillmissing!(M)
    for xi=1:size(M)[1]
        rys = 1:size(M)[2]
        rzs = M[xi,:]
        
        rys = [y for (y,z)=zip(rys, rzs) if !isnan(z)]
        rzs = [z for z=rzs if !isnan(z)]
        
        itp = LinearInterpolation((rys,), rzs, extrapolation_bc=Flat())
        
        for yi=1:size(M)[2]
            if isnan(M[xi,yi])
                M[xi,yi] = itp(yi)
            end
        end
    end
end

function vflip(M)
    M′ = similar(M)
    for (i, row) in enumerate(eachrow(M))
        M′[end-i+1,:] = row
    end
    M′
end

function varaouter()
    M = [
        0.84 0.73 0.60
        0.95 0.82 0.64
        1.14 0.96 0.71
        1.39 1.21 0.92
        1.73 1.59 1.40
        2.09 1.96 1.85
        2.61 2.45 2.35
        3.09 2.96 2.91
        3.56 3.58 3.39
        4.15 4.11 4.06
        4.63 4.64 4.66
        5.14 5.22 5.23
        5.70 5.78 5.81
        6.19 6.25 6.35
        6.64 6.79 6.89
        7.09 7.23 7.41
        7.50 7.65 7.85
        7.86 7.99 8.20
        8.10 8.24 8.46
        8.31 8.44 8.69
        8.46 8.60 8.72
    ]

    xs = 1:size(M)[1]
    ys = 1:size(M)[2]
    zs = M

    itp = LinearInterpolation((xs, ys), zs)
    xd = range(minimum(xs), maximum(xs), length=100)
    yd = range(minimum(ys), maximum(ys), length=100)
    zd = [itp(x, y) for x=xd, y=yd]

    heatmap(
        yd,
        xd,
        vflip(zd),
        ratio=1,
    )

    contour!(
        yd,
        xd,
        vflip(zd),
        color=:white,
    )

    savefig("output/var-A-outer-field.png")

    df(x, y) = (
        itp(size(M)[1] - y + 1, x + 1) - itp(size(M)[1] - y + 1, x),
        itp(size(M)[1] - y, x) - itp(size(M)[1] - y + 1, x)
    )

    quiver(
        reshape([y for x=xs[begin:end-1],y=ys[begin:end-1]], :),
        reshape([x for x=xs[begin:end-1],y=ys[begin:end-1]], :),
        quiver=df,
        ratio=1,
    )

    savefig("output/var-A-outer-vectors.png")

    Ex1 = (M[2,2] - M[2,1])/ustrip(u"m", 10u"mm")
    Ey1 = (M[1,1] - M[2,1])/ustrip(u"m", 10u"mm")

    Ex2 = (M[end,2] - M[end,1])/ustrip(u"m", 10u"mm")
    Ey2 = (M[end-1,1] - M[end,1])/ustrip(u"m", 10u"mm")

    println("E(1,2) = ($Ex1, $Ey1)")
    println("E(1,$(size(M)[1])) = ($Ex2, $Ey2)")
end

function varainner()
    M = [
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

    CSV.write("output/var-A-raw.csv", 
        DataFrame(M, ["V$i" for i=1:size(M)[2]])[:,[3,6,9,12,15]],
        delim=" & ",
        newline=" \\\\\n",
    )

    xs = 1:size(M)[1]
    ys = 1:size(M)[2]
    zs = M

    heatmap(
        ys,
        xs,
        vflip(zs),
        ratio=1,
    )

    savefig("output/var-A-raw-data.png")

    fillmissing!(zs)

    itp = LinearInterpolation((xs, ys), zs)
    xd = range(minimum(xs), maximum(xs), length=100)
    yd = range(minimum(ys), maximum(ys), length=100)
    zd = [itp(x, y) for x=xd, y=yd]

    heatmap(
        yd,
        xd,
        vflip(zd),
        ratio=1,
    )

    savefig("output/var-A-interp.png")

    contour!(
        yd,
        xd,
        vflip(zd),
        color=:white,
    )

    savefig("output/var-A-lines.png")

    df(x, y) = (
        itp(size(M)[1] - y + 1, x + 1) - itp(size(M)[1] - y + 1, x),
        itp(size(M)[1] - y, x) - itp(size(M)[1] - y + 1, x)
    )

    quiver(
        reshape([y for x=xs[begin:end-1],y=ys[begin:end-1]], :),
        reshape([x for x=xs[begin:end-1],y=ys[begin:end-1]], :),
        quiver=df,
        ratio=1,
    )

    savefig("output/var-A-vectors.png")

    dfV = DataFrame(:x=>[], :Vdosw=>[], :Vteor=>[])

    M = (M)u"V"
    XS = range(10u"mm", step=10u"mm", length=size(M)[1])
    d = 120u"mm"
    U = 10u"V"

    for (x, row) = zip(XS, eachrow(M))
        vavg = mean(row)
        push!(dfV, (x, vavg, x/d * U))
    end

    dfE = DataFrame(:x=>[], :Edosw=>[], :Eteor=>[])

    for i = 1:size(M)[1] - 1
        x = (XS[i] + XS[i+1])/2
        Edosw = (dfV[i+1, :Vdosw] - dfV[i, :Vdosw])/(XS[i+1] - XS[i])
        Edosw = uconvert(u"V/m", Edosw)
        Eteor = uconvert(u"V/m", U/d)
        push!(dfE, (x, Edosw, Eteor))
    end

    CSV.write("output/var-A-dfV.csv", 
        select(dfV,
            :x => (x -> ustrip.(u"mm", x)) => :x,
            :Vdosw => (x -> round.(ustrip.(u"V", x), digits=2)) => :Vdosw,
            :Vteor => (x -> round.(ustrip.(u"V", x), digits=2)) => :Vteor,
        ),
        delim=" & ",
        newline=" \\\\\n",
    )

    CSV.write("output/var-A-dfE.csv", 
        select(dfE,
            :x => (x -> ustrip.(u"mm", x)) => :x,
            :Edosw => (x -> round.(ustrip.(u"V/m", x), digits=2)) => :Edosw,
            :Eteor => (x -> round.(ustrip.(u"V/m", x), digits=2)) => :Eteor,
        ),
        delim=" & ",
        newline=" \\\\\n",
    )

    println("dfV")
    display(dfV)
    println("dfE")
    display(dfE)

    scatter(
        ustrip.(u"mm", dfV[:,:x]),
        ustrip.(u"V", dfV[:,:Vdosw]),
        xlabel="x [mm]",
        ylabel="V [V]",
        legend=false,
    )
    plot!(
        ustrip.(u"mm", dfV[:,:x]),
        ustrip.(u"V", dfV[:,:Vteor]),
    )

    savefig("output/var-A-V-vs-x.png")
    
    scatter(
        ustrip.(u"mm", dfE[:,:x]),
        ustrip.(u"V/m", dfE[:,:Edosw]),
        xlabel="x [mm]",
        ylabel="E [V/m]",
        ylims=(50,90),
        legend=false,
    )
    plot!(
        ustrip.(u"mm", dfE[:,:x]),
        ustrip.(u"V/m", dfE[:,:Eteor]),
    )

    savefig("output/var-A-E-vs-x.png")
    
    nothing
end

function varb()
    measurements = [
        7.96 8.11 7.82
        6.42 4.46 6.16
        5.28 5.23 5.03
        4.15 4.25 4.09
        3.32 3.41 3.24
        2.55 2.59 2.45
        1.86 1.94 1.86
        1.27 1.37 1.30
        0.74 0.79 0.69
    ]

    M = measurements

    CSV.write("output/var-B-raw.csv", 
        DataFrame(M, ["V$i" for i=1:size(M)[2]]),
        delim=" & ",
        newline=" \\\\\n",
    )

    M = hcat(M, M[:,1])

    xs = 1:size(M)[1]
    ys = 1:size(M)[2]
    zs = M

    itp = LinearInterpolation((xs, ys), zs)

    M = fill(NaN, 199, 199)
    cx, cy = 100, 100
    for x=1:199, y=1:199
        r = sqrt((x-cx)^2 + (y-cy)^2)
        if 35 <= r <= 83
            phi = atan(y, x)
            while phi < 0
                phi += 2pi
            end
            M[end-y+1, x] = itp(1 + (r - 35)/(83-35), 1 + phi/2pi * 3)
        end
    end

    heatmap(
        1:199,
        1:199,
        vflip(M),
        ratio=1,
    )

    savefig("output/var-B-interp.png")

    contour!(
        1:199,
        1:199,
        vflip(M),
        color=:white,
    )

    savefig("output/var-B-lines.png")

    dfx(x, y) = M[size(M)[1] - y + 1, x + 1] - M[size(M)[1] - y + 1, x]

    dfy(x, y) = M[size(M)[1] - y, x] - M[size(M)[1] - y + 1, x]

    quiver(
        reshape([round(Int, r*cos(phi) + 100) for r=36:10:83, phi=0:0.2:2pi], :),
        reshape([round(Int, r*sin(phi) + 100) for r=36:10:83, phi=0:0.2:2pi], :),
        quiver=(
            reshape([dfx(round(Int, r*cos(phi)) + 100, round(Int, r*sin(phi)) + 100) for r=36:10:83, phi=0:0.2:2pi], :),
            reshape([dfy(round(Int, r*cos(phi)) + 100, round(Int, r*sin(phi)) + 100) for r=36:10:83, phi=0:0.2:2pi], :)
        ),
        ratio=1,
    )

    savefig("output/var-B-vectors.png")

    dfV = DataFrame(:x=>[], :Vdosw=>[], :Vteor=>[])

    M = (measurements[:,1:3])u"V"
    XS = range(30u"mm", step=10u"mm", length=size(M)[1])
    rw = XS[1] - 10u"mm"
    rz = XS[end] + 10u"mm"
    U = 10u"V"

    for (x, row) = zip(XS, eachrow(M))
        r = x
        Vdosw = mean(row)
        Vteor = log(r/rz) * U / log(rz/rw)
        push!(dfV, (x, Vdosw, -Vteor))
    end

    dfE = DataFrame(:x=>[], :Edosw=>[], :Eteor=>[])

    for i = 2:size(M)[1]
        x = (XS[i] + XS[i-1])/2
        r = x
        Edosw = (dfV[i-1, :Vdosw] - dfV[i, :Vdosw])/(XS[i-1] - XS[i])
        Edosw = uconvert(u"V/m", Edosw)
        Eteor = uconvert(u"V/m", -U/(r*log(rz/rw)))
        push!(dfE, (x, -Edosw, -Eteor))
    end

    CSV.write("output/var-B-dfV.csv", 
        select(dfV,
            :x => (x -> ustrip.(u"mm", x)) => :x,
            :Vdosw => (x -> round.(ustrip.(u"V", x), digits=2)) => :Vdosw,
            :Vteor => (x -> round.(ustrip.(u"V", x), digits=2)) => :Vteor,
        ),
        delim=" & ",
        newline=" \\\\\n",
    )

    CSV.write("output/var-B-dfE.csv", 
        select(dfE,
            :x => (x -> ustrip.(u"mm", x)) => :x,
            :Edosw => (x -> round.(ustrip.(u"V/m", x), digits=2)) => :Edosw,
            :Eteor => (x -> round.(ustrip.(u"V/m", x), digits=2)) => :Eteor,
        ),
        delim=" & ",
        newline=" \\\\\n",
    )

    println("dfV")
    display(dfV)
    println("dfE")
    display(dfE)

    scatter(
        ustrip.(u"mm", dfV[:,:x]),
        ustrip.(u"V", dfV[:,:Vdosw]),
        xlabel="x [mm]",
        ylabel="V [V]",
        legend=false,
    )
    plot!(
        ustrip.(u"mm", dfV[:,:x]),
        ustrip.(u"V", dfV[:,:Vteor]),
    )

    savefig("output/var-B-V-vs-x.png")
    
    scatter(
        ustrip.(u"mm", dfE[:,:x]),
        ustrip.(u"V/m", dfE[:,:Edosw]),
        xlabel="x [mm]",
        ylabel="E [V/m]",
        legend=false,
    )
    plot!(
        ustrip.(u"mm", dfE[:,:x]),
        ustrip.(u"V/m", dfE[:,:Eteor]),
    )

    savefig("output/var-B-E-vs-x.png")

    nothing
end

function main()
    println("Variant A")
    varainner()
    varaouter()

    println("Variant B")
    varb()

    nothing
end

end
