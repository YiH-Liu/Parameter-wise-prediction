using Plots
using Measures
using CSV
using DataFrames
using FilePathsBase
using LaTeXStrings
gr()

# In this script, we create Figure S2 and Figure S3

# Get the path of the current script
path = dirname(@__FILE__)
t=150
# Read the CSV file: index for agent at initial condition
index_0 = DataFrame(CSV.File("$path\\index_0.csv")); 

# Read the CSV file: coordinate for agent at time t
index = DataFrame(CSV.File("$path\\index.csv")); 

# Read the CSV file: radial density of discrete model.
df_DiscreteDens = CSV.read("$path\\DiscreteDens.csv", DataFrame)
c_d =  df_DiscreteDens.u
r_d = df_DiscreteDens.r

# Read the CSV file: solution of PDE in radial coordinate.
df_ContinuumDens = CSV.read("$path\\ContinuumDens.csv", DataFrame)
c_c =  df_ContinuumDens.u
r_c = df_ContinuumDens.r

# Read the CSV file: solution of PDE in Cartesian coordinate.
df_continuum0 = CSV.read("$path\\A0Continuum.csv", DataFrame)
A0Continuum = Matrix(df_continuum0)
df_continuum = CSV.read("$path\\AContinuum.csv", DataFrame)
AContinuum = Matrix(df_continuum)

# Read the CSV file: density of discrete model.
df_discrete0 = CSV.read("$path\\A0Discrete.csv", DataFrame)
A0Discrete = Matrix(df_discrete0)
df_discrete = CSV.read("$path\\ADiscrete.csv", DataFrame)
ADiscrete = Matrix(df_discrete)

# difference between density of discrete model and solution of PDE in 2D
A0_diff = abs.(A0Continuum.-A0Discrete)
A_diff = abs.(AContinuum.-ADiscrete)




# Lattice size.
W = 199; H = 199; Δ = 1; δ=1


# Compare radial density
θ = LinRange(0, 2π, 500)
pf1 = scatter((index_0.i.-1).*Δ, (index_0.j.-1).*Δ, label=false, marker=:circle,
            xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
            yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),legend = false,markersize=0.5,
            grid=false,tickfontsize = tfs,aspect_ratio = :equal,
            color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
pf1 = annotate!(W*0.06, H*1.1, text(L"k = 0", :left, tfs))
pf1 = annotate!(W*0.5, -H*0.1, text(L"i", :left, tfs))
pf1 = annotate!(-W*0.18, H*0.5, text(L"j", :left, tfs,rotation=90))
pf1 = annotate!(-W*0.68, H*1.1, text(L"(\mathrm{a})", :left, tfs))

# Plot the snapshot for subpopulation 1 at time t.
pf2 = scatter((index.i.-1).*Δ, (index.j.-1).*Δ, label=false, marker=:circle,
            xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
            yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),legend = false,markersize=0.5,
            grid=false,tickfontsize = tfs,aspect_ratio = :equal,
            color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
# Plot each dashed circle
for r in range(0, 100; length = 22)
        pf2=plot!(100.5 .+ r * cos.(θ),100.5 .+ r * sin.(θ), line = (:dash, 1), color = :black, label = "")
end
pf2 = annotate!(W*0.06, H*1.1, text(L"k = 600", :left, tfs))
pf2 = annotate!(W*0.5, -H*0.1, text(L"i", :left, tfs))
pf2 = annotate!(-W*0.18, H*0.5, text(L"j", :left, tfs,rotation=90))
pf2 = annotate!(-W*0.645, H*1.1, text(L"(\mathrm{b})", :left, tfs))

p1=scatter(r_d,c_d, label=false, markersize=3, color=:black,margin=10mm,left_margin=30mm)
p1=plot!(r_c.-1,c_c, label=false, grid=false,
xlims=[0, 100], xticks=([0,25, 50, 75, 100], [latexstring("0"), "", latexstring("50"), "", latexstring("100")]),
ylims=[-0.01, 1], yticks=([ 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8, 0.9, 1], [ latexstring("0"), "", "", "", "", latexstring("0.5"), "", "", "", "", latexstring("1")]),
tickfontsize=tfs,color=:orange,markerstrokecolor=:orange,markersize=1,  legendfontsize=15,framestyle=:box,linewidth=2)
#= p1 = annotate!(200*0.05, 2*1.1, text(L"k = 100", :left, 15)) =#
p1 = annotate!(100*0.5, -1*0.23, text(L"r ", :left, tfs))
p1 = annotate!(-100*0.25, 0.3, text(L"c(r,600)", :left, tfs,rotation=90))
p1 = annotate!(-W*0.2, 1.15, text(L"(\mathrm{c})", :left, tfs))

# First contour plot
xd = 1:1:200; yd = 1:1:200
# Create the contour plot
lv = [0.25,0.75]

p2 = contour(xd, yd, ADiscrete, fill=false, framestyle=:box,levels=lv, aspect_ratio=:equal, colorbar=false,
             xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
             yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]), linewidth=2,
             margin=10mm, tickfontsize=tfs,color=:black)
p2 = contour!(xd, yd, AContinuum, fill=false, framestyle=:box,levels=lv, aspect_ratio=:equal, colorbar=false,
             xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
             yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]), linewidth=2,
             margin=10mm, tickfontsize=tfs,color=:orange)

p2 = annotate!(100, -40, text(L"i", :left, tfs))
p2 = annotate!(-60, 100, text(L"j", :left, tfs))
p2 = annotate!(-W*0.6, 200*1.15, text(L"(\mathrm{d})", :left, tfs))
p = plot(pf1,pf2,p1,p2,layout = @layout[a b;c d], link = :x, top_margin = 6mm,bottom_margin=4mm, right_margin = 10mm,left_margin = 10mm,size=(580,400))  # Here we set the top margin
savefig(p,"$path\\FigureS2.pdf")
display(p)




# Heat Map and difference
# Define the levels and color scheme
lvs = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
# Define LaTeX formatted tick labels for the color bar
colorbar_tick_labels = [latexstring("0.01"), latexstring("0.02"), latexstring("0.04"), latexstring("0.06"), latexstring("0.08")]

color_scheme = cgrad([:white,:blue],lvs)
# Define LaTeX formatted tick labels for the color bar
colorbar_tick_labels = [latexstring("0.01"), latexstring("0.02"), latexstring("0.04"), latexstring("0.06"), latexstring("0.08")]

x_range = 1:1:200
y_range = 1:1:200


s1a1 = heatmap(x_range, y_range, A0Discrete, fill=true, levels=lvs, framestyle=:box, color=color_scheme, aspect_ratio=:equal, colorbar=false,
xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]), linewidth=1, linecolor=:transparent,
margin=10mm, tickfontsize=tfs,clims=(0,1),grid=true)
s1a1 = annotate!(100, -50, text(L"i ", :left, tfs))
s1a1 = annotate!(-70, 100, text(L"j", :left, tfs,rotation=90))
s1a1 = annotate!(-W*0.5, H*1.2, text(L"(\mathrm{a})", :left, tfs))
s1a1 = plot!([1, 200, 200, 1, 1], [1, 1, 200, 200, 1], seriestype=:shape, lw=0.5, label=false, legend=false,fillalpha=0,linecolor=:black)



s1a2 = heatmap(x_range, y_range, A0Continuum, fill=true, levels=lvs, framestyle=:box,  color=color_scheme, aspect_ratio=:equal, colorbar=false,
xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]), linewidth=1, linecolor=:transparent,
margin=10mm, tickfontsize=tfs,clims=(0,1),grid=true)
s1a2 = annotate!(100, -50, text(L"i ", :left, tfs))
s1a2 = annotate!(-70, 100, text(L"j", :left, tfs,rotation=90))
s1a2 = annotate!(-W*0.5, H*1.2, text(L"(\mathrm{b})", :left, tfs))
s1a2 = plot!([1, 200, 200, 1, 1], [1, 1, 200, 200, 1], seriestype=:shape, lw=0.5, label=false, legend=false,fillalpha=0,linecolor=:black)


s1a3 = heatmap(x_range, y_range, A0_diff, fill=true, levels=lvs, framestyle=:box, color=color_scheme, aspect_ratio=:equal, colorbar=false,
xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]), linewidth=1, linecolor=:transparent,
margin=10mm, tickfontsize=tfs,clims=(0,1))
s1a3 = annotate!(100, -50, text(L"i ", :left, tfs))
s1a3 = annotate!(-70, 100, text(L"j", :left, tfs,rotation=90))
s1a3 = annotate!(-W*0.5, H*1.2, text(L"(\mathrm{c})", :left, tfs))
s1a3 = plot!([1, 200, 200, 1, 1], [1, 1, 200, 200, 1], seriestype=:shape, lw=0.5, label=false, legend=false,fillalpha=0,linecolor=:black)


s1b1 = heatmap(x_range, y_range, ADiscrete, fill=true, levels=lvs, framestyle=:box, color=color_scheme, aspect_ratio=:equal, colorbar=false,
xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]), linewidth=1, linecolor=:transparent,
margin=10mm, tickfontsize=tfs,clims=(0,1),grid=true)
s1b1 = annotate!(100, -50, text(L"i ", :left, tfs))
s1b1 = annotate!(-70, 100, text(L"j", :left, tfs,rotation=90))
s1b1 = annotate!(-W*0.5, H*1.2, text(L"(\mathrm{d})", :left, tfs))
s1b1 = plot!([1, 200, 200, 1, 1], [1, 1, 200, 200, 1], seriestype=:shape, lw=0.5, label=false, legend=false,fillalpha=0,linecolor=:black)


s1b2 = heatmap(x_range, y_range, AContinuum, fill=true, levels=lvs, framestyle=:box,  color=color_scheme, aspect_ratio=:equal, colorbar=false,
xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]), linewidth=4,
margin=10mm, tickfontsize=tfs,clims=(0,1),grid=true)
s1b2 = annotate!(100, -50, text(L"i ", :left, tfs))
s1b2 = annotate!(-70, 100, text(L"j", :left, tfs,rotation=90))
s1b2 = annotate!(-W*0.5, H*1.2, text(L"(\mathrm{e})", :left, tfs))
s1b2 = plot!([1, 200, 200, 1, 1], [1, 1, 200, 200, 1], seriestype=:shape, lw=0.5, label=false, legend=false,fillalpha=0,linecolor=:black)


s1b3 = heatmap(x_range, y_range, A_diff, fill=true, levels=lvs, framestyle=:box, color=color_scheme, aspect_ratio=:equal, colorbar=false,
xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]), linewidth=1, linecolor=:transparent,
margin=10mm, tickfontsize=tfs,clims=(0,1))
s1b3 = annotate!(100, -50, text(L"i ", :left, tfs))
s1b3 = annotate!(-70, 100, text(L"j", :left, tfs,rotation=90))
s1b3 = annotate!(-W*0.5, H*1.2, text(L"(\mathrm{f})", :left, tfs))
s1b3 = plot!([1, 200, 200, 1, 1], [1, 1, 200, 200, 1], seriestype=:shape, lw=0.5, label=false, legend=false,fillalpha=0,linecolor=:black)


s = plot(s1a1, s1a2, s1a3,s1b1, s1b2, s1b3, layout = (2,3), link = :x, top_margin = 0mm, right_margin = 0mm,left_margin = 0mm,size=(650,500))  # Here we set the top margin


lvs = [1,100,200,300,400,500,600,700,800,900,1000]
# Define LaTeX formatted tick labels for the color bar

color_scheme = cgrad([:white,:blue],lvs)
dummy_data = reshape(1:1000, 1, 1000) # This is just to generate the color scale
p_colorbar = heatmap(dummy_data, color=color_scheme,
xticks=([1,100,200,300,400,500,600,700,800,900,1000], [latexstring("0"),latexstring("0.1"), latexstring("0.2"),latexstring("0.3"), latexstring("0.4"),latexstring("0.5"), latexstring("0.6"),latexstring("0.7"), latexstring("0.8"),latexstring("0.9"),latexstring("1")]),
tickfontsize=tfs, legend=false,size=[1800,100],bottom_margin=0mm,yticks=nothing)

c1 = plot(s,p_colorbar,layout = @layout[a{0.95h};b{0.05h}], link = :x,bottom_margin=5mm, top_margin = 5mm, right_margin = 6mm,left_margin =5mm,size=(590,450))  # Here we set the top margin (1800,1250)
savefig(c1,"$path\\FigureS3.pdf")
display(c1)

