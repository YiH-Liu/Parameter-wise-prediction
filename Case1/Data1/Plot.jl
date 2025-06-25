using Plots
using Measures
using CSV
using DataFrames
using FilePathsBase
using LaTeXStrings
# In this script, we create Figure 2 and Supplementary Information (SI) Figure S1

# Get the path of the current script
path = dirname(@__FILE__)
t=72
t1 = 24; t2 = 48; t3 = 72

# Read the CSV file: indices for each agent in subpopulations 1 and 2 at the initial condition.
index1_0 = DataFrame(CSV.File("$path\\index1_0.csv")); 
index2_0 = DataFrame(CSV.File("$path\\index2_0.csv")); 

# Read the CSV file: indices for each agent in subpopulations 1 and 2 at time t.
index1 = DataFrame(CSV.File("$path\\index1.csv")); 
index2 = DataFrame(CSV.File("$path\\index2.csv")); 


# Read the CSV file: count data for subpopulations 1 and 2 at the initial condition.
data_0 = DataFrame(CSV.File("$path\\data_0.csv")); 
C_1_0 = data_0.a;
C_2_0 = data_0.b;

# Read the CSV file: count data for subpopulations 1 and 2 at time t.
data = DataFrame(CSV.File("$path\\data.csv")); 
C_1 = data.a;
C_2 = data.b;

# Read the CSV file: discrete model density for subpopulations 1 and 2 at the initial condition.
DiscreteDens = CSV.read("$path\\Adata.csv", DataFrame)
C1 = DiscreteDens.a
C2 = DiscreteDens.b

# Read the CSV file: continuum model density for subpopulations 1 and 2 at the initial condition.
ContinuumDens = CSV.read("$path\\solution0.csv", DataFrame)
DensC10 = ContinuumDens.a
DensC20 = ContinuumDens.b


# Read the CSV file: continuum model density for subpopulations 1 and 2 at time t.
ContinuumDens = CSV.read("$path\\solution.csv", DataFrame)
DensC1 = ContinuumDens.a
DensC2 = ContinuumDens.b
xC = 1:0.5:200


# Lattice size.
W = 199; H = 19; Δ = 1; δ=1; tfs = 8;

# Plot the snapshot for subpopulation 1 at initial condition.
ms1=0.5
s1 = scatter((index1_0.i.-1).*Δ, (index1_0.j.-1).*Δ, label=false, marker=:circle,
            color=:red,markerstrokecolor=:red,markersize=ms1,framestyle=:box,alpha=0.6)
s1 = scatter!((index2_0.i.-1).*Δ, (index2_0.j.-1).*Δ, label=false, marker=:circle,
            xlims = [1,200],ylims = [1,20],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
            yticks=([1,10,20],[latexstring("0"),"",latexstring("20")]),legend = false,markersize=ms1,
            grid=false,tickfontsize = tfs,aspect_ratio = :equal,
            color=:green,markerstrokecolor=:green,framestyle=:box,alpha=0.6)
s1 = annotate!(W*0.03, H*1.4, text(L"k = 0", :left, tfs))
s1 = annotate!(W*0.5, -H*0.6, text(L"i", :left, tfs))
s1 = annotate!(-W*0.15, H*0.5, text(L"j", :left, tfs))
s1 = annotate!(-W*0.25, H*1.4, text(L"(\mathrm{a})", :left, tfs))


# Plot the snapshot for subpopulation 1 at t.
s2 = scatter((index1.i.-1).*Δ, (index1.j.-1).*Δ, label=false, marker=:circle,
            color=:red,markerstrokecolor=:red,markersize=ms1,framestyle=:box,alpha=0.6)
s2 = scatter!((index2.i.-1).*Δ, (index2.j.-1).*Δ, label=false, marker=:circle,
            xlims = [1,200],ylims = [1,20],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
            yticks=([1,10,20],[latexstring("0"),"",latexstring("20")]),legend = false,markersize=ms1,
            grid=false,tickfontsize = tfs,aspect_ratio = :equal,
            color=:green,markerstrokecolor=:green,framestyle=:box,alpha=0.6)
s2 = annotate!(W*0.03, H*1.4, text(L"k = 100", :left, tfs))
s2 = annotate!(W*0.5, -H*0.6, text(L"i", :left, tfs))
s2 = annotate!(-W*0.15, H*0.5, text(L"j", :left, tfs))
s2 = annotate!(-W*0.25, H*1.4, text(L"(\mathrm{b})", :left, tfs))




# Plot the count data for subpopulation 1 at the initial condition.
d1a = hline([0],legend=false,lw = 1.5, linestyle = :dash, color=:black)
d1a = scatter!(1:1:200,C_1_0, label=false, marker=:circle,
        xlims = [1,200],ylims = [-10,50],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20,30,40,50],[latexstring("0"),"","","","",latexstring("50")]),legend = false,markersize=1.5,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:red,markerstrokecolor=:red,framestyle=:box)
d1a = annotate!(200*0.03, 51*1.15, text(L"k = 0", :left, tfs))
d1a = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
d1a = annotate!(-200*0.15, 20, text(L"C_{1}^{\mathrm{o}}", :left, tfs,rotation=90))
d1a = annotate!(-W*0.25, 51*1.15, text(L"(\mathrm{c})", :left, tfs))

# Plot the count data for subpopulation 2 at the initial condition.
d1b = hline([0],legend=false,lw = 1.5, linestyle = :dash, color=:black)
d1b = scatter!(1:1:200,C_2_0, label=false, marker=:circle,
        xlims = [1,200],ylims = [-10,50],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20,30,40,50],[latexstring("0"),"","","","",latexstring("50")]),legend = false,markersize=1.5,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:green,markerstrokecolor=:green,framestyle=:box)
d1b = annotate!(200*0.03, 51*1.15, text(L"k = 0", :left, tfs))
d1b = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
d1b = annotate!(-200*0.15, 20, text(L"C_{2}^{\mathrm{o}}", :left, tfs,rotation=90))
d1b = annotate!(-W*0.25, 51*1.15, text(L"(\mathrm{d})", :left, tfs))

# Plot the count data for subpopulation 1 at t.
δ = 1
ms = 1.5
d2a = hline([0],legend=false,lw = 1.5, linestyle = :dash, color=:black)
d2a = scatter!(1:δ:200,C_1[1:δ:200], label=false, marker=:circle,
        xlims = [1,200],ylims = [-10,50],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20,30,40,50],[latexstring("0"),"","","","",latexstring("50")]),legend = false,markersize=ms,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:red,markerstrokecolor=:red,framestyle=:box)
d2a = annotate!(200*0.03, 51*1.15, text(L"k = 100", :left, tfs))
d2a = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
d2a = annotate!(-200*0.15, 20, text(L"C_{1}^{\mathrm{o}}", :left, tfs,rotation=90))
d2a = annotate!(-W*0.25, 51*1.15, text(L"(\mathrm{e})", :left, tfs))

# Plot the count data for subpopulation 2 at t.
d2b = hline([0],legend=false,lw = 1.5, linestyle = :dash, color=:black)
d2b = scatter!(1:δ:200,C_2[1:δ:200], label=false, marker=:circle,
        xlims = [1,200],ylims = [-10,50],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20,30,40,50],[latexstring("0"),"","","","",latexstring("50")]),legend = false,markersize=ms,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:green,markerstrokecolor=:green,framestyle=:box)
d2b = annotate!(200*0.03, 51*1.15, text(L"k = 100", :left, tfs))
d2b = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
d2b = annotate!(-200*0.15, 20, text(L"C_{2}^{\mathrm{o}}", :left, tfs,rotation=90))
d2b = annotate!(-W*0.25, 51*1.15, text(L"(\mathrm{f})", :left, tfs))


# Plot the subset of count data for subpopulation 1 at t.
δ = 18
ms = 1.5
d3a = hline([0],legend=false,lw = 2, linestyle = :dash, color=:black)
d3a = scatter!(1:δ:200,C_1[1:δ:200], label=false, marker=:circle,
        xlims = [1,200],ylims = [-10,50],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20,30,40,50],[latexstring("0"),"","","","",latexstring("50")]),legend = false,markersize=ms,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:red,markerstrokecolor=:red,framestyle=:box)
d3a = annotate!(200*0.03, 51*1.15, text(L"k = 100", :left, tfs))
d3a = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
d3a = annotate!(-200*0.15, 20, text(L"C_{1}^{\mathrm{o}}", :left, tfs,rotation=90))
d3a = annotate!(-W*0.25, 51*1.15, text(L"(\mathrm{g})", :left, tfs))

# Plot the subset of count data for subpopulation 2 at t.
d3b = hline([0],legend=false,lw = 2, linestyle = :dash, color=:black)
d3b = scatter!(1:δ:200,C_2[1:δ:200], label=false, marker=:circle,
        xlims = [1,200],ylims = [-10,50],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20,30,40,50],[latexstring("0"),"","","","",latexstring("50")]),legend = false,markersize=ms,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:green,markerstrokecolor=:green,framestyle=:box)
d3b = annotate!(200*0.03, 51*1.15, text(L"k = 100", :left, tfs))
d3b = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
d3b = annotate!(-200*0.15, 20, text(L"C_{2}^{\mathrm{o}}", :left, tfs,rotation=90))
d3b = annotate!(-W*0.25, 51*1.15, text(L"(\mathrm{h})", :left, tfs))

f = plot(s1, s2, d1a, d1b, d2a, d2b,d3a,d3b,layout = layout = @layout[a{0.0955h} b{0.0955h};c{0.301h} d{0.301h}; e{0.301h} f{0.301h};g{0.301h} h{0.301h}], link = :x, 
top_margin = 1.5mm, bottom_margin=1mm, right_margin = 0mm,left_margin = 0mm,size=(580,350))  # Here we set the top margin
savefig(f,"$path\\Figure2.pdf")
display(f)

# Plot Supplementary Information (SI) Figure S1.

tfs = 8;ms=2;lws=2
x = 1:1:200
p1a=scatter(x,C_1_0, label=false, markersize=ms, color=:black)
p1a=plot!(xC,DensC10.*20, label=false, grid=false,
xlims=[1, 200], xticks=([1, 50, 100, 150, 200], [latexstring("1"), "", latexstring("150"), "", latexstring("200")]),
ylims=[-0.01, 50], yticks=([0,12.5,25,37.5,50], [ latexstring("0"), "", latexstring("25"), "", latexstring("50")]),
tickfontsize=tfs,color=:orange,markerstrokecolor=:orange,  framestyle=:box,linewidth=lws)
p1a = annotate!(200*0.05, 50*1.1, text(L"k = 0", :left, tfs))
p1a = annotate!(200*0.5, -50*0.3, text(L"i ", :left, tfs))
p1a = annotate!(-200*0.2, 50*0.5, text(L"C_{1}^{\mathrm{o}}", :left, tfs,rotation=90))
p1a = annotate!(-200*0.25, 50*1.1, text(L"(\mathrm{c})", :left, tfs))

x = 1:1:200
p1b=scatter(x,C_2_0, label=false, markersize=ms, color=:black)
p1b=plot!(xC,DensC20.*20, label=false, grid=false,
xlims=[1, 200], xticks=([1, 50, 100, 150, 200], [latexstring("1"), "", latexstring("150"), "", latexstring("200")]),
ylims=[-0.01, 50], yticks=([0,12.5,25,37.5,50], [ latexstring("0"), "", latexstring("25"), "", latexstring("50")]),
tickfontsize=tfs,color=:orange,markerstrokecolor=:orange, framestyle=:box,linewidth=lws)
p1b = annotate!(200*0.05, 50*1.1, text(L"k = 0", :left, tfs))
p1b = annotate!(200*0.5, -50*0.3, text(L"i ", :left, tfs))
p1b = annotate!(-200*0.2, 50*0.5, text(L"C_{2}^{\mathrm{o}}", :left, tfs,rotation=90))
p1b = annotate!(-200*0.25, 50*1.1, text(L"(\mathrm{d})", :left, tfs))

x = 1:1:200
p2a=scatter(x,C1, label=false, markersize=ms, color=:black)
p2a=plot!(xC,DensC1.*20, label=false, grid=false,
xlims=[1, 200], xticks=([1, 50, 100, 150, 200], [latexstring("1"), "", latexstring("150"), "", latexstring("200")]),
ylims=[-0.01, 50], yticks=([0,12.5,25,37.5,50], [ latexstring("0"), "", latexstring("25"), "", latexstring("50")]),
tickfontsize=tfs,color=:orange,markerstrokecolor=:orange, framestyle=:box,linewidth=lws)
p2a = annotate!(200*0.05, 50*1.1, text(L"k = 100", :left, tfs))
p2a = annotate!(200*0.5, -50*0.3, text(L"i ", :left, tfs))
p2a = annotate!(-200*0.2, 50*0.5, text(L"C_{1}^{\mathrm{o}}", :left, tfs,rotation=90))
p2a = annotate!(-200*0.25, 50*1.1, text(L"(\mathrm{e})", :left, tfs))

x = 1:1:200
p2b=scatter(x,C2, label=false, markersize=ms, color=:black)
p2b=plot!(xC,DensC2.*20, label=false, grid=false,
xlims=[1, 200], xticks=([1, 50, 100, 150, 200], [latexstring("1"), "", latexstring("150"), "", latexstring("200")]),
ylims=[-0.01, 50], yticks=([0,12.5,25,37.5,50], [ latexstring("0"), "", latexstring("25"), "", latexstring("50")]),
tickfontsize=tfs,color=:orange,markerstrokecolor=:orange,framestyle=:box,linewidth=lws)
p2b = annotate!(200*0.05, 50*1.1, text(L"k = 100", :left, tfs))
p2b = annotate!(200*0.5, -50*0.3, text(L"i ", :left, tfs))
p2b = annotate!(-200*0.2, 50*0.5, text(L"C_{2}^{\mathrm{o}}", :left, tfs,rotation=90))
p2b = annotate!(-200*0.25, 50*1.1, text(L"(\mathrm{f})", :left, tfs))


p = plot(s1,s2,p1a, p1b,p2a, p2b,layout = @layout[a{0.1h} b{0.1h};c{0.45h} d{0.4h}; e{0.4h} f{0.4h}], link = :x, 
top_margin = 1.5mm, bottom_margin=2mm, right_margin = 6mm,left_margin = 6mm,size=(570,350)) 
savefig(p,"$path\\FigureS1.pdf")
display(p)
