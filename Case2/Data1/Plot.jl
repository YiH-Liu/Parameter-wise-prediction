using Plots
using Measures
using CSV
using DataFrames
using FilePathsBase
using LaTeXStrings

# In this script, we create Figure 3

# Get the path of the current script
path = dirname(@__FILE__)
# Read the CSV file: index for agent at initial condition
index_0 = DataFrame(CSV.File("$path\\index_0.csv")); 
# Read the CSV file: coordinate for agent at time t
index = DataFrame(CSV.File("$path\\index.csv")); 
index2 = DataFrame(CSV.File("$path\\Data2\\index.csv")); 

# Read the CSV file: count data for agent at time 100 and 600
count0 = DataFrame(CSV.File("$path\\data0.csv")); 
count1 = DataFrame(CSV.File("$path\\data.csv")); 
count2 = DataFrame(CSV.File("$path\\Data2\\data.csv")); 
count0 = count0.a
count1 = count1.a
count2 = count2.a


# Lattice size.
W = 199; H = 199; Δ = 1; δ=1

# Plot the snapshot for subpopulation 1 at initial condition.

x = 1
y = 91
width = 199
height = 19
            
# Define the corners of the rectangle
x_rect = [x, x+width, x+width, x, x]
y_rect = [y, y, y+height, y+height, y]

tfs=10

f1 = scatter((index_0.i.-1).*Δ, (index_0.j.-1).*Δ, label=false, marker=:circle,
            xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
            yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),legend = false,markersize=0.5,
            grid=false,tickfontsize = tfs,aspect_ratio = :equal,
            color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f1 = plot!(x_rect, y_rect, seriestype=:shape, lw=1.5, label=false, legend=false,fillalpha=0,linecolor=:green)
f1 = hline!([100.5],legend=false,lw = 1.5, linestyle = :dash, color=:purple)
f1 = annotate!(W*0.06, H*0.9, text(L"k = 0", :left, tfs))
f1 = annotate!(W*0.5, -H*0.1, text(L"i", :left, tfs))
f1 = annotate!(-W*0.18, H*0.5, text(L"j", :left, tfs,rotation=90))
f1 = annotate!(-W*0.3, H*1.1, text(L"(\mathrm{a})", :left, tfs))

# Plot the snapshot for subpopulation 1 at time t.
f2 = scatter((index.i.-1).*Δ, (index.j.-1).*Δ, label=false, marker=:circle,
            xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
            yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),legend = false,markersize=0.5,
            grid=false,tickfontsize = tfs,aspect_ratio = :equal,
            color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f2 = plot!(x_rect, y_rect, seriestype=:shape, lw=1.5, label=false, legend=false,fillalpha=0,linecolor=:green)
f2 = hline!([100.5],legend=false,lw = 1.5, linestyle = :dash, color=:purple)
f2 = annotate!(W*0.06, H*0.9, text(L"k = 100", :left, tfs))
f2 = annotate!(W*0.5, -H*0.1, text(L"i", :left, tfs))
f2 = annotate!(-W*0.18, H*0.5, text(L"j", :left, tfs,rotation=90))
f2 = annotate!(-W*0.3, H*1.1, text(L"(\mathrm{b})", :left, tfs))



# Plot the snapshot for subpopulation 1 at initial condition.
f3 = scatter((index_0.i.-1).*Δ, (index_0.j.-1).*Δ, label=false, marker=:circle,
            xlims = [1,200],ylims = [91,110],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),yticks=false,
           #=  yticks=([91,110], [latexstring("1"),latexstring("20")]), =#legend = false,markersize=0.5,
            grid=false,tickfontsize = tfs,aspect_ratio = :equal,
            color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f3 = plot!([1, 200, 200, 1, 1], [91, 91, 110, 110, 91], 
            label=false, lw=1.5, color=:green, legend=false)
f3 = annotate!(100, 73, text(L"i", :left, tfs))
f3 = annotate!(-W*0.3, 120, text(L"(\mathrm{d})", :left, tfs))
#= f3 = annotate!(W*0.03, 120, text(L"k = 0", :left, tfs)) =#




# Plot the snapshot for subpopulation 1 at initial condition.
f4 = scatter((index.i.-1).*Δ, (index.j.-1).*Δ, label=false, marker=:circle,
            xlims = [1,200],ylims = [91,110],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),yticks=false,
            #=  yticks=([91,110], [latexstring("1"),latexstring("20")]), =#legend = false,markersize=0.5,
            grid=false,tickfontsize = tfs,aspect_ratio = :equal,
            color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f4 = plot!([1, 200, 200, 1, 1], [91, 91, 110, 110, 91], 
label=false, lw=1.5, color=:green, legend=false)
f4 = annotate!(100, 73, text(L"i", :left, tfs))
f4 = annotate!(-W*0.3, 120, text(L"(\mathrm{e})", :left, tfs))
#= f4 = annotate!(W*0.03, 120, text(L"k = 100", :left, tfs)) =#


δ = 1;ms=1;ms2=1.5
f5 = hline([0],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f5 = hline!([20],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f5 = scatter!(1:δ:200,count0[1:δ:200], label=false, marker=:circle,
        xlims = [1,200],ylims = [-5,25],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20],[latexstring("0"),"",latexstring("20")]),legend = false,markersize=ms2,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f5 = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
f5 = annotate!(-200*0.21, 0, text(L"C^{\mathrm{o}}", :left, tfs,rotation=90))
f5 = annotate!(-W*0.3, 40, text(L"(\mathrm{g})", :left, tfs))
#= f5 = annotate!(W*0.03, 20*1.7, text(L"k = 0", :left, tfs)) =#

f6 = hline([0],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f6 = hline!([20],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f6 = scatter!(1:δ:200,count1[1:δ:200], label=false, marker=:circle,
        xlims = [1,200],ylims = [-5,25],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20],[latexstring("0"),"",latexstring("20")]),legend = false,markersize=ms2,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f6 = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
f6 = annotate!(-200*0.21, 0, text(L"C^{\mathrm{o}}", :left, tfs,rotation=90))
f6 = annotate!(-W*0.3, 40, text(L"(\mathrm{h})", :left, tfs))
#= f6 = annotate!(W*0.03, 20*1.7, text(L"k = 100", :left, tfs))
 =#
δ = 10;ms=1;ms2=1.5
f7 = hline([0],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f7= hline!([20],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f7 = scatter!(1:δ:200,count0[1:δ:200], label=false, marker=:circle,
        xlims = [1,200],ylims = [-5,25],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20],[latexstring("0"),"",latexstring("20")]),legend = false,markersize=ms2,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f7 = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
f7 = annotate!(-200*0.21, 0, text(L"C^{\mathrm{o}}", :left, tfs,rotation=90))
f7 = annotate!(-W*0.3, 40, text(L"(\mathrm{j})", :left, tfs))
#= f7 = annotate!(W*0.03, 20*1.7, text(L"k = 0", :left, tfs))
 =#
f8 = hline([0],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f8 = hline!([20],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f8 = scatter!(1:δ:200,count1[1:δ:200], label=false, marker=:circle,
        xlims = [1,200],ylims = [-5,25],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20],[latexstring("0"),"",latexstring("20")]),legend = false,markersize=ms2,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f8 = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
f8 = annotate!(-200*0.21, 0, text(L"C^{\mathrm{o}}", :left, tfs,rotation=90))
f8 = annotate!(-W*0.3, 40, text(L"(\mathrm{k})", :left, tfs))
#= f8 = annotate!(W*0.03, 20*1.7, text(L"k = 100", :left, tfs)) =#


# Plot the snapshot for subpopulation 1 at time t.
f9 = scatter((index2.i.-1).*Δ, (index2.j.-1).*Δ, label=false, marker=:circle,
            xlims = [1,200],ylims = [1,200],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
            yticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),legend = false,markersize=0.5,
            grid=false,tickfontsize = tfs,aspect_ratio = :equal,
            color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f9 = plot!(x_rect, y_rect, seriestype=:shape, lw=1.5, label=false, legend=false,fillalpha=0,linecolor=:green)
f9 = hline!([100.5],legend=false,lw = 1.5, linestyle = :dash, color=:purple)
f9 = annotate!(W*0.06, H*0.9, text(L"k = 600", :left, tfs))
f9 = annotate!(W*0.5, -H*0.1, text(L"i", :left, tfs))
f9 = annotate!(-W*0.18, H*0.5, text(L"j", :left, tfs,rotation=90))
f9 = annotate!(-W*0.3, H*1.1, text(L"(\mathrm{c})", :left, tfs))

# Plot the snapshot for subpopulation 1 at initial condition.
f10 = scatter((index2.i.-1).*Δ, (index2.j.-1).*Δ, label=false, marker=:circle,
            xlims = [1,200],ylims = [91,110],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),yticks=false,
            #=  yticks=([91,110], [latexstring("1"),latexstring("20")]), =#legend = false,markersize=0.5,
            grid=false,tickfontsize = tfs,aspect_ratio = :equal,
            color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f10 = plot!([1, 200, 200, 1, 1], [91, 91, 110, 110, 91], 
label=false, lw=1.5, color=:green, legend=false)
f10 = annotate!(100, 73, text(L"i", :left, tfs))
f10 = annotate!(-W*0.3, 120, text(L"(\mathrm{f})", :left, tfs))


δ = 1
f11 = hline([0],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f11 = hline!([20],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f11 = scatter!(1:δ:200,count2[1:δ:200], label=false, marker=:circle,
        xlims = [1,200],ylims = [-5,25],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20],[latexstring("0"),"",latexstring("20")]),legend = false,markersize=ms2,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f11 = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
f11 = annotate!(-200*0.21, 0, text(L"C^{\mathrm{o}}", :left, tfs,rotation=90))
f11 = annotate!(-W*0.3, 40, text(L"(\mathrm{i})", :left, tfs))

δ = 10
f12 = hline([0],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f12 = hline!([20],legend=false,lw = 1.5, linestyle = :dash, color=:black)
f12 = scatter!(1:δ:200,count2[1:δ:200], label=false, marker=:circle,
        xlims = [1,200],ylims = [-5,25],xticks=([1,100,200], [latexstring("1"),"",latexstring("200")]),
        yticks=([0,10,20],[latexstring("0"),"",latexstring("20")]),legend = false,markersize=ms2,
        grid=false,tickfontsize = tfs,aspect_ratio = :equal,
        color=:red,markerstrokecolor=:red,framestyle=:box,margin=0mm,left_margin=0mm)
f12 = annotate!(200*0.5, -21*1.1, text(L"i", :left, tfs))
f12 = annotate!(-200*0.21, 0, text(L"C^{\mathrm{o}}", :left, tfs,rotation=90))
f12 = annotate!(-W*0.3, 40, text(L"(\mathrm{l})", :left, tfs))


f = plot(f1,f2,f9,f3,f4,f10,f5,f6,f11,f7,f8,f12,
layout = @layout[a{0.473h} b{0.473h} c{0.473h};d{0.0764h} e{0.0764h} f{0.0764h}; g{0.1216h} h{0.1216h} i{0.1216h};j{0.1216h} k{0.1216h} l{0.1216h}], link = :x, 
top_margin = 0mm,bottom_margin=0mm, right_margin = 2mm,left_margin = 3mm,size=(680,560))  # Here we set the top margin
savefig(f,"$path\\Figure3.pdf")
display(f)
