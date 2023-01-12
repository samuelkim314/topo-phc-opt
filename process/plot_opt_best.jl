# Plot optimization curves, i.e., best value of objective found so far as a function of iteration.
using DelimitedFiles
using Printf

nin = 4
sgnum = 13
D = 3

# prefix = "opt-local"
prefix = "local2-4-579-best"
out_dir_top = "output/sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))"


obj_best = []
open(joinpath(out_dir_top,"$prefix/obj-best.txt")) do io
    for line in eachline(io)
        splitline = split(line, '\t')
        push!(obj_best, parse.(Float64, splitline))
    end
end


using PlotlyJS
traces = Vector{AbstractTrace}()
for obj_best_i in obj_best
    push!(traces, PlotlyJS.scatter(x=1:length(obj_best_i), y=obj_best_i))
end
p = PlotlyJS.plot(traces,
    Layout(xaxis_title_text="N",
        yaxis_title_text="Complete gap",
        # yaxis_range=[-0.5,-0.05]
        ))
# PlotlyJS.savefig(p, joinpath(out_dir_top, "$prefix-best.png"))
display(p)


# # Overlay markers on the curve, like Fig 3(b) in the paper
# using PyPlot
# PyPlot.figure()
# for line in obj_best
#     PyPlot.plot(line .* 100) 
# end
# # localbestx = [444, 1000, 1497, 1996, 2464]
# localbestx = [1, 233, 1180, 1438, 2102]
# PyPlot.plot(localbestx, obj_best[1][localbestx]*100, "o", c="r")
# PyPlot.ylabel("Complete gap [%]")
# PyPlot.xlabel("Iteration")
# PyPlot.display_figs()
# # PyPlot.savefig(joinpath(out_dir_top, "$prefix/best.svg"))