using MPBUtils, Crystalline, PlotlyJS
using Printf, DelimitedFiles
include("../src/bands.jl")
include("../src/write_ctl.jl")

# sgnum = 81
# D = 3
# calcname = "4384"
# doscalcname = "sg81-9-3098-4384"
# bandsbelow = 4
sgnum = 82
# calcname = "4803"
# doscalcname = "sg82-4-3886-$calcname"
# bandsbelow = 6
calcname = "4920-16bands"
doscalcname = "sg82-8-4882-4920"
bandsbelow = 7
# calcname = "3136-16bands"
# doscalcname = "sg82-2-4911-3136"
# calcname = "4997"
# doscalcname = "sg82-1-4286-4997"
# bandsbelow = 6


cntr = centering(sgnum, 3) # will be ‘C’ for base-centered

in_dir = "sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))/ff-direct/rand-local-8-4882-best"
out_dir = in_dir
in_dir_full = joinpath("output", in_dir)
out_dir_full = joinpath("output", out_dir)


completegap = calccompletebandgap(joinpath(out_dir, calcname), bandsbelow, dir="output")
gap_bounds = getgap(joinpath(out_dir, calcname), bandsbelow, dir="output")

# Extract lattice vectors to extract k-points
in_path = joinpath(out_dir_full, calcname * ".sh")
pRs, pflat, isoval, epsin, epsout, kvecs = lattice_from_mpbparams(in_path)
Rs = conventionalize(pRs, cntr)

pbands = plotbands(joinpath(out_dir, calcname), Rs, sg=sgnum, dim=3,
    mu1=bandsbelow,
    title="$doscalcname    Compgap: $(@sprintf("%.4f", completegap))"*
        "    Topological index: [1]    Bands below: $bandsbelow",
    nkpoints=300, 
    width=600, height=250,
    ylim=[0.46, 0.6], 
    range_fill=[gap_bounds[1], gap_bounds[2]], 
    hlinesy=[0.535],
    savename=joinpath(out_dir_full, "bands-zoom-$calcname.svg")
)

dosdata = readdlm(joinpath("fullbz/DOS-calculation/data", doscalcname*"-output.txt"))
dosfreq = dosdata[:, 1]
dos = dosdata[:, 2]

# println(pbands.plot.layout)
# asdf

pdos = PlotlyJS.scatter(x=dos, y=dosfreq, mode="lines", line_color="black")
p = PlotlyJS.plot(pdos, 
    # Layout(plot_bgcolor=layout[:plot_bgcolor],
    #     yaxis=layout[:yaxis],
    #     xaxis_showline=true,
    #     xaxis_ticks="outside",
    #     xaxis_linecolor="black",
    #     xaxis_range=[0, 65]))
    #     p = PlotlyJS.plot(pdos, 
    merge!(pbands.plot.layout, Layout(
        xaxis_showline=true,
        xaxis_ticks="outside",
        xaxis_linecolor="black",
        xaxis_range=[0, 65],
        xaxis_mirror=true,
        # xaxis_title="DOS",
        title="$doscalcname")))
PlotlyJS.savefig(p, joinpath("fullbz/DOS-calculation/data", doscalcname*"-dos.svg"), width=300, height=250, scale=2)

# println(fieldnames(pbands.plot))
# println(p.plot.layout)

# traces = append!(p.plot.data, [pdos])
# PlotlyJS.plot(traces, p.plot.layout)

# PlotlyJS.relayout!([PlotlyJS.plot(pdos), PlotlyJS.plot(pdos)])

# p = make_subplots(rows=1, cols=2, shared_yaxes=true)
# for tracei in pbands.plot.data
#     add_trace!(p, tracei, row=1, col=1)
# end
# relayout!(p, layout=pbands.plot.layout)
# add_trace!(p, pdos, row=1, col=2)

# p

# p = add_trace!(pbands, pdos)
# layout = pbands.plot.layout
# layout.subplots = Subplots(cols=2, shared_yaxes=true)
# dump(layout)
# relayout!(pbands, subplots.cols=2)
# add_trace!(pbands, pdos, row=1, col=2)
# layout[Symbol("xaxis4")] = copy(get(layout, :xaxis2, attr()))
# layout[Symbol("xaxis4_domain")] = [50, 100]
# relayout!(pbands, xaxis4= copy(get(layout, :xaxis2, attr())))
# relayout!(pbands, layout=layout)
# pbands
# pbands