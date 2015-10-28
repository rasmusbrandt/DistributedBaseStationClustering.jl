type MovieState
    num_events::Int
    event_types::Array{Symbol,1}
    descriptions::Array{ASCIIString,1}
    acting_BSs::Array{Int,1}
    other_involved_BSs::Array{Vector{Int},1}
    involved_BSs_votes::Array{BitArray,1}
    partition_evolution::Array{Partition,1}
    BS_utilities_evolution::Array{Vector{Float64},1}
end
MovieState() =
    MovieState(0, Array(Symbol, 0), Array(ASCIIString, 0), Array(Int, 0),
               Array(Vector{Int}, 0), Array(Vector{Bool}, 0), Array(Partition, 0),
               Array(Vector{Float64}, 0))

function push_event!(ms::MovieState, event_type, desc, acting_BS, other_involved_BSs, involved_BSs_votes, partition, BS_utilities)
    ms.num_events += 1
    push!(ms.event_types, event_type)
    push!(ms.descriptions, desc)
    push!(ms.acting_BSs, acting_BS)
    push!(ms.other_involved_BSs, other_involved_BSs)
    push!(ms.involved_BSs_votes, involved_BSs_votes)
    push!(ms.partition_evolution, partition)
    push!(ms.BS_utilities_evolution, BS_utilities)
end

function generate_movie(network, movie_state, dir, output_filename, frame_rate)
    I = get_num_BSs(network); K = get_num_MSs(network)

    PyPlot.rc("lines", linewidth=1., markersize=3, markeredgewidth=0.5)
    PyPlot.rc("font", size=8, family="serif", serif="Computer Modern Sans Serif")
    PyPlot.rc("text", usetex=true)
    PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
    PyPlot.rc("axes", linewidth=0.5, labelsize=8)
    PyPlot.rc("xtick", labelsize=8)
    PyPlot.rc("ytick", labelsize=8)
    PyPlot.rc("legend", fancybox=true, fontsize=8)
    PyPlot.rc("figure", figsize=(3.50,3.50), dpi=125)

    # Evolution statistics
    cell_rate_max = maximum([maximum(movie_state.BS_utilities_evolution[iter]) for iter = 1:movie_state.num_events])
    sum_rate_evol = [sum(movie_state.BS_utilities_evolution[iter]) for iter = 1:movie_state.num_events]
    text_cols = ASCIIString[]
    for iter = 1:movie_state.num_events
        if movie_state.event_types[iter] == :accepted_by_all
            push!(text_cols, "ForestGreen")
        elseif movie_state.event_types[iter] == :rejected_by_all
            push!(text_cols, "Crimson")
        elseif movie_state.event_types[iter] == :mixed_response
            push!(text_cols, "Crimson")
        else
            push!(text_cols, "Black")
        end
    end

    for iter = 1:movie_state.num_events
        fig, (ax_map, ax_text, ax_hist, ax_sum_throughput) = PyPlot.subplots(2, 2)
        fig[:set_size_inches](9.32,7)

        # Map
        local line_BS
        for i = 1:I
            pos = network.BSs[i].position
            col = "SteelBlue"
            if i == movie_state.acting_BSs[iter]
                ms = 10
            else
                ms = 6
            end
            line_BS = ax_map[:plot](pos.x, pos.y; marker="o", color=col, markeredgecolor="k", markerfacecolor=col, markersize=ms, linewidth=0, label="BS")
            ax_map[:text](pos.x + 20, pos.y + 20, "BS $i")
        end
        # local line_MS
        # for k = 1:K
        #     pos = network.MSs[k].position
        #     line_MS = ax_map[:plot](pos.x, pos.y; marker="o", color="r", markersize=2, linewidth=0, label="MS")
        # end
        ax_map[:set_xlabel]("x coordinate [m]")
        ax_map[:set_ylabel]("y coordinate [m]")
        #legend = ax_map[:legend](handles=[line_BS[1], line_MS[1]], loc="upper left", numpoints=1)
        #legend_frame = legend[:get_frame]()
        #PyPlot.setp(legend_frame, linewidth=0.5)
        ax_map[:set_aspect]("equal", "datalim") # square plot

        # Show cluster lines (replace by ellipses)
        for block in movie_state.partition_evolution[iter].blocks
            for i in block.elements
                for j in block.elements
                    if i == j; continue; end
                    pos1 = network.BSs[i].position; pos2 = network.BSs[j].position
                    ax_map[:plot]([pos1.x, pos2.x], [pos1.y, pos2.y]; color="DarkGray")
                end
            end
        end

        # Show action lines
        i = movie_state.acting_BSs[iter]
        if i != 0
            for (idx, j) in enumerate(movie_state.other_involved_BSs[iter])
                if movie_state.event_types[iter] == :accepted_by_all
                    col = "ForestGreen"
                elseif movie_state.event_types[iter] == :rejected_by_all
                    col = "Crimson"
                elseif movie_state.event_types[iter] == :mixed_response
                    if movie_state.involved_BSs_votes[iter][idx] == true
                        col = "ForestGreen"
                    else
                        col = "Crimson"
                    end
                else
                    col = "Black"
                end

                pos1 = network.BSs[i].position; pos2 = network.BSs[j].position
                ax_map[:plot]([pos1.x, pos2.x], [pos1.y, pos2.y]; color=col)
            end
        end

        # Utility histogram
        bar_width = 0.35
        ax_hist[:bar](1:I, movie_state.BS_utilities_evolution[iter], bar_width, color="SteelBlue")
        ax_hist[:set_xlabel]("Cell ")
        ax_hist[:set_xticks](1:I)
        ax_hist[:set_ylabel]("Cell throughput [bits/s/Hz]")
        ax_hist[:set_ylim](0, cell_rate_max)

        # Text
        ax_text[:set_xlim](1, 10)
        ax_text[:set_ylim](1, 10)
        for text_iter = max(1, iter-8):iter
            ax_text[:text](0, text_iter - iter + 9.5, movie_state.descriptions[text_iter], fontdict=Dict(:fontsize => 12, :color => text_cols[text_iter]))
        end
        ax_text[:axis]("off")

        # # Sum rate evolution
        ax_sum_throughput[:plot](1:iter, sum_rate_evol[1:iter], "SteelBlue")
        ax_sum_throughput[:set_xlabel]("Event")
        ax_sum_throughput[:set_ylabel]("Sum throughput [bits/s/Hz]")

        # Save as png
        fig[:savefig](joinpath(dir, @sprintf("%03d.png", iter)))
        PyPlot.close(fig)
    end

    # Generate movie with ffmpeg
    run(`ffmpeg -r $frame_rate -i $dir/%03d.png -pix_fmt yuv420p -vb 20M $output_filename.mp4`)
    run(`open $output_filename.mp4`)
end

