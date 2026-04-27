using DelimitedFiles
using ArgParse
using Printf
using Interpolations

function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--scale", "-s"
            help = "scale factor from inches"
            arg_type = Float64
            default = 1.0
        "input"
            help = "input PE0 file"
            arg_type = String
            required = true
        "output"
            help = "output VSP BEM file"
            arg_type = String
            required = true
    end

    parsed_args = parse_args(s)

    infile = parsed_args["input"]
    outfile = parsed_args["output"]
    scale = parsed_args["scale"]

    # Read the entire file as text
    content = read(infile, String)

    # Extract RADIUS value
    radius_match = match(r"RADIUS:\s+([\d.]+)", content)
    if radius_match === nothing
        error("Could not find RADIUS in PE0 file")
    end
    R = parse(Float64, radius_match.captures[1])

    # Extract BLADES value
    blades_match = match(r"BLADES:\s+(\d+)", content)
    if blades_match === nothing
        error("Could not find BLADES in PE0 file")
    end
    num_blade = parse(Int, blades_match.captures[1])

    # Extract the geometry table
    lines = split(content, '\n')

    # Find the line with "STATION" header (the line that starts the geometry table header)
    header_idx_local = 0
    for (i, line) in enumerate(lines)
        if contains(line, "STATION") && contains(line, "CHORD")
            header_idx_local = i
            break
        end
    end

    if header_idx_local == 0
        error("Could not find geometry table header (STATION/CHORD) in PE0 file")
    end

    println("✓ Found geometry table header at line $header_idx_local")

    # Skip the header lines and the next units line, then read data
    data_lines = []
    for i in (header_idx_local + 3):length(lines)
        line = strip(lines[i])
        # Stop when we hit a blank line or the summary section
        if isempty(line) || contains(line, "RADIUS:") || contains(line, "-----")
            break
        end
        if !isempty(line)
            push!(data_lines, line)
        end
    end

    if isempty(data_lines)
        error("No geometry data found in PE0 file. Check that data starts 2 lines after STATION header.")
    end

    println("✓ Found $(length(data_lines)) geometry stations")

    # Parse the data
    # Each line has: STATION, CHORD, PITCH(QUOTED), PITCH(LE-TE), PITCH(PRATHER), 
    #               SWEEP(Y), RAKE(Z), THICKNESS_RATIO, TWIST, MAX-THICK, CROSS-SECTION, ZHIGH, CGY, CGZ
    # Data is tab or space separated

    apc = zeros(length(data_lines), 14)
    for (i, line) in enumerate(data_lines)
        # Split by whitespace (handles both tabs and spaces)
        values = parse.(Float64, filter(!isempty, split(line)))
        
        if length(values) < 14
            error("Line $i has $(length(values)) values, expected at least 14. Content: $line")
        end
        
        apc[i, :] = values[1:14]
    end

    num_sections = size(apc, 1)

    println("✓ Successfully parsed PE0 file")
    println("  Radius: $R inches")
    println("  Number of sections: $num_sections")
    println("  Number of blades: $num_blade")
    println()

    # Extract columns from the apc data
    # STATION, CHORD, PITCH(QUOTED), PITCH(LE-TE), PITCH(PRATHER), SWEEP(Y), RAKE(Z), 
    # THICKNESS_RATIO, TWIST, MAX-THICK, CROSS-SECTION, ZHIGH, CGY, CGZ

    station = apc[:,1]          # STATION (IN)
    radius_R = station / R

    chord = apc[:,2]            # CHORD (IN)
    chord_R = chord / R

    rake = apc[:,7]             # RAKE(Z) (IN)
    axial = -rake/R

    sweep = apc[:,6]            # SWEEP(Y) (IN)
    # NOTE: SWEEP IS DEFINED WITH (MOLD) LE PARTING LINE.
    #skew_R = -sweep / R
    tangential = -sweep/R

    t_c = apc[:,8]              # THICKNESS RATIO
    twist_deg = apc[:,9]        # TWIST (DEG)

    sweep_deg = zeros(num_sections)
    CLi = zeros(num_sections)
    skew_R = zeros(num_sections)

    # tangential = zeros(num_sections)
    rake_R = zeros(num_sections)

    # Radius/R, Chord/R, Twist (deg), Rake/R, Skew/R, Sweep, t/c, CLi, Axial, Tangential
    bem = [radius_R chord_R twist_deg rake_R skew_R sweep_deg t_c CLi axial tangential]

    # Interpolate twist at 75% radius
    linterp(A, B, at) = interpolate((A,), B, Gridded(Linear()))[at]

    diameter = 2 * R * scale

    beta3_4 = 0.0

    try
        itp = interpolate((radius_R,), twist_deg, Gridded(Linear()))
        beta3_4 = itp(0.75)
    catch
        idx_75 = argmin(abs.(radius_R .- 0.75))
        beta3_4 = twist_deg[idx_75]
        println("⚠ Using closest section to 75% radius for Beta 3/4 twist (index $idx_75)")
    end

    feather = 0.0
    pre_cone = 0.0
    center = [0.0, 0.0, 0.0]
    normal = [-1.0, 0.0, 0.0]

    function writemat(io::IO, a::Matrix{<:AbstractFloat})
        lastc = last(axes(a, 2))
        for i = axes(a, 1)
            for j = axes(a, 2)
                @printf io "%.8f" a[i,j]
                j == lastc ? print(io, '\n') : print(io, ", ")
            end
        end
    end

    open(outfile, "w") do io
        println(io, "...BEM Propeller...")
        @printf io "Num_Sections: %i\n" num_sections
        @printf io "Num_Blade: %i\n" num_blade
        @printf io "Diameter: %.8f\n" diameter
        @printf io "Beta 3/4 (deg): %.8f\n" beta3_4
        @printf io "Feather (deg): %.8f\n" feather
        @printf io "Pre_Cone (deg): %.8f\n" pre_cone
        @printf io "Center: %.8f, %.8f, %.8f\n" center[1] center[2] center[3]
        @printf io "Normal: %.8f, %.8f, %.8f\n" normal[1] normal[2] normal[3]
        println(io)

        println(io, "Radius/R, Chord/R, Twist (deg), Rake/R, Skew/R, Sweep, t/c, CLi, Axial, Tangential")
        writemat(io, bem)
    end

    println("✓ Successfully converted PE0 to BEM format!")
    println("✓ Output file: $outfile")
    println()
    println("After importing BEM file into OpenVSP, set the following properties:")
    println("  Construction X/C: 0.000")
    println("  Feather Axis: 0.000")
    println()
    println("Summary:")
    @printf "  Diameter: %.2f inches (scaled to %.4f with factor %.5f)\n" (2*R) diameter scale
    println("  Number of sections: $num_sections")
    println("  Number of blades: $num_blade")
    @printf "  Beta 3/4: %.2f degrees\n" beta3_4
end


main()
