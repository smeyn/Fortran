struct onerun 
    splinetype
    kx
    ky
    smoothing_factor
    r_square
    error_flag
    nx
    ny
    knots_x
    knots_y
    spline_coefficients
    eval_results
    eval_x
    eval_y
end

function get(f, offset, parsetype, endMatch)
    
    m = match(endMatch, f, offset)
    if isnothing(m)
        return parse(parsetype,"0"), offset
    end
    @info offset, m.match
    s = strip(m.captures[1])
    return parse(parsetype, s), m.offset #+length(m.match)
end
"""same as get but parse n repetitions of type"""
function get_n(f, offset, parsetype, endMatch)
    m = match(endMatch, f, offset)
    @debug offset, m.match
    s = strip(m.captures[1])
    @debug offset, s
    parts = split(s, " ")
    result = [parse(parsetype, strip(item)) for item in parts if strip(item) != ""]
    return result, m.offset #+length(m.match)
end

function read_section(fortran, offset)
    """ assumes we are already past '0smoothing spline of degrees '"""
    @debug "Reading section" , offset
    m = match(r"^\d(.* spline) of degrees\s+(\d+)\s+(\d+)"m, fortran, offset)
    @debug m.offset, m.match
    splinetype = m.captures[1]
    kx = parse(Int, m.captures[2])
    ky =  parse(Int, m.captures[3])
    offset = m.offset #+ length(m.match) 

    smoothing_factor, offset = get(fortran, offset, Float64, r"smoothing factor s=\s+(\d+\S)")
    r_square, offset = get(fortran, offset, Float64, r"sum squared residuals =\s+(\S+)")
    error_flag, offset = get(fortran, offset, Float64, r"error flag=\s+(\S+)")
    nx, offset = get(fortran, offset, Int, r"total number of knots in the x-direction =\s+(\d+)")
    tx, offset = get_n(fortran, offset, Float64, r"position of the knots\s+\R(.*)\R")
    ny, offset = get(fortran, offset, Int, r"number of knots in the y-direction =\s+(\d+)")
    ty, offset = get_n(fortran, offset, Float64,  r"position of the knots\s+\R(.*)\R")
    @info offset
    coefficients, offset = get_n(fortran, offset, Float64, r"coefficients\s+\R([\r\n .0-9+-]*)\dspl"ms)
    x_coords, offset = get_n(fortran, offset, Float64, r"evaluation on a given grid\s+\R.x\s+([\r\n .0-9+-]*)"ms)
    
    matvals, offset = get_n(fortran, offset, Float64, r"y\s*\R([\r\n .0-9+-]*)\R"ms)

    #matvals, offset = get_n(fortran, offset, Float64, r"(.*)0 debug: b-spline coefficients")
    ycoords = matvals[1,:]
    eval_results = matvals[2:end,:]
    return onerun(
        splinetype,
        kx, ky, smoothing_factor, r_square, error_flag, 
        nx, ny, tx, ty,coefficients, eval_results,x_coords, ycoords
    ), offset
end

"""read the fortran results file"""
function readResults(path)::Vector
    @info "reading Fortran results $path"
    list = Vector()
    fortran = read(path, String)
    index = 1
          
    while index < length(fortran)
        @info "skipping to start of section"
        m = match(r"\w+\s+spline of degrees "s, fortran, index)
        @info m
        if isnothing(m)
            break
        end
        index = m.offset -50 

        if index >= length(fortran)
            @info "eof reached"
            @info "last read" m.match
            break
        end
        section, index = read_section(fortran, index)
        if !isnothing(section)
            push!(list, section)
        else
            break
        end
        # now skip to next section
    end
    return list
end




