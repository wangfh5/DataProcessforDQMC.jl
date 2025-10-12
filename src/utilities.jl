#=
Utility Functions (utilities.jl)

This file contains common utility functions shared across different modules.
=#

export find_closest_k_point, find_coordinate_index

# ---------------------------------------------------------------------------- #
#                        Coordinate matching utilities                         #
# ---------------------------------------------------------------------------- #

"""
    find_closest_k_point(k_points, target_k, tolerance)

Find the closest k-point in a set of k-points to a target k-point.

# Arguments
- `k_points::Vector{<:Tuple{<:Real,<:Real}}`: Vector of k-point coordinates [(kx1, ky1), (kx2, ky2), ...]
- `target_k::Tuple{<:Real,<:Real}`: Target k-point (kx, ky)
- `tolerance::Float64`: Tolerance for matching k-points

# Returns
- `Tuple{Real, Real}`: Closest k-point coordinates (kx, ky)
- `Bool`: Whether it's an exact match (distance <= tolerance)
- `Float64`: Distance to the target k-point

# Example
```julia
k_points = [(0.0, 0.0), (0.5, 0.0), (0.5, 0.5)]
target = (0.49, 0.51)
closest_k, exact_match, distance = find_closest_k_point(k_points, target, 0.1)
# Returns: ((0.5, 0.5), true, 0.014...)
```
"""
function find_closest_k_point(k_points::Vector{<:Tuple{<:Real,<:Real}}, target_k::Tuple{<:Real,<:Real}, tolerance::Float64)
    # Calculate distance from each k-point to target
    distances = [(kx, ky, sqrt((kx - target_k[1])^2 + (ky - target_k[2])^2)) for (kx, ky) in k_points]
    
    # Sort by distance
    sort!(distances, by = x -> x[3])
    
    # Check if closest point is within tolerance
    closest_k = (distances[1][1], distances[1][2])
    exact_match = distances[1][3] <= tolerance
    
    # Return closest k-point, exact match flag, and distance
    return closest_k, exact_match, distances[1][3]
end

"""
    find_coordinate_index(coords, target_coord, tolerance)

Find the index of a coordinate in a coordinate matrix.

This function handles both r-space (integer) and k-space (float) coordinates:
- For r-space (integer tuples): exact integer matching
- For k-space (float tuples): nearest neighbor matching with tolerance

# Arguments
- `coords::Matrix`: Matrix of coordinates (Nx2, first two columns of bin file)
- `target_coord::Union{Tuple{Int, Int}, Tuple{Float64, Float64}}`: Target coordinate
- `tolerance::Float64=1e-6`: Tolerance for k-space matching (ignored for r-space)

# Returns
- `Union{Int, Nothing}`: Index of the coordinate (1-based), or Nothing if not found
- `Bool`: Whether the match is exact (k-space) or found (r-space)

# Examples
```julia
# R-space example
coords = [1.0 1.0; 1.0 2.0; 2.0 1.0]
idx, exact = find_coordinate_index(coords, (1, 2), 1e-6)
# Returns: (2, true)

# K-space example
coords = [0.0 0.0; 0.5 0.0; 0.5 0.5]
idx, exact = find_coordinate_index(coords, (0.49, 0.51), 0.1)
# Returns: (3, true)  # distance = 0.014 < 0.1
```
"""
function find_coordinate_index(coords::Matrix, target_coord::Union{Tuple{Int, Int}, Tuple{Float64, Float64}}, tolerance::Float64=1e-6)
    # Check if coordinates are integers (r-space) or floats (k-space)
    if target_coord isa Tuple{Int, Int}
        # R-space: exact integer matching
        for i in 1:size(coords, 1)
            if Int(coords[i, 1]) == target_coord[1] && Int(coords[i, 2]) == target_coord[2]
                return i, true
            end
        end
        return nothing, false
    else
        # K-space: use find_closest_k_point for consistent behavior
        k_points = [(coords[i, 1], coords[i, 2]) for i in 1:size(coords, 1)]
        closest_k, exact_match, dist = find_closest_k_point(k_points, target_coord, tolerance)
        
        # Find the index of the closest k-point
        for i in 1:size(coords, 1)
            if coords[i, 1] == closest_k[1] && coords[i, 2] == closest_k[2]
                return i, exact_match
            end
        end
        
        return nothing, false
    end
end

