maxval = 1.0  # maximum absolute value of perturbations

# Large eddy size (m).
# Above this scale, the spectrum is
# white.  If you don't want to use
# this, make it larger than
# max(ni,nj)*dx.
large_eddy_scale = 1000000000.0

function ring_filter(Nx, Ny, Δx, Δy)
    # center point indices
    i_c = Nx/2 + 1
    j_c = Ny/2 + 1

    # Compute wavenumbers
    κ = zeros(Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        κ[i, j] = √( (2π*(i-i_c) / (Nx*Δx))^2 + (2π*(j-j_c) / (Ny*Δy))^2 )
    end

    # Generate spectrum
    A = 1
    c = zeros(ComplexF64, Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        r1 = rand()
        r2 = rand()
        c[i, j] = A * exp(-2π*im*Nx*r1) * exp(-2π*im*Ny*r2)
    end

    L₁ = 1000

    i_min = j_min = 0
    v_min = Inf
    for j in 1:Ny, i in 1:Nx
        val = abs(κ[i, j] - 2π/L₁)
        if val < v_min
            v_min = val
            i_min = i
            j_min = j
        end
    end

    threshold = max(
        abs(κ[i_min, j_min+1] - 2π/L₁),
        abs(κ[i_min, j_min-1] - 2π/L₁),
        abs(κ[i_min+1, j_min] - 2π/L₁),
        abs(κ[i_min-1, j_min] - 2π/L₁)
    )

    for j in 1:Ny, i in 1:Nx
        if abs(κ[i, j] - 2π/L₁) > threshold
            c[i, j] = 0
        end
    end

    spectrum_2d = deepcopy(c)

    c = circshift(c, (Nx/2, Ny/2))
    KE = ifft(c)
    c = circshift(fft(KE), (Nx/2, Ny/2))

    Nk = min(Nx, Ny)
    sp = zeros(Nk)
    np = zeros(Nk)

    for j in 1:Ny, i in 1:Nx
        if Nx > Ny
            ik = √( ((i-i_c)*Ny/Nx)^2 + (j-j_c)^2 )
            ik = round(Int, ik)
        elseif Ny > Nx
            ik = √( (i-i_c)^2 + ((j-j_c)*Nx/Ny)^2 )
            ik = round(Int, ik)
        else
            ik = √( (i-i_c)^2 + (j-j_c)^2 )
            ik = round(Int, ik)
        end

        if ik > 0 && ik < Nk/2
            sp[ik] += abs(ccc[i, j])^2
            np[ik] += 1
        end
    end

    spectrum_1d = deepcopy(sp)

    return abs2.(KE), spectrum_2d, spectrum_1d
end

Nx = Ny = 192
Δx = Δy = 125

x = range(0, Nx*Δx, length=Nx)
y = range(0, Ny*Δy, length=Ny)

KE, s2, s1 = ring_filter(Nx, Ny, Δx, Δy)

l = @layout [
    a{0.7w} [b{0.5h}
             c{0.5h}]
]

pert_plot = contourf(x/1000, y/1000, KE, color=:dense, linewidth=0, xticks=0:4:24, xlabel="x (km)", yticks=0:4:24, ylabel="y (km)", title="Perturbation field", colorbar=nothing)
s2_plot = contourf(abs2.(s2), color=:dense, linewidth=0, xticks=[], yticks=[], colorbar=nothing, framestyle=:box, title="2D spectrum", xlabel="kx", ylabel="ky")
s1_plot = plot(log10.(s1), linewidth=2, label="", xlims=(0, 96), xticks=0:32:96, xlabel="k", ylabel="log E(k)", title="1D spectrum")
plot(pert_plot, s2_plot, s1_plot, layout=l)
