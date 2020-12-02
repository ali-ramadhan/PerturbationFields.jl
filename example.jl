Nx = Ny = 192

Δx = Δy = 125.0

maxval = 1.0  # maximum absolute value of perturbations

# integer, parameter :: spectrum_type = 1
# ! 1 = white noise (kappa^0)
# ! 2 = kappa^(-5/3)
# ! 3 = kappa^(-3)

# κ = 0

# Large eddy size (m).
# Above this scale, the spectrum is
# white.  If you don't want to use
# this, make it larger than
# max(ni,nj)*dx.
large_eddy_scale = 1000000000.0

# integer, parameter :: filter_type = 2
# ! 0 = no filter
# ! 1 = remove all scales smaller than
# !     6 Delta
# ! 2 = disk:  keep all scales between
# !     L1 and L2
# ! 3 = ring:  keep only scale L1

# ! For filter_type 2 and 3 only:
# real, parameter :: L1 =  1000.0   ! smaller scale (m)
# real, parameter :: L2 =  8000.0   ! larger scale (m)

i_c = Nx/2 + 1
j_c = Ny/2 + 1

# do j=1,nj
# do i=1,ni
#   kappa(i,j)=sqrt( ( 2.0*pi*(float(i)-i_c)/(ni*dx))**2   &
#                   +( 2.0*pi*(float(j)-j_c)/(nj*dx))**2 )
# enddo
# enddo

κ = zeros(Nx, Ny)
for j in 1:Ny, i in 1:Nx
    κ[i, j] = √( (2π*(i-i_c) / (Nx*Δx))^2 + (2π*(j-j_c) / (Ny*Δy))^2 )
end

c = zeros(ComplexF64, Nx, Ny)
for j in 1:Ny, i in 1:Nx
    r1 = rand()
    r2 = rand()
    c[i, j] = exp(-2π*im*Nx*r1) * exp(-2π*im*Ny*r2)
end

# do j=1,nj
# do i=1,ni
#   val = abs(kappa(i,j)-(2.0*pi/L1))
#   if(val.lt.vmin)then
#     vmin = val
#     imin = i
#     jmin = j
#   endif
# enddo
# enddo

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

cc = circshift(c, (Nx/2, Ny/2))

KE = ifft(cc)

ccc = circshift(fft(KE), (Nx/2, Ny/2))

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
