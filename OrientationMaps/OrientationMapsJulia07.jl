# This Julia code impelments the long-range interation model introduced by Wolf (PRL 2005).
 # COPYRIGHT: Engelken, R. (2025) based on a matlab implementation by Manuel Schottdorf. 
 # GPLv3
using FFTW, Random, Statistics, PyPlot, LinearAlgebra

FFTW.set_num_threads(10)

function lrim(T, 𝒢, σ)
    # ─── Initialization ──────────────────────────────────────────────
    Random.seed!(1)   # fixed seed for reproducibility
    N   = 512         # grid size
    nₕ  = 22          # spatial scaling factor
    Δt  = 1/2^10       # initial time step
    t   = 0.0         # simulation time
    plot_interval = 30
    tol = 2e-1        # error tolerance
    σ₀ = 1e-2         # std dev of initial noise

    # Build the wavevector grid (KK) and meshgrid
    KK  = [(2*(k-1)/N - 1)*N/(2*nₕ) for k in 1:N]
    KKx = [xi for xi in KK, _ in KK]
    KKy = [yi for _ in KK, yi in KK]

    # Gaussian kernel + FFT
    s = σ*N/nₕ
    h = [exp(-(((x - N/2)^2 + (y - N/2)^2)/(2*s^2))) for x in 1:N, y in 1:N]
    h ./= sum(h)
    𝒢̂ = fft(fftshift(h))

    # Orientation preference map (OPM)
    z = σ₀ * randn(ComplexF64, N, N)
    R𝒦 = Matrix{Float64}(undef, N, N)  # shift kernel

    # ── Preallocate intermediate arrays to avoid repeated allocation ──
    absz²,z²,N₁,N₂,nonlocal_term,local_term,Z,Zs,z0,z_full,z_half = ntuple(_ -> similar(z), 11)

    # ─── Semi-implicit update step ───────────────────────────────────
    function step!(z, Δt_local)
        @inbounds for k in 1:N, l in 1:N
            R𝒦[k, l] = 1 / (1 - Δt_local*(0.01 - (1 - (KK[k]^2 + KK[l]^2))^2))
        end

        # absz² = |z|^2
        absz² .= abs2.(z)
        # N₁ = ifft( 𝒢̂ * fft(absz²) ) .* z
        N₁ .= ifft(𝒢̂ .* fft(absz²)) .* z
        # z² = z^2
        z² .= z .^ 2
        # N₂ = ifft( 𝒢̂ * fft(z²) ) .* conj(z) * 0.5
        N₂ .= ifft(𝒢̂ .* fft(z²)) .* conj.(z) .* 0.5

        nonlocal_term .= (𝒢 - 2) .* (N₁ .+ N₂) .* Δt_local
        local_term   .= (1 - 𝒢) .* absz² .* z .* Δt_local

        # Regularized shift in Fourier space
        Z .= fft(z)
        Zs .= fftshift(Z)
        Zs .*= R𝒦

        z .= ifft(ifftshift(Zs)) .+ local_term .+ nonlocal_term
    end

    # ─── Main adaptive time-stepping loop ─────────────────────────────
    stepcount = -1
    fig, axs = plt.subplots(1, 2; figsize=(12,5), constrained_layout=true)

    while t < T
        stepcount += 1

        # Step doubling for local error estimate
        z0 .= z
        z_full .= z0;  step!(z_full, Δt)
        z_half .= z0;  step!(z_half, Δt/2);  step!(z_half, Δt/2)
        err = norm(z_full .- z_half) / (norm(z_half) + 1e-10)

        if err <= tol+1e-13 
            z .= z_half
            t += Δt
            Δt *= min(1.2, (tol/err)^0.5)
        else
            Δt *= max(0.5, (tol/err)^0.5)
        end

        # Periodic plotting
        if stepcount % plot_interval == 0 || ((t<300) && stepcount % div(plot_interval,10) == 0 )
            ax_phase, ax_spec = axs
            ax_phase.clear()
            ax_spec.clear()

            c_phase = ax_phase.imshow(angle.(z), aspect="equal", cmap="hsv")
            ax_phase.set_title("OPM Phase")
            ax_phase.axis("off")
            (stepcount == 0) && fig.colorbar(c_phase, ax=ax_phase)

            spec = abs.(fftshift(fft(z .- mean(z)))).^2
            c_spec = ax_spec.pcolormesh(KKx, KKy, spec, cmap="gray_r", shading="auto")
            ax_spec.set_aspect("equal", "box")
            ax_spec.set_xlim(-1.5, 1.5)
            ax_spec.set_ylim(-1.5, 1.5)
            ax_spec.set_title("OPM Spectrum")
            (stepcount == 0) && fig.colorbar(c_spec, ax=ax_spec)

            fig.suptitle("t = $(round(t,digits=1))")
            plt.show()
        end

        println("t = $t, step = $stepcount, Δt = $Δt, err = $err, ⟨|z|⟩ = $(mean(abs.(z)))")
    end
end

# Example usage:
# lrim(1e5, 0.98, 1.7)

