import ModelParams: soil_moisture_Q0!, TriSolver


"""
    SolveSM_Bonan(st, ps, inf, kstep)

Bonan-Q0 implicit Crank-Nicolson step for BEPS soil moisture.

`inf` [cm/h]: infiltration from `update_surface_water!`
`kstep` [s]: external timestep

Replaces `solve_SM_beps` with one implicit tridiagonal solve per step.
"""
function SolveSM_Bonan(st::StateBEPS, ps::ParamBEPS, inf::Float64, kstep::Float64)
  (; θ, ψ, θ_prev, ψ_prev, ∂θ∂ψ, K, K₊ₕ, tri,
    ibeg, Δz_cm, Δz₊ₕ_cm) = st
  (; hydraulic) = ps

  N = Int(st.n_layer)   # Cint → Int64 (soil_moisture_Q0! requires Int)
  FT = Float64

  Q0 = FT(-inf)                          # negative = downward infiltration
  sink = @view(st.ETi[1:N]) .* FT(3.6e5) # ETi [m/s] → sink [cm/h]

  soil_moisture_Q0!(
    θ, ψ, θ_prev, ψ_prev, ∂θ∂ψ, K, K₊ₕ, tri,
    hydraulic, sink, Q0;
    ibeg, N, Δz_cm, Δz₊ₕ_cm, dt=FT(kstep)
  )
  # Keep the raw pressure head from ModelParams.
  # Under saturated conditions ψ may become positive, so do not force a sign here.
end
