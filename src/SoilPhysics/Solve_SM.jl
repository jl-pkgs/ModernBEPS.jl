import ModelParams: soil_moisture_Q0!, TriSolver
export UpdateSoilMoisture_Q0!


"""
    solve_SM_beps_Q0!(st, ps, inf, kstep)

Bonan-Q0 implicit Crank-Nicolson step for BEPS soil moisture.

`inf` [cm/h]: infiltration from `update_surface_water!`
`kstep` [s]: external timestep

Replaces `solve_SM_beps` with one implicit tridiagonal solve per step.
"""
function solve_SM_beps_Q0!(st::StateBEPS, ps::ParamBEPS, inf::Float64, kstep::Float64)
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


"""
    UpdateSoilMoisture_Q0!(st, ps, kstep; fix_sm=false)

Drop-in replacement for `UpdateSoilMoisture` using the Bonan-Q0 implicit solver.
Surface infiltration calculation is identical; only the sub-surface transport
method changes (one implicit step vs. many adaptive explicit steps).
"""
function UpdateSoilMoisture_Q0!(st::StateBEPS, ps::ParamBEPS, kstep::Float64; fix_sm::Bool=false)
  n = Int(st.n_layer)
  (; θ, θ_prev, ice_ratio) = st

  θ_prev .= θ
  inf = update_surface_water!(st, ps, kstep)
  fix_sm && return

  solve_SM_beps_Q0!(st, ps, inf, kstep)

  for i in 1:n
    θ[i] > 0 && (ice_ratio[i] = min(1.0, ice_ratio[i] * θ_prev[i] / θ[i]))
  end
end
