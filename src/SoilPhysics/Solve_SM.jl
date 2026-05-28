import ModelParams: soil_moisture_Q0!, TriSolver
export UpdateSoilMoisture_Q0!


"""
    solve_SM_beps_Q0!(st, ps, inf, kstep)

Bonan-Q0 implicit Crank-Nicolson step for BEPS soil moisture.

`inf` [m/s]: infiltration from `update_surface_water!`
`kstep` [s]: external timestep

Replaces `solve_SM_beps` with one implicit tridiagonal solve per step.
"""
function solve_SM_beps_Q0!(st::StateBEPS, ps::ParamBEPS, inf::Float64, kstep::Float64)
  n = Int(st.n_layer)   # Cint → Int64 (soil_moisture_Q0! requires Int)
  FT = Float64
  # inf [m/s] → Q0 [cm/h], negative = downward infiltration
  Q0 = FT(-inf * 3.6e5)
  # ETi [m/s] → sink [cm/h]
  sink = @view(st.ETi[1:n]) .* FT(3.6e5)

  tri = st.tri
  soil_moisture_Q0!(
    st.θ, st.ψ, st.θ_prev, st.ψ_prev, st.∂θ∂ψ, st.K, st.K₊ₕ, tri,
    ps.hydraulic, sink, Q0;
    ibeg=st.ibeg, N=n,
    Δz_cm=st.Δz_cm, Δz₊ₕ_cm=st.Δz₊ₕ_cm,
    dt=FT(kstep)
  )
  # soil_moisture_Q0! writes ψ in negative cm (Richards convention).
  # Convert back to positive m (BEPS convention) for downstream code.
  @inbounds for i in 1:n
    st.ψ[i] = -st.ψ[i] / 100.0
  end
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
