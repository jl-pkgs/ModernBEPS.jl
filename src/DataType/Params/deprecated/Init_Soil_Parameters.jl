"""
    Init_Soil_Parameters(p::Soil, VegType::Integer, SoilType::Integer, r_root_decay::Float64)

Initialize soil parameters

- `K_sat`         : saturated hydraulic conductivity
- `porosity`     : porosity
- `θ_vwp`        : wilt point
- `ψ_sat`        : water potential at saturate
- `κ_dry`        : thermal conductivity
"""
function Init_Soil_Parameters(soil::Soil, VegType::Integer, SoilType::Integer, r_root_decay::Float64)
  soil.n_layer = 5

  if VegType == 6 || VegType == 9 # DBF or EBF, low constaint threshold
    soil.ψ_min = 10.0 # ψ_min
    soil.alpha = 1.5
  else
    soil.ψ_min = 33.0 # ψ_min
    soil.alpha = 0.4
  end

  soil.dz[1:5] .= [0.05, 0.10, 0.20, 0.40, 1.25] # BEPS V2023
  # z = [0, 5, 15, 25, 35, 45, 55.0] ./ 100
  # z_mid = (z[1:end-1] .+ z[2:end]) ./ 2
  # dz = diff(z)
  # n = length(z)
  # p.dz[1:n-1] = dz
  soil.r_root_decay = r_root_decay
  UpdateRootFraction!(soil)

  idx = (1 <= SoilType <= 11) ? SoilType : 11
  par = SOIL_PARAMS[idx]

  soil.ρ_soil[1:5] .= [1300.0, 1500.0, 1517.0, 1517.0, 1517.0] # from flux tower.
  soil.V_SOM[1:5] .= [5, 2, 1, 1, 0.3] ./ 100 # volume fraction, 0-1; fixed 20260502

  soil.b[1:5] .= par.b
  soil.K_sat[1:5] .= par.K_sat
  soil.θ_sat[1:5] .= fill(par.θ_sat, 5)
  # soil.θ_vfc[1:5] .= fill(par.θ_vfc, 5)
  soil.θ_vwp[1:5] .= fill(par.θ_vwp, 5)
  soil.κ_dry[1:5] .= fill(par.κ_dry, 5)
  soil.ψ_sat[1:5] .= par.ψ_sat
  return soil
end


function Init_Soil_var(soil::AbstractSoil, state::Union{StateBEPS,Vector}, Ta::FT;
  VegType::Int=25, SoilType::Int=8,
  r_drainage::FT, r_root_decay::FT,
  Tsoil0::FT=Ta, θ0::FT=0.2, z_snow0::FT=0.0,
  ignored...
) where {FT<:Real}
  # r_drainage = param[27]
  # r_root_decay = param[28]
  soil.r_drainage = r_drainage
  Init_Soil_Parameters(soil, VegType, SoilType, r_root_decay)
  Init_Soil_T_θ!(soil, Tsoil0, Ta, θ0, z_snow0)
  InitState!(soil, state, Ta)
end


function Init_Soil_T_θ!(st::S, Tsoil::Float64, Tair::Float64, θ0::Float64, snowdepth::Float64) where {
  S<:Union{StateBEPS,Soil}}

  dT = clamp(Tsoil - Tair, -5.0, 5.0)
  st.z_snow = snowdepth
  # st.z_water = 0.0
  # st.r_rain_g = 0.0
  T_scale_factors = [0.4, 0.5, 1.0, 1.2, 1.4]
  θ_scale_factors = [0.8, 1.0, 1.05, 1.10, 1.15]

  for i in 1:st.n_layer
    st.Tsoil_c[i] = Tair + T_scale_factors[i] * dT
    st.Tsoil_p[i] = Tair + T_scale_factors[i] * dT
    st.θ[i] = θ_scale_factors[i] * θ0
    st.θ_prev[i] = θ_scale_factors[i] * θ0
    st.ice_ratio[i] = get_ice_ratio(st.Tsoil_c[i])
  end
end


# for (i=3;i<=8;i++)   var_o[i] = tem;
# for(i=9;i<=14;i++) var_o[i] = soil->Tsoil_p[i-9];
# for(i=21;i<=26;i++) var_o[i] = soil->θ_prev[i-21];
# for(i=27;i<=32;i++) var_o[i] = soil->ice_ratio[i-27];
function InitState!(soil::AbstractSoil, state::Vector, Ta)
  state .= 0
  for i = 1:6
    state[i+3] = Ta
    state[i+9] = soil.Tsoil_p[i]
    state[i+21] = soil.θ_prev[i]
    state[i+27] = soil.ice_ratio[i]
  end
  return nothing
end

function InitState!(soil::AbstractSoil, state::StateBEPS, Ta)
  n = 5
  state.Tsnow_c .= Ta
  # state.Ts_prev .= soil.Tsoil_p[1:5]
  state.θ_prev[1:n] .= soil.θ_prev[1:n]
  state.ice_ratio[1:n] .= soil.ice_ratio[1:n]
  return nothing
end


export Init_Soil_var
