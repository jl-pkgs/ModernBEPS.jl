# Function to compute soil water stress factor
function soil_water_factor_v2(st::S, ps::P) where {S<:Union{StateBEPS,Soil},P<:Union{ParamBEPS,Soil}}
  (; ψ_min, alpha) = ps
  (; θ_sat, b) = get_hydraulic(ps)
  ψ_sat_m = _get_ψ_sat_m(ps)   # positive [m], handles both ParamBEPS and Soil

  θ = st.θ
  n = st.n_layer

  t1 = -0.02
  t2 = 2.0

  if st.ψ[1] <= 0.000001
    for i in 1:n
      st.ψ[i] = cal_ψ(θ[i], θ_sat[i], ψ_sat_m[i], b[i])
    end
  end

  for i in 1:n
    # psi_sr in m H2O! He 2017 JGR-B, Eq. 4
    st.f_stress[i] = st.ψ[i] > ψ_min ? 1.0 / (1 + ((st.ψ[i] - ψ_min) / ψ_min)^alpha) : 1.0
    st.f_temp[i] = st.Tsoil_p[i] > 0.0 ? 1.0 - exp(t1 * st.Tsoil_p[i]^t2) : 0

    st.f_stress[i] *= st.f_temp[i]
    st.w_root[i] = FW_VERSION == 1 ? st.f_root[i] * st.f_stress[i] : st.f_root[i]
  end

  w_root_sum = sum(st.w_root) # 每层的土壤水分限制因子

  if w_root_sum < 0.000001
    st.f_soilwater = 0.1
  else
    f_stress_sum = 0.0
    for i in 1:n
      st.w_norm[i] = st.w_root[i] / w_root_sum
      f_stress_sum += st.f_stress[i] * st.w_norm[i]
    end
    st.f_soilwater = max(0.1, f_stress_sum)
  end
end
soil_water_factor_v2(p::Soil) = soil_water_factor_v2(p, p)
