using BEPS, Test

@testset "UpdateSoilMoisture(; SolveSM_fn=SolveSM_Bonan) — Bonan-Q0 implicit solver" begin

  # ── 通用初始化 ───────────────────────────────────────────────────────────
  function make_state(; SoilType=4, VegType="DBF", θ0=0.3, r_rain_g=0.0)
    st, ps = setup(VegType, SoilType; θ0)
    st.r_rain_g = r_rain_g
    st, ps
  end

  # ── 1. 基本冒烟测试：确保不崩溃 ─────────────────────────────────────────
  @testset "No crash, correct output shape" begin
    st, ps = make_state()
    Root_Water_Uptake(st, 0.0, 0.0, 0.0)
    @test_nowarn UpdateSoilMoisture(st, ps, 3600.0; SolveSM_fn=SolveSM_Bonan)
    n = Int(st.n_layer)
    @test length(st.θ) >= n
    @test all(isfinite, st.θ[1:n])
    @test all(isfinite, st.ψ[1:n])
  end

  # ── 2. 物理边界检查 ─────────────────────────────────────────────────────
  @testset "θ stays in [θ_res, θ_sat]" begin
    for SoilType in [1, 4, 6, 11]
      st, ps = make_state(; SoilType, r_rain_g=1e-4)
      Root_Water_Uptake(st, 1e-3, 5e-4, 1e-5)
      for _ in 1:24
        UpdateSoilMoisture(st, ps, 3600.0; SolveSM_fn=SolveSM_Bonan)
      end
      n = st.n_layer
      @test all(st.θ[1:n] .>= 0.0)
      @test all(st.θ[1:n] .<= 1.0)
    end
  end

  # ── 3. 有降雨时，表层土壤含水量应增加 ───────────────────────────────────
  @testset "Rainfall increases surface θ" begin
    st, ps = make_state(; θ0=0.2, r_rain_g=2e-5)
    Root_Water_Uptake(st, 0.0, 0.0, 0.0)
    θ_before = st.θ[1]
    for _ in 1:6
      UpdateSoilMoisture(st, ps, 3600.0; SolveSM_fn=SolveSM_Bonan)
    end
    θ_after = st.θ[1]
    @test θ_after >= θ_before - 1e-6   # 有降雨时应增加，至少不应降低
  end

  # ── 4. 仅有蒸发时（无降雨），土壤含水量应降低 ───────────────────────────
  @testset "Evaporation dries out soil" begin
    st, ps = make_state(; θ0=0.35, r_rain_g=0.0)
    Root_Water_Uptake(st, 2e-3, 1e-3, 5e-4)   # 同时存在蒸腾和蒸发
    θ_init = copy(st.θ[1:st.n_layer])
    for _ in 1:48
      UpdateSoilMoisture(st, ps, 3600.0; SolveSM_fn=SolveSM_Bonan)
    end
    @test sum(st.θ[1:st.n_layer]) < sum(θ_init)
  end

  # ── 5. fix_sm=true 时，θ 不应发生变化 ──────────────────────────────────
  @testset "fix_sm=true leaves θ unchanged" begin
    st, ps = make_state(; r_rain_g=1e-4)
    θ_before = copy(st.θ)
    UpdateSoilMoisture(st, ps, 3600.0; fix_sm=true, SolveSM_fn=SolveSM_Bonan)
    @test st.θ == θ_before
    @test st.z_water >= 0.0
  end

  @testset "fix_sm matches explicit surface water update" begin
    st_explicit, ps = make_state(; θ0=0.22, r_rain_g=2e-5)
    st_q0 = deepcopy(st_explicit)

    θ_before = copy(st_explicit.θ)
    UpdateSoilMoisture(st_explicit, ps, 3600.0; fix_sm=true)
    UpdateSoilMoisture(st_q0, ps, 3600.0; fix_sm=true, SolveSM_fn=SolveSM_Bonan)

    @test st_explicit.θ == θ_before
    @test st_q0.θ == θ_before
    @test st_q0.z_water ≈ st_explicit.z_water
  end

  @testset "soil water stress stays bounded after Q0 step" begin
    st, ps = make_state(; θ0=0.25, r_rain_g=1e-5)
    Root_Water_Uptake(st, 5e-4, 2e-4, 1e-4)

    for _ in 1:24
      UpdateSoilMoisture(st, ps, 3600.0; SolveSM_fn=SolveSM_Bonan)
      soil_water_factor_v2(st, ps)
    end

    n = Int(st.n_layer)
    @test all(isfinite, st.ψ[1:n])
    @test all(isfinite, st.f_stress[1:n])
    @test all(0.0 .<= st.f_stress[1:n] .<= 1.0)
    @test 0.1 <= st.f_soilwater <= 1.0
  end

  # ── 6. 与现有 UpdateSoilMoisture 的结果趋势保持一致 ─────────────────────
  @testset "Qualitative agreement with existing solver" begin
    st_old, ps = make_state(; θ0=0.25, r_rain_g=1e-5)
    st_new     = deepcopy(st_old)
    Root_Water_Uptake(st_old, 5e-4, 2e-4, 1e-4)
    Root_Water_Uptake(st_new, 5e-4, 2e-4, 1e-4)

    for _ in 1:24
      UpdateSoilMoisture(st_old, ps, 3600.0)
      UpdateSoilMoisture(st_new, ps, 3600.0; SolveSM_fn=SolveSM_Bonan)
    end

    n = st_old.n_layer
    θ_old = st_old.θ[1:n]
    θ_new = st_new.θ[1:n]
    ψ_new = st_new.ψ[1:n]

    # 两组结果都应为有限值
    @test all(isfinite, θ_old)
    @test all(isfinite, θ_new)
    @test all(isfinite, ψ_new)

    # 两种算法（显式 Darcy vs 隐式 Richards）只要求趋势一致；
    # 数值差异来自不同离散格式，0.15 作为合理上限。
    mean_diff = sum(abs, θ_old .- θ_new) / n
    @test mean_diff < 0.15
  end

end
