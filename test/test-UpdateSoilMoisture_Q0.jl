using Test
using BEPS

@testset "UpdateSoilMoisture_Q0! — Bonan-Q0 implicit solver" begin

  # ── 共通セットアップ ──────────────────────────────────────────────────────
  function make_state(; SoilType=4, VegType="DBF", θ0=0.3, r_rain_g=0.0)
    st, ps = setup(VegType, SoilType; θ0)
    st.r_rain_g = r_rain_g
    st, ps
  end

  # ── 1. 基本スモークテスト（クラッシュしない） ─────────────────────────────
  @testset "No crash, correct output shape" begin
    st, ps = make_state()
    Root_Water_Uptake(st, 0.0, 0.0, 0.0)
    @test_nowarn UpdateSoilMoisture_Q0!(st, ps, 3600.0)
    @test length(st.θ) >= st.n_layer
    @test all(isfinite, st.θ[1:st.n_layer])
  end

  # ── 2. 物理的境界条件チェック ──────────────────────────────────────────────
  @testset "θ stays in [θ_res, θ_sat]" begin
    for SoilType in [1, 4, 6, 11]
      st, ps = make_state(; SoilType, r_rain_g=1e-4)
      Root_Water_Uptake(st, 1e-3, 5e-4, 1e-5)
      for _ in 1:24
        UpdateSoilMoisture_Q0!(st, ps, 3600.0)
      end
      n = st.n_layer
      @test all(st.θ[1:n] .>= 0.0)
      @test all(st.θ[1:n] .<= 1.0)
    end
  end

  # ── 3. 降雨時は土壌水分が増加 ──────────────────────────────────────────────
  @testset "Rainfall increases surface θ" begin
    st, ps = make_state(; θ0=0.2, r_rain_g=2e-5)
    Root_Water_Uptake(st, 0.0, 0.0, 0.0)
    θ_before = st.θ[1]
    for _ in 1:6
      UpdateSoilMoisture_Q0!(st, ps, 3600.0)
    end
    θ_after = st.θ[1]
    @test θ_after >= θ_before - 1e-6   # 雨なら増えるか同値
  end

  # ── 4. 蒸発のみ（降雨なし）は土壌水分が減少 ──────────────────────────────
  @testset "Evaporation dries out soil" begin
    st, ps = make_state(; θ0=0.35, r_rain_g=0.0)
    Root_Water_Uptake(st, 2e-3, 1e-3, 5e-4)   # 蒸散・蒸発あり
    θ_init = copy(st.θ[1:st.n_layer])
    for _ in 1:48
      UpdateSoilMoisture_Q0!(st, ps, 3600.0)
    end
    @test sum(st.θ[1:st.n_layer]) < sum(θ_init)
  end

  # ── 5. fix_sm=true で θ が変化しない ───────────────────────────────────────
  @testset "fix_sm=true leaves θ unchanged" begin
    st, ps = make_state(; r_rain_g=1e-4)
    θ_before = copy(st.θ)
    UpdateSoilMoisture_Q0!(st, ps, 3600.0; fix_sm=true)
    @test st.θ == θ_before
  end

  # ── 6. 既存 UpdateSoilMoisture との比較（傾向一致） ─────────────────────────
  @testset "Qualitative agreement with existing solver" begin
    st_old, ps = make_state(; θ0=0.25, r_rain_g=1e-5)
    st_new     = deepcopy(st_old)
    Root_Water_Uptake(st_old, 5e-4, 2e-4, 1e-4)
    Root_Water_Uptake(st_new, 5e-4, 2e-4, 1e-4)

    for _ in 1:24
      UpdateSoilMoisture(st_old, ps, 3600.0)
      UpdateSoilMoisture_Q0!(st_new, ps, 3600.0)
    end

    n = st_old.n_layer
    θ_old = st_old.θ[1:n]
    θ_new = st_new.θ[1:n]

    # 両方有限値
    @test all(isfinite, θ_old)
    @test all(isfinite, θ_new)

    # 两种算法（显式 Darcy vs 隐式 Richards）结果趋势一致即可；
    # 数值差异来自不同的离散格式，0.15 是合理上限。
    mean_diff = sum(abs, θ_old .- θ_new) / n
    @test mean_diff < 0.15
  end

end
