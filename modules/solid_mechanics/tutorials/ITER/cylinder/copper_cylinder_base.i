### Base input file for ITER copper cylinder simulations
### Contains common mesh, materials, and material properties
### Include this file in specialized simulations

# ITER Geometry
cable_length = '${units 0.45 m -> mm}'
cable_radius = '${units 10 mm}'

# Temperature Parameters
starting_temperature = '${units 4.5 K}'
ending_temperature = '${units 300 K}'
stress_free_temperature = '${starting_temperature}'

# Simulation Parameters
simulation_time = '${units 10 s}'

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [circle]
    type = ConcentricCircleMeshGenerator
    num_sectors = 4
    has_outer_square = false
    preserve_volumes = true
    radii = '${cable_radius}'
    rings = '2'
    smoothing_max_it = 3
  []
  [name_blocks]
    type = RenameBlockGenerator
    input = circle
    old_block = '1'
    new_block = 'cable'
  []
  [stabilizer]
    type = AdvancedExtruderGenerator
    input = name_blocks
    direction = '0 0 1'
    heights = '${cable_length}'
    num_layers = '24'  # Refined to match radial element size (~6.25 mm vs. 5-7.85 mm radial)
    biases = '1'
    bottom_boundary = 'axial_start'
    top_boundary = 'axial_end'
  []
  [pin_bottom_0_deg]
    type = ExtraNodesetGenerator
    input = stabilizer
    new_boundary = pin_bottom_0_deg
    coord = '${cable_radius} 0 0'
    use_closest_node = true
  []
  [pin_bottom_90_deg]
    type = ExtraNodesetGenerator
    input = pin_bottom_0_deg
    new_boundary = pin_bottom_90_deg
    coord = '0 ${cable_radius} 0'
    use_closest_node = true
  []
  [pin_bottom_180_deg]
    type = ExtraNodesetGenerator
    input = pin_bottom_90_deg
    new_boundary = pin_bottom_180_deg
    coord = '${fparse -cable_radius} 0 0'
    use_closest_node = true
  []
  [pin_bottom_center]
    type = ExtraNodesetGenerator
    input = pin_bottom_180_deg
    new_boundary = pin_bottom_center
    coord = '0 0 0'
    use_closest_node = true
  []
[]

[Physics/SolidMechanics]
  [QuasiStatic]
    [all]
      strain = FINITE
      add_variables = true
      # decomposition_method = EIGENSOLUTION
      eigenstrain_names = 'thermal_expansion'
      generate_output = 'vonmises_stress strain_xx strain_yy strain_zz strain_xy strain_xz strain_yz stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz'
      material_output_family = LAGRANGE
      material_output_order = FIRST
      temperature = T
      use_automatic_differentiation = true
    []
  []
[]

[Materials]
  [stress]
    type = ADComputeFiniteStrainElasticStress
  []

  [youngs_modulus]
    type = ADParsedMaterial
    property_name = youngs_modulus
    coupled_variables = T
    # NIST correlation: E(T) = 1e3 * (137 - 1.27e-04 * T²) [MPa]
    # Source: Page 6-1, NIST Monograph 177
    expression = '1e3*(137-1.27e-04*T^2)'
  []

  [poissons_ratio]
    type = ADParsedMaterial
    property_name = poissons_ratio
    coupled_variables = T
    # NIST correlation: ν(T) = 0.339 + 7.03e-08 * T²
    # Source: Page 6-23, NIST Monograph 177
    expression = '0.339 + 7.03e-08*T^2'
  []

  [elasticity]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
  []

  [expansion]
    type = ADComputeInstantaneousThermalExpansionFunctionEigenstrain
    temperature = T
    thermal_expansion_function = thermal_expansion_func
    stress_free_temperature = ${stress_free_temperature}
    eigenstrain_name = thermal_expansion
    outputs = 'exodus'
  []
[]

[FunctorMaterials]
  [thermal_expansion_coeff_functor]
    # Expression duplicated from thermal_expansion_func for postprocessing with spatially-varying T
    type = ADParsedFunctorMaterial
    expression = '1e-6 * (10^(c0 + c1*log10(T) + c2*log10(T)^2 + c3*log10(T)^3 + c4*log10(T)^4 + c5*log10(T)^5 + c6*log10(T)^6))'
    functor_names = 'T -17.9081289 67.131914 -118.809316 109.9845997 -53.8696089 13.30247491 -1.30843441'
    functor_symbols = 'T c0 c1 c2 c3 c4 c5 c6'
    property_name = thermal_expansion_coeff_functor
    output_properties = thermal_expansion_coeff_functor
    outputs = 'exodus'
  []
[]

[AuxVariables]
  [T]
    family = LAGRANGE
    order = FIRST
    initial_condition = '${starting_temperature}'
  []
  [thermal_strain_xx_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [thermal_strain_yy_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [thermal_strain_zz_aux]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [temperature_ramp]
    type = FunctionAux
    variable = T
    function = 'temperature_func'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [thermal_strain_xx_extract]
    type = ADRankTwoAux
    variable = thermal_strain_xx_aux
    rank_two_tensor = thermal_expansion
    index_i = 0
    index_j = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [thermal_strain_yy_extract]
    type = ADRankTwoAux
    variable = thermal_strain_yy_aux
    rank_two_tensor = thermal_expansion
    index_i = 1
    index_j = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [thermal_strain_zz_extract]
    type = ADRankTwoAux
    variable = thermal_strain_zz_aux
    rank_two_tensor = thermal_expansion
    index_i = 2
    index_j = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Functions]
  [temperature_func]
    type = PiecewiseLinear
    x = '0 ${simulation_time}'
    y = '${starting_temperature} ${ending_temperature}'
  []

  [thermal_expansion_func]
    # 't' parameter represents temperature when used by thermal expansion material
    # NIST OFHC Copper thermal expansion coefficient correlation (4-300 K)
    # log₁₀(α [10⁻⁶/K]) = Σ cᵢ × log₁₀(T)ⁱ
    # Source: https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm
    type = ParsedFunction
    expression = '1e-6 * (10^(c0 + c1*log10(t) + c2*log10(t)^2 + c3*log10(t)^3 + c4*log10(t)^4 + c5*log10(t)^5 + c6*log10(t)^6))'
    symbol_names = 'c0 c1 c2 c3 c4 c5 c6'
    symbol_values = '-17.9081289 67.131914 -118.809316 109.9845997 -53.8696089 13.30247491 -1.30843441'
  []
[]


[Postprocessors]
  # Simulation parameters (for verification script)
  # CSV only - no need for terminal output
  [stress_free_temperature]
    type = ConstantPostprocessor
    value = ${stress_free_temperature}
    execute_on = 'INITIAL'
    outputs = 'csv'
  []
  [starting_temperature]
    type = ConstantPostprocessor
    value = ${starting_temperature}
    execute_on = 'INITIAL'
    outputs = 'csv'
  []
  [ending_temperature]
    type = ConstantPostprocessor
    value = ${ending_temperature}
    execute_on = 'INITIAL'
    outputs = 'csv'
  []
  [simulation_time]
    type = ConstantPostprocessor
    value = ${simulation_time}
    execute_on = 'INITIAL'
    outputs = 'csv'
  []

  # Temperature and material properties
  # Keep temperature_average on console for monitoring
  [temperature_average]
    type = ElementAverageValue
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [thermal_expansion_coeff_pp]
    type = ElementAverageFunctorPostprocessor
    functor = thermal_expansion_coeff_functor
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [thermal_strain_xx_average]
    type = ElementAverageValue
    variable = thermal_strain_xx_aux
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [thermal_strain_yy_average]
    type = ElementAverageValue
    variable = thermal_strain_yy_aux
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [thermal_strain_zz_average]
    type = ElementAverageValue
    variable = thermal_strain_zz_aux
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [strain_xy]
    type = ElementAverageValue
    variable = strain_xy
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [strain_xz]
    type = ElementAverageValue
    variable = strain_xz
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [strain_yz]
    type = ElementAverageValue
    variable = strain_yz
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [stress_xx]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [stress_yy]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [stress_zz]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [stress_xy]
    type = ElementAverageValue
    variable = stress_xy
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [stress_xz]
    type = ElementAverageValue
    variable = stress_xz
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  [stress_yz]
    type = ElementAverageValue
    variable = stress_yz
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  []
  # Keep vonmises_stress on console for monitoring
  [vonmises_stress]
    type = ElementAverageValue
    variable = vonmises_stress
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_gmres_restart
  #                      -pc_hypre_boomeramg_nodal_coarsen
  #                      -pc_hypre_boomeramg_vec_interp_variant'
  # petsc_options_value  = 'gmres     hypre    boomeramg      201
                        # 1
                        # 1'
  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'preonly   lu       mumps'
  # petsc_options_iname = '-pc_type'
  # petsc_options_value  = 'lu'
  end_time = ${simulation_time}
  num_steps = 10
  nl_rel_tol = 5e-9
  # nl_abs_tol = 1e-12
  nl_max_its = 50
  l_max_its = 100
[]
