### Mechanical Simulation of Copper Stabilizer in the same dimensions as the cable bundle in the ITER CS superconducting cable

# ITER Geometry
cable_length = '${units 0.45 m -> mm}'
# cable_radius = '${units 16.3 mm}'
cable_radius = '${units 10 mm}'

# EM Material Properties
vacuum_permeability = '${units 1.25663706e-6 N/A^2}'
current_density_z = '${units 1e6 A/m^2 -> A/mm^2}'
# current_density_z = 0 # Use to turn off Lorentz force

# Copper Material Properties (L.B. Freund, S. Suresh.  Thin Film Materials: Stress, Defect Formation, and Surface Evolution (2003)), Note: These are room temperature properties
youngs_modulus = '${units 130 GPa -> MPa}' # Room Temperature
poissons_ratio = 0.34 # Room Temperature
# thermal_expansion_coeff = '${units 1.657e-05 1/K}' # 293 K
# thermal_expansion_coeff = '${units 2.8e-09 1/K}' # 4.5 K
# thermal_expansion_coeff = '${units 4.473e-06 1/K}' # Effective coefficient

# Simulation Parameters

simulation_time = '${units 10 s}'

# Auxilliary
starting_temperature = '${units 4.5 K}'
ending_temperature = '${units 300 K}'
stress_free_temperature = '${starting_temperature}' # Starting with initial temperature for quench simulation


[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  # second_order = true
  [circle]
    type = ConcentricCircleMeshGenerator
    num_sectors = 8
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
    num_layers = '20'
    biases = '1'
    bottom_boundary = 'axial_start'
    top_boundary = 'axial_end'
    show_info = true
  []
  [pin_bottom_0_deg]
    type = ExtraNodesetGenerator
    input = stabilizer
    new_boundary = pin_bottom_0_deg
    coord = '${cable_radius} 0 0'
    use_closest_node=true
  []
  [pin_bottom_90_deg]
    type = ExtraNodesetGenerator
    input = pin_bottom_0_deg
    new_boundary = pin_bottom_90_deg
    coord = '0 ${cable_radius} 0'
    use_closest_node=true
  []
  [pin_bottom_180_deg]
    type = ExtraNodesetGenerator
    input = pin_bottom_90_deg
    new_boundary = pin_bottom_180_deg
    coord = '${fparse -cable_radius} 0 0'
    use_closest_node=true
  []
[]

[Physics/SolidMechanics]
  [QuasiStatic]
    [all]
      strain             = SMALL
      add_variables      = true
      eigenstrain_names = 'thermal_expansion'
      generate_output = 'vonmises_stress'
      # material_output_family = MONOMIAL
      # material_output_order = FIRST
      temperature = T
      use_automatic_differentiation = true
    []
  []
  [MaterialVectorBodyForce]
    [all]
      body_force = lorentz
    []
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${youngs_modulus}
    poissons_ratio = ${poissons_ratio}
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
  [expansion]
    type = ADComputeInstantaneousThermalExpansionFunctionEigenstrain
    temperature = T
    thermal_expansion_function = thermal_expansion_func
    stress_free_temperature = ${stress_free_temperature}
    eigenstrain_name = thermal_expansion
    outputs = 'exodus'
  []
  [lorentz]
    type = GenericFunctionVectorMaterial
    prop_names = lorentz
    prop_values = 'lorentz_x lorentz_y lorentz_z'
  []
  [thermal_expansion_coeff_functor]
    # Note: Expression duplicated from thermal_expansion_func to enable postprocessing with variable T
    # (Functions use 't' for time, but we need to evaluate with spatially-varying temperature field T)
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
  [lorentz_x_aux]
    family = LAGRANGE
    order = FIRST
  []
  [lorentz_y_aux]
    family = LAGRANGE
    order = FIRST
  []
  [lorentz_z_aux]
    family = LAGRANGE
    order = FIRST
  []
  [T]
    family = LAGRANGE
    order = FIRST
    initial_condition = '${starting_temperature}'
  []
  [thermal_strain_xx_aux]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [lorentz_x_extract]
    type = FunctionAux
    variable = lorentz_x_aux
    function = lorentz_x
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [lorentz_y_extract]
    type = FunctionAux
    variable = lorentz_y_aux
    function = lorentz_y
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [lorentz_z_extract]
    type = FunctionAux
    variable = lorentz_z_aux
    function = lorentz_z
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [temperature_ramp]
    type = FunctionAux
    variable = T
    function = 'temperature_func' #ramp up to final temperature over simulation time
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [thermal_strain_extract]
    type = ADRankTwoAux
    variable = thermal_strain_xx_aux
    rank_two_tensor = thermal_expansion
    index_i = 0
    index_j = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[BCs]
  [pin_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'pin_bottom_90_deg'
    value = 0
  []
  [pin_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'pin_bottom_0_deg pin_bottom_180_deg'
    value = 0
  []
  [pin_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'pin_bottom_0_deg pin_bottom_90_deg pin_bottom_180_deg'
    value = 0
  []
[]

[Functions]
  [jx]
    type = ParsedFunction
    expression = 0
  []
  [jy]
    type = ParsedFunction
    expression = 0
  []
  [jz]
    type = ParsedFunction
    expression = '${current_density_z}'
  []
  [bx]
    type = ParsedFunction
    expression = '-${vacuum_permeability}*jz*y/2'
    symbol_names = 'jz'
    symbol_values = 'jz'
  []
  [by]
    type = ParsedFunction
    expression = '${vacuum_permeability}*jz*x/2'
    symbol_names = 'jz'
    symbol_values = 'jz'
  []
  [bz]
    type = ParsedFunction
    expression = '0'
  []
  [lorentz_x]
    type = ParsedFunction
    expression = 'jy*bz - jz*by'
    symbol_names = 'jy bz jz by'
    symbol_values = 'jy bz jz by'
  []
  [lorentz_y]
    type = ParsedFunction
    expression = 'jz*bx - jx*bz'
    symbol_names = 'jz bx jx bz'
    symbol_values = 'jz bx jx bz'
  []
  [lorentz_z]
    type = ParsedFunction
    expression = 'jx*by - jy*bx'
    symbol_names = 'jx by jy bx'
    symbol_values = 'jx by jy bx'
  []

  [temperature_func]
    type = PiecewiseLinear
    x = '0 ${simulation_time}'
    y = '${starting_temperature} ${ending_temperature}'
  []

  [thermal_expansion_func] # t means Temperature when fed to Thermal Expansion Material
    type = ParsedFunction
    expression = '1e-6 * (10^(c0 + c1*log10(t) + c2*log10(t)^2 + c3*log10(t)^3 + c4*log10(t)^4 + c5*log10(t)^5 + c6*log10(t)^6))'
    symbol_names = 'c0 c1 c2 c3 c4 c5 c6'
    symbol_values = '-17.9081289 67.131914 -118.809316 109.9845997 -53.8696089 13.30247491 -1.30843441' # https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm
    # Note: Expression duplicated in thermal_expansion_coeff_functor material for postprocessing with variable T
  []
[]

[Postprocessors]
  [current_density_z]
    type = ConstantPostprocessor
    value = ${current_density_z}
    execute_on = 'INITIAL'
  []
  [vacuum_permeability]
    type = ConstantPostprocessor
    value = ${vacuum_permeability}
    execute_on = 'INITIAL'
  []
  [vonmises_max]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
  []

  [vonmises_average]
    type = ElementAverageValue
    variable = vonmises_stress
  []
  [temperature_average]
    type = ElementAverageValue
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [temperature_at_point]
    type = PointValue
    variable = T
    point = '${cable_radius} 0 ${fparse cable_length / 2}'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [thermal_expansion_coeff_pp]
    type = ElementAverageFunctorPostprocessor
    functor = thermal_expansion_coeff_functor
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [thermal_strain_xx_average]
    type = ElementAverageValue
    variable = thermal_strain_xx_aux
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[VectorPostprocessors]
  # Radial line at 0 degrees (along +x axis) at mid-height
  [line_sample_0deg]
    type = LineValueSampler
    start_point = '0 0 ${fparse cable_length / 2}'
    end_point = '${cable_radius} 0 ${fparse cable_length / 2}'
    num_points = 50
    variable = 'lorentz_x_aux lorentz_y_aux lorentz_z_aux'
    sort_by = id
  []
  # Radial line at 45 degrees at mid-height
  [line_sample_45deg]
    type = LineValueSampler
    start_point = '0 0 ${fparse cable_length / 2}'
    end_point = '${fparse cable_radius/sqrt(2)} ${fparse cable_radius/sqrt(2)} ${fparse cable_length / 2}'
    num_points = 50
    variable = 'lorentz_x_aux lorentz_y_aux lorentz_z_aux'
    sort_by = id
  []
  # Radial line at 90 degrees (along +y axis) at mid-height
  [line_sample_90deg]
    type = LineValueSampler
    start_point = '0 0 ${fparse cable_length / 2}'
    end_point = '0 ${cable_radius} ${fparse cable_length / 2}'
    num_points = 50
    variable = 'lorentz_x_aux lorentz_y_aux lorentz_z_aux'
    sort_by = id
  []
  # Axial line at constant radius to verify uniformity along z
  [line_sample_axial]
    type = LineValueSampler
    start_point = '${fparse cable_radius/2} 0 0'
    end_point = '${fparse cable_radius/2} 0 ${cable_length}'
    num_points = 50
    variable = 'lorentz_x_aux lorentz_y_aux lorentz_z_aux'
    sort_by = id
  []
[]

# consider all off-diagonal Jacobians for preconditioning
[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  end_time = ${simulation_time}
  num_steps = 20
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    file_base = 'data/copper_cylinder_out'
  []
[]
