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
thermal_expansion_coeff = '${units 1.657e-05 1/K}' # 293 K
# thermal_expansion_coeff = '${units 2.8e-09 1/K}' # 4.5 K
# thermal_expansion_coeff = '${units 4.473e-06 1/K}' # Effective coefficient

# Auxilliary
starting_temperature = '${units 293.15 K}'
stress_free_temperature = '${units 300 K}'


[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [circle]
    type = ConcentricCircleMeshGenerator
    num_sectors = 12
    has_outer_square = false
    preserve_volumes = true
    radii = '${cable_radius}'
    rings = '4'
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
      use_automatic_differentiation = true
      generate_output = 'vonmises_stress'
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
    type = ADComputeThermalExpansionEigenstrain
    temperature = T
    thermal_expansion_coeff = ${thermal_expansion_coeff}
    stress_free_temperature = ${stress_free_temperature}
    eigenstrain_name = thermal_expansion
  []
  [lorentz]
    type = GenericFunctionVectorMaterial
    prop_names = lorentz
    prop_values = 'lorentz_x lorentz_y lorentz_z'
  []
[]

[AuxVariables]
  [lorentz_x_aux]
    family = MONOMIAL
    order = FIRST
  []
  [lorentz_y_aux]
    family = MONOMIAL
    order = FIRST
  []
  [T]
    initial_condition = '${starting_temperature}'
  []
[]

[AuxKernels]
  [lorentz_x_extract]
    type = MaterialRealVectorValueAux
    variable = lorentz_x_aux
    property = lorentz
    component = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [lorentz_y_extract]
    type = MaterialRealVectorValueAux
    variable = lorentz_y_aux
    property = lorentz
    component = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [temperature_ramp]
    type = FunctionAux
    variable = T
    function = '${starting_temperature} + 5*t' #ramp up to 100 K over 10 seconds
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
    boundary = 'pin_bottom_0_deg pin_bottom_90_deg'
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
[]

[VectorPostprocessors]
  # Radial line at 0 degrees (along +x axis) at mid-height
  [line_sample_0deg]
    type = LineValueSampler
    start_point = '0 0 ${fparse cable_length / 2}'
    end_point = '${cable_radius} 0 ${fparse cable_length / 2}'
    num_points = 50
    variable = 'lorentz_x_aux lorentz_y_aux'
    sort_by = id
  []
  # Radial line at 45 degrees at mid-height
  [line_sample_45deg]
    type = LineValueSampler
    start_point = '0 0 ${fparse cable_length / 2}'
    end_point = '${fparse cable_radius/sqrt(2)} ${fparse cable_radius/sqrt(2)} ${fparse cable_length / 2}'
    num_points = 50
    variable = 'lorentz_x_aux lorentz_y_aux'
    sort_by = id
  []
  # Radial line at 90 degrees (along +y axis) at mid-height
  [line_sample_90deg]
    type = LineValueSampler
    start_point = '0 0 ${fparse cable_length / 2}'
    end_point = '0 ${cable_radius} ${fparse cable_length / 2}'
    num_points = 50
    variable = 'lorentz_x_aux lorentz_y_aux'
    sort_by = id
  []
  # Axial line at constant radius to verify uniformity along z
  [line_sample_axial]
    type = LineValueSampler
    start_point = '${fparse cable_radius/2} 0 0'
    end_point = '${fparse cable_radius/2} 0 ${cable_length}'
    num_points = 50
    variable = 'lorentz_x_aux lorentz_y_aux'
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
  end_time = 10
  dt = 1
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    file_base = 'data/copper_cylinder'
  []
[]
