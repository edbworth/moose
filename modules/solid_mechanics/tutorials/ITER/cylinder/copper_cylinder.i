### Full ITER Copper Cylinder Simulation
### Thermal expansion + Lorentz force with constrained BCs
### This is the general simulation with both physics active

!include copper_cylinder_base.i

# EM Parameters

# current_density_z = '${units 1e6 A/m^2 -> A/mm^2}'
current_density_z = 0

# EM Material Properties
# vacuum_permeability = '${units 1.25663706e-6 N/A^2}'
vacuum_permeability = 1.25663706e-6

# [Mesh]
#   # Add pin points for constrained boundary conditions
#   [pin_bottom_0_deg]
#     type = ExtraNodesetGenerator
#     input = stabilizer
#     new_boundary = pin_bottom_0_deg
#     coord = '${cable_radius} 0 0'
#     use_closest_node = true
#   []
#   [pin_bottom_90_deg]
#     type = ExtraNodesetGenerator
#     input = pin_bottom_0_deg
#     new_boundary = pin_bottom_90_deg
#     coord = '0 ${cable_radius} 0'
#     use_closest_node = true
#   []
#   [pin_bottom_180_deg]
#     type = ExtraNodesetGenerator
#     input = pin_bottom_90_deg
#     new_boundary = pin_bottom_180_deg
#     coord = '${fparse -cable_radius} 0 0'
#     use_closest_node = true
#   []
# []

[Kernels]
  [lorentz_x]
    type = ADBodyForce
    variable = disp_x
    function = lorentz_x
    use_displaced_mesh = true
  []
  [lorentz_y]
    type = ADBodyForce
    variable = disp_y
    function = lorentz_y
    use_displaced_mesh = true
  []
  [lorentz_z]
    type = ADBodyForce
    variable = disp_z
    function = lorentz_z
    use_displaced_mesh = true
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
[]

[BCs]
  [pin_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'pin_bottom_90_deg pin_bottom_center'
    value = 0
  []
  [pin_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'pin_bottom_180_deg pin_bottom_center'
    value = 0
  []
  [pin_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'pin_bottom_center pin_bottom_90_deg pin_bottom_180_deg'
    value = 0
  []
  # [pin_new]
  #   type = prese
[]

[Functions]
  # Electromagnetic field functions for analytical J×B force
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

[Postprocessors]
  [current_density_z]
    type = ConstantPostprocessor
    value = ${current_density_z}
    execute_on = 'INITIAL'
    outputs = 'csv'
  []
  [vacuum_permeability]
    type = ConstantPostprocessor
    value = ${vacuum_permeability}
    execute_on = 'INITIAL'
    outputs = 'csv'
  []
  [temperature_at_point]
    type = PointValue
    variable = T
    point = '${cable_radius} 0 ${fparse cable_length / 2}'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
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
  [bottom_face_disp]
    type = NodalValueSampler
    variable = 'disp_x disp_y disp_z'
    boundary = 'axial_start'
    sort_by = id
    outputs = 'csv'
  []
  [bottom_face_stress]
    type = NodalValueSampler
    variable = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz'
    boundary = 'axial_start'
    sort_by = id
    outputs = 'csv'
  []
[]

[Outputs]
  [exodus]
    type = Exodus
    file_base = 'data/copper_cylinder_out'
  []
  [csv]
    type = CSV
    file_base = 'data/copper_cylinder_out'
  []
[]
