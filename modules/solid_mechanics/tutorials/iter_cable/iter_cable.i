### Electro-thermal simulation of a section of ITER central solenoid cable

!include iter_cable.params

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [cable_2D]
    type = ConcentricCircleMeshGenerator
    has_outer_square = true
    num_sectors = 12
    preserve_volumes = true
    radii = '${channel_radius} ${cable_radius}'
    rings = '2 4 8'
    pitch = '${steel_jacket_side_length}'
    smoothing_max_it = 3
  []

  [name_blocks]
    type = RenameBlockGenerator
    input = cable_2D
    old_block = '1 2 3'
    new_block = 'channel cable jacket'
  []

  [delete_helium]
    type = BlockDeletionGenerator
    input = name_blocks
    block = 'channel'
    new_boundary = 'cable_helium_boundary'
  []

  # [cable_3D_coarse]
  #   type = AdvancedExtruderGenerator
  #   input = delete_helium
  #   direction = '0 0 1'
  #   heights = '${cable_length}'
  #   num_layers = '20'
  #   biases = '1'
  #   bottom_boundary = 'axial_start'
  #   top_boundary = 'axial_end'
  #   show_info = true
  # []

  [cable_3D]
  type = AdvancedExtruderGenerator
  input = delete_helium
  direction = '0 0 1'
  heights = '${fparse cable_length / 2} ${fparse cable_length / 2}'
  num_layers = '20 20'
  # Element sizes decrease toward the center, then increase away from it
  biases = '0.95 1.05'
  bottom_boundary = 'axial_start'
  top_boundary = 'axial_end'
  show_info = true
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain             = SMALL
    add_variables      = true
    use_automatic_differentiation = true
  []
[]

[BCs]
  [cable_joints_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'axial_start axial_end'
    value = 0
  []
  [cable_joints_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'axial_start axial_end'
    value = 0
  []
  [cable_joints_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'axial_start axial_end'
    value = 0
  []
  [Pressure]
    [left]
      boundary = 'left'
      function = 5e9*t
      use_automatic_differentiation = true
    []
  []
  # [fix_wall_x]
  #   type = ADDirichletBC
  #   variable = disp_x
  #   boundary = 'jacket_inner_wall'
  #   value = 0
  # []
  # [fix_wall_y]
  #   type = ADDirichletBC
  #   variable = disp_y
  #   boundary = 'jacket_inner_wall'
  #   value = 0
  # []
  # [fix_wall_z]
  #   type = ADDirichletBC
  #   variable = disp_z
  #   boundary = 'jacket_inner_wall'
  #   value = 0
  # []
[]

[Materials] # These material properties can be refined further for JK2LB steel, but are comparable to traditional stainless steel. Note: These are room temperature properties
  [elasticity] # https://www.sciencedirect.com/science/article/pii/S0011227508001835?via%3Dihub
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${youngs_modulus}
    poissons_ratio = ${poissons_ratio}
  []
  [stress]
    type = ADComputeLinearElasticStress
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
  type = Transient # Note that the model equations are not time dependent, but uses transient executioner to help solvers tackle strong nonlinearities
  # we chose a direct solver here
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  end_time = 5
  dt = 1
[]

[Outputs]
  exodus = true
[]

