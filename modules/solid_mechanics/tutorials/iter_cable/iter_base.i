### Mechanical simulation of a section of ITER cable steel jacket

!include iter_cable.params

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics]
  [QuasiStatic]
    [all]
      strain             = SMALL
      add_variables      = true
      use_automatic_differentiation = true
    []
  []
  [MaterialVectorBodyForce]
    [all]
      body_force = lorentz
    []
  []
[]

[BCs]
  [cable_joints_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'front back'
    value = 0
  []

[]

[Functions]
  [jx]
    type = ParsedFunction
    expression = 't*1e8'
  []
  [jy]
    type = ParsedFunction
    expression = 't*1e8'
  []
  [jz]
    type = ParsedFunction
    expression = '0'
  []
  [bx]
    type = ParsedFunction
    expression = '0'
  []
  [by]
    type = ParsedFunction
    expression = '0'
  []
  [bz]
    type = ParsedFunction
    expression = 't*1'
  []
  [lorentz_x]
    type = ParsedFunction
    expression = '${unit_scale}*(jy*bz - jz*by)'
    symbol_names = 'jy bz jz by'
    symbol_values = 'jy bz jz by'
  []
  [lorentz_y]
    type = ParsedFunction
    expression = '${unit_scale}*(jz*bx - jx*bz)'
    symbol_names = 'jz bx jx bz'
    symbol_values = 'jz bx jx bz'
  []
  [lorentz_z]
    type = ParsedFunction
    expression = '${unit_scale}*(jx*by - jy*bx)'
    symbol_names = 'jx by jy bx'
    symbol_values = 'jx by jy bx'
  []

[]

[Materials] # These material properties can be refined further for JK2LB steel, but are comparable to traditional stainless steel. Note: These are room temperature properties
  [elasticity] # https://www.sciencedirect.com/science/article/pii/S0011227508001835?via%3Dihub
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = '${youngs_modulus}'
    poissons_ratio = '${poissons_ratio}'
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
  [lorentz]
    type = GenericFunctionVectorMaterial
    prop_names = lorentz
    prop_values = 'lorentz_x lorentz_y lorentz_z'
  []
[]

# consider all off-diagonal Jacobians for preconditioning
[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

# [Debug]
#   show_material_props = true
#   show_parser = true
# []

[Executioner]
  type = Transient # Note that the model equations are not time dependent, but uses transient executioner to help solvers tackle strong nonlinearities
  # we chose a direct solver here
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  end_time = 1
  dt = 0.2
[]

[Outputs]
  exodus = true
[]
