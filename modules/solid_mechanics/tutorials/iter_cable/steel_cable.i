# Simulation of Solid Steel Cable: Iter Cable if inner cable and cooling channel are all JK2LB steel

!include iter_base.i

[Mesh]
  [full_steel_cable]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 100
    xmax = ${steel_jacket_side_length}
    ymax = ${steel_jacket_side_length}
    zmax = ${cable_length}
  []
[]

