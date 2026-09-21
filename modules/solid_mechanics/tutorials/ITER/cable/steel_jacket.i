# Mechanical Simulation of JK2LB steel jacket of ITER cable

!include iter_base.i

[Mesh]
  [jacket_tmp]
    type = ConcentricCircleMeshGenerator
    has_outer_square = true
    num_sectors = 6
    preserve_volumes = true
    radii = '${cable_radius}'
    rings = '1 2'
    pitch = '${steel_jacket_side_length}'
    smoothing_max_it = 3
  []

  [name_blocks]
    type = RenameBlockGenerator
    input = jacket_tmp
    old_block = '1 2'
    new_block = 'cable_void jacket'
  []

  [delete_cable_void]
    type = BlockDeletionGenerator
    input = name_blocks
    block = 'cable_void'
    new_boundary = 'jacket_inner_wall'
  []

  [extruded_mesh]
    type = MeshExtruderGenerator
    extrusion_vector = '0 0 ${cable_length}'
    input = delete_cable_void
    num_layers = 20
    bottom_sideset = 'bottom'
    top_sideset = 'top'
    show_info = true
  []
[]
