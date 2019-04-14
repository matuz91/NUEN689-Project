# This is a template for simulating nucleation in a phase-field model B
# Note that the fluctuations should still leave c conserved
# check this file for possible stress calc
## -> modules/combined/examples/phase_field-mechanics
### interfface_stress.i
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 64
  ny = 64
  xmin = 0
  xmax = 256
  ymin = 0
  ymax = 256
  uniform_refine = 2
  elem_type = QUAD4
  skip_partitioning = true
[]

# MIGHT USE THIS TO ADD NODESET AddExtraNodeset
#[MeshModifiers]
#  [./cnode]
#    type = AddExtraNodeset
#    coord = '0 0 0'
#    new_boundary = 100
#  [../]
#  [./anode]
#    type = AddExtraNodeset
#    coord = '0 10'
#    new_boundary = 101
#  [../]
#[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./c]
  [../]
  [./w]
  [../]
  [./eta]
  [../]
# added tensor mechanics
[./disp_x]
  order = FIRST
  family = LAGRANGE
[../]
[./disp_y]
  order = FIRST
  family = LAGRANGE
[../]
[]
# aux varaibles to track the free energy change (must decrease with time)
[AuxVariables]
  [./total_F]
    order = CONSTANT
    family = MONOMIAL
  [../]
# tensor mechanics
[./stress_xx]
  order = CONSTANT
  family = MONOMIAL
[../]
[./stress_yy]
  order = CONSTANT
  family = MONOMIAL
[../]
[./stress_zz]
  order = CONSTANT
  family = MONOMIAL
[../]
[./stress_xy]
  order = CONSTANT
  family = MONOMIAL
[../]
[./stress_yz]
  order = CONSTANT
  family = MONOMIAL
[../]
[./stress_zx]
  order = CONSTANT
  family = MONOMIAL
[../]
[]

[ICs]
  [./c_IC]
    type = ConstantIC
    variable = c
    value = 0.2
    # outside of the unstabel spinodal region (0.218<c<0.787) and inside the metastable region
  [../]
  [./IC_eta]
    x1 = 128
    y1 = 128
    radius = 64.0
    outvalue = 0.0
    variable = eta
    invalue = 1.0
    type = SmoothCircleIC
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
  [./left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

[Kernels]
  # Split form of Cahn-Hilliard equation
  # w is the chemical potential
  [./c_dot]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
  [./c_residual]
    type = SplitCHParsed
    variable = c
    f_name = F
    kappa_name = kappa_c
    w = w
  [../]
  [./w_residual]
    # args = 'c' in case the mobility is concentration dependent
    type = SplitCHWRes
    variable = w
    mob_name = M
  [../]
  [./conserved_langevin] # Adds fluctuations in the flux without changing the total concentration
    type = ConservedLangevinNoise
    amplitude = 0.02
    variable = w
    noise = normal_noise
  [../]
  # Allen-Cahn equation with c as coupled variable
  [./AC_bulk]
    type = AllenCahn
    variable = eta
    f_name = F
    args = c
    mob_name = L
  [../]
  [./AC_int]
    type = ACInterface
    variable = eta
    kappa_name = kappa_op
    mob_name = L
  [../]
  [./eta_dot]
    type = TimeDerivative
    variable = eta
  [../]
# tensor mechanics
[./TensorMechanics]
  use_displaced_mesh = true
[../]
[]

[AuxKernels]
  [./total_F]
    type = TotalFreeEnergy
    variable = total_F
    interfacial_vars = 'c eta'
    kappa_names = 'kappa_c kappa_op'
  [../]
# tensor mechanics
[./stress_xx]
  type = RankTwoAux
  rank_two_tensor = stress
  variable = stress_xx
  index_i = 0
  index_j = 0
[../]
[./stress_yy]
  type = RankTwoAux
  rank_two_tensor = stress
  variable = stress_yy
  index_i = 1
  index_j = 1
[../]
[./stress_zz]
  type = RankTwoAux
  rank_two_tensor = stress
  variable = stress_zz
  index_i = 2
  index_j = 2
[../]
[./stress_xy]
  type = RankTwoAux
  rank_two_tensor = stress
  variable = stress_xy
  index_i = 0
  index_j = 1
[../]
[./stress_yz]
  type = RankTwoAux
  rank_two_tensor = stress
  variable = stress_yz
  index_i = 1
  index_j = 2
[../]
[./stress_zx]
  type = RankTwoAux
  rank_two_tensor = stress
  variable = stress_zx
  index_i = 2
  index_j = 0
[../]
[]

[Materials]
  [./Bulk_Free_Eng]
    type = DerivativeParsedMaterial
    args ='eta c'
    f_name = F
    constant_names = 'A B D C_m C_p'
    # m:matrix, p: precipitate, C_p: Equilibrium concentration in precipitate
    # eta =1.0 inside the precipitate and eta=0.0 inside the matrix
    constant_expressions = '1.0 1.0 1.0 0.0 1.0'
    function = 'h:=(3.0*eta^2-2.0*eta^3);g_p:=h*A*(c-C_p)^2;g_m:=(1.0-h)*B*(c-C_m)^2;g_eta:=D*(eta^2*(eta-1.0)^2);g_m+g_p+g_eta'
    # h is an interpolation function
    derivative_order = 2
  [../]
  [./const]
    type = GenericConstantMaterial
    prop_names = 'kappa_c kappa_op L M'
    prop_values = '1.0 1.0 1.0 1.0'
  [../]
  [./precipitate_indicator]  # Returns 1 if inside the precipitate
    type = ParsedMaterial
    f_name = prec_indic
    args = c
    function = if(c>0.9,1.0,0)
  [../]
# tensor mechanics
[./elasticity_tensor]
  type = ComputeElasticityTensor
  block = '0'
  C_ijkl = '1.0e6  0.0   0.0 1.0e6  0.0  1.0e6 0.5e6 0.5e6 0.5e6'
  fill_method = symmetric9
[../]
[./strain]
  type = ComputeFiniteStrain
  block = '0'
  displacements = 'disp_x disp_y'
[../]
[./stress]
  type = ComputeFiniteStrainElasticStress
  block = '0'
[../]
[]

[Postprocessors]
  [./dofs]
    type = NumDOFs
  [../]
  [./_dt]
    # time step
    type = TimestepSize
  [../]
  [./nonlinear_its]
    # number of nonlinear iterations at each timestep
    type = NumNonlinearIterations
  [../]
  [./linear_its]
    # number of nonlinear iterations at each timestep
    type = NumLinearIterations
  [../]
  [./ElementInt_c]
    type = ElementIntegralVariablePostprocessor
    variable = c
  [../]
  [./total_F]
    type = ElementIntegralVariablePostprocessor
    variable = total_F
  [../]
  [./precipitate_area]      # Area of precipitate
    type = ElementIntegralMaterialProperty
    mat_prop = prec_indic
  [../]
  [./eta_precip_area] # since eta=1.0 inside the precipitate
    type = ElementIntegralVariablePostprocessor
    variable = 'eta'
  [../]
# tensor mechanics
[./stress_xx]
  type = ElementAverageValue
  variable = stress_xx
[../]
[./stress_yy]
  type = ElementAverageValue
  variable = stress_yy
[../]
[./stress_zz]
  type = ElementAverageValue
  variable = stress_zz
[../]
[./stress_xy]
  type = ElementAverageValue
  variable = stress_xy
[../]
[./stress_yz]
  type = ElementAverageValue
  variable = stress_yz
[../]
[./stress_zx]
  type = ElementAverageValue
  variable = stress_zx
[../]

[]

[Preconditioning]
  [./SMP] # to produce the complete perfect Jacobian
    type = SMP
    full = true
  [../]
[]
[UserObjects] # The object that creates the conserved fluctuation
  [./normal_noise]
    type = ConservedNormalNoise
  [../]
[]

[Executioner]
  type = Transient
  nl_max_its = 30
  scheme = bdf2
  solve_type = NEWTON
  petsc_options_iname = -pc_type
  petsc_options_value = asm
  l_max_its = 15
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-8
  start_time = 0.0
  num_steps = 200
  nl_abs_tol = 1e-9
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0
    growth_factor = 1.2
    cutback_factor = 0.75
    optimal_iterations = 6
  [../]
  [./Adaptivity]
    refine_fraction = 0.5
    coarsen_fraction = 0.01
    max_h_level = 3
    initial_adaptivity = 2
  [../]
[]

[Outputs]
  exodus = true
  csv = true
  interval = 1
  file_base = nucleation_model_b
[]
