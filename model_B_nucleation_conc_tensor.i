# This is a template for simulating nucleation in a phase-field model B
# Note that the fluctuations should still leave c conserved
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
[]

[Variables]
  [./c]
  [../]
  [./w]
  [../]
  [./eta]
  [../]
  # FROM ComputeConcentrationDependentElasticityTensor deck
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
##############################################################
[]
# aux varaibles to track the free energy change (must decrease with time)
[AuxVariables]
  [./total_F]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # FROM ComputeConcentrationDependentElasticityTensor deck
  [./C11_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s11_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
#############################################################
[]

[ICs]
# make these consistent, set equal to 0
# add body force and keep langevin noise the same
## DO THESE FIRST WITHOUT TensorMechanics
# elasticity at a later point and do concentration dependent
# Dr. Ahmed will send files at some point
  # [./c_IC]
  #   type = ConstantIC
  #   variable = c
  #   value = 0.2
  #   # outside of the unstabel spinodal region (0.218<c<0.787) and inside the metastable region
  # [../]
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
### newly added for tensor stuff
  [./left]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./bottom]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
#####################################
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
## newly added BodyForce
  [/BodyForce]
    type = BodyForce
    variable = c
    value = 10.0
  [../]
  [./TensorMechanics]
    displacements = 'disp_x disp_y' 
  [../]
[]

[AuxKernels]
  [./total_F]
    type = TotalFreeEnergy
    variable = total_F
    interfacial_vars = 'c eta'
    kappa_names = 'kappa_c kappa_op'
  [../]
  # FROM ComputeConcentrationDependentElasticityTensor deck
  [./matl_s11]
     type = RankTwoAux
     variable = s11_aux
     rank_two_tensor = stress
     index_i = 0
     index_j = 0
   [../]
   [./matl_C11]
     type = RankFourAux
     variable = C11_aux
     rank_four_tensor = elasticity_tensor
     index_l = 0
     index_j = 0
     index_k = 0
     index_i = 0
   [../]
###############################################################
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
# FROM ComputeConcentrationDependentElasticityTensor deck
  [./elasticity_tensor]
    type = ComputeConcentrationDependentElasticityTensor
    block = 0
    c = c
    C1_ijkl = '6 6'
    C0_ijkl = '1 1'
    fill_method1 = symmetric_isotropic
    fill_method0 = symmetric_isotropic
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 0
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
#######################################################
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
  [./stress_s11_aux]
    type = ElementAverageValue
    variable = s11_aux
  [../]
  [./elasticity_c11_aux]
    type = ElementAverageValue
    variable = C11_aux
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
  nl_max_its = 15
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
  file_base = nucleation_conc_tensor
[]
