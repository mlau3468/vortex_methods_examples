from dust.tools.DUST_tools import *
from time import time
if __name__=='__main__':
    # Create DUST model from VSP degenerate geometry
    geo_files, geo_names, mesh_stats = vspDegen2DustBasic('wing_DegenGeom.csv', 'Model/', -0.75, show=False, dup_chck='Full', scale=1/3.28)

    # ---------------------------------------------------------------------------------

    # Flight Conditions
    dt = 0.005
    t_steps = 100
    a = 0
    b = 0
    v = 30

    period = dt * t_steps
    vx = v*math.cos(math.radians(b))*math.cos(math.radians(a))
    vz = v*math.cos(math.radians(b))*math.sin(math.radians(a))
    vy = v*math.sin(math.radians(b))

    # ---------------------------------------------------------------------------------
    # Dust Settings
    a = Runner(out_dir='dust_out/', verbosity=1)

    #dust.in
    options = {}
    options['basename'] = './Output/alia'
    options['basename_debug'] = './Debug/alia'
    options['u_inf'] = [vx, vy, vz]
    options['tend'] = period
    options['dt'] = 0.005
    options['dt_out'] = 0.005
    options['dt_debug_out'] = 0.005
    options['n_wake_panels'] = 1
    options['n_wake_particles'] = 300000
    options['particles_box_min'] = [-10, -15, -15]
    options['particles_box_max'] = [50, 15, 15]
    options['FMM'] = False
    options['output_start'] = True
    options['LLdamp'] = 15
    options['BoxLength'] = 10
    options['NBox'] = [6, 3, 3]
    options['OctreeOrigin'] = [-10, -15, -15]
    options['FarFieldRatioSource'] = 1000
    options['FarFieldRatioDoublet'] = 1000
    options['debug_level'] = 50
    options['NOctreeLevels'] = 6
    options['MinOctreePart'] = 5
    options['MultipoleDegree'] = 2
    options['join_te'] = True
    options['Vortstretch'] = False
    options['Diffusion'] = False
    options['PenetrationAvoidance'] = False

    # Dust references.in
    new_ref = {'Reference_Tag': 'test', 'Parent_Tag' : '0', 'Origin': [0,0,0,], 'Orientation': [0.9961947,  0.0000000,  -0.0871557, 0.0000000,  1.0000000,  0.0000000, 0.0871557,  0.0000000,  0.9961947]}
    a.add_ref( new_ref, 'simple',)
    new_ref = {'Reference_Tag': 'test2', 'Parent_Tag' : 'test', 'Origin': [0,0,0,], 'Orientation': [1.0,  0.0,  0.0, 0.0,0.0,1.0, 0.0,  -1,  0.0]}
    a.add_ref(new_ref, 'simple', )

    # Components
    for i, g in enumerate(geo_files):
        a.add_comp(geo_names[i], g, 'test')

    # ---------------------------------------------------------------------------------
    # Convergence Settings
    conv_options = {}
    conv_options['conv_type'] = 'steady'
    conv_options['max_runs'] = 10
    conv_options['tol'] = 0.05
    conv_options['scheme'] = 'max'
    conv_options['tstep_num'] = 10
    conv_options['comp_axes'] = [[0,0,1], [0,0,1]]
    conv_options['comp_frames'] = ['local', 'local']

    # ---------------------------------------------------------------------------------
    # Run Dust
    a.set_dust_options(options)
    a.set_conv_options(conv_options)
    ans = a.run(conv_error=True)
    a.run_vis('./Postprocessing/alia')

    # ---------------------------------------------------------------------------------

    # Process Results
    case = ans['CaseObj']
    #case.show()
    print('Loads on Comp1:, reference frames: global, test (same as local), test2')
    print(case.components[1].get_loads(frame='global'))
    print(case.components[1].get_loads(frame='test'))
    print(case.components[1].get_loads(frame='test2'))

    print('Loads on all components')
    print(case.get_loads(comp_num=None, frame='global'))

    print('See output files. Spanwise loads on Comp1: reference frames: global, test (same as local) ')
    case.components[1].span_loads(span_locs='auto', loads_frame='global', filename='Analysis/span_loads_global')
    case.components[1].span_loads(span_locs='auto', loads_frame='local', filename='Analysis/span_loads_local')
    '''
    print('Probe a point to get data')
    print(case.components[1].probe([0,0,0], frame='local'))
    print('See output files. Data for all panels on this component')
    case.components[1].get_data(frame='global')
    case.components[1].get_data(frame='local')
    '''
    
    print('See output files. Chordwise pressure')
    case.components[1].p_line(2, axis='chord', filename='Analysis/chord_p')
    print('See output files. Chordwise Velocity')
    case.components[1].surf_vel(2, frame='local', proj_vec=None, axis='chord', filename='Analysis/chord_surfvel')
    
    print(case.get_comp_names())

