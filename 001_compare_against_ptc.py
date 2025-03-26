from cpymad.madx import Madx
import xtrack as xt

flag_fringe = 'true'
mad_dct = {}
tt_dct = {}
tt0_dct = {}
tw_dct = {}
tt_all_dct = {}
for flag_fringe in ['true', 'false']:
    mad = Madx()
    mad.input(f'''
    !********************************************
    ! TT43-TT41 (AWAKE e-line) model
    !
    ! E. Belli
    !********************************************

    TITLE, "AWAKE TT43 electron line";
    set, format="22.10e";

    /***************************************
    * Cleaning output files
    ***************************************/

    system, "rm -r output_files";
    system, "mkdir output_files";

    /**********************************************
    * Initial beam parameters
    **********************************************/
    !Beam energy (MeV)
    pc = 150.0e-3 ;

    ! Initial H&V beta, alpha, emittances
    bxi = 11;
    axi = -2.1;
    byi = 11;
    ayi = -2.1;
    emitx_n = 2e-6;
    emity_n = 2e-6;
    dpp = 0.002;




    /**********************************************
    * Beam and initial lattice parameters
    **********************************************/

    beam, particle=electron, PC = pc, EXN=emitx_n, EYN=emity_n;

    initial: BETA0, BETX = bxi, ALFX = axi, MUX = 0.0,
                    BETY = byi, ALFY = ayi, MUY = 0.0,
                    DX = 0.0, DPX = 0.0, DY = 0.0, DPY = 0.0,
                    X=0.0, PX = 0.0, Y=0.0, PY=0.0;


    /**********************************************
    * Sequence
    **********************************************/
    call, file = "general_tt43_python_no_foil_15bend.madx";
    use, sequence = TT43;


    !-------------------------------------- Tracking in the transfer line -----------------------!


    ptc_create_universe;
    ptc_create_layout, model=1, method=6, exact=True, NST=100;
    ptc_align;
    ptc_setswitch,fringe='''
    f'{flag_fringe}'
    ''';
    call, file = 'distribution_for_tracking_ptc.dat';
    ptc_observe,place='MERGE';
    ptc_track, icase=56, turns = 1, dump, onetable, element_by_element, recloss=true, closed_orbit=false, maxaper={0.025, 0.025, 0.025, 0.025, 1.0, 1},file="output_files/tracking_";

    ptc_twiss, icase=56, betx=1., bety=1., betz=1,x=-7.2583e-05, px=-5.3961e-05, y=-0.000178864, py=-3.3722e-05, t=0.0001142302, pt=0.0026680501099999995;
    ptc_track_end;

    ptc_end;
    ''')

    tw = xt.Table(mad.table.ptc_twiss)

    tt_all = xt.Table(mad.table.trackone, index='number')

    tt = tt_all.rows[15.3:15.35:'s']
    tt0 = tt_all.rows[0.0:0.05:'s']

    mad_dct[flag_fringe] = mad
    tt_dct[flag_fringe] = tt
    tt0_dct[flag_fringe] = tt0
    tw_dct[flag_fringe] = tw
    tt_all_dct[flag_fringe] = tt_all

tt_fringe = tt_dct['true']
tt_no_fringe = tt_dct['false']

tt0_fringe = tt0_dct['true']
tt0_no_fringe = tt0_dct['false']

tw_fringe = tw_dct['true']
tw_no_fringe = tw_dct['false']

tt_all_fringe = tt_all_dct['true']
tt_all_no_fringe = tt_all_dct['false']


tw_fringe.rows["merge:1"].cols["s y py"].show()
tt_all_fringe.rows[tt_all_fringe.number==8565,1].cols["s y py"].show()
tw_no_fringe.rows["merge:1"].cols["s y py"].show()
tt_all_no_fringe.rows[tt_all_no_fringe.number==8565,1].cols["s y py"].show()

tt43 = xt.Line.from_madx_sequence(mad_dct['true'].sequence.tt43,
                                  deferred_expressions=True)

env = tt43.env
line = tt43

tt = line.get_table()
tt_bend = tt.rows[(tt.element_type=='Bend') | (tt.element_type=='RBend')]
tt_quad = tt.rows[tt.element_type=='Quadrupole']
tt_sext = tt.rows[tt.element_type=='Sextupole']
tt_oct = tt.rows[tt.element_type=='Octupole']

line['on_edge_mult'] = 0 # Use variable to control edge effects

line.set(tt_bend, model='adaptive', integrator='uniform',
         num_multipole_kicks=100,
         edge_entry_model='full', edge_exit_model='full',
         edge_entry_active=True, edge_exit_active=True
         )
line.set(tt_quad + tt_sext + tt_oct,
         model='drift-kick-drift-exact', integrator='uniform',
         num_multipole_kicks=10,
         edge_entry_model='full', edge_exit_model='full',
         edge_entry_active='on_edge_mult', edge_exit_active='on_edge_mult')
line.config.XTRACK_USE_EXACT_DRIFTS = True
# Include aperture in the line for particle losses during tracking

import numpy as np

env.elements['aperture20'] = xt.LimitEllipse(a=2e-2, b=2e-2)
env.elements['aperture30'] = xt.LimitEllipse(a=3e-2, b=3e-2)
env.elements['aperture40'] = xt.LimitEllipse(a=4e-2, b=4e-2)
step = 0.001
s_apertures = np.arange(0, tt['s', -1], step)
insertions = []

tt = tt43.get_table()

for s in s_apertures:
    if ((s < tt.rows[tt.name==['mbawh.8']].s_end) and (s > tt.rows[tt.name==['mbawh.8']].s_start) or 
        (s < tt.rows[tt.name==['mbawh.3']].s_end) and (s > tt.rows[tt.name==['mbawh.3']].s_start)):
        insertions.append(env.place('aperture40', at=s)) # aperture in mm
    elif ((s < tt.rows[tt.name==['mqawd.9']].s_end) and (s > tt.rows[tt.name==['mqawd.9']].s_start) or 
          (s < tt.rows[tt.name==['mqawd.6']].s_end) and (s > tt.rows[tt.name==['mqawd.6']].s_start) or 
          (s < tt.rows[tt.name==['sd3']].s_end) and (s > tt.rows[tt.name==['sd3']].s_start) or 
          (s < tt.rows[tt.name==['sd1']].s_end) and (s > tt.rows[tt.name==['sd1']].s_start) or 
          (s < tt.rows[tt.name==['btv.4']].s_end) and (s > tt.rows[tt.name==['btv.4']].s_start)
          ):
        insertions.append(env.place('aperture30', at=s))
    else:
        insertions.append(env.place('aperture20', at=s))

tt43.insert(insertions)

line = tt43.select(end='merge')


p = xt.Particles(
    gamma0=mad.sequence.tt43.beam.gamma,
    mass0=xt.ELECTRON_MASS_EV,
    x=tt0.x, px=tt0.px, y=tt0.y, py=tt0.py,
    tau=tt0.t, ptau=tt0.pt)

p.to_json('particles_init.json')

p_no_fringe = p.copy()
line.track(p_no_fringe)

line['on_edge_mult'] = 1

p_fringe = p.copy()
line.track(p_fringe)

p_fringe_at_merge = p_fringe.filter(p_fringe.state > 0)
p_no_fringe_at_merge = p_no_fringe.filter(p_no_fringe.state > 0)



import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
ax1 = plt.subplot(111)
plt.plot(tt_fringe.y*1e3, tt_fringe.py*1e3, '.')
plt.plot(tt_no_fringe.y*1e3, tt_no_fringe.py*1e3, '.')
plt.xlabel('y [m]')
plt.ylabel('py [mrad]')
plt.title('PTC tracking')

plt.figure(2)
ax2 = plt.subplot(111, sharex=ax1, sharey=ax1)
plt.plot(p_fringe_at_merge.y*1e3, p_fringe_at_merge.py*1e3, '.')
plt.plot(p_no_fringe_at_merge.y*1e3, p_no_fringe_at_merge.py*1e3, '.')
plt.xlabel('y [m]')
plt.ylabel('py [mrad]')
plt.title('Xsuite tracking')

plt.show()

