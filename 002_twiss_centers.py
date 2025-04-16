import xtrack as xt
import numpy as np


from cpymad.madx import Madx
mad = Madx()
mad.input(f'''
    beam, particle=electron, PC = 150.0e-3;
    call, file = "general_tt43_python_no_foil_15bend.madx";
    use, sequence = TT43;
''')

ltt43 = xt.Line.from_madx_sequence(mad.sequence.tt43,
                                  deferred_expressions=True)
ltt43.particle_ref = xt.Particles(mass0=xt.ELECTRON_MASS_EV,
                                  gamma0=mad.sequence.tt43.beam.gamma)
env = ltt43.env

#############################
# Configure magnet modeling #
#############################

tt = ltt43.get_table()
tt_bend = tt.rows[(tt.element_type=='Bend') | (tt.element_type=='RBend')]
tt_quad = tt.rows[tt.element_type=='Quadrupole']
tt_sext = tt.rows[tt.element_type=='Sextupole']
tt_oct = tt.rows[tt.element_type=='Octupole']

ltt43['on_edge_mult'] = 0 # Use variable to control edge effects

# Bend model
ltt43.set(tt_bend, model='adaptive', integrator='uniform',
         num_multipole_kicks=100,
         edge_entry_model='full', edge_exit_model='full',
         edge_entry_active=True, edge_exit_active=True
         )
# Other magnets model
ltt43.set(tt_quad + tt_sext + tt_oct,
         model='drift-kick-drift-exact', integrator='yoshida4',
         num_multipole_kicks=20,
         edge_entry_model='full', edge_exit_model='full',
         edge_entry_active='on_edge_mult', edge_exit_active='on_edge_mult')
# Use exact drifts
ltt43.config.XTRACK_USE_EXACT_DRIFTS = True


tt_bpm = tt.rows['bpm.*']
tt_quad = tt.rows[tt.element_type=='Quadrupole']

# Insert markers at center of all quadrupoles and bpms
insertions = []
for nn in list(tt_quad.name) + list(tt_bpm.name):
    insertions.append(env.new(nn + '/center', xt.Marker, at=tt['s_center', nn]))

# You can do it on your line or on a shallow copy (shallow means that any strengths
# change on the main line is reflected on the copy)
line_with_centers = ltt43.copy(shallow=True)
line_with_centers.insert(insertions)

bxi = 11.0
axi = -2.1
byi = 11.0
ayi = -2.1
tw = line_with_centers.twiss(betx=bxi, alfx=axi, bety=byi, alfy=ayi)

tw.show()
# name                            s             x            px
# tt43$start                      0             0             0
# drift_0                         0             0             0
# lqm.0                       0.025             0             0
# drift_1                     0.025             0             0
# mqawd.0_entry                 0.1             0             0
# mqawd.0..entry_map            0.1             0             0
# mqawd.0..0                    0.1             0             0
# mqawd.0/center               0.25             0             0  # <-- quad center
# mqawd.0..1                   0.25             0             0
# mqawd.0..exit_map             0.4             0             0
# mqawd.0_exit                  0.4             0             0
# drift_2                       0.4             0             0
# lqp.0                       0.475             0             0
# drift_3                     0.475             0             0
# lbm.0                         0.5             0             0
# drift_4                       0.5             0             0
# bpm.0                         0.6             0             0
# bpm.0/center                  0.6             0             0 # <-- bpm center
# drift_5                       0.6             0             0
# ...

tw.rows['bpm.*/center']
# is:
# name                     s             x            px             y            py ...
# bpm.0/center           0.6             0             0             0             0
# bpm.1/center          2.19             0             0             0             0
# bpm.2/center         3.283             0             0             0             0
# bpm.3/center         5.306  -1.64398e-07  -3.00778e-07             0             0
# bpm.4/center         7.728  -2.16284e-07   1.02884e-07             0             0
# bpm.5/center       8.81391  -1.04561e-07   1.02883e-07             0             0
# bpm.6/center        10.124   3.34511e-08   1.09126e-07             0             0
# bpm.7/center       11.7175   2.07339e-07   1.09126e-07             0             0
# bpm.8/center        14.015   2.40379e-07  -3.05927e-07             0             0