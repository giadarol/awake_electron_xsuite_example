import xtrack as xt
import numpy as np

# env = xt.Environment()
# env.call('tt43_150MeV_Run2c.py')

from cpymad.madx import Madx
mad = Madx()
mad.input(f'''
    beam, particle=electron, PC = 150.0e-3;
    call, file = "general_tt43_python_no_foil_15bend.madx";
    use, sequence = TT43;
''')

ltt43 = xt.Line.from_madx_sequence(mad.sequence.tt43,
                                  deferred_expressions=True)
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

#####################
# Install apertures #
#####################

env.elements['aperture20'] = xt.LimitEllipse(a=2e-2, b=2e-2)
env.elements['aperture30'] = xt.LimitEllipse(a=3e-2, b=3e-2)
env.elements['aperture40'] = xt.LimitEllipse(a=4e-2, b=4e-2)
step = 0.001
s_apertures = np.arange(0, tt['s', -1], step)
insertions = []

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

ltt43.insert(insertions)

# Select part to track
line = ltt43.select(end='merge')

############
# Tracking #
############

p = xt.Particles.from_dict(xt.json.load('particles_init.json'))

p_no_fringe = p.copy()
print('Tracking without fringe fields')
line.track(p_no_fringe)
print('Done')

line['on_edge_mult'] = 1

p_fringe = p.copy()
print('Tracking with fringe fields')
line.track(p_fringe)
print('Done')

p_fringe_at_merge = p_fringe.filter(p_fringe.state > 0)
p_no_fringe_at_merge = p_no_fringe.filter(p_no_fringe.state > 0)

########
# Plot #
########

import matplotlib.pyplot as plt
plt.close('all')

plt.figure(2)
ax2 = plt.subplot(111)
plt.plot(p_fringe_at_merge.y*1e3, p_fringe_at_merge.py*1e3, '.')
plt.plot(p_no_fringe_at_merge.y*1e3, p_no_fringe_at_merge.py*1e3, '.')
plt.xlabel('y [m]')
plt.ylabel('py [mrad]')
plt.title('Xsuite tracking')

plt.show()