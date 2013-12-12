from graphics import tango

shared_style = dict(linestyle='-')

phi_pivot_style = dict(color=tango.red[1], **shared_style)
phi_nonpivot_style = dict(color=tango.red[1], zorder=-10, linestyle=':')

psi_pivot_style = dict(color=tango.blue[1], **shared_style)
psi_nonpivot_style = dict(color=tango.blue[1], zorder=-10, linestyle=':')

omega_style = dict(color=tango.brown[1], **shared_style)

style = {
        'phi0':   dict(**phi_pivot_style),
        'phi1':   dict(**phi_pivot_style),
        'phi2':   dict(**phi_nonpivot_style),
        'phi3':   dict(**phi_pivot_style),

        'psi0':   dict(**psi_pivot_style),
        'psi1':   dict(**psi_pivot_style),
        'psi2':   dict(**psi_nonpivot_style),
        'psi3':   dict(**psi_pivot_style),

        'omega0':   dict(**omega_style),
        'omega1':   dict(**omega_style),
        'omega2':   dict(**omega_style),
        'omega3':   dict(**omega_style) }




