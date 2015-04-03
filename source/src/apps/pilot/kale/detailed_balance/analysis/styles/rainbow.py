from graphics import tango

pivot_line = dict(linestyle='-')
nonpivot_line = dict(linestyle='-')

style = {
        'phi0':   dict(color=tango.red[1], **pivot_line),
        'psi3':   dict(color=tango.red[0], **pivot_line),

        'psi0':   dict(color=tango.orange[1], **pivot_line),
        'phi3':   dict(color=tango.orange[0], **pivot_line),

        'omega0':   dict(color=tango.green[2], **nonpivot_line),
        'omega2':   dict(color=tango.green[1], **nonpivot_line),

        'phi1':   dict(color=tango.blue[1], **pivot_line),
        'psi2':   dict(color=tango.blue[0], **nonpivot_line),

        'psi1':   dict(color=tango.purple[1], **pivot_line),
        'phi2':   dict(color=tango.purple[0], **nonpivot_line),

        'omega1':   dict(color=tango.brown[1], **nonpivot_line) }


