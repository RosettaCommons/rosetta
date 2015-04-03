from graphics import tango

pivot_line = dict(linestyle='-')
nonpivot_line = dict(linestyle='-')

style = {
        'phi0':   dict(color=tango.red[1], **pivot_line),
        'psi3':   dict(color=tango.red[1], **pivot_line),

        'psi0':   dict(color=tango.red[1], **pivot_line),
        'phi3':   dict(color=tango.red[1], **pivot_line),

        'phi1':   dict(color=tango.red[1], **pivot_line),
        'psi2':   dict(color=tango.red[1], **nonpivot_line),

        'psi1':   dict(color=tango.red[1], **pivot_line),
        'phi2':   dict(color=tango.red[1], **nonpivot_line) }


