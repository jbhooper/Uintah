
itcl_class Uintah_Visualization_ParticleEigenEvaluator {
    inherit Module

    constructor {config} {
	set name ParticleEigenEvaluator
	set_defaults
    }

    method set_defaults {} {
	global $this-eigenSelect
	set $this-eigenSelect 0
    }
    method ui {} {
	set w .ui[modname]
	if {[winfo exists $w]} {
	    raise $w
	    return;
	}
	toplevel $w
	set n "$this-c needexecute"

	frame $w.m 
	pack $w.m -padx 2 -pady 2 -fill x -expand yes

	make_labeled_radio $w.m.r "Which Eigen Value" $n left \
		$this-eigenSelect {{"Largest" 0} {"Middle" 1} {"Smallest" 2}}
	pack $w.m.r

	button $w.b -text Close -command "destroy $w"
	pack $w.b -side bottom -expand yes -fill x -padx 2 -pady 2
    }
}

