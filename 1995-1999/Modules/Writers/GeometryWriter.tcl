
catch {rename GeometryWriter ""}

itcl_class GeometryWriter {
    inherit Module
    constructor {config} {
	set name GeometryWriter
	set_defaults
    }
    method set_defaults {} {
	global $this-filetype
	set $this-filetype Binary
    }
    method ui {} {
	set w .ui$this
	if {[winfo exists $w]} {
	    raise $w
	    return;
	}
	toplevel $w

	make_labeled_radio $w.filetype "Format:" "" left $this-filetype \
		{Binary ASCII Gzip}
	pack $w.filetype
	entry $w.f -textvariable $this-filename -width 40 \
		-borderwidth 2 -relief sunken
	pack $w.f -side bottom
	bind $w.f <Return> "$this-c needexecute "
    }
}
