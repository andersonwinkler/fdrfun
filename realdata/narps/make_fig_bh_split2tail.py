import gl
gl.resetdefaults()
gl.bmpzoom(2)
gl.loadimage("MNI152_T1_2mm_brain.nii")
gl.overlayload("narps-4965_9U7M-hypo1_unthresh.nii")
gl.minmax(1, 7, 7)
gl.colorfromzero(1,1)
gl.opacity(1,70)
gl.colorname (1,"6warm")
gl.overlayload("narps-4965_9U7M-hypo1_unthresh.nii")
gl.minmax(2, -7, -1.91344)
gl.colorfromzero(1,1)
gl.opacity(1,70)
gl.colorname (2,"5winter")
gl.colorbarposition(0)
gl.mosaic("A L- H 0.15 -12 0 12 24")
gl.savebmp("./bh_split2tail.png")
exit()
