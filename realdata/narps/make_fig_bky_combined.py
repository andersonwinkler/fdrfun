import gl
gl.resetdefaults()
gl.backcolor(0, 0, 0)
gl.bmpzoom(2)
gl.loadimage("MNI152_T1_2mm_brain.nii")
gl.overlayloadsmooth(0)
gl.overlayload("narps-4965_9U7M-hypo1_unthresh.nii")
gl.minmax(1, 2.53531, 7)
gl.colorfromzero(1,1)
gl.opacity(1,100)
gl.colorname (1,"6warm")
gl.overlayloadsmooth(0)
gl.overlayload("narps-4965_9U7M-hypo1_unthresh.nii")
gl.minmax(2, -7, -2.53531)
gl.colorfromzero(2,1)
gl.opacity(2,100)
gl.colorname (2,"5winter")
gl.colorbarposition(0)
gl.mosaic("A L- H 0.15 -12 0 12 24")
gl.savebmp("./bky_combined.png")
exit()
