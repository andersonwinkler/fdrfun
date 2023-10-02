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
gl.overlaymaskwithbackground(0)
gl.mosaic("A L- V -0.1 -32 -16 0 8 16 32 40 50 S X R 0")
gl.savebmp("./bky_split2tail.png")
exit()
